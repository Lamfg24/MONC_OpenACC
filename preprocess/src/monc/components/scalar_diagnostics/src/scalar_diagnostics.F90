










module scalar_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type, &
       component_descriptor_type_v1
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use state_mod, only : FORWARD_STEPPING
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use science_constants_mod, only : cp, rlvap
  use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, LOG_DEBUG, log_master_log

  implicit none

  private

  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0,    &
       iqg=0
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: dz_rhon_fac
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ww_prime_res, uu_prime_res, &
       vv_prime_res, cloud_content
  real(kind=DEFAULT_PRECISION) :: qlcrit
  public initialisation_callback_scalar_diagnostics, timestep_callback_scalar_diagnostics

contains

  subroutine initialisation_callback_scalar_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, iter
    integer :: y_size_local, x_size_local
    real(kind=DEFAULT_PRECISION) :: realnum
    
    ! qlcrit declared in the config file
    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "qlcrit") then
        read(current_state%options_database_string(iter,2),*) realnum
        qlcrit = realnum
      end if
    end do

    y_size_local = current_state%local_grid%size(Y_INDEX)
    x_size_local = current_state%local_grid%size(X_INDEX)
    
    ! allocate scalar diagnostics as 2-D fields (horizontal slices) that 
    ! are published to the ioserver so that the ioserver can do the manipulation
    ! to obtain the scalar field
    allocate(current_state%wmax(y_size_local, x_size_local), current_state%wmin(y_size_local, x_size_local), &
         current_state%reske(y_size_local, x_size_local),                                      &
         current_state%senhf(y_size_local, x_size_local),current_state%lathf(y_size_local, x_size_local))
    ! allocate the 1d velocity arrays for the kinetic energy calc
    allocate(ww_prime_res(current_state%local_grid%size(Z_INDEX)), &
         uu_prime_res(current_state%local_grid%size(Z_INDEX)),     &
         vv_prime_res(current_state%local_grid%size(Z_INDEX)))
    ww_prime_res(:) = 0.0
    uu_prime_res(:) = 0.0
    vv_prime_res(:) = 0.0

    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
       iqv = current_state%water_vapour_mixing_ratio_index
       iql = current_state%liquid_water_mixing_ratio_index
       allocate(current_state%vwp(y_size_local, x_size_local), current_state%lwp(y_size_local, x_size_local), &
            current_state%qlmax(y_size_local, x_size_local), current_state%hqlmax(y_size_local, x_size_local), &
            current_state%cltop(y_size_local, x_size_local), current_state%clbas(y_size_local, x_size_local))
       allocate(cloud_content(current_state%local_grid%size(Z_INDEX)))
       ! allocate other hydrometeors. Allocation dependent on index being set in
       ! appropriate microphysics scheme (see casim component from example)
       if ((current_state%casim_enabled .eqv. .false.) .and. (current_state%simplecloud_enabled .eqv. .false.)) then
         call log_master_log(LOG_WARN, "CASIM or simplecloud should be activated to calculate hydro species scalar diagnostics")
       else
        if (current_state%rain_water_mixing_ratio_index > 0) then
            iqr = current_state%rain_water_mixing_ratio_index
            allocate(current_state%rwp(y_size_local, x_size_local))
        endif
        if (current_state%ice_water_mixing_ratio_index > 0) then
            iqi = current_state%ice_water_mixing_ratio_index
            allocate(current_state%iwp(y_size_local, x_size_local))
            allocate(current_state%tot_iwp(y_size_local, x_size_local))
        endif
        if (current_state%snow_water_mixing_ratio_index > 0) then
            iqs = current_state%snow_water_mixing_ratio_index
            allocate(current_state%swp(y_size_local, x_size_local))
        endif
        if (current_state%graupel_water_mixing_ratio_index > 0) then
            iqg = current_state%graupel_water_mixing_ratio_index
            allocate(current_state%gwp(y_size_local, x_size_local))
        endif
      end if
    endif
    
    ! Surface Flux of tracers
    if (current_state%n_tracers .gt. 0) then
      allocate(current_state%trsfflux(y_size_local, x_size_local, current_state%n_tracers))
    endif
    
    allocate(dz_rhon_fac(current_state%local_grid%size(Z_INDEX)))    
    do k=2, current_state%local_grid%size(Z_INDEX)
       ! used in the water path calculation
       dz_rhon_fac(k)=current_state%global_grid%configuration%vertical%dz(k)*&
            current_state%global_grid%configuration%vertical%rhon(k)
    end do    

  end subroutine initialisation_callback_scalar_diagnostics

  subroutine timestep_callback_scalar_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i, n
    integer :: current_y_index, current_x_index, target_x_index, target_y_index

    if ((current_state%modulo_number_0d .eq. 0) .or. (current_state%modulo_number_2d .eq. 0)) then

    current_y_index=current_state%column_local_y
    current_x_index=current_state%column_local_x
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    if (current_state%first_timestep_column) then
       ! maximum vertical velocity for each column
       current_state%wmax(:,:)=0.0
       ! minimum vertical velocity for each column
       current_state%wmin(:,:)=0.0
       ! resolved ke
       current_state%reske(:,:) = 0.0


        if (.not. current_state%passive_q) then! .and. current_state%number_q_fields .gt. 0) then
          ! maximum liquid water content in a column
          current_state%qlmax(:,:)=0.0
          ! the height of the maximum liquid water content in a column
          current_state%hqlmax(:,:)=0.0
          ! cloud top height where liqud water content is greater than qlcrit
          current_state%cltop(:,:)=0.0
          ! minimum cloud base where liquid water content is greater than qlcrit
          current_state%clbas(:,:)=0.0
          ! water vapour path for each column
          current_state%vwp(:,:)=0.0
          ! liquid water path for each column
          current_state%lwp(:,:)=0.0
          ! rain water path for each column
          if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
            if (current_state%rain_water_mixing_ratio_index > 0) current_state%rwp(:,:)=0.0
            ! ice water path for each column
            if (current_state%ice_water_mixing_ratio_index > 0) then
              current_state%iwp(:,:)=0.0
              ! total ice water path (iwp + swp + gwp) for each column
              current_state%tot_iwp(:,:)=0.0
            endif
            ! snow water path for each column
            if (current_state%snow_water_mixing_ratio_index > 0) current_state%swp(:,:)=0.0
            ! graupel water path for each column
            if (current_state%graupel_water_mixing_ratio_index > 0) current_state%gwp(:,:)=0.0
          end if
       end if
       ! surface sensible heat flux
       current_state%senhf(:,:)=0.0
       ! surface latent heat flux
       current_state%lathf(:,:)=0.0
       
       ! Surface Flux of radioactive tracers
       if (current_state%n_tracers .gt. 0) then
         current_state%trsfflux(:,:,:) = 0.0
       endif
       
    end if
    
    if (.not. current_state%halo_column) then
       ! maximum and minimum vertical velocity in each column
       current_state%wmax(target_y_index, target_x_index)=maxval(current_state%w%data(:, current_y_index, current_x_index))
       current_state%wmin(target_y_index, target_x_index)=minval(current_state%w%data(:, current_y_index, current_x_index))
               
       ! work out the column resolved ww, uu, vv
       ww_prime_res(:) = &
            (current_state%w%data(:,current_state%column_local_y,current_state%column_local_x)**2.)
       if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
          uu_prime_res(:) = &
               ((current_state%u%data(:,current_state%column_local_y,current_state%column_local_x) &
               - (current_state%global_grid%configuration%vertical%olubar(:) - current_state%ugal))**2.)
       else
          uu_prime_res(:) = 0.0
       endif
       if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
          vv_prime_res(:) = &
               ((current_state%v%data(:,current_state%column_local_y,current_state%column_local_x) &
               - (current_state%global_grid%configuration%vertical%olvbar(:) - current_state%vgal))**2.)
       else
          vv_prime_res(:) = 0.0
       endif
       ! use column resolved ww, uu, vv to derive total resolved KE for column
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          current_state%reske(target_y_index, target_x_index) = current_state%reske(target_y_index, target_x_index) + &
               (uu_prime_res(k) + uu_prime_res(k+1)                                       &
               + vv_prime_res(k) + vv_prime_res(k+1)                                      &
               + 2_DEFAULT_PRECISION * ww_prime_res(k))*0.25_DEFAULT_PRECISION*           &
               current_state%global_grid%configuration%vertical%dzn(k+1)*                 &
               current_state%global_grid%configuration%vertical%rho(k)
       enddo
       ! divide reske by altitude to make it column mean (as in the LEM)
       current_state%reske(target_y_index, target_x_index) = current_state%reske(target_y_index, target_x_index)/ &
            current_state%global_grid%configuration%vertical%z(current_state%local_grid%size(Z_INDEX))
       ! subke is derive in the subgrid_profile_diagnostics component
      if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then

          ! calculate the lwc maximum and height of lwc max for each column

          !if (current_state%liquid_water_mixing_ratio_index .gt. 0 .and. &
          !     current_state%number_q_fields .ge. current_state%liquid_water_mixing_ratio_index) then
            current_state%qlmax(target_y_index, target_x_index) = &
                  maxval(current_state%qv%data(:,current_y_index, current_x_index))
            !hqlmax(current_y_index, current_x_index) = &
            !     current_state%global_grid%configuration%vertical% &
            !     zn(maxloc(current_state%q(current_state%liquid_water_mixing_ratio_index)%data &
            !     (:,current_y_index, current_x_index)))


            ! calculate the cloud top maximum and minimum for each column
            !
            cloud_content(:) = current_state%ql%data(:,current_y_index, current_x_index)
            !> Include ice if present.
            if (current_state%ice_water_mixing_ratio_index .gt. 0)  &
                cloud_content(:) = cloud_content(:) + current_state%qi%data(:,current_y_index, current_x_index)
            do k = 2, current_state%local_grid%size(Z_INDEX)
              if (cloud_content(k) .gt. qlcrit) then
                current_state%cltop(target_y_index, target_x_index) = &
                      current_state%global_grid%configuration%vertical%zn(k)
              endif
              if (cloud_content(current_state%local_grid%size(Z_INDEX)+1-k) .gt. qlcrit) then
                  current_state%clbas(target_y_index, target_x_index) = &
                      current_state%global_grid%configuration%vertical%zn(current_state%local_grid%size(Z_INDEX)+1-k)
              end if
            enddo ! k loop over height
          !endif
          !
          ! calculate the vapour and liquid water path
          !
          !if (current_state%water_vapour_mixing_ratio_index .gt. 0 .and. &
          !     current_state%number_q_fields .ge. current_state%water_vapour_mixing_ratio_index) then
            do k = 2, current_state%local_grid%size(Z_INDEX)
                current_state%vwp(target_y_index, target_x_index)=current_state%vwp(target_y_index, target_x_index) &
                    +dz_rhon_fac(k)*current_state%qv%data(k, current_y_index, current_x_index)
                current_state%lwp(target_y_index, target_x_index)=current_state%lwp(target_y_index, target_x_index) &
                    +dz_rhon_fac(k)*current_state%ql%data(k, current_y_index, current_x_index)
            enddo
            if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
            !if (current_state%rain_water_mixing_ratio_index > 0) then
                do k = 2, current_state%local_grid%size(Z_INDEX)
                  current_state%rwp(target_y_index, target_x_index)=current_state%rwp(target_y_index, target_x_index) &
                        +dz_rhon_fac(k)*current_state%qr%data(k,current_y_index, current_x_index)
                enddo
            !endif
            !if (current_state%ice_water_mixing_ratio_index > 0) then
                do k = 2, current_state%local_grid%size(Z_INDEX)
                  current_state%iwp(target_y_index, target_x_index)=current_state%iwp(target_y_index, target_x_index) &
                        +dz_rhon_fac(k)*current_state%qi%data(k,current_y_index, current_x_index)
                enddo
                current_state%tot_iwp(target_y_index, target_x_index)=current_state%tot_iwp(target_y_index, target_x_index)+ &
                    current_state%iwp(target_y_index, target_x_index)
            !endif
            !if (current_state%snow_water_mixing_ratio_index > 0) then
                do k = 2, current_state%local_grid%size(Z_INDEX)
                  current_state%swp(target_y_index, target_x_index)=current_state%swp(target_y_index, target_x_index) &
                        +dz_rhon_fac(k)*current_state%qs%data(k,current_y_index, current_x_index)
                enddo
                current_state%tot_iwp(target_y_index, target_x_index)=current_state%tot_iwp(target_y_index, target_x_index)+ &
                    current_state%swp(target_y_index, target_x_index)
            !endif
            !if (current_state%graupel_water_mixing_ratio_index > 0) then
                do k = 2, current_state%local_grid%size(Z_INDEX)
                  current_state%gwp(target_y_index, target_x_index)=current_state%gwp(target_y_index, target_x_index) &
                        +dz_rhon_fac(k)*current_state%qg%data(k,current_y_index, current_x_index)
                enddo
                current_state%tot_iwp(target_y_index, target_x_index)=current_state%tot_iwp(target_y_index, target_x_index)+ &
                    current_state%gwp(target_y_index, target_x_index)
            !endif
          !end if
          end if
      end if

       ! surface flux diagnostics
       if (current_state%use_surface_boundary_conditions) then
          !if (current_state%water_vapour_mixing_ratio_index .gt. 0 .and. &
          !     current_state%number_q_fields .ge. current_state%water_vapour_mixing_ratio_index) then
               current_state%lathf(target_y_index, target_x_index)= &
               (current_state%diff_coefficient%data(1, current_y_index, current_x_index) *   &
                  current_state%global_grid%configuration%vertical%rdzn(2) *  & 
                  (current_state%qv%data(1,current_y_index,current_x_index) - &
                  current_state%qv%data(2,current_y_index,current_x_index))) &
                  * rlvap * current_state%global_grid%configuration%vertical%rhon(1)
          !endif

          if (current_state%th%active) then
              current_state%senhf(target_y_index, target_x_index)= &
              (current_state%diff_coefficient%data(1, current_y_index, current_x_index)  &
                   * current_state%global_grid%configuration%vertical%rdzn(2)  &            
                   * (current_state%th%data(1, current_y_index, current_x_index) &            
                   - current_state%th%data(2, current_y_index, current_x_index) &            
                   - current_state%global_grid%configuration%vertical%dthref(1))) &
                   * current_state%global_grid%configuration%vertical%rhon(1)*cp
          endif
          
          ! Surface Flux of tracers
          if (current_state%n_tracers .gt. 0) then
            if (current_state%scalar_stepping == FORWARD_STEPPING) then

              do n = 1, current_state%n_tracers
                current_state%trsfflux(target_y_index, target_x_index, n) = &
                   (current_state%diff_coefficient%data(1, current_y_index, current_x_index) *   &                   
                    current_state%global_grid%configuration%vertical%rhon(1) * &
                    current_state%global_grid%configuration%vertical%rdzn(2) * & 
                    (current_state%tracer(n)%data(1,current_y_index,current_x_index) - & 
                     current_state%tracer(n)%data(2,current_y_index,current_x_index))) 
              end do
              
            else

              do n = 1, current_state%n_tracers
                current_state%trsfflux(target_y_index, target_x_index, n) = &
                   (current_state%diff_coefficient%data(1, current_y_index, current_x_index) *   &                   
                    current_state%global_grid%configuration%vertical%rhon(1) * &
                    current_state%global_grid%configuration%vertical%rdzn(2) * & 
                    (current_state%ztracer(n)%data(1,current_y_index,current_x_index) - & 
                     current_state%ztracer(n)%data(2,current_y_index,current_x_index))) 
              end do
              
            endif
          endif

       endif      
    end if
    end if
  end subroutine timestep_callback_scalar_diagnostics

end module scalar_diagnostics_mod
