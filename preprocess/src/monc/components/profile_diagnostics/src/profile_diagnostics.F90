










module profile_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type, &
       component_descriptor_type_v1
  use grids_mod, only : vertical_grid_configuration_type, Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real, options_get_string, options_get_logical
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use q_indices_mod, only: get_q_index, standard_q_names
  use saturation_mod, only: qsaturation
  use logging_mod, only : LOG_ERROR, log_master_log  
  use def_tvd_diagnostic_terms, only: tvd_dgs_terms, allocate_tvd_diagnostic_terms, &
                                      deallocate_tvd_diagnostic_terms
  use conversions_mod, only : conv_to_uppercase
  use def_merge_atm, only: str_merge_atm

  implicit none

  private

  integer :: total_points, iqv=0, iql=0, iqr=0, iqi=0, iqs=0,    &
                           iqg=0
  integer ::  inl=0, inr=0, ini=0, ins=0,    &
              ing=0

  ! 3D and local profile total binary cloud masks
  character(len=STRING_LENGTH)                                :: cloud_mask_method
  logical :: l_partial_liq_ice

  real(kind=DEFAULT_PRECISION) :: qlcrit, qicrit, max_height_cloud
  ! character string to determine the advection scheme
  character(len=STRING_LENGTH) :: advect_theta, advect_q, advect_flow

  type(vertical_grid_configuration_type) :: vertical_grid

  public tvd_dgs_terms, initialisation_callback_profile_diagnostics, &
         timestep_callback_profile_diagnostics, finalisation_callback_profile_diagnostics

contains

  subroutine initialisation_callback_profile_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, iter
    real(kind=DEFAULT_PRECISION) :: realnum
    logical :: logicnum, l_cloud_mask

    vertical_grid=current_state%global_grid%configuration%vertical

    total_points=(current_state%local_grid%size(Y_INDEX) * current_state%local_grid%size(X_INDEX))
    
    ! allocate local arrays for the horizontal wind averages
    allocate(current_state%u_wind_tot(current_state%local_grid%size(Z_INDEX)) &
         , current_state%v_wind_tot(current_state%local_grid%size(Z_INDEX))   &
         , current_state%ww_tot(current_state%local_grid%size(Z_INDEX))       &
         , current_state%w_wind_tot(current_state%local_grid%size(Z_INDEX))   &
         , current_state%prefn(current_state%local_grid%size(Z_INDEX))        &
         , current_state%rho(current_state%local_grid%size(Z_INDEX))          &
         , current_state%rhon(current_state%local_grid%size(Z_INDEX))         &
         , current_state%uinit(current_state%local_grid%size(Z_INDEX))       &
         , current_state%vinit(current_state%local_grid%size(Z_INDEX)))
    allocate( current_state%uw_tot(current_state%local_grid%size(Z_INDEX))     &
         ,   current_state%vw_tot(current_state%local_grid%size(Z_INDEX))     &
         ,   current_state%uv_tot(current_state%local_grid%size(Z_INDEX))     &
         ,   current_state%wwww_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   current_state%www_tot(current_state%local_grid%size(Z_INDEX)) )
    
    if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
       allocate(current_state%uprime_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%uprime(current_state%local_grid%size(Z_INDEX)))
    end if
    
    if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
       allocate(current_state%vprime_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%vprime(current_state%local_grid%size(Z_INDEX)))
    end if

    if (allocated(current_state%vprime) .and. allocated(current_state%uprime)) then
       allocate(current_state%wke_tot(current_state%local_grid%size(Z_INDEX)))
    endif

    if (current_state%th%active) then
       allocate(current_state%theta_tot(current_state%local_grid%size(Z_INDEX)) &
            , current_state%thref(current_state%local_grid%size(Z_INDEX))       &
            , current_state%thinit(current_state%local_grid%size(Z_INDEX))      &
            , current_state%wtheta_ad_tot(current_state%local_grid%size(Z_INDEX)) &
            , current_state%wtheta_cn_tot(current_state%local_grid%size(Z_INDEX)) &
            , current_state%th2_tot(current_state%local_grid%size(Z_INDEX)) )
    endif
    
    ! determine advection scheme used for mom and scalars
    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "advection_flow_fields") then
        advect_flow = current_state%options_database_string(iter,2)
      else if (current_state%options_database_string(iter,1) .eq. "advection_theta_field") then
        advect_theta = current_state%options_database_string(iter,2)
      else if (current_state%options_database_string(iter,1) .eq. "advection_q_fields") then
        advect_q = current_state%options_database_string(iter,2)
      end if
    end do

    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then 
       iqv= 1 !get_q_index(standard_q_names%VAPOUR, 'profile_diags')
       iql= 2 !get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'profile_diags')

       do iter = 1,current_state%config_args
         if (current_state%options_database_string(iter,1) .eq. "qlcrit") then
           read(current_state%options_database_string(iter,2),*) realnum
           qlcrit = realnum
         else if (current_state%options_database_string(iter,1) .eq. "qicrit") then
           read(current_state%options_database_string(iter,2),*) realnum
           qicrit = realnum
         else if (current_state%options_database_string(iter,1) .eq. "max_height_cloud") then
           read(current_state%options_database_string(iter,2),*) realnum
           max_height_cloud = realnum
         else if (current_state%options_database_string(iter,1) .eq. "l_partial_liq_ice") then
           read(current_state%options_database_string(iter,2),*) logicnum
           l_partial_liq_ice = logicnum
         else if (current_state%options_database_string(iter,1) .eq. "l_cloud_mask") then
           read(current_state%options_database_string(iter,2),*) logicnum
           l_cloud_mask = logicnum
         else if (current_state%options_database_string(iter,1) .eq. "cloud_mask_method") then
           cloud_mask_method = conv_to_uppercase(current_state%options_database_string(iter,2))
         end if
       end do


       allocate(current_state%qv_tot(current_state%local_grid%size(Z_INDEX))  &
            , current_state%ql_tot(current_state%local_grid%size(Z_INDEX)),   &
            current_state%q_temp(current_state%local_grid%size(Z_INDEX)),   &
            current_state%wqv_cn_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%wqv_ad_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%wql_cn_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%wql_ad_tot(current_state%local_grid%size(Z_INDEX)))
       if (current_state%th%active) &
            allocate(current_state%rh_tot(current_state%local_grid%size(Z_INDEX)))

       ! allocate other hydrometeors. Allocation dependent on index being set in
       ! appropriate microphysics scheme (see casim component from example)
       if (current_state%rain_water_mixing_ratio_index > 0) then
          iqr = current_state%rain_water_mixing_ratio_index
          allocate(current_state%qr_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(current_state%wqr_cn_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(current_state%wqr_ad_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       if (current_state%ice_water_mixing_ratio_index > 0) then
          iqi = current_state%ice_water_mixing_ratio_index 
          allocate(current_state%qi_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(current_state%wqi_cn_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(current_state%wqi_ad_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       if (current_state%snow_water_mixing_ratio_index > 0) then
          iqs = current_state%snow_water_mixing_ratio_index
          allocate(current_state%qs_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(current_state%wqs_cn_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(current_state%wqs_ad_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       if (current_state%graupel_water_mixing_ratio_index > 0) then
          iqg = current_state%graupel_water_mixing_ratio_index
          allocate(current_state%qg_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(current_state%wqg_cn_tot(current_state%local_grid%size(Z_INDEX)))
          allocate(current_state%wqg_ad_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       ! allocate number concentrations for the hydrometeors
       if (current_state%liquid_water_nc_index > 0) then
          inl = current_state%liquid_water_nc_index
          allocate(current_state%nl_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       if (current_state%rain_water_nc_index > 0) then
          inr = current_state%rain_water_nc_index
          allocate(current_state%nr_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       if (current_state%ice_water_nc_index > 0) then
          ini = current_state%ice_water_nc_index 
          allocate(current_state%ni_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       if (current_state%snow_water_nc_index > 0) then
          ins = current_state%snow_water_nc_index
          allocate(current_state%ns_tot(current_state%local_grid%size(Z_INDEX)))
       endif
       if (current_state%graupel_water_nc_index > 0) then
          ing = current_state%graupel_water_nc_index
          allocate(current_state%ng_tot(current_state%local_grid%size(Z_INDEX)))
       endif

       ! arrange and allocate cloud fraction diagnostics...3d mask is optional
       if (.not. (cloud_mask_method == "DEFAULT"  .or. &
                  cloud_mask_method == "SOCRATES" .or. &
                  cloud_mask_method == "RCEMIP")) then
         call log_master_log(LOG_ERROR,                                            &
          "Requested cloud_mask_method is invalid.  Check profile_diagnostics.F90") 
       end if ! cloud_mask_method validity check
       if (l_cloud_mask) &
       allocate(current_state%cloud_mask(current_state%local_grid%size(Z_INDEX), &
                             current_state%local_grid%size(Y_INDEX),               &
                             current_state%local_grid%size(X_INDEX)))
       allocate(current_state%cloud_mask_tot(current_state%local_grid%size(Z_INDEX)))
       allocate(current_state%cloud_liq_mask_tot(current_state%local_grid%size(Z_INDEX)))
       allocate(current_state%cloud_ice_mask_tot(current_state%local_grid%size(Z_INDEX)))
    endif

    call allocate_tvd_diagnostic_terms(current_state, tvd_dgs_terms)

  end subroutine initialisation_callback_profile_diagnostics

  subroutine timestep_callback_profile_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i, iq_tmp, icol, jcol
    real(kind=DEFAULT_PRECISION) :: cltop_col, clbas_col, qv, qc, TdegK, Pmb &
         , q_s, exner
    real(kind=DEFAULT_PRECISION) :: uprime_w_local, vprime_w_local &
         , thprime_w_local, qprime_w_local 

    icol=current_state%column_local_x
    jcol=current_state%column_local_y

    if (current_state%modulo_number_1d .ne. 0) return

    if (current_state%first_timestep_column) then
       current_state%u_wind_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(current_state%uprime_tot)) current_state%uprime_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(current_state%uprime)) current_state%uprime(:) = 0.0_DEFAULT_PRECISION
       current_state%v_wind_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(current_state%vprime_tot)) current_state%vprime_tot(:) = 0.0_DEFAULT_PRECISION
       if (allocated(current_state%vprime)) current_state%vprime(:) = 0.0_DEFAULT_PRECISION
       if (allocated(current_state%wke_tot)) current_state%wke_tot(:) =  0.0_DEFAULT_PRECISION
       
       current_state%w_wind_tot(:) = 0.0_DEFAULT_PRECISION
       current_state%ww_tot(:) = 0.0_DEFAULT_PRECISION
       current_state%wwww_tot(:) = 0.0_DEFAULT_PRECISION
       current_state%www_tot(:) = 0.0_DEFAULT_PRECISION
       current_state%uw_tot(:)     = 0.0_DEFAULT_PRECISION
       current_state%vw_tot(:)     = 0.0_DEFAULT_PRECISION
       current_state%uv_tot(:)     = 0.0_DEFAULT_PRECISION
       
       if (current_state%th%active) then
          current_state%wtheta_ad_tot(:) = 0.0_DEFAULT_PRECISION
          current_state%wtheta_cn_tot(:) = 0.0_DEFAULT_PRECISION
          current_state%th2_tot(:)    = 0.0_DEFAULT_PRECISION
          current_state%theta_tot(:)=0.0_DEFAULT_PRECISION
       endif
       if (.not. current_state%passive_q .and. &
            current_state%number_q_fields .gt. 0) then 
          current_state%q_temp(:)=0.0_DEFAULT_PRECISION
          current_state%qv_tot(:)=0.0_DEFAULT_PRECISION
          current_state%ql_tot(:)=0.0_DEFAULT_PRECISION
          current_state%wqv_cn_tot(:)=0.0_DEFAULT_PRECISION
          current_state%wqv_ad_tot(:)=0.0_DEFAULT_PRECISION
          current_state%wql_cn_tot(:)=0.0_DEFAULT_PRECISION
          current_state%wql_ad_tot(:)=0.0_DEFAULT_PRECISION
          if (current_state%th%active) &
               current_state%rh_tot(:) = 0.0_DEFAULT_PRECISION
          if (iqr > 0) then 
             current_state%qr_tot(:) = 0.0_DEFAULT_PRECISION
             current_state%wqr_cn_tot(:) = 0.0_DEFAULT_PRECISION
             current_state%wqr_ad_tot(:) = 0.0_DEFAULT_PRECISION
          endif
          if (iqi > 0) then
             current_state%qi_tot(:) = 0.0_DEFAULT_PRECISION
             current_state%wqi_cn_tot(:) = 0.0_DEFAULT_PRECISION
             current_state%wqi_ad_tot(:) = 0.0_DEFAULT_PRECISION
          endif
          if (iqs > 0) then
             current_state%qs_tot(:) = 0.0_DEFAULT_PRECISION
             current_state%wqs_cn_tot(:) = 0.0_DEFAULT_PRECISION
             current_state%wqs_ad_tot(:) = 0.0_DEFAULT_PRECISION
          endif
          if (iqg > 0) then 
             current_state%qg_tot(:) = 0.0_DEFAULT_PRECISION
             current_state%wqg_cn_tot(:) = 0.0_DEFAULT_PRECISION
             current_state%wqg_ad_tot(:) = 0.0_DEFAULT_PRECISION
          endif

          if (inl > 0) current_state%nl_tot(:) = 0.0_DEFAULT_PRECISION
          if (inr > 0) current_state%nr_tot(:) = 0.0_DEFAULT_PRECISION
          if (ini > 0) current_state%ni_tot(:) = 0.0_DEFAULT_PRECISION
          if (ins > 0) current_state%ns_tot(:) = 0.0_DEFAULT_PRECISION
          if (ing > 0) current_state%ng_tot(:) = 0.0_DEFAULT_PRECISION
    
          if (allocated(current_state%cloud_mask)) current_state%cloud_mask(:,:,:) = 0.0_DEFAULT_PRECISION
          current_state%cloud_mask_tot(:) = 0.0_DEFAULT_PRECISION
          current_state%cloud_liq_mask_tot(:) = 0.0_DEFAULT_PRECISION
          current_state%cloud_ice_mask_tot(:) = 0.0_DEFAULT_PRECISION
       endif
    end if
    !

    if (.not. current_state%halo_column) then
    ! work out the sum of u and v wind over local domain
    do k=1, current_state%local_grid%size(Z_INDEX)
       current_state%u_wind_tot(k) = current_state%u_wind_tot(k) + &
            (current_state%u%data(k,jcol,icol)  &
            + current_state%ugal)
       if (allocated(current_state%uprime_tot)) then
          current_state%uprime(k) =  &
               (current_state%u%data(k, jcol, icol) &
               - (current_state%global_grid%configuration%vertical%olubar(k) - current_state%ugal))
          current_state%uprime_tot(k) = current_state%uprime_tot(k) + current_state%uprime(k)**2.0_DEFAULT_precision
       end if
       current_state%v_wind_tot(k) = current_state%v_wind_tot(k) + &
            (current_state%v%data(k,jcol,icol)  &
            + current_state%vgal)
       if (allocated(current_state%vprime_tot)) then
          current_state%vprime(k) = &
          (current_state%v%data(k, jcol, icol) &
               - (current_state%global_grid%configuration%vertical%olvbar(k) - current_state%vgal))
          current_state%vprime_tot(k) = current_state%vprime_tot(k) + current_state%vprime(k)**2.0_DEFAULT_precision
       end if
    enddo
    ! Note: If TVD advection used, a conversion for the grid (from the LEM (RESDGS.653,654))
    !       is used but it is unclear whether that is correct. Also, it is worth noting that 
    !       www and wwww have no offset relating to the advection scheme - AH 21/09/2017
    if (trim(advect_flow) .eq. "pw") then 
       current_state%ww_tot(:) = current_state%ww_tot(:) + &
            (current_state%w%data(:,jcol,icol)**2.)
    else if (trim(advect_flow) .eq. "tvd") then
       do k=2, current_state%local_grid%size(Z_INDEX)
          current_state%ww_tot(k) = current_state%ww_tot(k) + 0.5 * &
               ((current_state%w%data(k-1,jcol,icol) + current_state%w%data(k,jcol,icol)) * &
               tvd_dgs_terms%adv_w_dgs(k,jcol,icol))
       enddo
    endif
    
    do k=2, current_state%local_grid%size(Z_INDEX)       
       current_state%www_tot(k) = current_state%www_tot(k) + &
            (current_state%w%data(k,jcol,icol)**3.)
       current_state%wwww_tot(k) = current_state%wwww_tot(k) + &
            (current_state%w%data(k,jcol,icol)**4.)
       current_state%w_wind_tot(k) = current_state%w_wind_tot(k) + &
            (current_state%w%data(k,jcol,icol))
    enddo

!      <u'w'> and <v'w'> are on w-points, so we interpolate u and v both horizontally and vertically.
!      NOTE: all "prime_w" values are the prognostic interpolated on to the w levels
    if (trim(advect_flow) .eq. "pw") then 
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          uprime_w_local =  &
               0.25 * ( current_state%u%data(k,jcol,icol)   + &
               current_state%u%data(k,jcol,icol-1) + &
               current_state%u%data(k+1,jcol,icol) + &
               current_state%u%data(k+1,jcol,icol-1) ) + &
               current_state%ugal
          if (allocated(current_state%global_grid%configuration%vertical%olubar)) &
               uprime_w_local = uprime_w_local - &
               0.5  * ( current_state%global_grid%configuration%vertical%olubar(k) + &
                  current_state%global_grid%configuration%vertical%olubar(k+1) )
          current_state%uw_tot(k) = current_state%uw_tot(k) + uprime_w_local * &
               current_state%w%data(k,jcol,icol)
       enddo
    else if (trim(advect_flow) .eq. "tvd") then
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          current_state%uw_tot(k) = current_state%uw_tot(k) + 0.5 * (      &
               (current_state%w%data(k,jcol,icol) + current_state%w%data(k,jcol,icol+1)) * &
               tvd_dgs_terms%adv_u_dgs(k+1,jcol,icol))
       enddo
    endif
    
    if (trim(advect_flow) .eq. "pw") then 
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          vprime_w_local = &
               0.25 * ( current_state%v%data(k,jcol,icol)   + &
               current_state%v%data(k,jcol-1,icol) + &
               current_state%v%data(k+1,jcol,icol) + &
               current_state%v%data(k+1,jcol-1,icol) ) + &
               current_state%vgal
          if (allocated(current_state%global_grid%configuration%vertical%olvbar)) &
                  vprime_w_local = vprime_w_local - &
                  0.5  * ( current_state%global_grid%configuration%vertical%olvbar(k) + &
                  current_state%global_grid%configuration%vertical%olvbar(k+1) )             
          current_state%vw_tot(k) = current_state%vw_tot(k) + vprime_w_local * &
               current_state%w%data(k,jcol,icol)
       enddo
    else if (trim(advect_flow) .eq. "tvd") then
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          current_state%vw_tot(k) = current_state%vw_tot(k) + 0.5 * (      &
               (current_state%w%data(k,jcol,icol) + current_state%w%data(k,jcol+1,icol)) * &
                  tvd_dgs_terms%adv_v_dgs(k+1,jcol,icol))
       enddo
    endif
    
    if (allocated(current_state%global_grid%configuration%vertical%olvbar) .and.               &
         allocated(current_state%global_grid%configuration%vertical%olubar)) then
       ! LEM equivalent code RESDGS.784, 785, comment from LEM
       ! No attempt has been made to provide a calculation of this term which
       ! is consistent with TVD advection 
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          current_state%uv_tot(k) = current_state%uv_tot(k) + 0.25*current_state%uprime(k)*(current_state%vprime(k) &
               + (current_state%v%data(k,jcol,icol+1))   &
               ! v primed term (not squared) at y - 1 , do whole calc for vprime at y-1
               + (current_state%v%data(k,jcol-1,icol)   &
               - (current_state%global_grid%configuration%vertical%olvbar(k) - current_state%vgal))   &
               + (current_state%v%data(k,jcol-1,icol+1)))
       enddo
    endif
       
    if (allocated(current_state%wke_tot)) then
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          uprime_w_local =  &
               0.25_DEFAULT_PRECISION * ( current_state%u%data(k,jcol,icol)   + &
               current_state%u%data(k,jcol,icol-1) + &
               current_state%u%data(k+1,jcol,icol) + &
               current_state%u%data(k+1,jcol,icol-1) ) + &
               current_state%ugal
          if (allocated(current_state%global_grid%configuration%vertical%olubar)) &
               uprime_w_local = uprime_w_local - &
               0.5_DEFAULT_PRECISION * ( current_state%global_grid%configuration%vertical%olubar(k) + &
                  current_state%global_grid%configuration%vertical%olubar(k+1) )
          vprime_w_local = &
               0.25_DEFAULT_PRECISION * ( current_state%v%data(k,jcol,icol)   + &
               current_state%v%data(k,jcol-1,icol) + &
               current_state%v%data(k+1,jcol,icol) + &
               current_state%v%data(k+1,jcol-1,icol) ) + &
               current_state%vgal
          if (allocated(current_state%global_grid%configuration%vertical%olvbar)) &
               vprime_w_local = vprime_w_local - &
               0.5_DEFAULT_PRECISION * ( current_state%global_grid%configuration%vertical%olvbar(k) + &
               current_state%global_grid%configuration%vertical%olvbar(k+1) )

          current_state%wke_tot(k) = current_state%wke_tot(k) + 0.5_DEFAULT_PRECISION *                     &
               current_state%global_grid%configuration%vertical%rhon(k) *       &
               current_state%w%data(k,jcol,icol) *                              &
               ( uprime_w_local * uprime_w_local +                              &
                 vprime_w_local * vprime_w_local +                              &
                 current_state%w%data(k,jcol,icol) * current_state%w%data(k,jcol,icol) )
       enddo
    endif
       
    if (current_state%th%active) then
       do k=1, current_state%local_grid%size(Z_INDEX)
          current_state%theta_tot(k) = current_state%theta_tot(k) + &
               (current_state%th%data(k,jcol,icol) &
               + current_state%global_grid%configuration%vertical%thref(k))
          current_state%th2_tot(k) = current_state%th2_tot(k) + &
               (current_state%th%data(k,jcol,icol) - &
               current_state%global_grid%configuration%vertical%olthbar(k) )**2
       enddo
       !       <w'theta'> is on w-levels, so theta is interpolated to w-levels.
       !       Set wtheta_ad_tot to the actual model TH flux.
!       Set wtheta_cn_tot slot to flux seen by w times the w-momentum equation. 
       !       This is the flux relevant for the tke balance. 
       !       The two are equal if the centred scheme is used but differ if TVD advection is used on scalars). 
       if (trim(advect_theta) .eq. "pw") then 
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             thprime_w_local = 0.5 * (                       &
                  current_state%th%data(k,jcol,icol) +       &
                  current_state%th%data(k+1,jcol,icol) )
             current_state%wtheta_ad_tot(k) = current_state%wtheta_ad_tot(k) +           &
                  (current_state%w%data(k,jcol,icol) *        &
                  thprime_w_local)
             current_state%wtheta_cn_tot(k) = current_state%wtheta_cn_tot(k)
          enddo
       else if (trim(advect_theta) .eq. "tvd") then
          do k=2, current_state%local_grid%size(Z_INDEX)-1
             thprime_w_local = 0.5 * (                        &
                  current_state%th%data(k,jcol,icol) +       &
                  current_state%th%data(k+1,jcol,icol) )
             current_state%wtheta_cn_tot(k) = current_state%wtheta_cn_tot(k) + &
                  current_state%w%data(k,jcol,icol) * &
                  thprime_w_local
             current_state%wtheta_ad_tot(k) = current_state%wtheta_ad_tot(k) +           &
                  current_state%w%data(k,jcol,icol) *         &
                  tvd_dgs_terms%adv_th_dgs(k+1,jcol,icol)
          enddo
       endif
    endif
    
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
       current_state%qv_tot(:) = current_state%qv_tot(:) + (current_state%qv%data(:,jcol,icol))
       current_state%ql_tot(:) = current_state%ql_tot(:) + (current_state%ql%data(:,jcol,icol))
       !if (current_state%th%active) then
       !   ! temporary code for RH calculation
       !   do k=1, current_state%local_grid%size(Z_INDEX)
       !      exner = current_state%global_grid%configuration%vertical%rprefrcp(k)
       !      Pmb   = (current_state%global_grid%configuration%vertical%prefn(k)/100.)
       !      qv    = current_state%qv%data(k, jcol,icol)
       !      TdegK = (current_state%th%data(k,jcol,icol) &
       !           + current_state%global_grid%configuration%vertical%thref(k))*exner
       !      q_s = qsaturation(TdegK, Pmb)
       !      current_state%rh_tot(k) = current_state%rh_tot(k) + (qv/q_s)
       !   enddo
       !endif
       ! hydrometeor mass profiles
       if (iqr > 0) &
            current_state%qr_tot(:) = current_state%qr_tot(:) + (current_state%qr%data(:,jcol,icol))
       if (iqi > 0) &
            current_state%qi_tot(:) = current_state%qi_tot(:) + (current_state%qi%data(:,jcol,icol))
       if (iqs > 0) &
            current_state%qs_tot(:) = current_state%qs_tot(:) + (current_state%qs%data(:,jcol,icol))
       if (iqg > 0) &
            current_state%qg_tot(:) = current_state%qg_tot(:) + (current_state%qg%data(:,jcol,icol))
       ! hydrometeor number concentration
       if (inl > 0) &
            current_state%nl_tot(:) = current_state%nl_tot(:) + (current_state%nl%data(:,jcol,icol))
       if (inr > 0) &
            current_state%nr_tot(:) = current_state%nr_tot(:) + (current_state%nr%data(:,jcol,icol))
       if (ini > 0) &
            current_state%ni_tot(:) = current_state%ni_tot(:) + (current_state%ni%data(:,jcol,icol))
       if (ins > 0) &
            current_state%ns_tot(:) = current_state%ns_tot(:) + (current_state%ns%data(:,jcol,icol))
       if (ing > 0) &
            current_state%ng_tot(:) = current_state%ng_tot(:) + (current_state%ng%data(:,jcol,icol))
       !
       ! moisture field fluxes
       ! vapour
       call calculate_wq(current_state, jcol, icol, current_state%qv%data, current_state%wqv_cn_tot, current_state%wqv_ad_tot, &
                        advect_q, tvd_dgs_terms%adv_qv_dgs)
       ! cloud liquid
       call calculate_wq(current_state, jcol, icol, current_state%ql%data, current_state%wql_cn_tot, &
                         current_state%wql_ad_tot, advect_q, tvd_dgs_terms%adv_ql_dgs)
       ! rain mass
       if (iqr > 0) &
            call calculate_wq(current_state, jcol, icol, current_state%qr%data, current_state%wqr_cn_tot, &
                              current_state%wqr_ad_tot, advect_q, tvd_dgs_terms%adv_qr_dgs)
       ! ice mass
       if (iqi > 0) &
            call calculate_wq(current_state, jcol, icol, current_state%qi%data, current_state%wqi_cn_tot, &
                              current_state%wqi_ad_tot, advect_q, tvd_dgs_terms%adv_qi_dgs)
       ! snow mass
       if (iqs > 0) & 
            call calculate_wq(current_state, jcol, icol, current_state%qs%data, current_state%wqs_cn_tot, &
                              current_state%wqs_ad_tot, advect_q, tvd_dgs_terms%adv_qs_dgs)
          ! graupel mass
       if (iqg > 0) &
            call calculate_wq(current_state, jcol, icol, current_state%qg%data, &
                              current_state%wqg_cn_tot, current_state%wqg_ad_tot, advect_q, tvd_dgs_terms%adv_qg_dgs)

       ! cloud mask / cloud fraction
       call calculate_cloud_mask(current_state, jcol, icol)
    endif ! end q_passive and number q field test
   endif
  end subroutine timestep_callback_profile_diagnostics

!  !> Frees up the memory associated with the advection
!  !! @param current_state The current model state
  subroutine finalisation_callback_profile_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call deallocate_tvd_diagnostic_terms(current_state, tvd_dgs_terms)

  end subroutine finalisation_callback_profile_diagnostics

  subroutine calculate_wq(current_state, jcol, icol, q, wq_cn, wq_ad, advect_q, adv_q_dgs)
    
    type(model_state_type), target, intent(inout) :: current_state
    
    real(kind=DEFAULT_PRECISION), intent(inout) :: q(:,:,:), wq_cn(:), wq_ad(:), adv_q_dgs(:,:,:)
    character(len=*), intent(in) :: advect_q
    integer, intent(in) :: jcol, icol!, iq
    
    integer :: k
  
    if (trim(advect_q) .eq. "pw") then 
       do k=1, current_state%local_grid%size(Z_INDEX)-1
          wq_ad(k) = wq_ad(k) + ( 0.5 * ( & 
               q(k,jcol,icol) + &
               q(k+1,jcol,icol)) * &
               current_state%w%data(k,jcol,icol))
          wq_cn(k) = wq_ad(k)
       enddo
    else if (trim(advect_q) .eq. "tvd") then
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          wq_ad(k) = wq_ad(k) + (current_state%w%data(k,jcol,icol) & 
               *adv_q_dgs(k+1, jcol, icol) )
               !* tvd_dgs_terms%adv_q_dgs(k+1, jcol, icol, iq) )
          wq_cn(k) = wq_cn(k) +  ( 0.5 * (  & 
               q(k,jcol,icol) + &
               q(k+1,jcol,icol)) * &
               current_state%w%data(k,jcol,icol) )
       enddo
    endif

  end subroutine calculate_wq
    
  !---------------------------------------------------------------------
  !> This routine calculates:
  !    cloud_mask              a binary 3D total cloud mask
  !                            optional: l_cloud_mask
  !    cloud_mask_tot          the total cloud mask profile local sum
  !    cloud_liq_mask_tot      the liquid cloud mask profile local sum
  !    cloud_ice_mask_tot      the ice cloud mask profile local sum
  !  The latter 3 are intended to be transformed into
  !  cloud fraction profile diagnostics via xml processing.
  !
  !
  !> SOCRATES method
  !  The definition used here is consistent with that within 
  !  SOCRATES, where i_cloud_representation = 2, that is, liquid and
  !  ice cloud are treated separately, with their individual fractions
  !  within a cell summing to 1 based their mixing ratios relative to
  !  the cell total.  
  !  Each element of the SOCRATES calculation is reproduced
  !  here because the current formulation within SOCRATES will
  !  not be available here when that component is not enabled.
  !    See: def_merge_atm.F90
  !         merge_atm_data.F90
  !  A better solution might be to have both work from 
  !  the same source for easy consistency in the case where
  !  definitions of cloudy cells were to change.
  !  We do not apply special consideration to the cases where
  !  MONC is run with SOCRATES and i_cloud_representation != 2.
  !  That is, values are calculated here even when cloud is 
  !  off in SOCRATES, and no stratiform/convective distinction
  !  is made. 
  !
  !
  !> DEFAULT method
  !  Cloud fraction is based on species exceeding qlcrit and qicrit.
  !  This definition is NOT consistent with that used in the conditional
  !  diagnostics routine (condition "AC") because: 
  !         Ice cloud considers the combination of ice and snow species.
  !  It does not include consideration of rain and graupel fields.  
  !  Liquid and ice cloud are treated separately, with their individual
  !  fractions within a cell summing to 1 based their mixing ratios 
  !  relative to the cell total.
  !  Ice cloud considered as the combination of ice and snow species.
  subroutine calculate_cloud_mask(current_state, jcol, icol)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: jcol, icol
    integer :: k, target_y_index, target_x_index
    logical :: l_prepare_3d_mask, cloud_present

    ! Create a local copy of the precipitation multiplication factors
    ! for the SOCRATES cloud fraction calculation from the str_merge_atm
    ! type structure.
    ! The factors below were derived as part of J. Petch PhD
    ! These are used in the SOCRATES method.
    type (str_merge_atm) :: merge_fields

    ! Local temporary terms
    real(kind=DEFAULT_PRECISION) :: templ, tempi, tempt, qsat_thresh

    target_y_index = jcol - current_state%local_grid%halo_size(Y_INDEX)
    target_x_index = icol - current_state%local_grid%halo_size(X_INDEX)

    templ = 0.0_DEFAULT_PRECISION
    tempi = 0.0_DEFAULT_PRECISION    
    tempt = 0.0_DEFAULT_PRECISION

    l_prepare_3d_mask = allocated(current_state%cloud_mask)

    do k=2, current_state%local_grid%size(Z_INDEX)

      if (current_state%global_grid%configuration%vertical%zn(k) > max_height_cloud) exit

      !> Collect available condensate amounts
      if (iql > 0) &
        templ = current_state%ql%data(k, jcol, icol)
      if (iqi > 0) &
        tempi = current_state%qi%data(k, jcol, icol)

      !> Check cloud_mask_method and modify as needed
      !> The SOCRATES method considers rain, snow, and graupel.
      if (cloud_mask_method == "SOCRATES") then
        if (iqr > 0) &
          templ = templ + merge_fields%rainfac  * current_state%qr%data(k, jcol, icol)
        if (iqs > 0) &
          tempi = tempi + merge_fields%snowfac  * current_state%qs%data(k, jcol, icol)
        if (iqg > 0) &
          tempi = tempi + merge_fields%graupfac * current_state%qg%data(k, jcol, icol)
      end if ! check cloud_mask_method

      !> Consider the snow species as cloud for DEFAULT and RCEMIP
      if ((cloud_mask_method == "DEFAULT" .or. &
           cloud_mask_method == "RCEMIP") .and. iqs > 0) then
        tempi = tempi + current_state%qs%data(k, jcol, icol)
      end if ! check cloud_mask_method

      !> The RCEMIP method
      !> 1e-5 g/g, or 1 % of the saturation mixing ratio over water, whichever is smaller
      if (cloud_mask_method == "RCEMIP") then
        qsat_thresh = min(1e-5_DEFAULT_PRECISION,&
                          qsaturation( (current_state%th%data(k,jcol,icol) &
                            + current_state%global_grid%configuration%vertical%thref(k))*&
                            current_state%global_grid%configuration%vertical%rprefrcp(k), &
                          current_state%global_grid%configuration%vertical%prefn(k)/100.) &
                          * 0.01_DEFAULT_PRECISION)
      end if

      !> Work out cloud fractions
      tempt = templ + tempi
      if (cloud_mask_method == "SOCRATES") then
        cloud_present = (tempt > EPSILON(tempt))
      !else if (cloud_mask_method == "RCEMIP") then
      !  cloud_present = (tempt > qsat_thresh)
      else ! DEFAULT
        cloud_present = (templ > qlcrit .or. tempi > qicrit .or. (templ+tempi) > qlcrit)
      end if 

      if (cloud_present) then
        if (l_prepare_3d_mask)    &
           current_state%cloud_mask(k, target_y_index, target_x_index) = 1.0_DEFAULT_PRECISION
        current_state%cloud_mask_tot(k) = current_state%cloud_mask_tot(k) + 1.0_DEFAULT_PRECISION

        if (l_partial_liq_ice) then ! separated cloud, partial coverage
          current_state%cloud_liq_mask_tot(k) = current_state%cloud_liq_mask_tot(k) + (templ / tempt)
          current_state%cloud_ice_mask_tot(k) = current_state%cloud_ice_mask_tot(k) + (tempi / tempt)
        else ! homogeneous cloud, full coverage for both types of condensate
          if (templ > qlcrit)         &
             current_state%cloud_liq_mask_tot(k) = current_state%cloud_liq_mask_tot(k) + 1.0_DEFAULT_PRECISION
          if (tempi > qicrit)         &
             current_state%cloud_ice_mask_tot(k) = current_state%cloud_ice_mask_tot(k) + 1.0_DEFAULT_PRECISION
        end if ! check liquid/ice partition method
      end if ! check cloud_present

    end do ! k loop over vertical model levels

  end subroutine calculate_cloud_mask


end module profile_diagnostics_mod
