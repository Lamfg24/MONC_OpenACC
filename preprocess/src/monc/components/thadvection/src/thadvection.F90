










!> Specific theta advection, which involves the vertical advection of reference state and advection of mean baroclinicity
module thadvection_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type, &
      component_descriptor_type_v1
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use science_constants_mod, only : G
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_master_log
  use optionsdatabase_mod, only : options_get_real, options_get_integer, options_get_logical
implicit none

  private

  logical :: baroclinicity_use_geostrophic_shear
  logical :: l_advect_mean_baroclinicity
  real(kind=DEFAULT_PRECISION) :: fcoriol, fcoriol_over_G, rate_change_geostrophic_wind_x, rate_change_geostrophic_wind_y, &
       multiplicative_factor_x, multiplicative_factor_y

  ! Local tendency diagnostic variables for this component
  logical :: l_tend_3d_th, l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  logical :: l_tend_pr_tot_th, l_tend_pr_tot_tabs

  public initialisation_callback_thadvection, timestep_callback_thadvection, &
         finalisation_callback_thadvection
contains

  !> Initialisation callback to set up the variables and data needed by the component
  !! @param current_state The current model state
  subroutine initialisation_callback_thadvection(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: iter
    integer :: alloc_z, alloc_y, alloc_x
    real(kind=DEFAULT_PRECISION) :: realnum
    logical :: logicnum

    !AH (08/02/18) - by default l_advect_mean_baroclinicity is false. It
    !                is switch to true if passive_q true and baroclinicity_use_geostrophic_shear is true
    !                If passive q is false code will run but a warning is issued.
    !                Mean advection of baroclinicity only works for dry case
    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "baroclinicity_use_geostrophic_shear") then
        read(current_state%options_database_string(iter,2),*) logicnum
        baroclinicity_use_geostrophic_shear = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "fcoriol") then
        read(current_state%options_database_string(iter,2),*) realnum
        fcoriol = realnum
      else if (current_state%options_database_string(iter,1) .eq. "rate_change_geostrophic_wind_x") then
        read(current_state%options_database_string(iter,2),*) realnum
        rate_change_geostrophic_wind_x = realnum
      else if (current_state%options_database_string(iter,1) .eq. "rate_change_geostrophic_wind_y") then
        read(current_state%options_database_string(iter,2),*) realnum
        rate_change_geostrophic_wind_y = realnum
      end if
    end do

    if (baroclinicity_use_geostrophic_shear) then
       if (current_state%passive_q) then 
          l_advect_mean_baroclinicity = .true.
       else
          call log_master_log(LOG_WARN, "The combination of baroclinicity and active q is not allowed, code will run but"// & 
               " no advection of mean baroclinicity")
          l_advect_mean_baroclinicity = .false.
       endif
    endif


    fcoriol_over_G = fcoriol/G
    multiplicative_factor_x=rate_change_geostrophic_wind_x*current_state%thref0*fcoriol_over_G
    multiplicative_factor_y=rate_change_geostrophic_wind_y*current_state%thref0*fcoriol_over_G

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available
    alloc_z=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    alloc_y=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    alloc_x=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

    l_tend_pr_tot_th  = current_state%th%active 
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_tabs = l_tend_3d_th

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_th) then
      !allocate( current_state%tend_3d_th_thad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_th_thad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_tabs) then
      !allocate( current_state%tend_3d_tabs_thad(current_state%local_grid%size(Z_INDEX),  &
      !                       current_state%local_grid%size(Y_INDEX),  &
      !                       current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_tabs_thad(alloc_z, alloc_y, alloc_x))
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_th) then
      !allocate( current_state%tend_pr_tot_th_thad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_th_thad(alloc_z))
    endif
    if (l_tend_pr_tot_tabs) then
      !allocate( current_state%tend_pr_tot_tabs_thad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_tabs_thad(alloc_z))
    endif

  end subroutine initialisation_callback_thadvection


  subroutine finalisation_callback_thadvection(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(current_state%tend_3d_th_thad)) deallocate(current_state%tend_3d_th_thad)
    if (allocated(current_state%tend_3d_tabs_thad)) deallocate(current_state%tend_3d_tabs_thad)

    if (allocated(current_state%tend_pr_tot_th_thad)) deallocate(current_state%tend_pr_tot_th_thad)
    if (allocated(current_state%tend_pr_tot_tabs_thad)) deallocate(current_state%tend_pr_tot_tabs_thad)

  end subroutine finalisation_callback_thadvection


  !> Timestep callback, will call the two separate procedures to do their advection if needed
  !! @param current_state The current model state
  !! @param target_(x/y)_index This is the index with the halos subtracted. This is needed so that diagnostic does
  !!                           not include halos and to prevent array out-of-bounds
  subroutine timestep_callback_thadvection(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: current_x_index, current_y_index, target_x_index, target_y_index
    logical :: calculate_diagnostics

    !calculate_diagnostics = current_state%diagnostic_sample_timestep &
    !                        .and. .not. current_state%halo_column


    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_th) then
        current_state%tend_pr_tot_th_thad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_tabs) then
        current_state%tend_pr_tot_tabs_thad(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    if ((current_state%modulo_number_3d .eq. 0) &
        .and. (.not. current_state%halo_column)) then
        call save_precomponent_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)
    end if

    call vertical_advection_of_reference_state(current_state, current_state%column_local_y, current_state%column_local_x)
    call advection_of_mean_baroclinicity(current_state, current_state%column_local_y, current_state%column_local_x)

    if ((current_state%modulo_number_3d .eq. 0) &
        .and. (.not. current_state%halo_column)) then
        call compute_component_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback_thadvection

  !> Vertical advection of the reference state.  It doesn't seem consistent to do the advection in this way if
  !! TVD advection of the deviation from the reference state has been selected. Separate vertical advection of the reference
  !! state was introduced to improve energy conservation when carrying out idealized gravity wave simulations in a deep, dry
  !! isothermal layer, for which the difference in potential temp between top and bottom was of order 100K. In less extreme cases
  !! the benefits are unlikely to be significant and with TVD advection energy conservation has been compromised so the best
  !! way forward might be to recombine the reference state into l_th
  !! @param current_state The current model state
  !! @param local_y The local Y of the column
  !! @param local_x The local X of the column
  subroutine vertical_advection_of_reference_state(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: k
    real(kind=DEFAULT_PRECISION) :: sctmp1, sctmp2

    if (current_state%use_anelastic_equations) then
      ! This code only needs to be executed if anelastic, otherwise THREF is constant and the increment calculated here is zero
      do k=2, current_state%local_grid%size(Z_INDEX)
        sctmp1=current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%vertical%dthref(k-1)
        sctmp2=current_state%global_grid%configuration%vertical%tzc2(k)*2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%vertical%dthref(k)
        current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)-(sctmp1*&
             current_state%w%data(k-1, local_y, local_x) + sctmp2*current_state%w%data(k, local_y, local_x))
      end do
    end if
  end subroutine vertical_advection_of_reference_state

  !> Performs advection of the mean baroclinicity if appropriate
  !! @param current_state The current model state
  !! @param local_y The local Y of the column
  !! @param local_x The local X of the column
  subroutine advection_of_mean_baroclinicity(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: k

    if (l_advect_mean_baroclinicity) then
       if (current_state%use_anelastic_equations) then
          do k=2, current_state%local_grid%size(Z_INDEX)
             current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)+&
                  current_state%global_grid%configuration%vertical%thref(k)*fcoriol_over_G*&
                  ((current_state%v%data(k, local_y, local_x) + current_state%vgal) * rate_change_geostrophic_wind_x-&
                  (current_state%u%data(k, local_y, local_x) + current_state%ugal) * rate_change_geostrophic_wind_y)
          end do
       else
          do k=2, current_state%local_grid%size(Z_INDEX)
             current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)+&
                  ((current_state%v%data(k, local_y, local_x) + current_state%vgal) * multiplicative_factor_x-&
                  (current_state%u%data(k, local_y, local_x) + current_state%ugal) * multiplicative_factor_y)                    
          end do
       end if
    end if
    
  end subroutine advection_of_mean_baroclinicity  


   !> Save the 3d tendencies coming into the component.
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine save_precomponent_tendencies(current_state, cxn, cyn)!, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  cxn, cyn!, txn, tyn

    ! Save 3d tendency fields upon request (of 3d or profiles) and availability
    if (l_tend_3d_th) then
      current_state%tend_3d_th_thad(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      current_state%tend_3d_tabs_thad(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn) * &
      current_state%global_grid%configuration%vertical%rprefrcp(:)
    endif

  end subroutine save_precomponent_tendencies


   !> Computation of component tendencies
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine compute_component_tendencies(current_state, cxn, cyn)!, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  cxn, cyn!, txn, tyn

    ! Calculate change in tendency due to component
    if (l_tend_3d_th) then
      current_state%tend_3d_th_thad(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn) - &
      current_state%tend_3d_th_thad(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      current_state%tend_3d_tabs_thad(:,cyn,cxn)=   &
       current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)   &
        - current_state%tend_3d_tabs_thad(:,cyn,cxn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_th) then
      current_state%tend_pr_tot_th_thad(:)=current_state%tend_pr_tot_th_thad(:) + &
      current_state%tend_3d_th_thad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_tabs) then
      current_state%tend_pr_tot_tabs_thad(:)=current_state%tend_pr_tot_tabs_thad(:) + &
      current_state%tend_3d_tabs_thad(:,cyn,cxn)
    endif

  end subroutine compute_component_tendencies

end module thadvection_mod
