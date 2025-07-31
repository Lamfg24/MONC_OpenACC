!> This calculates the coriolis and mean pressure gradient terms which impact su and sv fields
module coriolis_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type, component_descriptor_type_v1
  use optionsdatabase_mod, only : options_get_logical, options_get_real, options_get_integer
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: baroclinicity_use_geostrophic_shear
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: geostrophic_wind_x, geostrophic_wind_y
  real(kind=DEFAULT_PRECISION) :: fcoriol

  ! Local tendency diagnostic variables for this component
  logical :: l_tend_3d_u, l_tend_3d_v
  ! Local mean tendency profile fields and logicals for their use
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v

  public initialisation_callback_coriolis, timestep_callback_coriolis, &
         finalisation_callback_coriolis

  contains

  !> Initialisation call back which will read in the coriolis configuration and set up the geostrophic winds
  !! @param current_state The current model state
  subroutine initialisation_callback_coriolis(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, iter
    integer :: alloc_z, alloc_y, alloc_x
    logical:: logicnum
    real(kind=DEFAULT_PRECISION) :: realnum

    !baroclinicity_use_geostrophic_shear=options_get_logical(current_state%options_database, "baroclinicity_use_geostrophic_shear")
    !fcoriol=options_get_real(current_state%options_database, "fcoriol")
    !current_state%geostrophic_wind_rate_of_change_in_x=options_get_real(current_state%options_database, &
    !     "geostrophic_wind_rate_of_change_in_x")
    !current_state%geostrophic_wind_rate_of_change_in_y=options_get_real(current_state%options_database, &
    !     "geostrophic_wind_rate_of_change_in_y")
    !current_state%surface_geostrophic_wind_x=options_get_real(current_state%options_database, "surface_geostrophic_wind_x")
    !current_state%surface_geostrophic_wind_y=options_get_real(current_state%options_database, "surface_geostrophic_wind_y")
    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "baroclinicity_use_geostrophic_shear") then
        read(current_state%options_database_string(iter,2),*) logicnum
        baroclinicity_use_geostrophic_shear = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "fcoriol") then
        read(current_state%options_database_string(iter,2),*) realnum
        fcoriol = realnum
      else if (current_state%options_database_string(iter,1) .eq. "geostrophic_wind_rate_of_change_in_x") then
        read(current_state%options_database_string(iter,2),*) realnum
        current_state%geostrophic_wind_rate_of_change_in_x = realnum
      else if (current_state%options_database_string(iter,1) .eq. "geostrophic_wind_rate_of_change_in_y") then
        read(current_state%options_database_string(iter,2),*) realnum
        current_state%geostrophic_wind_rate_of_change_in_y = realnum
      else if (current_state%options_database_string(iter,1) .eq. "surface_geostrophic_wind_x") then
        read(current_state%options_database_string(iter,2),*) realnum
        current_state%surface_geostrophic_wind_x = realnum
      else if (current_state%options_database_string(iter,1) .eq. "surface_geostrophic_wind_y") then
        read(current_state%options_database_string(iter,2),*) realnum
        current_state%surface_geostrophic_wind_y = realnum
      end if
    end do

    allocate(geostrophic_wind_x(current_state%local_grid%size(Z_INDEX)), &
         geostrophic_wind_y(current_state%local_grid%size(Z_INDEX)))

    do k=1,current_state%local_grid%size(Z_INDEX)
      geostrophic_wind_x(k)=current_state%surface_geostrophic_wind_x
      geostrophic_wind_y(k)=current_state%surface_geostrophic_wind_y
      if (baroclinicity_use_geostrophic_shear) then
        geostrophic_wind_x(k)=geostrophic_wind_x(k)+current_state%geostrophic_wind_rate_of_change_in_x*&
             current_state%global_grid%configuration%vertical%zn(k)
        geostrophic_wind_y(k)=geostrophic_wind_y(k)+current_state%geostrophic_wind_rate_of_change_in_y*&
             current_state%global_grid%configuration%vertical%zn(k)
      end if
    end do

    alloc_z=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    alloc_y=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    alloc_x=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

    ! Tendency Logicals
    l_tend_pr_tot_u   = current_state%u%active
    l_tend_pr_tot_v   = current_state%v%active

    l_tend_3d_u   = current_state%u%active .or. l_tend_pr_tot_u
    l_tend_3d_v   = current_state%v%active .or. l_tend_pr_tot_v

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_u) then
      !allocate( current_state%tend_3d_u_corio(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_u_corio(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_v) then
      !allocate( current_state%tend_3d_v_corio(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_v_corio(alloc_z, alloc_y, alloc_x))
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_u) then
      !allocate( current_state%tend_pr_tot_u_corio(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_u_corio(alloc_z))
    endif
    if (l_tend_pr_tot_v) then
      !allocate( current_state%tend_pr_tot_v_corio(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_v_corio(alloc_z))
    endif

  end subroutine initialisation_callback_coriolis


  subroutine finalisation_callback_coriolis(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(current_state%tend_3d_u_corio)) deallocate(current_state%tend_3d_u_corio)
    if (allocated(current_state%tend_3d_v_corio)) deallocate(current_state%tend_3d_v_corio)

    if (allocated(current_state%tend_pr_tot_u_corio)) deallocate(current_state%tend_pr_tot_u_corio)
    if (allocated(current_state%tend_pr_tot_v_corio)) deallocate(current_state%tend_pr_tot_v_corio)

  end subroutine finalisation_callback_coriolis

  !> For each none halo cell this will calculate the coriolis terms for su and sv fields
  !! @param current_state The current model state
  subroutine timestep_callback_coriolis(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: local_y, local_x, k, target_x_index, target_y_index
    logical :: calculate_diagnostics

    !calculate_diagnostics = current_state%diagnostic_sample_timestep &
    !                        .and. .not. current_state%halo_column

    local_y=current_state%column_local_y
    local_x=current_state%column_local_x
    target_y_index=local_y-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=local_x-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_u) then
        current_state%tend_pr_tot_u_corio(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_v) then
        current_state%tend_pr_tot_v_corio(:)= 0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    if (current_state%halo_column) then
      if (.not. ((current_state%column_local_y == current_state%local_grid%halo_size(Y_INDEX) .and. &
           current_state%column_local_x .le. current_state%local_grid%local_domain_end_index(X_INDEX) .and. &
           current_state%column_local_x .ge. current_state%local_grid%local_domain_start_index(X_INDEX)-1) .or. &
           (current_state%column_local_x == current_state%local_grid%halo_size(X_INDEX) .and. &
           current_state%column_local_y .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) &
           .and. current_state%column_local_y .le. current_state%local_grid%local_domain_end_index(Y_INDEX)) )) return
    end if

    if ((current_state%modulo_number_3d .eq. 0) &
        .and. (.not. current_state%halo_column)) then
      call save_precomponent_tendencies(current_state, local_x, local_y)!, target_x_index, target_y_index)
    end if
    
    do k=2,current_state%local_grid%size(Z_INDEX)
#if defined(U_ACTIVE) && defined(V_ACTIVE)
      current_state%su%data(k, current_state%column_local_y, current_state%column_local_x)=&
           current_state%su%data(k, current_state%column_local_y, current_state%column_local_x)+fcoriol*&
           (0.25_DEFAULT_PRECISION*(current_state%v%data(k, current_state%column_local_y, current_state%column_local_x)+&
           current_state%v%data(k, current_state%column_local_y, current_state%column_local_x+1)+&
           current_state%v%data(k, current_state%column_local_y-1, current_state%column_local_x)+&
           current_state%v%data(k, current_state%column_local_y-1, current_state%column_local_x+1))+current_state%vgal-&
           geostrophic_wind_y(k))

      current_state%sv%data(k, current_state%column_local_y, current_state%column_local_x)=&
           current_state%sv%data(k, current_state%column_local_y, current_state%column_local_x)-fcoriol*&
           (0.25_DEFAULT_PRECISION*(current_state%u%data(k, current_state%column_local_y, current_state%column_local_x)+&
           current_state%u%data(k, current_state%column_local_y, current_state%column_local_x-1)+&
           current_state%u%data(k, current_state%column_local_y+1, current_state%column_local_x)+&
           current_state%u%data(k, current_state%column_local_y+1, current_state%column_local_x-1))+current_state%ugal-&
           geostrophic_wind_x(k))

#endif
    end do

    if ((current_state%modulo_number_3d .eq. 0) &
        .and. (.not. current_state%halo_column)) then
      call compute_component_tendencies(current_state, local_x, local_y)!, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback_coriolis


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
    if (l_tend_3d_u) then
      !current_state%tend_3d_u_corio(:,tyn,txn)=current_state%su%data(:,cyn,cxn)
      current_state%tend_3d_u_corio(:,cyn,cxn)=current_state%su%data(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      !current_state%tend_3d_v_corio(:,tyn,txn)=current_state%sv%data(:,cyn,cxn)
      current_state%tend_3d_v_corio(:,cyn,cxn)=current_state%sv%data(:,cyn,cxn)
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
    if (l_tend_3d_u) then
      current_state%tend_3d_u_corio(:,cyn,cxn)=current_state%su%data(:,cyn,cxn) - &
             current_state%tend_3d_u_corio(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      current_state%tend_3d_v_corio(:,cyn,cxn)=current_state%sv%data(:,cyn,cxn) - &
             current_state%tend_3d_v_corio(:,cyn,cxn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_u) then
      current_state%tend_pr_tot_u_corio(:)=current_state%tend_pr_tot_u_corio(:) + &
            current_state%tend_3d_u_corio(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_v) then
      current_state%tend_pr_tot_v_corio(:)=current_state%tend_pr_tot_v_corio(:) + &
            current_state%tend_3d_v_corio(:,cyn,cxn)
    endif

  end subroutine compute_component_tendencies


end module coriolis_mod
