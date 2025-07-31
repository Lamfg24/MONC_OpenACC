!> Implements TVD advection for prognostic fields
module tvdadvection_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type, &
      component_descriptor_type_v1
  use stencil_mod, only : grid_stencil_type, interpolate_to_dual, create_stencil, free_stencil
  use state_mod, only : model_state_type, parallel_state_type, FORWARD_STEPPING
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use prognostics_mod, only : prognostic_field_type, prognostic_field_ptr_type
  use ultimateflux_mod, only : ultflx
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use optionsdatabase_mod, only : options_get_string, options_get_integer
  use collections_mod, only : map_type
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUS_IGNORE
  use q_indices_mod, only: get_q_index, standard_q_names
  ! Some tvd diagnostic terms 
  use def_tvd_diagnostic_terms, only: tvd_dgs_terms 

  implicit none

#ifndef TEST_MODE
  private
#endif

  type(grid_stencil_type), save :: star_stencil
  integer, save :: u_index=0, v_index=0, w_index=0
  logical :: advect_flow, advect_th, advect_q, advect_tracer
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: flux_x, flux_y, flux_z, u_advection, v_advection, &
       w_advection, th_advection, qv_advection, ql_advection, qr_advection, qi_advection, qs_advection, &
       qg_advection, qAitkenSolMass_advection, qAccumSolMass_advection, qAccumInsolMass_advection, &
       qCoarseSolMass_advection, qCoarseDustMass_advection, &
       nl_advection, nr_advection, ni_advection, ns_advection, ng_advection, &
       nAitkenSolNumber_advection, nAccumSolNumber_advection, nAccumInsolNumber_advection, &
       nCoarseSolNumber_advection, nCoarseDustnumber_advection
  !real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_advection
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: tracer_advection
  
  type(prognostic_field_type), dimension(:), allocatable :: interpolated_fields

  ! Local tendency diagnostic variables for this component
  logical :: l_tend_3d_u, l_tend_3d_v, l_tend_3d_w, l_tend_3d_th,l_tend_3d_qv,       &
             l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,       &
             l_tend_3d_tabs
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w,l_tend_pr_tot_th,l_tend_pr_tot_qv,        &
             l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,l_tend_pr_tot_qs,l_tend_pr_tot_qg,       &
             l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0

  public initialisation_callback_tvdadvection, timestep_callback_tvdadvection, &
         finalisation_callback_tvdadvection

contains

  !> Sets up the stencil_mod (used in interpolation) and allocates data for the flux fields
  !! @param current_state The current model state_mod
  subroutine initialisation_callback_tvdadvection(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    logical :: l_qdiag

    type(prognostic_field_ptr_type), dimension(3) :: fields
    integer, dimension(3, 2) :: sizes
    integer :: num_fields, iter
    integer :: alloc_z, alloc_y, alloc_x
    logical :: xdim, ydim

    xdim=.false.
    ydim=.false.
    num_fields=0
#ifdef U_ACTIVE    
    xdim=.true.
    num_fields = num_fields + 1
    fields(num_fields)%ptr => current_state%u
    sizes(num_fields,:) = (/ 2, 2 /) ! need um2 therefore -2 (applies to all interpolations)
    u_index = num_fields    
#endif

#ifdef V_ACTIVE
    ydim=.true.
    num_fields = num_fields + 1     
    fields(num_fields)%ptr => current_state%v
    sizes(num_fields,:) = (/ 1, 1 /)
    v_index=num_fields
#endif

#ifdef W_ACTIVE    
    num_fields = num_fields + 1
    fields(num_fields)%ptr => current_state%w
    sizes(num_fields,:) = (/ 1, 1 /)
    w_index=num_fields
#endif
    ! Allocate from 0, as any inactive dimensions will issue 0 to the ultimate flux which ignores the field
    allocate(interpolated_fields(0:num_fields))
#ifdef U_ACTIVE
    allocate(interpolated_fields(u_index)%data(current_state%global_grid%size(Z_INDEX), -1:3, -1:3))
    interpolated_fields(u_index)%active=.true.
#endif
#ifdef V_ACTIVE
    allocate(interpolated_fields(v_index)%data(current_state%global_grid%size(Z_INDEX), 0:2, 0:2))
    interpolated_fields(v_index)%active=.true.
#endif
#ifdef W_ACTIVE
    allocate(interpolated_fields(w_index)%data(current_state%global_grid%size(Z_INDEX), 0:2, 0:2))
    interpolated_fields(w_index)%active=.true.
#endif

    star_stencil = create_stencil(current_state%local_grid, fields, num_fields, 3, sizes, xdim, ydim)
    allocate(flux_y(current_state%global_grid%size(Z_INDEX)))
    allocate(flux_z(current_state%global_grid%size(Z_INDEX)))
    allocate(flux_x(current_state%global_grid%size(Z_INDEX)))
    allocate(u_advection(current_state%global_grid%size(Z_INDEX)), v_advection(current_state%global_grid%size(Z_INDEX)), &
         w_advection(current_state%global_grid%size(Z_INDEX)), th_advection(current_state%global_grid%size(Z_INDEX)), &
         qv_advection(current_state%global_grid%size(Z_INDEX)), ql_advection(current_state%global_grid%size(Z_INDEX)), &
         qr_advection(current_state%global_grid%size(Z_INDEX)), qi_advection(current_state%global_grid%size(Z_INDEX)), &
         qs_advection(current_state%global_grid%size(Z_INDEX)), qg_advection(current_state%global_grid%size(Z_INDEX)), &
         qAitkenSolMass_advection(current_state%global_grid%size(Z_INDEX)), &
         qAccumSolMass_advection(current_state%global_grid%size(Z_INDEX)), &
         qAccumInsolMass_advection(current_state%global_grid%size(Z_INDEX)), &
         qCoarseSolMass_advection(current_state%global_grid%size(Z_INDEX)), &
         qCoarseDustMass_advection(current_state%global_grid%size(Z_INDEX)), &
         nl_advection(current_state%global_grid%size(Z_INDEX)), nr_advection(current_state%global_grid%size(Z_INDEX)), &
         ni_advection(current_state%global_grid%size(Z_INDEX)), ns_advection(current_state%global_grid%size(Z_INDEX)), &
         ng_advection(current_state%global_grid%size(Z_INDEX)), &
         nAitkenSolNumber_advection(current_state%global_grid%size(Z_INDEX)), &
         nAccumSolNumber_advection(current_state%global_grid%size(Z_INDEX)), &
         nAccumInsolNumber_advection(current_state%global_grid%size(Z_INDEX)), &
         nCoarseSolNumber_advection(current_state%global_grid%size(Z_INDEX)), &
         nCoarseDustnumber_advection(current_state%global_grid%size(Z_INDEX)), &
         !q_advection(current_state%global_grid%size(Z_INDEX), current_state%number_q_fields), &
         tracer_advection(current_state%global_grid%size(Z_INDEX), current_state%n_tracers))

    ! Initialise terms
    flux_y(:) = 0.0_DEFAULT_PRECISION
    flux_x(:) = 0.0_DEFAULT_PRECISION
    flux_z(:) = 0.0_DEFAULT_PRECISION

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "advection_flow_fields") then
        advect_flow=determine_if_advection_here(current_state%options_database_string(iter,2))
      else if (current_state%options_database_string(iter,1) .eq. "advection_theta_field") then
        advect_th=determine_if_advection_here(current_state%options_database_string(iter,2))
        advect_tracer=determine_if_advection_here(current_state%options_database_string(iter,2))
      else if (current_state%options_database_string(iter,1) .eq. "advection_q_fields") then
        advect_q=determine_if_advection_here(current_state%options_database_string(iter,2))
      end if
    end do

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available

    alloc_z=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    alloc_y=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    alloc_x=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

    l_qdiag =  (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0)! .and. advect_q

    l_tend_pr_tot_u   = current_state%u%active ! .and. advect_flow
    l_tend_pr_tot_v   = current_state%v%active !.and. advect_flow
    l_tend_pr_tot_w   = current_state%w%active !.and. advect_flow
    l_tend_pr_tot_th  = current_state%th%active !.and. advect_th
    l_tend_pr_tot_qv  = l_qdiag .and. current_state%water_vapour_mixing_ratio_index > 0
    l_tend_pr_tot_ql  = l_qdiag .and. current_state%liquid_water_mixing_ratio_index > 0
    l_tend_pr_tot_qi  = l_qdiag .and. current_state%ice_water_mixing_ratio_index > 0
    l_tend_pr_tot_qr  = l_qdiag .and. current_state%rain_water_mixing_ratio_index > 0
    l_tend_pr_tot_qs  = l_qdiag .and. current_state%snow_water_mixing_ratio_index > 0
    l_tend_pr_tot_qg  = l_qdiag .and. current_state%graupel_water_mixing_ratio_index > 0
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_u   = current_state%u%active .or. l_tend_pr_tot_u!(current_state%u%active .and. advect_flow) .or. l_tend_pr_tot_u
    l_tend_3d_v   = current_state%v%active .or. l_tend_pr_tot_v!(current_state%v%active .and. advect_flow) .or. l_tend_pr_tot_v
    l_tend_3d_w   = current_state%w%active .or. l_tend_pr_tot_w!(current_state%w%active .and. advect_flow) .or. l_tend_pr_tot_w
    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th!(current_state%th%active .and. advect_th)  .or. l_tend_pr_tot_th
    l_tend_3d_qv  = (l_qdiag .and. current_state%water_vapour_mixing_ratio_index > 0) .or. l_tend_pr_tot_qv
    l_tend_3d_ql  = (l_qdiag .and. current_state%liquid_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_ql
    l_tend_3d_qi  = (l_qdiag .and. current_state%ice_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qi
    l_tend_3d_qr  = (l_qdiag .and. current_state%rain_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qr
    l_tend_3d_qs  = (l_qdiag .and. current_state%snow_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qs
    l_tend_3d_qg  = (l_qdiag .and. current_state%graupel_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qg
    l_tend_3d_tabs = l_tend_3d_th


    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_u) then
      !allocate( current_state%tend_3d_u_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_u_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_v) then
      !allocate( current_state%tend_3d_v_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_v_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_w) then
      !allocate( current_state%tend_3d_w_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_w_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_th) then
      !allocate( current_state%tend_3d_th_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_th_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qv) then
      !iqv= 1 !get_q_index(standard_q_names%VAPOUR, 'tvd_advection')
      !allocate( current_state%tend_3d_qv_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qv_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_ql) then
      !iql= 2 !get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'tvd_advection')
      !allocate( current_state%tend_3d_ql_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_ql_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qi) then
      !iqi= 3 !get_q_index(standard_q_names%ICE_MASS, 'tvd_advection')
      !allocate( current_state%tend_3d_qi_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qi_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qr) then
      !iqr= 4 !get_q_index(standard_q_names%RAIN_MASS, 'tvd_advection')
      !allocate( current_state%tend_3d_qr_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qr_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qs) then
      !iqs= 5 !get_q_index(standard_q_names%SNOW_MASS, 'tvd_advection')
      !allocate( current_state%tend_3d_qs_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qs_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qg) then
      !iqg= 6 !get_q_index(standard_q_names%GRAUPEL_MASS, 'tvd_advection')
      !allocate( current_state%tend_3d_qg_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qg_tvad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_tabs) then
      !allocate( current_state%tend_3d_tabs_tvad(current_state%local_grid%size(Z_INDEX),  &
      !                       current_state%local_grid%size(Y_INDEX),  &
      !                       current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_tabs_tvad(alloc_z, alloc_y, alloc_x))
    endif
    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_u) then
      !allocate( current_state%tend_pr_tot_u_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_u_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_v) then
      !allocate( current_state%tend_pr_tot_v_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_v_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_w) then
      !allocate( current_state%tend_pr_tot_w_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_w_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_th) then
      !allocate( current_state%tend_pr_tot_th_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_th_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_qv) then
      !allocate( current_state%tend_pr_tot_qv_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qv_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_ql) then
      !allocate( current_state%tend_pr_tot_ql_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_ql_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_qi) then
      !allocate( current_state%tend_pr_tot_qi_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qi_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_qr) then
      !allocate( current_state%tend_pr_tot_qr_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qr_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_qs) then
      !allocate( current_state%tend_pr_tot_qs_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qs_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_qg) then
      !allocate( current_state%tend_pr_tot_qg_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qg_tvad(alloc_z))
    endif
    if (l_tend_pr_tot_tabs) then
      !allocate( current_state%tend_pr_tot_tabs_tvad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_tabs_tvad(alloc_z))
    endif

  end subroutine initialisation_callback_tvdadvection

  !> Frees up the memory associated with the advection
  !! @param current_state The current model state
  subroutine finalisation_callback_tvdadvection(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call free_stencil(star_stencil)
    if (allocated(flux_x)) deallocate(flux_x)
    if (allocated(flux_y)) deallocate(flux_y)
    if (allocated(flux_z)) deallocate(flux_z)
    if (allocated(interpolated_fields)) deallocate(interpolated_fields)
    if (allocated(u_advection)) deallocate(u_advection)
    if (allocated(v_advection)) deallocate(v_advection)
    if (allocated(w_advection)) deallocate(w_advection)
    if (allocated(th_advection)) deallocate(th_advection)
    !if (allocated(q_advection)) deallocate(q_advection)
    if (allocated(tracer_advection)) deallocate(tracer_advection)

    if (allocated(current_state%tend_3d_u_tvad)) deallocate(current_state%tend_3d_u_tvad)
    if (allocated(current_state%tend_3d_v_tvad)) deallocate(current_state%tend_3d_v_tvad)
    if (allocated(current_state%tend_3d_w_tvad)) deallocate(current_state%tend_3d_w_tvad)
    if (allocated(current_state%tend_3d_th_tvad)) deallocate(current_state%tend_3d_th_tvad)
    if (allocated(current_state%tend_3d_qv_tvad)) deallocate(current_state%tend_3d_qv_tvad)
    if (allocated(current_state%tend_3d_ql_tvad)) deallocate(current_state%tend_3d_ql_tvad)
    if (allocated(current_state%tend_3d_qi_tvad)) deallocate(current_state%tend_3d_qi_tvad)
    if (allocated(current_state%tend_3d_qr_tvad)) deallocate(current_state%tend_3d_qr_tvad)
    if (allocated(current_state%tend_3d_qs_tvad)) deallocate(current_state%tend_3d_qs_tvad)
    if (allocated(current_state%tend_3d_qg_tvad)) deallocate(current_state%tend_3d_qg_tvad)
    if (allocated(current_state%tend_3d_tabs_tvad)) deallocate(current_state%tend_3d_tabs_tvad)

    if (allocated(current_state%tend_pr_tot_u_tvad)) deallocate(current_state%tend_pr_tot_u_tvad)
    if (allocated(current_state%tend_pr_tot_v_tvad)) deallocate(current_state%tend_pr_tot_v_tvad)
    if (allocated(current_state%tend_pr_tot_w_tvad)) deallocate(current_state%tend_pr_tot_w_tvad)
    if (allocated(current_state%tend_pr_tot_th_tvad)) deallocate(current_state%tend_pr_tot_th_tvad)
    if (allocated(current_state%tend_pr_tot_qv_tvad)) deallocate(current_state%tend_pr_tot_qv_tvad)
    if (allocated(current_state%tend_pr_tot_ql_tvad)) deallocate(current_state%tend_pr_tot_ql_tvad)
    if (allocated(current_state%tend_pr_tot_qi_tvad)) deallocate(current_state%tend_pr_tot_qi_tvad)
    if (allocated(current_state%tend_pr_tot_qr_tvad)) deallocate(current_state%tend_pr_tot_qr_tvad)
    if (allocated(current_state%tend_pr_tot_qs_tvad)) deallocate(current_state%tend_pr_tot_qs_tvad)
    if (allocated(current_state%tend_pr_tot_qg_tvad)) deallocate(current_state%tend_pr_tot_qg_tvad)
    if (allocated(current_state%tend_pr_tot_tabs_tvad)) deallocate(current_state%tend_pr_tot_tabs_tvad)

  end subroutine finalisation_callback_tvdadvection

  !> Timestep callback hook which performs the TVD advection for each prognostic field
  !! @param current_state The current model state_mod
  subroutine timestep_callback_tvdadvection(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: current_x_index, current_y_index, target_x_index, target_y_index, iter
    logical :: calculate_diagnostics, profile_diagnostics_enabled, logicnum

    !calculate_diagnostics = current_state%diagnostic_sample_timestep &
    !                        .and. .not. current_state%halo_column

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_u) then
        current_state%tend_pr_tot_u_tvad(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_v) then
        current_state%tend_pr_tot_v_tvad(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_w) then
        current_state%tend_pr_tot_w_tvad(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_th) then
        current_state%tend_pr_tot_th_tvad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qv) then
        current_state%tend_pr_tot_qv_tvad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_ql) then
        current_state%tend_pr_tot_ql_tvad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qi) then
        current_state%tend_pr_tot_qi_tvad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qr) then
        current_state%tend_pr_tot_qr_tvad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qs) then
        current_state%tend_pr_tot_qs_tvad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qg) then
        current_state%tend_pr_tot_qg_tvad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_tabs) then
        current_state%tend_pr_tot_tabs_tvad(:)=0.0_DEFAULT_PRECISION
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
        call save_precomponent_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)
    end if

    do iter = 1,current_state%config_args ! to review
      if (current_state%options_database_string(iter,1) .eq. "profile_diagnostics_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        profile_diagnostics_enabled = logicnum
      end if
    end do

    if (advect_flow) call advect_flow_fields(current_state, profile_diagnostics_enabled)
    if (advect_th) call advect_theta(current_state, profile_diagnostics_enabled)

    if (advect_q) call advect_q_fields(current_state, profile_diagnostics_enabled)
    if (advect_tracer) call advect_tracer_fields(current_state)

    if ((current_state%modulo_number_3d .eq. 0) &
        .and. (.not. current_state%halo_column)) then
        call compute_component_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback_tvdadvection

  !> Will advect the flow fields
  !! @param current_state The current model state_mod
  subroutine advect_flow_fields(current_state, profile_diagnostics_enabled)
    type(model_state_type), intent(inout) :: current_state
    logical, intent(in) :: profile_diagnostics_enabled

    real(kind=DEFAULT_PRECISION) :: dtm

    dtm = current_state%dtm*2.0_DEFAULT_PRECISION
    if (current_state%momentum_stepping == FORWARD_STEPPING) dtm=current_state%dtm

#ifdef U_ACTIVE
    call advect_u(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
         current_state%v, current_state%w, current_state%zu, current_state%su, current_state%global_grid, &
         current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
    if (profile_diagnostics_enabled) then ! à revoir
       ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument 
       !       list in advect_scalar_field.
       tvd_dgs_terms%adv_u_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
            flux_z(:)
    endif
#endif

#ifdef V_ACTIVE
    call advect_v(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
         current_state%v, current_state%w, current_state%zv, current_state%sv, current_state%global_grid, &
         current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
    if (profile_diagnostics_enabled) then !  a revoir
       ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument 
       !       list in advect_scalar_field.
       tvd_dgs_terms%adv_v_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
            flux_z(:)
    endif
#endif

#ifdef W_ACTIVE
    call advect_w(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
         current_state%v, current_state%w, current_state%zw, current_state%sw, current_state%global_grid,&
         current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
    if (profile_diagnostics_enabled) then ! a revoir
       ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument 
       !       list in advect_scalar_field.
       tvd_dgs_terms%adv_w_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
            flux_z(:)
    endif
#endif
  end subroutine advect_flow_fields

  !> Advects the Q fields
  !! @param current_state The current model state_mod
  subroutine advect_q_fields(current_state, profile_diagnostics_enabled)
    type(model_state_type), intent(inout) :: current_state
    logical, intent(in) :: profile_diagnostics_enabled

    integer :: i, iter, intnum, q_active
    real(kind=DEFAULT_PRECISION) :: dtm

    dtm = current_state%dtm*2.0_DEFAULT_PRECISION
    if (current_state%scalar_stepping == FORWARD_STEPPING) dtm=current_state%dtm

    do iter = 1,current_state%config_args ! to review
      if (current_state%options_database_string(iter,1) .eq. "number_q_fields") then
        read(current_state%options_database_string(iter,2),*) intnum
        q_active = intnum
      end if
    end do

    if (q_active > 0) then
      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqv, current_state%qv, current_state%sqv, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qv_advection=current_state%sqv%data(:, current_state%column_local_y, current_state%column_local_x)
      if (profile_diagnostics_enabled) then
        tvd_dgs_terms%adv_qv_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
              flux_z(:)
      end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zql, current_state%ql, current_state%sql, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      ql_advection=current_state%sql%data(:, current_state%column_local_y, current_state%column_local_x)
      if (profile_diagnostics_enabled) then
        tvd_dgs_terms%adv_ql_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
              flux_z(:)
      end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqr, current_state%qr, current_state%sqr, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qr_advection=current_state%sqr%data(:, current_state%column_local_y, current_state%column_local_x)
      if (profile_diagnostics_enabled) then
        tvd_dgs_terms%adv_qr_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
              flux_z(:)
      end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqi, current_state%qi, current_state%sqi, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qi_advection=current_state%sqi%data(:, current_state%column_local_y, current_state%column_local_x)
      if (profile_diagnostics_enabled) then
        tvd_dgs_terms%adv_qi_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
              flux_z(:)
      end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqs, current_state%qs, current_state%sqs, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qs_advection=current_state%sqs%data(:, current_state%column_local_y, current_state%column_local_x)
      if (profile_diagnostics_enabled) then
        tvd_dgs_terms%adv_qs_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
              flux_z(:)
      end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqg, current_state%qg, current_state%sqg, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qg_advection=current_state%sqg%data(:, current_state%column_local_y, current_state%column_local_x)
      if (profile_diagnostics_enabled) then
        tvd_dgs_terms%adv_qg_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
              flux_z(:)
      end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqAitkenSolMass, &
            current_state%qAitkenSolMass, current_state%sqAitkenSolMass, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qAitkenSolMass_advection=current_state%sqAitkenSolMass%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qAitkenSolMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqAccumSolMass, &
            current_state%qAccumSolMass, current_state%sqAccumSolMass, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qAccumSolMass_advection=current_state%sqAccumSolMass%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qAccumSolMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqAccumInsolMass, &
            current_state%qAccumInsolMass, current_state%sqAccumInsolMass, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qAccumInsolMass_advection=current_state%sqAccumInsolMass%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qAccumInsolMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqCoarseSolMass, &
            current_state%qCoarseSolMass, current_state%sqCoarseSolMass, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qCoarseSolMass_advection=current_state%sqCoarseSolMass%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseSolMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zqCoarseDustMass, &
            current_state%qCoarseDustMass, current_state%sqCoarseDustMass, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      qCoarseDustMass_advection=current_state%sqCoarseDustMass%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if



      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%znl, &
            current_state%nl, current_state%snl, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      nl_advection=current_state%snl%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%znr, &
            current_state%nr, current_state%snr, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      nr_advection=current_state%snr%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zni, &
            current_state%ni, current_state%sni, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      ni_advection=current_state%sni%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zns, &
            current_state%ns, current_state%sns, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      ns_advection=current_state%sns%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%zng, &
            current_state%ng, current_state%sng, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      ng_advection=current_state%sng%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%znAitkenSolNumber, &
            current_state%nAitkenSolNumber, current_state%snAitkenSolNumber, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      nAitkenSolNumber_advection=current_state%snAitkenSolNumber%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%znAccumSolNumber, &
            current_state%nAccumSolNumber, current_state%snAccumSolNumber, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      nAccumSolNumber_advection=current_state%snAccumSolNumber%data(:, current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%znAccumInsolNumber, &
            current_state%nAccumInsolNumber, current_state%snAccumInsolNumber, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      nAccumInsolNumber_advection=current_state%snAccumInsolNumber%data(:, &
      current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%znCoarseSolNumber, &
            current_state%nCoarseSolNumber, current_state%snCoarseSolNumber, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      nCoarseSolNumber_advection=current_state%snCoarseSolNumber%data(:, &
      current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if

      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
            current_state%v, current_state%w, current_state%znCoarseDustnumber, &
            current_state%nCoarseDustnumber, current_state%snCoarseDustnumber, &
            current_state%global_grid, current_state%local_grid, current_state%parallel, &
            current_state%halo_column, current_state%field_stepping)
      nCoarseDustnumber_advection=current_state%snCoarseDustnumber%data(:, &
      current_state%column_local_y, current_state%column_local_x)
      !if (profile_diagnostics_enabled) then
      !  tvd_dgs_terms%adv_qCoarseDustMass_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
      !        flux_z(:)
      !end if
    end if

    !do i=1,current_state%number_q_fields
    !  if (current_state%q(i)%active) then
    !    call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
    !         current_state%v, current_state%w, current_state%zq(i), current_state%q(i), current_state%sq(i), &
    !         current_state%global_grid, current_state%local_grid, current_state%parallel, &
    !         current_state%halo_column, current_state%field_stepping)
    !    q_advection(:,i)=current_state%sq(i)%data(:, current_state%column_local_y, current_state%column_local_x)
        !if (is_component_enabled(current_state%options_database, "profile_diagnostics")) then !
        !   ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument
        !   !       list in advect_scalar_field.
        !   tvd_dgs_terms%adv_q_dgs(:, current_state%column_local_y, current_state%column_local_x, i) =  &
        !        flux_z(:)
        !endif
      !end if
    !end do
  end subroutine advect_q_fields

  !> Advects the tracer fields
  !! @param current_state The current model state_mod
  subroutine advect_tracer_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i
    real(kind=DEFAULT_PRECISION) :: dtm

    dtm = current_state%dtm*2.0_DEFAULT_PRECISION
    if (current_state%scalar_stepping == FORWARD_STEPPING) dtm=current_state%dtm

    do i=1,current_state%n_tracers
      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u, &
           current_state%v, current_state%w, current_state%ztracer(i), current_state%tracer(i), current_state%stracer(i), &
           current_state%global_grid, current_state%local_grid, current_state%parallel, &
           current_state%halo_column, current_state%field_stepping)
      tracer_advection(:,i)=current_state%stracer(i)%data(:, current_state%column_local_y, current_state%column_local_x)          
    end do
  end subroutine advect_tracer_fields

  !> Advects the theta field if it is active
  !! @param current_state The current model state_mod
  subroutine advect_theta(current_state, profile_diagnostics_enabled)
    type(model_state_type), intent(inout) :: current_state
    logical, intent(in) :: profile_diagnostics_enabled

    real(kind=DEFAULT_PRECISION) :: dtm

    dtm = current_state%dtm*2.0_DEFAULT_PRECISION
    if (current_state%scalar_stepping == FORWARD_STEPPING) dtm=current_state%dtm

    if (current_state%th%active) then
      call advect_scalar_field(current_state%column_local_y, current_state%column_local_x, dtm, current_state%u,&
           current_state%v, current_state%w, current_state%zth, current_state%th, current_state%sth, current_state%global_grid,&
           current_state%local_grid, current_state%parallel, current_state%halo_column, current_state%field_stepping)
      th_advection=current_state%sth%data(:, current_state%column_local_y, current_state%column_local_x)
      if (profile_diagnostics_enabled) then !  a revoir
           ! NOTE: flux_z is declared at the top of module and then passed into ultflx, through argument
           !       list in advect_scalar_field.
         tvd_dgs_terms%adv_th_dgs(:, current_state%column_local_y, current_state%column_local_x) =  &
             flux_z(:)
      endif
    end if
  end subroutine advect_theta

  !> Advects a single scalar field
  subroutine advect_scalar_field(y_local_index, x_local_index, dt, u, v, w, z_q_field, q_field, source_field, &
       global_grid, local_grid, parallel, halo_column, field_stepping)

    integer, intent(in) ::y_local_index, x_local_index, field_stepping
    real(kind=DEFAULT_PRECISION), intent(in) ::dt   ! timestep (s)
    logical, intent(in) :: halo_column
    type(prognostic_field_type), intent(inout) :: u, w, v, z_q_field, q_field, source_field
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    type(parallel_state_type), intent(inout) :: parallel

    if (.not. allocated(q_field%flux_previous_x)) allocate(q_field%flux_previous_x(local_grid%size(Z_INDEX), &
         local_grid%size(Y_INDEX)+4))
    if (.not. allocated(q_field%flux_previous_y)) allocate(q_field%flux_previous_y(local_grid%size(Z_INDEX)))

    call register_y_flux_wrap_recv_if_required(y_local_index, q_field, parallel, local_grid)

    if (field_stepping == FORWARD_STEPPING) then
      call ultflx(y_local_index, x_local_index, u, v, w, y_local_index, x_local_index, q_field, local_grid, &
           global_grid%configuration, parallel, 0, dt, &
           flux_y, flux_z, flux_x, q_field%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz,&
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    else
      call ultflx(y_local_index, x_local_index, u, v, w, y_local_index, x_local_index, z_q_field, local_grid, &
           global_grid%configuration, parallel, 0, dt, &
           flux_y, flux_z, flux_x, q_field%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz,&
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    end if
    call complete_y_flux_wrap_recv_if_required(y_local_index, q_field, parallel, local_grid)

    if (.not. halo_column) then
      call differentiate_face_values(y_local_index, x_local_index, u, v, w, y_local_index, x_local_index, source_field, &
           local_grid, global_grid, q_field%flux_previous_y, q_field%flux_previous_x(:,y_local_index), &
           global_grid%configuration%vertical%tzc1, global_grid%configuration%vertical%tzc2, .false.)
    end if

    if (y_local_index .lt. local_grid%local_domain_end_index(Y_INDEX)) q_field%flux_previous_y(:) = flux_y(:)
    call register_y_flux_wrap_send_if_required(y_local_index, q_field, parallel, local_grid)
    call complete_y_flux_wrap_send_if_required(y_local_index, q_field, parallel, local_grid)
  end subroutine advect_scalar_field

  !> Advects the U flow field
  subroutine advect_u(y_local_index, x_local_index, dt, u, v, w, zf, su, global_grid, local_grid, parallel, &
       halo_column, field_stepping)

    integer, intent(in) :: y_local_index, x_local_index, field_stepping
    real(kind=DEFAULT_PRECISION), intent(in) ::dt   ! timestep (s)
    logical, intent(in) :: halo_column
    type(prognostic_field_type), intent(inout) :: u, w, v, zf, su
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    type(parallel_state_type), intent(inout) :: parallel

    if (.not. allocated(u%flux_previous_x)) allocate(u%flux_previous_x(local_grid%size(Z_INDEX), local_grid%size(Y_INDEX)+4))
    if (.not. allocated(u%flux_previous_y)) allocate(u%flux_previous_y(local_grid%size(Z_INDEX)))

    call register_y_flux_wrap_recv_if_required(y_local_index, u, parallel, local_grid)

    call interpolate_to_dual(local_grid, u, star_stencil, x_local_index, y_local_index, interpolated_fields, u_index)

    if (field_stepping == FORWARD_STEPPING) then
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, u, local_grid, global_grid%configuration, parallel, 0, &
           dt, flux_y, flux_z, flux_x, &
           u%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz, &
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    else
      ! Centred stepping
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, zf, local_grid, global_grid%configuration, parallel, 0, &
           dt, flux_y, flux_z, flux_x, &
           u%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz, &
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    end if

    call complete_y_flux_wrap_recv_if_required(y_local_index, u, parallel, local_grid)

    if (.not. halo_column) then
      call differentiate_face_values(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, su, local_grid, global_grid, &
           u%flux_previous_y, u%flux_previous_x(:,y_local_index), &
           global_grid%configuration%vertical%tzc1, global_grid%configuration%vertical%tzc2, .true.)
      u_advection=su%data(:, y_local_index, x_local_index)
    end if

    if (y_local_index .lt. local_grid%local_domain_end_index(Y_INDEX)) u%flux_previous_y(:) = flux_y(:)
    call register_y_flux_wrap_send_if_required(y_local_index, u, parallel, local_grid)
    call complete_y_flux_wrap_send_if_required(y_local_index, u, parallel, local_grid)
  end subroutine advect_u

  !> Advects the V flow field
  subroutine advect_v(y_local_index, x_local_index, dt, u, v, w, zf, sv, global_grid, local_grid, parallel, &
       halo_column, field_stepping)

    integer, intent(in) ::y_local_index, x_local_index, field_stepping
    real(kind=DEFAULT_PRECISION), intent(in) ::dt   ! timestep (s)
    logical, intent(in) :: halo_column
    type(prognostic_field_type), intent(inout) :: u, w, v, zf, sv
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    type(parallel_state_type), intent(inout) :: parallel

    if (.not. allocated(v%flux_previous_x)) allocate(v%flux_previous_x(local_grid%size(Z_INDEX), local_grid%size(Y_INDEX)+4))
    if (.not. allocated(v%flux_previous_y)) allocate(v%flux_previous_y(local_grid%size(Z_INDEX)))

    call register_y_flux_wrap_recv_if_required(y_local_index, v, parallel, local_grid)

    call interpolate_to_dual(local_grid, v, star_stencil, x_local_index, y_local_index, interpolated_fields, v_index)

    if (field_stepping == FORWARD_STEPPING) then
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, v, local_grid, global_grid%configuration, parallel, 0, &
           dt, flux_y, flux_z, flux_x, &
           v%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz, &
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    else
      ! Centred stepping
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, zf, local_grid, global_grid%configuration, parallel, 0, &
           dt, flux_y, flux_z, flux_x, &
           v%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdz, &
           global_grid%configuration%vertical%rdzn, global_grid%configuration%vertical%dzn, 2, local_grid%size(Z_INDEX))
    end if

    call complete_y_flux_wrap_recv_if_required(y_local_index, v, parallel, local_grid)

    if (.not. halo_column) then
      call differentiate_face_values(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, sv, local_grid, global_grid, &
           v%flux_previous_y, v%flux_previous_x(:,y_local_index), &           
           global_grid%configuration%vertical%tzc1, global_grid%configuration%vertical%tzc2, .true.)
      v_advection=sv%data(:, y_local_index, x_local_index)
    end if

    if (y_local_index .lt. local_grid%local_domain_end_index(Y_INDEX)) v%flux_previous_y(:) = flux_y(:)
    call register_y_flux_wrap_send_if_required(y_local_index, v, parallel, local_grid)
    call complete_y_flux_wrap_send_if_required(y_local_index, v, parallel, local_grid)
  end subroutine advect_v

  !> Advects the W flow field
  subroutine advect_w(y_local_index, x_local_index, dt, u, v, w, zf, sw, global_grid, local_grid, parallel, &
       halo_column, field_stepping)

    integer, intent(in) ::y_local_index, x_local_index, field_stepping
    real(kind=DEFAULT_PRECISION), intent(in) ::dt   ! timestep (s)
    logical, intent(in) :: halo_column
    type(prognostic_field_type), intent(inout) :: u, w, v, zf, sw
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    type(parallel_state_type), intent(inout) :: parallel

    if (.not. allocated(w%flux_previous_x)) allocate(w%flux_previous_x(local_grid%size(Z_INDEX), local_grid%size(Y_INDEX)+4))
    if (.not. allocated(w%flux_previous_y)) allocate(w%flux_previous_y(local_grid%size(Z_INDEX)))

    call register_y_flux_wrap_recv_if_required(y_local_index, w, parallel, local_grid)

    call interpolate_to_dual(local_grid, w, star_stencil, x_local_index, y_local_index, interpolated_fields, w_index)

    if (field_stepping == FORWARD_STEPPING) then
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, w, local_grid, global_grid%configuration, parallel, 1, &
           dt, flux_y, flux_z, flux_x,&
           w%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdzn, &
           global_grid%configuration%vertical%rdz, global_grid%configuration%vertical%dz, 1, local_grid%size(Z_INDEX)-1)
    else
      ! Centred stepping
      call ultflx(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), &
           interpolated_fields(w_index), y_local_index, x_local_index, zf, local_grid, global_grid%configuration, parallel, 1, &
           dt, flux_y, flux_z, flux_x,&
           w%flux_previous_x(:,y_local_index), global_grid%configuration%vertical%rdzn, &
           global_grid%configuration%vertical%rdz, global_grid%configuration%vertical%dz, 1, local_grid%size(Z_INDEX)-1)
    end if

    call complete_y_flux_wrap_recv_if_required(y_local_index, w, parallel, local_grid)

    if (.not. halo_column) then
      call differentiate_face_values(1, 1, interpolated_fields(u_index), interpolated_fields(v_index),&
           interpolated_fields(w_index), y_local_index, x_local_index, sw, local_grid, global_grid, &
           w%flux_previous_y, w%flux_previous_x(:,y_local_index),&
           global_grid%configuration%vertical%tzd1, global_grid%configuration%vertical%tzd2, .false.)
      w_advection=sw%data(:, y_local_index, x_local_index)
    end if
    if (y_local_index .lt. local_grid%local_domain_end_index(Y_INDEX)) w%flux_previous_y(:) = flux_y(:)
    call register_y_flux_wrap_send_if_required(y_local_index, w, parallel, local_grid)
    call complete_y_flux_wrap_send_if_required(y_local_index, w, parallel, local_grid)
  end subroutine advect_w

  !> Differentiates face values to update the source field
  subroutine differentiate_face_values(y_flow_index, x_flow_index, u, v, w, y_source_index, x_source_index, source_field, &
       local_grid, global_grid, flux_y_previous, flux_x_previous, tzc1, tzc2, differentiate_top)

    integer, intent(in) :: y_flow_index, x_flow_index, y_source_index, x_source_index
    logical, intent(in) :: differentiate_top
    real(kind=DEFAULT_PRECISION), intent(in), dimension(*) :: tzc1, tzc2
    type(prognostic_field_type), intent(inout) :: u, w, v
    type(prognostic_field_type), intent(inout) :: source_field
    type(global_grid_type), intent(inout) :: global_grid
    type(local_grid_type), intent(inout) :: local_grid
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: flux_y_previous, flux_x_previous

    integer :: k

    do k=2,local_grid%size(Z_INDEX)-1
#ifdef V_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (v%data(k, y_flow_index-1, x_flow_index)* flux_y_previous(k) - v%data(k, y_flow_index, x_flow_index)*flux_y(k))*&
           global_grid%configuration%horizontal%cy
#endif
#ifdef W_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           4.0_DEFAULT_PRECISION*(w%data(k-1, y_flow_index, x_flow_index)* flux_z(k)*tzc1(k) - &
           w%data(k, y_flow_index, x_flow_index)*flux_z(k+1)*tzc2(k))
#endif
#ifdef U_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (u%data(k, y_flow_index, x_flow_index-1)* flux_x(k) - u%data(k, y_flow_index, x_flow_index)*flux_x_previous(k))*&
           global_grid%configuration%horizontal%cx
#endif
    end do
    if (differentiate_top) then
      k=local_grid%size(Z_INDEX)
      source_field%data(k, y_source_index, x_source_index)=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (v%data(k, y_flow_index-1, x_flow_index)* flux_y_previous(k) - v%data(k, y_flow_index, x_flow_index)*flux_y(k))*&
           global_grid%configuration%horizontal%cy
#endif
#ifdef W_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           4.0_DEFAULT_PRECISION*tzc1(k)* w%data(k-1, y_flow_index, x_flow_index)*flux_z(k)
#endif
#ifdef U_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (u%data(k, y_flow_index, x_flow_index-1)* flux_x(k) -u%data(k, y_flow_index, x_flow_index)*flux_x_previous(k))*&
           global_grid%configuration%horizontal%cx
#endif
    end if
  end subroutine differentiate_face_values

  !> Completes the Y flux MPI asynchronous send if required
  !! @param y_local_index The local index in Y
  !! @param field The prognostic field
  !! @param parallel Parallel system description
  !! @param local_grid The local grid description
  subroutine complete_y_flux_wrap_send_if_required(y_local_index, field, parallel, local_grid)
    integer, intent(in) :: y_local_index
    type(prognostic_field_type), intent(inout) :: field
    type(parallel_state_type), intent(inout) :: parallel
    type(local_grid_type), intent(inout) :: local_grid

    integer :: ierr

    if (y_local_index == local_grid%local_domain_end_index(Y_INDEX) .and. parallel%my_coords(Y_INDEX) == 0 .and. &
         field%async_flux_handle .ne. MPI_REQUEST_NULL) then
      call mpi_wait(field%async_flux_handle, MPI_STATUS_IGNORE, ierr)
    end if
  end subroutine complete_y_flux_wrap_send_if_required

  !> Registers an asynchronous send for the Y flux if required.
  !!
  !! This is done after the second y is computed and we have until the entire Y dimension is completed
  !! until the communication must be complete. If the wrap around process is the same (one process in Y dimension)
  !! then just issues a local copy to the buffer.
  !! @param y_local_index The local index in Y
  !! @param field The prognostic field
  !! @param parallel Parallel system description
  !! @param local_grid The local grid description
  subroutine register_y_flux_wrap_send_if_required(y_local_index, field, parallel, local_grid)
    integer, intent(in) :: y_local_index
    type(prognostic_field_type), intent(inout) :: field
    type(parallel_state_type), intent(inout) :: parallel
    type(local_grid_type), intent(inout) :: local_grid

    integer :: ierr

    if (y_local_index == local_grid%local_domain_start_index(Y_INDEX)-1 .and. parallel%my_coords(Y_INDEX) == 0) then
      if (.not. allocated(field%flux_y_buffer)) allocate(field%flux_y_buffer(local_grid%size(Z_INDEX)))
      field%flux_y_buffer(:) = flux_y(:)
      if (parallel%my_rank .ne. local_grid%neighbours(Y_INDEX,1)) then      
        call mpi_isend(field%flux_y_buffer, local_grid%size(Z_INDEX), PRECISION_TYPE, local_grid%neighbours(Y_INDEX,1), 0, &
             parallel%neighbour_comm, field%async_flux_handle, ierr)
      end if
    end if
  end subroutine register_y_flux_wrap_send_if_required

  !> Completes the Y flux MPI asynchronous recieve if required. If the wrap around process is the same (one process
  !! in the y dimension) then just issues a local copy
  !! @param y_local_index The local index in Y
  !! @param field The prognostic field
  !! @param parallel Parallel system description
  !! @param local_grid The local grid description
  subroutine complete_y_flux_wrap_recv_if_required(y_local_index, field, parallel, local_grid)
    integer, intent(in) :: y_local_index
    type(prognostic_field_type), intent(inout) :: field
    type(parallel_state_type), intent(inout) :: parallel
    type(local_grid_type), intent(inout) :: local_grid

    integer :: ierr

    if (y_local_index == local_grid%local_domain_end_index(Y_INDEX) .and. parallel%my_coords(Y_INDEX) == &
         parallel%dim_sizes(Y_INDEX)-1) then
      if (field%async_flux_handle .ne. MPI_REQUEST_NULL) then
        call mpi_wait(field%async_flux_handle, MPI_STATUS_IGNORE, ierr)
      end if
      flux_y(:) = field%flux_y_buffer(:)
    end if
  end subroutine complete_y_flux_wrap_recv_if_required

  !> Registers an MPI asynchronous receive for the flux if required.
  !!
  !! This is registered at the start and we have until the last column in Y until it must be completed. No 
  !! communication is registered if this is a local operation
  !! @param y_local_index The local index in Y
  !! @param field The prognostic field
  !! @param parallel Parallel system description
  !! @param local_grid The local grid description
  subroutine register_y_flux_wrap_recv_if_required(y_local_index, field, parallel, local_grid)
    integer, intent(in) :: y_local_index
    type(prognostic_field_type), intent(inout) :: field
    type(parallel_state_type), intent(inout) :: parallel
    type(local_grid_type), intent(inout) :: local_grid

    integer :: ierr

    if (y_local_index == local_grid%local_domain_start_index(Y_INDEX) .and. parallel%my_coords(Y_INDEX) == &
         parallel%dim_sizes(Y_INDEX)-1) then
      if (parallel%my_rank .ne. local_grid%neighbours(Y_INDEX,3)) then
        if (.not. allocated(field%flux_y_buffer)) allocate(field%flux_y_buffer(local_grid%size(Z_INDEX)))
        call mpi_irecv(field%flux_y_buffer, local_grid%size(Z_INDEX), PRECISION_TYPE, local_grid%neighbours(Y_INDEX,3), 0, &
             parallel%neighbour_comm, field%async_flux_handle, ierr)
      end if
    end if
  end subroutine register_y_flux_wrap_recv_if_required


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
      current_state%tend_3d_u_tvad(:,cyn,cxn)=current_state%su%data(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      current_state%tend_3d_v_tvad(:,cyn,cxn)=current_state%sv%data(:,cyn,cxn)
    endif
    if (l_tend_3d_w) then
      current_state%tend_3d_w_tvad(:,cyn,cxn)=current_state%sw%data(:,cyn,cxn)
    endif
    if (l_tend_3d_th) then
      current_state%tend_3d_th_tvad(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      current_state%tend_3d_qv_tvad(:,cyn,cxn)=current_state%sqv%data(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      current_state%tend_3d_ql_tvad(:,cyn,cxn)=current_state%sql%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qi) then
      current_state%tend_3d_qi_tvad(:,cyn,cxn)=current_state%sqi%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qr) then
      current_state%tend_3d_qr_tvad(:,cyn,cxn)=current_state%sqr%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qs) then
      current_state%tend_3d_qs_tvad(:,cyn,cxn)=current_state%sqs%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qg) then
      current_state%tend_3d_qg_tvad(:,cyn,cxn)=current_state%sqg%data(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      current_state%tend_3d_tabs_tvad(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn) * &
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
    if (l_tend_3d_u) then
      current_state%tend_3d_u_tvad(:,cyn,cxn)=current_state%su%data(:,cyn,cxn) - &
      current_state%tend_3d_u_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      current_state%tend_3d_v_tvad(:,cyn,cxn)=current_state%sv%data(:,cyn,cxn) - &
      current_state%tend_3d_v_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_w) then
      current_state%tend_3d_w_tvad(:,cyn,cxn)=current_state%sw%data(:,cyn,cxn) - &
      current_state%tend_3d_w_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_th) then
      current_state%tend_3d_th_tvad(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn) - &
      current_state%tend_3d_th_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      current_state%tend_3d_qv_tvad(:,cyn,cxn)=current_state%sqv%data(:,cyn,cxn) - &
      current_state%tend_3d_qv_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      current_state%tend_3d_ql_tvad(:,cyn,cxn)=current_state%sql%data(:,cyn,cxn) - &
      current_state%tend_3d_ql_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_qi) then
      current_state%tend_3d_qi_tvad(:,cyn,cxn)=current_state%sqi%data(:,cyn,cxn) - &
      current_state%tend_3d_qi_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_qr) then
      current_state%tend_3d_qr_tvad(:,cyn,cxn)=current_state%sqr%data(:,cyn,cxn) - &
      current_state%tend_3d_qr_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_qs) then
      current_state%tend_3d_qs_tvad(:,cyn,cxn)=current_state%sqs%data(:,cyn,cxn) - &
      current_state%tend_3d_qs_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_qg) then
      current_state%tend_3d_qg_tvad(:,cyn,cxn)=current_state%sqg%data(:,cyn,cxn) - &
      current_state%tend_3d_qg_tvad(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      current_state%tend_3d_tabs_tvad(:,cyn,cxn)=   &
       current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)   &
        - current_state%tend_3d_tabs_tvad(:,cyn,cxn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_u) then
      current_state%tend_pr_tot_u_tvad(:)=current_state%tend_pr_tot_u_tvad(:) + &
      current_state%tend_3d_u_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_v) then
      current_state%tend_pr_tot_v_tvad(:)=current_state%tend_pr_tot_v_tvad(:) + &
      current_state%tend_3d_v_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_w) then
      current_state%tend_pr_tot_w_tvad(:)=current_state%tend_pr_tot_w_tvad(:) + &
      current_state%tend_3d_w_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_th) then
      current_state%tend_pr_tot_th_tvad(:)=current_state%tend_pr_tot_th_tvad(:) + &
      current_state%tend_3d_th_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qv) then
      current_state%tend_pr_tot_qv_tvad(:)=current_state%tend_pr_tot_qv_tvad(:) + &
      current_state%tend_3d_qv_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_ql) then
      current_state%tend_pr_tot_ql_tvad(:)=current_state%tend_pr_tot_ql_tvad(:) + &
      current_state%tend_3d_ql_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qi) then
      current_state%tend_pr_tot_qi_tvad(:)=current_state%tend_pr_tot_qi_tvad(:) + &
      current_state%tend_3d_qi_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qr) then
      current_state%tend_pr_tot_qr_tvad(:)=current_state%tend_pr_tot_qr_tvad(:) + &
      current_state%tend_3d_qr_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qs) then
      current_state%tend_pr_tot_qs_tvad(:)=current_state%tend_pr_tot_qs_tvad(:) + &
      current_state%tend_3d_qs_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qg) then
      current_state%tend_pr_tot_qg_tvad(:)=current_state%tend_pr_tot_qg_tvad(:) + &
      current_state%tend_3d_qg_tvad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_tabs) then
      current_state%tend_pr_tot_tabs_tvad(:)=current_state%tend_pr_tot_tabs_tvad(:) + &
      current_state%tend_3d_tabs_tvad(:,cyn,cxn)
    endif

  end subroutine compute_component_tendencies


  !> Parses a field string (read in from the configuration file) and determines whether this algorithm should be used
  !! for advecting that field
  !! @param field The string configuration of field advection
  !! @returns Whether or not the field is advected here
  logical function determine_if_advection_here(field)
    character(len=*), intent(in) :: field

    if (len_trim(field) .ne. 0) then
      if (trim(field) .eq. "tvd" .or. trim(field) .eq. "any") then
        determine_if_advection_here=.true.
      else
        determine_if_advection_here=.false.
      end if
    else
      determine_if_advection_here=.true.
    end if
  end function determine_if_advection_here
end module tvdadvection_mod
