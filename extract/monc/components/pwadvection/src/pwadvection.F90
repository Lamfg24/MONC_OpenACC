!> Piacsek-Williams advection scheme
module pwadvection_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type, &
      component_descriptor_type_v1
  use optionsdatabase_mod, only : options_get_string, options_get_integer
  use collections_mod, only : map_type
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use q_indices_mod, only: get_q_index, standard_q_names

implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: advect_flow, advect_th, advect_q

  logical :: l_toplevel=.true.

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  logical :: l_tend_3d_u, l_tend_3d_v, l_tend_3d_w, l_tend_3d_th,l_tend_3d_qv,       &
             l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,       &
             l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w, l_tend_pr_tot_th,l_tend_pr_tot_qv,       &
             l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,l_tend_pr_tot_qs,l_tend_pr_tot_qg,       &
             l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=1, iql=2, iqr=3, iqi=4, iqs=5, iqg=6

 public initialisation_callback_pwadvection, timestep_callback_pwadvection, &
        finalisation_callback_pwadvection
contains


  !> Initialisation callback, will set up the configuration of this advection scheme
  !! @param current_state The current model state
  subroutine initialisation_callback_pwadvection(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    logical :: l_qdiag
    integer :: iter
    integer :: alloc_z, alloc_y, alloc_x

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "advection_flow_fields") then
        advect_flow=determine_if_advection_here(current_state%options_database_string(iter,2))
      else if (current_state%options_database_string(iter,1) .eq. "advection_theta_field") then
        advect_th=determine_if_advection_here(current_state%options_database_string(iter,2))
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

    l_qdiag =  (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) !.and. advect_q

    l_tend_pr_tot_u   = current_state%u%active !.and. advect_flow
    l_tend_pr_tot_v   = current_state%v%active !.and. advect_flow
    l_tend_pr_tot_w   = current_state%w%active !.and. advect_flow
    l_tend_pr_tot_th  = current_state%th%active! .and. advect_th
    l_tend_pr_tot_qv  = l_qdiag !.and. current_state%water_vapour_mixing_ratio_index > 0
    l_tend_pr_tot_ql  = l_qdiag !.and. current_state%liquid_water_mixing_ratio_index > 0
    l_tend_pr_tot_qi  = l_qdiag !.and. current_state%ice_water_mixing_ratio_index > 0
    l_tend_pr_tot_qr  = l_qdiag !.and. current_state%rain_water_mixing_ratio_index > 0
    l_tend_pr_tot_qs  = l_qdiag !.and. current_state%snow_water_mixing_ratio_index > 0
    l_tend_pr_tot_qg  = l_qdiag !.and. current_state%graupel_water_mixing_ratio_index > 0
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_u   = (current_state%u%active .and. advect_flow) .or. l_tend_pr_tot_u 
    l_tend_3d_v   = (current_state%v%active .and. advect_flow) .or. l_tend_pr_tot_v
    l_tend_3d_w   = (current_state%w%active .and. advect_flow) .or. l_tend_pr_tot_w
    l_tend_3d_th  = (current_state%th%active .and. advect_th)  .or. l_tend_pr_tot_th
    l_tend_3d_qv  = (l_qdiag .and. current_state%water_vapour_mixing_ratio_index > 0) .or. l_tend_pr_tot_qv
    l_tend_3d_ql  = (l_qdiag .and. current_state%liquid_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_ql
    l_tend_3d_qi  = (l_qdiag .and. current_state%ice_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qi
    l_tend_3d_qr  = (l_qdiag .and. current_state%rain_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qr
    l_tend_3d_qs  = (l_qdiag .and. current_state%snow_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qs
    l_tend_3d_qg  = (l_qdiag .and. current_state%graupel_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qg
    l_tend_3d_tabs = l_tend_3d_th

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_u) then
      !allocate( current_state%tend_3d_u_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_u_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_v) then
      !allocate( current_state%tend_3d_v_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_v_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_w) then
      !allocate( current_state%tend_3d_w_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_w_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_th) then
      !allocate( current_state%tend_3d_th_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_th_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qv) then
      !iqv=get_q_index(standard_q_names%VAPOUR, 'pw_advection')
      !allocate( current_state%tend_3d_qv_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qv_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_ql) then
      !iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'pw_advection')
      !allocate( current_state%tend_3d_ql_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_ql_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qi) then
      !iqi=get_q_index(standard_q_names%ICE_MASS, 'pw_advection')
      !allocate( current_state%tend_3d_qi_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qi_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qr) then
      !iqr=get_q_index(standard_q_names%RAIN_MASS, 'pw_advection')
      !allocate( current_state%tend_3d_qr_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qr_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qs) then
      !iqs=get_q_index(standard_q_names%SNOW_MASS, 'pw_advection')
      !allocate( current_state%tend_3d_qs_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qs_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_qg) then
      !iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'pw_advection')
      !allocate( current_state%tend_3d_qg_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qg_pwad(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_tabs) then
      !allocate( current_state%tend_3d_tabs_pwad(current_state%local_grid%size(Z_INDEX),  &
      !                       current_state%local_grid%size(Y_INDEX),  &
      !                       current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_tabs_pwad(alloc_z, alloc_y, alloc_x))
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_u) then
      !allocate( current_state%tend_pr_tot_u_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_u_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_v) then
      !allocate( current_state%tend_pr_tot_v_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_v_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_w) then
      !allocate( current_state%tend_pr_tot_w_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_w_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_th) then
      !allocate( current_state%tend_pr_tot_th_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_th_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_qv) then
      !allocate( current_state%tend_pr_tot_qv_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qv_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_ql) then
      !allocate( current_state%tend_pr_tot_ql_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_ql_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_qi) then
      !allocate( current_state%tend_pr_tot_qi_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qi_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_qr) then
      !allocate( current_state%tend_pr_tot_qr_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qr_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_qs) then
      !allocate( current_state%tend_pr_tot_qs_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qs_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_qg) then
      !allocate( current_state%tend_pr_tot_qg_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qg_pwad(alloc_z))
    endif
    if (l_tend_pr_tot_tabs) then
      !allocate( current_state%tend_pr_tot_tabs_pwad(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_tabs_pwad(alloc_z))
    endif

  end subroutine initialisation_callback_pwadvection


  subroutine finalisation_callback_pwadvection(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(current_state%tend_3d_u_pwad)) deallocate(current_state%tend_3d_u_pwad)
    if (allocated(current_state%tend_3d_v_pwad)) deallocate(current_state%tend_3d_v_pwad)
    if (allocated(current_state%tend_3d_w_pwad)) deallocate(current_state%tend_3d_w_pwad)
    if (allocated(current_state%tend_3d_th_pwad)) deallocate(current_state%tend_3d_th_pwad)
    if (allocated(current_state%tend_3d_qv_pwad)) deallocate(current_state%tend_3d_qv_pwad)
    if (allocated(current_state%tend_3d_ql_pwad)) deallocate(current_state%tend_3d_ql_pwad)
    if (allocated(current_state%tend_3d_qi_pwad)) deallocate(current_state%tend_3d_qi_pwad)
    if (allocated(current_state%tend_3d_qr_pwad)) deallocate(current_state%tend_3d_qr_pwad)
    if (allocated(current_state%tend_3d_qs_pwad)) deallocate(current_state%tend_3d_qs_pwad)
    if (allocated(current_state%tend_3d_qg_pwad)) deallocate(current_state%tend_3d_qg_pwad)
    if (allocated(current_state%tend_3d_tabs_pwad)) deallocate(current_state%tend_3d_tabs_pwad)

    if (allocated(current_state%tend_pr_tot_u_pwad)) deallocate(current_state%tend_pr_tot_u_pwad)
    if (allocated(current_state%tend_pr_tot_v_pwad)) deallocate(current_state%tend_pr_tot_v_pwad)
    if (allocated(current_state%tend_pr_tot_w_pwad)) deallocate(current_state%tend_pr_tot_w_pwad)
    if (allocated(current_state%tend_pr_tot_th_pwad)) deallocate(current_state%tend_pr_tot_th_pwad)
    if (allocated(current_state%tend_pr_tot_qv_pwad)) deallocate(current_state%tend_pr_tot_qv_pwad)
    if (allocated(current_state%tend_pr_tot_ql_pwad)) deallocate(current_state%tend_pr_tot_ql_pwad)
    if (allocated(current_state%tend_pr_tot_qi_pwad)) deallocate(current_state%tend_pr_tot_qi_pwad)
    if (allocated(current_state%tend_pr_tot_qr_pwad)) deallocate(current_state%tend_pr_tot_qr_pwad)
    if (allocated(current_state%tend_pr_tot_qs_pwad)) deallocate(current_state%tend_pr_tot_qs_pwad)
    if (allocated(current_state%tend_pr_tot_qg_pwad)) deallocate(current_state%tend_pr_tot_qg_pwad)
    if (allocated(current_state%tend_pr_tot_tabs_pwad)) deallocate(current_state%tend_pr_tot_tabs_pwad)

  end subroutine finalisation_callback_pwadvection


  !> Called per column of data, this will perform Piacsek-Williams advection on the applicable fields for non halo data
  !! @param current_state The current model state
  !! @param target_(x/y)_index This is the index with the halos subtracted. This is needed so that diagnostic does 
  !!                           not include halos and to prevent array out-of-bounds
  subroutine timestep_callback_pwadvection(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: current_x_index, current_y_index, target_x_index, target_y_index
    logical :: calculate_diagnostics

    !calculate_diagnostics = current_state%diagnostic_sample_timestep

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)


    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_u) then
        current_state%tend_pr_tot_u_pwad(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_v) then
        current_state%tend_pr_tot_v_pwad(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_w) then
        current_state%tend_pr_tot_w_pwad(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_th) then
        current_state%tend_pr_tot_th_pwad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qv) then
        current_state%tend_pr_tot_qv_pwad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_ql) then
        current_state%tend_pr_tot_ql_pwad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qi) then
        current_state%tend_pr_tot_qi_pwad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qr) then
        current_state%tend_pr_tot_qr_pwad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qs) then
        current_state%tend_pr_tot_qs_pwad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qg) then
        current_state%tend_pr_tot_qg_pwad(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_tabs) then
        current_state%tend_pr_tot_tabs_pwad(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals


    if (current_state%halo_column) return
    
    if (current_state%modulo_number_3d .eq. 0) &
        call save_precomponent_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)
    if (advect_flow) call advect_flow_fields(current_state, current_x_index, current_y_index)
    if (advect_th) call advect_th_field(current_state, current_x_index, current_y_index)
    if (advect_q) call advect_q_field(current_state, current_x_index, current_y_index)

    if (current_state%modulo_number_3d .eq. 0) &
        call compute_component_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)

  end subroutine timestep_callback_pwadvection

  !> Advects the q fields in a column
  !! @param current_state The current model state
  !! @param current_x_index The current slice, x, index
  !! @param current_y_index The current column, y, index
  subroutine advect_q_field(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: k, n

    !do n=1,current_state%number_q_fields
    do k=2,current_state%local_grid%size(Z_INDEX)-1
#ifdef U_ACTIVE
      current_state%sqv%data(k, current_y_index, current_x_index)=&
            current_state%sqv%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qv%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qv%data(k, current_y_index, current_x_index+1))
      current_state%sql%data(k, current_y_index, current_x_index)=&
            current_state%sql%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%ql%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%ql%data(k, current_y_index, current_x_index+1))
      current_state%sqr%data(k, current_y_index, current_x_index)=&
            current_state%sqr%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qr%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qr%data(k, current_y_index, current_x_index+1))
      current_state%sqi%data(k, current_y_index, current_x_index)=&
            current_state%sqi%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qi%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qi%data(k, current_y_index, current_x_index+1))
      current_state%sqs%data(k, current_y_index, current_x_index)=&
            current_state%sqs%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qs%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qs%data(k, current_y_index, current_x_index+1))
      current_state%sqg%data(k, current_y_index, current_x_index)=&
            current_state%sqg%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qg%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qg%data(k, current_y_index, current_x_index+1))
      current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)+ &
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qAitkenSolMass%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qAitkenSolMass%data(k, current_y_index, current_x_index+1))
      current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)+ &
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qAccumSolMass%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qAccumSolMass%data(k, current_y_index, current_x_index+1))
      current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)+ &
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qAccumInsolMass%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qAccumInsolMass%data(k, current_y_index, current_x_index+1))
      current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)+ &
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qCoarseSolMass%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qCoarseSolMass%data(k, current_y_index, current_x_index+1))
      current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)=&
            current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%qCoarseDustMass%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%qCoarseDustMass%data(k, current_y_index, current_x_index+1))

      current_state%snl%data(k, current_y_index, current_x_index)=&
            current_state%snl%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%nl%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%nl%data(k, current_y_index, current_x_index+1))
      current_state%snr%data(k, current_y_index, current_x_index)=&
            current_state%snr%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%nr%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%nr%data(k, current_y_index, current_x_index+1))
      current_state%sni%data(k, current_y_index, current_x_index)=&
            current_state%sni%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%ni%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%ni%data(k, current_y_index, current_x_index+1))
      current_state%sns%data(k, current_y_index, current_x_index)=&
            current_state%sns%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%ns%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%ns%data(k, current_y_index, current_x_index+1))
      current_state%sng%data(k, current_y_index, current_x_index)=&
            current_state%sng%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%ng%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%ng%data(k, current_y_index, current_x_index+1))
      current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%nAitkenSolNumber%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%nAitkenSolNumber%data(k, current_y_index, current_x_index+1))
      current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%nAccumSolNumber%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%nAccumSolNumber%data(k, current_y_index, current_x_index+1))
      current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%nAccumInsolNumber%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%nAccumInsolNumber%data(k, current_y_index, current_x_index+1))
      current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%nCoarseSolNumber%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%nCoarseSolNumber%data(k, current_y_index, current_x_index+1))
      current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)=&
            current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cx*&
            0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
            current_state%nCoarseDustnumber%data(k, current_y_index, current_x_index-1)-&
            current_state%u%data(k, current_y_index, current_x_index)*&
            current_state%nCoarseDustnumber%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
      current_state%sqv%data(k, current_y_index, current_x_index)=&
            current_state%sqv%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qv%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qv%data(k, current_y_index+1, current_x_index))
      current_state%sql%data(k, current_y_index, current_x_index)=&
            current_state%sql%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%ql%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%ql%data(k, current_y_index+1, current_x_index))
      current_state%sqr%data(k, current_y_index, current_x_index)=&
            current_state%sqr%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qr%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qr%data(k, current_y_index+1, current_x_index))
      current_state%sqi%data(k, current_y_index, current_x_index)=&
            current_state%sqi%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qi%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qi%data(k, current_y_index+1, current_x_index))
      current_state%sqs%data(k, current_y_index, current_x_index)=&
            current_state%sqs%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qs%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qs%data(k, current_y_index+1, current_x_index))
      current_state%sqg%data(k, current_y_index, current_x_index)=&
            current_state%sqg%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qg%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qg%data(k, current_y_index+1, current_x_index))
      current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qAitkenSolMass%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qAitkenSolMass%data(k, current_y_index+1, current_x_index))
      current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qAccumSolMass%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qAccumSolMass%data(k, current_y_index+1, current_x_index))
      current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qAccumInsolMass%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qAccumInsolMass%data(k, current_y_index+1, current_x_index))
      current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qCoarseSolMass%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qCoarseSolMass%data(k, current_y_index+1, current_x_index))
      current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)=&
            current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%qCoarseDustMass%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%qCoarseDustMass%data(k, current_y_index+1, current_x_index))


      current_state%snl%data(k, current_y_index, current_x_index)=&
            current_state%snl%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%nl%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%nl%data(k, current_y_index+1, current_x_index))
      current_state%snr%data(k, current_y_index, current_x_index)=&
            current_state%snr%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%nr%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%nr%data(k, current_y_index+1, current_x_index))
      current_state%sni%data(k, current_y_index, current_x_index)=&
            current_state%sni%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%ni%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%ni%data(k, current_y_index+1, current_x_index))
      current_state%sns%data(k, current_y_index, current_x_index)=&
            current_state%sns%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%ns%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%ns%data(k, current_y_index+1, current_x_index))
      current_state%sng%data(k, current_y_index, current_x_index)=&
            current_state%sng%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%ng%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%ng%data(k, current_y_index+1, current_x_index))
      current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%nAitkenSolNumber%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%nAitkenSolNumber%data(k, current_y_index+1, current_x_index))
      current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%nAccumSolNumber%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%nAccumSolNumber%data(k, current_y_index+1, current_x_index))
      current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%nAccumInsolNumber%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%nAccumInsolNumber%data(k, current_y_index+1, current_x_index))
      current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%nCoarseSolNumber%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%nCoarseSolNumber%data(k, current_y_index+1, current_x_index))
      current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)=&
            current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)+&
            current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
            (current_state%v%data(k, current_y_index-1, current_x_index)*&
            current_state%nCoarseDustnumber%data(k, current_y_index-1, current_x_index)-&
            current_state%v%data(k, current_y_index, current_x_index)*&
            current_state%nCoarseDustnumber%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
      current_state%sqv%data(k, current_y_index, current_x_index)=&
            current_state%sqv%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qv%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qv%data(k+1, current_y_index, current_x_index))
      current_state%sql%data(k, current_y_index, current_x_index)=&
            current_state%sql%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%ql%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%ql%data(k+1, current_y_index, current_x_index))
      current_state%sqr%data(k, current_y_index, current_x_index)=&
            current_state%sqr%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qr%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qr%data(k+1, current_y_index, current_x_index))
      current_state%sqi%data(k, current_y_index, current_x_index)=&
            current_state%sqi%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qi%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qi%data(k+1, current_y_index, current_x_index))
      current_state%sqs%data(k, current_y_index, current_x_index)=&
            current_state%sqs%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qs%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qs%data(k+1, current_y_index, current_x_index))
      current_state%sqg%data(k, current_y_index, current_x_index)=&
            current_state%sqg%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qg%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qg%data(k+1, current_y_index, current_x_index))
      current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qAitkenSolMass%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qAitkenSolMass%data(k+1, current_y_index, current_x_index))
      current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qAccumSolMass%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qAccumSolMass%data(k+1, current_y_index, current_x_index))
      current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qAccumInsolMass%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qAccumInsolMass%data(k+1, current_y_index, current_x_index))
      current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)=&
            current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qCoarseSolMass%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qCoarseSolMass%data(k+1, current_y_index, current_x_index))
      current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)=&
            current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%qCoarseDustMass%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%qCoarseDustMass%data(k+1, current_y_index, current_x_index))


      current_state%snl%data(k, current_y_index, current_x_index)=&
            current_state%snl%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%nl%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%nl%data(k+1, current_y_index, current_x_index))
      current_state%snr%data(k, current_y_index, current_x_index)=&
            current_state%snr%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%nr%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%nr%data(k+1, current_y_index, current_x_index))
      current_state%sni%data(k, current_y_index, current_x_index)=&
            current_state%sni%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%ni%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%ni%data(k+1, current_y_index, current_x_index))
      current_state%sns%data(k, current_y_index, current_x_index)=&
            current_state%sns%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%ns%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%ns%data(k+1, current_y_index, current_x_index))
      current_state%sng%data(k, current_y_index, current_x_index)=&
            current_state%sng%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%ng%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%ng%data(k+1, current_y_index, current_x_index))
      current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%nAitkenSolNumber%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%nAitkenSolNumber%data(k+1, current_y_index, current_x_index))
      current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%nAccumSolNumber%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%nAccumSolNumber%data(k+1, current_y_index, current_x_index))
      current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%nAccumInsolNumber%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%nAccumInsolNumber%data(k+1, current_y_index, current_x_index))
      current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)=&
            current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%nCoarseSolNumber%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%nCoarseSolNumber%data(k+1, current_y_index, current_x_index))
      current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)=&
            current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)+&
            2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
            current_state%w%data(k-1, current_y_index, current_x_index)*&
            current_state%nCoarseDustnumber%data(k-1, current_y_index, current_x_index)-&
            current_state%global_grid%configuration%vertical%tzc2(k)*&
            current_state%w%data(k, current_y_index, current_x_index)*&
            current_state%nCoarseDustnumber%data(k+1, current_y_index, current_x_index))
#endif
    end do

    if (l_toplevel)then
    k=current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
    current_state%sqv%data(k, current_y_index, current_x_index)=&
          current_state%sqv%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qv%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qv%data(k, current_y_index, current_x_index+1))
    current_state%sql%data(k, current_y_index, current_x_index)=&
          current_state%sql%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%ql%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%ql%data(k, current_y_index, current_x_index+1))
    current_state%sqr%data(k, current_y_index, current_x_index)=&
          current_state%sqr%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qr%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qr%data(k, current_y_index, current_x_index+1))
    current_state%sqi%data(k, current_y_index, current_x_index)=&
          current_state%sqi%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qi%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qi%data(k, current_y_index, current_x_index+1))
    current_state%sqs%data(k, current_y_index, current_x_index)=&
          current_state%sqs%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qs%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qs%data(k, current_y_index, current_x_index+1))
    current_state%sqg%data(k, current_y_index, current_x_index)=&
          current_state%sqg%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qg%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qg%data(k, current_y_index, current_x_index+1))
    current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qAitkenSolMass%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qAitkenSolMass%data(k, current_y_index, current_x_index+1))
    current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qAccumSolMass%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qAccumSolMass%data(k, current_y_index, current_x_index+1))
    current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qAccumInsolMass%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qAccumInsolMass%data(k, current_y_index, current_x_index+1))
    current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qCoarseSolMass%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qCoarseSolMass%data(k, current_y_index, current_x_index+1))
    current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)=&
          current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%qCoarseDustMass%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%qCoarseDustMass%data(k, current_y_index, current_x_index+1))


    current_state%snl%data(k, current_y_index, current_x_index)=&
          current_state%snl%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%nl%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%nl%data(k, current_y_index, current_x_index+1))
    current_state%snr%data(k, current_y_index, current_x_index)=&
          current_state%snr%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%nr%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%nr%data(k, current_y_index, current_x_index+1))
    current_state%sni%data(k, current_y_index, current_x_index)=&
          current_state%sni%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%ni%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%ni%data(k, current_y_index, current_x_index+1))
    current_state%sns%data(k, current_y_index, current_x_index)=&
          current_state%sns%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%ns%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%ns%data(k, current_y_index, current_x_index+1))
    current_state%sng%data(k, current_y_index, current_x_index)=&
          current_state%sng%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%ng%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%ng%data(k, current_y_index, current_x_index+1))
    current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%nAitkenSolNumber%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%nAitkenSolNumber%data(k, current_y_index, current_x_index+1))
    current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%nAccumSolNumber%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%nAccumSolNumber%data(k, current_y_index, current_x_index+1))
    current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%nAccumInsolNumber%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%nAccumInsolNumber%data(k, current_y_index, current_x_index+1))
    current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%nCoarseSolNumber%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%nCoarseSolNumber%data(k, current_y_index, current_x_index+1))
    current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)=&
          current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%horizontal%cx*&
          0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
          current_state%nCoarseDustnumber%data(k, current_y_index, current_x_index-1)-&
          current_state%u%data(k, current_y_index, current_x_index)*&
          current_state%nCoarseDustnumber%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
    current_state%sqv%data(k, current_y_index, current_x_index)=&
          current_state%sqv%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qv%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qv%data(k, current_y_index+1, current_x_index))
    current_state%sql%data(k, current_y_index, current_x_index)=&
          current_state%sql%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%ql%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%ql%data(k, current_y_index+1, current_x_index))
    current_state%sqr%data(k, current_y_index, current_x_index)=&
          current_state%sqr%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qr%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qr%data(k, current_y_index+1, current_x_index))
    current_state%sqi%data(k, current_y_index, current_x_index)=&
          current_state%sqi%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qi%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qi%data(k, current_y_index+1, current_x_index))
    current_state%sqs%data(k, current_y_index, current_x_index)=&
          current_state%sqs%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qs%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qs%data(k, current_y_index+1, current_x_index))
    current_state%sqg%data(k, current_y_index, current_x_index)=&
          current_state%sqg%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qg%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qg%data(k, current_y_index+1, current_x_index))
    current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qAitkenSolMass%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qAitkenSolMass%data(k, current_y_index+1, current_x_index))
    current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qAccumSolMass%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qAccumSolMass%data(k, current_y_index+1, current_x_index))
    current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qAccumInsolMass%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qAccumInsolMass%data(k, current_y_index+1, current_x_index))
    current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qCoarseSolMass%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qCoarseSolMass%data(k, current_y_index+1, current_x_index))
    current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)=&
          current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%qCoarseDustMass%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%qCoarseDustMass%data(k, current_y_index+1, current_x_index))


    current_state%snl%data(k, current_y_index, current_x_index)=&
          current_state%snl%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%nl%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%nl%data(k, current_y_index+1, current_x_index))
    current_state%snr%data(k, current_y_index, current_x_index)=&
          current_state%snr%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%nr%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%nr%data(k, current_y_index+1, current_x_index))
    current_state%sni%data(k, current_y_index, current_x_index)=&
          current_state%sni%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%ni%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%ni%data(k, current_y_index+1, current_x_index))
    current_state%sns%data(k, current_y_index, current_x_index)=&
          current_state%sns%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%ns%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%ns%data(k, current_y_index+1, current_x_index))
    current_state%sng%data(k, current_y_index, current_x_index)=&
          current_state%sng%data(k, current_y_index, current_x_index)+current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%ng%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%ng%data(k, current_y_index+1, current_x_index))
    current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%nAitkenSolNumber%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%nAitkenSolNumber%data(k, current_y_index+1, current_x_index))
    current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%nAccumSolNumber%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%nAccumSolNumber%data(k, current_y_index+1, current_x_index))
    current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%nAccumInsolNumber%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%nAccumInsolNumber%data(k, current_y_index+1, current_x_index))
    current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%nCoarseSolNumber%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%nCoarseSolNumber%data(k, current_y_index+1, current_x_index))
    current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)=&
          current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)+ &
          current_state%global_grid%configuration%horizontal%cy*&
          0.5_DEFAULT_PRECISION*(current_state%v%data(k, current_y_index-1, current_x_index)*&
          current_state%nCoarseDustnumber%data(k, current_y_index-1, current_x_index)-&
          current_state%v%data(k, current_y_index, current_x_index)*&
          current_state%nCoarseDustnumber%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
    current_state%sqv%data(k, current_y_index, current_x_index)=&
          current_state%sqv%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qv%data(k-1, current_y_index, current_x_index)
    current_state%sql%data(k, current_y_index, current_x_index)=&
          current_state%sql%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%ql%data(k-1, current_y_index, current_x_index)
    current_state%sqr%data(k, current_y_index, current_x_index)=&
          current_state%sqr%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qr%data(k-1, current_y_index, current_x_index)
    current_state%sqi%data(k, current_y_index, current_x_index)=&
          current_state%sqi%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qi%data(k-1, current_y_index, current_x_index)
    current_state%sqs%data(k, current_y_index, current_x_index)=&
          current_state%sqs%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qs%data(k-1, current_y_index, current_x_index)
    current_state%sqg%data(k, current_y_index, current_x_index)=&
          current_state%sqg%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qg%data(k-1, current_y_index, current_x_index)
    current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqAitkenSolMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qAitkenSolMass%data(k-1, current_y_index, current_x_index)
    current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqAccumSolMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qAccumSolMass%data(k-1, current_y_index, current_x_index)
    current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqAccumInsolMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qAccumInsolMass%data(k-1, current_y_index, current_x_index)
    current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)=&
          current_state%sqCoarseSolMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qCoarseSolMass%data(k-1, current_y_index, current_x_index)
    current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)=&
          current_state%sqCoarseDustMass%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%qCoarseDustMass%data(k-1, current_y_index, current_x_index)


    current_state%snl%data(k, current_y_index, current_x_index)=&
          current_state%snl%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%nl%data(k-1, current_y_index, current_x_index)
    current_state%snr%data(k, current_y_index, current_x_index)=&
          current_state%snr%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%nr%data(k-1, current_y_index, current_x_index)
    current_state%sni%data(k, current_y_index, current_x_index)=&
          current_state%sni%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%ni%data(k-1, current_y_index, current_x_index)
    current_state%sns%data(k, current_y_index, current_x_index)=&
          current_state%sns%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%ns%data(k-1, current_y_index, current_x_index)
    current_state%sng%data(k, current_y_index, current_x_index)=&
          current_state%sng%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%ng%data(k-1, current_y_index, current_x_index)
    current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snAitkenSolNumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%nAitkenSolNumber%data(k-1, current_y_index, current_x_index)
    current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snAccumSolNumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%nAccumSolNumber%data(k-1, current_y_index, current_x_index)
    current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snAccumInsolNumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%nAccumInsolNumber%data(k-1, current_y_index, current_x_index)
    current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)=&
          current_state%snCoarseSolNumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%nCoarseSolNumber%data(k-1, current_y_index, current_x_index)
    current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)=&
          current_state%snCoarseDustnumber%data(k, current_y_index, current_x_index)+&
          current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
          current_state%w%data(k-1, current_y_index, current_x_index)*&
          current_state%nCoarseDustnumber%data(k-1, current_y_index, current_x_index)
#endif
    end if
    !end do
  end subroutine advect_q_field

  !> Advects the theta field in a column
  !! @param current_state The current model state
  !! @param current_x_index The current slice, x, index
  !! @param current_y_index The current column, y, index
  subroutine advect_th_field(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: k

    if (current_state%th%active) then
      do k=2,current_state%local_grid%size(Z_INDEX)-1
#ifdef U_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)= current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cx*&
             0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
             current_state%th%data(k, current_y_index, current_x_index-1)-&
             current_state%u%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
             (current_state%v%data(k, current_y_index-1, current_x_index)*&
             current_state%th%data(k, current_y_index-1, current_x_index)-&
             current_state%v%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             2.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc1(k)*&
             current_state%w%data(k-1, current_y_index, current_x_index)*&
             current_state%th%data(k-1, current_y_index, current_x_index)-&
             current_state%global_grid%configuration%vertical%tzc2(k)*&
             current_state%w%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k+1, current_y_index, current_x_index))
#endif
      end do

      if (l_toplevel)then

        k=current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)= current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cx*&
             0.5_DEFAULT_PRECISION*(current_state%u%data(k, current_y_index, current_x_index-1)*&
             current_state%th%data(k, current_y_index, current_x_index-1)-&
             current_state%u%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index, current_x_index+1))
#endif
#ifdef V_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%horizontal%cy*0.5_DEFAULT_PRECISION*&
             (current_state%v%data(k, current_y_index-1, current_x_index)*&
             current_state%th%data(k, current_y_index-1, current_x_index)-&
             current_state%v%data(k, current_y_index, current_x_index)*&
             current_state%th%data(k, current_y_index+1, current_x_index))
#endif
#ifdef W_ACTIVE
        current_state%sth%data(k, current_y_index, current_x_index)=current_state%sth%data(k, current_y_index, current_x_index)+&
             current_state%global_grid%configuration%vertical%tzc1(k)*2.0_DEFAULT_PRECISION*&
             current_state%w%data(k-1, current_y_index, current_x_index)*current_state%th%data(k-1, current_y_index, &
             current_x_index)
#endif
      end if
   end if
  end subroutine advect_th_field

  !> Advects the flow fields depending upon which fields are active in the model in a column
  !! @param current_state The current model state
  !! @param current_x_index The current slice, x, index
  !! @param current_y_index The current column, y, index
  subroutine advect_flow_fields(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: k

    do k=2,current_state%local_grid%size(Z_INDEX)-1
#ifdef U_ACTIVE
      current_state%su%data(k, current_y_index, current_x_index)=&
           current_state%global_grid%configuration%horizontal%tcx*(current_state%u%data(k, current_y_index, current_x_index-1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k, current_y_index, current_x_index-1))-&
           current_state%u%data(k, current_y_index, current_x_index+1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k, current_y_index, current_x_index+1)))
#ifdef V_ACTIVE
      current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcy*(current_state%u%data(k, current_y_index-1, current_x_index)*&
           (current_state%v%data(k, current_y_index-1, current_x_index)+&
           current_state%v%data(k, current_y_index-1, current_x_index+1))-&
           current_state%u%data(k, current_y_index+1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k, current_y_index, current_x_index+1)))
#endif
#ifdef W_ACTIVE
      current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
           (current_state%global_grid%configuration%vertical%tzc1(k)*current_state%u%data(k-1, current_y_index, current_x_index)*&
           (current_state%w%data(k-1, current_y_index, current_x_index)+&
           current_state%w%data(k-1, current_y_index, current_x_index+1))-&
           current_state%global_grid%configuration%vertical%tzc2(k)*current_state%u%data(k+1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k, current_y_index, current_x_index+1)))
#endif
#endif

#ifdef V_ACTIVE
      current_state%sv%data(k, current_y_index, current_x_index)=current_state%global_grid%configuration%horizontal%tcy*(&
           current_state%v%data(k, current_y_index-1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k, current_y_index-1, current_x_index))-&
           current_state%v%data(k, current_y_index+1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k, current_y_index+1, current_x_index)))
#ifdef U_ACTIVE
      current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcx*(current_state%v%data(k, current_y_index, current_x_index-1)*&
           (current_state%u%data(k, current_y_index, current_x_index-1)+&
           current_state%u%data(k, current_y_index+1, current_x_index-1))-&
           current_state%v%data(k, current_y_index, current_x_index+1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k, current_y_index+1, current_x_index)))
#endif
#ifdef W_ACTIVE
      current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
           (current_state%global_grid%configuration%vertical%tzc1(k)*current_state%v%data(k-1, current_y_index, current_x_index)*&
           (current_state%w%data(k-1, current_y_index, current_x_index)+&
           current_state%w%data(k-1, current_y_index+1, current_x_index))-&
           current_state%global_grid%configuration%vertical%tzc2(k)*current_state%v%data(k+1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k, current_y_index+1, current_x_index)))
#endif
#endif

#ifdef W_ACTIVE
      current_state%sw%data(k, current_y_index, current_x_index)=(current_state%global_grid%configuration%vertical%tzd1(k)*&
           current_state%w%data(k-1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k-1, current_y_index, current_x_index))-&
           current_state%global_grid%configuration%vertical%tzd2(k)*current_state%w%data(k+1, current_y_index, current_x_index)*&
           (current_state%w%data(k, current_y_index, current_x_index)+&
           current_state%w%data(k+1, current_y_index, current_x_index)))
#ifdef U_ACTIVE
      current_state%sw%data(k, current_y_index, current_x_index)=current_state%sw%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcx*(current_state%w%data(k, current_y_index, current_x_index-1)*&
           (current_state%u%data(k, current_y_index, current_x_index-1)+&
           current_state%u%data(k+1, current_y_index, current_x_index-1))-&
           current_state%w%data(k, current_y_index, current_x_index+1)*&
           (current_state%u%data(k, current_y_index, current_x_index)+&
           current_state%u%data(k+1, current_y_index, current_x_index)))
#endif
#ifdef V_ACTIVE
      current_state%sw%data(k, current_y_index, current_x_index)=current_state%sw%data(k, current_y_index, current_x_index)+&
           current_state%global_grid%configuration%horizontal%tcy*(current_state%w%data(k, current_y_index-1, current_x_index)*&
           (current_state%v%data(k, current_y_index-1, current_x_index)+&
           current_state%v%data(k+1, current_y_index-1, current_x_index))-&
           current_state%w%data(k, current_y_index+1, current_x_index)*&
           (current_state%v%data(k, current_y_index, current_x_index)+&
           current_state%v%data(k+1, current_y_index, current_x_index)))
#endif
#endif
    end do

    if (l_toplevel)then
    k=current_state%local_grid%size(Z_INDEX)
#ifdef U_ACTIVE
    current_state%su%data(k, current_y_index, current_x_index)=current_state%global_grid%configuration%horizontal%tcx*&
         (current_state%u%data(k, current_y_index, current_x_index-1)*&
         (current_state%u%data(k, current_y_index, current_x_index)+&
         current_state%u%data(k, current_y_index, current_x_index-1))-&
         current_state%u%data(k, current_y_index, current_x_index+1)*&
         (current_state%u%data(k, current_y_index, current_x_index)+&
         current_state%u%data(k, current_y_index, current_x_index+1)))
#ifdef V_ACTIVE
    current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%horizontal%tcy*(current_state%u%data(k, current_y_index-1, current_x_index)*&
         (current_state%v%data(k, current_y_index-1, current_x_index)+&
         current_state%v%data(k, current_y_index-1, current_x_index+1))-&
         current_state%u%data(k, current_y_index+1, current_x_index)*&
         (current_state%v%data(k, current_y_index, current_x_index)+&
         current_state%v%data(k, current_y_index, current_x_index+1)))
#endif
#ifdef W_ACTIVE
    current_state%su%data(k, current_y_index, current_x_index)=current_state%su%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%vertical%tzc1(k)*current_state%u%data(k-1, current_y_index, current_x_index)*&
         (current_state%w%data(k-1, current_y_index, current_x_index)+&
         current_state%w%data(k-1, current_y_index, current_x_index+1))
#endif
#endif

#ifdef V_ACTIVE
    current_state%sv%data(k, current_y_index, current_x_index)=current_state%global_grid%configuration%horizontal%tcy*&
         (current_state%v%data(k, current_y_index-1, current_x_index)*&
         (current_state%v%data(k, current_y_index, current_x_index)+&
         current_state%v%data(k, current_y_index-1, current_x_index))-&
         current_state%v%data(k, current_y_index+1, current_x_index)*&
         (current_state%v%data(k, current_y_index, current_x_index)+&
         current_state%v%data(k, current_y_index+1, current_x_index)))
#ifdef U_ACTIVE
    current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%horizontal%tcx*(current_state%v%data(k, current_y_index, current_x_index-1)*&
         (current_state%u%data(k, current_y_index, current_x_index-1)+&
         current_state%u%data(k, current_y_index+1, current_x_index-1))-&
         current_state%v%data(k, current_y_index, current_x_index+1)*&
         (current_state%u%data(k, current_y_index, current_x_index)+&
         current_state%u%data(k, current_y_index+1, current_x_index)))
#endif
#ifdef W_ACTIVE
    current_state%sv%data(k, current_y_index, current_x_index)=current_state%sv%data(k, current_y_index, current_x_index)+&
         current_state%global_grid%configuration%vertical%tzc1(k)*current_state%v%data(k-1, current_y_index, current_x_index)*&
         (current_state%w%data(k-1, current_y_index, current_x_index)+&
         current_state%w%data(k-1, current_y_index+1, current_x_index))
#endif
#endif
 end if
  end subroutine advect_flow_fields
  

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
      current_state%tend_3d_u_pwad(:,cyn,cxn)=current_state%su%data(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      current_state%tend_3d_v_pwad(:,cyn,cxn)=current_state%sv%data(:,cyn,cxn)
    endif
    if (l_tend_3d_w) then
      current_state%tend_3d_w_pwad(:,cyn,cxn)=current_state%sw%data(:,cyn,cxn)
    endif
    if (l_tend_3d_th) then
      current_state%tend_3d_th_pwad(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      current_state%tend_3d_qv_pwad(:,cyn,cxn)=current_state%sqv%data(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      current_state%tend_3d_ql_pwad(:,cyn,cxn)=current_state%sql%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qi) then
      current_state%tend_3d_qi_pwad(:,cyn,cxn)=current_state%sqi%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qr) then
      current_state%tend_3d_qr_pwad(:,cyn,cxn)=current_state%sqr%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qs) then
      current_state%tend_3d_qs_pwad(:,cyn,cxn)=current_state%sqs%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qg) then
      current_state%tend_3d_qg_pwad(:,cyn,cxn)=current_state%sqg%data(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      current_state%tend_3d_tabs_pwad(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn)* &
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
      current_state%tend_3d_u_pwad(:,cyn,cxn)=current_state%su%data(:,cyn,cxn) - &
      current_state%tend_3d_u_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      current_state%tend_3d_v_pwad(:,cyn,cxn)=current_state%sv%data(:,cyn,cxn) - &
      current_state%tend_3d_v_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_w) then
      current_state%tend_3d_w_pwad(:,cyn,cxn)=current_state%sw%data(:,cyn,cxn)  - &
      current_state%tend_3d_w_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_th) then
      current_state%tend_3d_th_pwad(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn) - &
      current_state%tend_3d_th_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      current_state%tend_3d_qv_pwad(:,cyn,cxn)=current_state%sqv%data(:,cyn,cxn) - &
      current_state%tend_3d_qv_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      current_state%tend_3d_ql_pwad(:,cyn,cxn)=current_state%sql%data(:,cyn,cxn) - &
      current_state%tend_3d_ql_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_qi) then
      current_state%tend_3d_qi_pwad(:,cyn,cxn)=current_state%sqi%data(:,cyn,cxn) - &
      current_state%tend_3d_qi_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_qr) then
      current_state%tend_3d_qr_pwad(:,cyn,cxn)=current_state%sqr%data(:,cyn,cxn) - &
      current_state%tend_3d_qr_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_qs) then
      current_state%tend_3d_qs_pwad(:,cyn,cxn)=current_state%sqs%data(:,cyn,cxn) - &
      current_state%tend_3d_qs_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_qg) then
      current_state%tend_3d_qg_pwad(:,cyn,cxn)=current_state%sqg%data(:,cyn,cxn) - &
      current_state%tend_3d_qg_pwad(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      current_state%tend_3d_tabs_pwad(:,cyn,cxn)=   &
       current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)   &
        - current_state%tend_3d_tabs_pwad(:,cyn,cxn)
    endif

   ! Add local tendency fields to the profile total 
    if (l_tend_pr_tot_u) then
      current_state%tend_pr_tot_u_pwad(:)=current_state%tend_pr_tot_u_pwad(:) + &
      current_state%tend_3d_u_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_v) then
      current_state%tend_pr_tot_v_pwad(:)=current_state%tend_pr_tot_v_pwad(:) + &
      current_state%tend_3d_v_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_w) then
      current_state%tend_pr_tot_w_pwad(:)=current_state%tend_pr_tot_w_pwad(:) + &
      current_state%tend_3d_w_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_th) then
      current_state%tend_pr_tot_th_pwad(:)=current_state%tend_pr_tot_th_pwad(:) + &
      current_state%tend_3d_th_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qv) then
      current_state%tend_pr_tot_qv_pwad(:)=current_state%tend_pr_tot_qv_pwad(:) + &
      current_state%tend_3d_qv_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_ql) then
      current_state%tend_pr_tot_ql_pwad(:)=current_state%tend_pr_tot_ql_pwad(:) + &
      current_state%tend_3d_ql_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qi) then
      current_state%tend_pr_tot_qi_pwad(:)=current_state%tend_pr_tot_qi_pwad(:) + &
      current_state%tend_3d_qi_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qr) then
      current_state%tend_pr_tot_qr_pwad(:)=current_state%tend_pr_tot_qr_pwad(:) + &
      current_state%tend_3d_qr_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qs) then
      current_state%tend_pr_tot_qs_pwad(:)=current_state%tend_pr_tot_qs_pwad(:) + &
      current_state%tend_3d_qs_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qg) then
      current_state%tend_pr_tot_qg_pwad(:)=current_state%tend_pr_tot_qg_pwad(:) + &
      current_state%tend_3d_qg_pwad(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_tabs) then
      current_state%tend_pr_tot_tabs_pwad(:)=current_state%tend_pr_tot_tabs_pwad(:) + &
      current_state%tend_3d_tabs_pwad(:,cyn,cxn)
    endif

  end subroutine compute_component_tendencies
  
  
  !> Parses a field string (read in from the configuration file) and determines whether this algorithm should be used
  !! for advecting that field
  !! @param field The string configuration of field advection
  !! @returns Whether or not the field is advected here
  logical function determine_if_advection_here(field)
    character(len=*), intent(in) :: field

    if (len_trim(field) .ne. 0) then
      if (trim(field) .eq. "pw" .or. trim(field) .eq. "any") then
        determine_if_advection_here=.true.
      else
        determine_if_advection_here=.false.
      end if
    else
      determine_if_advection_here=.true.
    end if
  end function determine_if_advection_here

end module pwadvection_mod

