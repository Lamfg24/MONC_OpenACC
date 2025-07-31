!> Diffusion on the TH and Q fields
module diffusion_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type, &
      component_descriptor_type_v1
  use optionsdatabase_mod, only : options_get_integer
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : Z_INDEX, X_INDEX, Y_INDEX
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, perform_local_data_copy_for_field, complete_nonblocking_halo_swap, &
       copy_buffer_to_corner
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_master_log
  use q_indices_mod, only: get_q_index, standard_q_names

  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: th_diffusion, qv_diffusion, ql_diffusion, qr_diffusion, &
                               qi_diffusion, qs_diffusion, qg_diffusion, qAitkenSolMass_diffusion, qAccumSolMass_diffusion, &
                               qAccumInsolMass_diffusion, qCoarseSolMass_diffusion, qCoarseDustMass_diffusion, &
                               nl_diffusion, nr_diffusion, ni_diffusion, ns_diffusion, ng_diffusion, &
                               nAitkenSolNumber_diffusion, nAccumSolNumber_diffusion, nAccumInsolNumber_diffusion, &
                               nCoarseSolNumber_diffusion, nCoarseDustnumber_diffusion
  !real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_diffusion
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: tracer_diffusion

  ! Local tendency diagnostic variables for this component
  logical :: l_tend_3d_th,l_tend_3d_qv,l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,l_tend_3d_tabs
  logical :: l_tend_pr_tot_th,l_tend_pr_tot_qv,l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,  &
             l_tend_pr_tot_qs,l_tend_pr_tot_qg,l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0

  public initialisation_callback_diffusion, timestep_callback_diffusion, &
         finalisation_callback_diffusion

contains

  !> Sets up the stencil_mod (used in interpolation) and allocates data for the flux fields
  !! @param current_state The current model state_mod
  subroutine initialisation_callback_diffusion(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: z_size, y_size, x_size
    integer :: alloc_z, alloc_y, alloc_x
    logical :: l_qdiag

    z_size=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    y_size=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    x_size=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    allocate(current_state%diff_coefficient%data(z_size, y_size, x_size))

    z_size=current_state%global_grid%size(Z_INDEX)
    allocate(th_diffusion(z_size))
    allocate(qv_diffusion(z_size))
    allocate(ql_diffusion(z_size))
    allocate(qr_diffusion(z_size))
    allocate(qi_diffusion(z_size))
    allocate(qs_diffusion(z_size))
    allocate(qg_diffusion(z_size))
    allocate(qAitkenSolMass_diffusion(z_size))
    allocate(qAccumSolMass_diffusion(z_size))
    allocate(qAccumInsolMass_diffusion(z_size))
    allocate(qCoarseSolMass_diffusion(z_size))
    allocate(qCoarseDustMass_diffusion(z_size))
    allocate(nl_diffusion(z_size))
    allocate(nr_diffusion(z_size))
    allocate(ni_diffusion(z_size))
    allocate(ns_diffusion(z_size))
    allocate(ng_diffusion(z_size))
    allocate(nAitkenSolNumber_diffusion(z_size))
    allocate(nAccumSolNumber_diffusion(z_size))
    allocate(nAccumInsolNumber_diffusion(z_size))
    allocate(nCoarseSolNumber_diffusion(z_size))
    allocate(nCoarseDustnumber_diffusion(z_size))

    !allocate(q_diffusion(z_size, current_state%number_q_fields))
    if (current_state%n_tracers > 0) allocate(tracer_diffusion(z_size, current_state%n_tracers))

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available
    l_qdiag =  (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0)

    l_tend_pr_tot_th  = current_state%th%active 
    l_tend_pr_tot_qv  = l_qdiag .and. current_state%water_vapour_mixing_ratio_index > 0
    l_tend_pr_tot_ql  = l_qdiag .and. current_state%liquid_water_mixing_ratio_index > 0
    l_tend_pr_tot_qi  = l_qdiag .and. current_state%ice_water_mixing_ratio_index > 0
    l_tend_pr_tot_qr  = l_qdiag .and. current_state%rain_water_mixing_ratio_index > 0
    l_tend_pr_tot_qs  = l_qdiag .and. current_state%snow_water_mixing_ratio_index > 0
    l_tend_pr_tot_qg  = l_qdiag .and. current_state%graupel_water_mixing_ratio_index > 0
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_qv  = (l_qdiag .and. current_state%water_vapour_mixing_ratio_index > 0) .or. l_tend_pr_tot_qv
    l_tend_3d_ql  = (l_qdiag .and. current_state%liquid_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_ql
    l_tend_3d_qi  = (l_qdiag .and. current_state%ice_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qi
    l_tend_3d_qr  = (l_qdiag .and. current_state%rain_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qr
    l_tend_3d_qs  = (l_qdiag .and. current_state%snow_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qs
    l_tend_3d_qg  = (l_qdiag .and. current_state%graupel_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qg
    l_tend_3d_tabs = l_tend_3d_th

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_th) then
      !allocate( current_state%tend_3d_th_diff(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_th_diff(z_size, y_size, x_size))
    endif
    if (l_tend_3d_qv) then
      iqv=1!get_q_index(standard_q_names%VAPOUR, 'diffusion')
      !allocate( current_state%tend_3d_qv_diff(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qv_diff(z_size, y_size, x_size))
    endif
    if (l_tend_3d_ql) then
      iql=2!get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'diffusion')
      !allocate( current_state%tend_3d_ql_diff(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_ql_diff(z_size, y_size, x_size))
    endif
    if (l_tend_3d_qi) then
      iqi=3!get_q_index(standard_q_names%ICE_MASS, 'diffusion')
      !allocate( current_state%tend_3d_qi_diff(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qi_diff(z_size, y_size, x_size))
    endif
    if (l_tend_3d_qr) then
      iqr=4!get_q_index(standard_q_names%RAIN_MASS, 'diffusion')
      !allocate( current_state%tend_3d_qr_diff(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qr_diff(z_size, y_size, x_size))
    endif
    if (l_tend_3d_qs) then
      iqs=5!get_q_index(standard_q_names%SNOW_MASS, 'diffusion')
      !allocate( current_state%tend_3d_qs_diff(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qs_diff(z_size, y_size, x_size))
    endif
    if (l_tend_3d_qg) then
      iqg=6!get_q_index(standard_q_names%GRAUPEL_MASS, 'diffusion')
      !allocate( current_state%tend_3d_qg_diff(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_qg_diff(z_size, y_size, x_size))
    endif
    if (l_tend_3d_tabs) then
      !allocate( current_state%tend_3d_tabs_diff(current_state%local_grid%size(Z_INDEX),  &
      !                       current_state%local_grid%size(Y_INDEX),  &
      !                       current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_tabs_diff(z_size, y_size, x_size))
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_th) then
      !allocate( current_state%tend_pr_tot_th_diff(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_th_diff(z_size))
    endif
    if (l_tend_pr_tot_qv) then
      !allocate( current_state%tend_pr_tot_qv_diff(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qv_diff(z_size))
    endif
    if (l_tend_pr_tot_ql) then
      !allocate( current_state%tend_pr_tot_ql_diff(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_ql_diff(z_size))
    endif
    if (l_tend_pr_tot_qi) then
      !allocate( current_state%tend_pr_tot_qi_diff(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qi_diff(z_size))
    endif
    if (l_tend_pr_tot_qr) then
      !allocate( current_state%tend_pr_tot_qr_diff(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qr_diff(z_size))
    endif
    if (l_tend_pr_tot_qs) then
      !allocate( current_state%tend_pr_tot_qs_diff(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qs_diff(z_size))
    endif
    if (l_tend_pr_tot_qg) then
      !allocate( current_state%tend_pr_tot_qg_diff(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_qg_diff(z_size))
    endif
    if (l_tend_pr_tot_tabs) then
      !allocate( current_state%tend_pr_tot_tabs_diff(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_tabs_diff(z_size))
    endif

  end subroutine initialisation_callback_diffusion


  subroutine finalisation_callback_diffusion(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(th_diffusion)) deallocate(th_diffusion)
    !if (allocated(q_diffusion)) deallocate(q_diffusion)
    if (allocated(tracer_diffusion)) deallocate(tracer_diffusion)

    if (allocated(current_state%tend_3d_th_diff)) deallocate(current_state%tend_3d_th_diff)
    if (allocated(current_state%tend_3d_qv_diff)) deallocate(current_state%tend_3d_qv_diff)
    if (allocated(current_state%tend_3d_ql_diff)) deallocate(current_state%tend_3d_ql_diff)
    if (allocated(current_state%tend_3d_qi_diff)) deallocate(current_state%tend_3d_qi_diff)
    if (allocated(current_state%tend_3d_qr_diff)) deallocate(current_state%tend_3d_qr_diff)
    if (allocated(current_state%tend_3d_qs_diff)) deallocate(current_state%tend_3d_qs_diff)
    if (allocated(current_state%tend_3d_qg_diff)) deallocate(current_state%tend_3d_qg_diff)
    if (allocated(current_state%tend_3d_tabs_diff)) deallocate(current_state%tend_3d_tabs_diff)

    if (allocated(current_state%tend_pr_tot_th_diff)) deallocate(current_state%tend_pr_tot_th_diff)
    if (allocated(current_state%tend_pr_tot_qv_diff)) deallocate(current_state%tend_pr_tot_qv_diff)
    if (allocated(current_state%tend_pr_tot_ql_diff)) deallocate(current_state%tend_pr_tot_ql_diff)
    if (allocated(current_state%tend_pr_tot_qi_diff)) deallocate(current_state%tend_pr_tot_qi_diff)
    if (allocated(current_state%tend_pr_tot_qr_diff)) deallocate(current_state%tend_pr_tot_qr_diff)
    if (allocated(current_state%tend_pr_tot_qs_diff)) deallocate(current_state%tend_pr_tot_qs_diff)
    if (allocated(current_state%tend_pr_tot_qg_diff)) deallocate(current_state%tend_pr_tot_qg_diff)
    if (allocated(current_state%tend_pr_tot_tabs_diff)) deallocate(current_state%tend_pr_tot_tabs_diff)

  end subroutine finalisation_callback_diffusion

  !> At each timestep will compute the diffusion source terms for TH and Q fields per column if these fields are active
  !! @param current_state The current model state
  subroutine timestep_callback_diffusion(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: local_y, local_x, target_x_index, target_y_index
    logical :: calculate_diagnostics

    !calculate_diagnostics = current_state%diagnostic_sample_timestep

    local_y=current_state%column_local_y
    local_x=current_state%column_local_x
    target_y_index=local_y-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=local_x-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_th) then
        current_state%tend_pr_tot_th_diff(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qv) then
        current_state%tend_pr_tot_qv_diff(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_ql) then
        current_state%tend_pr_tot_ql_diff(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qi) then
        current_state%tend_pr_tot_qi_diff(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qr) then
        current_state%tend_pr_tot_qr_diff(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qs) then
        current_state%tend_pr_tot_qs_diff(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qg) then
        current_state%tend_pr_tot_qg_diff(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_tabs) then
        current_state%tend_pr_tot_tabs_diff(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    if (.not. current_state%use_viscosity_and_diffusion .or. current_state%halo_column) return

    if (current_state%diffusion_halo_swap_state%swap_in_progress) then
      ! If there is a diffusion halo swap in progress then complete it
      call complete_nonblocking_halo_swap(current_state, current_state%diffusion_halo_swap_state, &
           perform_local_data_copy_for_diff, copy_halo_buffer_to_diff, copy_halo_buffer_to_diff_corners)
    end if

    if (current_state%modulo_number_3d .eq. 0) &
      call save_precomponent_tendencies(current_state, local_x, local_y)!, target_x_index, target_y_index)

    if (current_state%th%active) call perform_th_diffusion(current_state, local_y, local_x)
    if (current_state%number_q_fields .gt. 0) call perform_q_diffusion(current_state, local_y, local_x)
    if (current_state%n_radioactive_tracers .gt. 0) call perform_tracer_diffusion(current_state, local_y, local_x)

    if (current_state%modulo_number_3d .eq. 0) &
      call compute_component_tendencies(current_state, local_x, local_y)!, target_x_index, target_y_index)
  end subroutine timestep_callback_diffusion

  !> Computes the diffusion source terms for each tracer field
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine perform_tracer_diffusion(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: n

    do n=current_state%radioactive_tracer_index, current_state%radioactive_tracer_index + current_state%n_radioactive_tracers - 1 
      call general_diffusion(current_state, current_state%ztracer(n), current_state%stracer(n), local_y, local_x, &
        tracer_diffusion(:,n))
    end do
  end subroutine perform_tracer_diffusion  


  !> Computes the diffusion source terms for each Q field
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine perform_q_diffusion(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: n

    call general_diffusion(current_state, current_state%zqv, current_state%sqv, local_y, local_x, qv_diffusion)
    call general_diffusion(current_state, current_state%zql, current_state%sql, local_y, local_x, ql_diffusion)
    call general_diffusion(current_state, current_state%zqr, current_state%sqr, local_y, local_x, qr_diffusion)
    call general_diffusion(current_state, current_state%zqi, current_state%sqi, local_y, local_x, qi_diffusion)
    call general_diffusion(current_state, current_state%zqs, current_state%sqs, local_y, local_x, qs_diffusion)
    call general_diffusion(current_state, current_state%zqg, current_state%sqg, local_y, local_x, qg_diffusion)
    call general_diffusion(current_state, current_state%zqAitkenSolMass, current_state%sqAitkenSolMass, &
                           local_y, local_x, qAitkenSolMass_diffusion)
    call general_diffusion(current_state, current_state%zqAccumSolMass, current_state%sqAccumSolMass, &
                           local_y, local_x, qAccumSolMass_diffusion)
    call general_diffusion(current_state, current_state%zqAccumInsolMass, current_state%sqAccumInsolMass, &
                           local_y, local_x, qAccumInsolMass_diffusion)
    call general_diffusion(current_state, current_state%zqCoarseSolMass, current_state%sqCoarseSolMass, &
                           local_y, local_x, qCoarseSolMass_diffusion)
    call general_diffusion(current_state, current_state%zqCoarseDustMass, current_state%sqCoarseDustMass, &
                           local_y, local_x, qCoarseDustMass_diffusion)

    call general_diffusion(current_state, current_state%znl, current_state%snl, local_y, local_x, nl_diffusion)
    call general_diffusion(current_state, current_state%znr, current_state%snr, local_y, local_x, nr_diffusion)
    call general_diffusion(current_state, current_state%zni, current_state%sni, local_y, local_x, ni_diffusion)
    call general_diffusion(current_state, current_state%zns, current_state%sns, local_y, local_x, ns_diffusion)
    call general_diffusion(current_state, current_state%zng, current_state%sng, local_y, local_x, ng_diffusion)
    call general_diffusion(current_state, current_state%znAitkenSolNumber, current_state%snAitkenSolNumber, &
                           local_y, local_x, nAitkenSolNumber_diffusion)
    call general_diffusion(current_state, current_state%znAccumSolNumber, current_state%snAccumSolNumber, &
                           local_y, local_x, nAccumSolNumber_diffusion)
    call general_diffusion(current_state, current_state%znAccumInsolNumber, current_state%snAccumInsolNumber, &
                           local_y, local_x, nAccumInsolNumber_diffusion)
    call general_diffusion(current_state, current_state%znCoarseSolNumber, current_state%snCoarseSolNumber, &
                           local_y, local_x, nCoarseSolNumber_diffusion)
    call general_diffusion(current_state, current_state%znCoarseDustnumber, current_state%snCoarseDustnumber, &
                           local_y, local_x, nCoarseDustnumber_diffusion)
    !do n=1, current_state%number_q_fields
    !  call general_diffusion(current_state, current_state%zq(n), current_state%sq(n), local_y, local_x, q_diffusion(:,n))
    !end do
  end subroutine perform_q_diffusion  

  !> Computes the diffusion source terms for the theta field
  !! @param current_state The current model state
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine perform_th_diffusion(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_y, local_x

    integer :: k
    real(kind=DEFAULT_PRECISION) :: adjustment

    call general_diffusion(current_state, current_state%zth, current_state%sth, local_y, local_x, th_diffusion)

    if (current_state%use_anelastic_equations) then
     ! This code only needs to be executed if anelastic, otherwise THREF is constant and the increment calculated here is zero
      do k=2, current_state%local_grid%size(Z_INDEX)
        adjustment=(current_state%global_grid%configuration%vertical%cza(k)*&
             current_state%global_grid%configuration%vertical%dthref(k)*&
             current_state%diff_coefficient%data(k, local_y, local_x) - current_state%global_grid%configuration%vertical%czb(k)*&
             current_state%global_grid%configuration%vertical%dthref(k-1)*&
             current_state%diff_coefficient%data(k-1, local_y, local_x))
        current_state%sth%data(k, local_y, local_x)=current_state%sth%data(k, local_y, local_x)+adjustment
        th_diffusion(k)=th_diffusion(k)+adjustment
      end do
    end if

  end subroutine perform_th_diffusion

  !> General diffusion computation for any field which is provided as arguments. Works in a column
  !! @param current_state The current model state
  !! @param field The field to take values from, typically zth or zq(n)
  !! @param source_field The source target field to update, typically sth or sq(n)
  !! @param local_y Local Y index
  !! @param local_x Local X index
  subroutine general_diffusion(current_state, field, source_field, local_y, local_x, diagnostics)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: field, source_field
    integer, intent(in) :: local_y, local_x
    real(kind=DEFAULT_PRECISION), dimension(:), intent(inout), optional :: diagnostics

    real(kind=DEFAULT_PRECISION) :: term
    integer :: k

    do k=2, current_state%local_grid%size(Z_INDEX)
      term=current_state%global_grid%configuration%horizontal%cx2*0.25_DEFAULT_PRECISION*&
           (((current_state%diff_coefficient%data(k, local_y, local_x)+&
           current_state%diff_coefficient%data(k, local_y, local_x-1))+&
           (current_state%diff_coefficient%data(k-1, local_y, local_x)+&
           current_state%diff_coefficient%data(k-1, local_y, local_x-1)))&
           *(field%data(k, local_y, local_x-1)-field%data(k, local_y, local_x)) -&
           ((current_state%diff_coefficient%data(k, local_y, local_x+1)+&
           current_state%diff_coefficient%data(k, local_y, local_x))+&
           (current_state%diff_coefficient%data(k-1, local_y, local_x+1)+&
           current_state%diff_coefficient%data(k-1, local_y, local_x)))&
           *(field%data(k, local_y, local_x)-field%data(k, local_y, local_x+1)) )&
           +current_state%global_grid%configuration%horizontal%cy2*0.25_DEFAULT_PRECISION*(&
           ((current_state%diff_coefficient%data(k, local_y, local_x)+&
           current_state%diff_coefficient%data(k, local_y-1, local_x))+&
           (current_state%diff_coefficient%data(k-1, local_y, local_x)+&
           current_state%diff_coefficient%data(k-1, local_y-1, local_x)))&
           *(field%data(k, local_y-1, local_x)-field%data(k, local_y, local_x)) -&
           ((current_state%diff_coefficient%data(k, local_y+1, local_x)+&
           current_state%diff_coefficient%data(k, local_y, local_x))+&
           (current_state%diff_coefficient%data(k-1, local_y+1, local_x)+&
           current_state%diff_coefficient%data(k-1, local_y, local_x)))&
           *(field%data(k, local_y, local_x)-field%data(k, local_y+1, local_x)) )&
           +( current_state%global_grid%configuration%vertical%czb(k)*&
           current_state%diff_coefficient%data(k-1, local_y, local_x)*&
           (field%data(k-1, local_y, local_x)-field%data(k, local_y, local_x)))
      if (k .lt. current_state%local_grid%size(Z_INDEX)) then
        term=term - current_state%global_grid%configuration%vertical%cza(k)*&
             current_state%diff_coefficient%data(k, local_y, local_x)*&
             (field%data(k, local_y, local_x)-field%data(k+1, local_y, local_x))
      end if
      source_field%data(k, local_y, local_x)=source_field%data(k, local_y, local_x)+term
      if (present(diagnostics)) diagnostics(k)=term
    end do
  end subroutine general_diffusion

    !> Does local data copying for diffusion coefficient variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_diff(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%diff_coefficient%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_diff

  !> Copies the halo buffer to halo location for the diffusion coefficient field
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_diff(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%diff_coefficient%data, dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_diff

  !> Copies the corner halo buffer to the diffusion coefficient field corners
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param corner_loc The corner location
  !! @param x_target_index The X target index for the dimension we are receiving for
  !! @param y_target_index The Y target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_diff_corners(current_state, neighbour_description, corner_loc, x_target_index, &
       y_target_index, neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, x_target_index, y_target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%diff_coefficient%data, corner_loc, x_target_index, y_target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_diff_corners  

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
      current_state%tend_3d_th_diff(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      current_state%tend_3d_qv_diff(:,cyn,cxn)=current_state%sqv%data(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      current_state%tend_3d_ql_diff(:,cyn,cxn)=current_state%sql%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qi) then
      current_state%tend_3d_qi_diff(:,cyn,cxn)=current_state%sqi%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qr) then
      current_state%tend_3d_qr_diff(:,cyn,cxn)=current_state%sqr%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qs) then
      current_state%tend_3d_qs_diff(:,cyn,cxn)=current_state%sqs%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qg) then
      current_state%tend_3d_qg_diff(:,cyn,cxn)=current_state%sqg%data(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      current_state%tend_3d_tabs_diff(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn) * &
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
      current_state%tend_3d_th_diff(:,cyn,cxn)=current_state%sth%data(:,cyn,cxn) - &
      current_state%tend_3d_th_diff(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      current_state%tend_3d_qv_diff(:,cyn,cxn)=current_state%sqv%data(:,cyn,cxn) - &
      current_state%tend_3d_qv_diff(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      current_state%tend_3d_ql_diff(:,cyn,cxn)=current_state%sql%data(:,cyn,cxn) - &
      current_state%tend_3d_ql_diff(:,cyn,cxn)
    endif
    if (l_tend_3d_qi) then
      current_state%tend_3d_qi_diff(:,cyn,cxn)=current_state%sqi%data(:,cyn,cxn) - &
      current_state%tend_3d_qi_diff(:,cyn,cxn)
    endif
    if (l_tend_3d_qr) then
      current_state%tend_3d_qr_diff(:,cyn,cxn)=current_state%sqr%data(:,cyn,cxn) - &
      current_state%tend_3d_qr_diff(:,cyn,cxn)
    endif
    if (l_tend_3d_qs) then
      current_state%tend_3d_qs_diff(:,cyn,cxn)=current_state%sqs%data(:,cyn,cxn) - &
      current_state%tend_3d_qs_diff(:,cyn,cxn)
    endif
    if (l_tend_3d_qg) then
      current_state%tend_3d_qg_diff(:,cyn,cxn)=current_state%sqg%data(:,cyn,cxn) - &
      current_state%tend_3d_qg_diff(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      current_state%tend_3d_tabs_diff(:,cyn,cxn)=   &
       current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)   &
        - current_state%tend_3d_tabs_diff(:,cyn,cxn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_th) then
      current_state%tend_pr_tot_th_diff(:)=current_state%tend_pr_tot_th_diff(:) + &
      current_state%tend_3d_th_diff(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qv) then
      current_state%tend_pr_tot_qv_diff(:)=current_state%tend_pr_tot_qv_diff(:) + &
      current_state%tend_3d_qv_diff(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_ql) then
      current_state%tend_pr_tot_ql_diff(:)=current_state%tend_pr_tot_ql_diff(:) + &
      current_state%tend_3d_ql_diff(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qi) then
      current_state%tend_pr_tot_qi_diff(:)=current_state%tend_pr_tot_qi_diff(:) + &
      current_state%tend_3d_qi_diff(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qr) then
      current_state%tend_pr_tot_qr_diff(:)=current_state%tend_pr_tot_qr_diff(:) + &
      current_state%tend_3d_qr_diff(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qs) then
      current_state%tend_pr_tot_qs_diff(:)=current_state%tend_pr_tot_qs_diff(:) + &
      current_state%tend_3d_qs_diff(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_qg) then
      current_state%tend_pr_tot_qg_diff(:)=current_state%tend_pr_tot_qg_diff(:) + &
      current_state%tend_3d_qg_diff(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_tabs) then
      current_state%tend_pr_tot_tabs_diff(:)=current_state%tend_pr_tot_tabs_diff(:) + &
      current_state%tend_3d_tabs_diff(:,cyn,cxn)
    endif

  end subroutine compute_component_tendencies

end module diffusion_mod
