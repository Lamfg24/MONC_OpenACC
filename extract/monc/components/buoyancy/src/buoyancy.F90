!> Calculates buoyancy terms for the SW field
module buoyancy_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type, &
      component_descriptor_type_v1
  use optionsdatabase_mod, only : options_has_key, options_get_real_array, options_get_integer
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, X_INDEX, Y_INDEX
  use science_constants_mod
  use q_indices_mod, only: get_q_index, standard_q_names
implicit none

#ifndef TEST_MODE
  private
#endif
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: w_buoyancy

  real(kind=DEFAULT_PRECISION) :: G_over_2

  integer :: iqv ! Index for water vapour

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  logical :: l_tend_3d_w
  ! Local mean tendency profile fields and logicals for their use
  logical :: l_tend_pr_tot_w

  public initialisation_callback_buoyancy, timestep_callback_buoyancy, &
         finalisation_callback_buoyancy

contains


  !> The initialisation callback sets up the buoyancy coefficient
  !! @param current_state The current model state
  subroutine initialisation_callback_buoyancy(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: z_size
    integer :: alloc_z, alloc_y, alloc_x, i, iter

    if (.not. current_state%passive_q .and. current_state%number_q_fields > 0)then
      if (.not. allocated(current_state%cq))then
        allocate(current_state%cq(21))
        current_state%cq=0.0_DEFAULT_PRECISION
      end if
      !iqv = get_q_index(standard_q_names%VAPOUR, 'buoyancy')
      !current_state%cq(iqv) = ratio_mol_wts-1.0
      current_state%cq(1) = ratio_mol_wts-1.0
    end if

    G_over_2 = 0.5_DEFAULT_PRECISION*G
    z_size=current_state%global_grid%size(Z_INDEX)
    allocate(w_buoyancy(z_size))

    w_buoyancy = 0.0_DEFAULT_PRECISION

    ! Tendency logicals
    l_tend_pr_tot_w   = current_state%w%active
    l_tend_3d_w   = current_state%w%active .or. l_tend_pr_tot_w

    alloc_z=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    alloc_y=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    alloc_x=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_w) then
      !allocate( current_state%tend_3d_w_buoy(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_w_buoy(alloc_z, alloc_y, alloc_x))
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_w) then
      !allocate( current_state%tend_pr_tot_w_buoy(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_w_buoy(alloc_z))
    endif

  end subroutine initialisation_callback_buoyancy


  subroutine finalisation_callback_buoyancy(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(w_buoyancy)) deallocate(w_buoyancy)
    if (allocated(current_state%tend_3d_w_buoy)) deallocate(current_state%tend_3d_w_buoy)
    if (allocated(current_state%tend_pr_tot_w_buoy)) deallocate(current_state%tend_pr_tot_w_buoy)

  end subroutine finalisation_callback_buoyancy


  !> Called for each column per timestep this will calculate the buoyancy terms for the SW field
  !! @param current_state The current model state
  !! @param target_(x/y)_index This is the index with the halos subtracted. This is needed so that diagnostic does
  !!                           not include halos and to prevent array out-of-bounds
  subroutine timestep_callback_buoyancy(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, n, i, j
    integer :: current_x_index, current_y_index, target_x_index, target_y_index
    logical :: calculate_diagnostics

    !calculate_diagnostics = current_state%diagnostic_sample_timestep
    i=current_state%column_local_x
    j=current_state%column_local_y
    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_w) then
        current_state%tend_pr_tot_w_buoy(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    !if (calculate_diagnostics) &
    if ((current_state%modulo_number_3d .eq. 0) &
        .and. (.not. current_state%halo_column)) then
        call save_precomponent_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)
    end if
    
#ifdef W_ACTIVE
    if (.not. current_state%passive_th .and. current_state%th%active) then
      if (current_state%immersed%ib_enabled .eqv. .true.) then
        if (current_state%immersed%ib_col(j,i) .eqv. .true.) then
          do k=2,current_state%local_grid%size(Z_INDEX)-1
            if (current_state%immersed%indic_w(k,j,i) .eq. 0 &
            .and. current_state%immersed%indic_s(k,j,i) .eq. 0 &
            .and. current_state%immersed%indic_s(k+1,j,i) .eq. 0) then
              current_state%w_buoyancy(k)=(0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                (current_state%th%data(k,j,i)+current_state%th%data(k+1,j,i))
                  current_state%sw%data(k,j,i)=current_state%sw%data(k,j,i)+current_state%w_buoyancy(k)
            end if
          end do
        end if
      else
        do k=2,current_state%local_grid%size(Z_INDEX)-1
          w_buoyancy(k)=(0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
              (current_state%th%data(k, current_state%column_local_y, current_state%column_local_x)&
              +current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
              current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+w_buoyancy(k)
        end do
      end if
    end if
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
      if (current_state%use_anelastic_equations) then
        do k=2,current_state%local_grid%size(Z_INDEX)-1
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(1)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qv%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qv%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(2)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%ql%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%ql%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(3)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qr%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qr%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(4)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qi%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qi%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(5)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qs%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qs%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(6)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qg%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qg%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(7)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qAitkenSolMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qAitkenSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(8)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qAccumSolMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qAccumSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(9)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qAccumInsolMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qAccumInsolMass%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(10)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qCoarseSolMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qCoarseSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(11)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%qCoarseDustMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%qCoarseDustMass%data(k+1, current_state%column_local_y, current_state%column_local_x))

          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(12)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%nl%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%nl%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(13)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%nr%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%nr%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(14)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%ni%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%ni%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(15)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%ns%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%ns%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(16)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%ng%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%ng%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(17)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%nAitkenSolNumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%nAitkenSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(18)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%nAccumSolNumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%nAccumSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(19)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%nAccumInsolNumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%nAccumInsolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(20)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%nCoarseSolNumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%nCoarseSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k))*&
                current_state%cq(21)* (current_state%global_grid%configuration%vertical%thref(k)*&
                current_state%nCoarseDustnumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                current_state%global_grid%configuration%vertical%thref(k+1)*&
                current_state%nCoarseDustnumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
        end do
      else
        do k=2,current_state%local_grid%size(Z_INDEX)-1
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(1)*&
                 (current_state%qv%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qv%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(2)*&
                 (current_state%ql%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%ql%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(3)*&
                 (current_state%qr%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qr%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(4)*&
                 (current_state%qi%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qi%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(5)*&
                 (current_state%qs%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qs%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(6)*&
                 (current_state%qg%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qg%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(7)*&
                 (current_state%qAitkenSolMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qAitkenSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(8)*&
                 (current_state%qAccumSolMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qAccumSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(9)*&
                 (current_state%qAccumInsolMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qAccumInsolMass%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(10)*&
                 (current_state%qCoarseSolMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qCoarseSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(11)*&
                 (current_state%qCoarseDustMass%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%qCoarseDustMass%data(k+1, current_state%column_local_y, current_state%column_local_x))

          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(12)*&
                 (current_state%nl%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%nl%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(13)*&
                 (current_state%nr%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%nr%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(14)*&
                 (current_state%ni%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%ni%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(15)*&
                 (current_state%ns%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%ns%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(16)*&
                 (current_state%ng%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%ng%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(11)*&
                 (current_state%nAitkenSolNumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%nAitkenSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(11)*&
                 (current_state%nAccumSolNumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%nAccumSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(11)*&
                 (current_state%nAccumInsolNumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%nAccumInsolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(11)*&
                 (current_state%nCoarseSolNumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%nCoarseSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
          current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)=&
                 current_state%sw%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 G_over_2*current_state%cq(11)*&
                 (current_state%nCoarseDustnumber%data(k, current_state%column_local_y, current_state%column_local_x)+&
                 current_state%nCoarseDustnumber%data(k+1, current_state%column_local_y, current_state%column_local_x))
        end do
      end if
    end if
#endif

    !if (calculate_diagnostics) &
    if ((current_state%modulo_number_3d .eq. 0.0) &
        .and. (.not. current_state%halo_column)) then
        call compute_component_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback_buoyancy


   !> Save the 3d tendencies coming into the component.
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  !subroutine save_precomponent_tendencies(current_state, cxn, cyn, txn, tyn)
  subroutine save_precomponent_tendencies(current_state, cxn, cyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  cxn, cyn!, txn, tyn

    ! Save 3d tendency fields upon request (of 3d or profiles) and availability
    if (l_tend_3d_w) then
      current_state%tend_3d_w_buoy(:,cyn,cxn)=current_state%sw%data(:,cyn,cxn)
      !print*,"current_state%tend_3d_w_buoy(:,cyn,cxn) = ",current_state%tend_3d_w_buoy(:,cyn,cxn)
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
    if (l_tend_3d_w) then
      current_state%tend_3d_w_buoy(:,cyn,cxn)=current_state%sw%data(:,cyn,cxn) - current_state%tend_3d_w_buoy(:,cyn,cxn)
      !print*,"current_state%tend_3d_w_buoy(:,cyn,cxn) = ",current_state%tend_3d_w_buoy(:,cyn,cxn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_w) then
      current_state%tend_pr_tot_w_buoy(:)=current_state%tend_pr_tot_w_buoy(:) + current_state%tend_3d_w_buoy(:,cyn,cxn)
      !print*,"tend_pr_tot_w_buoy(:) = ",current_state%tend_pr_tot_w_buoy(:)
    endif

  end subroutine compute_component_tendencies

end module buoyancy_mod
