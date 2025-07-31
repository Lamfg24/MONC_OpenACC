










!> Does the swapping and smoothing which is called for each column as part of the pressure-terms group of components
!! Note that this does not currently implement smoothing for mean profiles
module swapsmooth_mod
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use state_mod, only : model_state_type, FORWARD_STEPPING
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  use saturation_mod, only: qsaturation, qisaturation
  implicit none

  private

  logical :: mean_profiles_active=.false. !< Whether or not mean profiles need smoothing

  public initialisation_callback_swapsmooth, timestep_callback_swapsmooth

contains

  !> Initialises the swap and smooth component
  !! @param current_state The current model state
  subroutine initialisation_callback_swapsmooth(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olubar)
    return
    mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olvbar)
    return
    if (current_state%th%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olthbar)   
      return
    end if

    if (current_state%qv%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqvbar)
      return
    end if

    if (current_state%ql%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqlbar)
      return
    end if

    if (current_state%qr%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqrbar)
      return
    end if

    if (current_state%qi%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqibar)
      return
    end if

    if (current_state%qs%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqsbar)
      return
    end if

    if (current_state%qg%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqgbar)
      return
    end if

    if (current_state%qAitkenSolMass%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqAitkenSolMassbar)
      return
    end if

    if (current_state%qAccumSolMass%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqAccumSolMassbar)
      return
    end if

    if (current_state%qAccumInsolMass%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqAccumInsolMassbar)
      return
    end if

    if (current_state%qCoarseSolMass%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqCoarseSolMassbar)
      return
    end if

    if (current_state%qCoarseDustMass%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olqCoarseDustMassbar)
      return
    end if

    !! NUMBER
    if (current_state%nl%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olnlbar)
      return
    end if

    if (current_state%nr%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olnrbar)
      return
    end if

    if (current_state%ni%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olnibar)
      return
    end if

    if (current_state%ns%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olnsbar)
      return
    end if

    if (current_state%ng%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olngbar)
      return
    end if

    if (current_state%nAitkenSolNumber%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar)
      return
    end if

    if (current_state%nAccumSolNumber%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olnAccumSolNumberbar)
      return
    end if

    if (current_state%nAccumInsolNumber%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar)
      return
    end if

    if (current_state%nCoarseSolNumber%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar)
      return
    end if

    if (current_state%nCoarseDustnumber%active) then
      mean_profiles_active=allocated(current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar)
      return
    end if

  end subroutine initialisation_callback_swapsmooth

  !> Called for each non halo timestep column and will perform swapping and smoothing as required on that column
  !! @param current_state The current model state_mod
  subroutine timestep_callback_swapsmooth(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (.not. current_state%halo_column) then
      if (current_state%field_stepping == FORWARD_STEPPING) then
        call swap_and_smooth_classic(current_state, .false.)
      else
        ! Centred stepping
        call swap_and_smooth_robert_filter(current_state)
      end if
    end if

    if (mean_profiles_active .and. current_state%last_timestep_column) then
      if (current_state%field_stepping == FORWARD_STEPPING) then
        call classic_for_average_profiles(current_state, .false.)
      else
        call robert_filter_for_average_profiles(current_state)
      end if
    end if
  end subroutine timestep_callback_swapsmooth

  !> Classic swap and smooth based upon the old or no smoothing
  !! @param current_state The current model state_mod
  !! @param old_smoother Whether to use the old smoother or not (not means no smoothing)
  subroutine swap_and_smooth_classic(current_state, old_smoother)
    type(model_state_type), intent(inout) :: current_state
    logical, intent(in) :: old_smoother

    integer :: y_index, x_index, k, n
    real(kind=DEFAULT_PRECISION) :: c1, c2, existing_value

    if (old_smoother) then
      c1 = 1.0_DEFAULT_PRECISION - current_state%tsmth                                               
      c2 = 2.0_DEFAULT_PRECISION * current_state%tsmth - 1.0_DEFAULT_PRECISION
    else
      c1 = 1.0_DEFAULT_PRECISION
      c2 = -1.0_DEFAULT_PRECISION
    end if

    x_index=current_state%column_local_x
    y_index=current_state%column_local_y
    do k=1,current_state%global_grid%size(Z_INDEX)
      existing_value = current_state%u%data(k,y_index,x_index) + current_state%zu%data(k,y_index,x_index)
      current_state%u%data(k,y_index,x_index)=existing_value * c1 + current_state%u%data(k,y_index,x_index) * c2
      current_state%zu%data(k,y_index,x_index)=existing_value - current_state%u%data(k,y_index,x_index)
      existing_value = current_state%v%data(k,y_index,x_index) + current_state%zv%data(k,y_index,x_index)
      current_state%v%data(k,y_index,x_index)=existing_value * c1 + current_state%v%data(k,y_index,x_index) * c2
      current_state%zv%data(k,y_index,x_index)=existing_value - current_state%v%data(k,y_index,x_index)
      existing_value = current_state%w%data(k,y_index,x_index) + current_state%zw%data(k,y_index,x_index)
      current_state%w%data(k,y_index,x_index)=existing_value * c1 + current_state%w%data(k,y_index,x_index) * c2
      current_state%zw%data(k,y_index,x_index)=existing_value - current_state%w%data(k,y_index,x_index)
      if (current_state%th%active) then
        existing_value = current_state%th%data(k,y_index,x_index) + current_state%zth%data(k,y_index,x_index)
        current_state%th%data(k,y_index,x_index)=existing_value * c1 + current_state%th%data(k,y_index,x_index) * c2
        current_state%zth%data(k,y_index,x_index)=existing_value - current_state%th%data(k,y_index,x_index)
        current_state%Tk%data(k,y_index, x_index) = (current_state%th%data(k,y_index, x_index)       &
            + current_state%global_grid%configuration%vertical%thref(k)) * &
            current_state%global_grid%configuration%vertical%rprefrcp(k)
        current_state%qv_saturation%data(k,y_index, x_index) = qsaturation(current_state%Tk%data(k,y_index, x_index), &
                            current_state%global_grid%configuration%vertical%prefn(k)/100.0)
        current_state%qi_saturation%data(k,y_index, x_index)= qisaturation(current_state%Tk%data(k,y_index, x_index), &
                            current_state%global_grid%configuration%vertical%prefn(k)/100.0)
      end if

      if (current_state%number_q_fields .gt. 0) then
        existing_value = current_state%qv%data(k,y_index,x_index) + current_state%zqv%data(k,y_index,x_index)
        current_state%qv%data(k,y_index,x_index)=existing_value * c1 + current_state%qv%data(k,y_index,x_index) * c2
        current_state%zqv%data(k,y_index,x_index)=existing_value - current_state%qv%data(k,y_index,x_index)

        current_state%RH%data(k,y_index, x_index) = current_state%qv%data(k,y_index,x_index) &
                                                    *100/current_state%qv_saturation%data(k,y_index, x_index)

        existing_value = current_state%ql%data(k,y_index,x_index) + current_state%zql%data(k,y_index,x_index)
        current_state%ql%data(k,y_index,x_index)=existing_value * c1 + current_state%ql%data(k,y_index,x_index) * c2
        current_state%zql%data(k,y_index,x_index)=existing_value - current_state%ql%data(k,y_index,x_index)

        existing_value = current_state%qr%data(k,y_index,x_index) + current_state%zqr%data(k,y_index,x_index)
        current_state%qr%data(k,y_index,x_index)=existing_value * c1 + current_state%qr%data(k,y_index,x_index) * c2
        current_state%zqr%data(k,y_index,x_index)=existing_value - current_state%qr%data(k,y_index,x_index)

        existing_value = current_state%qi%data(k,y_index,x_index) + current_state%zqi%data(k,y_index,x_index)
        current_state%qi%data(k,y_index,x_index)=existing_value * c1 + current_state%qi%data(k,y_index,x_index) * c2
        current_state%zqi%data(k,y_index,x_index)=existing_value - current_state%qi%data(k,y_index,x_index)

        current_state%RI%data(k,y_index, x_index) = current_state%qv%data(k,y_index,x_index) &
                                                    *100/current_state%qi_saturation%data(k,y_index, x_index)

        existing_value = current_state%qs%data(k,y_index,x_index) + current_state%zqs%data(k,y_index,x_index)
        current_state%qs%data(k,y_index,x_index)=existing_value * c1 + current_state%qs%data(k,y_index,x_index) * c2
        current_state%zqs%data(k,y_index,x_index)=existing_value - current_state%qs%data(k,y_index,x_index)

        existing_value = current_state%qg%data(k,y_index,x_index) + current_state%zqg%data(k,y_index,x_index)
        current_state%qg%data(k,y_index,x_index)=existing_value * c1 + current_state%qg%data(k,y_index,x_index) * c2
        current_state%zqg%data(k,y_index,x_index)=existing_value - current_state%qg%data(k,y_index,x_index)

        existing_value = current_state%qAitkenSolMass%data(k,y_index,x_index) + &
        current_state%zqAitkenSolMass%data(k,y_index,x_index)
        current_state%qAitkenSolMass%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%qAitkenSolMass%data(k,y_index,x_index) * c2
        current_state%zqAitkenSolMass%data(k,y_index,x_index)=existing_value - &
        current_state%qAitkenSolMass%data(k,y_index,x_index)

        existing_value = current_state%qAccumSolMass%data(k,y_index,x_index) + &
        current_state%zqAccumSolMass%data(k,y_index,x_index)
        current_state%qAccumSolMass%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%qAccumSolMass%data(k,y_index,x_index) * c2
        current_state%zqAccumSolMass%data(k,y_index,x_index)=existing_value - &
        current_state%qAccumSolMass%data(k,y_index,x_index)

        existing_value = current_state%qAccumInsolMass%data(k,y_index,x_index) + &
        current_state%zqAccumInsolMass%data(k,y_index,x_index)
        current_state%qAccumInsolMass%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%qAccumInsolMass%data(k,y_index,x_index) * c2
        current_state%zqAccumInsolMass%data(k,y_index,x_index)=existing_value - &
        current_state%qAccumInsolMass%data(k,y_index,x_index)

        existing_value = current_state%qCoarseSolMass%data(k,y_index,x_index) + &
        current_state%zqCoarseSolMass%data(k,y_index,x_index)
        current_state%qCoarseSolMass%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%qCoarseSolMass%data(k,y_index,x_index) * c2
        current_state%zqCoarseSolMass%data(k,y_index,x_index)=existing_value - &
        current_state%qCoarseSolMass%data(k,y_index,x_index)

        existing_value = current_state%qCoarseDustMass%data(k,y_index,x_index) + &
        current_state%zqCoarseDustMass%data(k,y_index,x_index)
        current_state%qCoarseDustMass%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%qCoarseDustMass%data(k,y_index,x_index) * c2
        current_state%zqCoarseDustMass%data(k,y_index,x_index)=existing_value - &
        current_state%qCoarseDustMass%data(k,y_index,x_index)


        existing_value = current_state%nl%data(k,y_index,x_index) + current_state%znl%data(k,y_index,x_index)
        current_state%nl%data(k,y_index,x_index)=existing_value * c1 + current_state%nl%data(k,y_index,x_index) * c2
        current_state%znl%data(k,y_index,x_index)=existing_value - current_state%nl%data(k,y_index,x_index)

        existing_value = current_state%nr%data(k,y_index,x_index) + current_state%znr%data(k,y_index,x_index)
        current_state%nr%data(k,y_index,x_index)=existing_value * c1 + current_state%nr%data(k,y_index,x_index) * c2
        current_state%znr%data(k,y_index,x_index)=existing_value - current_state%nr%data(k,y_index,x_index)

        existing_value = current_state%ni%data(k,y_index,x_index) + current_state%zni%data(k,y_index,x_index)
        current_state%ni%data(k,y_index,x_index)=existing_value * c1 + current_state%ni%data(k,y_index,x_index) * c2
        current_state%zni%data(k,y_index,x_index)=existing_value - current_state%ni%data(k,y_index,x_index)

        existing_value = current_state%ns%data(k,y_index,x_index) + current_state%zns%data(k,y_index,x_index)
        current_state%ns%data(k,y_index,x_index)=existing_value * c1 + current_state%ns%data(k,y_index,x_index) * c2
        current_state%zns%data(k,y_index,x_index)=existing_value - current_state%ns%data(k,y_index,x_index)

        existing_value = current_state%ng%data(k,y_index,x_index) + current_state%zng%data(k,y_index,x_index)
        current_state%ng%data(k,y_index,x_index)=existing_value * c1 + current_state%ng%data(k,y_index,x_index) * c2
        current_state%zng%data(k,y_index,x_index)=existing_value - current_state%ng%data(k,y_index,x_index)

        existing_value = current_state%nAitkenSolNumber%data(k,y_index,x_index) + &
        current_state%znAitkenSolNumber%data(k,y_index,x_index)
        current_state%nAitkenSolNumber%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%nAitkenSolNumber%data(k,y_index,x_index) * c2
        current_state%znAitkenSolNumber%data(k,y_index,x_index)=existing_value - &
        current_state%nAitkenSolNumber%data(k,y_index,x_index)

        existing_value = current_state%nAccumSolNumber%data(k,y_index,x_index) + &
        current_state%znAccumSolNumber%data(k,y_index,x_index)
        current_state%nAccumSolNumber%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%nAccumSolNumber%data(k,y_index,x_index) * c2
        current_state%znAccumSolNumber%data(k,y_index,x_index)=existing_value - &
        current_state%nAccumSolNumber%data(k,y_index,x_index)

        existing_value = current_state%nAccumInsolNumber%data(k,y_index,x_index) + &
        current_state%znAccumInsolNumber%data(k,y_index,x_index)
        current_state%nAccumInsolNumber%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%nAccumInsolNumber%data(k,y_index,x_index) * c2
        current_state%znAccumInsolNumber%data(k,y_index,x_index)=existing_value - &
        current_state%nAccumInsolNumber%data(k,y_index,x_index)

        existing_value = current_state%nCoarseSolNumber%data(k,y_index,x_index) + &
        current_state%znCoarseSolNumber%data(k,y_index,x_index)
        current_state%nCoarseSolNumber%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%nCoarseSolNumber%data(k,y_index,x_index) * c2
        current_state%znCoarseSolNumber%data(k,y_index,x_index)=existing_value - &
        current_state%nCoarseSolNumber%data(k,y_index,x_index)

        existing_value = current_state%nCoarseDustnumber%data(k,y_index,x_index) + &
        current_state%znCoarseDustnumber%data(k,y_index,x_index)
        current_state%nCoarseDustnumber%data(k,y_index,x_index)=existing_value *c1 +&
        current_state%nCoarseDustnumber%data(k,y_index,x_index) * c2
        current_state%znCoarseDustnumber%data(k,y_index,x_index)=existing_value - &
        current_state%nCoarseDustnumber%data(k,y_index,x_index)
      end if
    end do
  end subroutine swap_and_smooth_classic

  !> Swap and smooth with a Robert filter
  !! @param current_state The current model state
  subroutine swap_and_smooth_robert_filter(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: y_index, x_index, k, n
    real(kind=DEFAULT_PRECISION) :: c1, c2, existing_value

    x_index=current_state%column_local_x
    y_index=current_state%column_local_y

    c1 = 1.0_DEFAULT_PRECISION - 2.0_DEFAULT_PRECISION*current_state%tsmth
    c2 = current_state%tsmth

    do k=1,current_state%global_grid%size(Z_INDEX)
      existing_value = current_state%u%data(k,y_index,x_index)
      current_state%u%data(k,y_index,x_index)=current_state%zu%data(k,y_index,x_index)
      current_state%zu%data(k,y_index,x_index)=c1*existing_value+c2*(current_state%u%data(k, y_index, x_index)+&
           current_state%savu%data(k,y_index,x_index) -current_state%ugal)
      existing_value = current_state%v%data(k,y_index,x_index)
      current_state%v%data(k,y_index,x_index)=current_state%zv%data(k,y_index,x_index)
      current_state%zv%data(k,y_index,x_index)=c1*existing_value+c2*(current_state%v%data(k, y_index, x_index)+&
           current_state%savv%data(k,y_index,x_index)-current_state%vgal)
      existing_value = current_state%w%data(k,y_index,x_index)
      current_state%w%data(k,y_index,x_index)=current_state%zw%data(k,y_index,x_index)
      current_state%zw%data(k,y_index,x_index)=c1*existing_value+c2*(current_state%w%data(k, y_index, x_index)+&
           current_state%savw%data(k,y_index,x_index))
      if (current_state%th%active) then
        ! Uses the partial smooth of theta from stepfields
        existing_value = current_state%zth%data(k,y_index,x_index)
        current_state%zth%data(k,y_index,x_index)=current_state%th%data(k,y_index,x_index) + current_state%tsmth * existing_value
        current_state%th%data(k,y_index,x_index)=existing_value
        current_state%Tk%data(k,y_index, x_index) = (current_state%th%data(k,y_index, x_index)       &
            + current_state%global_grid%configuration%vertical%thref(k)) * &
            current_state%global_grid%configuration%vertical%rprefrcp(k)
        current_state%qv_saturation%data(k,y_index, x_index) = qsaturation(current_state%Tk%data(k,y_index, x_index), &
                            current_state%global_grid%configuration%vertical%prefn(k)/100.0)
        current_state%qi_saturation%data(k,y_index, x_index)= qisaturation(current_state%Tk%data(k,y_index, x_index), &
                            current_state%global_grid%configuration%vertical%prefn(k)/100.0)
      end if
      if (current_state%number_q_fields .gt. 0) then
        existing_value = current_state%zqv%data(k,y_index,x_index)
        current_state%zqv%data(k,y_index,x_index)=current_state%qv%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qv%data(k,y_index,x_index)=existing_value

        current_state%RH%data(k,y_index, x_index) = current_state%qv%data(k,y_index,x_index) &
                                                    *100/current_state%qv_saturation%data(k,y_index, x_index)

        existing_value = current_state%zql%data(k,y_index,x_index)
        current_state%zql%data(k,y_index,x_index)=current_state%ql%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%ql%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zqr%data(k,y_index,x_index)
        current_state%zqr%data(k,y_index,x_index)=current_state%qr%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qr%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zqi%data(k,y_index,x_index)
        current_state%zqi%data(k,y_index,x_index)=current_state%qi%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qi%data(k,y_index,x_index)=existing_value

        current_state%RI%data(k,y_index, x_index) = current_state%qv%data(k,y_index,x_index) &
                                                    *100/current_state%qi_saturation%data(k,y_index, x_index)

        existing_value = current_state%zqs%data(k,y_index,x_index)
        current_state%zqs%data(k,y_index,x_index)=current_state%qs%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qs%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zqg%data(k,y_index,x_index)
        current_state%zqg%data(k,y_index,x_index)=current_state%qg%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qg%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zqAitkenSolMass%data(k,y_index,x_index)
        current_state%zqAitkenSolMass%data(k,y_index,x_index)= &
           current_state%qAitkenSolMass%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qAitkenSolMass%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zqAccumSolMass%data(k,y_index,x_index)
        current_state%zqAccumSolMass%data(k,y_index,x_index)= &
           current_state%qAccumSolMass%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qAccumSolMass%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zqAccumInsolMass%data(k,y_index,x_index)
        current_state%zqAccumInsolMass%data(k,y_index,x_index)= &
           current_state%qAccumInsolMass%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qAccumInsolMass%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zqCoarseSolMass%data(k,y_index,x_index)
        current_state%zqCoarseSolMass%data(k,y_index,x_index)= &
           current_state%qCoarseSolMass%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qCoarseSolMass%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zqCoarseDustMass%data(k,y_index,x_index)
        current_state%zqCoarseDustMass%data(k,y_index,x_index)= &
           current_state%qCoarseDustMass%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%qCoarseDustMass%data(k,y_index,x_index)=existing_value

        existing_value = current_state%znl%data(k,y_index,x_index)
        current_state%znl%data(k,y_index,x_index)=current_state%nl%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%nl%data(k,y_index,x_index)=existing_value

        existing_value = current_state%znr%data(k,y_index,x_index)
        current_state%znr%data(k,y_index,x_index)=current_state%nr%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%nr%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zni%data(k,y_index,x_index)
        current_state%zni%data(k,y_index,x_index)=current_state%ni%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%ni%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zns%data(k,y_index,x_index)
        current_state%zns%data(k,y_index,x_index)=current_state%ns%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%ns%data(k,y_index,x_index)=existing_value

        existing_value = current_state%zng%data(k,y_index,x_index)
        current_state%zng%data(k,y_index,x_index)=current_state%ng%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%ng%data(k,y_index,x_index)=existing_value

        existing_value = current_state%znAitkenSolNumber%data(k,y_index,x_index)
        current_state%znAitkenSolNumber%data(k,y_index,x_index)= &
           current_state%nAitkenSolNumber%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%nAitkenSolNumber%data(k,y_index,x_index)=existing_value

        existing_value = current_state%znAccumSolNumber%data(k,y_index,x_index)
        current_state%znAccumSolNumber%data(k,y_index,x_index)= &
           current_state%nAccumSolNumber%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%nAccumSolNumber%data(k,y_index,x_index)=existing_value

        existing_value = current_state%znAccumInsolNumber%data(k,y_index,x_index)
        current_state%znAccumInsolNumber%data(k,y_index,x_index)= &
           current_state%nAccumInsolNumber%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%nAccumInsolNumber%data(k,y_index,x_index)=existing_value

        existing_value = current_state%znCoarseSolNumber%data(k,y_index,x_index)
        current_state%znCoarseSolNumber%data(k,y_index,x_index)= &
           current_state%nCoarseSolNumber%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%nCoarseSolNumber%data(k,y_index,x_index)=existing_value

        existing_value = current_state%znCoarseDustnumber%data(k,y_index,x_index)
        current_state%znCoarseDustnumber%data(k,y_index,x_index)= &
           current_state%nCoarseDustnumber%data(k,y_index,x_index)+&
              current_state%tsmth * existing_value
        current_state%nCoarseDustnumber%data(k,y_index,x_index)=existing_value
      end if
    end do
  end subroutine swap_and_smooth_robert_filter

  !> Does swapping and smoothing (using classic algorithm) for the average profiles (the bars)
  !! @param current_state The current model state
  subroutine classic_for_average_profiles(current_state, old_smoother)
    type(model_state_type), intent(inout) :: current_state
    logical, intent(in) :: old_smoother

    integer :: k, n
    real(kind=DEFAULT_PRECISION) :: c1, c2

    if (old_smoother) then
      c1 = 1.0_DEFAULT_PRECISION - current_state%tsmth
      c2 = 2.0_DEFAULT_PRECISION * current_state%tsmth - 1.0_DEFAULT_PRECISION
    else
      c1 = 1.0_DEFAULT_PRECISION
      c2 = -1.0_DEFAULT_PRECISION
    end if

    do k=1,current_state%global_grid%size(Z_INDEX)
      current_state%global_grid%configuration%vertical%olzubar(k)=current_state%global_grid%configuration%vertical%olubar(k) +&
           current_state%global_grid%configuration%vertical%olzubar(k)
      current_state%global_grid%configuration%vertical%olubar(k)=current_state%global_grid%configuration%vertical%olzubar(k) *&
           c1 + current_state%global_grid%configuration%vertical%olubar(k) * c2
      current_state%global_grid%configuration%vertical%olzubar(k)=current_state%global_grid%configuration%vertical%olzubar(k) -&
           current_state%global_grid%configuration%vertical%olubar(k)
      current_state%global_grid%configuration%vertical%olzvbar(k)=current_state%global_grid%configuration%vertical%olvbar(k) +&
           current_state%global_grid%configuration%vertical%olzvbar(k)
      current_state%global_grid%configuration%vertical%olvbar(k)=current_state%global_grid%configuration%vertical%olzvbar(k) *&
           c1 + current_state%global_grid%configuration%vertical%olvbar(k) * c2
      current_state%global_grid%configuration%vertical%olzvbar(k)=current_state%global_grid%configuration%vertical%olzvbar(k) -&
           current_state%global_grid%configuration%vertical%olvbar(k)
      if (current_state%th%active) then
        current_state%global_grid%configuration%vertical%olzthbar(k)=current_state%global_grid%configuration%vertical%olthbar(k)+&
             current_state%global_grid%configuration%vertical%olzthbar(k)
        current_state%global_grid%configuration%vertical%olthbar(k)=current_state%global_grid%configuration%vertical%olzthbar(k)*&
             c1 + current_state%global_grid%configuration%vertical%olthbar(k) * c2
        current_state%global_grid%configuration%vertical%olzthbar(k)=&
             current_state%global_grid%configuration%vertical%olzthbar(k)-&
             current_state%global_grid%configuration%vertical%olthbar(k)
      end if
      if (current_state%number_q_fields .gt. 0) then
        current_state%global_grid%configuration%vertical%olzqvbar(k)=&
              current_state%global_grid%configuration%vertical%olqvbar(k)+&
              current_state%global_grid%configuration%vertical%olzqvbar(k)
        current_state%global_grid%configuration%vertical%olqvbar(k)=&
              current_state%global_grid%configuration%vertical%olzqvbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqvbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqvbar(k)=&
              current_state%global_grid%configuration%vertical%olzqvbar(k)-&
              current_state%global_grid%configuration%vertical%olqvbar(k)

        current_state%global_grid%configuration%vertical%olzqlbar(k)=&
              current_state%global_grid%configuration%vertical%olqlbar(k)+&
              current_state%global_grid%configuration%vertical%olzqlbar(k)
        current_state%global_grid%configuration%vertical%olqlbar(k)=&
              current_state%global_grid%configuration%vertical%olzqlbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqlbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqlbar(k)=&
              current_state%global_grid%configuration%vertical%olzqlbar(k)-&
              current_state%global_grid%configuration%vertical%olqlbar(k)

        current_state%global_grid%configuration%vertical%olzqrbar(k)=&
              current_state%global_grid%configuration%vertical%olqrbar(k)+&
              current_state%global_grid%configuration%vertical%olzqrbar(k)
        current_state%global_grid%configuration%vertical%olqrbar(k)=&
              current_state%global_grid%configuration%vertical%olzqrbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqrbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqrbar(k)=&
              current_state%global_grid%configuration%vertical%olzqrbar(k)-&
              current_state%global_grid%configuration%vertical%olqrbar(k)

        current_state%global_grid%configuration%vertical%olzqibar(k)=&
              current_state%global_grid%configuration%vertical%olqibar(k)+&
              current_state%global_grid%configuration%vertical%olzqibar(k)
        current_state%global_grid%configuration%vertical%olqibar(k)=&
              current_state%global_grid%configuration%vertical%olzqibar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqibar(k) * c2
        current_state%global_grid%configuration%vertical%olzqibar(k)=&
              current_state%global_grid%configuration%vertical%olzqibar(k)-&
              current_state%global_grid%configuration%vertical%olqibar(k)

        current_state%global_grid%configuration%vertical%olzqsbar(k)=&
              current_state%global_grid%configuration%vertical%olqsbar(k)+&
              current_state%global_grid%configuration%vertical%olzqsbar(k)
        current_state%global_grid%configuration%vertical%olqsbar(k)=&
              current_state%global_grid%configuration%vertical%olzqsbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqsbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqsbar(k)=&
              current_state%global_grid%configuration%vertical%olzqsbar(k)-&
              current_state%global_grid%configuration%vertical%olqsbar(k)

        current_state%global_grid%configuration%vertical%olzqgbar(k)=&
              current_state%global_grid%configuration%vertical%olqgbar(k)+&
              current_state%global_grid%configuration%vertical%olzqgbar(k)
        current_state%global_grid%configuration%vertical%olqgbar(k)=&
              current_state%global_grid%configuration%vertical%olzqgbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqgbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqgbar(k)=&
              current_state%global_grid%configuration%vertical%olzqgbar(k)-&
              current_state%global_grid%configuration%vertical%olqgbar(k)

        current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k)+&
              current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(k)
        current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(k)-&
              current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k)

        current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k)+&
              current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(k)
        current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(k)-&
              current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k)

        current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k)+&
              current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(k)
        current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(k)-&
              current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k)

        current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k)+&
              current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(k)
        current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(k)-&
              current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k)

        current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(k)=&
              current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k)+&
              current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(k)
        current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k) * c2
        current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(k)=&
              current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(k)-&
              current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k)


        current_state%global_grid%configuration%vertical%olznlbar(k)=&
              current_state%global_grid%configuration%vertical%olnlbar(k)+&
              current_state%global_grid%configuration%vertical%olznlbar(k)
        current_state%global_grid%configuration%vertical%olnlbar(k)=&
              current_state%global_grid%configuration%vertical%olznlbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olnlbar(k) * c2
        current_state%global_grid%configuration%vertical%olznlbar(k)=&
              current_state%global_grid%configuration%vertical%olznlbar(k)-&
              current_state%global_grid%configuration%vertical%olnlbar(k)

        current_state%global_grid%configuration%vertical%olznrbar(k)=&
              current_state%global_grid%configuration%vertical%olnrbar(k)+&
              current_state%global_grid%configuration%vertical%olznrbar(k)
        current_state%global_grid%configuration%vertical%olnrbar(k)=&
              current_state%global_grid%configuration%vertical%olznrbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olnrbar(k) * c2
        current_state%global_grid%configuration%vertical%olznrbar(k)=&
              current_state%global_grid%configuration%vertical%olznrbar(k)-&
              current_state%global_grid%configuration%vertical%olnrbar(k)

        current_state%global_grid%configuration%vertical%olznibar(k)=&
              current_state%global_grid%configuration%vertical%olnibar(k)+&
              current_state%global_grid%configuration%vertical%olznibar(k)
        current_state%global_grid%configuration%vertical%olnibar(k)=&
              current_state%global_grid%configuration%vertical%olznibar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olnibar(k) * c2
        current_state%global_grid%configuration%vertical%olznibar(k)=&
              current_state%global_grid%configuration%vertical%olznibar(k)-&
              current_state%global_grid%configuration%vertical%olnibar(k)

        current_state%global_grid%configuration%vertical%olznsbar(k)=&
              current_state%global_grid%configuration%vertical%olnsbar(k)+&
              current_state%global_grid%configuration%vertical%olznsbar(k)
        current_state%global_grid%configuration%vertical%olnsbar(k)=&
              current_state%global_grid%configuration%vertical%olznsbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olnsbar(k) * c2
        current_state%global_grid%configuration%vertical%olznsbar(k)=&
              current_state%global_grid%configuration%vertical%olznsbar(k)-&
              current_state%global_grid%configuration%vertical%olnsbar(k)

        current_state%global_grid%configuration%vertical%olzngbar(k)=&
              current_state%global_grid%configuration%vertical%olngbar(k)+&
              current_state%global_grid%configuration%vertical%olzngbar(k)
        current_state%global_grid%configuration%vertical%olngbar(k)=&
              current_state%global_grid%configuration%vertical%olzngbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olngbar(k) * c2
        current_state%global_grid%configuration%vertical%olzngbar(k)=&
              current_state%global_grid%configuration%vertical%olzngbar(k)-&
              current_state%global_grid%configuration%vertical%olngbar(k)

        current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k)+&
              current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(k)
        current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k) * c2
        current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(k)-&
              current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k)

        current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k)+&
              current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(k)
        current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k) * c2
        current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(k)-&
              current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k)

        current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k)+&
              current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(k)
        current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k) * c2
        current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(k)-&
              current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k)

        current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k)+&
              current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(k)
        current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k) * c2
        current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(k)-&
              current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k)

        current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(k)=&
              current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k)+&
              current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(k)
        current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(k)*&
              c1 + current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k) * c2
        current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(k)=&
              current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(k)-&
              current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k)
      end if
    end do
  end subroutine classic_for_average_profiles

  !> Does swapping and smoothing (using robert filter algorithm) for the average profiles (the bars)
  !! @param current_state The current model state
  subroutine robert_filter_for_average_profiles(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k, n
    real(kind=DEFAULT_PRECISION) :: c1, c2, existing_value

    c1 = 1.0_DEFAULT_PRECISION - 2.0_DEFAULT_PRECISION*current_state%tsmth
    c2 = current_state%tsmth

    do k=1,current_state%global_grid%size(Z_INDEX)
      existing_value=current_state%global_grid%configuration%vertical%olubar(k)
      current_state%global_grid%configuration%vertical%olubar(k)=current_state%global_grid%configuration%vertical%olzubar(k)
      current_state%global_grid%configuration%vertical%olzubar(k)=c1*existing_value+c2*&
           (current_state%global_grid%configuration%vertical%olubar(k) + &
           current_state%global_grid%configuration%vertical%savolubar(k))
      existing_value=current_state%global_grid%configuration%vertical%olvbar(k)
      current_state%global_grid%configuration%vertical%olvbar(k)=current_state%global_grid%configuration%vertical%olzvbar(k)
      current_state%global_grid%configuration%vertical%olzvbar(k)=c1*existing_value+c2*&
           (current_state%global_grid%configuration%vertical%olvbar(k) + &
           current_state%global_grid%configuration%vertical%savolvbar(k))
      if (current_state%th%active) then
        existing_value=current_state%global_grid%configuration%vertical%olzthbar(k)
        current_state%global_grid%configuration%vertical%olzthbar(k)=&
             current_state%global_grid%configuration%vertical%olthbar(k) + current_state%tsmth * existing_value
        current_state%global_grid%configuration%vertical%olthbar(k)=existing_value
      end if
      if (current_state%number_q_fields .gt. 0) then
        existing_value=current_state%global_grid%configuration%vertical%olzqvbar(k)
          current_state%global_grid%configuration%vertical%olzqvbar(k)=&
               current_state%global_grid%configuration%vertical%olqvbar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqvbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqlbar(k)
          current_state%global_grid%configuration%vertical%olzqlbar(k)=&
               current_state%global_grid%configuration%vertical%olqlbar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqlbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqrbar(k)
          current_state%global_grid%configuration%vertical%olzqrbar(k)=&
               current_state%global_grid%configuration%vertical%olqrbar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqrbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqibar(k)
          current_state%global_grid%configuration%vertical%olzqibar(k)=&
               current_state%global_grid%configuration%vertical%olqibar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqibar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqsbar(k)
          current_state%global_grid%configuration%vertical%olzqsbar(k)=&
               current_state%global_grid%configuration%vertical%olqsbar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqsbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqgbar(k)
          current_state%global_grid%configuration%vertical%olzqgbar(k)=&
               current_state%global_grid%configuration%vertical%olqgbar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqgbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(k)
          current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(k)=&
               current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(k)
          current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(k)=&
               current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(k)
          current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(k)=&
               current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(k)
          current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(k)=&
               current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(k)
          current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(k)=&
               current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k)=existing_value


        existing_value=current_state%global_grid%configuration%vertical%olznlbar(k)
          current_state%global_grid%configuration%vertical%olznlbar(k)=&
               current_state%global_grid%configuration%vertical%olnlbar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olnlbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olznrbar(k)
          current_state%global_grid%configuration%vertical%olznrbar(k)=&
               current_state%global_grid%configuration%vertical%olnrbar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olnrbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olznibar(k)
          current_state%global_grid%configuration%vertical%olznibar(k)=&
               current_state%global_grid%configuration%vertical%olnibar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olnibar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olznsbar(k)
          current_state%global_grid%configuration%vertical%olznsbar(k)=&
               current_state%global_grid%configuration%vertical%olnsbar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olnsbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olzngbar(k)
          current_state%global_grid%configuration%vertical%olzngbar(k)=&
               current_state%global_grid%configuration%vertical%olngbar(k) + current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olngbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(k)
          current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(k)=&
               current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(k)
          current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(k)=&
               current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(k)
          current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(k)=&
               current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(k)
          current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(k)=&
               current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k)=existing_value
        existing_value=current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(k)
          current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(k)=&
               current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k) + &
                        current_state%tsmth * existing_value
          current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k)=existing_value
      end if
    end do
  end subroutine robert_filter_for_average_profiles
end module swapsmooth_mod
