!> Clears the source terms at the start of a timestep to then be populated by items in the dynamics group
module clearsourceterms_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use state_mod, only : model_state_type
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type, component_descriptor_type_v1
implicit none

#ifndef TEST_MODE
  private
#endif

  public timestep_callback_clearsourceterms, no_component
contains
  subroutine no_component(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    !print*,"no_component"
    continue
  end subroutine no_component

  !> Timestep callback which simply clears the source terms
  !! @param current_state The current model state
  subroutine timestep_callback_clearsourceterms(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: i

#ifdef U_ACTIVE
    current_state%su%data=0.0_DEFAULT_PRECISION
#endif
#ifdef V_ACTIVE
    current_state%sv%data=0.0_DEFAULT_PRECISION
#endif
#ifdef W_ACTIVE
    current_state%sw%data=0.0_DEFAULT_PRECISION
#endif

    if (current_state%sth%active) then
      current_state%sth%data=0.0_DEFAULT_PRECISION
    end if

    if (current_state%number_q_fields .gt. 0) then
      current_state%sqv%data=0.0_DEFAULT_PRECISION
      current_state%sql%data=0.0_DEFAULT_PRECISION
      current_state%sqr%data=0.0_DEFAULT_PRECISION
      current_state%sqi%data=0.0_DEFAULT_PRECISION
      current_state%sqs%data=0.0_DEFAULT_PRECISION
      current_state%sqg%data=0.0_DEFAULT_PRECISION
      current_state%snl%data=0.0_DEFAULT_PRECISION
      current_state%snr%data=0.0_DEFAULT_PRECISION
      current_state%sni%data=0.0_DEFAULT_PRECISION
      current_state%sns%data=0.0_DEFAULT_PRECISION
      current_state%sng%data=0.0_DEFAULT_PRECISION
      current_state%sqAitkenSolMass%data=0.0_DEFAULT_PRECISION
      current_state%sqAccumSolMass%data=0.0_DEFAULT_PRECISION
      current_state%sqAccumInsolMass%data=0.0_DEFAULT_PRECISION
      current_state%sqCoarseSolMass%data=0.0_DEFAULT_PRECISION
      current_state%sqCoarseDustMass%data=0.0_DEFAULT_PRECISION
      current_state%snAitkenSolNumber%data=0.0_DEFAULT_PRECISION
      current_state%snAccumSolNumber%data=0.0_DEFAULT_PRECISION
      current_state%snAccumInsolNumber%data=0.0_DEFAULT_PRECISION
      current_state%snCoarseSolNumber%data=0.0_DEFAULT_PRECISION
      current_state%snCoarseDustnumber%data=0.0_DEFAULT_PRECISION
    end if
    !do i=1, current_state%number_q_fields
    !  current_state%sq(i)%data=0.0_DEFAULT_PRECISION
    !end do

    if (current_state%n_tracers .gt. 0) then
      do i=1, current_state%n_tracers
        current_state%stracer(i)%data=0.0_DEFAULT_PRECISION
      end do
      if (current_state%traj_tracer_index >0 .and. current_state%reinit_tracer) then
        current_state%reinit_tracer=.false.
      end if
    end if
  end subroutine timestep_callback_clearsourceterms
end module clearsourceterms_mod
