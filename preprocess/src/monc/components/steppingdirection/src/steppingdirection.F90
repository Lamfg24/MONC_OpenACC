










!> Determines the current stepping direction, which can be either forward or centred. This is mainly for field stepping,
!! which is u, v, w fields but also scalars as well which is th and q.
module steppingdirection_mod
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use state_mod, only : model_state_type, CENTRED_STEPPING, FORWARD_STEPPING
  implicit none

  private

  public initialisation_callback_steppingdirection, timestep_callback_steppingdirection

contains

  !> Sets the scalar stepping on initialisation. This does not change throughout the model run so we can safely set it here
  !! @param current_state The current model state
  subroutine initialisation_callback_steppingdirection(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    ! TODO - allow configuration to use forward stepping for scalars (th and q) and momentum (u,v,w)
    current_state%scalar_stepping = CENTRED_STEPPING
    current_state%momentum_stepping = CENTRED_STEPPING
  end subroutine initialisation_callback_steppingdirection

  !> Determines whether we are forward or centre stepping
  !!
  !! This is important as forward stepping will effectively ignore the Z field
  !! and performs no smoothing on the field. Centre stepping uses the Z field and
  !! the Robert filter for smoothing.
  !! @param current_state The current model state
  subroutine timestep_callback_steppingdirection(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    
    current_state%field_stepping = CENTRED_STEPPING
    if (current_state%timestep .eq. current_state%start_timestep) then
      current_state%field_stepping = FORWARD_STEPPING
      if (current_state%reconfig_run) current_state%field_stepping = CENTRED_STEPPING
    end if
  end subroutine timestep_callback_steppingdirection
end module steppingdirection_mod
