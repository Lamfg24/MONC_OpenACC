!> Dummy stub when not compiling with CASIM microphysics
module casim_mod
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  implicit none

#ifndef TEST_MODE
  private
#endif

  public casim_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type_v1) function casim_get_descriptor()
    casim_get_descriptor%name="casim"
    casim_get_descriptor%version=0.0
    casim_get_descriptor%published_fields_on_off = .false.
  end function casim_get_descriptor

end module casim_mod
