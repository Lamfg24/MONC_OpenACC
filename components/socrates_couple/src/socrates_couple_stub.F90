!> Dummy stub when not compiling in socrates
module socrates_couple_mod
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  implicit none

  public socrates_couple_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type_v1) function socrates_couple_get_descriptor()
    socrates_couple_get_descriptor%name="socrates_couple"
    socrates_couple_get_descriptor%version=0.1
    socrates_couple_get_descriptor%published_fields_on_off = .false.

  end function socrates_couple_get_descriptor
    
end module socrates_couple_mod


    
    

    
