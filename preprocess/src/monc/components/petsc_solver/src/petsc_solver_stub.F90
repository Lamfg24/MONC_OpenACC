










!> Dummy stub when not compiling with PETSc iterative solver
module petsc_solver_mod
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  implicit none

  private

  public petsc_solver_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type_v1) function petsc_solver_get_descriptor()
    petsc_solver_get_descriptor%name="petsc_solver"
    petsc_solver_get_descriptor%version=0.0
    petsc_solver_get_descriptor%published_fields_on_off = .false.
  end function petsc_solver_get_descriptor
end module petsc_solver_mod
