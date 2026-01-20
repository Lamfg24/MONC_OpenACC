










!> Calculates the mean profiles of prognostic variables which are then used in smoothing and other areas
module meanprofiles_mod
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use mpi, only : MPI_SUM, MPI_IN_PLACE
  use netcdf, only : nf90_noerr, nf90_global, nf90_nowrite, nf90_inquire_attribute, nf90_open, nf90_strerror, &
       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_inquire, nf90_close, nf90_get_att
  implicit none

  private

  integer :: start_x, end_x, start_y, end_y, bar_fields

  real(kind=DEFAULT_PRECISION) :: rnhpts
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: bartmp
  ! Immersed boundary variables
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: rnspts_u, rnspts_v, rnspts_s

  public init_callback_meanprofiles, timestep_callback_meanprofiles, &
         finalisation_callback_meanprofiles
contains

  !> Called on MONC initialisation, will allocate appropriate data structures
  !! @param current_state The current model state
  subroutine init_callback_meanprofiles(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION) :: tmp
    integer :: k, ierr, myproc
    integer :: ncdf_id, variable_id


    myproc = current_state%parallel%my_rank

    if (current_state%continuation_run .eqv. .true.) then
      call check(nf90_open(current_state%options_database_string(2,2), nf90_nowrite, ncdf_id))
    end if

    bar_fields=0

    rnhpts=1.0_DEFAULT_PRECISION/real(current_state%global_grid%size(X_INDEX)*current_state%global_grid%size(Y_INDEX))

    start_x=current_state%local_grid%local_domain_start_index(X_INDEX)
    end_x=current_state%local_grid%local_domain_end_index(X_INDEX)
    start_y=current_state%local_grid%local_domain_start_index(Y_INDEX)
    end_y=current_state%local_grid%local_domain_end_index(Y_INDEX)

    if (.not. current_state%continuation_run) then
      allocate(current_state%global_grid%configuration%vertical%olubar(current_state%local_grid%size(Z_INDEX)),&
           current_state%global_grid%configuration%vertical%olzubar(current_state%local_grid%size(Z_INDEX)))
    else
      allocate(current_state%global_grid%configuration%vertical%olubar(current_state%local_grid%size(Z_INDEX)),&
           current_state%global_grid%configuration%vertical%olzubar(current_state%local_grid%size(Z_INDEX)))
      call check(nf90_inq_varid(ncdf_id, "olubar", variable_id))
      call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olubar( &
          current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
          start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      call check(nf90_inq_varid(ncdf_id, "olzubar", variable_id))
      call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzubar( &
          current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
          start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
    end if
    bar_fields=bar_fields+2
    if (.not. current_state%continuation_run) then
      allocate(current_state%global_grid%configuration%vertical%olvbar(current_state%local_grid%size(Z_INDEX)),&
           current_state%global_grid%configuration%vertical%olzvbar(current_state%local_grid%size(Z_INDEX)))
    else
      allocate(current_state%global_grid%configuration%vertical%olvbar(current_state%local_grid%size(Z_INDEX)),&
           current_state%global_grid%configuration%vertical%olzvbar(current_state%local_grid%size(Z_INDEX)))
      call check(nf90_inq_varid(ncdf_id, "olvbar", variable_id))
      call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olvbar( &
          current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
          start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      call check(nf90_inq_varid(ncdf_id, "olzvbar", variable_id))
      call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzvbar( &
          current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
          start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
    end if
    bar_fields=bar_fields+2
    if (current_state%th%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olthbar(current_state%local_grid%size(Z_INDEX)),&
             current_state%global_grid%configuration%vertical%olzthbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olthbar(current_state%local_grid%size(Z_INDEX)),&
             current_state%global_grid%configuration%vertical%olzthbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olthbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olthbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzthbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzthbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if

    !! MASS
    if (current_state%qv%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqvbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqvbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqvbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqvbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqvbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olqvbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqvbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzqvbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%ql%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqlbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqlbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqlbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqlbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqlbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olqlbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqlbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzqlbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%qr%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqrbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqrbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqrbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqrbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqrbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olqrbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqrbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzqrbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%qi%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqibar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqibar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqibar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqibar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqibar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olqibar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqibar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzqibar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%qs%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqsbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqsbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqsbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqsbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqsbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olqsbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqsbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzqsbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%qg%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqgbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqgbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqgbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqgbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqgbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olqgbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqgbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzqgbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%qAitkenSolMass%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqAitkenSolMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olqAitkenSolMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqAitkenSolMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%qAccumSolMass%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqAccumSolMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqAccumSolMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqAccumSolMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olqAccumSolMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqAccumSolMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olzqAccumSolMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%qAccumInsolMass%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqAccumInsolMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olqAccumInsolMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqAccumInsolMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%qCoarseSolMass%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqCoarseSolMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olqCoarseSolMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqCoarseSolMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%qCoarseDustMass%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olqCoarseDustMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olqCoarseDustMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzqCoarseDustMassbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if

    !! NUMBER
    if (current_state%nl%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olnlbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznlbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olnlbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznlbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olnlbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olnlbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olznlbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olznlbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%nr%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olnrbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznrbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olnrbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznrbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olnrbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olnrbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olznrbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olznrbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%ni%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olnibar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznibar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olnibar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznibar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olnibar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olnibar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olznibar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olznibar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%ns%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olnsbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznsbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olnsbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznsbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olnsbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olnsbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olznsbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olznsbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%ng%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olngbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzngbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olngbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olzngbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olngbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olngbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olzngbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, current_state%global_grid%configuration%vertical%olzngbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%nAitkenSolNumber%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olnAitkenSolNumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olznAitkenSolNumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%nAccumSolNumber%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olnAccumSolNumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olnAccumSolNumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olznAccumSolNumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olznAccumSolNumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%nAccumInsolNumber%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olnAccumInsolNumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olznAccumInsolNumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%nCoarseSolNumber%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olnCoarseSolNumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olznCoarseSolNumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%nCoarseDustnumber%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(current_state%local_grid%size(Z_INDEX)))
      else
        allocate(current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(current_state%local_grid%size(Z_INDEX)), &
                 current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(current_state%local_grid%size(Z_INDEX)))
        call check(nf90_inq_varid(ncdf_id, "olnCoarseDustnumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
        call check(nf90_inq_varid(ncdf_id, "olznCoarseDustnumberbar", variable_id))
        call check(nf90_get_var(ncdf_id, variable_id, &
            current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar( &
            current_state%local_grid%local_domain_start_index(Z_INDEX):current_state%local_grid%local_domain_end_index(Z_INDEX)),&
            start=[1, current_state%local_grid%start(Z_INDEX)], count=[1, current_state%local_grid%size(Z_INDEX)]))
      end if
      bar_fields=bar_fields+2
    end if
    allocate(bartmp(current_state%local_grid%size(Z_INDEX), bar_fields))


    ! Immersed boundary
    !if(current_state%immersed%ib_enabled) then ! à revoir

    !  tmp = real(current_state%local_grid%size(X_INDEX)*current_state%local_grid%size(Y_INDEX))

    !  allocate(rnspts_u(current_state%local_grid%size(Z_INDEX)))
    !  rnspts_u = 1.0_DEFAULT_PRECISION
    !  do k=1,current_state%local_grid%size(Z_INDEX)
    !    rnspts_u(k) = tmp-sum(current_state%immersed%indic_u(k, start_y:end_y, start_x:end_x))
    !  end do
    !  call mpi_allreduce(MPI_IN_PLACE, rnspts_u, current_state%local_grid%size(Z_INDEX), PRECISION_TYPE, MPI_SUM, &
    !                     current_state%parallel%monc_communicator, ierr)
    !  do k=1,current_state%local_grid%size(Z_INDEX)
    !    if(rnspts_u(k).gt.0.0)rnspts_u(k) = 1.0_DEFAULT_PRECISION/rnspts_u(k)
    !  end do

    !  allocate(rnspts_v(current_state%local_grid%size(Z_INDEX)))
    !  rnspts_v = 1.0_DEFAULT_PRECISION
    !  do k=1,current_state%local_grid%size(Z_INDEX)
    !    rnspts_v(k) = tmp-sum(current_state%immersed%indic_v(k, start_y:end_y, start_x:end_x))
    !  end do
    !  call mpi_allreduce(MPI_IN_PLACE, rnspts_v, current_state%local_grid%size(Z_INDEX), PRECISION_TYPE, MPI_SUM, &
    !                     current_state%parallel%monc_communicator, ierr)
    !  do k=1,current_state%local_grid%size(Z_INDEX)
    !    if(rnspts_v(k).gt.0.0)rnspts_v(k) = 1.0_DEFAULT_PRECISION/rnspts_v(k)
    !  end do

    !  if(current_state%th%active.or.(current_state%number_q_fields.gt.0)) then
    !    allocate(rnspts_s(current_state%local_grid%size(Z_INDEX)))
    !    rnspts_s = 1.0_DEFAULT_PRECISION
    !    do k=1,current_state%local_grid%size(Z_INDEX)
    !      rnspts_s(k) = tmp-sum(current_state%immersed%indic_s(k, start_y:end_y, start_x:end_x))
    !    end do
    !    call mpi_allreduce(MPI_IN_PLACE, rnspts_s, current_state%local_grid%size(Z_INDEX), PRECISION_TYPE, MPI_SUM, &
    !                       current_state%parallel%monc_communicator, ierr)
    !    do k=1,current_state%local_grid%size(Z_INDEX)
    !      if(rnspts_s(k).gt.0.0)rnspts_s(k) = 1.0_DEFAULT_PRECISION/rnspts_s(k)
    !    end do
    !  end if

    !end if


    ! Do the initial calculation for the first timestep
    if (.not. current_state%continuation_run) call calculate_mean_profiles(current_state)
  end subroutine init_callback_meanprofiles

  !> Will recalculate the mean profiles of each prognostic when called (for the entire local domain)
  !! @param current_state The current model state
  subroutine timestep_callback_meanprofiles(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call calculate_mean_profiles(current_state)
    
  end subroutine timestep_callback_meanprofiles

  !> Frees up the temporary data for the bars
  !! @param current_state The current model state
  subroutine finalisation_callback_meanprofiles(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(bartmp)) deallocate(bartmp)    
  end subroutine finalisation_callback_meanprofiles

  !> Calculates the global mean profiles and stores these in the ol bar arrays
  !! @param current_state The current model state
  subroutine calculate_mean_profiles(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: bar_index, i, k
    
    call calculate_sum_profiles(current_state)
    
    if(current_state%immersed%ib_enabled)then ! à revoir

    !  bar_index=1
    !  current_state%global_grid%configuration%vertical%olubar(:)=bartmp(:, bar_index)*rnspts_u(:)
    !  current_state%global_grid%configuration%vertical%olzubar(:)=bartmp(:, bar_index+1)*rnspts_u(:)
    !  bar_index=bar_index+2
    !  current_state%global_grid%configuration%vertical%olvbar(:)=bartmp(:, bar_index)*rnspts_v(:)
    !  current_state%global_grid%configuration%vertical%olzvbar(:)=bartmp(:, bar_index+1)*rnspts_v(:)
    !  bar_index=bar_index+2
    !  if (current_state%th%active) then
    !    current_state%global_grid%configuration%vertical%olthbar(:)=bartmp(:, bar_index)*rnspts_s(:)
    !    current_state%global_grid%configuration%vertical%olzthbar(:)=bartmp(:, bar_index+1)*rnspts_s(:)
    !    bar_index=bar_index+2
    !  end if
    !  do i=1,current_state%number_q_fields
    !    if (current_state%qv%active) then
    !      current_state%global_grid%configuration%vertical%olqvbar(:)=bartmp(:, bar_index)*rnspts_s(:)
    !      current_state%global_grid%configuration%vertical%olzqvbar(:)=bartmp(:, bar_index+1)*rnspts_s(:)
    !      bar_index=bar_index+2
    !    end if
    !  end do

    else

      bar_index=1
      current_state%global_grid%configuration%vertical%olubar(:)=bartmp(:, bar_index)*rnhpts
      current_state%global_grid%configuration%vertical%olzubar(:)=bartmp(:, bar_index+1)*rnhpts
      bar_index=bar_index+2
      current_state%global_grid%configuration%vertical%olvbar(:)=bartmp(:, bar_index)*rnhpts
      current_state%global_grid%configuration%vertical%olzvbar(:)=bartmp(:, bar_index+1)*rnhpts
      bar_index=bar_index+2
      if (current_state%th%active) then
        current_state%global_grid%configuration%vertical%olthbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzthbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      !! MASS
      if (current_state%qv%active) then
        current_state%global_grid%configuration%vertical%olqvbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqvbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%ql%active) then
        current_state%global_grid%configuration%vertical%olqlbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqlbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%qr%active) then
        current_state%global_grid%configuration%vertical%olqrbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqrbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%qi%active) then
        current_state%global_grid%configuration%vertical%olqibar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqibar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%qs%active) then
        current_state%global_grid%configuration%vertical%olqsbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqsbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%qg%active) then
        current_state%global_grid%configuration%vertical%olqgbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqgbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%qAitkenSolMass%active) then
        current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqAitkenSolMassbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%qAccumSolMass%active) then
        current_state%global_grid%configuration%vertical%olqAccumSolMassbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqAccumSolMassbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%qAccumInsolMass%active) then
        current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqAccumInsolMassbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%qCoarseSolMass%active) then
        current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqCoarseSolMassbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%qCoarseDustMass%active) then
        current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzqCoarseDustMassbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      !! NUMBER
      if (current_state%nl%active) then
        current_state%global_grid%configuration%vertical%olnlbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olznlbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%nr%active) then
        current_state%global_grid%configuration%vertical%olnrbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olznrbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%ni%active) then
        current_state%global_grid%configuration%vertical%olnibar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olznibar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%ns%active) then
        current_state%global_grid%configuration%vertical%olnsbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olznsbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%ng%active) then
        current_state%global_grid%configuration%vertical%olngbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzngbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%nAitkenSolNumber%active) then
        current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olznAitkenSolNumberbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%nAccumSolNumber%active) then
        current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%nAccumInsolNumber%active) then
        current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%nCoarseSolNumber%active) then
        current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      if (current_state%nCoarseDustnumber%active) then
        current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if

      !do i=1,current_state%number_q_fields
      !  if (current_state%qv%active) then
      !    current_state%global_grid%configuration%vertical%olqvbar(:)=bartmp(:, bar_index)*rnhpts
      !    current_state%global_grid%configuration%vertical%olzqvbar(:)=bartmp(:, bar_index+1)*rnhpts
      !    bar_index=bar_index+2
      !  end if
      !end do
    end if

  end subroutine calculate_mean_profiles    

  !> Calculates the sum profiles for the bars for each level globally
  !! @param current_state The current model state_mod
  subroutine calculate_sum_profiles(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k, n, bar_index, ierr    

    if(current_state%immersed%ib_enabled)then
      !do k=current_state%local_grid%local_domain_start_index(Z_INDEX), current_state%local_grid%local_domain_end_index(Z_INDEX)
      !  bar_index=1
      !  bartmp(k, bar_index)=sum(current_state%u%data(k, start_y:end_y, start_x:end_x),&
      !                           mask=current_state%immersed%indic_u(k, start_y:end_y,start_x:end_x).eq.0 )
      !  bartmp(k, bar_index+1)=sum(current_state%zu%data(k, start_y:end_y, start_x:end_x),&
      !                             mask=current_state%immersed%indic_u(k, start_y:end_y,start_x:end_x).eq.0 )
      !  bar_index=bar_index+2
      !  bartmp(k, bar_index)=sum(current_state%v%data(k, start_y:end_y, start_x:end_x),&
      !                           mask=current_state%immersed%indic_v(k, start_y:end_y,start_x:end_x).eq.0 )
      !  bartmp(k, bar_index+1)=sum(current_state%zv%data(k, start_y:end_y, start_x:end_x),&
      !                           mask=current_state%immersed%indic_v(k, start_y:end_y,start_x:end_x).eq.0 )
      !  bar_index=bar_index+2
      !  if (current_state%th%active) then
      !    bartmp(k, bar_index)=sum(current_state%th%data(k, start_y:end_y, start_x:end_x),&
      !                             mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
      !    bartmp(k, bar_index+1)=sum(current_state%zth%data(k, start_y:end_y, start_x:end_x),&
      !                               mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
      !    bar_index=bar_index+2
      !  end if
      !  if (current_state%qv%active) then
      !    bartmp(k, bar_index)=sum(current_state%qv%data(k, start_y:end_y, start_x:end_x),&
      !                             mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
      !    bartmp(k, bar_index+1)=sum(current_state%zqv%data(k, start_y:end_y, start_x:end_x),&
      !                               mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
      !    bar_index=bar_index+2
      !  end if
        !do n=1,current_state%number_q_fields
        !  if (current_state%q(n)%active) then
        !    bartmp(k, bar_index)=sum(current_state%q(n)%data(k, start_y:end_y, start_x:end_x),&
        !                             mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
        !    bartmp(k, bar_index+1)=sum(current_state%zq(n)%data(k, start_y:end_y, start_x:end_x),&
        !                               mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
        !    bar_index=bar_index+2
        !  end if
        !end do
      !end do

    else

     do k=current_state%local_grid%local_domain_start_index(Z_INDEX), current_state%local_grid%local_domain_end_index(Z_INDEX)
        bar_index=1
        bartmp(k, bar_index)=sum(current_state%u%data(k, start_y:end_y, start_x:end_x))
        bartmp(k, bar_index+1)=sum(current_state%zu%data(k, start_y:end_y, start_x:end_x))
        bar_index=bar_index+2
        bartmp(k, bar_index)=sum(current_state%v%data(k, start_y:end_y, start_x:end_x))
        bartmp(k, bar_index+1)=sum(current_state%zv%data(k, start_y:end_y, start_x:end_x))
        bar_index=bar_index+2
        if (current_state%th%active) then
          bartmp(k, bar_index)=sum(current_state%th%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zth%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        !! MASS
        if (current_state%qv%active) then
          bartmp(k, bar_index)=sum(current_state%qv%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqv%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%ql%active) then
          bartmp(k, bar_index)=sum(current_state%ql%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zql%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%qr%active) then
          bartmp(k, bar_index)=sum(current_state%qr%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqr%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%qi%active) then
          bartmp(k, bar_index)=sum(current_state%qi%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqi%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%qs%active) then
          bartmp(k, bar_index)=sum(current_state%qs%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqs%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%qg%active) then
          bartmp(k, bar_index)=sum(current_state%qg%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqg%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%qAitkenSolMass%active) then
          bartmp(k, bar_index)=sum(current_state%qAitkenSolMass%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqAitkenSolMass%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%qAccumSolMass%active) then
          bartmp(k, bar_index)=sum(current_state%qAccumSolMass%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqAccumSolMass%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%qAccumInsolMass%active) then
          bartmp(k, bar_index)=sum(current_state%qAccumInsolMass%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqAccumInsolMass%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%qCoarseSolMass%active) then
          bartmp(k, bar_index)=sum(current_state%qCoarseSolMass%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqCoarseSolMass%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%qCoarseDustMass%active) then
          bartmp(k, bar_index)=sum(current_state%qCoarseDustMass%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zqCoarseDustMass%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

      !! NUMBER
        if (current_state%nl%active) then
          bartmp(k, bar_index)=sum(current_state%nl%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%znl%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%nr%active) then
          bartmp(k, bar_index)=sum(current_state%nr%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%znr%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%ni%active) then
          bartmp(k, bar_index)=sum(current_state%ni%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zni%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%ns%active) then
          bartmp(k, bar_index)=sum(current_state%ns%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zns%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%ng%active) then
          bartmp(k, bar_index)=sum(current_state%ng%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zng%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%nAitkenSolNumber%active) then
          bartmp(k, bar_index)=sum(current_state%nAitkenSolNumber%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%znAitkenSolNumber%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%nAccumSolNumber%active) then
          bartmp(k, bar_index)=sum(current_state%nAccumSolNumber%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%znAccumSolNumber%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%nAccumInsolNumber%active) then
          bartmp(k, bar_index)=sum(current_state%nAccumInsolNumber%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%znAccumInsolNumber%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%nCoarseSolNumber%active) then
          bartmp(k, bar_index)=sum(current_state%nCoarseSolNumber%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%znCoarseSolNumber%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if

        if (current_state%nCoarseDustnumber%active) then
          bartmp(k, bar_index)=sum(current_state%nCoarseDustnumber%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%znCoarseDustnumber%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if
      end do
    end if

    call mpi_allreduce(MPI_IN_PLACE, bartmp, bar_fields*current_state%local_grid%size(Z_INDEX), PRECISION_TYPE, MPI_SUM, &
         current_state%parallel%monc_communicator, ierr)
  end subroutine calculate_sum_profiles

  subroutine check(status)
    integer(4), intent ( in) :: status
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine
end module meanprofiles_mod
