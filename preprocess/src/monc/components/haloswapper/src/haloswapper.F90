










!> Performs halo swapping. In the parallel case this is between neighbouring processes
!! and in the serial case it still needs to wrap the halos around for the boundary conditions. This module
!! determines the policy of halo swapping (i.e. the fields to communicate) and the halo communication
!! module is used to provide the actual mechanism.
module haloswapper_mod
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use state_mod, only : model_state_type
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_DEBUG, log_get_logging_level, log_log
  use conversions_mod, only : conv_to_string
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, &
       field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, &
       perform_local_data_copy_for_field,&
       init_halo_communication, finalise_halo_communication, blocking_halo_swap, &
       copy_corner_to_buffer, copy_buffer_to_corner
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
  use optionsdatabase_mod, only : options_get_integer
  implicit none

  private
  !< First call tracked for use in displaying debugging information only once
  logical, private :: first_call = .true. 
  type(halo_communication_type), save :: halo_swap_state

  public initialisation_callback_haloswapper, timestep_callback_haloswapper, &
         finalisation_callback_haloswapper

contains

  !> Initialisation callback hook which will set up the halo swapping state and cache some 
  !! precalculated data for fast(er) halo swapping at each timestep
  !! @param current_state The current model state
  subroutine initialisation_callback_haloswapper(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: halo_depth, iter, intnum

    ! get halo_depth and pass it to the halo_swapping routines
    !halo_depth = options_get_integer(current_state%options_database, "halo_depth")
    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "halo_depth") then
        read(current_state%options_database_string(iter,2),*) intnum
        halo_depth = intnum
      end if
    end do
    call init_halo_communication(current_state, get_fields_per_halo_cell, halo_swap_state, &
         halo_depth, .true.)
  end subroutine initialisation_callback_haloswapper

  !> Timestep callback hook which performs the halo swapping for each prognostic field
  !!
  !! In parallel this is performed with MPI communication calls and wrapping around. In serial
  !! still need to wrap data around
  !! @param current_state The current model state_mod
  subroutine timestep_callback_haloswapper(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call blocking_halo_swap(current_state, halo_swap_state, copy_fields_to_halo_buffer, &
         perform_local_data_copy_for_all_prognostics, copy_halo_buffer_to_field, &
         copy_corners_to_halo_buffer, copy_halo_buffer_to_corners)
    call display_debugging_info_if_needed(halo_swap_state%number_distinct_neighbours, &
         halo_swap_state%fields_per_cell, current_state%parallel%my_rank)

  end subroutine timestep_callback_haloswapper
  
  !> The finalisation callback hook which will clean up and free the memory associated with the 
  !! halo swapping
  !! @param current_state The current model state
  subroutine finalisation_callback_haloswapper(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call finalise_halo_communication(halo_swap_state)
  end subroutine finalisation_callback_haloswapper

  !> Copies the halo buffer to halo location in a corner for a halo cell/column and corner 
  !! location.
  !! The copies are performed for each prognostic field.
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are 
  !!                              accessing the buffer of
  !! @param corner_loc The location of the corner
  !! @param x_target_index The target index for the x dimension we are receiving for
  !! @param y_target_index The target index for the y dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been
  !!                     read and copied already)
  !! @param source_data Optional source data which is read from into send buffers and written
  !!                    into by receieve buffers
  subroutine copy_halo_buffer_to_corners(current_state, neighbour_description, corner_loc, &
       x_target_index, y_target_index, neighbour_location, current_page, source_data)

    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, x_target_index, y_target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: page_bookmark, i

    page_bookmark = current_page(neighbour_location)
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer,&
         current_state%u%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%zu%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%v%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%zv%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%w%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
         current_state%zw%data, corner_loc, x_target_index, y_target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    if (current_state%th%active) then
      call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
           current_state%th%data, corner_loc, x_target_index, y_target_index, page_bookmark)
      page_bookmark=page_bookmark+1
      call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
           current_state%zth%data, corner_loc, x_target_index, y_target_index, page_bookmark)
      page_bookmark=page_bookmark+1
    end if
    if (current_state%number_q_fields .gt. 0) then
      if (current_state%qv%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qv%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqv%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ql%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%ql%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zql%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qr%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qr%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqr%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qi%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qi%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqi%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qs%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qs%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqs%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qg%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qg%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqg%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAitkenSolMass%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qAitkenSolMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqAitkenSolMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAccumSolMass%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qAccumSolMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqAccumSolMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAccumInsolMass%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qAccumInsolMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqAccumInsolMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qCoarseSolMass%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qCoarseSolMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqCoarseSolMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qCoarseDustMass%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%qCoarseDustMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zqCoarseDustMass%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nl%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%nl%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%znl%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nr%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%nr%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%znr%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ni%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%ni%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zni%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ns%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%ns%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zns%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ng%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%ng%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%zng%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAitkenSolNumber%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%nAitkenSolNumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%znAitkenSolNumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAccumSolNumber%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%nAccumSolNumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%znAccumSolNumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAccumInsolNumber%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%nAccumInsolNumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%znAccumInsolNumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nCoarseSolNumber%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%nCoarseSolNumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%znCoarseSolNumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nCoarseDustnumber%active) then
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%nCoarseDustnumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
                current_state%znCoarseDustnumber%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
    end if
    if (current_state%n_tracers .gt. 0) then
      do i=1, current_state%n_tracers
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
             current_state%tracer(i)%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_corner(current_state%local_grid, neighbour_description%recv_corner_buffer, &
             current_state%ztracer(i)%data, corner_loc, x_target_index, y_target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end do
    end if

    current_page(neighbour_location)=page_bookmark
  end subroutine copy_halo_buffer_to_corners

  !> Copies the halo buffer to halo location in a field for a specific dimension and halo cell/column.
  !! The copies are performed for each prognostic field.
  !! @param current_state The current model state
  !! @param neighbour_description The halo swapping description of the neighbour we are accessing the buffer of
  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied
  !!        already)
  !! @param source_data Optional source data which is read from into send buffers and written into by receieve 
  !!        buffers
  subroutine copy_halo_buffer_to_field(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: page_bookmark, i

    page_bookmark = current_page(neighbour_location)
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%u%data,  dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%zu%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%v%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%zv%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%w%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
         current_state%zw%data, dim, target_index, page_bookmark)
    page_bookmark=page_bookmark+1
    if (current_state%th%active) then
       call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
            current_state%th%data, dim, target_index, page_bookmark)
       page_bookmark=page_bookmark+1
       call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
            current_state%zth%data, dim, target_index, page_bookmark)
       page_bookmark=page_bookmark+1
    end if
    if (current_state%number_q_fields .gt. 0) then
      if (current_state%qv%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qv%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqv%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ql%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%ql%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zql%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qr%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qr%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqr%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qi%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qi%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqi%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qs%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qs%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqs%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qg%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qg%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqg%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAitkenSolMass%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qAitkenSolMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqAitkenSolMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAccumSolMass%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qAccumSolMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqAccumSolMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAccumInsolMass%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qAccumInsolMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqAccumInsolMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qCoarseSolMass%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qCoarseSolMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqCoarseSolMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qCoarseDustMass%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%qCoarseDustMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zqCoarseDustMass%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nl%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%nl%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%znl%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nr%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%nr%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%znr%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ni%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%ni%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zni%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ns%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%ns%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zns%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ng%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%ng%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%zng%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAitkenSolNumber%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%nAitkenSolNumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%znAitkenSolNumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAccumSolNumber%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%nAccumSolNumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%znAccumSolNumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAccumInsolNumber%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%nAccumInsolNumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%znAccumInsolNumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nCoarseSolNumber%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%nCoarseSolNumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%znCoarseSolNumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nCoarseDustnumber%active) then
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%nCoarseDustnumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%znCoarseDustnumber%data, dim, target_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
    end if
    if (current_state%n_tracers .gt. 0) then
       do i=1, current_state%n_tracers
          call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%tracer(i)%data, dim, target_index, page_bookmark)
          page_bookmark=page_bookmark+1
          call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, &
               current_state%ztracer(i)%data, dim, target_index, page_bookmark)
          page_bookmark=page_bookmark+1
       end do
    end if
    current_page(neighbour_location)=page_bookmark
  end subroutine copy_halo_buffer_to_field

  !> Copies the prognostic field data to halo buffers for a specific process in a dimension and
  !! halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic
  !!        field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from into send buffers and written 
  !!        into by receieve buffers
  subroutine copy_fields_to_halo_buffer(current_state, neighbour_description, dim, source_index,&
       pid_location, current_page, source_data)

    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: page_bookmark, i

    page_bookmark = current_page(pid_location)

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer,&
         current_state%u%data, dim, source_index, page_bookmark)
    page_bookmark = page_bookmark + 1
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%zu%data, dim, source_index, page_bookmark)
    page_bookmark = page_bookmark + 1
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%v%data, dim, source_index, page_bookmark)
    page_bookmark = page_bookmark + 1
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%zv%data,  dim, source_index, page_bookmark)
    page_bookmark = page_bookmark + 1
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%w%data, dim, source_index, page_bookmark)
    page_bookmark = page_bookmark+1
    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, &
         current_state%zw%data,  dim, source_index, page_bookmark)
    page_bookmark = page_bookmark+1
    if (current_state%th%active) then
      call copy_field_to_buffer(current_state%local_grid, &
           neighbour_description%send_halo_buffer, current_state%th%data,  dim, source_index,&
           page_bookmark)
      page_bookmark = page_bookmark+1
      call copy_field_to_buffer(current_state%local_grid, &
           neighbour_description%send_halo_buffer, current_state%zth%data, dim, source_index, &
           page_bookmark)
      page_bookmark = page_bookmark+1
    end if
    if (current_state%number_q_fields .gt. 0) then
      if (current_state%qv%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qv%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqv%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ql%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%ql%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zql%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qr%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qr%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqr%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qi%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qi%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqi%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qs%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qs%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqs%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qg%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qg%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqg%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAitkenSolMass%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qAitkenSolMass%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqAitkenSolMass%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAccumSolMass%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qAccumSolMass%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqAccumSolMass%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAccumInsolMass%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qAccumInsolMass%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqAccumInsolMass%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qCoarseSolMass%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qCoarseSolMass%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqCoarseSolMass%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qCoarseDustMass%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%qCoarseDustMass%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zqCoarseDustMass%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nl%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%nl%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%znl%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nr%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%nr%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%znr%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ns%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%ns%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zns%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ni%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%ni%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zni%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ng%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%ng%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%zng%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAitkenSolNumber%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%nAitkenSolNumber%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%znAitkenSolNumber%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAccumSolNumber%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%nAccumSolNumber%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%znAccumSolNumber%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAccumInsolNumber%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%nAccumInsolNumber%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%znAccumInsolNumber%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nCoarseSolNumber%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%nCoarseSolNumber%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%znCoarseSolNumber%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nCoarseDustnumber%active) then
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, current_state%nCoarseDustnumber%data, dim, &
             source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%znCoarseDustnumber%data, dim, source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
    end if

    if (current_state%n_tracers .gt. 0) then
      do i=1, current_state%n_tracers
        call copy_field_to_buffer(current_state%local_grid, &
               neighbour_description%send_halo_buffer, current_state%tracer(i)%data, dim, &
               source_index, page_bookmark)
        page_bookmark = page_bookmark + 1
        call copy_field_to_buffer(current_state%local_grid, &
             neighbour_description%send_halo_buffer, &
             current_state%ztracer(i)%data, dim, source_index, page_bookmark)
        page_bookmark = page_bookmark + 1
      end do
    end if
    current_page(pid_location) = page_bookmark
  end subroutine copy_fields_to_halo_buffer

  !> Copies the prognostic corner field data to halo buffers for a specific process in a 
  !! dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param corner_loc Location of the corner
  !! @param x_source_index The X source index of the dimension we are reading from in the 
  !! prognostic field
  !! @param y_source_index The Y source index of the dimension we are reading from in the 
  !! prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from into send buffers and written
  !! into by receieve buffers
  subroutine copy_corners_to_halo_buffer(current_state, neighbour_description, corner_loc, &
       x_source_index, y_source_index, pid_location, current_page, source_data)

    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: corner_loc, pid_location, x_source_index, y_source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: page_bookmark, i

    page_bookmark = current_page(pid_location)

    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer,&
         current_state%u%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%zu%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%v%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%zv%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%w%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
    call copy_corner_to_buffer(current_state%local_grid, &
         neighbour_description%send_corner_buffer, &
         current_state%zw%data, corner_loc, x_source_index, y_source_index, page_bookmark)
    page_bookmark=page_bookmark+1
    if (current_state%th%active) then
      call copy_corner_to_buffer(current_state%local_grid, &
           neighbour_description%send_corner_buffer,&
           current_state%th%data,  corner_loc, x_source_index, y_source_index, page_bookmark)
      page_bookmark=page_bookmark+1
      call copy_corner_to_buffer(current_state%local_grid, &
           neighbour_description%send_corner_buffer, &
           current_state%zth%data, corner_loc, x_source_index, y_source_index, page_bookmark)
      page_bookmark=page_bookmark+1
    end if
    if (current_state%number_q_fields .gt. 0) then
      if (current_state%qv%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qv%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqv%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ql%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%ql%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zql%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qr%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qr%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqr%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qi%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qi%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqi%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qs%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qs%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqs%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qg%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qg%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqg%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAitkenSolMass%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qAitkenSolMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqAitkenSolMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAccumSolMass%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qAccumSolMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqAccumSolMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qAccumInsolMass%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qAccumInsolMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqAccumInsolMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qCoarseSolMass%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qCoarseSolMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqCoarseSolMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%qCoarseDustMass%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%qCoarseDustMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zqCoarseDustMass%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nl%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%nl%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%znl%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nr%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%nr%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%znr%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ni%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%ni%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zni%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ns%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%ns%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zns%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%ng%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%ng%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%zng%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAitkenSolNumber%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%nAitkenSolNumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%znAitkenSolNumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAccumSolNumber%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%nAccumSolNumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%znAccumSolNumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nAccumInsolNumber%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%nAccumInsolNumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%znAccumInsolNumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nCoarseSolNumber%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%nCoarseSolNumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%znCoarseSolNumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
      if (current_state%nCoarseDustnumber%active) then
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%nCoarseDustnumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%znCoarseDustnumber%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark + 1
      end if
    end if


    if (current_state%n_tracers .gt. 0) then
      do i=1, current_state%n_tracers
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%tracer(i)%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark+1
        call copy_corner_to_buffer(current_state%local_grid, &
             neighbour_description%send_corner_buffer, &
             current_state%ztracer(i)%data, corner_loc, x_source_index, y_source_index, page_bookmark)
        page_bookmark=page_bookmark+1
      end do
    end if
    current_page(pid_location)=page_bookmark
  end subroutine copy_corners_to_halo_buffer

  !> Deduces the number of fields per halo cell. This depends upon what fields are active in the model
  !! @param current_state The current model state
  integer function get_fields_per_halo_cell(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i

    get_fields_per_halo_cell=0

    get_fields_per_halo_cell=get_fields_per_halo_cell+2
    get_fields_per_halo_cell=get_fields_per_halo_cell+2
    get_fields_per_halo_cell=get_fields_per_halo_cell+2
    if (current_state%th%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qv%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%ql%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qr%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qi%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qs%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qg%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qAitkenSolMass%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qAccumSolMass%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qAccumInsolMass%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qCoarseSolMass%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%qCoarseDustMass%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    !! NUMBER
    if (current_state%nl%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%nr%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%ni%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%ns%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%ng%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%nAitkenSolNumber%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%nAccumSolNumber%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%nAccumInsolNumber%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%nCoarseSolNumber%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if

    if (current_state%nCoarseDustnumber%active) then
      get_fields_per_halo_cell=get_fields_per_halo_cell+2
    end if


    if (current_state%n_tracers .gt. 0) then
      do i=1, current_state%n_tracers
        get_fields_per_halo_cell=get_fields_per_halo_cell+2
      end do
    end if
  end function get_fields_per_halo_cell

  !> Displays some debugging information about who is sending what if that logging_mod level is selected
  !! @param neighbour_counts Number of neighbours that I have
  !! @param included_fields Number of prognostic fields that we include in this halo swap
  !! @param rank My current rank
  subroutine display_debugging_info_if_needed(neighbour_counts, included_fields, rank)
    integer, intent(in) :: neighbour_counts, rank, included_fields

    if (.not. first_call .or. log_get_logging_level() .lt. LOG_DEBUG) return
    call log_log(LOG_DEBUG, "Rank "//trim(conv_to_string(rank))//": "//trim(conv_to_string(neighbour_counts))//&
         " neighbours per timestep over "//trim(conv_to_string(included_fields))//" fields")
    first_call = .false.
  end subroutine display_debugging_info_if_needed

  !> Will do any local copying of data required for the boundary conditions. I.e. if all columns in a slice
  !! are on a process then it will copy in the y dimension
  !! @param current_state The current model state_mod
  !! @param copy_counts Number of local copies performed
  !! @param source_data Optional source data which is read from into send buffers and written into by receieve
  !!                     buffers
  subroutine perform_local_data_copy_for_all_prognostics(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    integer :: i

    call perform_local_data_copy_for_field(current_state%u%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
    call perform_local_data_copy_for_field(current_state%zu%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
    call perform_local_data_copy_for_field(current_state%v%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
    call perform_local_data_copy_for_field(current_state%zv%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
    call perform_local_data_copy_for_field(current_state%w%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
    call perform_local_data_copy_for_field(current_state%zw%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
    if (current_state%th%active) then
      call perform_local_data_copy_for_field(current_state%th%data, current_state%local_grid, &
           current_state%parallel%my_rank, halo_depth, involve_corners)
      call perform_local_data_copy_for_field(current_state%zth%data, current_state%local_grid, &
           current_state%parallel%my_rank, halo_depth, involve_corners)
    end if
    if (current_state%number_q_fields .gt. 0) then
      if (current_state%qv%active) then
        call perform_local_data_copy_for_field(current_state%qv%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqv%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%ql%active) then
        call perform_local_data_copy_for_field(current_state%ql%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zql%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%qr%active) then
        call perform_local_data_copy_for_field(current_state%qr%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqr%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%qi%active) then
        call perform_local_data_copy_for_field(current_state%qi%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqi%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%qs%active) then
        call perform_local_data_copy_for_field(current_state%qs%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqs%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%qg%active) then
        call perform_local_data_copy_for_field(current_state%qg%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqg%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%qAitkenSolMass%active) then
        call perform_local_data_copy_for_field(current_state%qAitkenSolMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqAitkenSolMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%qAccumSolMass%active) then
        call perform_local_data_copy_for_field(current_state%qAccumSolMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqAccumSolMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%qAccumInsolMass%active) then
        call perform_local_data_copy_for_field(current_state%qAccumInsolMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqAccumInsolMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%qCoarseSolMass%active) then
        call perform_local_data_copy_for_field(current_state%qCoarseSolMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqCoarseSolMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%qCoarseDustMass%active) then
        call perform_local_data_copy_for_field(current_state%qCoarseDustMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zqCoarseDustMass%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if

      if (current_state%nl%active) then
        call perform_local_data_copy_for_field(current_state%nl%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%znl%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%nr%active) then
        call perform_local_data_copy_for_field(current_state%nr%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%znr%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%ni%active) then
        call perform_local_data_copy_for_field(current_state%ni%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zni%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%ns%active) then
        call perform_local_data_copy_for_field(current_state%ns%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zns%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%ng%active) then
        call perform_local_data_copy_for_field(current_state%ng%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%zng%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%nAitkenSolNumber%active) then
        call perform_local_data_copy_for_field(current_state%nAitkenSolNumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%znAitkenSolNumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%nAccumSolNumber%active) then
        call perform_local_data_copy_for_field(current_state%nAccumSolNumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%znAccumSolNumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%nAccumInsolNumber%active) then
        call perform_local_data_copy_for_field(current_state%nAccumInsolNumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%znAccumInsolNumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%nCoarseSolNumber%active) then
        call perform_local_data_copy_for_field(current_state%nCoarseSolNumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%znCoarseSolNumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
      if (current_state%nCoarseDustnumber%active) then
        call perform_local_data_copy_for_field(current_state%nCoarseDustnumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%znCoarseDustnumber%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end if
    end if

    if (current_state%n_tracers .gt. 0) then
      do i=1, current_state%n_tracers
        call perform_local_data_copy_for_field(current_state%tracer(i)%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
        call perform_local_data_copy_for_field(current_state%ztracer(i)%data, current_state%local_grid, &
             current_state%parallel%my_rank, halo_depth, involve_corners)
      end do
    end if
  end subroutine perform_local_data_copy_for_all_prognostics
end module haloswapper_mod
