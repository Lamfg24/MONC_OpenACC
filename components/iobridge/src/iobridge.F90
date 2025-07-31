!> Bridge between MONC and the IO server, this registers the current MONC process, will issue data dumps
!! and deregister MONCs when the model run is completed.
module iobridge_mod
  use monc_component_mod, only : COMPONENT_DOUBLE_DATA_TYPE, COMPONENT_INTEGER_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type, &
       component_descriptor_type_v1, component_descriptor_type_v1_array
  use collections_mod, only : hashmap_type, map_type, list_type, c_contains, c_get_generic, c_get_string, c_put_generic, &
       c_put_integer, c_size, c_key_at, c_free
  use conversions_mod, only : conv_to_string
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX, local_grid_type, vertical_grid_configuration_type
  use optionsdatabase_mod, only : options_size, options_get_logical, options_get_integer, options_get_string, options_get_real
  use prognostics_mod, only : prognostic_field_type
  use datadefn_mod, only : DEFAULT_PRECISION, SINGLE_PRECISION, DOUBLE_PRECISION, STRING_LENGTH, PRECISION_TYPE
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log, log_master_log
  !use registry_mod, only : get_all_component_published_fields, get_component_field_value, &
  !     get_component_field_information, is_component_enabled
  use io_server_client_mod, only : COMMAND_TAG, DATA_TAG, REGISTER_COMMAND, DEREGISTER_COMMAND, DATA_COMMAND_START, &
       ARRAY_FIELD_TYPE, SCALAR_FIELD_TYPE, MAP_FIELD_TYPE, INTEGER_DATA_TYPE, BOOLEAN_DATA_TYPE, STRING_DATA_TYPE, &
       FLOAT_DATA_TYPE, DOUBLE_DATA_TYPE, LOCAL_SIZES_KEY, LOCAL_START_POINTS_KEY, LOCAL_END_POINTS_KEY, NUMBER_Q_INDICIES_KEY, &
       data_sizing_description_type, definition_description_type, field_description_type, build_mpi_type_data_sizing_description,&
       build_mpi_type_field_description, build_mpi_type_definition_description, populate_mpi_type_extents, append_mpi_datatype, &
       get_mpi_datatype_from_internal_representation, pack_scalar_field, pack_array_field, pack_map_field
  use mpi, only : MPI_COMM_WORLD, MPI_INT, MPI_BYTE, MPI_REQUEST_NULL, MPI_STATUSES_IGNORE, MPI_STATUS_IGNORE, MPI_STATUS_SIZE, &
                  MPI_SUM, MPI_IN_PLACE, MPI_REAL, MPI_DOUBLE, MPI_LOGICAL, MPI_CHAR
  use q_indices_mod, only : q_metadata_type, get_max_number_q_indices, get_indices_descriptor, get_number_active_q_indices
  use tracers_mod, only : get_tracer_name, reinitialise_trajectories, traj_interval
  use conditional_diagnostics_column_mod, only : ncond, ndiag, cond_request, diag_request, cond_long, diag_long
  !use socrates_couple_mod, only : socrates_couple_get_descriptor


  implicit none

#ifndef TEST_MODE
  private
#endif

  type io_server_sendable_field_sizing
     integer :: number_dimensions, dimensions(4)
  end type io_server_sendable_field_sizing

  type io_configuration_field_type
     character(len=STRING_LENGTH) :: name
     integer :: field_type, data_type
     logical :: optional, enabled
  end type io_configuration_field_type

  type io_configuration_data_definition_type
     character(len=STRING_LENGTH) :: name
     logical :: send_on_terminate
     integer :: number_of_data_fields, frequency, mpi_datatype, command_data
     type(io_configuration_field_type), dimension(:), allocatable :: fields
     integer :: dump_requests(2) !< Dump non blocking send request handles
     character, dimension(:), allocatable :: send_buffer !< Send buffer which holds the model during a dump
  end type io_configuration_data_definition_type

  type(io_configuration_data_definition_type), dimension(:), allocatable :: data_definitions
  type(map_type) :: unique_field_names, sendable_fields, component_field_descriptions
  logical :: io_server_enabled, in_finalisation_callback, socrates_enabled
  integer :: radiation_interval, checkpoint_frequency, init_v = 0
  type(component_descriptor_type_v1) :: socrates_descriptor ! modification
  real(kind=DEFAULT_PRECISION) :: dtmmin
  character(len=200) :: checkpoint_path, diagnostic_path

  public iobridge_get_descriptor, init_callback_iobridge, timestep_callback_iobridge

contains 

  type(component_descriptor_type_v1) function iobridge_get_descriptor()
    iobridge_get_descriptor%name="iobridge"
    iobridge_get_descriptor%version=0.1
    iobridge_get_descriptor%published_fields_on_off = .false.
    !iobridge_get_descriptor%initialisation=>init_callback
    !iobridge_get_descriptor%timestep=>timestep_callback
    !iobridge_get_descriptor%finalisation=>finalisation_callback
  end function iobridge_get_descriptor

  !> Initialisation call back, called at the start of the model run
  !! @param current_state The current model state
  subroutine init_callback_iobridge(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: mpi_type_data_sizing_description, mpi_type_definition_description, mpi_type_field_description, &
               ierr, iter, intnum
    logical :: logicnum

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "enable_io_server") then
        read(current_state%options_database_string(iter,2),*) logicnum
        socrates_enabled = logicnum
        if (.not. socrates_enabled) then
          call log_master_log(LOG_WARN, &
          "Enabled IO bridge but missing IO server compilation, therefore ignoring IO bridge component")
          return
        end if
      end if
    end do
    !if (.not. options_get_logical(current_state%options_database, "enable_io_server")) then
    !  io_server_enabled=.false.
    !  call log_master_log(LOG_WARN, "Enabled IO bridge but missing IO server compilation, therefore ignoring IO bridge component")
    !  return
    !end if

    io_server_enabled=.true.
    in_finalisation_callback=.false.

    call read_and_check_timing_options(current_state)

    !call populate_sendable_fields(current_state)

    !mpi_type_data_sizing_description=build_mpi_type_data_sizing_description()
    !mpi_type_definition_description=build_mpi_type_definition_description()
    !mpi_type_field_description=build_mpi_type_field_description()
    
   !call register_with_io_server(current_state, mpi_type_definition_description, mpi_type_field_description, component_descriptions)
!     !call send_monc_specific_data_to_server(current_state, mpi_type_data_sizing_description)
    

    !call mpi_type_free(mpi_type_data_sizing_description, ierr)
    !call mpi_type_free(mpi_type_definition_description, ierr)
    !call mpi_type_free(mpi_type_field_description, ierr)

    !call build_mpi_data_types()

    !call setup_timing_parameters(current_state)

  end subroutine init_callback_iobridge

  subroutine read_and_check_timing_options(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: iter, intnum
    logical ::  logicnum


  ! Obtain logical switch for "as-needed" diagnostic calculations (only when sampling)
  !current_state%only_compute_on_sample_timestep = &
  !   options_get_logical(current_state%options_database, "only_compute_on_sample_timestep")

  ! Obtain logical switch to ensure that samples are sent on the requestes output_frequency
  ! time_basis=.true. does this automatically
  !current_state%force_output_on_interval = &
  !   options_get_logical(current_state%options_database, "force_output_on_interval")

  ! Obtain logical switch for time_basis handling
  !current_state%time_basis = options_get_logical(current_state%options_database, "time_basis")

    do iter = 1,current_state%config_args
      ! Obtain logical switch for "as-needed" diagnostic calculations (only when sampling)
      if (current_state%options_database_string(iter,1) .eq. "only_compute_on_sample_timestep") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%only_compute_on_sample_timestep = logicnum
      ! Obtain logical switch to ensure that samples are sent on the requestes output_frequency
      ! time_basis=.true. does this automatically
      else if (current_state%options_database_string(iter,1) .eq. "force_output_on_interval") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%force_output_on_interval = logicnum
      ! Obtain logical switch for time_basis handling
      else if (current_state%options_database_string(iter,1) .eq. "time_basis") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%time_basis = logicnum
      !else if (current_state%options_database_string(iter,1) .eq. "socrates_couple_enabled") then
      !  read(current_state%options_database_string(iter,2),*) logicnum
      !  socrates_enabled = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "rad_interval") then
        read(current_state%options_database_string(iter,2),*) intnum
        radiation_interval = intnum
      else if (current_state%options_database_string(iter,1) .eq. "checkpoint_path") then
        checkpoint_path = current_state%options_database_string(iter,2)
      else if (current_state%options_database_string(iter,1) .eq. "diagnostic_path") then
        diagnostic_path = current_state%options_database_string(iter,2)
      end if
    end do

    ! Send logical behaviour message
    if (current_state%force_output_on_interval .and. current_state%time_basis) &
      call log_master_log(LOG_WARN, "Both force_output_on_interval and time_basis are set to "//&
                                    ".true..  Behaviour defaults to that of time_basis.")
  !
!     ! Record SOCRATES information, if enabled
!     socrates_enabled = is_component_enabled(current_state%options_database, "socrates_couple")
!     radiation_interval = options_get_integer(current_state%options_database, "rad_interval")
!     if (socrates_enabled) then
!       socrates_descriptor = socrates_couple_get_descriptor()
!     end if
!
  end subroutine read_and_check_timing_options

  !> Model dump call back, called for each model dump phase
  !! @param current_state The current model state
  subroutine timestep_callback_iobridge(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    type(vertical_grid_configuration_type) :: vertical_grid
    integer :: command = 9999, ierr, rank_sender, size_prog
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: global_array_3d
    real(kind=DEFAULT_PRECISION) :: time_model
    integer :: global_grid_z_size, global_grid_y_size, global_grid_x_size, &
            local_grid_x_size, local_grid_y_size, local_grid_z_size, &
            global_grid_x_start, global_grid_y_start, global_grid_z_start, &
            global_array_size_3d, global_array_size_2d, global_array_size_1d, &
            iter_z, iter_j

    local_grid_z_size = current_state%local_grid%size(Z_INDEX)
    local_grid_y_size = current_state%local_grid%size(Y_INDEX)
    local_grid_x_size = current_state%local_grid%size(X_INDEX)

    global_grid_z_size = current_state%global_grid%size(Z_INDEX)
    global_grid_y_size = current_state%global_grid%size(Y_INDEX)
    global_grid_x_size = current_state%global_grid%size(X_INDEX)

    vertical_grid=current_state%global_grid%configuration%vertical

    if (current_state%parallel%my_global_rank == 1 .and. init_v .eq. 0) then
      call mpi_send(current_state%global_grid%size, 3, MPI_INT, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%termination_time, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%global_grid%resolution, 3, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%global_grid%bottom, 3, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%global_grid%top, 3, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%last_cfl_timestep, 1, MPI_INT, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%th%active, 1, MPI_LOGICAL, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%number_q_fields, 1, MPI_INT, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%checkpoint_frequency, 1, MPI_INT, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%diagnostic_file_0d_write_frequency, 1, MPI_INT, &
                current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%diagnostic_file_1d_write_frequency, 1, MPI_INT, &
                current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%diagnostic_file_2d_write_frequency, 1, MPI_INT, &
                current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%diagnostic_file_3d_write_frequency, 1, MPI_INT, &
                current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
      call mpi_send(vertical_grid%z, global_grid_z_size, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(vertical_grid%zn, global_grid_z_size, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)


      call mpi_send(current_state%buoyancy_enabled, 1, MPI_LOGICAL, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%forcing_enabled, 1, MPI_LOGICAL, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%profile_diagnostics_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%casim_profile_dgs_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%scalar_diagnostics_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%stepfields_enabled, 1, MPI_LOGICAL, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%mean_profiles_enabled, 1, MPI_LOGICAL, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%pstep_enabled, 1, MPI_LOGICAL, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%conditional_diagnostics_column_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%diagnostics_3d_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%conditional_diagnostics_whole_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%subgrid_profile_diagnostics_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%casim_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%simplecloud_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%coriolis_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%diffusion_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%pw_advection_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%th_advection_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%tvd_advection_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%viscosity_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%socrates_enabled, 1, MPI_LOGICAL, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(checkpoint_path, 200, MPI_CHAR, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(diagnostic_path, 200, MPI_CHAR, &
                    current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      init_v = init_v + 1
    end if


    global_grid_z_start = current_state%local_grid%start(Z_INDEX)
    global_grid_y_start = current_state%local_grid%start(Y_INDEX)
    global_grid_x_start = current_state%local_grid%start(X_INDEX)

    global_array_size_3d = global_grid_z_size*global_grid_y_size*global_grid_x_size
    global_array_size_2d = global_grid_y_size*global_grid_x_size
    global_array_size_1d = global_grid_z_size

    if (current_state%parallel%my_global_rank == 1) then
      call mpi_send(current_state%time, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                      2000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%dtm, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                      2000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%timestep, 1, MPI_INT, current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
    end if

    if (current_state%modulo_number_0d .eq. 0) then
      call diagnostic_file_0d_write(current_state, vertical_grid, global_array_size_2d, &
                                global_grid_y_size, global_grid_x_size, &
                                local_grid_y_size, local_grid_x_size, &
                                global_grid_y_start, global_grid_x_start)
    end if

    if (current_state%modulo_number_1d .eq. 0) then
      call diagnostic_file_1d_write(current_state, vertical_grid, global_array_size_1d, &
                                  global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                                  local_grid_z_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_z_start, global_grid_y_start, global_grid_x_start)
    end if

    if (current_state%modulo_number_2d .eq. 0) then
      call diagnostic_file_2d_write(current_state, vertical_grid, global_array_size_2d, &
                                global_grid_y_size, global_grid_x_size, &
                                local_grid_y_size, local_grid_x_size, &
                                global_grid_y_start, global_grid_x_start)
    end if

    if (current_state%modulo_number_3d .eq. 0) then
      call diagnostic_file_3d_write(current_state, vertical_grid, global_array_size_3d, &
                                  global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                                  local_grid_z_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_z_start, global_grid_y_start, global_grid_x_start)
    end if

    if ((current_state%modulo_number_check .eq. 0) .or. &
                                (current_state%time .ge. current_state%termination_time)) then
      call send_data_for_checkpoint(current_state, vertical_grid, global_array_size_3d, &
                                  global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                                  local_grid_z_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_z_start, global_grid_y_start, global_grid_x_start)
    end if
!
!     integer :: i, snc!, request, ierr
!     logical :: data_sent,transmit_on_off !! LAMBERT transmit_on_off
!
!     data_sent=.false.
!
!     if (.not. io_server_enabled) return
!
!     ! Send data definitions for active sampling intervals
!     do snc=1,size(current_state%sampling(:)%interval)
!       if (current_state%sampling(snc)%active) then
!         do i=1, size(data_definitions)
!           if (data_definitions(i)%frequency == current_state%sampling(snc)%interval) then
!             print*,"snc = ",snc
!             print*,"i = ",i
!             print*,"data_definitions(i)%frequency = ",data_definitions(i)%frequency
!             print*,"current_state%sampling(snc)%interval = ", current_state%sampling(snc)%interval
!             call send_data_to_io_server(current_state, i)
!             data_sent=.true.
!           end if
!         end do
!       end if
!     end do
!
!     ! Update time parameters if data was sent
!     if (data_sent) then
!       if (current_state%time_basis) then
!         ! Release the time_basis dtm lock and eliminate accumulated rounding error on 'time'
!         if (current_state%diagnostic_sample_timestep) then
!           current_state%normal_step=.true.
!           current_state%time =                                                            &
!                 real(nint(current_state%time + current_state%dtm),kind=DEFAULT_PRECISION) &
!                                              - current_state%dtm
!         end if
!         ! Increment next sample time if data was written for this interval
!         where(current_state%sampling(:)%active)                                        &
!           current_state%sampling(:)%next_time = current_state%sampling(:)%next_time +  &
!                                                 current_state%sampling(:)%interval
!
!       ! force_output_on_interval
!       else if (current_state%force_output_on_interval) then
!         do i=1,size(current_state%sampling(:))
!           if (current_state%timestep .eq. current_state%sampling(i)%next_step) then
!             ! Release the time_basis dtm lock and eliminate accumulated rounding error on 'time'
!             if (abs(current_state%time+current_state%dtm - current_state%sampling(i)%next_time) &
!                   .lt. dtmmin/10) then
!               current_state%normal_step=.true.
!               current_state%time =                                                          &
!                   real(nint(current_state%time + current_state%dtm),kind=DEFAULT_PRECISION) &
!                                                - current_state%dtm
!               ! Update next output time
!               current_state%sampling(i)%next_time = &
!                    minval(((nint(current_state%time + current_state%dtm)     &
!                             / current_state%sampling(i)%output(:)) + 1)      &
!                           * current_state%sampling(i)%output(:))
!             end if
!             !Update next sample step.
!             current_state%sampling(i)%next_step = current_state%sampling(i)%next_step &
!                                                   + current_state%sampling(i)%interval
!           end if
!         end do
!       end if ! time_basis or force_output_on_interval check
!     end if ! Update time parameters if data was sent
!
!     ! Adjust radiation timings under time_basis
!     if (current_state%radiation_timestep .and. current_state%time_basis) then
!       current_state%normal_step=.true.
!       current_state%time =                                                            &
!             real(nint(current_state%time + current_state%dtm),kind=DEFAULT_PRECISION) &
!                                          - current_state%dtm
!       where (current_state%sampling(:)%radiation)                                     &
!           current_state%sampling(:)%next_time = current_state%sampling(:)%next_time + &
!                                                 current_state%sampling(:)%interval
!     end if
!
  end subroutine timestep_callback_iobridge

  subroutine diagnostic_file_0d_write(current_state, vertical_grid, global_array_size, &
                                global_grid_y_size, global_grid_x_size, &
                                local_grid_y_size, local_grid_x_size, &
                                global_grid_y_start, global_grid_x_start)

    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_y_size, global_grid_x_size, &
                           local_grid_y_size, local_grid_x_size, &
                           global_grid_y_start, global_grid_x_start

    integer :: ierr

    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%rwp)
        call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%iwp)
        call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%swp)
        call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%gwp)
        call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%tot_iwp)
      end if
      call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%lathf)
      call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%lwp)
      call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%reske)
      call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%senhf)
      if (current_state%casim_enabled .eqv. .true.) then
        call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%surface_precip)
      end if
      call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%vwp)
      call gather2D_and_send_Min_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%wmin)
      call gather2D_and_send_Max_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%wmax)
    end if
    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%subke_2d)
    end if

    !print*,"send_data_for_diag_0D MONC"
  end subroutine diagnostic_file_0d_write

  subroutine diagnostic_file_1d_write(current_state, vertical_grid, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                                local_grid_z_size, local_grid_y_size, local_grid_x_size, &
                                global_grid_z_start, global_grid_y_start, global_grid_x_start)

    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           local_grid_z_size, local_grid_y_size, local_grid_x_size, &
                           global_grid_z_start, global_grid_y_start, global_grid_x_start

    integer :: ierr

    !if (current_state%parallel%my_global_rank == 1) then
    !  call mpi_send(current_state%dtm, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
    !                1000, MPI_COMM_WORLD, ierr)
    !end if
    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%uwsg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%vwsg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%uusg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%vvsg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wwsg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tkesg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wtsg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%th2sg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%sed_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                   global_array_size, local_grid_z_size, &
                                   global_grid_z_start, current_state%ssub_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dissipation_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%buoysg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wkesg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%theta_dis_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%vis_coef_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%diff_coef_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%richardson_number_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%richardson_squared_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wqv_sg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wql_sg_tot)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqr_sg_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqi_sg_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqs_sg_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqg_sg_tot)
      end if
    end if

    if (current_state%casim_profile_dgs_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqc_mphys_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqg_mphys_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqi_mphys_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqr_mphys_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqs_mphys_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqv_mphys_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dth_mphys_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dth_cond_evap_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqv_cond_evap_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%phomc_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pinuc_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pidep_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psdep_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%piacw_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psacw_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psacr_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pisub_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%pssub_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pimlt_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psmlt_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psaut_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psaci_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%praut_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pracw_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%prevp_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pgacw_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pgacs_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pgmlt_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pgsub_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psedi_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pseds_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psedr_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psedg_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%psedl_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%pcond_tot)
    end if

    if (current_state%profile_diagnostics_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%u_wind_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%uprime_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%v_wind_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%vprime_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wke_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%ww_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%www_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wwww_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%theta_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%w_wind_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%rh_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wtheta_ad_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wtheta_cn_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%uw_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%vw_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%uv_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%th2_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%thref)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%global_grid%configuration%vertical%rho)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%global_grid%configuration%vertical%rhon)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%thinit)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%qv_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%ql_tot)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%qr_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%qi_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%qs_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%qg_tot)
      end if
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wqv_cn_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wql_cn_tot)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqr_cn_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqi_cn_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqs_cn_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqg_cn_tot)
      end if
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wqv_ad_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%wql_ad_tot)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqr_ad_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqi_ad_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqs_ad_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%wqg_ad_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%cloud_mask_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%cloud_liq_mask_tot)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%cloud_ice_mask_tot)
      end if
    end if

    if (current_state%forcing_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%du_subs_profile_diag)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dv_subs_profile_diag)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dtheta_subs_profile_diag)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqv_subs_profile_diag)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dql_subs_profile_diag)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqr_subs_profile_diag)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqi_subs_profile_diag)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqs_subs_profile_diag)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%dqg_subs_profile_diag)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_u_forc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_v_forc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_th_forc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qv_forc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_ql_forc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qi_forc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qr_forc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qs_forc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qg_forc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_tabs_forc)
    end if

    if (current_state%diffusion_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_th_diff)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qv_diff)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_ql_diff)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qi_diff)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qr_diff)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qs_diff)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qg_diff)
      end if
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_tabs_diff)
    end if
    if (current_state%buoyancy_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_w_buoy)
    end if
    if (current_state%coriolis_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_u_corio)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_v_corio)
    end if
    if (current_state%pstep_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tendp_pr_tot_u_pt)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tendp_pr_tot_v_pt)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tendp_pr_tot_w_pt)
    end if
    if (current_state%pw_advection_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_u_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_v_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_w_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_th_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qv_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_ql_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qi_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qr_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qs_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qg_pwad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_tabs_pwad)
    end if
    if (current_state%stepfields_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_th_sf)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qv_sf)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_ql_sf)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qi_sf)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qr_sf)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qs_sf)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qg_sf)
      end if
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_tabs_sf)
    end if
    if (current_state%th_advection_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_th_thad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_tabs_thad)
    end if
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_u_tvad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_v_tvad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_w_tvad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_th_tvad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_qv_tvad)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_ql_tvad)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qi_tvad)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qr_tvad)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qs_tvad)
        call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                      global_array_size, local_grid_z_size, &
                                      global_grid_z_start, current_state%tend_pr_tot_qg_tvad)
      end if
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_tabs_tvad)
    end if
    if (current_state%viscosity_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_u_visc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_v_visc)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_w_visc)
    end if
    if (current_state%socrates_enabled .eqv. .true.) then
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%cloud_reff_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%longwave_hr_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%shortwave_hr_tot)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_th_lw)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_tabs_lw)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_th_sw)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_tabs_sw)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_th_total)
      call gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                    global_array_size, local_grid_z_size, &
                                    global_grid_z_start, current_state%tend_pr_tot_tabs_total)
    end if

    !print*,"send_data_for_diag_1D MONC"
  end subroutine diagnostic_file_1d_write

  subroutine diagnostic_file_2d_write(current_state, vertical_grid, global_array_size, &
                                global_grid_y_size, global_grid_x_size, &
                                local_grid_y_size, local_grid_x_size, &
                                global_grid_y_start, global_grid_x_start)

    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_y_size, global_grid_x_size, &
                           local_grid_y_size, local_grid_x_size, &
                           global_grid_y_start, global_grid_x_start

    integer :: ierr

    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%qlmax)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%hqlmax)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%cltop)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%clbas)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%vwp)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%lwp)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%wmax)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%wmin)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%reske)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%senhf)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%lathf)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%rwp)
        call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%iwp)
        call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%swp)
        call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%gwp)
        call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                    global_array_size, local_grid_y_size, local_grid_x_size, &
                                    global_grid_y_start, global_grid_x_start, current_state%tot_iwp)
      end if
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%subke_2d)
    end if

    if (current_state%casim_enabled .eqv. .true.) then
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%surface_precip)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%surface_cloudsed)
      call gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                  global_array_size, local_grid_y_size, local_grid_x_size, &
                                  global_grid_y_start, global_grid_x_start, current_state%surface_rainsed)
    end if
  end subroutine diagnostic_file_2d_write

  subroutine diagnostic_file_3d_write(current_state, vertical_grid, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                                local_grid_z_size, local_grid_y_size, local_grid_x_size, &
                                global_grid_z_start, global_grid_y_start, global_grid_x_start)

    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           local_grid_z_size, local_grid_y_size, local_grid_x_size, &
                           global_grid_z_start, global_grid_y_start, global_grid_x_start

    integer :: ierr

    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%u%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%v%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%w%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%th%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%p%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qv%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%ql%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qr%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qi%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qs%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qg%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qAitkenSolMass%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qAccumSolMass%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qAccumInsolMass%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qCoarseSolMass%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%qCoarseDustMass%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%nl%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%nr%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%ni%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%ns%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%ng%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%nAitkenSolNumber%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%nAccumSolNumber%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%nAccumInsolNumber%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%nCoarseSolNumber%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%nCoarseDustnumber%data)
    ! diag bouyancy
    if (current_state%buoyancy_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_w_buoy)
    end if
    ! diag coriolis
    if (current_state%coriolis_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_u_corio)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_v_corio)
    end if
    ! diag diffusion
    if (current_state%diffusion_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_th_diff)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qv_diff)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_ql_diff)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qr_diff)
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qi_diff)
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qs_diff)
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qg_diff)
      end if
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_tabs_diff)
    end if
    ! diag forcing
    if (current_state%forcing_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_u_forc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_v_forc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_th_forc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qv_forc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_ql_forc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qr_forc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qi_forc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qs_forc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qg_forc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_tabs_forc)
    end if
    ! diag pstep
    if (current_state%pstep_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tendp_3d_u_pt)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tendp_3d_v_pt)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tendp_3d_w_pt)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_u_pt)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_v_pt)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_w_pt)
    end if
    ! diag pwadvection
    if (current_state%pw_advection_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_u_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_v_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_w_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_th_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qv_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_ql_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qr_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qi_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qs_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qg_pwad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_tabs_pwad)
    end if
    ! diag stepfields
    if (current_state%stepfields_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_th_sf)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qv_sf)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_ql_sf)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qr_sf)
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qi_sf)
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qs_sf)
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qg_sf)
      end if
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_tabs_sf)
    end if
    ! diag thadvection
    if (current_state%th_advection_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_th_thad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_tabs_thad)
    end if
    ! diag tvadvection
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_u_tvad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_v_tvad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_w_tvad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_th_tvad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_qv_tvad)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_ql_tvad)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qr_tvad)
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qi_tvad)
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qs_tvad)
        call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                              local_grid_z_size, current_state%tend_3d_qg_tvad)
      end if
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_tabs_tvad)
    end if
    ! diag viscosity
    if (current_state%viscosity_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_u_visc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_v_visc)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%tend_3d_w_visc)
    end if
        !print*,"send_data_for_diag_3D MONC"
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%rdAitkenSol%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%rdAccumSol%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%rdCoarseSol%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_cloud%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_rain%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_ice%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_snow%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_graupel%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%RH%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%RI%data)
    ! SOCRATES
    if (current_state%socrates_enabled .eqv. .true.) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%cloud_reff%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%lwrad_hr%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%swrad_hr%data)
      call gather3D_and_send_to_IO_SOCRATES(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                            global_array_size, local_grid_z_size, current_state%tend_3d_tabs_lw)
      call gather3D_and_send_to_IO_SOCRATES(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                            global_array_size, local_grid_z_size, current_state%tend_3d_tabs_sw)
      call gather3D_and_send_to_IO_SOCRATES(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                            global_array_size, local_grid_z_size, current_state%tend_3d_tabs_total)
    end if
  end subroutine diagnostic_file_3d_write

  subroutine send_data_for_checkpoint(current_state, vertical_grid, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                                local_grid_z_size, local_grid_y_size, local_grid_x_size, &
                                global_grid_z_start, global_grid_y_start, global_grid_x_start)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           local_grid_z_size, local_grid_y_size, local_grid_x_size, &
                           global_grid_z_start, global_grid_y_start, global_grid_x_start

    integer :: ierr

    if (current_state%parallel%my_global_rank == 1) then
      call mpi_send(current_state%dtm_new, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%absolute_new_dtm, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%ugal, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%vgal, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%global_grid%configuration%vertical%thref, global_grid_z_size, MPI_DOUBLE, &
                    current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
      call mpi_send(current_state%global_grid%configuration%vertical%prefn, global_grid_z_size, MPI_DOUBLE, &
                    current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
      call mpi_send(vertical_grid%olubar, global_grid_z_size, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(vertical_grid%olvbar, global_grid_z_size, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      if (current_state%th%active) then
        call mpi_send(vertical_grid%olthbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
      end if
      if (current_state%number_q_fields .ne. 0) then
        call mpi_send(vertical_grid%olqvbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqlbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqrbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqibar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqsbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqgbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqAitkenSolMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqAccumSolMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqAccumInsolMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqCoarseSolMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olqCoarseDustMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olnlbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olnrbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olnibar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olnsbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olngbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, &
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olnAitkenSolNumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olnAccumSolNumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olnAccumInsolNumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olnCoarseSolNumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olnCoarseDustnumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
      end if
      call mpi_send(vertical_grid%olzubar, global_grid_z_size, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      call mpi_send(vertical_grid%olzvbar, global_grid_z_size, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process, &
                    1000, MPI_COMM_WORLD, ierr)
      if (current_state%th%active) then
        call mpi_send(vertical_grid%olzthbar, global_grid_z_size, MPI_DOUBLE, &
                    current_state%parallel%corresponding_io_server_process,&
                    1000, MPI_COMM_WORLD, ierr)
      end if
      if (current_state%number_q_fields .ne. 0) then
        call mpi_send(vertical_grid%olzqvbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqlbar, global_grid_z_size, &
                      MPI_DOUBLE, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqrbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqibar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqsbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqgbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqAitkenSolMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqAccumSolMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqAccumInsolMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqCoarseSolMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzqCoarseDustMassbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olznlbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olznrbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olznibar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olznsbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olzngbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olznAitkenSolNumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olznAccumSolNumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olznAccumInsolNumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olznCoarseSolNumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
        call mpi_send(vertical_grid%olznCoarseDustnumberbar, global_grid_z_size, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
      end if
    end if

    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%u%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%v%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%w%data)
    if (current_state%th%active) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%th%data)
    end if
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%p%data)
    if (current_state%number_q_fields .ne. 0) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qv%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%ql%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qr%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qi%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qs%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qg%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qAitkenSolMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qAccumSolMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qAccumInsolMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qCoarseSolMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%qCoarseDustMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%nl%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%nr%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%ni%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%ns%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%ng%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%nAitkenSolNumber%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%nAccumSolNumber%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%nAccumInsolNumber%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%nCoarseSolNumber%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%nCoarseDustnumber%data)
    end if
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%zu%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%zv%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%zw%data)
    if (current_state%th%active) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%zth%data)
    end if
    if (current_state%number_q_fields .ne. 0) then
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqv%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zql%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqr%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqi%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqs%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqg%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqAitkenSolMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqAccumSolMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqAccumInsolMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqCoarseSolMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zqCoarseDustMass%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%znl%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%znr%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zni%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zns%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%zng%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%znAitkenSolNumber%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%znAccumSolNumber%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%znAccumInsolNumber%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%znCoarseSolNumber%data)
      call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                            local_grid_z_size, current_state%znCoarseDustnumber%data)
    end if
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%rdAitkenSol%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%rdAccumSol%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%rdCoarseSol%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_cloud%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_rain%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_ice%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_snow%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%D0_graupel%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%RH%data)
    call gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                          local_grid_z_size, current_state%RI%data)
    !print*,"send_data_for_checkpoint MONC"
  end subroutine send_data_for_checkpoint

  subroutine gather1D_and_send_to_IO(current_state, global_grid_z_size, &
                                   global_array_size, local_grid_z_size, &
                                   global_grid_z_start, data_1D)


    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: global_grid_z_size, global_array_size, &
                           local_grid_z_size, global_grid_z_start
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: data_1D
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: global_array_1d
    integer :: iter_z, ierr

    allocate(global_array_1d(global_grid_z_size))
    global_array_1d = 0.0_DEFAULT_PRECISION

    do iter_z = 1, local_grid_z_size
      global_array_1d(global_grid_z_start+ iter_z -1) = &
                      data_1D(iter_z)
    end do

    if (current_state%parallel%my_rank == 0) then
        call mpi_reduce(MPI_IN_PLACE , global_array_1d, global_array_size, &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    else
      call mpi_reduce(global_array_1d, global_array_1d, global_array_size, &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    end if
    if (current_state%parallel%my_global_rank == 1) then
      call mpi_send(global_array_1d, global_array_size*2, MPI_REAL, current_state%parallel%corresponding_io_server_process,&
                    1000, MPI_COMM_WORLD, ierr)
      deallocate(global_array_1d)
    end if
  end subroutine gather1D_and_send_to_IO

  subroutine gather2D_and_send_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                   global_array_size, local_grid_y_size, local_grid_x_size, &
                                   global_grid_y_start, global_grid_x_start, data_2D)

    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: global_grid_y_size, global_grid_x_size, global_array_size, &
                           local_grid_y_size, local_grid_x_size, &
                           global_grid_y_start, global_grid_x_start
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: data_2D
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: global_array_2d
    integer :: iter_j , iter_i, ierr

    allocate(global_array_2d(global_grid_y_size, global_grid_x_size))
    global_array_2d = 0.0_DEFAULT_PRECISION


    do iter_j = 1, local_grid_y_size
      do iter_i = 1, local_grid_x_size
        global_array_2d(global_grid_y_start+ iter_j -1, global_grid_x_start+ iter_i -1) = &
                    data_2D(iter_j, iter_i)
      end do
    end do
    if (current_state%parallel%my_rank == 0) then
      call mpi_reduce(MPI_IN_PLACE , global_array_2d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    else
      call mpi_reduce(global_array_2d, global_array_2d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    end if
    if (current_state%parallel%my_global_rank == 1) then
      call mpi_send(global_array_2d, global_array_size*2, MPI_REAL, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      deallocate(global_array_2d)
    end if
  end subroutine gather2D_and_send_to_IO

  subroutine gather2D_and_send_Min_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                   global_array_size, local_grid_y_size, local_grid_x_size, &
                                   global_grid_y_start, global_grid_x_start, data_2D)

    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: global_grid_y_size, global_grid_x_size, global_array_size, &
                           local_grid_y_size, local_grid_x_size, &
                           global_grid_y_start, global_grid_x_start
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: data_2D
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: global_array_2d
    integer :: iter_j , iter_i, ierr

    allocate(global_array_2d(global_grid_y_size, global_grid_x_size))
    global_array_2d = 0.0_DEFAULT_PRECISION

    do iter_j = 1, local_grid_y_size
      do iter_i = 1, local_grid_x_size
        global_array_2d(global_grid_y_start+ iter_j -1, global_grid_x_start+ iter_i -1) = &
                    data_2D(iter_j, iter_i)
      end do
    end do
    if (current_state%parallel%my_rank == 0) then
      call mpi_reduce(MPI_IN_PLACE , global_array_2d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    else
      call mpi_reduce(global_array_2d, global_array_2d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    end if
    if (current_state%parallel%my_global_rank == 1) then
      call mpi_send(minval(global_array_2d), 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      deallocate(global_array_2d)
    end if
  end subroutine gather2D_and_send_Min_to_IO

  subroutine gather2D_and_send_Max_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                   global_array_size, local_grid_y_size, local_grid_x_size, &
                                   global_grid_y_start, global_grid_x_start, data_2D)

    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: global_grid_y_size, global_grid_x_size, global_array_size, &
                           local_grid_y_size, local_grid_x_size, &
                           global_grid_y_start, global_grid_x_start
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: data_2D
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: global_array_2d
    integer :: iter_j , iter_i, ierr

    allocate(global_array_2d(global_grid_y_size, global_grid_x_size))
    global_array_2d = 0.0_DEFAULT_PRECISION


    do iter_j = 1, local_grid_y_size
      do iter_i = 1, local_grid_x_size
        global_array_2d(global_grid_y_start+ iter_j -1, global_grid_x_start+ iter_i -1) = &
                    data_2D(iter_j, iter_i)
      end do
    end do
    if (current_state%parallel%my_rank == 0) then
      call mpi_reduce(MPI_IN_PLACE , global_array_2d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    else
      call mpi_reduce(global_array_2d, global_array_2d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    end if
    if (current_state%parallel%my_global_rank == 1) then
      call mpi_send(maxval(global_array_2d), 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      deallocate(global_array_2d)
    end if
  end subroutine gather2D_and_send_Max_to_IO

  subroutine gather2D_and_send_MMM_to_IO(current_state, global_grid_y_size, global_grid_x_size, &
                                   global_array_size, local_grid_y_size, local_grid_x_size, &
                                   global_grid_y_start, global_grid_x_start, data_2D)

    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: global_grid_y_size, global_grid_x_size, global_array_size, &
                           local_grid_y_size, local_grid_x_size, &
                           global_grid_y_start, global_grid_x_start
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: data_2D
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: global_array_2d
    integer :: iter_j , iter_i, ierr
    real :: min_value, mean_value, max_value

    allocate(global_array_2d(global_grid_y_size, global_grid_x_size))
    global_array_2d = 0.0_DEFAULT_PRECISION


    do iter_j = 1, local_grid_y_size
      do iter_i = 1, local_grid_x_size
        global_array_2d(global_grid_y_start+ iter_j -1, global_grid_x_start+ iter_i -1) = &
                    data_2D(iter_j, iter_i)
      end do
    end do
    if (current_state%parallel%my_rank == 0) then
      call mpi_reduce(MPI_IN_PLACE , global_array_2d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    else
      call mpi_reduce(global_array_2d, global_array_2d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    end if
    if (current_state%parallel%my_global_rank == 1) then
      min_value = minval(global_array_2d)
      mean_value = sum(global_array_2d)/size(global_array_2d)
      max_value = maxval(global_array_2d)
      call mpi_send(min_value, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      call mpi_send(mean_value, 1, MPI_DOUBLE, &
                      current_state%parallel%corresponding_io_server_process, 1000, MPI_COMM_WORLD, ierr)
      call mpi_send(max_value, 1, MPI_DOUBLE, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      deallocate(global_array_2d)
    end if
  end subroutine gather2D_and_send_MMM_to_IO

  subroutine gather3D_and_send_to_IO(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                                   global_array_size, local_grid_z_size, data_3D)

    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                           local_grid_z_size
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: data_3D

    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: global_array_3d
    integer :: iter_z, iter_j , iter_i, ierr, y_count = 0, x_count = 0


    allocate(global_array_3d(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    global_array_3d = 0.0_DEFAULT_PRECISION

    do iter_z = 1, local_grid_z_size
      y_count = current_state%global_grid%bottom(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX)
      do iter_j = (current_state%local_grid%start(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX)), &
                  (current_state%local_grid%end(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX))
        x_count = current_state%global_grid%bottom(X_INDEX) + current_state%local_grid%halo_size(X_INDEX)
        y_count = y_count + 1
        do iter_i = (current_state%local_grid%start(X_INDEX) + current_state%local_grid%halo_size(X_INDEX)), &
                    (current_state%local_grid%end(X_INDEX) + current_state%local_grid%halo_size(X_INDEX))
          x_count = x_count + 1
          global_array_3d(iter_z , iter_j - current_state%local_grid%halo_size(Y_INDEX), &
                         iter_i - current_state%local_grid%halo_size(Y_INDEX)) =  &
                         data_3D(iter_z,y_count,x_count)
        end do
      end do
    end do
    if (current_state%parallel%my_rank == 0) then
      call mpi_reduce(MPI_IN_PLACE , global_array_3d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    else
      call mpi_reduce(global_array_3d, global_array_3d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    end if
    if (current_state%parallel%my_global_rank == 1) then
      call mpi_send(global_array_3d, global_array_size*2, MPI_REAL, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      deallocate(global_array_3d)
    end if
  end subroutine gather3D_and_send_to_IO

  subroutine gather3D_and_send_to_IO_SOCRATES(current_state, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                                   global_array_size, local_grid_z_size, data_3D)

    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: global_grid_z_size, global_grid_y_size, global_grid_x_size, global_array_size, &
                           local_grid_z_size
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: data_3D

    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: global_array_3d
    integer :: iter_z, iter_j , iter_i, ierr, y_count = 0, x_count = 0


    allocate(global_array_3d(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    global_array_3d = 0.0_DEFAULT_PRECISION

    do iter_z = 1, local_grid_z_size
      y_count = current_state%global_grid%bottom(Y_INDEX)
      do iter_j = (current_state%local_grid%start(Y_INDEX)), &
                  (current_state%local_grid%end(Y_INDEX))
        x_count = current_state%global_grid%bottom(X_INDEX)
        y_count = y_count + 1
        do iter_i = (current_state%local_grid%start(X_INDEX)), &
                    (current_state%local_grid%end(X_INDEX))
          x_count = x_count + 1
          global_array_3d(iter_z , iter_j , &
                         iter_i) =  &
                         data_3D(iter_z,y_count,x_count)
        end do
      end do
    end do
    if (current_state%parallel%my_rank == 0) then
      call mpi_reduce(MPI_IN_PLACE , global_array_3d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    else
      call mpi_reduce(global_array_3d, global_array_3d, global_array_size, &
                    PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    end if
    if (current_state%parallel%my_global_rank == 1) then
      call mpi_send(global_array_3d, global_array_size*2, MPI_REAL, current_state%parallel%corresponding_io_server_process,&
                      1000, MPI_COMM_WORLD, ierr)
      deallocate(global_array_3d)
    end if
  end subroutine gather3D_and_send_to_IO_SOCRATES
!
!   !> Sends data to the IO server
!   !! @param current_state The current model state
!   !! @param data_index The specific data index to send over
!   subroutine send_data_to_io_server(current_state, data_index)
!     type(model_state_type), target, intent(inout) :: current_state
!     integer, intent(in) :: data_index
!
!     integer :: command_to_send, ierr, tag_close_check=0, tag_close_diag, tag_diag, &
!                tag_command, tag_begin, request, status(MPI_STATUS_SIZE), rank
!
!     if (data_definitions(data_index)%dump_requests(1) .ne. MPI_REQUEST_NULL .or. &
!          data_definitions(data_index)%dump_requests(2) .ne. MPI_REQUEST_NULL) then
!       ! Here wait for previous data dump to complete (consider extending to using buffers for performance)
!       call mpi_waitall(2, data_definitions(data_index)%dump_requests, MPI_STATUSES_IGNORE, ierr)
!     end if
!
!     ! Pack the send buffer and send it to the IO server
!     call pack_send_buffer(current_state, data_definitions(data_index))
!     print*,"send_data_to_io_server000"
!     data_definitions(data_index)%command_data=DATA_COMMAND_START+data_index
!     call mpi_issend(data_definitions(data_index)%command_data, 1, MPI_INT, &
!             current_state%parallel%corresponding_io_server_process, &
!          COMMAND_TAG, MPI_COMM_WORLD, data_definitions(data_index)%dump_requests(1), ierr)  ! LAMBERT envoie les commandes de diag et checkpoint
!     print*,"data_index = ",data_index
!     !print*,"data_definitions(data_index)%name = ", data_definitions(data_index)%name
!     print*,"data_definitions(data_index)%command_data = ", data_definitions(data_index)%command_data
!     !print*,"DATA_TAG+data_index = ",DATA_TAG+data_index
!     print*,"size(data_definitions(data_index)%send_buffer) = ",size(data_definitions(data_index)%send_buffer)
!     call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
!     print*,"rank = ",rank
!     call mpi_issend(data_definitions(data_index)%send_buffer, 1, data_definitions(data_index)%mpi_datatype, &
!          current_state%parallel%corresponding_io_server_process, DATA_TAG+data_index, MPI_COMM_WORLD, &
!          data_definitions(data_index)%dump_requests(2), ierr)  ! envoie les données de monc sous format byte à la fonction data_receive du fichiercommunication.F90 ! LAMBERT
!     print*,"send_buffer DATA_TAG+data_index = ",DATA_TAG+data_index
!     !print*,"test"
!     !****** LAMBERT MODIF *******
!     if (data_definitions(data_index)%command_data .eq. 20) then
!       print*,"before_mpi"
!       call mpi_recv(tag_close_check, 1, MPI_INT, 0, 520, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
!       !call mpi_irecv(tag_close_check, 1, MPI_INT, 0, 520, MPI_COMM_WORLD, request, ierr)
!       print*,"tag_close_check = ",tag_close_check
!       !call mpi_wait
!     else if (data_definitions(data_index)%command_data .eq. 1) then
!       call mpi_recv(tag_begin, 1, MPI_INT, 0, 501, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
!       print*,"tag_begin = ",tag_begin
!     else if (data_definitions(data_index)%command_data .ge. 5 .and. &
!              data_definitions(data_index)%command_data .le. 18) then
!       call mpi_recv(tag_diag, 1, MPI_INT, 0, 5518, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
!       print*,"tag_diag = ",tag_diag
!     else if (data_definitions(data_index)%command_data .eq. 19) then
!       call mpi_recv(tag_close_diag, 1, MPI_INT, 0, 519, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
!       print*,"tag_close_diag = ",tag_close_diag
!   !else
!      !call mpi_recv(tag_command, 1, MPI_INT, 0, 500, MPI_COMM_WORLD, ierr)
!      !print*,"tag_command = ",tag_command
!    end if
!    !****** LAMBERT MODIF *******
!    !print*,"data_definitions(data_index)%command_data = ", data_definitions(data_index)%command_data
!    !print*,"------------------------------------------------------------------------"
!    !print*,"------------------------------------------------------------------------"
!    print*,"send_data_to_io_server1111"
!   end subroutine send_data_to_io_server
!
!   !> Finalisation call back, called at the end of the model run
!   !! @param current_state The current model state
!   subroutine finalisation_callback(current_state)
!     type(model_state_type), target, intent(inout) :: current_state
!
!     integer :: ierr, i
!
!     if (.not. io_server_enabled) return
!     in_finalisation_callback=.true.
!
!     print*,"finalisation_callback0"
!     do i=1, size(data_definitions)
!       if (data_definitions(i)%send_on_terminate) then
!         print*,"finalisation_callbackaaa"
!         call send_data_to_io_server(current_state, i) ! LAMBERT
!         print*,"finalisation_callbackbbb"
!       end if
!       if (data_definitions(i)%dump_requests(1) .ne. MPI_REQUEST_NULL .or. &
!            data_definitions(i)%dump_requests(2) .ne. MPI_REQUEST_NULL) then
!            print*,"finalisation_callbackccc"
!         call mpi_waitall(2, data_definitions(i)%dump_requests, MPI_STATUSES_IGNORE, ierr)
!       end if
!       if (allocated(data_definitions(i)%send_buffer)) deallocate(data_definitions(i)%send_buffer)
!       call mpi_type_free(data_definitions(i)%mpi_datatype, ierr)
!     end do
!     print*,"finalisation_callback and command = ",DEREGISTER_COMMAND
!     call mpi_send(DEREGISTER_COMMAND, 1, MPI_INT, current_state%parallel%corresponding_io_server_process, &
!          COMMAND_TAG, MPI_COMM_WORLD, ierr)
!     print*,"finalisation_callback1"
!   end subroutine finalisation_callback
!
!   !> Builds the MPI data types that correspond to the field descriptions and sizings
!   subroutine build_mpi_data_types()
!     integer :: i, dump_send_buffer_size
!
!     do i=1, size(data_definitions)
!       dump_send_buffer_size=build_mpi_data_type_for_definition(data_definitions(i))
!       data_definitions(i)%dump_requests=(/MPI_REQUEST_NULL, MPI_REQUEST_NULL/)
!       allocate(data_definitions(i)%send_buffer(dump_send_buffer_size))
!     end do
!   end subroutine build_mpi_data_types
!
!   !> Builds the MPI data type for a specific definition with sizing information
!   !! @param specific_data_definition The data definition to build the type for
!   !! @returns The size (in bytes) that the send buffer needs to be to store the data for the MPI operation
!   integer function build_mpi_data_type_for_definition(specific_data_definition)
!     type(io_configuration_data_definition_type), intent(inout) :: specific_data_definition
!
!     integer :: type_extents(5), type_counts, i, j, tempsize, field_start, data_type, field_array_sizes, &
!          temp_size, prev_data_type, old_types(20), offsets(20), block_counts(20), ierr, field_ignores
!     logical :: field_found
!     type(io_server_sendable_field_sizing) :: field_size_info
!
!     type_extents=populate_mpi_type_extents()
!
!     field_start=1
!     data_type=0
!     type_counts=0
!     field_array_sizes=0
!     field_ignores=0
!     do i=1, specific_data_definition%number_of_data_fields
!       if (data_type == 0) then
!         prev_data_type=data_type
!         data_type=specific_data_definition%fields(i)%data_type
!       else
!         if (data_type .ne. specific_data_definition%fields(i)%data_type) then
!           ! For efficient type packing, combine multiple fields with the same type - therefore when the type changes work the previous one pack
!           call append_mpi_datatype(field_start, i-1-field_ignores, field_array_sizes, data_type, &
!                type_extents, prev_data_type, type_counts+1, old_types, offsets, block_counts)
!           field_start=i
!           field_array_sizes=0
!           field_ignores=0
!           prev_data_type=data_type
!           data_type=specific_data_definition%fields(i)%data_type
!           type_counts=type_counts+1
!         end if
!       end if
!
!       if (specific_data_definition%fields(i)%field_type .eq. ARRAY_FIELD_TYPE .or. &
!            specific_data_definition%fields(i)%field_type .eq. MAP_FIELD_TYPE) then
!         ! Grab the field info based upon the field name
!         field_size_info=get_sendable_field_sizing(specific_data_definition%fields(i)%name, field_found)
!         specific_data_definition%fields(i)%enabled=field_found
!         if (.not. field_found .or. field_size_info%number_dimensions == 0) then
!           ! If no field info, or the dimension is 0 then this MONC process is not sending that field - check it is optional
!           if (.not. specific_data_definition%fields(i)%optional) then
!             call log_log(LOG_ERROR, "Non optional field `"//trim(specific_data_definition%fields(i)%name)//&
!                  "' omitted from MONC IO server registration")
!           end if
!           field_ignores=field_ignores+1
!         else
!           ! If the field is specified then use the size data to assemble the field size and append to current size
!           temp_size=1
!           do j=1, field_size_info%number_dimensions
!             temp_size=temp_size*field_size_info%dimensions(j)
!           end do
!           if (specific_data_definition%fields(i)%field_type .eq. MAP_FIELD_TYPE) then
!             field_array_sizes=(field_array_sizes+temp_size*STRING_LENGTH*2)-1
!           else
!             field_array_sizes=(field_array_sizes+temp_size)-1
!           end if
!         end if
!       else
!         if (specific_data_definition%fields(i)%optional) then
!           field_size_info=get_sendable_field_sizing(specific_data_definition%fields(i)%name, field_found)
!           specific_data_definition%fields(i)%enabled=field_found
!           if (.not. field_found) field_ignores=field_ignores+1
!         end if
!       end if
!     end do
!     if (field_start .le. i-1) then
!       ! If there are outstanding fields to process then we do this here
!       call append_mpi_datatype(field_start, i-1, field_array_sizes, data_type, &
!                type_extents, prev_data_type, type_counts+1, old_types, offsets, block_counts)
!       type_counts=type_counts+1
!     end if
!
!     call mpi_type_struct(type_counts, block_counts, offsets, old_types, specific_data_definition%mpi_datatype, ierr)
!     call mpi_type_commit(specific_data_definition%mpi_datatype, ierr)
!     call mpi_type_size(specific_data_definition%mpi_datatype, tempsize, ierr)
!     build_mpi_data_type_for_definition=tempsize
!   end function build_mpi_data_type_for_definition
!
!   !> Populates the sendable field definitions with the field sizing information
!   !! @param current_state The current model state
!   subroutine populate_sendable_fields(current_state)
!     type(model_state_type), target, intent(inout) :: current_state
!
!     type(list_type) :: published_field_descriptors
!     integer :: i
!
!     call populate_globally_visible_sendable_fields(current_state)
!     published_field_descriptors=get_all_component_published_fields()
!     do i=1, c_size(published_field_descriptors)
!       call populate_component_public_field(current_state, c_get_string(published_field_descriptors, i))
!     end do
!   end subroutine populate_sendable_fields
!
!   !> Populates the field information for a specific publically available field offered by one of the components
!   !! @param field_visibility_descriptor The field descriptor which contains sizing information
!   subroutine populate_component_public_field(current_state, field_name)
!     type(model_state_type), target, intent(inout) :: current_state
!     character(len=*), intent(in) :: field_name
!
!     class(*), pointer :: generic_data
!     type(io_server_sendable_field_sizing), pointer :: field_sizing
!     type(component_field_information_type) :: field_information
!     type(component_field_information_type), pointer :: field_information_info_alloc
!
!
!     field_information=get_component_field_information(current_state, field_name)
!     if (field_information%enabled) then
!       allocate(field_information_info_alloc, source=field_information)
!       generic_data=>field_information_info_alloc
!       call c_put_generic(component_field_descriptions, field_name, generic_data, .false.)
!
!       allocate(field_sizing)
!       field_sizing%number_dimensions=field_information%number_dimensions
!       field_sizing%dimensions=field_information%dimension_sizes
!       generic_data=>field_sizing
!       call c_put_generic(sendable_fields, field_name, generic_data, .false.)
!     end if
!   end subroutine populate_component_public_field
!
!   !> Populates the globally visible sendable fields which is a key value pair mapping between name and description of that field
!   !! @param current_state The current model state
!  subroutine populate_globally_visible_sendable_fields(current_state)
!    type(model_state_type), target, intent(inout) :: current_state
  !
!    integer :: x_size, y_size, z_size
!     class(*), pointer :: raw_generic
!
!    z_size=current_state%local_grid%size(Z_INDEX)
!    y_size=current_state%local_grid%size(Y_INDEX)
!    x_size=current_state%local_grid%size(X_INDEX)
!
!     raw_generic=>generate_sendable_description(options_size(current_state%options_database))
!     call c_put_generic(sendable_fields, "options_database", raw_generic, .false.)
!     if (get_number_active_q_indices() .gt. 0) then
!       raw_generic=>generate_sendable_description(get_number_active_q_indices())
!       call c_put_generic(sendable_fields, "q_indicies", raw_generic, .false.)
!     end if
!     raw_generic=>generate_sendable_description(3)
!     call c_put_generic(sendable_fields, "local_grid_size", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "local_grid_start", raw_generic, .false.)
! #ifdef U_ACTIVE
!     raw_generic=>generate_sendable_description()
!     call c_put_generic(sendable_fields, "x_size", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "x_bottom", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "x_resolution", raw_generic, .false.)
!     raw_generic=>generate_sendable_description(z_size, y_size, x_size)
!     call c_put_generic(sendable_fields, "u", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "u_nogal", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "zu", raw_generic, .false.)
!     if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
!       raw_generic=>generate_sendable_description(z_size)
!       call c_put_generic(sendable_fields, "olubar", raw_generic, .false.)
!       call c_put_generic(sendable_fields, "olzubar", raw_generic, .false.)
!     end if
! #endif
! #ifdef V_ACTIVE
!     raw_generic=>generate_sendable_description()
!     call c_put_generic(sendable_fields, "y_size", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "y_bottom", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "y_resolution", raw_generic, .false.)
!     raw_generic=>generate_sendable_description(z_size, y_size, x_size)
!     call c_put_generic(sendable_fields, "v", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "v_nogal", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "zv", raw_generic, .false.)
!     if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
!       raw_generic=>generate_sendable_description(z_size)
!       call c_put_generic(sendable_fields, "olvbar", raw_generic, .false.)
!       call c_put_generic(sendable_fields, "olzvbar", raw_generic, .false.)
!     end if
! #endif
! #ifdef W_ACTIVE
!     raw_generic=>generate_sendable_description(z_size)
!     call c_put_generic(sendable_fields, "z", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "thref", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "prefn", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "rhon", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "rho", raw_generic, .false.)
!     raw_generic=>generate_sendable_description(z_size, y_size, x_size)
!     call c_put_generic(sendable_fields, "w", raw_generic, .false.)
!     call c_put_generic(sendable_fields, "zw", raw_generic, .false.)
! #endif
!     if (current_state%number_q_fields .gt. 0) then
!       raw_generic=>generate_sendable_description(z_size, y_size, x_size, current_state%number_q_fields)
!       call c_put_generic(sendable_fields, "q", raw_generic, .false.)
!       call c_put_generic(sendable_fields, "zq", raw_generic, .false.)
!       if (allocated(current_state%global_grid%configuration%vertical%olqbar)) then
!         raw_generic=>generate_sendable_description(z_size, current_state%number_q_fields)
!         call c_put_generic(sendable_fields, "olqbar", raw_generic, .false.)
!         call c_put_generic(sendable_fields, "olzqbar", raw_generic, .false.)
!       end if
!     end if
!     if (current_state%th%active) then
!       raw_generic=>generate_sendable_description(z_size, y_size, x_size)
!       call c_put_generic(sendable_fields, "th", raw_generic, .false.)
!       call c_put_generic(sendable_fields, "zth", raw_generic, .false.)
!       if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
!         raw_generic=>generate_sendable_description(z_size)
!         call c_put_generic(sendable_fields, "olthbar", raw_generic, .false.)
!         call c_put_generic(sendable_fields, "olzthbar", raw_generic, .false.)
!       end if
!     end if
!     if (current_state%p%active) then
!       raw_generic=>generate_sendable_description(z_size, y_size, x_size)
!       call c_put_generic(sendable_fields, "p", raw_generic, .false.)
!     end if
!     if (current_state%n_tracers .gt. 0) then
!       raw_generic=>generate_sendable_description(z_size, y_size, x_size, current_state%n_tracers)
!       call c_put_generic(sendable_fields, "tracer", raw_generic, .false.)
!       call c_put_generic(sendable_fields, "ztracer", raw_generic, .false.)
!     end if
!     ! need to dump heating rate tendency from socrates radiation
!     if (is_component_enabled(current_state%options_database, "socrates_couple")) then
!        raw_generic=>generate_sendable_description(z_size, y_size, x_size)
!        call c_put_generic(sendable_fields, "sth_lw", raw_generic, .false.)
!        call c_put_generic(sendable_fields, "sth_sw", raw_generic, .false.)
!     endif
!
!     if (is_component_enabled(current_state%options_database, "pdf_analysis")) then
!       if (allocated(current_state%global_grid%configuration%vertical%w_up)) then
!         raw_generic=>generate_sendable_description(z_size)
!         call c_put_generic(sendable_fields, "w_up", raw_generic, .false.)
!       end if
!
!       if (allocated(current_state%global_grid%configuration%vertical%w_dwn)) then
!         raw_generic=>generate_sendable_description(z_size)
!         call c_put_generic(sendable_fields, "w_dwn", raw_generic, .false.)
!       end if
!     end if
!
!   end subroutine populate_globally_visible_sendable_fields
!
!   !> Generates a sendable description based upon the dimension information supplied, missing arguments means that dimension
!   !! does not exist
!   !! @param dim1 Optional size of dimension one
!   !! @param dim2 Optional size of dimension two
!   !! @param dim3 Optional size of dimension three
!   !! @param dim4 Optional size of dimension four
!   !! @returns The corresponding sendable description for the field, with the number of dimensions and sizing of each of these
!   function generate_sendable_description(dim1, dim2, dim3, dim4)
!     integer, intent(in), optional :: dim1, dim2, dim3, dim4
!     class(*), pointer :: generate_sendable_description
!
!     type(io_server_sendable_field_sizing), pointer :: field
!     integer :: number_dimensions
!
!     allocate(field)
!     number_dimensions=0
!     field%dimensions=0
!     if (present(dim1)) then
!       field%dimensions(1)=dim1
!       number_dimensions=number_dimensions+1
!     end if
!     if (present(dim2)) then
!       field%dimensions(2)=dim2
!       number_dimensions=number_dimensions+1
!     end if
!     if (present(dim3)) then
!       field%dimensions(3)=dim3
!       number_dimensions=number_dimensions+1
!     end if
!     if (present(dim4)) then
!       field%dimensions(4)=dim4
!       number_dimensions=number_dimensions+1
!     end if
!     field%number_dimensions=number_dimensions
!     generate_sendable_description=>field
!   end function generate_sendable_description
!
!   !> Sends this MONC specific information to the IO server, which is field info (sizing & availability) as well as meta data
!   !! such as ZN field and Q field names
!   !! @param current_state The current model state
!   !! @param mpi_type_data_sizing_description MPI data type representing the sizing message
!   subroutine send_monc_specific_data_to_server(current_state, mpi_type_data_sizing_description)
!     type(model_state_type), target, intent(inout) :: current_state
!     integer, intent(in) :: mpi_type_data_sizing_description
!
!     type(data_sizing_description_type), dimension(:), allocatable :: data_description
!     character, dimension(:), allocatable :: buffer
!     integer :: number_unique_fields, buffer_size, request_handles(2), ierr
!     real(kind=DEFAULT_PRECISION) :: dreal
!
!     number_unique_fields=c_size(unique_field_names)
!     allocate(data_description(number_unique_fields+4))
!     request_handles(1)=send_data_field_sizes_to_server(current_state, mpi_type_data_sizing_description, &
!        data_description, number_unique_fields)
!     buffer_size=(kind(dreal)*current_state%local_grid%size(Z_INDEX))*2 + (STRING_LENGTH * current_state%number_q_fields &
!                  + STRING_LENGTH * current_state%n_tracers + 4*ncond*STRING_LENGTH + 2*ndiag*STRING_LENGTH )
!     allocate(buffer(buffer_size))
!     request_handles(2)=send_general_monc_information_to_server(current_state, buffer)
!     call mpi_waitall(2, request_handles, MPI_STATUSES_IGNORE, ierr)
!     deallocate(data_description)
!     deallocate(buffer)
!   end subroutine send_monc_specific_data_to_server
!
!   !> Assembles all the data field sizing information and sends this to the IO server
!   !! @param current_state The current model state
!   !! @param mpi_type_data_sizing_description MPI data type representing the sizing message
!   !! @param data_description Data descriptions which will be populated and then sent
!   !! @param number_unique_fields The number of unique fields that we are sending over
!   integer function send_data_field_sizes_to_server(current_state, mpi_type_data_sizing_description, &
!        data_description, number_unique_fields)
!     type(model_state_type), target, intent(inout) :: current_state
!     integer, intent(in) :: mpi_type_data_sizing_description, number_unique_fields
!     type(data_sizing_description_type), dimension(:), intent(inout) :: data_description
!
!     integer :: ierr, i, next_index, request_handle
!     character(len=STRING_LENGTH) :: field_name
!
!     call package_local_monc_decomposition_into_descriptions(current_state, data_description)
!     next_index=5
!     do i=1, number_unique_fields
!       field_name=c_key_at(unique_field_names, i)
!       if (c_contains(sendable_fields, field_name)) then
!         call assemble_individual_description(data_description, next_index, field_name, get_sendable_field_sizing(field_name))
!         next_index=next_index+1
!       end if
!     end do
!
!     call mpi_isend(data_description, next_index-1, mpi_type_data_sizing_description, &
!          current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, request_handle, ierr)
!     send_data_field_sizes_to_server=request_handle
!   end function send_data_field_sizes_to_server
!
!   !> Sends the general MONC information (ZN field and Q field names) to the IO server
!   !! @param current_state The current model state
!   !! @param buffer The communication buffer to use
!   !! @returns Handle to nonblocking send
!   integer function send_general_monc_information_to_server(current_state, buffer)
!     type(model_state_type), target, intent(inout) :: current_state
!     character, dimension(:), intent(inout) :: buffer
!
!     character(len=STRING_LENGTH) :: q_field_name, tracer_name, cd_field_name
!     type(q_metadata_type) :: q_meta_data
!     integer :: current_loc, n, ierr, request_handle
!
!     current_loc=1
!     current_loc=pack_array_field(buffer, current_loc, real_array_1d=current_state%global_grid%configuration%vertical%zn)
!     if (current_state%number_q_fields .gt. 0) then
!       do n=1, current_state%number_q_fields
!         q_meta_data=get_indices_descriptor(n)
!         if (q_meta_data%l_used) then
!           q_field_name=q_meta_data%name
!         else
!           q_field_name="qfield_"//trim(conv_to_string(n))
!         end if
!         current_loc=pack_scalar_field(buffer, current_loc, string_value=q_field_name)
!       end do
!     end if
!
!     if (current_state%n_tracers .gt. 0) then
!       do n=1, current_state%n_tracers
!         tracer_name=get_tracer_name(n, current_state%traj_tracer_index, &
!                       current_state%radioactive_tracer_index, &
!                       current_state%n_radioactive_tracers, current_state%n_tracers)
!         current_loc=pack_scalar_field(buffer, current_loc, string_value=tracer_name)
!       end do
!     end if
!
!     current_loc=pack_array_field(buffer, current_loc, real_array_1d=current_state%global_grid%configuration%vertical%z)
!
!     do n=1,ncond*2
!       if (n .le. ncond) then
!         cd_field_name = cond_request(n)
!         current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
!        cd_field_name = cond_long(n)
!         current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
!       else
!         cd_field_name = ".not. "//trim(cond_request(n-ncond))
!         current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
!         cd_field_name = ".not. "//trim(cond_long(n-ncond))
!         current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
!       end if
!     end do
!     do n=1,ndiag
!       cd_field_name = diag_request(n)
!       current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
!       cd_field_name = diag_long(n)
!       current_loc=pack_scalar_field(buffer, current_loc, string_value=cd_field_name)
!     end do
!
!
!     call mpi_isend(buffer, current_loc-1, MPI_BYTE, current_state%parallel%corresponding_io_server_process, &
!          DATA_TAG, MPI_COMM_WORLD, request_handle, ierr)
!     send_general_monc_information_to_server=request_handle
!   end function send_general_monc_information_to_server
!
!   !> Packages the local MONC decomposition information into descriptions for communication
!   !! @param current_state The current model state
!   !! @param data_description The data description to pack into
!   subroutine package_local_monc_decomposition_into_descriptions(current_state, data_description)
!     type(model_state_type), target, intent(inout) :: current_state
!     type(data_sizing_description_type), dimension(:), intent(inout) :: data_description
!
!     type(io_server_sendable_field_sizing) :: sizing_info
!
!     sizing_info%number_dimensions=3
!     sizing_info%dimensions(Z_INDEX)=current_state%local_grid%size(Z_INDEX)
!     sizing_info%dimensions(Y_INDEX)=current_state%local_grid%size(Y_INDEX)
!     sizing_info%dimensions(X_INDEX)=current_state%local_grid%size(X_INDEX)
!     call assemble_individual_description(data_description, 1, LOCAL_SIZES_KEY, sizing_info)
!     sizing_info%dimensions(Z_INDEX)=current_state%local_grid%start(Z_INDEX)
!     sizing_info%dimensions(Y_INDEX)=current_state%local_grid%start(Y_INDEX)
!     sizing_info%dimensions(X_INDEX)=current_state%local_grid%start(X_INDEX)
!     call assemble_individual_description(data_description, 2, LOCAL_START_POINTS_KEY, sizing_info)
!     sizing_info%dimensions(Z_INDEX)=current_state%local_grid%end(Z_INDEX)
!     sizing_info%dimensions(Y_INDEX)=current_state%local_grid%end(Y_INDEX)
!     sizing_info%dimensions(X_INDEX)=current_state%local_grid%end(X_INDEX)
!     call assemble_individual_description(data_description, 3, LOCAL_END_POINTS_KEY, sizing_info)
!     sizing_info%number_dimensions=1
!     sizing_info%dimensions(1)=get_number_active_q_indices()
!     call assemble_individual_description(data_description, 4, NUMBER_Q_INDICIES_KEY, sizing_info)
!   end subroutine package_local_monc_decomposition_into_descriptions
!
!   !> Retrieves the sizing information associated with a specific field
!   !! @param field_name The field name to look up
!   !! @param field_found Optional flag depicting whether the field was found or not
!   type(io_server_sendable_field_sizing) function get_sendable_field_sizing(field_name, field_found)
!     character(len=*), intent(in) :: field_name
!     logical, intent(out), optional :: field_found
!
!     class(*), pointer :: generic
!
!     if (present(field_found)) field_found=.false.
!     if (c_contains(sendable_fields, field_name)) then
!       generic=>c_get_generic(sendable_fields, field_name)
!       select type(generic)
!       type is (io_server_sendable_field_sizing)
!         get_sendable_field_sizing=generic
!         if (present(field_found)) field_found=.true.
!       end select
!     end if
!   end function get_sendable_field_sizing
!
!   !> Retrieves the descriptor associated with some component's field based upon the field name
!   !! @param field_name The field name
!   type(component_field_information_type) function get_component_field_descriptor(field_name)
!     character(len=*), intent(in) :: field_name
!
!     class(*), pointer :: generic
!     if (c_contains(component_field_descriptions, field_name)) then
!       generic=>c_get_generic(component_field_descriptions, field_name)
!       select type(generic)
!       type is (component_field_information_type)
!         get_component_field_descriptor=generic
!       end select
!     end if
!   end function get_component_field_descriptor
!
!   !> Will assemble an individual description of an array data field
!   !! @param data_description The data structure used to describe the size of arrays
!   !! @param index The index of this current field
!   !! @param field_name The corresponding field name that we are describing
!   !! @param dimensions The number of dimensions (zero means the field is inactive)
!   !! @param dim1 Optional size of dimension 1
!   !! @param dim2 Optional size of dimension 2
!   !! @param dim3 Optional size of dimension 3
!   !! @param dim4 Optional size of dimension 4
!   subroutine assemble_individual_description(data_description, index, field_name, field_sizing_description)
!     integer, intent(in) :: index
!     character(len=*), intent(in) :: field_name
!     type(io_server_sendable_field_sizing), intent(in) :: field_sizing_description
!     type(data_sizing_description_type), dimension(:), intent(inout) :: data_description
!
!     data_description(index)%field_name=field_name
!     data_description(index)%dimensions=field_sizing_description%number_dimensions
!     data_description(index)%dim_sizes=field_sizing_description%dimensions
!   end subroutine assemble_individual_description
!
!   !> Registers this MONC with the corresponding IO server. This will encapsulate the entire protocol, which is sending the
!   !! registration command, receiving the data and field definitions from the IO server and then sending back the sizing
!   !! for the fields that this MONC will contribute.
!   !! Additionally, this receives the unique sampling/output interval pairs from the IO server and
!   !! stores them in the current_state%sampling structure.
!   !! @param current_state The current model state
!   !! @param mpi_type_definition_description MPI data type for data definition message
!   !! @param mpi_type_field_description MPI data type for field definition message
!    subroutine register_with_io_server(current_state, mpi_type_definition_description, &
!               mpi_type_field_description, component_descriptions)
!      type(model_state_type), target, intent(inout) :: current_state
!      integer, intent(in) :: mpi_type_definition_description, mpi_type_field_description
!      type(component_descriptor_type_v1_array), intent(inout) :: component_descriptions

!     type(definition_description_type), dimension(:), allocatable :: definition_descriptions
!     type(field_description_type), dimension(:), allocatable :: field_descriptions
!     integer :: number_defns, number_fields, status(MPI_STATUS_SIZE), ierr, nvalues, i, psize, request_wait, &
!                nprocs, iter_rank, rank
!     integer, dimension(:,:), allocatable :: tmparr
!     logical, dimension(:), allocatable :: mask
!
     !call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr) ! modif
     !call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr) ! modif
     !print*,"nprocs = ",nprocs
     !print*,"register_with_io_server = REGISTER_COMMAND = ",REGISTER_COMMAND
     !do iter_rank = 1, nprocs ! modif
     !if (rank .eq. iter_rank) then ! modif
     !call mpi_send(REGISTER_COMMAND, 1, MPI_INT, current_state%parallel%corresponding_io_server_process, &
     !     COMMAND_TAG, MPI_COMM_WORLD, ierr)

     !call mpi_probe(current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, status, ierr)
     !call mpi_get_count(status, mpi_type_definition_description, number_defns, ierr)
!
!     allocate(definition_descriptions(number_defns))
!
!     call mpi_recv(definition_descriptions, number_defns, mpi_type_definition_description, &
!          current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) ! LAMBERT recoit de ioserver.F90/L661
!     number_fields=get_total_number_of_fields(definition_descriptions, number_defns)
!     !print*,"register_with_io_server_number_defns",number_defns
!     allocate(field_descriptions(number_fields))
!     call mpi_recv(field_descriptions, number_fields, mpi_type_field_description, &
!          current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) ! LAMBERT recoit de ioserver.F90/L664
!
!     call populate_data_definition_configuration(definition_descriptions, number_defns, field_descriptions, number_fields)
!     deallocate(definition_descriptions)
!
!     call mpi_probe(current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, status, ierr)
!     !print*,"register_with_io_server_status",status
!
!     call mpi_get_count(status, MPI_INT, nvalues, ierr)
!     allocate(tmparr(nvalues/2,2))
!     call mpi_recv(tmparr, nvalues, MPI_INT, &
!          current_state%parallel%corresponding_io_server_process, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) ! LAMBERT recoit de ioserver.F90/L667
!     !print*,"register_with_io_server_nvalues",nvalues
!
     !end if !modif
     !end do !modif
!
!     ! Store all unique non-zero user-selected sampling intervals (integers of time or timestep)
!     ! tmparr contains unique sample/output pairs
!     ! tmparr index(:,1)=sampling intervals
!     ! tmparr ineex(:,2)=output intervals
!     psize = 0
!     allocate(mask(nvalues/2))
!     mask(:) = .false.
!     ! Isolate unique sampling intervals
!     do i=1,nvalues/2
!       if ((count( tmparr(i,1) == tmparr(:,1)) .eq. 1)  &
!          .or. .not. any(tmparr(i,1) == tmparr(:,1) .and. mask)) &
!          mask(i) = .true.
!     end do
!
!     ! If radiation is enabled under time_basis and not required every timestep,
!     !   add a sampling entry to track the radiation calculation
!     if (socrates_enabled .and. current_state%time_basis .and. radiation_interval .gt. 0) then
!       psize = count(mask) + 1
!     else
!       psize = count(mask)
!     end if
!
!     ! Allocate the sampling structure and enter the interval values
!     allocate(current_state%sampling(psize))
!     current_state%sampling(1:count(mask))%interval = pack(tmparr(:,1), mask)
!
!     ! For each unique diagnostic sample interval,
!     !   find all non-zero, associated output intervals and store
!     do i=1,size(current_state%sampling(1:count(mask)))
!       mask(:) = tmparr(:,1) == current_state%sampling(i)%interval .and.     &
!                 tmparr(:,2) .gt. 0
!       nvalues=count(mask)
!       allocate(current_state%sampling(i)%output(nvalues))
!       current_state%sampling(i)%output(:) = pack(tmparr(:,2), mask)
!     end do
!     deallocate(mask)
!     deallocate(tmparr)
!
!     ! Populate the radiation interval if not already used
!     if (socrates_enabled .and. current_state%time_basis .and. radiation_interval .gt. 0) then
!       current_state%sampling(psize)%radiation = .true.
!       current_state%sampling(psize)%interval = radiation_interval
!       allocate(current_state%sampling(psize)%output(1))
!       current_state%sampling(psize)%output(1) = -9999  ! this is a dummy item in this case
!     end if
!
!     ! Enforce time_basis sampling regularity (prevents bugs)
!     if (current_state%time_basis) then
!       if (any(mod(current_state%sampling(:)%interval,minval(current_state%sampling(:)%interval)) &
!            .ne. 0)) then
!         call log_master_log(LOG_ERROR, "Under time_basis, all sampling intervals must be "//&
!                  "divisible by the smallest sampling interval ("//&
!                  trim(conv_to_string(minval(current_state%sampling(:)%interval)))//&
!                  " s).  This rule applies system-wide (diagnostics, radiation, checkpoints, etc.")
!
!       end if
!     end if
!
!     ! Check l_constant_dtm for consistency
!     if (options_get_logical(current_state%options_database, "l_constant_dtm")) then
!       if (current_state%time_basis) then
!         if (any(modulo(real( current_state%sampling(:)%interval), real(current_state%dtm)) .gt. 0)) then
!           call log_master_log(LOG_ERROR, "All sampling intervals must be a multiple of dtm "//&
!                                          "when l_constant_dtm=.true.")
!         end if
!       else if (current_state%force_output_on_interval) then
!         call log_master_log(LOG_ERROR, "Use of l_constant_dtm requires force_output_on_interval"//&
!                             "=.false. and time_basis=.false. or time_basis=.true. with all"//&
!                             " sampling intervals a multiple of dtm.")
!       end if
!     end if
! end subroutine register_with_io_server
!
!   !> Retrieve the total number of fields, which is all the fields in all the data definitions
!   !! @param definition_descriptions Data definition descriptions
!   !! @param number_defns The number of data definitions
!   !! @returns The total number of fields
!   integer function get_total_number_of_fields(definition_descriptions, number_defns)
!     type(definition_description_type), dimension(:), intent(inout) :: definition_descriptions
!     integer, intent(in) :: number_defns
!
!     integer :: i
!
!     get_total_number_of_fields=0
!     do i=1, number_defns
!       get_total_number_of_fields=get_total_number_of_fields+definition_descriptions(i)%number_fields
!     end do
!   end function get_total_number_of_fields
!
!   !> Based upon the received data and field definitions this will configure the IO bridge internal representation of these
!   !! facets, which is a structured tree, data defintions holding their own fields rather than the unstructured data
!   !! we get from the IO server
!   !! @param definition_descriptions data definitions received from the IO server
!   !! @param number_defns The number of data definitions
!   !! @param field_descriptions The field descriptions that we received from the IO server
!   !! @param number_fields The total number of field descriptions received
!   subroutine populate_data_definition_configuration(definition_descriptions, number_defns, field_descriptions, number_fields)
!     type(definition_description_type), dimension(:), intent(inout) :: definition_descriptions
!     type(field_description_type), dimension(:), intent(inout) :: field_descriptions
!     integer, intent(in) :: number_defns, number_fields
!
!     integer :: i, definition_index, field_index
!
!     allocate(data_definitions(number_defns))
!     do i=1, number_defns
!       data_definitions(i)%name=definition_descriptions(i)%definition_name
!       data_definitions(i)%send_on_terminate=definition_descriptions(i)%send_on_terminate
!       data_definitions(i)%number_of_data_fields=0 ! Will increment this for each field
!       data_definitions(i)%frequency=definition_descriptions(i)%frequency
!       allocate(data_definitions(i)%fields(definition_descriptions(i)%number_fields))
!     end do
!     do i=1, number_fields
!       definition_index=get_definition_index(field_descriptions(i)%definition_name)
!       field_index=data_definitions(definition_index)%number_of_data_fields+1
!       data_definitions(definition_index)%number_of_data_fields=field_index
!       data_definitions(definition_index)%fields(field_index)%name=field_descriptions(i)%field_name
!       data_definitions(definition_index)%fields(field_index)%field_type=field_descriptions(i)%field_type
!       data_definitions(definition_index)%fields(field_index)%data_type=field_descriptions(i)%data_type
!       data_definitions(definition_index)%fields(field_index)%optional=field_descriptions(i)%optional
!       if (field_descriptions(i)%optional .or. field_descriptions(i)%field_type == ARRAY_FIELD_TYPE .or. &
!            field_descriptions(i)%field_type == MAP_FIELD_TYPE) then
!         call c_put_integer(unique_field_names, field_descriptions(i)%field_name, 1)
!       end if
!       if (.not. field_descriptions(i)%optional) data_definitions(definition_index)%fields(field_index)%enabled=.true.
!
!       ! Verify valid configuration of SOCRATES radiation diagnostics
!       if (socrates_enabled .and. radiation_interval .gt. 0 .and. &
!           definition_descriptions(definition_index)%frequency .gt. 0) then
!         if (any(field_descriptions(i)%field_name .eq. socrates_descriptor%published_fields(:))) then
!           if (mod(definition_descriptions(definition_index)%frequency, radiation_interval) .ne. 0 .or. &
!               definition_descriptions(definition_index)%frequency .lt. radiation_interval) then
!             call log_master_log(LOG_ERROR, "To guarantee availability of radiation "//&
!              "diagnostics, the sampling interval (currently: "//&
!                 trim(conv_to_string(definition_descriptions(definition_index)%frequency))//&
!              ") of the diagnostic ("//trim(field_descriptions(i)%field_name)//&
!              ") must be a MULTIPLE of AND .GE. the radiation calculation interval "//&
!              "(currently: rad_interval="//trim(conv_to_string(radiation_interval))//&
!              ").  This means radiation diagnostics cannot be sampled between calculations."//&
!              "  Check variable's sampling interval (frequency=#) in the xml data-definition." )
!           end if !err
!         end if !compare
!       end if !enabled
!
!     end do
!   end subroutine populate_data_definition_configuration
!
!   !> Looks up a specific definition based upon its name and returns the index
!   !! @param name The defintion name to search for
!   !! @returns The index where the corresponding definition can be found or -1 if no definition is located with that name
!   integer function get_definition_index(name)
!     character(len=*), intent(in) :: name
!
!     integer :: i
!     do i=1, size(data_definitions)
!       if (data_definitions(i)%name .eq. name) then
!         get_definition_index=i
!         return
!       end if
!     end do
!     get_definition_index=-1
!   end function get_definition_index
!
!   !> Packs the current state into the send buffer. This iterates through each field in the data description and adds it to the
!   !! send buffer
!   !! @param current_state The current model state
!   !! @param data_definition The definition of the data which hold the send buffer and the fields
!   subroutine pack_send_buffer(current_state, data_definition)
!     type(model_state_type), target, intent(inout) :: current_state
!     type(io_configuration_data_definition_type), intent(inout) :: data_definition
!
!     integer :: current_buffer_point, i
!
!     current_buffer_point=1
!     do i=1, data_definition%number_of_data_fields
!       if (data_definition%fields(i)%enabled) then
!         if (data_definition%fields(i)%field_type == ARRAY_FIELD_TYPE) then
!           current_buffer_point=pack_array_into_send_buffer(current_state, data_definition, data_definition%fields(i), &
!                current_buffer_point)
!         else if (data_definition%fields(i)%field_type == MAP_FIELD_TYPE) then
!           current_buffer_point=pack_map_into_send_buffer(current_state, data_definition, data_definition%fields(i), &
!                current_buffer_point)
!         else if (data_definition%fields(i)%field_type == SCALAR_FIELD_TYPE) then
!           current_buffer_point=pack_scalar_into_send_buffer(current_state, data_definition, data_definition%fields(i), &
!                current_buffer_point)
!         end if
!       end if
!     end do
!
!     if (current_state%traj_tracer_index .gt. 0 .and. data_definition%name == "3d_tracer_data") then
!       if (mod(nint(current_state%time+current_state%dtm),traj_interval) .eq. 0) then
!         call reinitialise_trajectories(current_state)
!         current_state%reinit_tracer=.true.
!       end if
!     end if
!
!   end subroutine pack_send_buffer
!
!   !> Packs scalar fields into the send bufer
!   !! @param current_state The current model state
!   !! @param data_definition The data definition description
!   !! @param field The specific field we are looking up
!   !! @param current_buffer_point The current point in the buffer where this data will be entered
!   !! @returns The new current buffer point which is after the data addition has taken place
!   integer function pack_scalar_into_send_buffer(current_state, data_definition, field, current_buffer_point)
!     type(model_state_type), target, intent(inout) :: current_state
!     type(io_configuration_data_definition_type), intent(inout) :: data_definition
!     type(io_configuration_field_type), intent(in) :: field
!     integer, intent(in) :: current_buffer_point
!
!     integer :: normal_step_int
!
!     normal_step_int = 0
!
!     if (field%name .eq. "timestep") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            int_value=current_state%timestep)
!     else if (field%name .eq. "terminated") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            logical_value=.not. current_state%continue_timestep .and. in_finalisation_callback)
!     else if (field%name .eq. "z_size") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            int_value=current_state%global_grid%size(Z_INDEX))
!     else if (field%name .eq. "y_size") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            int_value=current_state%global_grid%size(Y_INDEX))
!     else if (field%name .eq. "y_bottom") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%global_grid%bottom(Y_INDEX))
!     else if (field%name .eq. "y_top") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%global_grid%top(Y_INDEX))
!     else if (field%name .eq. "y_resolution") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%global_grid%resolution(Y_INDEX))
!     else if (field%name .eq. "x_size") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            int_value=current_state%global_grid%size(X_INDEX))
!     else if (field%name .eq. "x_bottom") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%global_grid%bottom(X_INDEX))
!     else if (field%name .eq. "x_top") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%global_grid%top(X_INDEX))
!     else if (field%name .eq. "x_resolution") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%global_grid%resolution(X_INDEX))
!     else if (field%name .eq. "time") then
!       ! The time is incremented with dtm as the model was about to increment for the next step and this is needed for diagnostics
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%time+current_state%dtm)
!     else if (field%name .eq. "ugal") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%ugal)
!     else if (field%name .eq. "vgal") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%vgal)
!     else if (field%name .eq. "nqfields") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            int_value=current_state%number_q_fields)
!     else if (field%name .eq. "ntracers") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            int_value=current_state%n_tracers)
!     else if (field%name .eq. "nradtracers") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            int_value=current_state%n_radioactive_tracers)
!     else if (field%name .eq. "dtm") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%dtm)
!     else if (field%name .eq. "dtm_new") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%dtm_new)
!     else if (field%name .eq. "absolute_new_dtm") then
!       pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            real_value=current_state%absolute_new_dtm)
!     else if (field%name .eq. "normal_step") then
!        if (current_state%normal_step) normal_step_int = 1
!        pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!            int_value=normal_step_int)
!     else if (field%name .eq. "rad_last_time") then
!        pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!             real_value=current_state%rad_last_time)
!     else if (field%name .eq. "last_cfl_timestep") then
!        pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!             int_value=current_state%last_cfl_timestep)
!     else if (field%name .eq. "reconfig_timestep_offset") then
!        pack_scalar_into_send_buffer=pack_scalar_field(data_definition%send_buffer, current_buffer_point, &
!             int_value=current_state%reconfig_timestep_offset)
!     else
!       ! Handle component field here
!       pack_scalar_into_send_buffer=handle_component_field_scalar_packing_into_send_buffer(current_state, &
!            data_definition, field, current_buffer_point)
!     end if
!   end function pack_scalar_into_send_buffer
!
!   !> Packs a components field scalar into the send buffer, these are fields that are served up by components rather than
!   !! explicitly available
!   !! @param current_state The current model state
!   !! @param data_definition The data definition description
!   !! @param field The specific field we are looking up
!   !! @param current_buffer_point The current point in the buffer where this data will be entered
!   !! @returns The new current buffer point which is after the data addition has taken place
!   integer function handle_component_field_scalar_packing_into_send_buffer(current_state, data_definition, &
!        field, current_buffer_point)
!     type(model_state_type), target, intent(inout) :: current_state
!     type(io_configuration_data_definition_type), intent(inout) :: data_definition
!     type(io_configuration_field_type), intent(in) :: field
!     integer, intent(in) :: current_buffer_point
!
!     type(component_field_information_type) :: field_descriptor
!     type(component_field_value_type) :: published_value
!
!     field_descriptor=get_component_field_descriptor(field%name)
!     published_value=get_component_field_value(current_state, field%name)
!     if (field_descriptor%data_type == COMPONENT_DOUBLE_DATA_TYPE) then
!       handle_component_field_scalar_packing_into_send_buffer=pack_scalar_field(data_definition%send_buffer, &
!            current_buffer_point, real_value=published_value%scalar_real)
!     else if (field_descriptor%data_type == COMPONENT_INTEGER_DATA_TYPE) then
!       handle_component_field_scalar_packing_into_send_buffer=pack_scalar_field(data_definition%send_buffer, &
!            current_buffer_point, int_value=published_value%scalar_int)
!     end if
!   end function handle_component_field_scalar_packing_into_send_buffer
!
!   !> Packs map fields into the send buffer
!   !! @param current_state The current model state
!   !! @param data_definition The data definition description
!   !! @param field The specific field we are looking up
!   !! @param current_buffer_point The current point in the buffer where this data will be entered
!   !! @returns The new current buffer point which is after the data addition has taken place
!   integer function pack_map_into_send_buffer(current_state, data_definition, field, current_buffer_point)
!     type(model_state_type), target, intent(inout) :: current_state
!     type(io_configuration_data_definition_type), intent(inout) :: data_definition
!     type(io_configuration_field_type), intent(in) :: field
!     integer, intent(in) :: current_buffer_point
!
!     integer :: i
!     type(q_metadata_type) :: specific_q_data
!     type(hashmap_type) :: q_indicies_map
!
!     if (field%name .eq. "options_database") then
!       pack_map_into_send_buffer=pack_map_field(data_definition%send_buffer, current_buffer_point, current_state%options_database)
!     else if (field%name .eq. "q_indicies") then
!       do i=1, get_max_number_q_indices()
!         specific_q_data=get_indices_descriptor(i)
!         if (specific_q_data%l_used) then
!           call c_put_integer(q_indicies_map, specific_q_data%name, i)
!         end if
!       end do
!       pack_map_into_send_buffer=pack_map_field(data_definition%send_buffer, current_buffer_point, q_indicies_map)
!       call c_free(q_indicies_map)
!     end if
!   end function pack_map_into_send_buffer
!
!   !> Packs array fields into the send bufer
!   !! @param current_state The current model state
!   !! @param data_definition The data definition description
!   !! @param field The specific field we are looking up
!   !! @param current_buffer_point The current point in the buffer where this data will be entered
!   !! @returns The new current buffer point which is after the data addition has taken place
!   integer function pack_array_into_send_buffer(current_state, data_definition, field, current_buffer_point)
!     type(model_state_type), target, intent(inout) :: current_state
!     type(io_configuration_data_definition_type), intent(inout) :: data_definition
!     type(io_configuration_field_type), intent(in) :: field
!     integer, intent(in) :: current_buffer_point
!
!     if (field%name .eq. "local_grid_size") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            int_array=current_state%local_grid%size)
!     else if (field%name .eq. "local_grid_start") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            int_array=current_state%local_grid%start)
!     else if (field%name .eq. "z") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%z)
!     else if (field%name .eq. "olubar") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%olubar)
!     else if (field%name .eq. "olzubar") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%olzubar)
!     else if (field%name .eq. "olvbar") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%olvbar)
!     else if (field%name .eq. "olzvbar") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%olzvbar)
!     else if (field%name .eq. "olthbar") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%olthbar)
!     else if (field%name .eq. "olzthbar") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%olzthbar)
!     else if (field%name .eq. "olqbar") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_2d=current_state%global_grid%configuration%vertical%olqbar)
!     else if (field%name .eq. "olzqbar") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_2d=current_state%global_grid%configuration%vertical%olzqbar)
!     else if (field%name .eq. "thref") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%thref)
!     else if (field%name .eq. "prefn") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%prefn)
!     else if (field%name .eq. "rhon") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%rhon)
!     else if (field%name .eq. "rho") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%rho)
!     else if (field%name .eq. "u") then
!       current_state%u%data=current_state%u%data+current_state%ugal
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%u, &
!            current_buffer_point, current_state%local_grid)
!       current_state%u%data=current_state%u%data-current_state%ugal
!     else if (field%name .eq. "u_nogal") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%u, current_buffer_point, &
!            current_state%local_grid)
!     else if (field%name .eq. "zu") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%zu, current_buffer_point,&
!            current_state%local_grid)
!     else if (field%name .eq. "v") then
!       current_state%v%data=current_state%v%data+current_state%vgal
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%v, current_buffer_point, &
!            current_state%local_grid)
!       current_state%v%data=current_state%v%data-current_state%vgal
!     else if (field%name .eq. "v_nogal") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%v, current_buffer_point, &
!            current_state%local_grid)
!     else if (field%name .eq. "zv") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%zv, current_buffer_point,&
!            current_state%local_grid)
!     else if (field%name .eq. "w") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%w, current_buffer_point, &
!            current_state%local_grid)
!     else if (field%name .eq. "zw") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%zw, current_buffer_point,&
!            current_state%local_grid)
!     else if (field%name .eq. "q") then
!       pack_array_into_send_buffer=pack_q_fields(data_definition%send_buffer, current_state%q, current_state%number_q_fields, &
!            current_buffer_point, current_state%local_grid)
!     else if (field%name .eq. "zq") then
!       pack_array_into_send_buffer=pack_q_fields(data_definition%send_buffer, current_state%zq, current_state%number_q_fields, &
!            current_buffer_point, current_state%local_grid)
!     else if (field%name .eq. "tracer") then
!       pack_array_into_send_buffer=pack_q_fields(data_definition%send_buffer, current_state%tracer, current_state%n_tracers, &
!            current_buffer_point, current_state%local_grid)
!     else if (field%name .eq. "ztracer") then
!       pack_array_into_send_buffer=pack_q_fields(data_definition%send_buffer, current_state%ztracer, current_state%n_tracers, &
!            current_buffer_point, current_state%local_grid)
!     else if (field%name .eq. "th") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%th, current_buffer_point,&
!            current_state%local_grid)
!     else if (field%name .eq. "zth") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%zth, current_buffer_point,&
!            current_state%local_grid)
!     else if (field%name .eq. "p") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer, current_state%p, current_buffer_point, &
!            current_state%local_grid)
!     else if (field%name .eq. "sth_lw") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer,  &
!             current_state%sth_lw, current_buffer_point, current_state%local_grid)
!     else if (field%name .eq. "sth_sw") then
!       pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer,  &
!             current_state%sth_sw, current_buffer_point, current_state%local_grid)
!     else if (field%name .eq. "w_up") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%w_up)
!     else if (field%name .eq. "w_dwn") then
!       pack_array_into_send_buffer=pack_array_field(data_definition%send_buffer, current_buffer_point, &
!            real_array_1d=current_state%global_grid%configuration%vertical%w_dwn)
! !!$    else if (field%name .eq. "sw_down_surf") then
! !!$      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer,  &
! !!$            current_state%sth_sw, current_buffer_point, current_state%local_grid)
! !!$    else if (field%name .eq. "lww_down_surf") then
! !!$      pack_array_into_send_buffer=pack_prognostic_flow_field(data_definition%send_buffer,  &
! !!$            current_state%sth_sw, current_buffer_point, current_state%local_grid)
!     else
!        ! Handle component field here
!        pack_array_into_send_buffer=handle_component_field_array_packing_into_send_buffer(current_state, &
!             data_definition, field, current_buffer_point)
!     end if
!   end function pack_array_into_send_buffer
!
!   !> Packs a components field array into the send buffer, these are fields that are served up by components rather than
!   !! explicitly available
!   !! @param current_state The current model state
!   !! @param data_definition The data definition description
!   !! @param field The specific field we are looking up
!   !! @param current_buffer_point The current point in the buffer where this data will be entered
!   !! @returns The new current buffer point which is after the data addition has taken place
!   integer function handle_component_field_array_packing_into_send_buffer(current_state, data_definition, &
!        field, current_buffer_point)
!     type(model_state_type), target, intent(inout) :: current_state
!     type(io_configuration_data_definition_type), intent(inout) :: data_definition
!     type(io_configuration_field_type), intent(in) :: field
!     integer, intent(in) :: current_buffer_point
!
!     type(component_field_information_type) :: field_descriptor
!     type(component_field_value_type) :: published_value
!
!     field_descriptor=get_component_field_descriptor(field%name)
!     published_value=get_component_field_value(current_state, field%name)
!     if (field_descriptor%data_type == COMPONENT_DOUBLE_DATA_TYPE) then
!       if (field_descriptor%number_dimensions == 1) then
!         handle_component_field_array_packing_into_send_buffer=pack_array_field(data_definition%send_buffer, &
!              current_buffer_point, real_array_1d=published_value%real_1d_array)
!         deallocate(published_value%real_1d_array)
!       else if (field_descriptor%number_dimensions == 2) then
!         handle_component_field_array_packing_into_send_buffer=pack_array_field(data_definition%send_buffer, &
!              current_buffer_point, real_array_2d=published_value%real_2d_array)
!         deallocate(published_value%real_2d_array)
!       else if (field_descriptor%number_dimensions == 3) then
!         handle_component_field_array_packing_into_send_buffer=pack_array_field(data_definition%send_buffer, &
!              current_buffer_point, real_array_3d=published_value%real_3d_array)
!         deallocate(published_value%real_3d_array)
!       else if (field_descriptor%number_dimensions == 4) then
!         handle_component_field_array_packing_into_send_buffer=pack_array_field(data_definition%send_buffer, &
!              current_buffer_point, real_array_4d=published_value%real_4d_array)
!         deallocate(published_value%real_4d_array)
!       end if
!     end if
!   end function handle_component_field_array_packing_into_send_buffer
!
!   !> Packs the data of a specific prognostic field into a buffer
!   !! @param buffer The buffer to pack the field into
!   !! @param prognostic The prognostic field
!   !! @param start_offset The starting offset to write into the buffer
!   !! @param local_grid Description of the local grid
!   !! @returns The next location in the buffer to write to (next start offset)
!   integer function pack_prognostic_flow_field(buffer, prognostic, start_offset, local_grid)
!     character, dimension(:), allocatable, intent(inout) :: buffer
!     type(prognostic_field_type), intent(inout) :: prognostic
!     integer, intent(in) :: start_offset
!     type(local_grid_type), intent(inout) :: local_grid
!
!     integer :: target_end
!
!     target_end=start_offset + (local_grid%size(Z_INDEX)*local_grid%size(Y_INDEX)*local_grid%size(X_INDEX)*kind(prognostic%data)-1)
!
!     buffer(start_offset : target_end) = transfer(prognostic%data(&
!          local_grid%local_domain_start_index(Z_INDEX): local_grid%local_domain_end_index(Z_INDEX),&
!          local_grid%local_domain_start_index(Y_INDEX): local_grid%local_domain_end_index(Y_INDEX), &
!          local_grid%local_domain_start_index(X_INDEX): local_grid%local_domain_end_index(X_INDEX)), &
!          buffer(start_offset : target_end))
!     pack_prognostic_flow_field=target_end+1
!   end function pack_prognostic_flow_field
!
!   !> Packs the Q fields into the send buffer
!   !! @param buffer The send buffer to pack into
!   !! @param q_fields Q prognostic fields
!   !! @param number_q_fields The number of Q fields
!   !! @param start_offset Starting offset in the buffer to pack into
!   !! @param local_grid Local grid description
!   !! @returns Updated write location, which is the next location in the buffer to write to
!   integer function pack_q_fields(buffer, q_fields, number_q_fields, start_offset, local_grid)
!     character, dimension(:), allocatable, intent(inout) :: buffer
!     type(prognostic_field_type), dimension(:), intent(inout) :: q_fields
!     integer, intent(in) :: start_offset, number_q_fields
!     type(local_grid_type), intent(inout) :: local_grid
!
!     integer :: target_end, i, current_starting_index
!
!     current_starting_index=start_offset
!
!     do i=1,number_q_fields
!       target_end=current_starting_index + (local_grid%size(Z_INDEX)*local_grid%size(Y_INDEX)*&
!            local_grid%size(X_INDEX)*kind(q_fields(i)%data)-1)
!       buffer(current_starting_index : target_end) = transfer(q_fields(i)%data(&
!            local_grid%local_domain_start_index(Z_INDEX): local_grid%local_domain_end_index(Z_INDEX),&
!            local_grid%local_domain_start_index(Y_INDEX): local_grid%local_domain_end_index(Y_INDEX), &
!            local_grid%local_domain_start_index(X_INDEX): local_grid%local_domain_end_index(X_INDEX)), &
!            buffer(current_starting_index : target_end))
!       current_starting_index=target_end+1
!     end do
!     pack_q_fields=target_end+1
!   end function pack_q_fields
!
!
!   !> Reads in and checks the timing options

!
!
!   !> Initialises timimg paramters
!   !! @param current_state The current model state
!   subroutine setup_timing_parameters(current_state)
!     type(model_state_type), target, intent(inout) :: current_state
!     integer :: sample_nts, next_sample_time, i
!
!     ! Set up the sampling times for time_basis or force_output_on_interval
!     if (current_state%time_basis) then
!       dtmmin = options_get_real(current_state%options_database, "cfl_dtmmin")
!       current_state%sampling(:)%next_time = ((int(current_state%time + dtmmin)             &
!                                                / current_state%sampling(:)%interval) + 1)  &
!                                             * current_state%sampling(:)%interval
!
!     else if (current_state%force_output_on_interval) then
!       dtmmin = options_get_real(current_state%options_database, "cfl_dtmmin")
!       do i=1,size(current_state%sampling(:))
!         current_state%sampling(i)%next_time = minval(((int(current_state%time + dtmmin)    &
!                                                  / current_state%sampling(i)%output(:)) + 1)  &
!                                               * current_state%sampling(i)%output(:))
!         if (size(current_state%sampling(i)%output(:)) .eq. 0) then
!           ! There are no specified output intervals for this sampling interval, so we set the
!           ! next "output time" (which could change dtm) to be the largest possible integer.
!           ! This ensures that in these cases (possibly for a non-zero checkpoint_frequency or
!           ! a radiation interval) the samples simply occur on the timestep interval, and the
!           ! request has no impact on dtm changes.  In the case of a specified non-zero
!           ! checkpoint_frequency, though, checkpoints will be written at that sampling frequency
!           ! without any consideration for the model time.
!           ! Further, note that this only needs to happen here.  Updates to %next_time in this
!           ! module's timestep_callback only occur when writing at an existing %next_time, which
!           ! won't reasonably be reached in this case.
!           current_state%sampling(i)%next_time = huge(current_state%sampling(i)%next_time)
!         else
!           current_state%sampling(i)%next_time = minval(((int(current_state%time + dtmmin)         &
!                                                      / current_state%sampling(i)%output(:)) + 1)  &
!                                                            * current_state%sampling(i)%output(:) )
!         end if
!         current_state%sampling(i)%next_step = (current_state%timestep &
!                                                / current_state%sampling(i)%interval + 1) &
!                                               * current_state%sampling(i)%interval
!         if (mod(current_state%sampling(i)%interval,minval(current_state%sampling(:)%interval)) &
!             .ne. 0) then
!           call log_master_log(LOG_ERROR, "Use of force_output_on_interval requires that all"//&
!              " sampling intervals be evenly divisible by the smallest sampling interval. "//&
!              " Smallest: "//trim(conv_to_string(minval(current_state%sampling(:)%interval)))//&
!              " Conflicting: "//trim(conv_to_string(current_state%sampling(i)%interval)))
!         end if
!       end do
!     end if ! time_basis=.true. or force_output_on_interval=.true.
!
!     ! If we are restarting from a NON-normal_step under time_basis,
!     !   then the next sample steps need to be set now.
!     ! This is similar to the code in cfltest's evaluate_time_basis.
!     if (.not. current_state%normal_step .and. current_state%time_basis) then
!
!       next_sample_time = minval(current_state%sampling(:)%next_time)
!       sample_nts = nint((next_sample_time - current_state%time) / current_state%dtm_new) - 1
!
!       ! Record the next sampling step for intervals matching the next time
!       where(next_sample_time .eq. current_state%sampling(:)%next_time)        &
!         current_state%sampling(:)%next_step = current_state%timestep + sample_nts
!     end if
!
!
!     ! If this is a reconfig_run, it's possible that the previous run did not use time_basis, and it could
!     !   be the case that taking 1 timestep would put us beyond an expected sample time.
!     ! Check for and correct for this case, if needed.
!     !   Correction is to set the dtm to align with the next sample time and set the next_step to the current step.
!     if (current_state%reconfig_run .and. current_state%time_basis .and. current_state%normal_step) then
!       next_sample_time = minval(current_state%sampling(:)%next_time)
!       if (next_sample_time .lt. current_state%time + current_state%dtm) then
!         current_state%dtm = next_sample_time - current_state%time
!         current_state%normal_step = .false.
!         where(next_sample_time .eq. current_state%sampling(:)%next_time)        &
!           current_state%sampling(:)%next_step = current_state%timestep
!       end if ! check for passing next_sample_time in one step
!     end if ! check for reconfig_run with time_basis
!
!
!   end subroutine setup_timing_parameters

end module iobridge_mod
