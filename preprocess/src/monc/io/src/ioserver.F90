










!> The main IO server functionality which handles waiting for commands and data both of which are delt with.
!! The lower level details of the communication, configuration parsing etc are all held elsewhere. The server
!! can be thought of similar to a bus, with command and data channels. The command gives context to what is on
!! the data channel and not all commands require data (such as deregistration of MONC process)
module io_server_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH, LONG_STRING_LENGTH
  use configuration_parser_mod, only : DATA_SIZE_STRIDE, io_configuration_type, io_configuration_data_definition_type, &
       io_configuration_registered_monc_type, configuration_parse, extend_registered_moncs_array, retrieve_data_definition, &
       build_definition_description_type_from_configuration, build_field_description_type_from_configuration, get_monc_location, &
       get_io_xml, cond_request, diag_request, cond_long, diag_long, ncond, ndiag, l_thoff
  use state_mod, only : model_state_type
  use mpi_communication_mod, only : build_mpi_datatype, data_receive, test_for_command, register_command_receive, &
       cancel_requests, free_mpi_type, get_number_io_servers, get_my_io_rank, test_for_inter_io, lock_mpi, unlock_mpi, &
       waitall_for_mpi_requests, initialise_mpi_communication, pause_for_mpi_interleaving
  use diagnostic_federator_mod, only : initialise_diagnostic_federator, finalise_diagnostic_federator, &
       check_diagnostic_federator_for_completion, pass_fields_to_diagnostics_federator, determine_diagnostics_fields_available
  use writer_federator_mod, only : initialise_writer_federator, finalise_writer_federator, check_writer_for_trigger, &
       inform_writer_federator_fields_present, inform_writer_federator_time_point, provide_q_field_names_to_writer_federator, &
       provide_tracer_names_to_writer_federator, any_pending
  use writer_field_manager_mod, only : initialise_writer_field_manager, finalise_writer_field_manager, &
       provide_monc_data_to_writer_federator
  use collections_mod, only : hashset_type, hashmap_type, map_type, iterator_type, c_get_integer, c_put_integer, c_is_empty, &
       c_remove, c_add_string, c_integer_at, c_free, c_get_iterator, c_has_next, c_next_mapentry
  use conversions_mod, only : conv_to_string
  use string_utils_mod, only : replace_character
  use io_server_client_mod, only : REGISTER_COMMAND, DEREGISTER_COMMAND, INTER_IO_COMMUNICATION, DATA_COMMAND_START, DATA_TAG, &
       LOCAL_SIZES_KEY, LOCAL_START_POINTS_KEY, LOCAL_END_POINTS_KEY, NUMBER_Q_INDICIES_KEY, SCALAR_FIELD_TYPE, &
       data_sizing_description_type, definition_description_type, field_description_type, build_mpi_type_data_sizing_description,&
       get_data_description_from_name, build_mpi_type_field_description, build_mpi_type_definition_description
  use forthread_mod, only : forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_tryrdlock, &
       forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy, forthread_mutex_init, forthread_mutex_lock, &
       forthread_mutex_unlock, forthread_cond_wait, forthread_cond_signal, forthread_cond_init
  use threadpool_mod, only : threadpool_init, threadpool_finalise, threadpool_start_thread, check_thread_status, &
       threadpool_deactivate, threadpool_is_idle
  use global_callback_inter_io_mod, only : perform_global_callback
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log, initialise_logging
  use optionsdatabase_mod, only : options_get_logical, options_get_string
  use mpi, only : MPI_COMM_WORLD, MPI_STATUSES_IGNORE, MPI_BYTE, MPI_INT, MPI_STATUS_IGNORE, MPI_REAL, MPI_INFO_NULL, MPI_DOUBLE, &
        MPI_LOGICAL, MPI_CHAR
  use io_server_state_reader_mod, only : read_io_server_configuration
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX, local_grid_type, vertical_grid_configuration_type
  use netcdf
  implicit none

  private

  integer :: mpi_type_data_sizing_description, & !< The MPI type for field sizing (i.e. array size etc send when MONCs register)
       mpi_type_definition_description, & !< The MPI data type for data descriptions sent to MONCs
       mpi_type_field_description !< The MPI data type for field descriptions sent to MONCs
  type(io_configuration_type), volatile, save :: io_configuration !< Internal representation of the IO configuration
  logical, volatile :: continue_poll_messages, & !< Whether to continue waiting command messages from any MONC processes
       initialised_present_data
  logical, volatile :: continue_poll_interio_messages, already_registered_finishing_call
  type(field_description_type), dimension(:), allocatable :: registree_field_descriptions
  type(definition_description_type), dimension(:), allocatable :: registree_definition_descriptions
  integer, dimension(:,:), allocatable :: sample_output_pairs
  character(len=200) :: checkpoint_path, diagnostic_path
  integer, volatile :: monc_registration_lock
  integer :: global_grid_z_size, global_grid_y_size, global_grid_x_size
  integer :: time_id = 1

  public io_server_run
contains

  !> Called to start the IO server and once this subroutine returns then it indicates that the IO server has finished.
  !! The runtime is spent in here awaiting commands and then dealing with them. Termination occurs when all MONC processes
  !! have deregistered, note that to trigger this then at least one MONC process must first register
  !! @param io_communicator_arg The IO communicator containing just the IO servers
  !! @param io_xml_configuration Textual XML configuration that is used to set up the IO server
  subroutine io_server_run(current_state, options_database_real, options_database_string, io_communicator_arg, &
       provided_threading, total_global_processes, continuation_run, reconfig_initial_time, io_configuration_file, &
       my_global_rank)
    !type(hashmap_type), intent(inout) ::
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(1200,500), intent(inout) :: options_database_real
    character(len=STRING_LENGTH), dimension(1200,2), intent(inout) :: options_database_string
    integer, intent(in) :: io_communicator_arg, provided_threading, total_global_processes, my_global_rank
    logical, intent(in) :: continuation_run
    real(kind=DEFAULT_PRECISION), intent(in) :: reconfig_initial_time
    character(len=LONG_STRING_LENGTH), intent(in) :: io_configuration_file

    type(vertical_grid_configuration_type) :: vertical_grid
    integer :: number_of_io_servers, number_of_global_moncs, my_io_rank, ierr

    integer :: global_array_size_3d, global_array_size_2d, &
               ncdf_id, status, global_array_size_1d

    integer :: z_size_id, y_size_id, x_size_id, scalar_size_id, init_servor = 0

    logical ::  continue_servor = .true.

    real(kind=DEFAULT_PRECISION) :: time_model, mod_freq
    real :: modulo_number0d, modulo_number1d, modulo_number2d, modulo_number3d, modulo_number_check


    print*,"io_server_run"


    number_of_io_servers=get_number_io_servers(io_communicator_arg)
    number_of_global_moncs=total_global_processes-number_of_io_servers
    my_io_rank=get_my_io_rank(io_communicator_arg)
    call initialise_logging(my_io_rank)

    if (init_servor .eq. 0) then
      call mpi_recv(current_state%global_grid%size, 3, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%termination_time, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%global_grid%resolution, 3, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%global_grid%bottom, 3, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%global_grid%top, 3, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%last_cfl_timestep, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th%active, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%number_q_fields, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%checkpoint_frequency, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%diagnostic_file_0d_write_frequency, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%diagnostic_file_1d_write_frequency, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%diagnostic_file_2d_write_frequency, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%diagnostic_file_3d_write_frequency, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)


      global_grid_z_size = current_state%global_grid%size(Z_INDEX)
      global_grid_y_size = current_state%global_grid%size(Y_INDEX)
      global_grid_x_size = current_state%global_grid%size(X_INDEX)

      global_array_size_3d = global_grid_z_size*global_grid_y_size*global_grid_x_size
      global_array_size_2d = global_grid_y_size*global_grid_x_size
      global_array_size_1d = global_grid_z_size
      vertical_grid=current_state%global_grid%configuration%vertical

      allocate(vertical_grid%z(global_grid_z_size))
      allocate(vertical_grid%zn(global_grid_z_size))
      vertical_grid%z = 0.0_DEFAULT_PRECISION
      vertical_grid%zn = 0.0_DEFAULT_PRECISION
      call mpi_recv(vertical_grid%z, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%zn, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%buoyancy_enabled, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%forcing_enabled, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%profile_diagnostics_enabled, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%casim_profile_dgs_enabled, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%scalar_diagnostics_enabled, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%stepfields_enabled, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_profiles_enabled, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pstep_enabled, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%conditional_diagnostics_column_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%diagnostics_3d_enabled, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%conditional_diagnostics_whole_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%subgrid_profile_diagnostics_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%casim_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%simplecloud_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%coriolis_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%diffusion_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pw_advection_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_advection_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tvd_advection_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%viscosity_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%socrates_enabled, 1, MPI_LOGICAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(checkpoint_path, 200, MPI_CHAR, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(diagnostic_path, 200, MPI_CHAR, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      init_servor = init_servor + 1
    end if


    do while (continue_servor .eqv. .true.)
      call mpi_recv(current_state%time, 1, MPI_DOUBLE, 1, 2000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dtm, 1, MPI_DOUBLE, 1, 2000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%timestep, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      if (current_state%diagnostic_file_0d_write_frequency .eq. 0) then
        modulo_number0d = -1.0
      else
        if (current_state%time_frequency_enabled .eqv. .true.) then
          mod_freq = modulo(current_state%time, float(current_state%diagnostic_file_0d_write_frequency))
          if ((mod_freq - current_state%dtm) .le. 0.0) then
            modulo_number0d = int(mod_freq)
          else
            modulo_number0d = -1.0
          end if
        else
          modulo_number0d = modulo(current_state%timestep, current_state%diagnostic_file_0d_write_frequency)
        end if
      end if
      if (current_state%diagnostic_file_1d_write_frequency .eq. 0) then
        modulo_number1d = -1.0
      else
        if (current_state%time_frequency_enabled .eqv. .true.) then
          mod_freq = modulo(current_state%time, float(current_state%diagnostic_file_1d_write_frequency))
          if ((mod_freq - current_state%dtm) .le. 0.0) then
            modulo_number1d = int(mod_freq)
          else
            modulo_number1d = -1.0
          end if
        else
          modulo_number1d = modulo(current_state%timestep, current_state%diagnostic_file_1d_write_frequency)
        end if
      end if
      if (current_state%diagnostic_file_2d_write_frequency .eq. 0) then
        modulo_number2d = -1.0
      else
        if (current_state%time_frequency_enabled .eqv. .true.) then
          mod_freq = modulo(current_state%time, float(current_state%diagnostic_file_2d_write_frequency))
          if ((mod_freq - current_state%dtm) .le. 0.0) then
            modulo_number2d = int(mod_freq)
          else
            modulo_number2d = -1.0
          end if
        else
          modulo_number2d = modulo(current_state%timestep, current_state%diagnostic_file_2d_write_frequency)
        end if
      end if
      if (current_state%diagnostic_file_3d_write_frequency .eq. 0) then
        modulo_number3d = -1.0
      else
        if (current_state%time_frequency_enabled .eqv. .true.) then
          mod_freq = modulo(current_state%time, float(current_state%diagnostic_file_3d_write_frequency))
          if ((mod_freq - current_state%dtm) .le. 0.0) then
            modulo_number3d = int(mod_freq)
          else
            modulo_number3d = -1.0
          end if
        else
          modulo_number3d = modulo(current_state%timestep, current_state%diagnostic_file_3d_write_frequency)
        end if
      end if
      if (current_state%checkpoint_frequency .eq. 0) then
        modulo_number_check = -1.0
      else
        if (current_state%time_frequency_enabled .eqv. .true.) then
          mod_freq = modulo(current_state%time, float(current_state%checkpoint_frequency))
          if ((mod_freq - current_state%dtm) .le. 0.0) then
            modulo_number_check = int(mod_freq)
          else
            modulo_number_check = -1.0
          end if
        else
          modulo_number_check = modulo(current_state%timestep, current_state%checkpoint_frequency)
        end if
      end if


      if (modulo_number0d .eq. 0.0) then
        call diagnostic_file_0d_generation(current_state, vertical_grid, io_communicator_arg)
      end if

      if (modulo_number1d .eq. 0.0) then
        call diagnostic_file_1d_generation(current_state, vertical_grid, io_communicator_arg, global_array_size_1d, &
                    global_grid_z_size, global_grid_y_size, global_grid_x_size)
      end if

      if (modulo_number2d .eq. 0.0) then
        call diagnostic_file_2d_generation(current_state, vertical_grid, io_communicator_arg, global_array_size_2d, &
                                  global_grid_z_size, global_grid_y_size, global_grid_x_size)
      end if

      if (modulo_number3d .eq. 0.0) then
        call diagnostic_file_3d_generation(current_state, vertical_grid, io_communicator_arg, global_array_size_3d, &
                                  global_grid_z_size, global_grid_y_size, global_grid_x_size)
      end if

      if ((modulo_number_check .eq. 0.0) .or. &
                                (current_state%time .ge. current_state%termination_time)) then
        call checkpoint_file_generation(current_state, vertical_grid, io_communicator_arg, global_array_size_3d, &
                  global_grid_z_size, global_grid_y_size, global_grid_x_size, time_model)
      end if

      if (current_state%time .ge. current_state%termination_time) then
        continue_servor = .false.
      end if
    end do
    print*,"end of servor IO"

  end subroutine io_server_run

  subroutine diagnostic_file_0d_generation(current_state, vertical_grid, io_communicator_arg)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: io_communicator_arg

    integer :: scalar_size_id, time_id
    integer :: ncdf_id, ierr
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename
    real(kind=DEFAULT_PRECISION) :: gwp_min, gwp_mean, gwp_max, iwp_min, iwp_mean, iwp_max, lathf_min, lathf_mean, lathf_max, &
                                    lwp_min, lwp_mean, lwp_max, reske_min, reske_mean, reske_max, rwp_min, rwp_mean, rwp_max, &
                                    senhf_min, senhf_mean, senhf_max, subke_2d_min, subke_2d_mean, subke_2d_max, &
                                    surface_precip_min, surface_precip_mean, surface_precip_max, swp_min, swp_mean, swp_max, &
                                    tot_iwp_min, tot_iwp_mean, tot_iwp_max, vwp_min, vwp_mean, vwp_max, w_min, w_max

    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(rwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(rwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(rwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(iwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(iwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(iwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(swp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(swp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(swp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(gwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(gwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(gwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(tot_iwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(tot_iwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(tot_iwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(lathf_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lathf_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lathf_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(reske_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(reske_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(reske_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(senhf_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(senhf_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(senhf_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if (current_state%casim_enabled .eqv. .true.) then
        call mpi_recv(surface_precip_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(surface_precip_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(surface_precip_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(vwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(w_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(w_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(subke_2d_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(subke_2d_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(subke_2d_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ls1 = len_trim(diagnostic_path)
    ls2 = 0
    do i = 1,ls1
        if(diagnostic_path(i:i).ne.' ') then
          ls2 = ls2 + 1
        endif
    enddo
    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = diagnostic_path(:ls2)//"/full_diag_0d_time_"//trim(conv_to_string(int(current_state%time)))//".nc"
    else
      unique_filename = diagnostic_path(:ls2)//"/full_diag_0d_timestep_"//trim(conv_to_string(current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "time", 1, time_id))
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))

    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)
    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "rwp_min", "minimum rain water path", rwp_min, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "rwp_mean", "mean rain water path", rwp_mean, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "rwp_max", "maximum rain water path", rwp_max, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "iwp_min", "minimum ice water path", iwp_min, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "iwp_mean", "mean ice water path", iwp_mean, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "iwp_max", "maximum ice water path", iwp_max, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "swp_min", "minimum snow water path", swp_min, "g.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "swp_mean", "mean snow water path", swp_mean, "g.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "swp_max", "maximum snow water path", swp_max, "g.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "gwp_min", "minimum graupel water path", gwp_min, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "gwp_mean", "mean graupel water path", gwp_mean, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "gwp_max", "maximum graupel water path", gwp_max, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "tot_iwp_min", "minimum total ice water path", tot_iwp_min, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "tot_iwp_mean", "mean total ice water path", tot_iwp_mean, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "tot_iwp_max", "maximum total ice water path", tot_iwp_max, "kg.m-2", current_state%time)
      end if
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lathf_min", "minimum surface latent heat flux", lathf_min, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lathf_mean", "mean surface latent heat flux", lathf_mean, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lathf_max", "maximum surface latent heat flux", lathf_max, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lwp_min", "minimum liquid water path", lwp_min, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lwp_mean", "mean liquid water path", lwp_mean, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lwp_max", "maximum liquid water path", lwp_max, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "reske_min", "minimum resolved ke", reske_min, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "reske_mean", "mean resolved ke", reske_mean, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "reske_max", "maximum resolved ke", reske_max, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "senhf_min", "minimum surface sensible heat flux", senhf_min, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "senhf_mean", "mean surface sensible heat flux", senhf_mean, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "senhf_max", "maximum surface sensible heat flux", senhf_max, "W.m-2", current_state%time)
      if (current_state%casim_enabled .eqv. .true.) then
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "surface_precip_min", "minimum surface precipitation", surface_precip_min, "???", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "surface_precip_mean", "mean surface precipitation", surface_precip_mean, "???", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "surface_precip_max", "maximum surface precipitation", surface_precip_max, "???", current_state%time)
      end if
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "vwp_min", "minimum water vapour path", vwp_min, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "vwp_mean", "mean water vapour path", vwp_mean, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "vwp_max", "maximum water vapour path", vwp_max, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "w_min", "minimum vertical wind speed", w_min, "m.s-1", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "w_max", "maximum vertical wind speed", w_max, "m.s-1", current_state%time)
    end if
    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "subke_2d_min", "minimum subgrid ke", subke_2d_min, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "subke_2d_mean", "mean subgrid ke", subke_2d_mean, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "subke_2d_max", "maximum subgrid ke", subke_2d_max, "J or kg.m2.s-2", current_state%time)
    end if

    call check(nf90_close(ncdf_id))

  end subroutine diagnostic_file_0d_generation

  subroutine diagnostic_file_1d_generation(current_state, vertical_grid, io_communicator_arg, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           io_communicator_arg

    integer :: z_size_id, y_size_id, x_size_id, scalar_size_id, time_id
    integer :: ncdf_id, ierr, surface
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename

    surface = global_grid_x_size * global_grid_y_size

    allocate(current_state%diff_coef_tot(global_grid_z_size))
    allocate(current_state%dissipation_tot(global_grid_z_size))
    allocate(current_state%dqc_mphys_tot(global_grid_z_size))
    allocate(current_state%dqg_mphys_tot(global_grid_z_size))
    allocate(current_state%dqi_mphys_tot(global_grid_z_size))
    allocate(current_state%dqr_mphys_tot(global_grid_z_size))
    allocate(current_state%dqs_mphys_tot(global_grid_z_size))
    allocate(current_state%dqv_mphys_tot(global_grid_z_size))
    allocate(current_state%dth_mphys_tot(global_grid_z_size))
    allocate(current_state%dth_cond_evap_tot(global_grid_z_size))
    allocate(current_state%dqv_cond_evap_tot(global_grid_z_size))
    allocate(current_state%qs_tot(global_grid_z_size))
    allocate(current_state%qg_tot(global_grid_z_size))
    allocate(current_state%dqg_subs_profile_diag(global_grid_z_size))
    allocate(current_state%cloud_ice_mask_tot(global_grid_z_size))
    allocate(current_state%qi_tot(global_grid_z_size))
    allocate(current_state%dqi_subs_profile_diag(global_grid_z_size))
    allocate(current_state%cloud_liq_mask_tot(global_grid_z_size))
    allocate(current_state%ql_tot(global_grid_z_size))
    allocate(current_state%pcond_tot(global_grid_z_size))
    allocate(current_state%pgacs_tot(global_grid_z_size))
    allocate(current_state%pgacw_tot(global_grid_z_size))
    allocate(current_state%pgmlt_tot(global_grid_z_size))
    allocate(current_state%pgsub_tot(global_grid_z_size))
    allocate(current_state%phomc_tot(global_grid_z_size))
    allocate(current_state%piacw_tot(global_grid_z_size))
    allocate(current_state%pidep_tot(global_grid_z_size))
    allocate(current_state%pimlt_tot(global_grid_z_size))
    allocate(current_state%pinuc_tot(global_grid_z_size))
    allocate(current_state%pisub_tot(global_grid_z_size))
    allocate(current_state%pracw_tot(global_grid_z_size))
    allocate(current_state%praut_tot(global_grid_z_size))
    allocate(current_state%prevp_tot(global_grid_z_size))
    allocate(current_state%psaci_tot(global_grid_z_size))
    allocate(current_state%psacr_tot(global_grid_z_size))
    allocate(current_state%psacw_tot(global_grid_z_size))
    allocate(current_state%psaut_tot(global_grid_z_size))
    allocate(current_state%psdep_tot(global_grid_z_size))
    allocate(current_state%psedg_tot(global_grid_z_size))
    allocate(current_state%psedi_tot(global_grid_z_size))
    allocate(current_state%psedl_tot(global_grid_z_size))
    allocate(current_state%psedr_tot(global_grid_z_size))
    allocate(current_state%pseds_tot(global_grid_z_size))
    allocate(current_state%psmlt_tot(global_grid_z_size))
    allocate(current_state%pssub_tot(global_grid_z_size))
    allocate(current_state%qr_tot(global_grid_z_size))
    allocate(current_state%dqr_subs_profile_diag(global_grid_z_size))
    allocate(current_state%rh_tot(global_grid_z_size))
    allocate(current_state%global_grid%configuration%vertical%rho(global_grid_z_size))
    allocate(current_state%global_grid%configuration%vertical%rhon(global_grid_z_size))
    allocate(current_state%richardson_number_tot(global_grid_z_size))
    allocate(current_state%richardson_squared_tot(global_grid_z_size))
    allocate(current_state%dqs_subs_profile_diag(global_grid_z_size))
    allocate(current_state%buoysg_tot(global_grid_z_size))
    allocate(current_state%ssub_tot(global_grid_z_size))
    allocate(current_state%sed_tot(global_grid_z_size))
    allocate(current_state%tend_pr_tot_th_diff(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qv_diff(global_grid_z_size))
    allocate(current_state%tend_pr_tot_ql_diff(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qi_diff(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qr_diff(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qs_diff(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qg_diff(global_grid_z_size))
    allocate(current_state%tend_pr_tot_tabs_diff(global_grid_z_size))
    allocate(current_state%tend_pr_tot_w_buoy(global_grid_z_size))
    allocate(current_state%tend_pr_tot_u_corio(global_grid_z_size))
    allocate(current_state%tend_pr_tot_v_corio(global_grid_z_size))
    allocate(current_state%tend_pr_tot_u_forc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_v_forc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_th_forc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qv_forc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_ql_forc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qi_forc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qr_forc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qs_forc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qg_forc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_tabs_forc(global_grid_z_size))
    allocate(current_state%tendp_pr_tot_u_pt(global_grid_z_size))
    allocate(current_state%tendp_pr_tot_v_pt(global_grid_z_size))
    allocate(current_state%tendp_pr_tot_w_pt(global_grid_z_size))
    allocate(current_state%tend_pr_tot_u_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_v_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_w_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_th_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qv_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_ql_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qi_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qr_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qs_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qg_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_tabs_pwad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_th_sf(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qv_sf(global_grid_z_size))
    allocate(current_state%tend_pr_tot_ql_sf(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qi_sf(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qr_sf(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qs_sf(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qg_sf(global_grid_z_size))
    allocate(current_state%tend_pr_tot_tabs_sf(global_grid_z_size))
    allocate(current_state%tend_pr_tot_th_thad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_tabs_thad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_u_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_v_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_w_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_th_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qv_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_ql_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qi_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qr_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qs_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_qg_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_tabs_tvad(global_grid_z_size))
    allocate(current_state%tend_pr_tot_u_visc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_v_visc(global_grid_z_size))
    allocate(current_state%tend_pr_tot_w_visc(global_grid_z_size))
    allocate(current_state%th2_tot(global_grid_z_size))
    allocate(current_state%th2sg_tot(global_grid_z_size))
    allocate(current_state%dtheta_subs_profile_diag(global_grid_z_size))
    allocate(current_state%theta_dis_tot(global_grid_z_size))
    allocate(current_state%theta_tot(global_grid_z_size))
    allocate(current_state%thinit(global_grid_z_size))
    allocate(current_state%thref(global_grid_z_size))
    allocate(current_state%tkesg_tot(global_grid_z_size))
    allocate(current_state%cloud_mask_tot(global_grid_z_size))
    allocate(current_state%du_subs_profile_diag(global_grid_z_size))
    allocate(current_state%u_wind_tot(global_grid_z_size))
    allocate(current_state%uprime_tot(global_grid_z_size))
    allocate(current_state%uusg_tot(global_grid_z_size))
    allocate(current_state%uv_tot(global_grid_z_size))
    allocate(current_state%uw_tot(global_grid_z_size))
    allocate(current_state%uwsg_tot(global_grid_z_size))
    allocate(current_state%dv_subs_profile_diag(global_grid_z_size))
    allocate(current_state%v_wind_tot(global_grid_z_size))
    allocate(current_state%vprime_tot(global_grid_z_size))
    allocate(current_state%vvsg_tot(global_grid_z_size))
    allocate(current_state%vw_tot(global_grid_z_size))
    allocate(current_state%vwsg_tot(global_grid_z_size))
    allocate(current_state%qv_tot(global_grid_z_size))
    allocate(current_state%dqv_subs_profile_diag(global_grid_z_size))
    allocate(current_state%dql_subs_profile_diag(global_grid_z_size))
    allocate(current_state%vis_coef_tot(global_grid_z_size))
    allocate(current_state%w_wind_tot(global_grid_z_size))
    allocate(current_state%wke_tot(global_grid_z_size))
    allocate(current_state%wkesg_tot(global_grid_z_size))
    allocate(current_state%wtheta_ad_tot(global_grid_z_size))
    allocate(current_state%wqv_ad_tot(global_grid_z_size))
    allocate(current_state%wql_ad_tot(global_grid_z_size))
    allocate(current_state%wqr_ad_tot(global_grid_z_size))
    allocate(current_state%wqi_ad_tot(global_grid_z_size))
    allocate(current_state%wqs_ad_tot(global_grid_z_size))
    allocate(current_state%wqg_ad_tot(global_grid_z_size))
    allocate(current_state%wtheta_cn_tot(global_grid_z_size))
    allocate(current_state%wqv_cn_tot(global_grid_z_size))
    allocate(current_state%wql_cn_tot(global_grid_z_size))
    allocate(current_state%wqr_cn_tot(global_grid_z_size))
    allocate(current_state%wqi_cn_tot(global_grid_z_size))
    allocate(current_state%wqs_cn_tot(global_grid_z_size))
    allocate(current_state%wqg_cn_tot(global_grid_z_size))
    allocate(current_state%wqv_sg_tot(global_grid_z_size))
    allocate(current_state%wql_sg_tot(global_grid_z_size))
    allocate(current_state%wqr_sg_tot(global_grid_z_size))
    allocate(current_state%wqi_sg_tot(global_grid_z_size))
    allocate(current_state%wqs_sg_tot(global_grid_z_size))
    allocate(current_state%wqg_sg_tot(global_grid_z_size))
    allocate(current_state%wtsg_tot(global_grid_z_size))
    allocate(current_state%ww_tot(global_grid_z_size))
    allocate(current_state%wwsg_tot(global_grid_z_size))
    allocate(current_state%www_tot(global_grid_z_size))
    allocate(current_state%wwww_tot(global_grid_z_size))
    allocate(current_state%cloud_reff_tot(global_grid_z_size))
    allocate(current_state%longwave_hr_tot(global_grid_z_size))
    allocate(current_state%shortwave_hr_tot(global_grid_z_size))
    allocate(current_state%tend_pr_tot_th_lw(global_grid_z_size))
    allocate(current_state%tend_pr_tot_tabs_lw(global_grid_z_size))
    allocate(current_state%tend_pr_tot_th_sw(global_grid_z_size))
    allocate(current_state%tend_pr_tot_tabs_sw(global_grid_z_size))
    allocate(current_state%tend_pr_tot_th_total(global_grid_z_size))
    allocate(current_state%tend_pr_tot_tabs_total(global_grid_z_size))

    current_state%diff_coef_tot = 0.0_DEFAULT_PRECISION
    current_state%dissipation_tot = 0.0_DEFAULT_PRECISION
    current_state%dqc_mphys_tot = 0.0_DEFAULT_PRECISION
    current_state%dqg_mphys_tot = 0.0_DEFAULT_PRECISION
    current_state%dqi_mphys_tot = 0.0_DEFAULT_PRECISION
    current_state%dqr_mphys_tot = 0.0_DEFAULT_PRECISION
    current_state%dqs_mphys_tot = 0.0_DEFAULT_PRECISION
    current_state%dqv_mphys_tot = 0.0_DEFAULT_PRECISION
    current_state%dth_mphys_tot = 0.0_DEFAULT_PRECISION
    current_state%dth_cond_evap_tot = 0.0_DEFAULT_PRECISION
    current_state%dqv_cond_evap_tot = 0.0_DEFAULT_PRECISION
    current_state%qs_tot = 0.0_DEFAULT_PRECISION
    current_state%qg_tot = 0.0_DEFAULT_PRECISION
    current_state%dqg_subs_profile_diag = 0.0_DEFAULT_PRECISION
    current_state%cloud_ice_mask_tot = 0.0_DEFAULT_PRECISION
    current_state%qi_tot = 0.0_DEFAULT_PRECISION
    current_state%dqi_subs_profile_diag = 0.0_DEFAULT_PRECISION
    current_state%cloud_liq_mask_tot = 0.0_DEFAULT_PRECISION
    current_state%ql_tot = 0.0_DEFAULT_PRECISION
    current_state%pcond_tot = 0.0_DEFAULT_PRECISION
    current_state%pgacs_tot = 0.0_DEFAULT_PRECISION
    current_state%pgacw_tot = 0.0_DEFAULT_PRECISION
    current_state%pgmlt_tot = 0.0_DEFAULT_PRECISION
    current_state%pgsub_tot = 0.0_DEFAULT_PRECISION
    current_state%phomc_tot = 0.0_DEFAULT_PRECISION
    current_state%piacw_tot = 0.0_DEFAULT_PRECISION
    current_state%pidep_tot = 0.0_DEFAULT_PRECISION
    current_state%pimlt_tot = 0.0_DEFAULT_PRECISION
    current_state%pinuc_tot = 0.0_DEFAULT_PRECISION
    current_state%pisub_tot = 0.0_DEFAULT_PRECISION
    current_state%pracw_tot = 0.0_DEFAULT_PRECISION
    current_state%praut_tot = 0.0_DEFAULT_PRECISION
    current_state%prevp_tot = 0.0_DEFAULT_PRECISION
    current_state%psaci_tot = 0.0_DEFAULT_PRECISION
    current_state%psacr_tot = 0.0_DEFAULT_PRECISION
    current_state%psacw_tot = 0.0_DEFAULT_PRECISION
    current_state%psaut_tot = 0.0_DEFAULT_PRECISION
    current_state%psdep_tot = 0.0_DEFAULT_PRECISION
    current_state%psedg_tot = 0.0_DEFAULT_PRECISION
    current_state%psedi_tot = 0.0_DEFAULT_PRECISION
    current_state%psedl_tot = 0.0_DEFAULT_PRECISION
    current_state%psedr_tot = 0.0_DEFAULT_PRECISION
    current_state%pseds_tot = 0.0_DEFAULT_PRECISION
    current_state%psmlt_tot = 0.0_DEFAULT_PRECISION
    current_state%pssub_tot = 0.0_DEFAULT_PRECISION
    current_state%qr_tot = 0.0_DEFAULT_PRECISION
    current_state%dqr_subs_profile_diag = 0.0_DEFAULT_PRECISION
    current_state%rh_tot = 0.0_DEFAULT_PRECISION
    current_state%global_grid%configuration%vertical%rho = 0.0_DEFAULT_PRECISION
    current_state%global_grid%configuration%vertical%rhon = 0.0_DEFAULT_PRECISION
    current_state%richardson_number_tot = 0.0_DEFAULT_PRECISION
    current_state%richardson_squared_tot = 0.0_DEFAULT_PRECISION
    current_state%dqs_subs_profile_diag = 0.0_DEFAULT_PRECISION
    current_state%buoysg_tot = 0.0_DEFAULT_PRECISION
    current_state%ssub_tot = 0.0_DEFAULT_PRECISION
    current_state%sed_tot = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_th_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qv_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_ql_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qi_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qr_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qs_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qg_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_tabs_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_w_buoy = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_u_corio = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_v_corio = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_u_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_v_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_th_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qv_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_ql_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qi_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qr_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qs_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qg_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_tabs_forc = 0.0_DEFAULT_PRECISION
    current_state%tendp_pr_tot_u_pt = 0.0_DEFAULT_PRECISION
    current_state%tendp_pr_tot_v_pt = 0.0_DEFAULT_PRECISION
    current_state%tendp_pr_tot_w_pt = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_u_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_v_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_w_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_th_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qv_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_ql_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qi_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qr_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qs_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qg_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_tabs_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_th_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qv_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_ql_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qi_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qr_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qs_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qg_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_tabs_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_th_thad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_tabs_thad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_u_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_v_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_w_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_th_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qv_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_ql_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qi_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qr_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qs_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_qg_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_tabs_tvad = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_u_visc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_v_visc = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_w_visc = 0.0_DEFAULT_PRECISION
    current_state%th2_tot = 0.0_DEFAULT_PRECISION
    current_state%th2sg_tot = 0.0_DEFAULT_PRECISION
    current_state%dtheta_subs_profile_diag = 0.0_DEFAULT_PRECISION
    current_state%theta_dis_tot = 0.0_DEFAULT_PRECISION
    current_state%theta_tot = 0.0_DEFAULT_PRECISION
    current_state%thinit = 0.0_DEFAULT_PRECISION
    current_state%thref = 0.0_DEFAULT_PRECISION
    current_state%tkesg_tot = 0.0_DEFAULT_PRECISION
    current_state%cloud_mask_tot = 0.0_DEFAULT_PRECISION
    current_state%du_subs_profile_diag = 0.0_DEFAULT_PRECISION
    current_state%u_wind_tot = 0.0_DEFAULT_PRECISION
    current_state%uprime_tot = 0.0_DEFAULT_PRECISION
    current_state%uusg_tot = 0.0_DEFAULT_PRECISION
    current_state%uv_tot = 0.0_DEFAULT_PRECISION
    current_state%uw_tot = 0.0_DEFAULT_PRECISION
    current_state%uwsg_tot = 0.0_DEFAULT_PRECISION
    current_state%dv_subs_profile_diag = 0.0_DEFAULT_PRECISION
    current_state%v_wind_tot = 0.0_DEFAULT_PRECISION
    current_state%vprime_tot = 0.0_DEFAULT_PRECISION
    current_state%vvsg_tot = 0.0_DEFAULT_PRECISION
    current_state%vw_tot = 0.0_DEFAULT_PRECISION
    current_state%vwsg_tot = 0.0_DEFAULT_PRECISION
    current_state%qv_tot = 0.0_DEFAULT_PRECISION
    current_state%dqv_subs_profile_diag = 0.0_DEFAULT_PRECISION
    current_state%dql_subs_profile_diag = 0.0_DEFAULT_PRECISION
    current_state%vis_coef_tot = 0.0_DEFAULT_PRECISION
    current_state%w_wind_tot = 0.0_DEFAULT_PRECISION
    current_state%wke_tot = 0.0_DEFAULT_PRECISION
    current_state%wkesg_tot = 0.0_DEFAULT_PRECISION
    current_state%wtheta_ad_tot = 0.0_DEFAULT_PRECISION
    current_state%wqv_ad_tot = 0.0_DEFAULT_PRECISION
    current_state%wql_ad_tot = 0.0_DEFAULT_PRECISION
    current_state%wqr_ad_tot = 0.0_DEFAULT_PRECISION
    current_state%wqi_ad_tot = 0.0_DEFAULT_PRECISION
    current_state%wqs_ad_tot = 0.0_DEFAULT_PRECISION
    current_state%wqg_ad_tot = 0.0_DEFAULT_PRECISION
    current_state%wqv_sg_tot = 0.0_DEFAULT_PRECISION
    current_state%wql_sg_tot = 0.0_DEFAULT_PRECISION
    current_state%wqr_sg_tot = 0.0_DEFAULT_PRECISION
    current_state%wqi_sg_tot = 0.0_DEFAULT_PRECISION
    current_state%wqs_sg_tot = 0.0_DEFAULT_PRECISION
    current_state%wqg_sg_tot = 0.0_DEFAULT_PRECISION
    current_state%wtsg_tot = 0.0_DEFAULT_PRECISION
    current_state%ww_tot = 0.0_DEFAULT_PRECISION
    current_state%wwsg_tot = 0.0_DEFAULT_PRECISION
    current_state%www_tot = 0.0_DEFAULT_PRECISION
    current_state%wwww_tot = 0.0_DEFAULT_PRECISION
    current_state%cloud_reff_tot = 0.0_DEFAULT_PRECISION
    current_state%longwave_hr_tot = 0.0_DEFAULT_PRECISION
    current_state%shortwave_hr_tot = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_th_lw = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_tabs_lw = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_th_sw = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_tabs_sw = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_th_total = 0.0_DEFAULT_PRECISION
    current_state%tend_pr_tot_tabs_total = 0.0_DEFAULT_PRECISION


    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(current_state%uwsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vwsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uusg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vvsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wwsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tkesg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wtsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th2sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%sed_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ssub_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dissipation_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%buoysg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wkesg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%theta_dis_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vis_coef_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%diff_coef_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%richardson_number_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%richardson_squared_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wqv_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wql_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%wqr_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqi_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqs_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqg_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
    end if

    if (current_state%casim_profile_dgs_enabled .eqv. .true.) then
      call mpi_recv(current_state%dqc_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqg_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqi_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqr_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqs_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqv_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dth_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dth_cond_evap_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqv_cond_evap_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%phomc_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pinuc_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pidep_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psdep_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%piacw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psacw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psacr_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pisub_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pssub_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pimlt_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psmlt_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psaut_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psaci_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%praut_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pracw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%prevp_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pgacw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pgacs_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pgmlt_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pgsub_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psedi_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pseds_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psedr_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psedg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psedl_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pcond_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%profile_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(current_state%u_wind_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uprime_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%v_wind_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vprime_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wke_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ww_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%www_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wwww_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%theta_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_wind_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%rh_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wtheta_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wtheta_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uv_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th2_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%thref, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%global_grid%configuration%vertical%rho, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%global_grid%configuration%vertical%rhon, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%thinit, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qv_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ql_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%qr_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%qi_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%qs_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%qg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%wqv_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wql_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%wqr_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqi_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqs_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqg_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%wqv_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wql_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%wqr_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqi_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqs_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqg_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%cloud_mask_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%cloud_liq_mask_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%cloud_ice_mask_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
    end if

    if (current_state%forcing_enabled .eqv. .true.) then
      call mpi_recv(current_state%du_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dv_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dtheta_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqv_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dql_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqr_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqi_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqs_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqg_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_u_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qi_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qr_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qs_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qg_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%diffusion_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_th_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_pr_tot_qi_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qr_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qs_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qg_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_pr_tot_tabs_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%buoyancy_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_w_buoy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%coriolis_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_u_corio, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_corio, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%pstep_enabled .eqv. .true.) then
      call mpi_recv(current_state%tendp_pr_tot_u_pt, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tendp_pr_tot_v_pt, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tendp_pr_tot_w_pt, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%pw_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_u_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_w_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qi_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qr_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qs_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qg_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%stepfields_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_th_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_pr_tot_qi_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qr_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qs_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qg_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_pr_tot_tabs_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%th_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_th_thad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_thad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_u_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_w_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_pr_tot_qi_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qr_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qs_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qg_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_pr_tot_tabs_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%viscosity_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_u_visc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_visc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_w_visc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%socrates_enabled .eqv. .true.) then
      call mpi_recv(current_state%cloud_reff_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%longwave_hr_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%shortwave_hr_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_lw, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_lw, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_sw, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_sw, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_total, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_total, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ls1 = len_trim(diagnostic_path)
      ls2 = 0
      do i = 1,ls1
          if(diagnostic_path(i:i).ne.' ') then
            ls2 = ls2 + 1
          endif
      enddo
    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = diagnostic_path(:ls2)//"/full_diag_1d_time_"//trim(conv_to_string(int(current_state%time)))//".nc"
    else
      unique_filename = diagnostic_path(:ls2)//"/full_diag_1d_timestep_"//trim(conv_to_string(current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "t", 1, time_id))
    call check(nf90_def_dim(ncdf_id, "z_size", global_grid_z_size, z_size_id))
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))


    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "zh", "heights at w levels (m)", vertical_grid%z, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "zhn", "heights at pressure levels (m)", vertical_grid%zn, "m", current_state%time)
    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uwsg_mean", "uwsg_mean", &
                                          current_state%uwsg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vwsg_tot", "vwsg_mean", &
                                          current_state%vwsg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uusg_mean", "uusg_mean", &
                                          current_state%uusg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vvsg_tot", "vvsg_mean", &
                                          current_state%vvsg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wwsg_mean", "wwsg_mean", &
                                          current_state%wwsg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tkesg_mean", "tkesg_mean", &
                                          current_state%tkesg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wtsg_mean", "wtsg_mean", &
                                          current_state%wtsg_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th2sg_mean", "th2sg_mean", &
                                          current_state%th2sg_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "sed_mean", "mean subgrid turbulent transport", &
                                          current_state%sed_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                         "ssub_mean", "mean subgrid shear stress", &
                                          current_state%ssub_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dissipation_mean", "dissipation", &
                                          current_state%dissipation_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "buoysg_mean", "mean subgrid buoyant production", &
                                          current_state%buoysg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wkesg_mean", "wkesg_mean", &
                                          current_state%wkesg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "theta_dis_mean", "theta_dis_mean", current_state%theta_dis_tot, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vis_coef_tot", "viscosity_coef_mean", &
                                          current_state%vis_coef_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "diff_coef_mean", "mean diffusion coefficient",&
                                          current_state%diff_coef_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                   "richardson_number_mean", "mean Richardson number", &
                                    current_state%richardson_number_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "richardson_squared_mean", "mean squared Richardson number", &
                                      current_state%richardson_squared_tot/surface, "K2", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqv_sg_mean", "wqv_sg_mean", &
                                            current_state%wqv_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wql_sg_mean", "wql_sg_mean", &
                                            current_state%wql_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqr_sg_mean", "wqr_sg_mean", &
                                            current_state%wqr_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqi_sg_mean", "wqi_sg_mean", &
                                            current_state%wqi_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqs_sg_mean", "wqs_sg_mean", &
                                            current_state%wqs_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqg_sg_mean", "wqg_sg_mean", &
                                            current_state%wqg_sg_tot/surface, "K.m.s-1", current_state%time)
      end if
    end if

    if (current_state%casim_profile_dgs_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqc_mphys_mean", "mean qv rate", &
                                          current_state%dqc_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqg_mphys_mean", "mean qg rate", &
                                          current_state%dqg_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqi_mphys_mean", "mean qi rate", &
                                          current_state%dqi_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqr_mphys_mean", "mean qr rate", &
                                          current_state%dqr_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqs_mphys_mean", "mean qs rate", &
                                          current_state%dqs_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqv_mphys_mean", " mean qv rate", &
                                          current_state%dqv_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dth_mphys_mean", "mean th rate", &
                                          current_state%dth_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dth_cond_evap_mean", "mean condensation/evaporation rate due to th", &
                                          current_state%dth_cond_evap_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqv_cond_evap_mean", "mean condensation/evaporation rate due to qv", &
                                          current_state%dqv_cond_evap_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "phomc_mean", "mean mass homogeneous freezing of cloud droplet rate", &
                                          current_state%phomc_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pinuc_mean", "mean mass ice nucleation rate", &
                                          current_state%pinuc_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pidep_mean", "mean mass ice deposition rate", &
                                          current_state%pidep_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psdep_mean", "mean mass snow deposition rate", &
                                          current_state%psdep_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "piacw_mean", "mean mass ice -> cloud -> ice accretion rate", &
                                          current_state%piacw_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psacw_mean", "mean mass snow -> cloud -> snow accretion rate", &
                                          current_state%psacw_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                            "psacr_mean", "mean mass snow -> rain -> graupel and snow -> rain -> snow accretion rate", &
                                          current_state%psacr_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pisub_mean", "mean mass ice sublimation rate", &
                                          current_state%pisub_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pssub_mean", "mean mass snow sublimation rate", &
                                          current_state%pssub_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pimlt_mean", "mean mass ice melting rate", &
                                          current_state%pimlt_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psmlt_mean", "mean mass snow melting rate", &
                                          current_state%psmlt_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psaut_mean", "mean mass autoconversion to snow rate", &
                                          current_state%psaut_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psaci_mean", "mean mass snow -> ice -> snow accretion rate", &
                                          current_state%psaci_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "praut_mean", "mean mass autoconversion to rain rate", &
                                          current_state%praut_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pracw_mean", "mean mass rain accreting cloud rate", &
                                          current_state%pracw_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "prevp_mean", "mean mass evaporation of rain rate", &
                                          current_state%prevp_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pgacw_mean", "mean mass graupel -> cloud -> graupel accretion rate", &
                                          current_state%pgacw_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pgacs_mean", "mean mass graupel -> snow accretion rate", &
                                          current_state%pgacs_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pgmlt_mean", "mean mass graupel melting rate", &
                                          current_state%pgmlt_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pgsub_mean", "mean mass graupel sublimation rate", &
                                          current_state%pgsub_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psedi_mean", "mean mass ice sedimentation rate", &
                                          current_state%psedi_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pseds_mean", "mean mass snow sedimentation rate", &
                                          current_state%pseds_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psedr_mean", "mean mass rain sedimentation rate", &
                                          current_state%psedr_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psedg_mean", "mean mass graupel sedimentation rate", &
                                          current_state%psedg_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psedl_mean", "mean mass liquid/cloud sedimentation rate", &
                                          current_state%psedl_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pcond_mean", "mean mass condensation rate", &
                                          current_state%pcond_tot/surface, "kg.kg-1.s-1", current_state%time)
    end if

    if (current_state%profile_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "u_wind_mean", "u_wind_mean", &
                                          current_state%u_wind_tot/surface, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uu_mean", "uu_mean", &
                                          current_state%uprime_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "v_wind_tot", "v_wind_mean", &
                                          current_state%v_wind_tot/surface, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vprime_tot", "vv_mean", &
                                          current_state%vprime_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wke_mean", "wke_mean", &
                                          current_state%wke_tot/surface, "J.m-2.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ww_mean", "ww_mean", &
                                          current_state%ww_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "www_mean", "www_mean", &
                                          current_state%www_tot/surface, "m3.s-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wwww_mean", "wwww_mean", &
                                          current_state%wwww_tot/surface, "m4.s-4", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "theta_mean", "theta_mean", &
                                          current_state%theta_tot/surface, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "w_wind_mean", "w_wind_mean", &
                                          current_state%w_wind_tot/surface, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "rh_mean", "mean Relative Humidity", &
                                          current_state%rh_tot/surface, "%", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wtheta_ad_mean", "wtheta_ad_mean", &
                                          current_state%wtheta_ad_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wtheta_cn_mean", "wtheta_cn_mean", &
                                          current_state%wtheta_cn_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uw_mean", "uw_mean", &
                                          current_state%uw_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vw_tot", "vw_mean", &
                                          current_state%vw_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uv_mean", "uv_mean", &
                                          current_state%uv_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th2_mean", "th2_mean", &
                                          current_state%th2_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "thref", "thref", current_state%thref, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "rho", "density", &
                                          current_state%global_grid%configuration%vertical%rho/surface, &
                                          "kg.m-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "rhon", "density (altitude following pressure)", &
                                          current_state%global_grid%configuration%vertical%rhon/surface, &
                                          "kg.m-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "thinit", "thinit", current_state%thinit, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qv_tot", "mean vapour mass", &
                                          current_state%qv_tot/surface, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ql_mean", "mean liquid mass", current_state%ql_tot/surface, &
                                          "kg.kg-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "qr_mean", "mean rain mass", current_state%qr_tot/surface, &
                                            "kg.kg-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "qi_mean", "mean ice mass", &
                                            current_state%qi_tot/surface, "kg.kg-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "qs_mean", "snow_mmr_mean", current_state%qs_tot/surface, &
                                            "kg.kg-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "qg_mean", "graupel_mmr_mean", current_state%qg_tot/surface, &
                                            "kg.kg-1", current_state%time)
      end if
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wqv_cn_mean", "wqv_cn_mean", &
                                          current_state%wqv_cn_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wql_cn_mean", "wql_cn_mean", &
                                          current_state%wql_cn_tot/surface, "K.m.s-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqr_cn_mean", "wqr_cn_mean", &
                                            current_state%wqr_cn_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqi_cn_mean", "wqi_cn_mean", &
                                            current_state%wqi_cn_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqs_cn_mean", "wqs_cn_mean", &
                                            current_state%wqs_cn_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqg_cn_mean", "wqg_cn_mean", &
                                            current_state%wqg_cn_tot/surface, "K.m.s-1", current_state%time)
      end if
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wqv_ad_mean", "wqv_ad_mean", &
                                          current_state%wqv_ad_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wql_ad_mean", "wql_ad_mean", &
                                          current_state%wql_ad_tot/surface, "K.m.s-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqr_ad_mean", "wqr_ad_mean", &
                                            current_state%wqr_ad_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqi_ad_mean", "wqi_ad_mean", &
                                            current_state%wqi_ad_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqs_ad_mean", "wqs_ad_mean", &
                                            current_state%wqs_ad_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqg_ad_mean", "wqg_ad_mean", &
                                            current_state%wqg_ad_tot/surface, "K.m.s-1", current_state%time)

        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "cloud_mask_tot", "total_cloud_fraction", &
                                            current_state%cloud_mask_tot/surface, "no unit", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "cloud_liq_mask_mean", "liquid_cloud_fraction", &
                                            current_state%cloud_liq_mask_tot/surface, "no unit", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "cloud_ice_mask_mean", "ice_cloud_fraction", &
                                            current_state%cloud_ice_mask_tot/surface, "no unit", current_state%time)
      end if
    end if

    if (current_state%forcing_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "du_subs_profile_diag", "u_subsidence_mean", current_state%du_subs_profile_diag, &
                                          "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dv_subs_profile_diag", "v_subsidence_mean", current_state%dv_subs_profile_diag, &
                                          "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                            "dtheta_subs_profile_diag", "th_subsidence_mean", current_state%dtheta_subs_profile_diag, &
                            "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqv_subs_profile_diag", "vapour_mmr_subsidence_mean", &
                                          current_state%dqv_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dql_subs_profile_diag", "liquid_mmr_subsidence_mean", &
                                          current_state%dql_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqr_subs_profile_diag", "mean qr rate (subsidence)", &
                                    current_state%dqr_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqi_subs_profile_diag", "mean qi rate (subsidence)",&
                                    current_state%dqi_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqs_subs_profile_diag", "mean qs rate (subsidence)", &
                                      current_state%dqs_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqg_subs_profile_diag", "dqg_subs_profile_diag", &
                                    current_state%dqg_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_forc_mean", "mean tendancy u from forcing", &
                                      current_state%tend_pr_tot_u_forc/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_forc_mean", "mean tendancy v from forcing", &
                                      current_state%tend_pr_tot_v_forc/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_forc_mean", "mean tendancy th from forcing", &
                                      current_state%tend_pr_tot_th_forc/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_forc_mean", "mean tendancy qv from forcing", &
                                      current_state%tend_pr_tot_qv_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_forc_mean", "mean tendancy ql from forcing", &
                                      current_state%tend_pr_tot_ql_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qi_forc_mean", "mean tendancy qi from forcing", &
                                      current_state%tend_pr_tot_qi_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qr_forc_mean", "mean tendancy qr from forcing", &
                                      current_state%tend_pr_tot_qr_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qs_forc_mean", "mean tendancy qs from forcing", &
                                      current_state%tend_pr_tot_qs_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qg_forc_mean", "mean tendancy qg from forcing", &
                                      current_state%tend_pr_tot_qg_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "tend_tabs_forc_mean", "mean tendancy tabs from forcing", &
                                    current_state%tend_pr_tot_tabs_forc/surface, "???", current_state%time)
    end if

    if (current_state%diffusion_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_th_diff_mean", "mean tendancy th from diffusion", &
                                          current_state%tend_pr_tot_th_diff/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_diff_mean", "mean tendancy qv from diffusion", &
                                      current_state%tend_pr_tot_qv_diff/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_diff_mean", "mean tendancy ql from diffusion", &
                                      current_state%tend_pr_tot_ql_diff/surface, "kg.kg-1.s-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qi_diff_mean", "mean tendancy qi from diffusion", &
                                        current_state%tend_pr_tot_qi_diff/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qr_diff_mean", "mean tendancy qr from diffusion", &
                                        current_state%tend_pr_tot_qr_diff/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qs_diff_mean", "mean tendancy qs from diffusion", &
                                        current_state%tend_pr_tot_qs_diff/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qg_diff_mean", "mean tendancy qg from diffusion", &
                                        current_state%tend_pr_tot_qg_diff/surface, "kg.kg-1.s-1", current_state%time)
      end if
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_diff_mean", "mean tendancy tabs from diffusion", &
                                      current_state%tend_pr_tot_tabs_diff/surface, "???", current_state%time)
    end if
    if (current_state%buoyancy_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_pr_tot_w_buoy", "mean tendancy w from bouyancy", &
                                      current_state%tend_pr_tot_w_buoy/surface, "m.s-2", current_state%time)
    end if
    if (current_state%coriolis_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_corio_mean", "mean tendancy u from coriolis", &
                                      current_state%tend_pr_tot_u_corio/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_corio_mean", "mean tendancy v from coriolis", &
                                      current_state%tend_pr_tot_v_corio/surface, "m.s-2", current_state%time)
    end if
    if (current_state%pstep_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_pt_mean", "mean tendancy u from pstep", &
                                      current_state%tendp_pr_tot_u_pt/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_pt_mean", "mean tendancy v from pstep", &
                                      current_state%tendp_pr_tot_v_pt/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_w_pt_mean", "mean tendancy w from pstep", &
                                      current_state%tendp_pr_tot_w_pt/surface, "m.s-2", current_state%time)
    end if
    if (current_state%pw_advection_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_pwad_mean", "mean tendancy u from pwadvection", &
                                      current_state%tend_pr_tot_u_pwad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_pwad_mean", "mean tendancy v from pwadvection", &
                                      current_state%tend_pr_tot_v_pwad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_w_pwad_mean", "mean tendancy w from pwadvection", &
                                      current_state%tend_pr_tot_w_pwad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_pwad_mean", "mean tendancy th from pwadvection", &
                                      current_state%tend_pr_tot_th_pwad/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_pwad_mean", "mean tendancy qv from pwadvection", &
                                      current_state%tend_pr_tot_qv_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_pwad_mean", "mean tendancy ql from pwadvection", &
                                      current_state%tend_pr_tot_ql_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qi_pwad_mean", "mean tendancy qi from pwadvection", &
                                      current_state%tend_pr_tot_qi_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qr_pwad_mean", "mean tendancy qr from pwadvection", &
                                      current_state%tend_pr_tot_qr_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qs_pwad_mean", "mean tendancy qs from pwadvection", &
                                      current_state%tend_pr_tot_qs_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qg_pwad_mean", "mean tendancy qg from pwadvection", &
                                      current_state%tend_pr_tot_qg_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_pwad_mean", "mean tendancy tabs from pwadvection", &
                                      current_state%tend_pr_tot_tabs_pwad/surface, "???", current_state%time)
    end if
    if (current_state%stepfields_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_sf_mean", "mean tendancy th from stepfield", &
                                      current_state%tend_pr_tot_th_sf/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_sf_mean", "mean tendancy qv from stepfield", &
                                      current_state%tend_pr_tot_qv_sf/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_sf_mean", "mean tendancy ql from stepfield", &
                                      current_state%tend_pr_tot_ql_sf/surface, "kg.kg-1.s-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qi_sf_mean", "mean tendancy qi from stepfield", &
                                        current_state%tend_pr_tot_qi_sf/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qr_sf_mean", "mean tendancy qr from stepfield", &
                                        current_state%tend_pr_tot_qr_sf/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qs_sf_mean", "mean tendancy qs from stepfield", &
                                        current_state%tend_pr_tot_qs_sf/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qg_sf_mean", "mean tendancy qg from stepfield", &
                                        current_state%tend_pr_tot_qg_sf/surface, "kg.kg-1.s-1", current_state%time)
      end if
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_sf_mean", "mean tendancy tabs from stepfield", &
                                      current_state%tend_pr_tot_tabs_sf/surface, "???", current_state%time)
    end if
    if (current_state%th_advection_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_thad_mean", "mean tendancy th from thadvection", &
                                      current_state%tend_pr_tot_th_thad/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_thad_mean", "mean tendancy tabs from thadvection", &
                                      current_state%tend_pr_tot_tabs_thad/surface, "???", current_state%time)
    end if
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_tvad_mean", "mean tendancy u from tvdadvection", &
                                      current_state%tend_pr_tot_u_tvad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_tvad_mean", "mean tendancy v from tvdadvection", &
                                      current_state%tend_pr_tot_v_tvad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_w_tvad_mean", "mean tendancy w from tvdadvection", &
                                      current_state%tend_pr_tot_w_tvad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_tvad_mean", "mean tendancy th from tvdadvection", &
                                      current_state%tend_pr_tot_th_tvad/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_tvad_mean", "mean tendancy qv from tvdadvection", &
                                      current_state%tend_pr_tot_qv_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_tvad_mean", "mean tendancy ql from tvdadvection", &
                                      current_state%tend_pr_tot_ql_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qi_tvad_mean", "mean tendancy qi from tvdadvection", &
                                      current_state%tend_pr_tot_qi_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qr_tvad_mean", "mean tendancy qr from tvdadvection", &
                                      current_state%tend_pr_tot_qr_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qs_tvad_mean", "mean tendancy qs from tvdadvection", &
                                      current_state%tend_pr_tot_qs_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qg_tvad_mean", "mean tendancy qg from tvdadvection", &
                                      current_state%tend_pr_tot_qs_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_tvad_mean", "mean tendancy tabs from tvdadvection", &
                                      current_state%tend_pr_tot_tabs_tvad/surface, "???", current_state%time)
    end if
    if (current_state%viscosity_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_u_visc_mean", "mean tendancy u from viscosity", &
                                          current_state%tend_pr_tot_u_visc/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_v_visc_mean", "mean tendancy v from viscosity", &
                                          current_state%tend_pr_tot_v_visc/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_w_visc_mean", "mean tendancy w from viscosity", &
                                          current_state%tend_pr_tot_w_visc/surface, "m.s-2", current_state%time)
    end if
    if (current_state%socrates_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "cloud_reff_mean", "mean cloud reflectivity", &
                                          current_state%cloud_reff_tot/surface, " ", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "longwave_hr_mean", "mean lw heating rate", &
                                          current_state%longwave_hr_tot/surface, "K.d-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "shortwave_hr_mean", "mean sw heating rate", &
                                          current_state%shortwave_hr_tot/surface, "K.d-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_th_lw", "mean tendency th from lw heating rate", &
                                          current_state%tend_pr_tot_th_lw/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_tabs_lw", "mean tendency tabs from lw heating rate", &
                                          current_state%tend_pr_tot_tabs_lw/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_th_sw", "mean tendency th from sw heating rate", &
                                          current_state%tend_pr_tot_th_sw/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_tabs_sw", "mean tendency tabs from sw heating rate", &
                                          current_state%tend_pr_tot_tabs_sw/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_th_total", "mean tendency th from total heating rate", &
                                          current_state%tend_pr_tot_th_total/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_tabs_total", "mean tendency th from tabs heating rate", &
                                          current_state%tend_pr_tot_tabs_total/surface, "K.s-1", current_state%time)
      call check(nf90_close(ncdf_id))
    end if


    deallocate(current_state%dql_subs_profile_diag, current_state%qg_tot, current_state%diff_coef_tot, &
                current_state%dqc_mphys_tot, current_state%dqg_mphys_tot, current_state%dqi_mphys_tot, &
                current_state%dqr_mphys_tot, current_state%dqs_mphys_tot, current_state%dqv_mphys_tot, &
                current_state%dth_mphys_tot, current_state%dth_cond_evap_tot, current_state%dqv_cond_evap_tot, &
                current_state%qs_tot, current_state%dqg_subs_profile_diag, current_state%cloud_ice_mask_tot, &
                current_state%qi_tot, current_state%dqi_subs_profile_diag, current_state%cloud_liq_mask_tot, &
                current_state%ql_tot, current_state%pcond_tot, current_state%pgacs_tot, current_state%pgacw_tot, &
                current_state%pgmlt_tot, current_state%pgsub_tot, current_state%phomc_tot, current_state%piacw_tot, &
                current_state%pidep_tot, current_state%pimlt_tot, current_state%pinuc_tot, current_state%pisub_tot, &
                current_state%pracw_tot, current_state%praut_tot, current_state%prevp_tot, current_state%psaci_tot, &
                current_state%psacr_tot, current_state%psacw_tot, current_state%psaut_tot, current_state%psdep_tot, &
                current_state%psedg_tot, current_state%psedi_tot, current_state%psedl_tot, current_state%psedr_tot, &
                current_state%pseds_tot, current_state%psmlt_tot, current_state%pssub_tot, current_state%qr_tot, &
                current_state%dqr_subs_profile_diag, current_state%rh_tot, current_state%global_grid%configuration%vertical%rho, &
                current_state%global_grid%configuration%vertical%rhon, current_state%richardson_number_tot, &
                current_state%richardson_squared_tot, current_state%dqs_subs_profile_diag, current_state%buoysg_tot, &
                current_state%ssub_tot, current_state%sed_tot, current_state%tend_pr_tot_th_diff, &
                current_state%tend_pr_tot_qv_diff, current_state%tend_pr_tot_ql_diff, current_state%tend_pr_tot_qi_diff, &
                current_state%tend_pr_tot_qr_diff, current_state%tend_pr_tot_qs_diff, current_state%tend_pr_tot_qg_diff, &
                current_state%tend_pr_tot_tabs_diff, current_state%tend_pr_tot_w_buoy, current_state%tend_pr_tot_u_corio, &
                current_state%tend_pr_tot_v_corio, current_state%tend_pr_tot_u_forc, current_state%tend_pr_tot_v_forc, &
                current_state%tend_pr_tot_th_forc, current_state%tend_pr_tot_qv_forc, current_state%tend_pr_tot_ql_forc, &
                current_state%tend_pr_tot_qi_forc, current_state%tend_pr_tot_qr_forc, current_state%tend_pr_tot_qs_forc, &
                current_state%tend_pr_tot_qg_forc, current_state%tend_pr_tot_tabs_forc, current_state%tendp_pr_tot_u_pt, &
                current_state%tendp_pr_tot_v_pt, current_state%tendp_pr_tot_w_pt, current_state%tend_pr_tot_u_pwad, &
                current_state%tend_pr_tot_v_pwad, current_state%tend_pr_tot_w_pwad, current_state%tend_pr_tot_th_pwad, &
                current_state%tend_pr_tot_qv_pwad, current_state%tend_pr_tot_ql_pwad, current_state%tend_pr_tot_qi_pwad, &
                current_state%tend_pr_tot_qr_pwad, current_state%tend_pr_tot_qs_pwad, current_state%tend_pr_tot_qg_pwad, &
                current_state%tend_pr_tot_tabs_pwad, current_state%tend_pr_tot_th_sf, current_state%tend_pr_tot_qv_sf, &
                current_state%tend_pr_tot_ql_sf, current_state%tend_pr_tot_qi_sf, current_state%tend_pr_tot_qr_sf, &
                current_state%tend_pr_tot_qs_sf, current_state%tend_pr_tot_qg_sf, current_state%tend_pr_tot_tabs_sf, &
                current_state%tend_pr_tot_th_thad, current_state%tend_pr_tot_tabs_thad, current_state%tend_pr_tot_u_tvad, &
                current_state%tend_pr_tot_v_tvad, current_state%tend_pr_tot_w_tvad, current_state%tend_pr_tot_th_tvad, &
                current_state%tend_pr_tot_qv_tvad, current_state%tend_pr_tot_ql_tvad, current_state%tend_pr_tot_qi_tvad, &
                current_state%tend_pr_tot_qr_tvad, current_state%tend_pr_tot_qs_tvad, current_state%tend_pr_tot_qg_tvad, &
                current_state%tend_pr_tot_tabs_tvad, current_state%tend_pr_tot_u_visc, current_state%tend_pr_tot_v_visc, &
                current_state%tend_pr_tot_w_visc, current_state%th2_tot, current_state%th2sg_tot, &
                current_state%dtheta_subs_profile_diag, current_state%theta_dis_tot, current_state%theta_tot, &
                current_state%thinit, current_state%thref, current_state%tkesg_tot, current_state%cloud_mask_tot, &
                current_state%du_subs_profile_diag, current_state%u_wind_tot, current_state%uprime_tot, current_state%uusg_tot, &
                current_state%uv_tot, current_state%uw_tot, current_state%uwsg_tot, current_state%dv_subs_profile_diag, &
                current_state%vw_tot, current_state%vvsg_tot, current_state%qv_tot, current_state%dqv_subs_profile_diag, &
                current_state%vis_coef_tot, current_state%v_wind_tot, current_state%vprime_tot, current_state%vwsg_tot, &
                current_state%w_wind_tot, current_state%wke_tot, current_state%wkesg_tot, current_state%wtheta_ad_tot, &
                current_state%wqv_ad_tot, current_state%wql_ad_tot, current_state%wqr_ad_tot, current_state%wqi_ad_tot, &
                current_state%wqs_ad_tot, current_state%wqg_ad_tot, current_state%wtheta_cn_tot, current_state%wqv_cn_tot, &
                current_state%wql_cn_tot, current_state%wqr_cn_tot, current_state%wqi_cn_tot, current_state%wqs_cn_tot, &
                current_state%wqg_cn_tot, current_state%wqv_sg_tot, current_state%wql_sg_tot, current_state%wqr_sg_tot, &
                current_state%wqi_sg_tot, current_state%wqs_sg_tot, current_state%wqg_sg_tot, current_state%wtsg_tot, &
                current_state%ww_tot, current_state%wwsg_tot, current_state%www_tot, current_state%wwww_tot, &
                current_state%dissipation_tot, current_state%cloud_reff_tot, current_state%longwave_hr_tot, &
                current_state%shortwave_hr_tot, current_state%tend_pr_tot_th_lw, current_state%tend_pr_tot_tabs_lw, &
                current_state%tend_pr_tot_th_sw, current_state%tend_pr_tot_tabs_sw, &
                current_state%tend_pr_tot_th_total, current_state%tend_pr_tot_tabs_total)
    !print*,"send_data_for_diag_1D MONC"
  end subroutine diagnostic_file_1d_generation

  subroutine diagnostic_file_2d_generation(current_state, vertical_grid, io_communicator_arg, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           io_communicator_arg

    integer :: z_size_id, y_size_id, x_size_id, scalar_size_id, time_id
    integer :: ncdf_id, ierr
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename

    allocate(current_state%qlmax(global_grid_y_size, global_grid_x_size))
    allocate(current_state%hqlmax(global_grid_y_size, global_grid_x_size))
    allocate(current_state%cltop(global_grid_y_size, global_grid_x_size))
    allocate(current_state%clbas(global_grid_y_size, global_grid_x_size))
    allocate(current_state%vwp(global_grid_y_size, global_grid_x_size))
    allocate(current_state%lwp(global_grid_y_size, global_grid_x_size))
    allocate(current_state%wmax(global_grid_y_size, global_grid_x_size))
    allocate(current_state%wmin(global_grid_y_size, global_grid_x_size))
    allocate(current_state%reske(global_grid_y_size, global_grid_x_size))
    allocate(current_state%senhf(global_grid_y_size, global_grid_x_size))
    allocate(current_state%lathf(global_grid_y_size, global_grid_x_size))
    allocate(current_state%rwp(global_grid_y_size, global_grid_x_size))
    allocate(current_state%iwp(global_grid_y_size, global_grid_x_size))
    allocate(current_state%swp(global_grid_y_size, global_grid_x_size))
    allocate(current_state%gwp(global_grid_y_size, global_grid_x_size))
    allocate(current_state%tot_iwp(global_grid_y_size, global_grid_x_size))
    allocate(current_state%subke_2d(global_grid_y_size, global_grid_x_size))
    allocate(current_state%surface_precip(global_grid_y_size, global_grid_x_size))
    allocate(current_state%surface_cloudsed(global_grid_y_size, global_grid_x_size))
    allocate(current_state%surface_rainsed(global_grid_y_size, global_grid_x_size))
    current_state%qlmax = 0.0_DEFAULT_PRECISION
    current_state%hqlmax = 0.0_DEFAULT_PRECISION
    current_state%cltop = 0.0_DEFAULT_PRECISION
    current_state%clbas = 0.0_DEFAULT_PRECISION
    current_state%vwp = 0.0_DEFAULT_PRECISION
    current_state%lwp = 0.0_DEFAULT_PRECISION
    current_state%wmax = 0.0_DEFAULT_PRECISION
    current_state%wmin = 0.0_DEFAULT_PRECISION
    current_state%reske = 0.0_DEFAULT_PRECISION
    current_state%senhf = 0.0_DEFAULT_PRECISION
    current_state%lathf = 0.0_DEFAULT_PRECISION
    current_state%rwp = 0.0_DEFAULT_PRECISION
    current_state%iwp = 0.0_DEFAULT_PRECISION
    current_state%swp = 0.0_DEFAULT_PRECISION
    current_state%gwp = 0.0_DEFAULT_PRECISION
    current_state%tot_iwp = 0.0_DEFAULT_PRECISION
    current_state%surface_precip = 0.0_DEFAULT_PRECISION
    current_state%surface_cloudsed = 0.0_DEFAULT_PRECISION
    current_state%surface_rainsed = 0.0_DEFAULT_PRECISION
    current_state%subke_2d = 0.0_DEFAULT_PRECISION
    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(current_state%qlmax, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%hqlmax, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%cltop, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%clbas, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%lwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wmax, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wmin, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%reske, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%senhf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%lathf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%rwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%iwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%swp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%gwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tot_iwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
    end if
    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(current_state%subke_2d, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%casim_enabled .eqv. .true.) then
      call mpi_recv(current_state%surface_precip, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%surface_cloudsed, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%surface_rainsed, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ls1 = len_trim(diagnostic_path)
    ls2 = 0
    do i = 1,ls1
        if(diagnostic_path(i:i).ne.' ') then
          ls2 = ls2 + 1
        endif
    enddo
    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = diagnostic_path(:ls2)//"/full_diag_2d_time_"//trim(conv_to_string(int(current_state%time)))//".nc"
    else
      unique_filename = diagnostic_path(:ls2)//"/full_diag_2d_timestep_"//trim(conv_to_string(current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "t", 1, time_id))
    call check(nf90_def_dim(ncdf_id, "z", global_grid_z_size, z_size_id))
    call check(nf90_def_dim(ncdf_id, "y", global_grid_y_size, y_size_id))
    call check(nf90_def_dim(ncdf_id, "x", global_grid_x_size, x_size_id))
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))


    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)
    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "zh", "heights at w levels (m)", vertical_grid%z, "m", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "zhn", "heights at pressure levels (m)", vertical_grid%zn, "m", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "qlmax", "maximum liquid water content in a column", current_state%qlmax, &
                                        "kg.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                      "hqlmax", "height of the maximum liquid water content in a column", current_state%hqlmax, &
                      "m", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "cltop", "cloud top height", current_state%cltop, "m.s-1", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "clbas", "cloud base height", current_state%clbas, "m.s-1", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "vwp", "water vapour path for each column", current_state%vwp, "kg.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "lwp", "liquid water path for each column", current_state%lwp, "kg.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "wmax", "maximum vertical velocity for each column", current_state%wmax, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "wmin", "minimum vertical velocity for each column", current_state%wmin, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "reske", "resolved ke", current_state%reske, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "senhf", "surface sensible heat flux", current_state%senhf, "W.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "lathf", "surface latent heat flux", current_state%lathf, "W.m-2", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "rwp", "rain water path for each column", current_state%rwp, "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "iwp", "ice water path for each column", current_state%iwp, "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "swp", "snow water path for each column", current_state%swp, "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "gwp", "graupel water path for each column", &
                                          current_state%gwp, "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                        "tot_iwp", "total ice water path (iwp + swp + gwp) for each column", current_state%tot_iwp, &
                        "kg.m-2", current_state%time)
      end if
    end if
    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "subke_2d", "subke_2d", current_state%subke_2d, "J or kg.m2.s-2", current_state%time)
    end if
    if (current_state%casim_enabled .eqv. .true.) then
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "surface_precip", "surface_precip", current_state%surface_precip, &
                                        "???", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "surface_cloudsed", "surface_cloudsed", current_state%surface_cloudsed, &
                                        "???", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "surface_rainsed", "surface_rainsed", current_state%surface_rainsed, &
                                        "???", current_state%time)
    end if

    call check(nf90_close(ncdf_id))

    deallocate(current_state%qlmax, current_state%hqlmax, current_state%cltop, current_state%clbas, current_state%vwp, &
                current_state%lwp, current_state%wmax, current_state%wmin, current_state%reske, current_state%senhf, &
                current_state%lathf, current_state%rwp, current_state%iwp, current_state%swp, current_state%gwp, &
                current_state%tot_iwp, current_state%surface_precip, current_state%surface_cloudsed, &
                current_state%surface_rainsed, current_state%subke_2d)
    !print*,"send_data_for_diag_2D MONC"
  end subroutine diagnostic_file_2d_generation

  subroutine diagnostic_file_3d_generation(current_state, vertical_grid, io_communicator_arg, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           io_communicator_arg

    integer :: z_size_id, y_size_id, x_size_id, scalar_size_id, time_id
    integer :: ncdf_id, ierr
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename

    allocate(current_state%u%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%v%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%w%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%th%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%p%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qv%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ql%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qi%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qs%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qg%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAitkenSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAccumSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAccumInsolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qCoarseSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qCoarseDustMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nl%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ni%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ns%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ng%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAitkenSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAccumSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAccumInsolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nCoarseSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nCoarseDustnumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag bouyancy
    allocate(current_state%tend_3d_w_buoy(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag coriolis
    allocate(current_state%tend_3d_u_corio(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_v_corio(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag diffusion
    allocate(current_state%tend_3d_th_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qv_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_ql_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qr_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qi_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qs_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qg_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_tabs_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag forcing
    allocate(current_state%tend_3d_u_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_v_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_th_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qv_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_ql_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qr_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qi_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qs_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qg_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_tabs_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag pstep
    allocate(current_state%tendp_3d_u_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tendp_3d_v_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tendp_3d_w_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_u_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_v_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_w_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag pwadvection
    allocate(current_state%tend_3d_u_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_v_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_w_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_th_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qv_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_ql_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qr_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qi_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qs_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qg_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_tabs_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag stepfields
    allocate(current_state%tend_3d_th_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qv_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_ql_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qr_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qi_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qs_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qg_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_tabs_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag thadvection
    allocate(current_state%tend_3d_th_thad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_tabs_thad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag tvadvection
    allocate(current_state%tend_3d_u_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_v_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_w_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_th_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qv_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_ql_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qr_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qi_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qs_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_qg_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_tabs_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag viscosity
    allocate(current_state%tend_3d_u_visc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_v_visc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_w_visc(global_grid_z_size, global_grid_y_size, global_grid_x_size))

    allocate(current_state%rdAitkenSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdAccumSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdCoarseSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_cloud%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_rain%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_ice%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_snow%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_graupel%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%RH%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%RI%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))

    ! SOCRATES
    allocate(current_state%cloud_reff%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%lwrad_hr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%swrad_hr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_tabs_lw(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_tabs_sw(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%tend_3d_tabs_total(global_grid_z_size, global_grid_y_size, global_grid_x_size))

    current_state%u%data = 0.0_DEFAULT_PRECISION
    current_state%v%data = 0.0_DEFAULT_PRECISION
    current_state%w%data = 0.0_DEFAULT_PRECISION
    current_state%th%data = 0.0_DEFAULT_PRECISION
    current_state%p%data = 0.0_DEFAULT_PRECISION
    current_state%qv%data = 0.0_DEFAULT_PRECISION
    current_state%ql%data = 0.0_DEFAULT_PRECISION
    current_state%qr%data = 0.0_DEFAULT_PRECISION
    current_state%qi%data = 0.0_DEFAULT_PRECISION
    current_state%qs%data = 0.0_DEFAULT_PRECISION
    current_state%qg%data = 0.0_DEFAULT_PRECISION
    current_state%qAitkenSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qAccumSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qAccumInsolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qCoarseSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qCoarseDustMass%data = 0.0_DEFAULT_PRECISION
    current_state%nl%data = 0.0_DEFAULT_PRECISION
    current_state%nr%data = 0.0_DEFAULT_PRECISION
    current_state%ni%data = 0.0_DEFAULT_PRECISION
    current_state%ns%data = 0.0_DEFAULT_PRECISION
    current_state%ng%data = 0.0_DEFAULT_PRECISION
    current_state%nAitkenSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nAccumSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nAccumInsolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nCoarseSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nCoarseDustnumber%data = 0.0_DEFAULT_PRECISION
    ! diag bouyancy
    current_state%tend_3d_w_buoy = 0.0_DEFAULT_PRECISION
    ! diag coriolis
    current_state%tend_3d_u_corio = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_v_corio = 0.0_DEFAULT_PRECISION
    ! diag diffusion
    current_state%tend_3d_th_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qv_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_ql_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qr_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qi_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qs_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qg_diff = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_tabs_diff = 0.0_DEFAULT_PRECISION
    ! diag forcing
    current_state%tend_3d_u_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_v_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_th_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qv_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_ql_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qr_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qi_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qs_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qg_forc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_tabs_forc = 0.0_DEFAULT_PRECISION
    ! diag pstep
    current_state%tendp_3d_u_pt = 0.0_DEFAULT_PRECISION
    current_state%tendp_3d_v_pt = 0.0_DEFAULT_PRECISION
    current_state%tendp_3d_w_pt = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_u_pt = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_v_pt = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_w_pt = 0.0_DEFAULT_PRECISION
    ! diag pwadvection
    current_state%tend_3d_u_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_v_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_w_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_th_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qv_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_ql_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qr_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qi_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qs_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qg_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_tabs_pwad = 0.0_DEFAULT_PRECISION
    ! diag stepfields
    current_state%tend_3d_th_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qv_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_ql_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qr_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qi_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qs_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qg_sf = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_tabs_sf = 0.0_DEFAULT_PRECISION
    ! diag thadvection
    current_state%tend_3d_th_thad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_tabs_thad = 0.0_DEFAULT_PRECISION
    ! diag tvadvection
    current_state%tend_3d_u_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_v_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_w_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_th_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qv_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_ql_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qr_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qi_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qs_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_qg_pwad = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_tabs_pwad = 0.0_DEFAULT_PRECISION
    ! diag viscosity
    current_state%tend_3d_u_visc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_v_visc = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_w_visc = 0.0_DEFAULT_PRECISION

    current_state%rdAitkenSol%data = 0.0_DEFAULT_PRECISION
    current_state%rdAccumSol%data = 0.0_DEFAULT_PRECISION
    current_state%rdCoarseSol%data = 0.0_DEFAULT_PRECISION
    current_state%D0_cloud%data = 0.0_DEFAULT_PRECISION
    current_state%D0_rain%data = 0.0_DEFAULT_PRECISION
    current_state%D0_ice%data = 0.0_DEFAULT_PRECISION
    current_state%D0_snow%data = 0.0_DEFAULT_PRECISION
    current_state%D0_graupel%data = 0.0_DEFAULT_PRECISION
    current_state%RH%data = 0.0_DEFAULT_PRECISION
    current_state%RI%data = 0.0_DEFAULT_PRECISION

    ! SOCRATES
    current_state%cloud_reff%data = 0.0_DEFAULT_PRECISION
    current_state%lwrad_hr%data = 0.0_DEFAULT_PRECISION
    current_state%swrad_hr%data = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_tabs_lw = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_tabs_sw = 0.0_DEFAULT_PRECISION
    current_state%tend_3d_tabs_total = 0.0_DEFAULT_PRECISION

    call mpi_recv(current_state%u%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%v%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%w%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%th%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%p%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qv%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ql%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qi%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qs%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qg%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qAitkenSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qAccumSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qAccumInsolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qCoarseSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qCoarseDustMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nl%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ni%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ns%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ng%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nAitkenSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nAccumSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nAccumInsolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nCoarseSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nCoarseDustnumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    ! diag bouyancy
    if (current_state%buoyancy_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_w_buoy, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag coriolis
    if (current_state%coriolis_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_corio, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_corio, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag diffusion
    if (current_state%diffusion_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_th_diff, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_diff, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_diff, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_3d_qr_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qi_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qs_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qg_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_3d_tabs_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag forcing
    if (current_state%forcing_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_th_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qr_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qi_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qs_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qg_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag pstep
    if (current_state%pstep_enabled .eqv. .true.) then
      call mpi_recv(current_state%tendp_3d_u_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tendp_3d_v_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tendp_3d_w_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_u_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_w_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag pwadvection
    if (current_state%pw_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_w_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_th_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qr_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qi_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qs_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qg_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag stepfields
    if (current_state%stepfields_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_th_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_3d_qr_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qi_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qs_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qg_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_3d_tabs_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag thadvection
    if (current_state%th_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_th_thad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_thad, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag tvadvection
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_w_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_th_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_3d_qr_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qi_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qs_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qg_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_3d_tabs_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag viscosity
    if (current_state%viscosity_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_visc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_visc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_w_visc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    call mpi_recv(current_state%rdAitkenSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%rdAccumSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%rdCoarseSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_cloud%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_rain%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_ice%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_snow%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_graupel%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%RH%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%RI%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

    ! SOCRATES
    if (current_state%socrates_enabled .eqv. .true.) then
      call mpi_recv(current_state%cloud_reff%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%lwrad_hr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%swrad_hr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_lw, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_sw, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_total, global_array_size*2, MPI_REAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ls1 = len_trim(diagnostic_path)
    ls2 = 0
    do i = 1,ls1
        if(diagnostic_path(i:i).ne.' ') then
          ls2 = ls2 + 1
        endif
    enddo
    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = diagnostic_path(:ls2)//"/full_diag_3d_time_"//trim(conv_to_string(int(current_state%time)))//".nc"
    else
      unique_filename = diagnostic_path(:ls2)//"/full_diag_3d_timestep_"//trim(conv_to_string(current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "t", 1, time_id))
    call check(nf90_def_dim(ncdf_id, "z", global_grid_z_size, z_size_id))
    call check(nf90_def_dim(ncdf_id, "y", global_grid_y_size, y_size_id))
    call check(nf90_def_dim(ncdf_id, "x", global_grid_x_size, x_size_id))
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))


    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "zh", "heights at w levels", vertical_grid%z, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "zhn", "heights at pressure levels", vertical_grid%zn, "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_bottom", "x bottom coordinate", current_state%global_grid%bottom(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_top", "x top coordinate", current_state%global_grid%top(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_resolution", "resolution along x axis", current_state%global_grid%resolution(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_bottom", "y bottom coordinate", current_state%global_grid%bottom(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_top", "y top coordinate", current_state%global_grid%top(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_resolution", "resolution along y axis", current_state%global_grid%resolution(2), "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "u", "u wind component", current_state%u%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                     "v", "v wind component", current_state%v%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "w", "w wind component", current_state%w%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "th", "potentiel temperature perturbations", current_state%th%data, "K", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "p", "pressure", current_state%p%data, "Pa", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qv", "water vapour mass mixing ratio", current_state%qv%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ql", "liquid water mass mixing ratio", current_state%ql%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qr", "rain water mass mixing ratio", current_state%qr%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qi", "ice water mass mixing ratio", current_state%qi%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qs", "snow water mass mixing ratio", current_state%qs%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qg", "graupel water mass mixing ratio", current_state%qg%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAitkenSolMass", "Aitken soluble aerosol mass mixing ratio", &
                                      current_state%qAitkenSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAccumSolMass", "accumulation soluble aerosol mass mixing ratio",&
                                      current_state%qAccumSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAccumInsolMass", "accumulation insoluble aerosol mass mixing ratio",&
                                      current_state%qAccumInsolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qCoarseSolMass", "coarse soluble aerosol mixing ratio",&
                                      current_state%qCoarseSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qCoarseDustMass", "dust coarse aerosol mixing ratio",&
                                      current_state%qCoarseDustMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nl", "liquid water number concentration", current_state%nl%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nr", "rain water number concentration", current_state%nr%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ni", "ice water number concentration", current_state%ni%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ns", "snow water number concentration", current_state%ns%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ng", "graupel water number concentration", current_state%ng%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAitkenSolNumber", "Aitken soluble aerosol number concentration", &
                                      current_state%nAitkenSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAccumSolNumber", "accumulation soluble aerosol number concentration",&
                                      current_state%nAccumSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAccumInsolNumber", "accumulation insoluble aerosol number concentration",&
                                      current_state%nAccumInsolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nCoarseSolNumber", "coarse soluble aerosol number concentration",&
                                      current_state%nCoarseSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nCoarseDustnumber", "dust coarse aerosol number concentration",&
                                      current_state%nCoarseDustnumber%data, "kg-1", current_state%time)
    ! diag bouyancy
    if (current_state%buoyancy_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_buoy", "tend_3d_w_buoy", current_state%tend_3d_w_buoy, &
                                        "m.s-1", current_state%time)
    end if
    ! diag coriolis
    if (current_state%coriolis_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_corio", "3d tendency u from coriolis", current_state%tend_3d_u_corio, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_corio", "3d tendency v from coriolis", current_state%tend_3d_v_corio, &
                                        "m.s-1", current_state%time)
    end if
    ! diag diffusion
    if (current_state%diffusion_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_diff", "3d tendency th from diffusion", current_state%tend_3d_th_diff, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_diff", "3d tendency qv from diffusion", current_state%tend_3d_qv_diff, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_diff", "3d tendency ql from diffusion", current_state%tend_3d_ql_diff, &
                                        "kg.kg-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qr_diff", "3d tendency qr from diffusion", current_state%tend_3d_qr_diff, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qi_diff", "3d tendency qi from diffusion", current_state%tend_3d_qi_diff, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qs_diff", "3d tendency qs from diffusion", current_state%tend_3d_qs_diff, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qg_diff", "3d tendency qg from diffusion", current_state%tend_3d_qg_diff, &
                                          "kg.kg-1", current_state%time)
      end if
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_diff", "3d tendency tabs from diffusion", &
                                        current_state%tend_3d_tabs_diff, "???", current_state%time)
    end if
    ! diag forcing
    if (current_state%forcing_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_forc", "3d tendency u from forcing", current_state%tend_3d_u_forc, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_forc", "3d tendency v from forcing", current_state%tend_3d_v_forc, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_forc", "3d tendency th from forcing", current_state%tend_3d_th_forc, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_forc", "3d tendency qv from forcing", current_state%tend_3d_qv_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_forc", "3d tendency ql from forcing", current_state%tend_3d_ql_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qr_forc", "3d tendency qr from forcing", current_state%tend_3d_qr_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qi_forc", "3d tendency qi from forcing", current_state%tend_3d_qi_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qs_forc", "3d tendency qs from forcing", current_state%tend_3d_qs_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qg_forc", "3d tendency qg from forcing", current_state%tend_3d_qg_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                "tend_3d_tabs_forc", "3d tendency tabs from forcing", current_state%tend_3d_tabs_forc, &
                                "???", current_state%time)
    end if
    ! diag pstep
    if (current_state%pstep_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tendp_3d_u_pt", "3d previous tendency u from pstep", current_state%tendp_3d_u_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tendp_3d_v_pt", "3d previous tendency v from pstep", current_state%tendp_3d_v_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tendp_3d_w_pt", "3d previous tendency w from pstep", current_state%tendp_3d_w_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_pt", "3d tendency u from pstep", current_state%tend_3d_u_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_pt", "3d tendency v from pstep", current_state%tend_3d_v_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_pt", "3d tendency w from pstep", current_state%tend_3d_w_pt, &
                                        "m.s-1", current_state%time)
    end if
    ! diag pwadvection
    if (current_state%pw_advection_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_pwad", "3d tendency u from pwadvection", current_state%tend_3d_u_pwad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_pwad", "3d tendency v from pwadvection", current_state%tend_3d_v_pwad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_pwad", "3d tendency w from pwadvection", current_state%tend_3d_w_pwad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_pwad", "3d tendency th from pwadvection", current_state%tend_3d_th_pwad, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_pwad", "3d tendency qv from pwadvection", &
                                        current_state%tend_3d_qv_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_pwad", "3d tendency ql from pwadvection", &
                                        current_state%tend_3d_ql_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qr_pwad", "3d tendency qr from pwadvection", &
                                        current_state%tend_3d_qr_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qi_pwad", "3d tendency qi from pwadvection", &
                                        current_state%tend_3d_qi_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qs_pwad", "3d tendency qs from pwadvection", &
                                        current_state%tend_3d_qs_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qg_pwad", "3d tendency qg from pwadvection", &
                                        current_state%tend_3d_qg_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_pwad", "3d tendency tabs from pwadvection", &
                                        current_state%tend_3d_tabs_pwad, "???", current_state%time)
    end if
    ! diag stepfields
    if (current_state%stepfields_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_sf", "3d tendency th from stepfield", current_state%tend_3d_th_sf, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_sf", "3d tendency qv from stepfield", current_state%tend_3d_qv_sf, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_sf", "3d tendency ql from stepfield", current_state%tend_3d_ql_sf, &
                                        "kg.kg-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qr_sf", "3d tendency qr from stepfield", current_state%tend_3d_qr_sf, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qi_sf", "3d tendency qi from stepfield", current_state%tend_3d_qi_sf, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qs_sf", "3d tendency qs from stepfield", current_state%tend_3d_qs_sf, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qg_sf", "3d tendency qg from stepfield", current_state%tend_3d_qg_sf, &
                                          "kg.kg-1", current_state%time)
      end if
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_sf", "3d tendency tabs from stepfield", current_state%tend_3d_tabs_sf, &
                                        "???", current_state%time)
    end if
    ! diag thadvection
    if (current_state%th_advection_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_thad", "3d tendency th from thadvection", current_state%tend_3d_th_thad, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_thad", "3d tendency th tabs thadvection", &
                                        current_state%tend_3d_tabs_thad, "???", current_state%time)
    end if
    ! diag tvadvection
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_tvad", "3d tendency u from tvdadvection", current_state%tend_3d_u_tvad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_tvad", "3d tendency v from tvdadvection", current_state%tend_3d_v_tvad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_tvad", "3d tendency w from tvdadvection", current_state%tend_3d_w_tvad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_tvad", "3d tendency th from tvdadvection", current_state%tend_3d_th_tvad, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_tvad", "3d tendency qv from tvdadvection", &
                                        current_state%tend_3d_qv_tvad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_tvad", "3d tendency ql from tvdadvection", &
                                        current_state%tend_3d_ql_tvad, "kg.kg-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qr_tvad", "3d tendency qr from tvdadvection", &
                                          current_state%tend_3d_qr_tvad, "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qi_tvad", "3d tendency qi from tvdadvection", &
                                          current_state%tend_3d_qi_tvad, "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qs_tvad", "3d tendency qs from tvdadvection", &
                                          current_state%tend_3d_qs_tvad, "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qg_tvad", "3d tendency qg from tvdadvection", &
                                          current_state%tend_3d_qg_tvad, "kg.kg-1", current_state%time)
      end if
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_tvad", "3d tendency tabs from tvdadvection", &
                                        current_state%tend_3d_tabs_tvad, "???", current_state%time)
    end if
    ! diag viscosity
    if (current_state%viscosity_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_visc", "3d tendency u from viscosity", current_state%tend_3d_u_visc, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_visc", "3d tendency v from viscosity", current_state%tend_3d_v_visc, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_visc", "3d tendency w from viscosity", current_state%tend_3d_w_visc, &
                                        "m.s-1", current_state%time)
    end if
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdAitkenSol", "Aitken soluble aerosol radius", &
                                      current_state%rdAitkenSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdAccumSol", "accumSol soluble aerosol radius", &
                                      current_state%rdAccumSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdCoarseSol", "coarseSol soluble aerosol radius", &
                                      current_state%rdCoarseSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_cloud", "cloud droplet diameter", &
                                      current_state%D0_cloud%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_rain", "rain droplet diameter", &
                                      current_state%D0_rain%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_ice", "ice particle diameter", &
                                      current_state%D0_ice%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_snow", "snow particle diameter", &
                                      current_state%D0_snow%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_graupel", "graupel particle diameter", &
                                      current_state%D0_graupel%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "RH", "relative Humidity", &
                                      current_state%RH%data, "%", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "RHI", "relative Himidity with respect to Ice", &
                                      current_state%RI%data, "%", current_state%time)

    ! SOCRATES cloud_reff, lwrad_hr, swrad_hr
    if (current_state%socrates_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "cloud_reff", "cloud reflectivity", &
                                        current_state%cloud_reff%data, " ", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "lwrad_hr", "longwave heating rate", &
                                        current_state%lwrad_hr%data, "K.d-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "swrad_hr", "shortwave heating rate", &
                                        current_state%swrad_hr%data, "K.d-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_lw", "3d tendency th from Socrates (lw)", &
                                        current_state%tend_3d_tabs_lw, "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_sw", "3d tendency th from Socrates (sw)", &
                                        current_state%tend_3d_tabs_sw, "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_total", "3d tendency th from Socrates", &
                                        current_state%tend_3d_tabs_total, "K", current_state%time)
      call check(nf90_close(ncdf_id))
    end if

    deallocate(current_state%u%data, current_state%v%data, current_state%w%data, current_state%th%data, &
              current_state%p%data, current_state%qv%data, current_state%ql%data, current_state%qr%data, &
              current_state%qi%data, current_state%qs%data, current_state%qg%data, &
              current_state%qAitkenSolMass%data, current_state%qAccumSolMass%data, current_state%qAccumInsolMass%data, &
              current_state%qCoarseSolMass%data, current_state%qCoarseDustMass%data,&
              current_state%nl%data, current_state%nr%data, current_state%ni%data, &
              current_state%ns%data, current_state%ng%data, &
              current_state%nAitkenSolNumber%data, current_state%nAccumSolNumber%data, current_state%nAccumInsolNumber%data, &
              current_state%nCoarseSolNumber%data, current_state%nCoarseDustnumber%data,&
              current_state%tend_3d_w_buoy, current_state%tend_3d_u_corio, current_state%tend_3d_v_corio, &
              current_state%tend_3d_th_diff, current_state%tend_3d_qv_diff, current_state%tend_3d_ql_diff, &
              current_state%tend_3d_qr_diff, current_state%tend_3d_qi_diff, current_state%tend_3d_qs_diff, &
              current_state%tend_3d_qg_diff, current_state%tend_3d_tabs_diff, &
              current_state%tend_3d_u_forc, current_state%tend_3d_v_forc, current_state%tend_3d_th_forc, &
              current_state%tend_3d_qv_forc, current_state%tend_3d_ql_forc, current_state%tend_3d_qr_forc, &
              current_state%tend_3d_qi_forc, current_state%tend_3d_qs_forc, current_state%tend_3d_qg_forc, &
              current_state%tend_3d_tabs_forc, current_state%tendp_3d_u_pt, current_state%tendp_3d_v_pt,&
              current_state%tendp_3d_w_pt, current_state%tend_3d_u_pt, current_state%tend_3d_v_pt, current_state%tend_3d_w_pt, &
              current_state%tend_3d_u_pwad, current_state%tend_3d_v_pwad, current_state%tend_3d_w_pwad, &
              current_state%tend_3d_th_pwad, current_state%tend_3d_qv_pwad, current_state%tend_3d_ql_pwad, &
              current_state%tend_3d_qr_pwad, current_state%tend_3d_qi_pwad, current_state%tend_3d_qs_pwad,&
              current_state%tend_3d_qg_pwad, current_state%tend_3d_tabs_pwad, &
              current_state%tend_3d_th_sf, current_state%tend_3d_qv_sf, current_state%tend_3d_ql_sf, &
              current_state%tend_3d_qr_sf, current_state%tend_3d_qi_sf, current_state%tend_3d_qs_sf, &
              current_state%tend_3d_qg_sf, current_state%tend_3d_tabs_sf, &
              current_state%tend_3d_th_thad, current_state%tend_3d_tabs_thad, &
              current_state%tend_3d_u_tvad, current_state%tend_3d_v_tvad, current_state%tend_3d_w_tvad, &
              current_state%tend_3d_th_tvad, current_state%tend_3d_qv_tvad, current_state%tend_3d_ql_tvad, &
              current_state%tend_3d_qr_tvad, current_state%tend_3d_qi_tvad, current_state%tend_3d_qs_tvad,&
              current_state%tend_3d_qg_tvad, current_state%tend_3d_tabs_tvad, &
              current_state%tend_3d_u_visc, current_state%tend_3d_v_visc, current_state%tend_3d_w_visc, &
              current_state%rdAitkenSol%data, current_state%rdAccumSol%data, current_state%rdCoarseSol%data, &
              current_state%D0_cloud%data, current_state%D0_rain%data, current_state%D0_ice%data, &
              current_state%D0_snow%data, current_state%D0_graupel%data, &
              current_state%RH%data, current_state%RI%data, &
              current_state%cloud_reff%data, current_state%lwrad_hr%data, current_state%swrad_hr%data, &
              current_state%tend_3d_tabs_lw, current_state%tend_3d_tabs_sw, current_state%tend_3d_tabs_total)


  end subroutine diagnostic_file_3d_generation

  subroutine checkpoint_file_generation(current_state, vertical_grid, io_communicator_arg, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size, time_model)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           io_communicator_arg
    real(kind=DEFAULT_PRECISION), intent(out):: time_model

    integer :: z_size_id, y_size_id, x_size_id, scalar_size_id, time_id
    integer :: ncdf_id, ierr
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename

    call mpi_recv(current_state%dtm_new, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%absolute_new_dtm, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ugal, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%vgal, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    time_model = current_state%time! + current_state%dtm

    allocate(current_state%global_grid%configuration%vertical%thref(global_grid_z_size))
    allocate(vertical_grid%prefn(global_grid_z_size))
    allocate(vertical_grid%olubar(global_grid_z_size))
    allocate(vertical_grid%olvbar(global_grid_z_size))
    allocate(vertical_grid%olthbar(global_grid_z_size))
    allocate(vertical_grid%olqvbar(global_grid_z_size))
    allocate(vertical_grid%olqlbar(global_grid_z_size))
    allocate(vertical_grid%olqrbar(global_grid_z_size))
    allocate(vertical_grid%olqibar(global_grid_z_size))
    allocate(vertical_grid%olqsbar(global_grid_z_size))
    allocate(vertical_grid%olqgbar(global_grid_z_size))
    allocate(vertical_grid%olqAitkenSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olqAccumSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olqAccumInsolMassbar(global_grid_z_size))
    allocate(vertical_grid%olqCoarseSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olqCoarseDustMassbar(global_grid_z_size))
    allocate(vertical_grid%olnlbar(global_grid_z_size))
    allocate(vertical_grid%olnrbar(global_grid_z_size))
    allocate(vertical_grid%olnibar(global_grid_z_size))
    allocate(vertical_grid%olnsbar(global_grid_z_size))
    allocate(vertical_grid%olngbar(global_grid_z_size))
    allocate(vertical_grid%olnAitkenSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olnAccumSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olnAccumInsolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olnCoarseSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olnCoarseDustnumberbar(global_grid_z_size))
    allocate(vertical_grid%olzubar(global_grid_z_size))
    allocate(vertical_grid%olzvbar(global_grid_z_size))
    allocate(vertical_grid%olzthbar(global_grid_z_size))
    allocate(vertical_grid%olzqvbar(global_grid_z_size))
    allocate(vertical_grid%olzqlbar(global_grid_z_size))
    allocate(vertical_grid%olzqrbar(global_grid_z_size))
    allocate(vertical_grid%olzqibar(global_grid_z_size))
    allocate(vertical_grid%olzqsbar(global_grid_z_size))
    allocate(vertical_grid%olzqgbar(global_grid_z_size))
    allocate(vertical_grid%olzqAitkenSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olzqAccumSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olzqAccumInsolMassbar(global_grid_z_size))
    allocate(vertical_grid%olzqCoarseSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olzqCoarseDustMassbar(global_grid_z_size))
    allocate(vertical_grid%olznlbar(global_grid_z_size))
    allocate(vertical_grid%olznrbar(global_grid_z_size))
    allocate(vertical_grid%olznibar(global_grid_z_size))
    allocate(vertical_grid%olznsbar(global_grid_z_size))
    allocate(vertical_grid%olzngbar(global_grid_z_size))
    allocate(vertical_grid%olznAitkenSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olznAccumSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olznAccumInsolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olznCoarseSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olznCoarseDustnumberbar(global_grid_z_size))
    current_state%global_grid%configuration%vertical%thref = 0.0_DEFAULT_PRECISION
    vertical_grid%prefn = 0.0_DEFAULT_PRECISION
    vertical_grid%olubar = 0.0_DEFAULT_PRECISION
    vertical_grid%olvbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olthbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqvbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqlbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqrbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqibar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqsbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqgbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqAitkenSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqAccumSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqAccumInsolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqCoarseSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqCoarseDustMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnlbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnrbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnibar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnsbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olngbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnAitkenSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnAccumSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnAccumInsolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnCoarseSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnCoarseDustnumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzubar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzvbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzthbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqvbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqlbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqrbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqibar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqsbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqgbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqAitkenSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqAccumSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqAccumInsolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqCoarseSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqCoarseDustMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznlbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznrbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznibar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznsbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzngbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznAitkenSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznAccumSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznAccumInsolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznCoarseSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznCoarseDustnumberbar = 0.0_DEFAULT_PRECISION
    call mpi_recv(current_state%global_grid%configuration%vertical%thref, &
        global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(vertical_grid%prefn, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(vertical_grid%olubar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(vertical_grid%olvbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%th%active) then
      call mpi_recv(vertical_grid%olthbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%number_q_fields .ne. 0) then
      call mpi_recv(vertical_grid%olqvbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqlbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqrbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqibar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqsbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqgbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqAitkenSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqAccumSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqAccumInsolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqCoarseSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqCoarseDustMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnlbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnrbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnibar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnsbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olngbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnAitkenSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnAccumSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnAccumInsolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnCoarseSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnCoarseDustnumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    call mpi_recv(vertical_grid%olzubar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(vertical_grid%olzvbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%th%active) then
      call mpi_recv(vertical_grid%olzthbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%number_q_fields .ne. 0) then
      call mpi_recv(vertical_grid%olzqvbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqlbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqrbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqibar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqsbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqgbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqAitkenSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqAccumSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqAccumInsolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqCoarseSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqCoarseDustMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznlbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznrbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznibar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznsbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzngbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznAitkenSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznAccumSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznAccumInsolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznCoarseSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznCoarseDustnumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if


    allocate(current_state%u%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%v%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%w%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%th%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%p%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qv%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ql%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qi%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qs%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qg%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAitkenSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAccumSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAccumInsolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qCoarseSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qCoarseDustMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nl%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ni%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ns%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ng%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAitkenSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAccumSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAccumInsolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nCoarseSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nCoarseDustnumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zu%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zv%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zw%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zth%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqv%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zql%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqi%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqs%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqg%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqAitkenSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqAccumSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqAccumInsolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqCoarseSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqCoarseDustMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znl%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zni%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zns%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zng%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znAitkenSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znAccumSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znAccumInsolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znCoarseSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znCoarseDustnumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdAitkenSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdAccumSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdCoarseSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_cloud%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_rain%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_ice%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_snow%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_graupel%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%RH%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%RI%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    current_state%u%data = 0.0_DEFAULT_PRECISION
    current_state%v%data = 0.0_DEFAULT_PRECISION
    current_state%w%data = 0.0_DEFAULT_PRECISION
    current_state%th%data = 0.0_DEFAULT_PRECISION
    current_state%p%data = 0.0_DEFAULT_PRECISION
    current_state%qv%data = 0.0_DEFAULT_PRECISION
    current_state%ql%data = 0.0_DEFAULT_PRECISION
    current_state%qr%data = 0.0_DEFAULT_PRECISION
    current_state%qi%data = 0.0_DEFAULT_PRECISION
    current_state%qs%data = 0.0_DEFAULT_PRECISION
    current_state%qg%data = 0.0_DEFAULT_PRECISION
    current_state%qAitkenSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qAccumSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qAccumInsolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qCoarseSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qCoarseDustMass%data = 0.0_DEFAULT_PRECISION
    current_state%nl%data = 0.0_DEFAULT_PRECISION
    current_state%nr%data = 0.0_DEFAULT_PRECISION
    current_state%ni%data = 0.0_DEFAULT_PRECISION
    current_state%ns%data = 0.0_DEFAULT_PRECISION
    current_state%ng%data = 0.0_DEFAULT_PRECISION
    current_state%nAitkenSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nAccumSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nAccumInsolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nCoarseSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nCoarseDustnumber%data = 0.0_DEFAULT_PRECISION
    current_state%zu%data = 0.0_DEFAULT_PRECISION
    current_state%zv%data = 0.0_DEFAULT_PRECISION
    current_state%zw%data = 0.0_DEFAULT_PRECISION
    current_state%zth%data = 0.0_DEFAULT_PRECISION
    current_state%zqv%data = 0.0_DEFAULT_PRECISION
    current_state%zql%data = 0.0_DEFAULT_PRECISION
    current_state%zqr%data = 0.0_DEFAULT_PRECISION
    current_state%zqi%data = 0.0_DEFAULT_PRECISION
    current_state%zqs%data = 0.0_DEFAULT_PRECISION
    current_state%zqg%data = 0.0_DEFAULT_PRECISION
    current_state%zqAitkenSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%zqAccumSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%zqAccumInsolMass%data = 0.0_DEFAULT_PRECISION
    current_state%zqCoarseSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%zqCoarseDustMass%data = 0.0_DEFAULT_PRECISION
    current_state%znl%data = 0.0_DEFAULT_PRECISION
    current_state%znr%data = 0.0_DEFAULT_PRECISION
    current_state%zni%data = 0.0_DEFAULT_PRECISION
    current_state%zns%data = 0.0_DEFAULT_PRECISION
    current_state%zng%data = 0.0_DEFAULT_PRECISION
    current_state%znAitkenSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%znAccumSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%znAccumInsolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%znCoarseSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%znCoarseDustnumber%data = 0.0_DEFAULT_PRECISION
    current_state%rdAitkenSol%data = 0.0_DEFAULT_PRECISION
    current_state%rdAccumSol%data = 0.0_DEFAULT_PRECISION
    current_state%rdCoarseSol%data = 0.0_DEFAULT_PRECISION
    current_state%D0_cloud%data = 0.0_DEFAULT_PRECISION
    current_state%D0_rain%data = 0.0_DEFAULT_PRECISION
    current_state%D0_ice%data = 0.0_DEFAULT_PRECISION
    current_state%D0_snow%data = 0.0_DEFAULT_PRECISION
    current_state%D0_graupel%data = 0.0_DEFAULT_PRECISION
    current_state%RH%data = 0.0_DEFAULT_PRECISION
    current_state%RI%data = 0.0_DEFAULT_PRECISION
    call mpi_recv(current_state%u%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%v%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%w%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%th%active) then
      call mpi_recv(current_state%th%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    call mpi_recv(current_state%p%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%number_q_fields .ne. 0) then
      call mpi_recv(current_state%qv%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ql%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qi%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qs%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qg%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAitkenSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumInsolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseDustMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nl%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ni%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ns%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ng%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAitkenSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumInsolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseDustnumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    call mpi_recv(current_state%zu%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%zv%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%zw%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%th%active) then
      call mpi_recv(current_state%zth%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%number_q_fields .ne. 0) then
      call mpi_recv(current_state%zqv%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zql%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqi%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqs%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqg%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqAitkenSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqAccumSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqAccumInsolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqCoarseSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqCoarseDustMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znl%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zni%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zns%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zng%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znAitkenSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znAccumSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znAccumInsolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znCoarseSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znCoarseDustnumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    call mpi_recv(current_state%rdAitkenSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%rdAccumSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%rdCoarseSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_cloud%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_rain%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_ice%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_snow%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_graupel%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%RH%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%RI%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

    ls1 = len_trim(checkpoint_path)
    ls2 = 0
    do i = 1,ls1
        if(checkpoint_path(i:i).ne.' ') then
          ls2 = ls2 + 1
        endif
    enddo
    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = checkpoint_path(:ls2)//"/checkpt_run_time_"//trim(conv_to_string(int(current_state%time)))//".nc"
    else
      unique_filename = checkpoint_path(:ls2)//"/checkpt_run_timestep_"//trim(conv_to_string(current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "t", 1, time_id))
    call check(nf90_def_dim(ncdf_id, "z", global_grid_z_size, z_size_id))
    call check(nf90_def_dim(ncdf_id, "y", global_grid_y_size, y_size_id))
    call check(nf90_def_dim(ncdf_id, "x", global_grid_x_size, x_size_id))
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))


    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "last_cfl_timestep", "last cfl timestep", current_state%last_cfl_timestep, "no unit", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "zh", "heights at w levels", vertical_grid%z, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "zhn", "heights at pressure levels", vertical_grid%zn, "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_bottom", "x bottom coordinate", current_state%global_grid%bottom(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_top", "x top coordinate", current_state%global_grid%top(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_resolution", "resolution along x axis", current_state%global_grid%resolution(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_bottom", "y bottom coordinate", current_state%global_grid%bottom(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_top", "y top coordinate", current_state%global_grid%top(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_resolution", "resolution along y axis", current_state%global_grid%resolution(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "dtm", "time timestep", current_state%dtm, "s", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "dtm_new", "new time timestep (cfltest activated)", current_state%absolute_new_dtm, "s", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "absolute_new_dtm", "absolute time timestep", current_state%absolute_new_dtm, "s", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "ugal", "galilean u wind component", current_state%ugal, "m.s-1", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "vgal", "galilean v wind component", current_state%vgal, "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "thref", "reference potential temperature profile", &
                                      current_state%global_grid%configuration%vertical%thref, "K", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "prefn", "reference pressure profil", vertical_grid%prefn, "Pa", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olubar", "u mean along z axis", vertical_grid%olubar, "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olvbar", "v mean along z axis", vertical_grid%olvbar, "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olthbar", "th mean along z axis", vertical_grid%olthbar, "K", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqvbar", "qv mean along z axis", vertical_grid%olqvbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqlbar", "ql mean along z axis", vertical_grid%olqlbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqrbar", "qr mean along z axis", vertical_grid%olqrbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqibar", "qi mean along z axis", vertical_grid%olqibar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqsbar", "qs mean along z axis", vertical_grid%olqsbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqgbar", "qg mean along z axis", vertical_grid%olqgbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqAitkenSolMassbar", "qAitkenSolMass mean along z axis", vertical_grid%olqAitkenSolMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqAccumSolMassbar", "qAccumSolMass mean along z axis", vertical_grid%olqAccumSolMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqAccumInsolMassbar", "qAccumInsolMass mean along z axis", vertical_grid%olqAccumInsolMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqCoarseSolMassbar", "qCoarseSolMass mean along z axis", vertical_grid%olqCoarseSolMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqCoarseDustMassbar", "qCoarseDustMass mean along z axis", vertical_grid%olqCoarseDustMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olnlbar", "nl mean along z axis", vertical_grid%olnlbar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olnrbar", "nr mean along z axis", vertical_grid%olnrbar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olnibar", "ni mean along z axis", vertical_grid%olnibar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olnsbar", "ns mean along z axis", vertical_grid%olnsbar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olngbar", "ng mean along z axis", vertical_grid%olngbar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnAitkenSolNumberbar", "nAitkenSolNumber mean along z axis", vertical_grid%olnAitkenSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnAccumSolNumberbar", "nAccumSolNumber mean along z axis", vertical_grid%olnAccumSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnAccumInsolNumberbar", "nAccumInsolNumber mean along z axis", vertical_grid%olnAccumInsolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnCoarseSolNumberbar", "nCoarseSolNumber mean along z axis", vertical_grid%olnCoarseSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnCoarseDustnumberbar", "nCoarseDustnumber mean along z axis", vertical_grid%olnCoarseDustnumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzubar", "previous timestep u mean along z axis", vertical_grid%olzubar, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzvbar", "previous timestep v mean along z axis", vertical_grid%olzvbar, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzthbar", "previous timestep th mean along z axis", vertical_grid%olzthbar, &
                                      "K", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqvbar", "previous timestep qv mean along z axis", vertical_grid%olzqvbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqlbar", "previous timestep ql mean along z axis", vertical_grid%olzqlbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqrbar", "previous timestep qr mean along z axis", vertical_grid%olzqrbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqibar", "previous timestep qi mean along z axis", vertical_grid%olzqibar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqsbar", "previous timestep qs mean along z axis", vertical_grid%olzqsbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqgbar", "previous timestep qg mean along z axis", vertical_grid%olzqgbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqAitkenSolMassbar", "previous timestep qAitkenSolMass mean along z axis", &
                                                      vertical_grid%olzqAitkenSolMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqAccumSolMassbar", "previous timestep qAccumSolMass mean along z axis", &
                                                    vertical_grid%olzqAccumSolMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqAccumInsolMassbar", "previous timestep qAccumInsolMass mean along z axis", &
                                                      vertical_grid%olzqAccumInsolMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqCoarseSolMassbar", "previous timestep qCoarseSolMass mean along z axis", &
                                                     vertical_grid%olzqCoarseSolMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqCoarseDustMassbar", "previous timestep qCoarseDustMass mean along z axis", &
                                                      vertical_grid%olzqCoarseDustMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olznlbar", "previous timestep nl mean along z axis", vertical_grid%olznlbar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olznrbar", "previous timestep nr mean along z axis", vertical_grid%olznrbar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olznibar", "previous timestep ni mean along z axis", vertical_grid%olznibar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olznsbar", "previous timestep ns mean along z axis", vertical_grid%olznsbar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzngbar", "previous timestep ng mean along z axis", vertical_grid%olzngbar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id,"olznAitkenSolNumberbar", &
                    "previous timestep nAitkenSolNumber mean along z axis", vertical_grid%olznAitkenSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, "olznAccumSolNumberbar", &
                    "previous timestep nAccumSolNumber mean along z axis", vertical_grid%olznAccumSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, "olznAccumInsolNumberbar", &
                    "previous timestep nAccumInsolNumber mean along z axis", vertical_grid%olznAccumInsolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, "olznCoarseSolNumberbar", &
                    "previous timestep nCoarseSolNumber mean along z axis", vertical_grid%olznCoarseSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, "olznCoarseDustnumberbar", &
                    "previous timestep nCoarseDustnumber mean along z axis", vertical_grid%olznCoarseDustnumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "u", "u wind component", current_state%u%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "v", "v wind component", current_state%v%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "w", "w wind component", current_state%w%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "th", "potentiel temperature perturbations", current_state%th%data, "K", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "p", "pressure perturbations", current_state%p%data, "Pa", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qv", "water vapour mass mixing ratio", current_state%qv%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ql", "liquid water mass mixing ratio", current_state%ql%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qr", "rain water mass mixing ratio", current_state%qr%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qi", "ice water mass mixing ratio", current_state%qi%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qs", "snow water mass mixing ratio", current_state%qs%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qg", "graupel water mass mixing ratio", current_state%qg%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAitkenSolMass", "Aitken soluble aerosol mass mixing ratio", &
                                      current_state%qAitkenSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAccumSolMass", "accumulation soluble aerosol mass mixing ratio",&
                                      current_state%qAccumSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAccumInsolMass", "accumulation insoluble aerosol mass mixing ratio",&
                                      current_state%qAccumInsolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qCoarseSolMass", "coarse soluble aerosol mixing ratio",&
                                      current_state%qCoarseSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qCoarseDustMass", "dust coarse aerosol mixing ratio",&
                                      current_state%qCoarseDustMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nl", "liquid water number concentration", current_state%nl%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nr", "rain water number concentration", current_state%nr%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ni", "ice water number concentration", current_state%ni%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ns", "snow water number concentration", current_state%ns%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ng", "graupel water number concentration", current_state%ng%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAitkenSolNumber", "Aitken soluble aerosol number concentration", &
                                      current_state%nAitkenSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAccumSolNumber", "accumulation soluble aerosol number concentration",&
                                      current_state%nAccumSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAccumInsolNumber", "accumulation insoluble aerosol number concentration",&
                                      current_state%nAccumInsolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nCoarseSolNumber", "coarse soluble aerosol number concentration",&
                                      current_state%nCoarseSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nCoarseDustnumber", "dust coarse aerosol number concentration",&
                                      current_state%nCoarseDustnumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zu", "previous timestep u wind component", current_state%zu%data, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zv", "previous timestep v wind component", current_state%zv%data, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zw", "previous timestep w wind component", current_state%zw%data, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zth", "previous timestep potentiel temperature", current_state%zth%data, &
                                      "K", current_state%time)

    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqv", "previous timestep water vapour mass mixing ratio", current_state%zqv%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zql", "previous timestep liquid water mass mixing ratio", current_state%zql%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqr", "previous timestep rain water mass mixing ratio", current_state%zqr%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqi", "previous timestep ice water mass mixing ratio", current_state%zqi%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqs", "previous timestep snow water mass mixing ratio", current_state%zqs%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqg", "previous timestep graupel water mass mixing ratio", current_state%zqg%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqAitkenSolMass", "previous timestep Aitken soluble aerosol mass mixing ratio", &
                                      current_state%zqAitkenSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqAccumSolMass", "previous timestep accumulation soluble aerosol mass mixing ratio",&
                                      current_state%zqAccumSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqAccumInsolMass", "previous timestep accumulation insoluble aerosol mass mixing ratio",&
                                      current_state%zqAccumInsolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqCoarseSolMass", "previous timestep coarse soluble aerosol mixing ratio",&
                                      current_state%zqCoarseSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqCoarseDustMass", "previous timestep dust coarse aerosol mixing ratio",&
                                      current_state%zqCoarseDustMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znl", "previous timestep liquid water number concentration", current_state%znl%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znr", "previous timestep rain water number concentration", current_state%znr%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zni", "previous timestep ice water number concentration", current_state%zni%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zns", "previous timestep snow water number concentration", current_state%zns%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zng", "previous timestep graupel water number concentration", current_state%zng%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znAitkenSolNumber", "previous timestep Aitken soluble aerosol number concentration", &
                                      current_state%znAitkenSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znAccumSolNumber", "previous timestep accumulation soluble aerosol number concentration",&
                                      current_state%znAccumSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                  "znAccumInsolNumber", "previous timestep accumulation insoluble aerosol number concentration",&
                                  current_state%znAccumInsolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znCoarseSolNumber", "previous timestep coarse soluble aerosol number concentration",&
                                      current_state%znCoarseSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znCoarseDustnumber", "previous timestep dust coarse aerosol number concentration",&
                                      current_state%znCoarseDustnumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdAitkenSol", "Aitken soluble aerosol radius", &
                                      current_state%rdAitkenSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdAccumSol", "accumSol soluble aerosol radius", &
                                      current_state%rdAccumSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdCoarseSol", "coarseSol soluble aerosol radius", &
                                      current_state%rdCoarseSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_cloud", "cloud droplet diameter", &
                                      current_state%D0_cloud%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_rain", "rain droplet diameter", &
                                      current_state%D0_rain%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_ice", "ice particle diameter", &
                                      current_state%D0_ice%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_snow", "snow particle diameter", &
                                      current_state%D0_snow%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_graupel", "graupel particle diameter", &
                                      current_state%D0_graupel%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "RH", "relative Humidity", &
                                      current_state%RH%data, "%", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "RHI", "relative Himidity with respect to Ice", &
                                      current_state%RI%data, "%", current_state%time)

    call check(nf90_close(ncdf_id))

    deallocate( current_state%global_grid%configuration%vertical%thref, vertical_grid%prefn, &
                vertical_grid%olubar, vertical_grid%olvbar, vertical_grid%olthbar, &
                vertical_grid%olqvbar, vertical_grid%olqlbar, vertical_grid%olqrbar, &
                vertical_grid%olqibar, vertical_grid%olqsbar, vertical_grid%olqgbar, &
                vertical_grid%olqAitkenSolMassbar, vertical_grid%olqAccumSolMassbar, vertical_grid%olqAccumInsolMassbar, &
                vertical_grid%olqCoarseSolMassbar, vertical_grid%olqCoarseDustMassbar, &
                vertical_grid%olnlbar, vertical_grid%olnrbar, vertical_grid%olnibar,&
                vertical_grid%olnsbar, vertical_grid%olngbar, &
                vertical_grid%olnAitkenSolNumberbar, vertical_grid%olnAccumSolNumberbar, vertical_grid%olnAccumInsolNumberbar, &
                vertical_grid%olnCoarseSolNumberbar, vertical_grid%olnCoarseDustnumberbar, &
                vertical_grid%olzubar, vertical_grid%olzvbar, vertical_grid%olzthbar, &
                vertical_grid%olzqvbar, vertical_grid%olzqlbar, vertical_grid%olzqrbar, &
                vertical_grid%olzqibar, vertical_grid%olzqsbar, vertical_grid%olzqgbar, &
                vertical_grid%olzqAitkenSolMassbar, vertical_grid%olzqAccumSolMassbar, vertical_grid%olzqAccumInsolMassbar, &
                vertical_grid%olzqCoarseSolMassbar, vertical_grid%olzqCoarseDustMassbar, &
                vertical_grid%olznlbar, vertical_grid%olznrbar, vertical_grid%olznibar,&
                vertical_grid%olznsbar, vertical_grid%olzngbar, &
                vertical_grid%olznAitkenSolNumberbar, vertical_grid%olznAccumSolNumberbar, &
                vertical_grid%olznAccumInsolNumberbar, &
                vertical_grid%olznCoarseSolNumberbar, vertical_grid%olznCoarseDustnumberbar, &
              )
    deallocate(current_state%u%data, current_state%v%data, current_state%w%data, current_state%th%data, &
              current_state%p%data, current_state%qv%data, current_state%ql%data, current_state%qr%data, &
              current_state%qi%data, current_state%qs%data, current_state%qg%data, &
              current_state%qAitkenSolMass%data, current_state%qAccumSolMass%data, current_state%qAccumInsolMass%data, &
              current_state%qCoarseSolMass%data, current_state%qCoarseDustMass%data,&
              current_state%nl%data, current_state%nr%data, current_state%ni%data, &
              current_state%ns%data, current_state%ng%data, &
              current_state%nAitkenSolNumber%data, current_state%nAccumSolNumber%data, current_state%nAccumInsolNumber%data, &
              current_state%nCoarseSolNumber%data, current_state%nCoarseDustnumber%data,&
              current_state%zu%data, current_state%zv%data, current_state%zw%data, current_state%zth%data, &
              current_state%zqv%data, current_state%zql%data, current_state%zqr%data, &
              current_state%zqi%data, current_state%zqs%data, current_state%zqg%data, &
              current_state%zqAitkenSolMass%data, current_state%zqAccumSolMass%data, current_state%zqAccumInsolMass%data, &
              current_state%zqCoarseSolMass%data, current_state%zqCoarseDustMass%data,&
              current_state%znl%data, current_state%znr%data, current_state%zni%data, &
              current_state%zns%data, current_state%zng%data, &
              current_state%znAitkenSolNumber%data, current_state%znAccumSolNumber%data, current_state%znAccumInsolNumber%data, &
              current_state%znCoarseSolNumber%data, current_state%znCoarseDustnumber%data,&
              current_state%rdAitkenSol%data, current_state%rdAccumSol%data, current_state%rdCoarseSol%data, &
              current_state%D0_cloud%data, current_state%D0_rain%data, current_state%D0_ice%data, &
              current_state%D0_snow%data, current_state%D0_graupel%data, &
              current_state%RH%data, current_state%RI%data)
    !print*,"checkpoint_file_generation IO"
  end subroutine checkpoint_file_generation

  subroutine define_and_write_variable_integer_scalar(ncdf_id, name_var, standard_name, data_scalar, unit, time)
    integer, intent(in) :: ncdf_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    integer, intent(in) :: data_scalar

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,1) :: data_scalar_time
    data_scalar_time(1,1) = data_scalar
    call check(nf90_def_var(ncdf_id, name_var, NF90_INT, (/ time_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_scalar_time))
  end subroutine define_and_write_variable_integer_scalar

  subroutine define_and_write_variable_real_scalar(ncdf_id, name_var, standard_name, data_scalar, unit, time)
    integer, intent(in) :: ncdf_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    real(kind=DEFAULT_PRECISION), intent(in) :: data_scalar

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,1) :: data_scalar_time
    data_scalar_time(1,1) = data_scalar

    call check(nf90_def_var(ncdf_id, name_var, NF90_DOUBLE, (/ time_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_scalar_time))
  end subroutine define_and_write_variable_real_scalar

  subroutine define_and_write_variable_1D(ncdf_id, z_size_id, name_var, standard_name, data_1D, unit, time)
    integer, intent(in) :: ncdf_id, z_size_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: data_1D

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,global_grid_z_size) :: data_1D_time
    data_1D_time(1,:) = data_1D

    call check(nf90_def_var(ncdf_id, name_var, NF90_DOUBLE, (/ time_id, z_size_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_1D_time))
  end subroutine define_and_write_variable_1D

  subroutine define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, name_var, standard_name, data_2D, unit, time)
    integer, intent(in) :: ncdf_id, y_size_id, x_size_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: data_2D

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,global_grid_y_size,global_grid_x_size) :: data_2D_time
    data_2D_time(1,:,:) = data_2D

    call check(nf90_def_var(ncdf_id, name_var, NF90_DOUBLE, (/ time_id, y_size_id, x_size_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_2D_time))
  end subroutine define_and_write_variable_2D

  subroutine define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, name_var, standard_name, data_3D, unit, time)
    integer, intent(in) :: ncdf_id, z_size_id, y_size_id, x_size_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: data_3D

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,global_grid_z_size,global_grid_y_size,global_grid_x_size) :: data_3D_time
    data_3D_time(1,:,:,:) = data_3D

    call check(nf90_def_var(ncdf_id, name_var, NF90_DOUBLE, (/ time_id, z_size_id, y_size_id, x_size_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_3D_time))
  end subroutine define_and_write_variable_3D

  subroutine check(status)
    integer(4), intent ( in) :: status
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine

  !> Handle potential conditional diagnostics conflict
  !! Provides a more helpful error in the case where conditional diagnostics are requested as output,
  !! but their components are not enabled.
  !! We check this by searching the io_xml_configuration.
  !! @param raw_contents, intended to be the io_xml_configuration character array
  !! @param options_database
  subroutine check_for_condi_conflict(raw_contents, options_database_string)
    character, dimension(:), intent(in) :: raw_contents ! here io_xml_configuration
    !type(hashmap_type), intent(inout) :: options_database
    character(len=STRING_LENGTH), dimension(1200,2), intent(inout) :: options_database_string
    character(len=size(raw_contents)) :: string_to_process
    integer :: i, iter
    logical :: diag_enabled = .false.

    do iter = 51,350
      if (options_database_string(iter,1) .eq. "conditional_diagnostics_column_enabled") then
        diag_enabled = .true.
      end if
    end do
    if (diag_enabled .eqv. .false.) then
      do i=1, size(raw_contents)
        string_to_process(i:i)=raw_contents(i)
      end do
      !print*,"string_to_process = ",string_to_process
      if (index(string_to_process,"CondDiags_") .ne. 0) then
        call log_log(LOG_ERROR, &
            "Conditional diagnostics are DISABLED but requested via xml.  Enable or remove request to resolve.")
      end if
    end if
  end subroutine check_for_condi_conflict

  !> Awaits a command or shutdown from MONC processes and other IO servers
  !! @param command The command received is output
  !! @param source The source process received is output
  !! @returns Whether to continue polling for commands (and whether to process the current output)
  !logical function await_command(command, source, data_buffer)
  !logical function await_command(current_state, command)!, source, data_buffer)
  !  integer, intent(out) :: command!, source
  !  type(model_state_type), target, intent(inout) :: current_state
  !  !character, dimension(:), allocatable :: data_buffer

  !  logical :: completed!, inter_io_complete
  !  integer :: ierr

  !  completed=.false.
  !  await_command=.false.
  !  call mpi_recv(command, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  !  if (command .eq. 9999) then
  !    await_command = .true.
  !    call exit
  !  end if
  !  print*,"await_command = ",await_command

    !do while(.not. completed)
    !  if (.not. continue_poll_messages .and. .not. continue_poll_interio_messages) return
    !  if (continue_poll_messages) then
    !    if (test_for_command(command, source)) then
    !      await_command=.true.
    !      return
    !    end if
    !  end if
    !  if (continue_poll_interio_messages .and. allocated(io_configuration%inter_io_communications)) then
    !    inter_io_complete=test_for_inter_io(io_configuration%inter_io_communications, &
    !         io_configuration%number_inter_io_communications, io_configuration%io_communicator, command, source, data_buffer)
    !    if (inter_io_complete) then
    !      await_command=.true.
    !      return
    !    end if
    !  end if
    !  if (.not. continue_poll_messages .and. .not. already_registered_finishing_call) then
    !    if (check_diagnostic_federator_for_completion(io_configuration) .and. &
    !        (.not. any_pending()) .and. threadpool_is_idle()) then
    !      already_registered_finishing_call=.true.
    !      call perform_global_callback(io_configuration, "termination", 1, termination_callback)
    !    end if
    !  end if
    !  if (.not. completed) call pause_for_mpi_interleaving()
    !end do
  !end function await_command

  !> This is the termination callback which is called once all MONCs have deregistered, no sends are active by inter IO
  !! communications and all threads are idle. This shuts down the inter IO listening and kickstarts finalisation and closure
  !! @param io_configuration The IO server configuration
  !! @param values Values (ignored)
  !! @param field_name Field name identifier
  !! @param timestep Timestep identifier
  subroutine termination_callback(io_configuration, values, field_name, timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(DEFAULT_PRECISION), dimension(:) :: values
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep

    continue_poll_interio_messages=.false.
  end subroutine termination_callback

  !> Called to handle a specific command that has been received
  !! @param command The command which has been received from some process
  !! @param source The PID of the source (MONC) process
!   subroutine handle_command_message(command, source, data_buffer)
!     integer, intent(in) :: command, source
!     character, dimension(:), allocatable, intent(inout) :: data_buffer
!
!     if (command == REGISTER_COMMAND) then
!       if (l_thoff) then
!         call handle_monc_registration((/ source /))
!       else
!         !call threadpool_start_thread(handle_monc_registration, (/ source /)) ! original
!         call handle_monc_registration((/ source /)) ! modif
!       end if
!     else if (command == DEREGISTER_COMMAND) then
!       if (l_thoff) then
!         call handle_deregistration_command((/ source /))
!       else
!         print*,"DEREGISTER_COMMAND0"
!         call threadpool_start_thread(handle_deregistration_command, (/ source /))
!         print*,"DEREGISTER_COMMAND1"
!       end if
!     else if (command == INTER_IO_COMMUNICATION) then
!       if (l_thoff) then
!         call handle_inter_io_communication_command((/ source /), data_buffer=data_buffer)
!       else
!         call threadpool_start_thread(handle_inter_io_communication_command, (/ source /), data_buffer=data_buffer)
!       end if
!       deallocate(data_buffer)
!     else if (command .ge. DATA_COMMAND_START) then
!       print*,"handle_command_message_command000 and command =", command
!       call pull_back_data_message_and_handle(source, command-DATA_COMMAND_START)
!       print*,"handle_command_message_command111 and command =", command
!     end if
!   end subroutine handle_command_message

  !> Handles inter IO server communications
  !! @param arguments The thread based arguments, this is the index of the inter IO server description
!   subroutine handle_inter_io_communication_command(arguments, data_buffer)
!     integer, dimension(:), intent(in) :: arguments
!     character, dimension(:), allocatable, intent(inout), optional :: data_buffer
!
!     integer :: source
!
!     source=arguments(1)
!
!     call io_configuration%inter_io_communications(source)%handling_procedure(io_configuration, data_buffer, source)
!   end subroutine handle_inter_io_communication_command
!
!   !> Frees up the memory associated with individual registered MONCs. This is done at the end for all MONCs as we can't
!   !! deallocate dynamically in a threaded environment without excessive ordering and locking in case some data processing
!   !! is queued or in progress
!   subroutine free_individual_registered_monc_aspects()
!     integer :: i, specific_monc_data_type
!     type(iterator_type) :: types_iterator
!
!     do i=1, size(io_configuration%registered_moncs)
!       types_iterator=c_get_iterator(io_configuration%registered_moncs(i)%registered_monc_types)
!       do while (c_has_next(types_iterator))
!         specific_monc_data_type=c_get_integer(c_next_mapentry(types_iterator))
!         call free_mpi_type(specific_monc_data_type)
!       end do
!       if (allocated(io_configuration%registered_moncs(i)%field_start_locations)) &
!            deallocate(io_configuration%registered_moncs(i)%field_start_locations)
!       if (allocated(io_configuration%registered_moncs(i)%field_end_locations)) &
!            deallocate(io_configuration%registered_moncs(i)%field_end_locations)
!       if (allocated(io_configuration%registered_moncs(i)%definition_names)) &
!            deallocate(io_configuration%registered_moncs(i)%definition_names)
!       if (allocated(io_configuration%registered_moncs(i)%dimensions)) deallocate(io_configuration%registered_moncs(i)%dimensions)
!     end do
!   end subroutine free_individual_registered_monc_aspects
!
!   !> Deregisteres a specific MONC source process
!   !! @param source The MONC process PID that we are deregistering
!   subroutine handle_deregistration_command(arguments, data_buffer)
!     integer, dimension(:), intent(in) :: arguments
!     character, dimension(:), allocatable, intent(inout), optional :: data_buffer
!
!     integer :: monc_location, source
!
!     source=arguments(1)
!     monc_location=get_monc_location(io_configuration, source)
!     call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
!     do while (io_configuration%registered_moncs(monc_location)%active_threads .gt. 0)
!       call check_thread_status(forthread_cond_wait(io_configuration%registered_moncs(monc_location)%deactivate_condition_variable,&
!              io_configuration%registered_moncs(monc_location)%active_mutex))
!     end do
!     call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))
!     call check_thread_status(forthread_rwlock_wrlock(monc_registration_lock))
!     io_configuration%active_moncs=io_configuration%active_moncs-1
!     if (io_configuration%active_moncs==0) continue_poll_messages=.false.
!     call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))
!   end subroutine handle_deregistration_command
!
!   !> Retrieves the message from MONC off the data channel and throws this to a thread in the thread pool to actually process
!   !! We do it this way to enforce ordering between the command (including the data set ID) and the raw data itself
!   !! @param source Source PID of the MONC process
!   !! @param data_set ID of the data set being communicated
!   subroutine pull_back_data_message_and_handle(source, data_set)
!     integer, intent(in) :: source, data_set
!
!     integer :: specific_monc_data_type, specific_monc_buffer_size, recv_count, monc_location, matched_datadefn_index
!     character, dimension(:), allocatable :: data_buffer
!
!     call check_thread_status(forthread_rwlock_rdlock(monc_registration_lock))
!     monc_location=get_monc_location(io_configuration, source)
!
!     specific_monc_data_type=c_get_integer(io_configuration%registered_moncs(monc_location)%registered_monc_types, &
!          conv_to_string(data_set))
!     specific_monc_buffer_size=c_get_integer(io_configuration%registered_moncs(monc_location)%registered_monc_buffer_sizes, &
!          conv_to_string(data_set))
!
!     allocate(data_buffer(specific_monc_buffer_size))
!     recv_count=data_receive(specific_monc_data_type, 1, source, dump_data=data_buffer, data_dump_id=data_set)
!
!
!     ! This call is not handled by threading...should aid in ensuring that all time points are listed sequentially
!     matched_datadefn_index=retrieve_data_definition(io_configuration, &
!          io_configuration%registered_moncs(monc_location)%definition_names(data_set))
!     if (matched_datadefn_index .gt. 0) then
!       call inform_writer_federator_time_point(io_configuration, source, data_set, data_buffer)
!     end if
!
!     call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))
!
!     if (l_thoff) then
!       call handle_data_message((/ source,  data_set /), data_buffer=data_buffer)
!     else
!       !call threadpool_start_thread(handle_data_message, (/ source,  data_set /), data_buffer=data_buffer) ! LAMBERT
!       call handle_data_message((/ source,  data_set /), data_buffer=data_buffer)
!     end if
!
!     deallocate(data_buffer)
!   end subroutine pull_back_data_message_and_handle
!
!   !> Handles the command for data download from a specific process. This will allocate the receive buffer
!   !! and then call to get the data. Once it has been received then the data is run against handling rules
!   !! @param arguments, element 1 is the source & element 2 is the data_set
!   !! @param data_buffer The actual data from MONC read from the data channel
!   subroutine handle_data_message(arguments, data_buffer)
!     integer, dimension(:), intent(in) :: arguments
!     character, dimension(:), allocatable, intent(inout), optional :: data_buffer
!
!     integer :: monc_location, data_set, source, matched_datadefn_index
!
!     source=arguments(1)
!     data_set=arguments(2)
!     print*,"handle_data_message0 and command-4 = ",data_set
!     call check_thread_status(forthread_rwlock_rdlock(monc_registration_lock))
!     monc_location=get_monc_location(io_configuration, source)
!
!     call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
!     io_configuration%registered_moncs(monc_location)%active_threads=&
!          io_configuration%registered_moncs(monc_location)%active_threads+1
!     call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))
!     !print*,"handle_data_message_aaaa"
!     matched_datadefn_index=retrieve_data_definition(io_configuration, &
!          io_configuration%registered_moncs(monc_location)%definition_names(data_set))
!     if (matched_datadefn_index .gt. 0) then
!       !print*,"handle_data_message_ccccc"
!       call pass_fields_to_diagnostics_federator(io_configuration, source, data_set, data_buffer)
!       print*,"handle_data_message_ddddd"
!       call provide_monc_data_to_writer_federator(io_configuration, source, data_set, data_buffer) ! LAMBERT recoit et convertit les données de bytes en interger/float ---- essaie de passer les commandes writer
!       print*,"handle_data_message_aaa and command-4 = ",data_set
!       call check_writer_for_trigger(io_configuration, source, data_set, data_buffer) ! ordonne l'écriture LAMBERT
!       !print*,"handle_data_message_fffff"
!     else
!       call log_log(LOG_WARN, "IO server can not find matching data definition with name "&
!            //io_configuration%registered_moncs(monc_location)%definition_names(data_set))
!     end if
!
!     call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
!     io_configuration%registered_moncs(monc_location)%active_threads=&
!          io_configuration%registered_moncs(monc_location)%active_threads-1
!
!     call check_thread_status(forthread_cond_signal(io_configuration%registered_moncs(monc_location)%deactivate_condition_variable))
!     call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))
!     call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))
!     print*,"handle_data_message1 and command-4 = ",data_set
!   end subroutine handle_data_message
!
!   !> Handles registration from some MONC process. The source process sends some data description to this IO server which
!   !! basically tells the IO server the size of the array datas (which might be different on different processes in the case
!   !! of uneven decomposition.) Based upon this a communication (MPI) data type is constructed and the data size in bytes determined
!   !! @param source The PID of the MONC process that is registering itself
!   subroutine handle_monc_registration(arguments, data_buffer)
!     integer, dimension(:), intent(in) :: arguments
!     character, dimension(:), allocatable, intent(inout), optional :: data_buffer
!
!     integer :: configuration_send_request(3), ierr, number_data_definitions, this_monc_index, source
!
!     source=arguments(1)
!     configuration_send_request=send_configuration_to_registree(source)
!     number_data_definitions=io_configuration%number_of_data_definitions
!
!     call check_thread_status(forthread_rwlock_wrlock(monc_registration_lock))
!
!     io_configuration%number_of_moncs=io_configuration%number_of_moncs+1
!     this_monc_index=io_configuration%number_of_moncs
!     if (io_configuration%number_of_moncs .gt. size(io_configuration%registered_moncs)) then
!       call log_log(LOG_ERROR, "You have a high ratio of computational cores to IO servers, the limit is currently 100")
!       ! The extension of the MONC registration array is broken as the pointers involved in the map does not get copied across
!       ! we could manually do this, but that is for another day! If you need to extend these limits either increase the constants
!       ! or fix the extension, I don't think it will be too hard to fix the extension bit (copy the maps manually)
!       call extend_registered_moncs_array(io_configuration)
!     end if
!
!     io_configuration%active_moncs=io_configuration%active_moncs+1
!     call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))
!
!     call c_put_integer(io_configuration%monc_to_index, conv_to_string(source), this_monc_index)
!
!     call check_thread_status(forthread_mutex_init(io_configuration%registered_moncs(this_monc_index)%active_mutex, -1))
!     call check_thread_status(forthread_cond_init(&
!          io_configuration%registered_moncs(this_monc_index)%deactivate_condition_variable, -1))
!     io_configuration%registered_moncs(this_monc_index)%active_threads=0
!     io_configuration%registered_moncs(this_monc_index)%source_id=source
!
!     allocate(io_configuration%registered_moncs(this_monc_index)%field_start_locations(number_data_definitions), &
!          io_configuration%registered_moncs(this_monc_index)%field_end_locations(number_data_definitions), &
!          io_configuration%registered_moncs(this_monc_index)%definition_names(number_data_definitions), &
!          io_configuration%registered_moncs(this_monc_index)%dimensions(number_data_definitions))
!
!     ! Wait for configuration to have been sent to registree
!     call waitall_for_mpi_requests(configuration_send_request, 3)
!     call init_data_definition(source, io_configuration%registered_moncs(this_monc_index))
!   end subroutine handle_monc_registration
!
!   !> Sends the data and field descriptions to the MONC process that just registered with the IO server
!   !! @param source The MPI rank (MPI_COMM_WORLD) of the registree
!   !! @returns The nonblocking send request handles which can be waited for completion later (overlap compute and communication)
!   function send_configuration_to_registree(source)
!     integer, intent(in) :: source
!     integer :: send_configuration_to_registree(3)
!
!     integer :: ierr, srequest(3)
!
!     !**** a voir si l'asynchrone marche
!     call lock_mpi()
!     call mpi_isend(registree_definition_descriptions, size(registree_definition_descriptions), mpi_type_definition_description, &
!          source, DATA_TAG, MPI_COMM_WORLD, srequest(1), ierr) ! LAMBERT envoie dans iobridge.F90/L772
!          !print*,"size(registree_definition_descriptions) = ",size(registree_definition_descriptions)
!     call mpi_isend(registree_field_descriptions, size(registree_field_descriptions), mpi_type_field_description, &
!          source, DATA_TAG, MPI_COMM_WORLD, srequest(2), ierr) ! LAMBERT envoie dans iobridge.F90/L777
!          !print*,"size(registree_field_descriptions) = ",size(registree_field_descriptions)
!     call mpi_isend(sample_output_pairs, size(sample_output_pairs), MPI_INT, &
!          source, DATA_TAG, MPI_COMM_WORLD, srequest(3), ierr) ! LAMBERT envoie dans iobridge.F90/L790
!          !print*,"size(sample_output_pairs) = ",size(sample_output_pairs)
!     call unlock_mpi()
!     call unlock_mpi()
!
!     send_configuration_to_registree=srequest
!   end function send_configuration_to_registree
!
!   !> Initialise the sizing of data definitions from a MONC process. The IO server determines, from configuration, the
!   !! structure of each data definition but the size of the arrays depends upon the MONC process (due to uneven distribution
!   !! of data etc...) This receives the sizing message and then builds the MPI datatype for each data definition that the IO
!   !! server will receive from that specific MONC process. The field sizings are for all fields in every data definition, and
!   !! these are applied to each data definition which will simply ignore non matching fields
!   !! @param source The source MONC PID
!   !! @param monc_defn The corresponding MONC definition data structure
!   subroutine init_data_definition(source, monc_defn)
!     integer, intent(in) :: source
!     type(io_configuration_registered_monc_type), intent(inout) :: monc_defn
!
!     type(data_sizing_description_type) :: data_description(io_configuration%number_of_distinct_data_fields+4)
!     integer :: created_mpi_type, data_size, recv_count, i
!     type(data_sizing_description_type) :: field_description
!     logical :: field_found
!
!     recv_count=data_receive(mpi_type_data_sizing_description, io_configuration%number_of_distinct_data_fields+4, &
!          source, description_data=data_description)
!
!     call handle_monc_dimension_information(data_description, monc_defn)
!
!     do i=1, io_configuration%number_of_data_definitions
!       created_mpi_type=build_mpi_datatype(io_configuration%data_definitions(i), data_description, data_size, &
!            monc_defn%field_start_locations(i), monc_defn%field_end_locations(i), monc_defn%dimensions(i))
!
!       call c_put_integer(monc_defn%registered_monc_types, conv_to_string(i), created_mpi_type)
!       call c_put_integer(monc_defn%registered_monc_buffer_sizes, conv_to_string(i), data_size)
!
!       monc_defn%definition_names(i)=io_configuration%data_definitions(i)%name
!     end do
!     if (.not. initialised_present_data) then
!       initialised_present_data=.true.
!       field_found=get_data_description_from_name(data_description, NUMBER_Q_INDICIES_KEY, field_description)
!       call c_put_integer(io_configuration%dimension_sizing, "active_q_indicies", field_description%dim_sizes(1))
!       call register_present_field_names_to_federators(data_description, recv_count)
!     end if
!     call get_monc_information_data(source)
!   end subroutine init_data_definition
!
!   !> Retrieves MONC information data, this is sent by MONC (and received) regardless, but only actioned if the data has not
!   !! already been set
!   !! @param source MONC source process
!   subroutine get_monc_information_data(source)
!     integer, intent(in) :: source
!
!     character, dimension(:), allocatable :: buffer
!     character(len=STRING_LENGTH) :: q_field_name, tracer_name, cd_field_name
!     integer :: buffer_size, z_size, num_q_fields, num_tracers, n, current_point, recv_count
!     type(data_sizing_description_type) :: field_description
!     real(kind=DEFAULT_PRECISION) :: dreal
!     logical :: field_found
!
!
!     z_size=c_get_integer(io_configuration%dimension_sizing, "z")
!     num_q_fields=c_get_integer(io_configuration%dimension_sizing, "qfields")
!     num_tracers=c_get_integer(io_configuration%dimension_sizing, "tfields")
!
!     buffer_size=(kind(dreal)*z_size)*2 + (STRING_LENGTH * num_q_fields) + STRING_LENGTH * num_tracers &
!                  + 2*ncond*STRING_LENGTH + 2*ndiag*STRING_LENGTH
!     allocate(buffer(buffer_size))
!     recv_count=data_receive(MPI_BYTE, buffer_size, source, buffer)
!     if (.not. io_configuration%general_info_set) then
!       call check_thread_status(forthread_mutex_lock(io_configuration%general_info_mutex))
!       if (.not. io_configuration%general_info_set) then
!         io_configuration%general_info_set=.true.
!         allocate(io_configuration%zn_field(z_size))
!         allocate(io_configuration%z_field(z_size))
!         io_configuration%zn_field=transfer(buffer(1:kind(dreal)*z_size), io_configuration%zn_field)
!         current_point=(kind(dreal)*z_size)
!         if (num_q_fields .gt. 0) then
!           do n=1, num_q_fields
!             q_field_name=transfer(buffer(current_point+1:current_point+STRING_LENGTH), q_field_name)
!             current_point=current_point+STRING_LENGTH
!             call replace_character(q_field_name, " ", "_")
!             call c_add_string(io_configuration%q_field_names, q_field_name)
!           end do
!         end if
!
!         if (num_tracers .gt. 0) then
!           do n=1, num_tracers
!             tracer_name=transfer(buffer(current_point+1:current_point+STRING_LENGTH), tracer_name)
!             current_point=current_point+STRING_LENGTH
!             call c_add_string(io_configuration%tracer_names, tracer_name)
!           end do
!         end if
!
!         io_configuration%z_field=transfer(buffer(current_point+1:current_point+kind(dreal)*z_size), &
!                                           io_configuration%z_field)
!         current_point=current_point+(kind(dreal)*z_size)
!
!         do n=1,ncond
!           cond_request(n)=transfer(buffer(current_point+1:current_point+STRING_LENGTH), cd_field_name)
!           current_point=current_point+STRING_LENGTH
!           cond_long(n)=transfer(buffer(current_point+1:current_point+STRING_LENGTH), cd_field_name)
!           current_point=current_point+STRING_LENGTH
!         end do
!
!         do n=1,ndiag
!           diag_request(n)=transfer(buffer(current_point+1:current_point+STRING_LENGTH), cd_field_name)
!           current_point=current_point+STRING_LENGTH
!           diag_long(n)=transfer(buffer(current_point+1:current_point+STRING_LENGTH), cd_field_name)
!           current_point=current_point+STRING_LENGTH
!         end do
!
!       end if
!       call provide_q_field_names_to_writer_federator(io_configuration%q_field_names)
!       call provide_tracer_names_to_writer_federator(io_configuration%tracer_names)
!       call check_thread_status(forthread_mutex_unlock(io_configuration%general_info_mutex))
!     end if
!     deallocate(buffer)
!   end subroutine get_monc_information_data
!
!   !> Registers with the writer federator the set of fields (prognostic and diagnostic) that are available, this is based on
!   !! the array/optional fields present from MONC and the non-optional scalars. This is quite an expensive operation, so only
!   !! done once
!   !! @param data_description Array of data descriptions from MONC
!   !! @param recv_count Number of data descriptions
!   subroutine register_present_field_names_to_federators(data_description, recv_count)
!     type(data_sizing_description_type), dimension(:), intent(in) :: data_description
!     integer, intent(in) :: recv_count
!
!     type(hashset_type) :: present_field_names
!     type(hashmap_type) :: diagnostics_field_names_and_roots
!     integer :: i, j
!
!     do i=1, recv_count
!       call c_add_string(present_field_names, data_description(i)%field_name)
!     end do
!     do i=1, io_configuration%number_of_data_definitions
!       do j=1, io_configuration%data_definitions(i)%number_of_data_fields
!         if (io_configuration%data_definitions(i)%fields(j)%field_type == SCALAR_FIELD_TYPE .and. .not. &
!              io_configuration%data_definitions(i)%fields(j)%optional) then
!           call c_add_string(present_field_names, io_configuration%data_definitions(i)%fields(j)%name)
!         end if
!       end do
!     end do
!     call c_add_string(present_field_names, "time")
!     call c_add_string(present_field_names, "timestep")
!     call inform_writer_federator_fields_present(io_configuration, present_field_names)
!     diagnostics_field_names_and_roots=determine_diagnostics_fields_available(present_field_names)
!     call inform_writer_federator_fields_present(io_configuration, diag_field_names_and_roots=diagnostics_field_names_and_roots)
!     call c_free(present_field_names)
!     call c_free(diagnostics_field_names_and_roots)
!   end subroutine register_present_field_names_to_federators
!
!   !> Handles the provided local MONC dimension and data layout information
!   !! @param data_description The data descriptions sent over from MONC
!   !! @param monc_defn The corresponding MONC definition data structure
!   subroutine handle_monc_dimension_information(data_description, monc_defn)
!     type(io_configuration_registered_monc_type), intent(inout) :: monc_defn
!     type(data_sizing_description_type), dimension(:) :: data_description
!
!     type(data_sizing_description_type) :: field_description
!     integer :: i
!     logical :: field_found
!
!     field_found=get_data_description_from_name(data_description, LOCAL_SIZES_KEY, field_description)
!     if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local size information")
!     do i=1,3
!       monc_defn%local_dim_sizes(i)=field_description%dim_sizes(i)
!     end do
!     field_found=get_data_description_from_name(data_description, LOCAL_START_POINTS_KEY, field_description)
!     if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local start point information")
!     do i=1,3
!       monc_defn%local_dim_starts(i)=field_description%dim_sizes(i)
!     end do
!     field_found=get_data_description_from_name(data_description, LOCAL_END_POINTS_KEY, field_description)
!     if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local end point information")
!     do i=1,3
!       monc_defn%local_dim_ends(i)=field_description%dim_sizes(i)
!     end do
!   end subroutine handle_monc_dimension_information
end module io_server_mod
