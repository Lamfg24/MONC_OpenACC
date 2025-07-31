










!> This component will check for termination conditions at stages of the model run and terminate that
!! specific stage if the parameters have been met
module terminationcheck_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use state_mod, only : model_state_type, TIME_TERMINATION_REASON, TIMESTEP_TERMINATION_REASON, MESSAGE_TERMINATION_REASON, &
       WALLTIME_TERMINATION_REASON
  use conversions_mod, only : conv_single_real_to_double, conv_to_integer, conv_to_string,&
                              conv_to_lowercase
  use optionsdatabase_mod, only : options_get_integer, options_has_key, options_get_real, options_add, options_get_string
  use logging_mod, only : LOG_WARN, log_master_log
  use mpi, only : MPI_INT, MPI_LOGICAL, MPI_IN_PLACE, MPI_LOR, mpi_wtime
  implicit none

  private

  integer, parameter :: FILE_LINE_LEN=100, FILE_UNIT=10
  integer :: max_timesteps, check_messages_file_frequency, check_walltime_frequency, max_walltime_secs
  !real(kind=DEFAULT_PRECISION) :: termination_time
  character(len=STRING_LENGTH) :: messages_file_name
  logical :: check_for_walltime

  public terminationcheck_get_descriptor, init_callback_terminationcheck, timestep_callback_terminationcheck

  contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type_v1) function terminationcheck_get_descriptor()
    terminationcheck_get_descriptor%name="termination_check"
    terminationcheck_get_descriptor%version=0.1
    terminationcheck_get_descriptor%published_fields_on_off = .false.
    !terminationcheck_get_descriptor%initialisation=>init_callback
    !terminationcheck_get_descriptor%timestep=>timestep_callback
  end function terminationcheck_get_descriptor

  !> Called upon model initialisation. Will basically read from the options database and set options in
  !! the database that are appropriate
  !! @param current_state The current model state_mod
  subroutine init_callback_terminationcheck(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: i, idx, pidx, walltime_secs, walltime_mins, walltime_hours, intnum, iter
    character(len=STRING_LENGTH) :: walltime_string
    logical :: mangled
    real :: realnum

    !max_timesteps=options_get_integer(current_state%options_database, "nn_timesteps")
    !termination_time=options_get_real(current_state%options_database, "termination_time")
    !check_messages_file_frequency=options_get_integer(current_state%options_database, "check_msg_frequency")
    !messages_file_name=options_get_string(current_state%options_database, "msg_filename")
    !check_walltime_frequency=options_get_integer(current_state%options_database, "check_walltime_frequency")
    !walltime_string=options_get_string(current_state%options_database, "walltime_limit")
     do iter = 1,current_state%config_args
       if (current_state%options_database_string(iter,1) .eq. "nn_timesteps") then
         read(current_state%options_database_string(iter,2),*) intnum
         max_timesteps = intnum
       else if (current_state%options_database_string(iter,1) .eq. "termination_time") then
         read(current_state%options_database_string(iter,2),*) realnum
         current_state%termination_time = realnum
       else if (current_state%options_database_string(iter,1) .eq. "check_msg_frequency") then
         read(current_state%options_database_string(iter,2),*) intnum
         check_messages_file_frequency = intnum
       else if (current_state%options_database_string(iter,1) .eq. "msg_filename") then
         messages_file_name = current_state%options_database_string(iter,2)
       else if (current_state%options_database_string(iter,1) .eq. "check_walltime_frequency") then
         read(current_state%options_database_string(iter,2),*) intnum
         check_walltime_frequency = intnum
       else if (current_state%options_database_string(iter,1) .eq. "walltime_limit") then
         walltime_string = current_state%options_database_string(iter,2)
       end if
     end do
     check_for_walltime=conv_to_lowercase(trim(walltime_string)) /= "none"
     if (check_for_walltime) then
       pidx=1
       mangled=.false.
       do i=1, 2
         idx=index(walltime_string(pidx:), ":")
         if (idx .gt. 0) then
           read(walltime_string(pidx:pidx+idx-2),*) intnum
           if (i==1) walltime_hours=intnum
           if (i==2) walltime_mins=intnum
           pidx=pidx+idx
         else
           call log_master_log(LOG_WARN, "Walltime limit of  `"//trim(walltime_string)//&
                "` does not contains hh:mm:ss, defaulting to no limit")
           check_for_walltime=.false.
           exit
         end if
       end do
       if (check_for_walltime) then
          read(walltime_string(pidx:),*) intnum
          walltime_secs=intnum
         if (walltime_mins .lt. 0 .or. walltime_mins .gt. 59) then
           walltime_mins=0
           mangled=.true.
         end if
         if (walltime_secs .lt. 0 .or. walltime_secs .gt. 59) then
           walltime_secs=0
           mangled=.true.
         end if
         if (mangled) then
           call log_master_log(LOG_WARN, "Walltime limit of `"//trim(walltime_string)//"` mangled, defaulting to "//&
                trim(conv_to_string(walltime_hours))//":"//trim(conv_to_string(walltime_mins))//":"//&
                trim(conv_to_string(walltime_secs)))
         end if
         max_walltime_secs=(walltime_hours*60*60)+(walltime_mins*60)+walltime_secs
       end if
     end if
  end subroutine init_callback_terminationcheck

  !> Timestep hook which is called at each timestep to determine whether or not to terminate timestep iterations
  !! @param current_state The current model state_mod
  subroutine timestep_callback_terminationcheck(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: ierr, file_message_status

    if (max_timesteps .gt. 0) then
      current_state%continue_timestep=mod(current_state%timestep, max_timesteps) /= 0
      if (.not. current_state%continue_timestep) current_state%termination_reason=TIMESTEP_TERMINATION_REASON
    else
      current_state%continue_timestep=.true.
    end if
    if (current_state%continue_timestep) then
      current_state%continue_timestep = current_state%time .lt. current_state%termination_time
      if (.not. current_state%continue_timestep) current_state%termination_reason=TIME_TERMINATION_REASON
    end if
    !print*,"current_state%time = ",current_state%time
    !print*,"termination_time = ",termination_time
    !print*,"termination_check current_state%continue_timestepaaa = ",current_state%continue_timestep
    !print*,"termination_check check_for_walltime = ",check_for_walltime
    !print*,"termination_check mod(current_state%timestep, check_walltime_frequency) = ", &
    !        mod(current_state%timestep, check_walltime_frequency)
    !print*,"termination_check int(mpi_wtime()) = ",int(mpi_wtime())
    !print*,"termination_check int(current_state%model_start_wtime) = ",int(current_state%model_start_wtime)
    !print*,"termination_check max_walltime_secs = ",max_walltime_secs

    if (current_state%time .eq. current_state%termination_time) then
      print*,"TERMINER MONC"
      current_state%continue_timestep=.false.
    else
      continue
    end if
    !if (current_state%continue_timestep .and. check_for_walltime .and. & ! LAMBERT ESSAI
    !     mod(current_state%timestep, check_walltime_frequency) == 0) then
    !  current_state%continue_timestep=int(mpi_wtime() - current_state%model_start_wtime) .lt. max_walltime_secs ! LAMBERT
    !  print*,"termination_check current_state%continue_timestepaaawwwwwww = ",current_state%continue_timestep
    !  call mpi_allreduce(MPI_IN_PLACE, current_state%continue_timestep, 1, MPI_LOGICAL, MPI_LOR, &
    !       current_state%parallel%monc_communicator, ierr)
    !  if (.not. current_state%continue_timestep) current_state%termination_reason=WALLTIME_TERMINATION_REASON
    !end if
    !print*,"termination_check current_state%continue_timestebbb = ",current_state%continue_timestep
    !if (current_state%continue_timestep .and. mod(current_state%timestep, check_messages_file_frequency) == 0) then
    !  if (current_state%parallel%my_rank == 0) then
    !    file_message_status=check_messages_file(current_state)
    !    call mpi_bcast(file_message_status, 1, MPI_INT, 0, current_state%parallel%monc_communicator, ierr)
    !  else
    !    call mpi_bcast(file_message_status, 1, MPI_INT, 0, current_state%parallel%monc_communicator, ierr)
    !    if (file_message_status == 1) current_state%continue_timestep=.false.
    !  end if
    !  if (.not. current_state%continue_timestep) current_state%termination_reason=MESSAGE_TERMINATION_REASON
    !end if ! LAMBERT ESSAI
    !print*,"termination_check current_state%continue_timestepccc = ",current_state%continue_timestep
    !print*,"termination_check111"
  end subroutine timestep_callback_terminationcheck

  !> Checks the messages file for commands which determine user control of the model
  !! @param current_state The current model state
  integer function check_messages_file(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    character(len=FILE_LINE_LEN) :: msg_line
    integer :: ierr

    check_messages_file=0
    open (unit=FILE_UNIT, file=messages_file_name, status='OLD', iostat=ierr)
    if (ierr == 0) then
      read(FILE_UNIT,"(A)",iostat=ierr) msg_line
      if (ierr == 0) then
        if (trim(msg_line) .eq. "terminate") then
          check_messages_file=1
          current_state%continue_timestep=.false.
        end if
      end if      
    end if    
    close(FILE_UNIT)
  end function check_messages_file  
end module terminationcheck_mod
