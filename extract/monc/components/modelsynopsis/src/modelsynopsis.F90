!> Displays information about the current state_mod of the model run
module modelsynopsis_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use state_mod, only : model_state_type, TIME_TERMINATION_REASON, TIMESTEP_TERMINATION_REASON, MESSAGE_TERMINATION_REASON, &
       WALLTIME_TERMINATION_REASON
  use logging_mod, only : LOG_INFO, log_log, log_is_master, log_newline
  use conversions_mod, only : conv_to_string
  use optionsdatabase_mod, only : options_get_integer, options_get_real, options_get_string
  use mpi, only : mpi_wtime
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: previous_ts, reporting_frequency
  double precision :: start_time

  public initialisation_callback_modelsynopsis, timestep_callback_modelsynopsis, &
         finalisation_callback_modelsynopsis

contains

  subroutine initialisation_callback_modelsynopsis(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: iter, intnum

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "display_synopsis_frequency") then
        read(current_state%options_database_string(iter,2),*) intnum
        reporting_frequency = intnum
      end if
    end do
    previous_ts=current_state%timestep
    start_time=mpi_wtime()
  end subroutine initialisation_callback_modelsynopsis

  !> Timestep callback hook which performs the halo swapping for each prognostic field
  !!
  !! In parallel this is performed with MPI communication calls and wrapping around. In serial still
  !! need to wrap data around
  !! @param current_state The current model state_mod
  subroutine timestep_callback_modelsynopsis(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    double precision :: end_time

    if (mod(current_state%timestep, reporting_frequency)==0 .and. log_is_master()) then
      end_time=mpi_wtime()
      call log_newline()
      call log_log(LOG_INFO, "Number of completed timesteps "//conv_to_string(current_state%timestep))
      call log_log(LOG_INFO, "Completed "//trim(conv_to_string((current_state%timestep-previous_ts)))//&
           " timesteps in "//trim(conv_to_string(int((end_time-start_time) * 1000)))//"ms")
      call log_log(LOG_INFO, "Model time "//trim(conv_to_string(current_state%time, 5))//" seconds; dtm="//&
           trim(conv_to_string(current_state%dtm, 5)))
      previous_ts=current_state%timestep
      start_time=mpi_wtime()
    end if
  end subroutine timestep_callback_modelsynopsis

  !> Called at the end of the MONC run, will log the reason why the model is terminating
  !!@param current_state The current model state
  subroutine finalisation_callback_modelsynopsis(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: iter
    real :: realnum, termination_time_final
    character(len=STRING_LENGTH) :: walltime_limit_final

    if (log_is_master()) then
      call log_newline()
      do iter = 1,current_state%config_args
        if (current_state%options_database_string(iter,1) .eq. "termination_time") then
          read(current_state%options_database_string(iter,2),*) realnum
          termination_time_final = realnum
        end if
        if (current_state%options_database_string(iter,1) .eq. "walltime_limit") then
          walltime_limit_final = current_state%options_database_string(iter,2)
        end if
      end do
      if (current_state%termination_reason == TIME_TERMINATION_REASON) then
        call log_log(LOG_INFO, "Model run complete due to model time "//&
             trim(conv_to_string(current_state%time, 2))//" exceeding limit of "//&
             trim(conv_to_string(termination_time_final)))
             !trim(conv_to_string(options_get_real(current_state%options_database, "termination_time"), 2)))
      else if (current_state%termination_reason == TIMESTEP_TERMINATION_REASON) then
        call log_log(LOG_INFO, "Model run complete due to timestep completion, model time is "//&
             trim(conv_to_string(current_state%time, 2)))
      else if (current_state%termination_reason == MESSAGE_TERMINATION_REASON) then
        call log_log(LOG_INFO, "Model run complete due to messages file containing termination command, model time is "//&
             trim(conv_to_string(current_state%time, 2)))
      else if (current_state%termination_reason == WALLTIME_TERMINATION_REASON) then
        call log_log(LOG_INFO, "Model run complete due to walltime limit of '"//&
             trim(walltime_limit_final)//"' reached, model time is "//&
             !trim(options_get_string(current_state%options_database, "walltime_limit"))//"' reached, model time is "//&
             trim(conv_to_string(current_state%time, 2)))
      else
        call log_log(LOG_INFO, "Model run complete due to unknown reason, model time is "//&
             trim(conv_to_string(current_state%time, 2)))
      end if
    end if
  end subroutine finalisation_callback_modelsynopsis
end module modelsynopsis_mod
