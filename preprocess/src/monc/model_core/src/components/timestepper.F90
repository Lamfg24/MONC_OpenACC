










!> Performs the actual time stepping over groups of components. Each group can be the whole (which is one call per
!! component per timestep) or column, which calls components for each column of the timestep. Groups are executed
!! sequentially in the order that they have been configured (which is already set up in the registry)
module timestepper_mod
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  use optionsdatabase_mod, only : options_get_integer
  use clearsourceterms_mod ! also contains the no_component function
  use steppingdirection_mod
  use haloswapper_mod
  use setfluxlook_mod
  use lowerbc_mod
  use smagorinsky_mod
  use pwadvection_mod
  use tvdadvection_mod
  use thadvection_mod
  use diffusion_mod
  use viscosity_mod
  use coriolis_mod
  use buoyancy_mod
  use damping_mod
  use forcing_mod
  use set_consistent_lowbc_mod
  use simplecloud_mod
  use casim_mod
  use subgrid_profile_diagnostics_mod
  use diverr_mod
  use pressuresource_mod
  use diagnostics_3d_mod
  use profile_diagnostics_mod
  use casim_profile_dgs_mod
  use scalar_diagnostics_mod
  use stepfields_mod
  use pdf_analysis_mod
  use meanprofiles_mod
  use fftsolver_mod
  use iterativesolver_mod
  use cfltest_mod
  use pstep_mod
  use swapsmooth_mod
  use conditional_diagnostics_column_mod
  use conditional_diagnostics_whole_mod
  use checkpointer_mod
  use modelsynopsis_mod
  use terminationcheck_mod
  use tank_experiments_mod
  use iobridge_mod
  use socrates_couple_mod
  implicit none

  private

  integer :: radiation_interval
  logical :: socrates_enabled

  abstract interface
    subroutine component_func(current_state)
        use state_mod, only : model_state_type
        type(model_state_type), target, intent(inout) :: current_state
    end subroutine component_func
  end interface

  procedure (component_func), pointer :: timestep_callback_clearsourceterms_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_steppingdirection_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_haloswapper_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_setfluxlook_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_lowerbc_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_smagorinsky_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_pwadvection_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_tvdadvection_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_thadvection_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_diffusion_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_viscosity_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_coriolis_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_buoyancy_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_damping_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_forcing_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_set_consistent_lowbc_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_simplecloud_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_casim_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_subgrid_profile_diagnostics_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_diverr_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_pressuresource_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_profile_diagnostics_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_casim_profile_dgs_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_scalar_diagnostics_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_stepfields_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_pdf_analysis_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_meanprofiles_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_fftsolver_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_iterativesolver_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_cfltest_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_pstep_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_swapsmooth_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_diagnostics_column_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_diagnostics_3d_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_conditional_diagnostics_whole_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_tank_experiments_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_modelsynopsis_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_terminationcheck_ptr => no_component
  procedure (component_func), pointer :: timestep_callbtank_experimentsack_tank_experiments_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_iobridge_ptr => no_component
  procedure (component_func), pointer :: timestep_callback_socrates_ptr => no_component

  public init_timestepper, timestep, finalise_timestepper
contains

  !> Initialises the timestepper by prefetching the groups in the order that they will be executed, this is for optimised
  !! execution in the timestep calls
  !! @param current_state The current model state
  subroutine init_timestepper(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: iter, intnum
    logical :: logicnum

    !call get_ordered_groups(group_descriptors)
    !radiation_interval=options_get_integer(current_state%options_database, "rad_interval")
    !socrates_enabled=is_component_enabled(current_state%options_database, "socrates_couple")
    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "rad_interval") then
        read(current_state%options_database_string(iter,2),*) intnum
        radiation_interval = intnum
      else if (current_state%options_database_string(iter,1) .eq. "tracking_variables_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%tracking_variables_enabled = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "socrates_couple_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%socrates_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_socrates_ptr => timestep_callback_socrates
      else if (current_state%options_database_string(iter,1) .eq. "clearsourceterms_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_clearsourceterms_ptr => timestep_callback_clearsourceterms
      else if (current_state%options_database_string(iter,1) .eq. "stepping_direction_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_steppingdirection_ptr => timestep_callback_steppingdirection
      else if (current_state%options_database_string(iter,1) .eq. "halo_swapper_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_haloswapper_ptr => timestep_callback_haloswapper
      else if (current_state%options_database_string(iter,1) .eq. "setfluxlook_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_setfluxlook_ptr => timestep_callback_setfluxlook
      else if (current_state%options_database_string(iter,1) .eq. "lower_bc_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_lowerbc_ptr => timestep_callback_lowerbc
      else if (current_state%options_database_string(iter,1) .eq. "smagorinsky_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_smagorinsky_ptr => timestep_callback_smagorinsky
      else if (current_state%options_database_string(iter,1) .eq. "pw_advection_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%pw_advection_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_pwadvection_ptr => timestep_callback_pwadvection
      else if (current_state%options_database_string(iter,1) .eq. "tvd_advection_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%tvd_advection_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_tvdadvection_ptr => timestep_callback_tvdadvection
      else if (current_state%options_database_string(iter,1) .eq. "th_advection_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%th_advection_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_thadvection_ptr => timestep_callback_thadvection
      else if (current_state%options_database_string(iter,1) .eq. "diffusion_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%diffusion_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_diffusion_ptr => timestep_callback_diffusion
      else if (current_state%options_database_string(iter,1) .eq. "viscosity_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%viscosity_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_viscosity_ptr => timestep_callback_viscosity
      else if (current_state%options_database_string(iter,1) .eq. "coriolis_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%coriolis_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_coriolis_ptr => timestep_callback_coriolis
      else if (current_state%options_database_string(iter,1) .eq. "buoyancy_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%buoyancy_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_buoyancy_ptr => timestep_callback_buoyancy
      else if (current_state%options_database_string(iter,1) .eq. "damping_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_damping_ptr => timestep_callback_damping
      else if (current_state%options_database_string(iter,1) .eq. "forcing_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%forcing_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_forcing_ptr => timestep_callback_forcing
      else if (current_state%options_database_string(iter,1) .eq. "set_consistent_lowbc_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_set_consistent_lowbc_ptr => timestep_callback_set_consistent_lowbc
      else if (current_state%options_database_string(iter,1) .eq. "simplecloud_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%simplecloud_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_simplecloud_ptr => timestep_callback_simplecloud
      else if (current_state%options_database_string(iter,1) .eq. "casim_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%casim_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_casim_ptr => timestep_callback_casim
      else if (current_state%options_database_string(iter,1) .eq. "subgrid_profile_diagnostics_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%subgrid_profile_diagnostics_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_subgrid_profile_diagnostics_ptr => &
                                             timestep_callback_subgrid_profile_diagnostics
      else if (current_state%options_database_string(iter,1) .eq. "diverr_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_diverr_ptr => timestep_callback_diverr
      else if (current_state%options_database_string(iter,1) .eq. "diverr_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_diverr_ptr => timestep_callback_diverr
      else if (current_state%options_database_string(iter,1) .eq. "pressure_source_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_pressuresource_ptr => timestep_callback_pressuresource
      else if (current_state%options_database_string(iter,1) .eq. "profile_diagnostics_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%profile_diagnostics_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_profile_diagnostics_ptr => timestep_callback_profile_diagnostics
      else if (current_state%options_database_string(iter,1) .eq. "casim_profile_dgs_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%casim_profile_dgs_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_casim_profile_dgs_ptr => timestep_callback_casim_profile_dgs
      else if (current_state%options_database_string(iter,1) .eq. "scalar_diagnostics_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%scalar_diagnostics_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_scalar_diagnostics_ptr => timestep_callback_scalar_diagnostics
      else if (current_state%options_database_string(iter,1) .eq. "stepfields_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%stepfields_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_stepfields_ptr => timestep_callback_stepfields
      else if (current_state%options_database_string(iter,1) .eq. "pdf_analysis_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_pdf_analysis_ptr => timestep_callback_pdf_analysis
      else if (current_state%options_database_string(iter,1) .eq. "mean_profiles_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%mean_profiles_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_meanprofiles_ptr => timestep_callback_meanprofiles
      else if (current_state%options_database_string(iter,1) .eq. "fftsolver_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_fftsolver_ptr => timestep_callback_fftsolver
      else if (current_state%options_database_string(iter,1) .eq. "iterativesolver_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_iterativesolver_ptr => timestep_callback_iterativesolver
      else if (current_state%options_database_string(iter,1) .eq. "cfltest_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_cfltest_ptr => timestep_callback_cfltest
      else if (current_state%options_database_string(iter,1) .eq. "pstep_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%pstep_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_pstep_ptr => timestep_callback_pstep
      else if (current_state%options_database_string(iter,1) .eq. "swap_smooth_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_swapsmooth_ptr => timestep_callback_swapsmooth
      else if (current_state%options_database_string(iter,1) .eq. "conditional_diagnostics_column_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%conditional_diagnostics_column_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_diagnostics_column_ptr => timestep_callback_diagnostics_column
      else if (current_state%options_database_string(iter,1) .eq. "diagnostics_3d_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%diagnostics_3d_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_diagnostics_3d_ptr => timestep_callback_diagnostics_3d
      else if (current_state%options_database_string(iter,1) .eq. "conditional_diagnostics_whole_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%conditional_diagnostics_whole_enabled = logicnum
        if (logicnum .eqv. .true.) timestep_callback_conditional_diagnostics_whole_ptr => &
                                            timestep_callback_conditional_diagnostics_whole
      else if (current_state%options_database_string(iter,1) .eq. "tank_experiments_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_tank_experiments_ptr => timestep_callback_tank_experiments
      else if (current_state%options_database_string(iter,1) .eq. "model_synopsis_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_modelsynopsis_ptr => timestep_callback_modelsynopsis
      else if (current_state%options_database_string(iter,1) .eq. "termination_check_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_terminationcheck_ptr => timestep_callback_terminationcheck
      else if (current_state%options_database_string(iter,1) .eq. "iobridge_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) timestep_callback_iobridge_ptr => timestep_callback_iobridge
      end if
    end do
  end subroutine init_timestepper

  !> Performs a timestep, which is comprised of executing each group of components in the order that they have been configured
  !! in. The components in a group can be called, depending on the type, just once per timestep (WHOLE) or per column (COLUMN).
  !! @param current_state The current model state
  subroutine timestep(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i
    real(kind=DEFAULT_PRECISION) :: mod_freq

    if (current_state%diagnostic_file_0d_write_frequency .eq. 0) then
      current_state%modulo_number_0d = -1.0
    else
      if (current_state%time_frequency_enabled .eqv. .true.) then
        mod_freq = modulo(current_state%time, float(current_state%diagnostic_file_0d_write_frequency))
        if ((mod_freq - current_state%dtm) .le. 0.0) then
          current_state%modulo_number_0d = int(mod_freq)
        else
          current_state%modulo_number_0d = -1.0
        end if
      else
        current_state%modulo_number_0d = modulo(current_state%timestep, current_state%diagnostic_file_0d_write_frequency)
      end if
    end if
    if (current_state%diagnostic_file_1d_write_frequency .eq. 0) then
      current_state%modulo_number_1d = -1.0
    else
      if (current_state%time_frequency_enabled .eqv. .true.) then
        mod_freq = modulo(current_state%time, float(current_state%diagnostic_file_1d_write_frequency))
        if ((mod_freq - current_state%dtm) .le. 0.0) then
          current_state%modulo_number_1d = int(mod_freq)
        else
          current_state%modulo_number_1d = -1.0
        end if
      else
        current_state%modulo_number_1d = modulo(current_state%timestep, current_state%diagnostic_file_1d_write_frequency)
      end if
    end if
    if (current_state%diagnostic_file_2d_write_frequency .eq. 0) then
      current_state%modulo_number_2d = -1.0
    else
      if (current_state%time_frequency_enabled .eqv. .true.) then
        mod_freq = modulo(current_state%time, float(current_state%diagnostic_file_2d_write_frequency))
        if ((mod_freq - current_state%dtm) .le. 0.0) then
          current_state%modulo_number_2d = int(mod_freq)
        else
          current_state%modulo_number_2d = -1.0
        end if
      else
        current_state%modulo_number_2d = modulo(current_state%timestep, current_state%diagnostic_file_2d_write_frequency)
      end if
    end if
    if (current_state%diagnostic_file_3d_write_frequency .eq. 0) then
      current_state%modulo_number_3d = -1.0
    else
      if (current_state%time_frequency_enabled .eqv. .true.) then
        mod_freq = modulo(current_state%time, float(current_state%diagnostic_file_3d_write_frequency))
        if ((mod_freq - current_state%dtm) .le. 0.0) then
          current_state%modulo_number_3d = int(mod_freq)
        else
          current_state%modulo_number_3d = -1.0
        end if
      else
        current_state%modulo_number_3d = modulo(current_state%timestep, current_state%diagnostic_file_3d_write_frequency)
      end if
    end if
    if (current_state%checkpoint_frequency .eq. 0) then
      current_state%modulo_number_check = -1.0
    else
      if (current_state%time_frequency_enabled .eqv. .true.) then
        mod_freq = modulo(current_state%time, float(current_state%checkpoint_frequency))
        if ((mod_freq - current_state%dtm) .le. 0.0) then
          current_state%modulo_number_check = int(mod_freq)
        else
          current_state%modulo_number_check = -1.0
        end if
      else
        current_state%modulo_number_check = modulo(current_state%timestep, current_state%checkpoint_frequency)
      end if
    end if


    call timestep_whole_start_group_contents(current_state)
    if (current_state%tracking_variables_enabled .eqv. .true.) then
      print*,"mean(p)_ts1 = ",sum(current_state%p%data)/size(current_state%p%data)
      print*,"       "
      print*,"mean(su)_ts1 = ",sum(current_state%su%data)/size(current_state%su%data)
      print*,"mean(u)_ts1 = ",sum(current_state%u%data)/size(current_state%u%data)
      print*,"mean(zu)_ts1 = ",sum(current_state%zu%data)/size(current_state%zu%data)
      print*,"       "
      print*,"mean(sv)_ts1 = ",sum(current_state%sv%data)/size(current_state%sv%data)
      print*,"mean(v)_ts1 = ",sum(current_state%v%data)/size(current_state%v%data)
      print*,"mean(zv)_ts1 = ",sum(current_state%zv%data)/size(current_state%zv%data)
      print*,"       "
      print*,"mean(sw)_ts1 = ",sum(current_state%sw%data)/size(current_state%sw%data)
      print*,"mean(w)_ts1 = ",sum(current_state%w%data)/size(current_state%w%data)
      print*,"mean(zw)_ts1 = ",sum(current_state%zw%data)/size(current_state%zw%data)
      print*,"       "
      print*,"mean(sth)_ts1 = ",sum(current_state%sth%data)/size(current_state%sth%data)
      print*,"mean(th)_ts1 = ",sum(current_state%th%data)/size(current_state%th%data)
      print*,"mean(zth)_ts1 = ",sum(current_state%zth%data)/size(current_state%zth%data)
      print*,"       "
      print*,"mean(sqv)_ts1 = ",sum(current_state%sqv%data)/size(current_state%sqv%data)
      print*,"mean(qv)_ts1 = ",sum(current_state%qv%data)/size(current_state%qv%data)
      print*,"mean(zqv)_ts1 = ",sum(current_state%zqv%data)/size(current_state%zqv%data)
      print*,"       "
      print*,"############################################################################"
      print*,"       "
    end if

    call timestep_column_subgrid_group_contents(current_state)
    if (current_state%tracking_variables_enabled .eqv. .true.) then
      print*,"mean(p)_ts2 = ",sum(current_state%p%data)/size(current_state%p%data)
      print*,"       "
      print*,"mean(su)_ts2 = ",sum(current_state%su%data)/size(current_state%su%data)
      print*,"mean(u)_ts2 = ",sum(current_state%u%data)/size(current_state%u%data)
      print*,"mean(zu)_ts2 = ",sum(current_state%zu%data)/size(current_state%zu%data)
      print*,"       "
      print*,"mean(sv)_ts2 = ",sum(current_state%sv%data)/size(current_state%sv%data)
      print*,"mean(v)_ts2 = ",sum(current_state%v%data)/size(current_state%v%data)
      print*,"mean(zv)_ts2 = ",sum(current_state%zv%data)/size(current_state%zv%data)
      print*,"       "
      print*,"mean(sw)_ts2 = ",sum(current_state%sw%data)/size(current_state%sw%data)
      print*,"mean(w)_ts2 = ",sum(current_state%w%data)/size(current_state%w%data)
      print*,"mean(zw)_ts2 = ",sum(current_state%zw%data)/size(current_state%zw%data)
      print*,"       "
      print*,"mean(sth)_ts2 = ",sum(current_state%sth%data)/size(current_state%sth%data)
      print*,"mean(th)_ts2 = ",sum(current_state%th%data)/size(current_state%th%data)
      print*,"mean(zth)_ts2 = ",sum(current_state%zth%data)/size(current_state%zth%data)
      print*,"       "
      print*,"mean(sqv)_ts2 = ",sum(current_state%sqv%data)/size(current_state%sqv%data)
      print*,"mean(qv)_ts2 = ",sum(current_state%qv%data)/size(current_state%qv%data)
      print*,"mean(zqv)_ts2 = ",sum(current_state%zqv%data)/size(current_state%zqv%data)
      print*,"       "
      print*,"############################################################################"
      print*,"       "
    end if

    call timestep_column_dynamics_group_contents(current_state)
    if (current_state%tracking_variables_enabled .eqv. .true.) then
      print*,"mean(p)_ts3 = ",sum(current_state%p%data)/size(current_state%p%data)
      print*,"       "
      print*,"mean(su)_ts3 = ",sum(current_state%su%data)/size(current_state%su%data)
      print*,"mean(u)_ts3 = ",sum(current_state%u%data)/size(current_state%u%data)
      print*,"mean(zu)_ts3 = ",sum(current_state%zu%data)/size(current_state%zu%data)
      print*,"       "
      print*,"mean(sv)_ts3 = ",sum(current_state%sv%data)/size(current_state%sv%data)
      print*,"mean(v)_ts3 = ",sum(current_state%v%data)/size(current_state%v%data)
      print*,"mean(zv)_ts3 = ",sum(current_state%zv%data)/size(current_state%zv%data)
      print*,"       "
      print*,"mean(sw)_ts3 = ",sum(current_state%sw%data)/size(current_state%sw%data)
      print*,"mean(w)_ts3 = ",sum(current_state%w%data)/size(current_state%w%data)
      print*,"mean(zw)_ts3 = ",sum(current_state%zw%data)/size(current_state%zw%data)
      print*,"       "
      print*,"mean(sth)_ts3 = ",sum(current_state%sth%data)/size(current_state%sth%data)
      print*,"mean(th)_ts3 = ",sum(current_state%th%data)/size(current_state%th%data)
      print*,"mean(zth)_ts3 = ",sum(current_state%zth%data)/size(current_state%zth%data)
      print*,"       "
      print*,"mean(sqv)_ts3 = ",sum(current_state%sqv%data)/size(current_state%sqv%data)
      print*,"mean(qv)_ts3 = ",sum(current_state%qv%data)/size(current_state%qv%data)
      print*,"mean(zqv)_ts3 = ",sum(current_state%zqv%data)/size(current_state%zqv%data)
      print*,"       "
      print*,"############################################################################"
      print*,"       "
    end if

    call timestep_whole_solver_group_contents(current_state)
    if (current_state%tracking_variables_enabled .eqv. .true.) then
      print*,"mean(p)_ts4 = ",sum(current_state%p%data)/size(current_state%p%data)
      print*,"       "
      print*,"mean(su)_ts4 = ",sum(current_state%su%data)/size(current_state%su%data)
      print*,"mean(u)_ts4 = ",sum(current_state%u%data)/size(current_state%u%data)
      print*,"mean(zu)_ts4 = ",sum(current_state%zu%data)/size(current_state%zu%data)
      print*,"       "
      print*,"mean(sv)_ts4 = ",sum(current_state%sv%data)/size(current_state%sv%data)
      print*,"mean(v)_ts4 = ",sum(current_state%v%data)/size(current_state%v%data)
      print*,"mean(zv)_ts4 = ",sum(current_state%zv%data)/size(current_state%zv%data)
      print*,"       "
      print*,"mean(sw)_ts4 = ",sum(current_state%sw%data)/size(current_state%sw%data)
      print*,"mean(w)_ts4 = ",sum(current_state%w%data)/size(current_state%w%data)
      print*,"mean(zw)_ts4 = ",sum(current_state%zw%data)/size(current_state%zw%data)
      print*,"       "
      print*,"mean(sth)_ts4 = ",sum(current_state%sth%data)/size(current_state%sth%data)
      print*,"mean(th)_ts4 = ",sum(current_state%th%data)/size(current_state%th%data)
      print*,"mean(zth)_ts4 = ",sum(current_state%zth%data)/size(current_state%zth%data)
      print*,"       "
      print*,"mean(sqv)_ts4 = ",sum(current_state%sqv%data)/size(current_state%sqv%data)
      print*,"mean(qv)_ts4 = ",sum(current_state%qv%data)/size(current_state%qv%data)
      print*,"mean(zqv)_ts4 = ",sum(current_state%zqv%data)/size(current_state%zqv%data)
      print*,"       "
      print*,"############################################################################"
      print*,"       "
    end if

    call timestep_column_pressure_terms_group_contents(current_state)
    if (current_state%tracking_variables_enabled .eqv. .true.) then
      print*,"mean(p)_ts5 = ",sum(current_state%p%data)/size(current_state%p%data)
      print*,"       "
      print*,"mean(su)_ts5 = ",sum(current_state%su%data)/size(current_state%su%data)
      print*,"mean(u)_ts5 = ",sum(current_state%u%data)/size(current_state%u%data)
      print*,"mean(zu)_ts5 = ",sum(current_state%zu%data)/size(current_state%zu%data)
      print*,"       "
      print*,"mean(sv)_ts5 = ",sum(current_state%sv%data)/size(current_state%sv%data)
      print*,"mean(v)_ts5 = ",sum(current_state%v%data)/size(current_state%v%data)
      print*,"mean(zv)_ts5 = ",sum(current_state%zv%data)/size(current_state%zv%data)
      print*,"       "
      print*,"mean(sw)_ts5 = ",sum(current_state%sw%data)/size(current_state%sw%data)
      print*,"mean(w)_ts5 = ",sum(current_state%w%data)/size(current_state%w%data)
      print*,"mean(zw)_ts5 = ",sum(current_state%zw%data)/size(current_state%zw%data)
      print*,"       "
      print*,"mean(sth)_ts5 = ",sum(current_state%sth%data)/size(current_state%sth%data)
      print*,"mean(th)_ts5 = ",sum(current_state%th%data)/size(current_state%th%data)
      print*,"mean(zth)_ts5 = ",sum(current_state%zth%data)/size(current_state%zth%data)
      print*,"       "
      print*,"mean(sqv)_ts5 = ",sum(current_state%sqv%data)/size(current_state%sqv%data)
      print*,"mean(qv)_ts5 = ",sum(current_state%qv%data)/size(current_state%qv%data)
      print*,"mean(zqv)_ts5 = ",sum(current_state%zqv%data)/size(current_state%zqv%data)
      print*,"       "
      print*,"############################################################################"
      print*,"       "
    end if

    call timestep_whole_last_group_contents(current_state)
    if (current_state%tracking_variables_enabled .eqv. .true.) then
      print*,"mean(p)_ts6 = ",sum(current_state%p%data)/size(current_state%p%data)
      print*,"       "
      print*,"mean(su)_ts6 = ",sum(current_state%su%data)/size(current_state%su%data)
      print*,"mean(u)_ts6 = ",sum(current_state%u%data)/size(current_state%u%data)
      print*,"mean(zu)_ts6 = ",sum(current_state%zu%data)/size(current_state%zu%data)
      print*,"       "
      print*,"mean(sv)_ts6 = ",sum(current_state%sv%data)/size(current_state%sv%data)
      print*,"mean(v)_ts6 = ",sum(current_state%v%data)/size(current_state%v%data)
      print*,"mean(zv)_ts6 = ",sum(current_state%zv%data)/size(current_state%zv%data)
      print*,"       "
      print*,"mean(sw)_ts6 = ",sum(current_state%sw%data)/size(current_state%sw%data)
      print*,"mean(w)_ts6 = ",sum(current_state%w%data)/size(current_state%w%data)
      print*,"mean(zw)_ts6 = ",sum(current_state%zw%data)/size(current_state%zw%data)
      print*,"       "
      print*,"mean(sth)_ts6 = ",sum(current_state%sth%data)/size(current_state%sth%data)
      print*,"mean(th)_ts6 = ",sum(current_state%th%data)/size(current_state%th%data)
      print*,"mean(zth)_ts6 = ",sum(current_state%zth%data)/size(current_state%zth%data)
      print*,"       "
      print*,"mean(sqv)_ts6 = ",sum(current_state%sqv%data)/size(current_state%sqv%data)
      print*,"mean(qv)_ts6 = ",sum(current_state%qv%data)/size(current_state%qv%data)
      print*,"mean(zqv)_ts6 = ",sum(current_state%zqv%data)/size(current_state%zqv%data)
      print*,"       "
      print*,"############################################################################"
      print*,"############################################################################"
      print*,"       "
    end if

  end subroutine timestep

  !> Finalises the timestepper by cleaning up allocated memory
  subroutine finalise_timestepper()
  end subroutine finalise_timestepper

  !> Executes a timestep for components in a group which are designed to be executed once per timestep
  !! @param current_state The current model state
  !! @param group_descriptor Description of the group of components to execute
  subroutine timestep_whole_start_group_contents(current_state)!, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    !type(group_descriptor_type), intent(in) :: group_descriptor

    ! For print_debug_data, the column_global fields must match the requested coordinate.
    ! This is already handled for the timestep_column, but needs to be specially set for timestep_whole.
    !  This shouldù not affect update_state_sitation_flags, as that is only used for timestep_column.
    if (current_state%print_debug_data) then
      current_state%column_global_x = current_state%pdd_x
      current_state%column_global_y = current_state%pdd_y
    end if
    call timestep_callback_clearsourceterms_ptr(current_state)
    call timestep_callback_steppingdirection_ptr(current_state)
    call timestep_callback_haloswapper_ptr(current_state)
    call timestep_callback_setfluxlook_ptr(current_state)

    !call execute_timestep_callbacks(current_state, group_descriptor%id)

  end subroutine timestep_whole_start_group_contents

  subroutine timestep_whole_solver_group_contents(current_state)!, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    !type(group_descriptor_type), intent(in) :: group_descriptor

    ! For print_debug_data, the column_global fields must match the requested coordinate.
    ! This is already handled for the timestep_column, but needs to be specially set for timestep_whole.
    !  This should not affect update_state_sitation_flags, as that is only used for timestep_column.
    if (current_state%print_debug_data) then
      current_state%column_global_x = current_state%pdd_x
      current_state%column_global_y = current_state%pdd_y
    end if
    call timestep_callback_pdf_analysis_ptr(current_state)
    call timestep_callback_meanprofiles_ptr(current_state)
    call timestep_callback_fftsolver_ptr(current_state)
    call timestep_callback_iterativesolver_ptr(current_state)
    call timestep_callback_cfltest_ptr(current_state)
    !call execute_timestep_callbacks(current_state, group_descriptor%id)

  end subroutine timestep_whole_solver_group_contents

  subroutine timestep_whole_last_group_contents(current_state)!, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    !type(group_descriptor_type), intent(in) :: group_descriptor

    ! For print_debug_data, the column_global fields must match the requested coordinate.
    ! This is already handled for the timestep_column, but needs to be specially set for timestep_whole.
    !  This should not affect update_state_sitation_flags, as that is only used for timestep_column.
    if (current_state%print_debug_data) then
      current_state%column_global_x = current_state%pdd_x
      current_state%column_global_y = current_state%pdd_y
    end if
    call timestep_callback_conditional_diagnostics_whole_ptr(current_state)
    !!call timestep_callback_checkpointer(current_state)
    call timestep_callback_modelsynopsis_ptr(current_state)
    call timestep_callback_terminationcheck_ptr(current_state)
    call timestep_callback_iobridge_ptr(current_state)
    !call execute_timestep_callbacks(current_state, group_descriptor%id)s
  end subroutine timestep_whole_last_group_contents


    !> Performs timestepping for a group of components on a per column basis. Each component in the group is executed
  !! for every column.
  !! @param current_state The current model state
  !! @param group_descriptor Description of the group of components to execu
  subroutine timestep_column_subgrid_group_contents(current_state)!, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    !type(group_descriptor_type), intent(in) :: group_descriptor
    current_state%column_global_x=current_state%local_grid%start(X_INDEX) - current_state%local_grid%halo_size(X_INDEX)
    current_state%column_local_x=1
    do while (current_state%column_global_x .le. &
         current_state%local_grid%end(X_INDEX)+current_state%local_grid%halo_size(X_INDEX))
      current_state%column_global_y = current_state%local_grid%start(Y_INDEX) - current_state%local_grid%halo_size(Y_INDEX)
      current_state%column_local_y=1
      do while (current_state%column_global_y .le. &
           current_state%local_grid%end(Y_INDEX)+current_state%local_grid%halo_size(Y_INDEX))
        call update_state_sitation_flags(current_state)
        !call execute_timestep_callbacks(current_state, group_descriptor%id)
        call timestep_callback_lowerbc_ptr(current_state)
        call timestep_callback_smagorinsky_ptr(current_state)
        current_state%column_global_y = current_state%column_global_y + 1
        current_state%column_local_y = current_state%column_local_y + 1
      end do
      current_state%column_global_x = current_state%column_global_x + 1
      current_state%column_local_x = current_state%column_local_x + 1
    end do
  end subroutine timestep_column_subgrid_group_contents

  subroutine timestep_column_dynamics_group_contents(current_state)!, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    !type(group_descriptor_type), intent(in) :: group_descriptor
    current_state%column_global_x=current_state%local_grid%start(X_INDEX) - current_state%local_grid%halo_size(X_INDEX)
    current_state%column_local_x=1
    do while (current_state%column_global_x .le. &
         current_state%local_grid%end(X_INDEX)+current_state%local_grid%halo_size(X_INDEX))
      current_state%column_global_y = current_state%local_grid%start(Y_INDEX) - current_state%local_grid%halo_size(Y_INDEX)
      current_state%column_local_y=1
      do while (current_state%column_global_y .le. &
           current_state%local_grid%end(Y_INDEX)+current_state%local_grid%halo_size(Y_INDEX))
        call update_state_sitation_flags(current_state)
        !call execute_timestep_callbacks(current_state, group_descriptor%id)

        call timestep_callback_pwadvection_ptr(current_state)
        call timestep_callback_tvdadvection_ptr(current_state)
        call timestep_callback_thadvection_ptr(current_state)
        call timestep_callback_diffusion_ptr(current_state)
        call timestep_callback_viscosity_ptr(current_state)

        call timestep_callback_coriolis_ptr(current_state)
        call timestep_callback_buoyancy_ptr(current_state)
        call timestep_callback_damping_ptr(current_state)
        call timestep_callback_forcing_ptr(current_state)
        call timestep_callback_socrates_ptr(current_state)
        call timestep_callback_set_consistent_lowbc_ptr(current_state)
        call timestep_callback_simplecloud_ptr(current_state)
        call timestep_callback_tank_experiments_ptr(current_state)
        call timestep_callback_casim_ptr(current_state)

        call timestep_callback_diverr_ptr(current_state)
        call timestep_callback_pressuresource_ptr(current_state)!!!! attention de la position du calcul de pression
        call timestep_callback_profile_diagnostics_ptr(current_state)
        call timestep_callback_subgrid_profile_diagnostics_ptr(current_state)
        call timestep_callback_casim_profile_dgs_ptr(current_state)
        call timestep_callback_scalar_diagnostics_ptr(current_state)
        call timestep_callback_stepfields_ptr(current_state)
        if (current_state%timestep .lt. 2) then
        end if
        current_state%column_global_y = current_state%column_global_y + 1
        current_state%column_local_y = current_state%column_local_y + 1
      end do
      current_state%column_global_x = current_state%column_global_x + 1
      current_state%column_local_x = current_state%column_local_x + 1
    end do
  end subroutine timestep_column_dynamics_group_contents

  subroutine timestep_column_pressure_terms_group_contents(current_state)!, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    !type(group_descriptor_type), intent(in) :: group_descriptor
    current_state%column_global_x=current_state%local_grid%start(X_INDEX) - current_state%local_grid%halo_size(X_INDEX)
    current_state%column_local_x=1
    do while (current_state%column_global_x .le. &
         current_state%local_grid%end(X_INDEX)+current_state%local_grid%halo_size(X_INDEX))
      current_state%column_global_y = current_state%local_grid%start(Y_INDEX) - current_state%local_grid%halo_size(Y_INDEX)
      current_state%column_local_y=1
      do while (current_state%column_global_y .le. &
           current_state%local_grid%end(Y_INDEX)+current_state%local_grid%halo_size(Y_INDEX))
        call update_state_sitation_flags(current_state)
        !call execute_timestep_callbacks(current_state, group_descriptor%id)
        call timestep_callback_pstep_ptr(current_state)
        call timestep_callback_swapsmooth_ptr(current_state)
        call timestep_callback_diagnostics_column_ptr(current_state)
        call timestep_callback_diagnostics_3d_ptr(current_state)
        current_state%column_global_y = current_state%column_global_y + 1
        current_state%column_local_y = current_state%column_local_y + 1
      end do
      current_state%column_global_x = current_state%column_global_x + 1
      current_state%column_local_x = current_state%column_local_x + 1
    end do
  end subroutine timestep_column_pressure_terms_group_contents

  !> Updates the states situation flags for easy retrieval in the components that are
  !! run per timestep
  !! @param state The current model state
  subroutine update_state_sitation_flags(current_state)
    type(model_state_type), intent(inout) :: current_state

    current_state%first_timestep_column = (current_state%column_local_x == 1 .and. current_state%column_local_y == 1)
    current_state%last_timestep_column = (current_state%column_global_x == &
         current_state%local_grid%end(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) .and. &
         current_state%column_global_y == current_state%local_grid%end(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX))

    current_state%first_nonhalo_timestep_column = (current_state%column_local_x == current_state%local_grid%halo_size(X_INDEX)+1 &
       .and. current_state%column_local_y == current_state%local_grid%halo_size(Y_INDEX)+1)

    current_state%halo_column = current_state%column_local_y .lt. current_state%local_grid%local_domain_start_index(Y_INDEX) .or.&
         current_state%column_local_x .lt. current_state%local_grid%local_domain_start_index(X_INDEX) .or.&
         current_state%column_local_y .gt. current_state%local_grid%local_domain_end_index(Y_INDEX) .or.&
         current_state%column_local_x .gt. current_state%local_grid%local_domain_end_index(X_INDEX)
  end subroutine update_state_sitation_flags

  !> Updates the diagnostic sampling flag for the new timestep
  !! @param state The current model state
  subroutine handle_sampling(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i

    current_state%diagnostic_sample_timestep = .false.
    current_state%sampling(:)%active = .false.
    current_state%radiation_timestep = .false.  ! for computation timing under time_basis

    if (.not. current_state%only_compute_on_sample_timestep) then
      ! always compute the diagnostic in this case
      current_state%diagnostic_sample_timestep = .true.
    end if

    ! The following three cases will only compute diangnostics at requested intervals.
    ! However, it does ALL diagnostics regardless of specific request.
    ! MONC isn't STATSH-smart...though radiation diagnostics come close
    if (current_state%time_basis) then
      ! enable calculations and sampling at specified step only
      ! (at sampling time interval, which is also an output or write interval)
      do i=1, size(current_state%sampling)
        if (current_state%timestep .eq. current_state%sampling(i)%next_step) then
          if (current_state%sampling(i)%radiation) then
            ! Only possible when socrates_enabled and radiation_interval .gt. 0 (iobridge)
            ! Permits radiation without needing to do all diagnostics
            ! Never set %active for the %radiation case - does not denote a iob data_definition
            current_state%radiation_timestep = .true.
          else
            current_state%diagnostic_sample_timestep = .true.
            current_state%sampling(i)%active = .true.
          end if
        end if
      end do
    else ! timestep basis or force_output_on_interval
      ! enable radiation calculation
      if (current_state%socrates_enabled .and. radiation_interval .gt. 0) then
        if (mod(current_state%timestep, radiation_interval) == 0) &
             current_state%radiation_timestep = .true.
      end if
      ! enable diagnostic calculation and sampling on the sampling timestep interval.
      !do i=1,size(current_state%sampling)
      !  if (mod(current_state%timestep, current_state%sampling(i)%interval) == 0) then
      !    current_state%diagnostic_sample_timestep = .true.
      !    current_state%sampling(i)%active = .true.
      !  end if
      !end do
    end if
  end subroutine handle_sampling

end module timestepper_mod
