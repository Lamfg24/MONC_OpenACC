!> MONC component registry
!!
!! Supports management of components. Each stage is called via the registry which will execute the
!! installed callback hooks for that stage in order.
module registry_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use collections_mod, only : list_type, hashmap_type, map_type, iterator_type, c_size, c_generic_at, c_key_at, c_get_integer, &
       c_get_string, c_get_generic, c_remove, c_put_generic, c_put_string, c_put_integer, c_put_real, c_is_empty, &
       c_contains, c_add_generic, c_add_string, c_free, c_get_iterator, c_has_next, c_next_mapentry, mapentry_type
  use monc_component_mod, only : component_descriptor_type, component_field_value_type, component_field_information_type, &
       FINALISATION_PRIORITY_INDEX, INIT_PRIORITY_INDEX, TIMESTEP_PRIORITY_INDEX, &
       pointer_wrapper_value_type, pointer_wrapper_info_type, component_descriptor_type_v1_array
  use conversions_mod, only : conv_to_string
  use state_mod, only : model_state_type
  use optionsdatabase_mod, only : options_has_key, options_get_string, options_get_logical, options_get_array_size
  use logging_mod, only : LOG_INFO, LOG_ERROR, LOG_WARN, log_master_log, log_is_master
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use clearsourceterms_mod ! contains the no_component function
  use decomposition_mod
  use checkpointer_mod
  use simplesetup_mod
  use gridmanager_mod
  use meanprofiles_mod
  use swapsmooth_mod
  use terminationcheck_mod
  use simplecloud_mod
  use casim_mod
  use coriolis_mod
  use buoyancy_mod
  use damping_mod
  use cfltest_mod
  use diverr_mod
  use haloswapper_mod
  use fftsolver_mod
  use iterativesolver_mod
  use setfluxlook_mod
  use lowerbc_mod
  use pressuresource_mod
  use pwadvection_mod
  use diffusion_mod
  use set_consistent_lowbc_mod
  use viscosity_mod
  use smagorinsky_mod
  use stepfields_mod
  use steppingdirection_mod
  use tvdadvection_mod
  use modelsynopsis_mod
  use thadvection_mod
  use randomnoise_mod
  use forcing_mod
  use diagnostics_3d_mod
  use profile_diagnostics_mod
  use casim_profile_dgs_mod
  use conditional_diagnostics_column_mod
  use conditional_diagnostics_whole_mod
  use pdf_analysis_mod
  use subgrid_profile_diagnostics_mod
  use scalar_diagnostics_mod
  use pstep_mod
  use tank_experiments_mod
  use drybl_mod
  use iobridge_mod
  use socrates_couple_mod
  implicit none

#ifndef TEST_MODE
  private
#endif

  abstract interface
    subroutine component_init(current_state)
        use state_mod, only : model_state_type
        type(model_state_type), target, intent(inout) :: current_state
    end subroutine component_init
  end interface

  procedure (component_init), pointer :: init_callback_decomposition_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_simplesetup_ptr => no_component
  procedure (component_init), pointer :: initialise_callback_gridmanager_ptr => no_component
  procedure (component_init), pointer :: init_callback_meanprofiles_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_swapsmooth_ptr => no_component
  procedure (component_init), pointer :: init_callback_terminationcheck_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_simplecloud_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_casim_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_coriolis_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_buoyancy_ptr => no_component
  procedure (component_init), pointer :: init_callback_damping_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_cfltest_ptr => no_component
  procedure (component_init), pointer :: init_callback_diverr_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_haloswapper_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_fftsolver_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_iterativesolver_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_setfluxlook_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_lowerbc_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_pressuresource_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_pwadvection_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_diffusion_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_set_consistent_lowbc_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_viscosity_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_smagorinsky_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_stepfields_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_steppingdirection_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_tvdadvection_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_modelsynopsis_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_thadvection_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_randomnoise_ptr => no_component
  procedure (component_init), pointer :: init_callback_forcing_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_diagnostics_3d_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_profile_diagnostics_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_casim_profile_dgs_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_conditional_diagnostics_column_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_conditional_diagnostics_whole_ptr => no_component
  procedure (component_init), pointer :: init_callback_pdf_analysis_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_subgrid_profile_diagnostics_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_scalar_diagnostics_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_pstep_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_tank_experiments_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_drybl_ptr => no_component
  procedure (component_init), pointer :: init_callback_iobridge_ptr => no_component
  procedure (component_init), pointer :: initialisation_callback_socrates_ptr => no_component

  abstract interface
    subroutine component_final(current_state)
        use state_mod, only : model_state_type
        type(model_state_type), target, intent(inout) :: current_state
    end subroutine component_final
  end interface

  procedure (component_final), pointer :: finalise_callback_gridmanager_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_meanprofiles_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_coriolis_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_simplecloud_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_buoyancy_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_damping_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_diverr_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_haloswapper_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_fftsolver_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_iterativesolver_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_lowerbc_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_pwadvection_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_diffusion_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_viscosity_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_smagorinsky_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_stepfields_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_tvdadvection_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_modelsynopsis_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_thadvection_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_forcing_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_profile_diagnostics_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_conditional_diagnostics_column_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_conditional_diagnostics_whole_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_pdf_analysis_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_pstep_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_tank_experiments_ptr => no_component
  procedure (component_final), pointer :: finalisation_callback_socrates_ptr => no_component


  public register_component, execute_initialisation_callbacks, execute_finalisation_callbacks
contains

  !> Initialises the registry with the provided configuration file
  !! @param configurationFileName The filename of the configuration file to parse
  subroutine init_registry(options_database_real, options_database_string)
    !type(hashmap_type), intent(inout) :: options_database
    real(kind=DEFAULT_PRECISION), dimension(1200,500), intent(inout) :: options_database_real
    character(len=STRING_LENGTH), dimension(1200,2), intent(inout) :: options_database_string

    !call read_group_configurations(options_database_real, options_database_string)
    !call read_initialisation_and_finalisation_orders(options_database)
    !allocate(timestep_callbacks(c_size(group_descriptors)))
  end subroutine init_registry

  !> Will deregister all components and free up the registry data structures. This can either be called
  !! at the end of execution to clean memory up or used to clear the registry

  !> Will register a component and install the nescesary callback hooks
  !! @param descriptor The component descriptor and a separate copy of this it stored as reference
  !subroutine register_component(options_database, descriptor)
  subroutine register_component(current_state)
    type(model_state_type), intent(inout) :: current_state
    integer :: iter
    logical :: logicnum

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "decomposition_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) init_callback_decomposition_ptr => init_callback_decomposition
      else if (current_state%options_database_string(iter,1) .eq. "simplesetup_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_simplesetup_ptr => initialisation_callback_simplesetup
      else if (current_state%options_database_string(iter,1) .eq. "gridmanager_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialise_callback_gridmanager_ptr => initialise_callback_gridmanager
          finalise_callback_gridmanager_ptr => finalise_callback_gridmanager
        end if
      else if (current_state%options_database_string(iter,1) .eq. "mean_profiles_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          init_callback_meanprofiles_ptr => init_callback_meanprofiles
          finalisation_callback_meanprofiles_ptr => finalisation_callback_meanprofiles
        end if
      else if (current_state%options_database_string(iter,1) .eq. "swap_smooth_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_swapsmooth_ptr => initialisation_callback_swapsmooth
      else if (current_state%options_database_string(iter,1) .eq. "termination_check_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) init_callback_terminationcheck_ptr => init_callback_terminationcheck
      else if (current_state%options_database_string(iter,1) .eq. "simplecloud_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_simplecloud_ptr => initialisation_callback_simplecloud
          finalisation_callback_simplecloud_ptr => finalisation_callback_simplecloud
        end if
      else if (current_state%options_database_string(iter,1) .eq. "casim_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_casim_ptr => initialisation_callback_casim
      else if (current_state%options_database_string(iter,1) .eq. "coriolis_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_coriolis_ptr => initialisation_callback_coriolis
          finalisation_callback_coriolis_ptr => finalisation_callback_coriolis
        end if
      else if (current_state%options_database_string(iter,1) .eq. "buoyancy_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_buoyancy_ptr => initialisation_callback_buoyancy
          finalisation_callback_buoyancy_ptr => finalisation_callback_buoyancy
        end if
      else if (current_state%options_database_string(iter,1) .eq. "damping_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          init_callback_damping_ptr => init_callback_damping
          finalisation_callback_damping_ptr => finalisation_callback_damping
        end if
      else if (current_state%options_database_string(iter,1) .eq. "cfltest_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_cfltest_ptr => initialisation_callback_cfltest
      else if (current_state%options_database_string(iter,1) .eq. "diverr_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          init_callback_diverr_ptr => init_callback_diverr
          finalisation_callback_diverr_ptr => finalisation_callback_diverr
        end if
      else if (current_state%options_database_string(iter,1) .eq. "halo_swapper_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_haloswapper_ptr => initialisation_callback_haloswapper
          finalisation_callback_haloswapper_ptr => finalisation_callback_haloswapper
        end if
      else if (current_state%options_database_string(iter,1) .eq. "fftsolver_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_fftsolver_ptr => initialisation_callback_fftsolver
          finalisation_callback_fftsolver_ptr => finalisation_callback_fftsolver
        end if
      else if (current_state%options_database_string(iter,1) .eq. "iterativesolver_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_iterativesolver_ptr => initialisation_callback_iterativesolver
          finalisation_callback_iterativesolver_ptr => finalisation_callback_iterativesolver
        end if
      else if (current_state%options_database_string(iter,1) .eq. "setfluxlook_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_setfluxlook_ptr => initialisation_callback_setfluxlook
      else if (current_state%options_database_string(iter,1) .eq. "lower_bc_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_lowerbc_ptr => initialisation_callback_lowerbc
          finalisation_callback_lowerbc_ptr => finalisation_callback_lowerbc
        end if
      else if (current_state%options_database_string(iter,1) .eq. "pressure_source_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_pressuresource_ptr => initialisation_callback_pressuresource
      else if (current_state%options_database_string(iter,1) .eq. "pw_advection_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_pwadvection_ptr => initialisation_callback_pwadvection
          finalisation_callback_pwadvection_ptr => finalisation_callback_pwadvection
        end if
      else if (current_state%options_database_string(iter,1) .eq. "diffusion_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_diffusion_ptr => initialisation_callback_diffusion
          finalisation_callback_diffusion_ptr => finalisation_callback_diffusion
        end if
      else if (current_state%options_database_string(iter,1) .eq. "set_consistent_lowbc_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_set_consistent_lowbc_ptr => initialisation_callback_set_consistent_lowbc
      else if (current_state%options_database_string(iter,1) .eq. "viscosity_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_viscosity_ptr => initialisation_callback_viscosity
          finalisation_callback_viscosity_ptr => finalisation_callback_viscosity
        end if
      else if (current_state%options_database_string(iter,1) .eq. "smagorinsky_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_smagorinsky_ptr => initialisation_callback_smagorinsky
          finalisation_callback_smagorinsky_ptr => finalisation_callback_smagorinsky
        end if
      else if (current_state%options_database_string(iter,1) .eq. "stepfields_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_stepfields_ptr => initialisation_callback_stepfields
          finalisation_callback_stepfields_ptr => finalisation_callback_stepfields
        end if
      else if (current_state%options_database_string(iter,1) .eq. "stepping_direction_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_steppingdirection_ptr => initialisation_callback_steppingdirection
      else if (current_state%options_database_string(iter,1) .eq. "tvd_advection_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_tvdadvection_ptr => initialisation_callback_tvdadvection
          finalisation_callback_tvdadvection_ptr => finalisation_callback_tvdadvection
        end if
      else if (current_state%options_database_string(iter,1) .eq. "model_synopsis_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_modelsynopsis_ptr => initialisation_callback_modelsynopsis
          finalisation_callback_modelsynopsis_ptr => finalisation_callback_modelsynopsis
        end if
      else if (current_state%options_database_string(iter,1) .eq. "th_advection_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_thadvection_ptr => initialisation_callback_thadvection
          finalisation_callback_thadvection_ptr => finalisation_callback_thadvection
        end if
      else if (current_state%options_database_string(iter,1) .eq. "randomnoise_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_randomnoise_ptr => initialisation_callback_randomnoise
      else if (current_state%options_database_string(iter,1) .eq. "forcing_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          init_callback_forcing_ptr => init_callback_forcing
          finalisation_callback_forcing_ptr => finalisation_callback_forcing
        end if
      else if (current_state%options_database_string(iter,1) .eq. "diagnostics_3d_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_diagnostics_3d_ptr => initialisation_callback_diagnostics_3d
      else if (current_state%options_database_string(iter,1) .eq. "profile_diagnostics_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_profile_diagnostics_ptr => initialisation_callback_profile_diagnostics
          finalisation_callback_profile_diagnostics_ptr => finalisation_callback_profile_diagnostics
        end if
      else if (current_state%options_database_string(iter,1) .eq. "casim_profile_dgs_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_casim_profile_dgs_ptr => initialisation_callback_casim_profile_dgs
      else if (current_state%options_database_string(iter,1) .eq. "conditional_diagnostics_column_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_conditional_diagnostics_column_ptr => initialisation_callback_conditional_diagnostics_column
          finalisation_callback_conditional_diagnostics_column_ptr => finalisation_callback_conditional_diagnostics_column
        end if
      else if (current_state%options_database_string(iter,1) .eq. "conditional_diagnostics_whole_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_conditional_diagnostics_whole_ptr => &
                        initialisation_callback_conditional_diagnostics_whole
          finalisation_callback_conditional_diagnostics_whole_ptr => finalisation_callback_conditional_diagnostics_whole
        end if
      else if (current_state%options_database_string(iter,1) .eq. "pdf_analysis_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          init_callback_pdf_analysis_ptr => init_callback_pdf_analysis
          finalisation_callback_pdf_analysis_ptr => finalisation_callback_pdf_analysis
        end if
      else if (current_state%options_database_string(iter,1) .eq. "subgrid_profile_diagnostics_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_subgrid_profile_diagnostics_ptr => &
                                                               initialisation_callback_subgrid_profile_diagnostics
      else if (current_state%options_database_string(iter,1) .eq. "scalar_diagnostics_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) initialisation_callback_scalar_diagnostics_ptr => initialisation_callback_scalar_diagnostics
      else if (current_state%options_database_string(iter,1) .eq. "pstep_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_pstep_ptr => initialisation_callback_pstep
          finalisation_callback_pstep_ptr => finalisation_callback_pstep
        end if
      else if (current_state%options_database_string(iter,1) .eq. "tank_experiments_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_tank_experiments_ptr => initialisation_callback_tank_experiments
          finalisation_callback_tank_experiments_ptr => finalisation_callback_tank_experiments
        end if
      else if (current_state%options_database_string(iter,1) .eq. "drybl_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_drybl_ptr => initialisation_callback_drybl
        end if
      else if (current_state%options_database_string(iter,1) .eq. "socrates_couple_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          initialisation_callback_socrates_ptr => initialisation_callback_socrates
          finalisation_callback_socrates_ptr => finalisation_callback_socrates
        end if
      else if (current_state%options_database_string(iter,1) .eq. "iobridge_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) init_callback_iobridge_ptr => init_callback_iobridge
      end if
    end do

!     allocate(registry_descriptor, source=descriptor) ! Make copy of the descriptor (which registry might modify)
!
!     if (options_has_key(options_database, trim(descriptor%name)//"_enabled")) then
!       component_enabled=options_get_logical(options_database, trim(descriptor%name)//"_enabled")
!     else
!       component_enabled=.false.
!       call log_master_log(LOG_WARN, "No enabled configuration for component "//trim(descriptor%name)//" therefore disabling this")
!     end if
!
!     if (component_enabled) then
!       if (c_contains(component_groups, trim(descriptor%name))) then
!         group_name=c_get_string(component_groups, trim(descriptor%name))
!         call load_callback_hooks(registry_descriptor, group_name)
!       else
!         call load_callback_hooks(registry_descriptor)
!       end if
!       call load_published_fields(descriptor)
!     end if
!
!     description_data => registry_descriptor
!     call c_put_generic(component_descriptions, descriptor%name, description_data, .false.)
  end subroutine register_component

  !> Calls all initialisation callbacks with the specified state
  !! @param currentState The current model state which may (and often is) modified
  subroutine execute_initialisation_callbacks(current_state)
    type(model_state_type), intent(inout) :: current_state
    !type(component_descriptor_type_v1_array), intent(inout) :: component_descriptions

    call init_callback_decomposition_ptr(current_state) ! check
    !call initialisation_callback_checkpointer(current_state) ! check
    call initialisation_callback_simplesetup_ptr(current_state) ! check
    call initialise_callback_gridmanager_ptr(current_state) ! check
    call init_callback_meanprofiles_ptr(current_state) ! check
    call initialisation_callback_swapsmooth_ptr(current_state) ! check
    call init_callback_terminationcheck_ptr(current_state) ! check
    call initialisation_callback_simplecloud_ptr(current_state)
    call initialisation_callback_casim_ptr(current_state) ! check
    call initialisation_callback_coriolis_ptr(current_state) ! check
    call initialisation_callback_buoyancy_ptr(current_state) ! check
    call init_callback_damping_ptr(current_state) ! check
    call initialisation_callback_cfltest_ptr(current_state) ! check
    call init_callback_diverr_ptr(current_state) ! check
    call initialisation_callback_haloswapper_ptr(current_state) ! check
    call initialisation_callback_fftsolver_ptr(current_state)
    call initialisation_callback_iterativesolver_ptr(current_state) ! check
    call initialisation_callback_setfluxlook_ptr(current_state) ! check
    call initialisation_callback_lowerbc_ptr(current_state) ! check
    call initialisation_callback_pressuresource_ptr(current_state) ! check
    call initialisation_callback_pwadvection_ptr(current_state) ! check
    call initialisation_callback_diffusion_ptr(current_state) ! check
    call initialisation_callback_set_consistent_lowbc_ptr(current_state) ! check
    call initialisation_callback_viscosity_ptr(current_state) ! check
    call initialisation_callback_smagorinsky_ptr(current_state) ! check
    call initialisation_callback_stepfields_ptr(current_state) ! check
    call initialisation_callback_steppingdirection_ptr(current_state) ! check
    call initialisation_callback_tvdadvection_ptr(current_state) ! check
    call initialisation_callback_modelsynopsis_ptr(current_state) ! check
    call initialisation_callback_thadvection_ptr(current_state) ! check
    call initialisation_callback_randomnoise_ptr(current_state) ! check
    call init_callback_forcing_ptr(current_state)
    call initialisation_callback_diagnostics_3d_ptr(current_state) ! check
    call initialisation_callback_profile_diagnostics_ptr(current_state) ! check
    call initialisation_callback_casim_profile_dgs_ptr(current_state) ! check
    call initialisation_callback_conditional_diagnostics_column_ptr(current_state) ! check
    call initialisation_callback_conditional_diagnostics_whole_ptr(current_state) ! check
    call init_callback_pdf_analysis_ptr(current_state) ! check
    call initialisation_callback_subgrid_profile_diagnostics_ptr(current_state) ! check
    call initialisation_callback_scalar_diagnostics_ptr(current_state) ! check
    call initialisation_callback_pstep_ptr(current_state) ! check
    call initialisation_callback_tank_experiments_ptr(current_state)
    call initialisation_callback_drybl_ptr(current_state)
    call init_callback_iobridge_ptr(current_state)
    call initialisation_callback_socrates_ptr(current_state)
    !!call execute_callbacks(init_callbacks, current_state, "initialisation_callback")
  end subroutine execute_initialisation_callbacks


  !> Calls all finalisation callbacks with the specified state
  !! @param currentState The current model state which may (and often is) modified
  subroutine execute_finalisation_callbacks(current_state)
    type(model_state_type), intent(inout) :: current_state
    !call finalisation_callback_checkpointer(current_state)
    call finalise_callback_gridmanager_ptr(current_state)
    call finalisation_callback_meanprofiles_ptr(current_state)
    call finalisation_callback_coriolis_ptr(current_state)
    call finalisation_callback_buoyancy_ptr(current_state)
    call finalisation_callback_damping_ptr(current_state)
    call finalisation_callback_diverr_ptr(current_state)
    call finalisation_callback_haloswapper_ptr(current_state)
    call finalisation_callback_fftsolver_ptr(current_state)
    call finalisation_callback_iterativesolver_ptr(current_state)
    call finalisation_callback_lowerbc_ptr(current_state)
    call finalisation_callback_pwadvection_ptr(current_state)
    call finalisation_callback_diffusion_ptr(current_state)
    call finalisation_callback_viscosity_ptr(current_state)
    call finalisation_callback_smagorinsky_ptr(current_state)
    call finalisation_callback_stepfields_ptr(current_state)
    call finalisation_callback_tvdadvection_ptr(current_state)
    call finalisation_callback_modelsynopsis_ptr(current_state)
    call finalisation_callback_thadvection_ptr(current_state)
    call finalisation_callback_forcing_ptr(current_state)
    call finalisation_callback_profile_diagnostics_ptr(current_state)
    call finalisation_callback_conditional_diagnostics_column_ptr(current_state)
    call finalisation_callback_conditional_diagnostics_whole_ptr(current_state)
    call finalisation_callback_pdf_analysis_ptr(current_state)
    call finalisation_callback_pstep_ptr(current_state)
    call finalisation_callback_tank_experiments_ptr(current_state)
    call finalisation_callback_socrates_ptr(current_state)

    !!call execute_callbacks(finalisation_callbacks, current_state, "finalisation_callback")
  end subroutine execute_finalisation_callbacks


  
end module registry_mod
