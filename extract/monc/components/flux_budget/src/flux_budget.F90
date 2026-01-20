!> Flux budget component which produces diagnostic data for the flux aspects of the model
module flux_budget_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use state_mod, only : model_state_type
  use science_constants_mod, only: G, ratio_mol_wts
  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: some_theta_flux_diagnostics_enabled, some_q_flux_diagnostics_enabled, some_uw_vw_diagnostics_enabled, &
       some_prognostic_budget_diagnostics_enabled, some_thetal_diagnostics_enabled, some_mse_diagnostics_enabled, &
       some_qt_diagnostics_enabled, some_tke_diagnostics_enabled
       
  public initialisation_callback_flux_budget, timestep_callback_flux_budget, finalisation_callback_flux_budget

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The flux budget component descriptor
!   type(component_descriptor_type_v1) function flux_budget_get_descriptor()
!     type(iterator_type) :: iterator
!     type(mapentry_type) :: mapentry
!     integer :: current_index, total_number_published_fields
!
!     flux_budget_get_descriptor%name="flux_budget"
!     flux_budget_get_descriptor%version=0.1
!     !flux_budget_get_descriptor%initialisation=>initialisation_callback
!     !flux_budget_get_descriptor%timestep=>timestep_callback
!     !flux_budget_get_descriptor%finalisation=>finalisation_callback
!
!     !flux_budget_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
!     !flux_budget_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
!     call populate_field_names()
!     total_number_published_fields=c_size(heat_flux_fields) + c_size(q_flux_fields) + c_size(uw_vw_fields) +       &
!                                   c_size(prognostic_budget_fields) + c_size(thetal_fields) + c_size(mse_fields) + &
!                                   c_size(qt_fields) + c_size(scalar_fields) + c_size(tke_fields)
!     flux_budget_get_descriptor%published_fields_on_off = .true.
!     allocate(flux_budget_get_descriptor%published_fields(total_number_published_fields))
!
!     current_index=1
!     iterator=c_get_iterator(heat_flux_fields)
!     do while (c_has_next(iterator))
!       mapentry=c_next_mapentry(iterator)
!       flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
!       current_index=current_index+1
!     end do
!     iterator=c_get_iterator(q_flux_fields)
!     do while (c_has_next(iterator))
!       mapentry=c_next_mapentry(iterator)
!       flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
!       current_index=current_index+1
!     end do
!     iterator=c_get_iterator(uw_vw_fields)
!     do while (c_has_next(iterator))
!       mapentry=c_next_mapentry(iterator)
!       flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
!       current_index=current_index+1
!     end do
!     iterator=c_get_iterator(tke_fields)
!     do while (c_has_next(iterator))
!       mapentry=c_next_mapentry(iterator)
!       flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
!       current_index=current_index+1
!     end do
!     iterator=c_get_iterator(prognostic_budget_fields)
!     do while (c_has_next(iterator))
!       mapentry=c_next_mapentry(iterator)
!       flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
!       current_index=current_index+1
!     end do
!     iterator=c_get_iterator(thetal_fields)
!     do while (c_has_next(iterator))
!       mapentry=c_next_mapentry(iterator)
!       flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
!       current_index=current_index+1
!     end do
!     iterator=c_get_iterator(mse_fields)
!     do while (c_has_next(iterator))
!       mapentry=c_next_mapentry(iterator)
!       flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
!       current_index=current_index+1
!     end do
!     iterator=c_get_iterator(qt_fields)
!     do while (c_has_next(iterator))
!       mapentry=c_next_mapentry(iterator)
!       flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
!       current_index=current_index+1
!     end do
!     iterator=c_get_iterator(scalar_fields)
!     do while (c_has_next(iterator))
!       mapentry=c_next_mapentry(iterator)
!       flux_budget_get_descriptor%published_fields(current_index)=mapentry%key
!       current_index=current_index+1
!     end do
!   end function flux_budget_get_descriptor

  !> Initialisation call back
  !! @param current_state Current model state
  subroutine initialisation_callback_flux_budget(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call initialise_theta_flux_diagnostics(current_state)
    call initialise_q_flux_diagnostics(current_state)
    call initialise_uw_vw_diagnostics(current_state)
    call initialise_prognostic_budget_diagnostics(current_state)
    call initialise_thetal_diagnostics(current_state)
    call initialise_scalar_diagnostics(current_state)
    call initialise_tke_diagnostics(current_state)
    !call initialise_mse_diagnostics(current_state)
    !call initialise_qt_diagnostics(current_state)
  end subroutine initialisation_callback_flux_budget

  !> Initialises the heat flux diagnostic areas and enabled flags depending upon the configuration of the model
  !! @param current_state The current model state
  subroutine initialise_theta_flux_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    logical :: th_flux_term_enabled, th_tendency_term_enabled, th_diff_enabled, th_gradient_term_enabled, th_buoyancy_enabled

    th_flux_term_enabled = current_state%th%active .and. current_state%w%active
    th_tendency_term_enabled = current_state%th%active .and. current_state%w%active
    th_gradient_term_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) .and. &
                              current_state%w%active .and. current_state%th%active
    th_diff_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
                      .and. current_state%w%active .and. current_state%th%active
    th_buoyancy_enabled = current_state%buoyancy_enabled .and. current_state%th%active
    some_theta_flux_diagnostics_enabled = th_flux_term_enabled .or. th_buoyancy_enabled

    if (th_flux_term_enabled) then
      !call set_published_field_enabled_state(heat_flux_fields, "heat_flux_transport_local", .true.)
      allocate(current_state%th_flux_values(current_state%local_grid%size(Z_INDEX)))
    end if
    if (th_tendency_term_enabled) then
      !call set_published_field_enabled_state(heat_flux_fields, "heat_flux_tendency_local", .true.)
      allocate(current_state%th_tendency(current_state%local_grid%size(Z_INDEX)))
    end if
    if (th_diff_enabled) then
      !call set_published_field_enabled_state(heat_flux_fields, "heat_flux_dissipation_local", .true.)
      allocate(current_state%th_diff(current_state%local_grid%size(Z_INDEX)))
    end if
    if (th_gradient_term_enabled) then
      !call set_published_field_enabled_state(heat_flux_fields, "heat_flux_gradient_local", .true.)
      allocate(current_state%th_gradient(current_state%local_grid%size(Z_INDEX)))
    end if
    if (th_buoyancy_enabled) then
      !call set_published_field_enabled_state(heat_flux_fields, "heat_flux_buoyancy_local", .true.)
      allocate(current_state%th_buoyancy(current_state%local_grid%size(Z_INDEX)))
    end if
  end subroutine initialise_theta_flux_diagnostics

  !> Initialises the Q field flux diagnostic areas and enabled flags depending upon the configuration of the model
  !! @param current_state The current model state
  subroutine initialise_q_flux_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    logical :: q_flux_term_enabled, q_tendency_term_enabled, q_gradient_term_enabled, q_diff_enabled, &
       q_buoyancy_enabled

    q_flux_term_enabled = current_state%w%active
    q_tendency_term_enabled = current_state%w%active
    q_gradient_term_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) .and. &
                              current_state%w%active
    q_diff_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
                    .and. current_state%w%active
    q_buoyancy_enabled = current_state%buoyancy_enabled
    some_q_flux_diagnostics_enabled = q_flux_term_enabled .or. q_buoyancy_enabled

    if (q_flux_term_enabled) then
      !call set_published_field_enabled_state(q_flux_fields, "q_flux_transport_local", .true.)
      allocate(current_state%qv_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ql_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qr_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qi_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qg_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qs_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAitkenSolMass_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumSolMass_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumInsolMass_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseSolMass_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseDustMass_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nl_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nr_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ni_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ng_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ns_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAitkenSolNumber_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumSolNumber_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumInsolNumber_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseSolNumber_flux_values(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseDustnumber_flux_values(current_state%local_grid%size(Z_INDEX)))
    end if
    if (q_tendency_term_enabled) then
      !call set_published_field_enabled_state(q_flux_fields, "q_flux_tendency_local", .true.)
      allocate(current_state%qv_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ql_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qr_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qi_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qs_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qg_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAitkenSolMass_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumSolMass_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumInsolMass_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseSolMass_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseDustMass_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nl_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nr_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ni_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ns_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ng_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAitkenSolNumber_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumSolNumber_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumInsolNumber_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseSolNumber_tendency(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseDustnumber_tendency(current_state%local_grid%size(Z_INDEX)))
    end if
    if (q_gradient_term_enabled) then
      !call set_published_field_enabled_state(q_flux_fields, "q_flux_gradient_local", .true.)
      allocate(current_state%qv_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ql_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qr_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qi_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qs_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qg_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAitkenSolMass_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumSolMass_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumInsolMass_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseSolMass_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseDustMass_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nl_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nr_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ni_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ng_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ns_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAitkenSolNumber_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumSolNumber_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumInsolNumber_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseSolNumber_gradient(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseDustnumber_gradient(current_state%local_grid%size(Z_INDEX)))
    end if
    if (q_diff_enabled) then
      !call set_published_field_enabled_state(q_flux_fields, "q_flux_dissipation_local", .true.)
      allocate(current_state%qv_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ql_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qr_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qi_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qs_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qg_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAitkenSolMass_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumSolMass_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumInsolMass_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseSolMass_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseDustMass_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nl_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nr_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ni_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ng_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ns_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAitkenSolNumber_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumSolNumber_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumInsolNumber_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseSolNumber_diff(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseDustnumber_diff(current_state%local_grid%size(Z_INDEX)))
    end if
    if (q_buoyancy_enabled) then
      !call set_published_field_enabled_state(q_flux_fields, "q_flux_buoyancy_local", .true.)
      allocate(current_state%qv_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ql_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qr_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qi_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qs_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qg_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAitkenSolMass_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumSolMass_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qAccumInsolMass_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseSolMass_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%qCoarseDustMass_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nl_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nr_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ni_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ng_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%ns_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAitkenSolNumber_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumSolNumber_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nAccumInsolNumber_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseSolNumber_buoyancy(current_state%local_grid%size(Z_INDEX)))
      allocate(current_state%nCoarseDustnumber_buoyancy(current_state%local_grid%size(Z_INDEX)))
    end if
  end subroutine initialise_q_flux_diagnostics

  !> Initialises the UW and VW diagnostics
  !! @param current_state Current model state
  subroutine initialise_uw_vw_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: column_size
    logical :: uw_advection_term_enabled, vw_advection_term_enabled, uw_viscosity_term_enabled, &
       vw_viscosity_term_enabled, uw_buoyancy_term_enabled, vw_buoyancy_term_enabled, uw_tendency_term_enabled, &
       vw_tendency_term_enabled, uw_w_term_enabled, vw_w_term_enabled

    uw_advection_term_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
                                .and. current_state%w%active
    vw_advection_term_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
                                .and. current_state%w%active
    uw_viscosity_term_enabled = current_state%viscosity_enabled .and. current_state%w%active
    vw_viscosity_term_enabled = current_state%viscosity_enabled .and. current_state%w%active
    uw_buoyancy_term_enabled = current_state%buoyancy_enabled
    vw_buoyancy_term_enabled = current_state%buoyancy_enabled
    uw_tendency_term_enabled = current_state%w%active .and. current_state%u%active
    vw_tendency_term_enabled = current_state%w%active .and. current_state%v%active
    uw_w_term_enabled = current_state%w%active .and. current_state%u%active
    vw_w_term_enabled = current_state%w%active .and. current_state%v%active
    some_uw_vw_diagnostics_enabled = uw_buoyancy_term_enabled .or. vw_buoyancy_term_enabled .or. uw_w_term_enabled &
      .or. vw_w_term_enabled .or. uw_advection_term_enabled .or. vw_advection_term_enabled .or. uw_viscosity_term_enabled .or. &
      vw_viscosity_term_enabled

    column_size=current_state%local_grid%size(Z_INDEX)

    if (uw_advection_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "uw_advection_local", .true.)
      allocate(current_state%uw_advection(column_size))
    end if
    if (vw_advection_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "vw_advection_local", .true.)
      allocate(current_state%vw_advection(column_size))
    end if
    if (uw_viscosity_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "uw_viscosity_local", .true.)
      allocate(current_state%uw_viscosity(column_size))
    end if
    if (vw_viscosity_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "vw_viscosity_local", .true.)
      allocate(current_state%vw_viscosity(column_size))
    end if
    if (uw_buoyancy_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "uw_buoyancy_local", .true.)
      allocate(current_state%uw_buoyancy(column_size))
    end if
    if (vw_buoyancy_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "vw_buoyancy_local", .true.)
      allocate(current_state%vw_buoyancy(column_size))
    end if
    if (uw_tendency_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "uw_tendency_local", .true.)
      allocate(current_state%uw_tendency(column_size))
    end if
    if (vw_tendency_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "vw_tendency_local", .true.)
      allocate(current_state%vw_tendency(column_size))
    end if
    if (uw_w_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "uw_w_local", .true.)
      allocate(current_state%uw_w(column_size))
    end if
    if (vw_w_term_enabled) then
      !call set_published_field_enabled_state(uw_vw_fields, "vw_w_local", .true.)
      allocate(current_state%vw_w(column_size))
    end if
  end subroutine initialise_uw_vw_diagnostics

  !> Initialises the prognostic (uu, vv, ww) budget diagnostics
  !! @param current_state The current model state
  subroutine initialise_prognostic_budget_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: column_size
    logical :: tu_su_enabled, uu_advection_enabled, uu_viscosity_enabled, wu_u_enabled, tv_sv_enabled, vv_advection_enabled, &
         vv_viscosity_enabled, wv_v_enabled, tw_sw_enabled, ww_advection_enabled, ww_viscosity_enabled, ww_buoyancy_enabled

    tu_su_enabled = current_state%u%active
    uu_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
                            .and. current_state%u%active
    uu_viscosity_enabled = current_state%viscosity_enabled .and. current_state%u%active
    wu_u_enabled = current_state%u%active .and. current_state%w%active
    tv_sv_enabled = current_state%v%active
    vv_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
                            .and. current_state%v%active
    vv_viscosity_enabled = current_state%viscosity_enabled .and. current_state%v%active
    wv_v_enabled = current_state%v%active .and. current_state%w%active
    tw_sw_enabled = current_state%w%active
    ww_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
                            .and. current_state%w%active
    ww_viscosity_enabled = current_state%viscosity_enabled .and. current_state%w%active
    ww_buoyancy_enabled = current_state%buoyancy_enabled .and. current_state%w%active
    some_prognostic_budget_diagnostics_enabled=tu_su_enabled .or. tv_sv_enabled .or. tw_sw_enabled

    column_size=current_state%local_grid%size(Z_INDEX)

    if (tu_su_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "tu_su_local", .true.)
      allocate(current_state%tu_su(column_size))
    end if
    if (uu_advection_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "uu_advection_local", .true.)
      allocate(current_state%uu_advection(column_size))
    end if
    if (uu_viscosity_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "uu_viscosity_local", .true.)
      allocate(current_state%uu_viscosity(column_size))
    end if
    if (wu_u_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "wu_u_local", .true.)
      allocate(current_state%wu_u(column_size))
    end if
    if (tv_sv_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "tv_sv_local", .true.)
      allocate(current_state%tv_sv(column_size))
    end if
    if (vv_advection_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "vv_advection_local", .true.)
      allocate(current_state%vv_advection(column_size))
    end if
    if (vv_viscosity_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "vv_viscosity_local", .true.)
      allocate(current_state%vv_viscosity(column_size))
    end if
    if (wv_v_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "wv_v_local", .true.)
      allocate(current_state%wv_v(column_size))
    end if
    if (tw_sw_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "tw_sw_local", .true.)
      allocate(current_state%tw_sw(column_size))
    end if
    if (ww_advection_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "ww_advection_local", .true.)
      allocate(current_state%ww_advection(column_size))
    end if
    if (ww_viscosity_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "ww_viscosity_local", .true.)
      allocate(current_state%ww_viscosity(column_size))
    end if
    if (ww_buoyancy_enabled) then
      !call set_published_field_enabled_state(prognostic_budget_fields, "ww_buoyancy_local", .true.)
      allocate(current_state%ww_buoyancy(column_size))
    end if
  end subroutine initialise_prognostic_budget_diagnostics


  !> Initialises the thetal diagnostics
  !! @param current_state The current model state
  subroutine initialise_thetal_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: column_size
    logical :: u_thetal_enabled, us_thetal_enabled, u_thetal_advection_enabled, u_thetal_viscosity_diffusion_enabled, &
         wu_thetal_enabled, v_thetal_enabled, vs_thetal_enabled, v_thetal_advection_enabled, &
         v_thetal_viscosity_diffusion_enabled, wv_thetal_enabled, w_thetal_enabled, ws_thetal_enabled, &
         w_thetal_advection_enabled, w_thetal_viscosity_diffusion_enabled, w_thetal_buoyancy_enabled, ww_thetal_enabled, &
         thetal_thetal_enabled, sthetal_thetal_enabled, thetal_thetal_advection_enabled, &
         thetal_thetal_diffusion_enabled, wthetal_thetal_enabled

    column_size=current_state%local_grid%size(Z_INDEX)

    u_thetal_enabled = current_state%u%active .and. current_state%th%active
    us_thetal_enabled = u_thetal_enabled
    u_thetal_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
                                  .and. u_thetal_enabled
    u_thetal_viscosity_diffusion_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
                                          .and. u_thetal_enabled
    wu_thetal_enabled = current_state%w%active .and. u_thetal_enabled
    some_thetal_diagnostics_enabled = u_thetal_enabled

    if (u_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "u_thetal_local", .true.)
      allocate(current_state%u_thetal(column_size))
    end if
    if (us_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "us_thetal_local", .true.)
      allocate(current_state%us_thetal(column_size))
    end if
    if (u_thetal_advection_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "u_thetal_advection_local", .true.)
      allocate(current_state%u_thetal_advection(column_size))
    end if
    if (u_thetal_viscosity_diffusion_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "u_thetal_viscosity_diffusion_local", .true.)
      allocate(current_state%u_thetal_viscosity_diffusion(column_size))
    end if
    if (wu_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "wu_thetal_local", .true.)
      allocate(current_state%wu_thetal(column_size))
    end if

    v_thetal_enabled = current_state%v%active .and. current_state%th%active
    vs_thetal_enabled = v_thetal_enabled
    v_thetal_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
                                  .and. v_thetal_enabled
    v_thetal_viscosity_diffusion_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
                                          .and. v_thetal_enabled
    wv_thetal_enabled = current_state%w%active .and. v_thetal_enabled
    some_thetal_diagnostics_enabled = some_thetal_diagnostics_enabled .or. v_thetal_enabled

    if (v_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "v_thetal_local", .true.)
      allocate(current_state%v_thetal(column_size))
    end if
    if (vs_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "vs_thetal_local", .true.)
      allocate(current_state%vs_thetal(column_size))
    end if
    if (v_thetal_advection_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "v_thetal_advection_local", .true.)
      allocate(current_state%v_thetal_advection(column_size))
    end if
    if (v_thetal_viscosity_diffusion_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "v_thetal_viscosity_diffusion_local", .true.)
      allocate(current_state%v_thetal_viscosity_diffusion(column_size))
    end if
    if (wv_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "wv_thetal_local", .true.)
      allocate(current_state%wv_thetal(column_size))
    end if

    w_thetal_enabled = current_state%w%active .and. current_state%th%active
    ws_thetal_enabled = w_thetal_enabled
    w_thetal_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
                                .and. w_thetal_enabled
    w_thetal_viscosity_diffusion_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
                                          .and. w_thetal_enabled
    w_thetal_buoyancy_enabled = current_state%th%active .and. current_state%buoyancy_enabled
    ww_thetal_enabled = w_thetal_enabled
    some_thetal_diagnostics_enabled=some_thetal_diagnostics_enabled .or. w_thetal_enabled

    if (w_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "w_thetal_local", .true.)
      allocate(current_state%w_thetal(column_size))
    end if
    if (ws_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "ws_thetal_local", .true.)
      allocate(current_state%ws_thetal(column_size))
    end if
    if (w_thetal_advection_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "w_thetal_advection_local", .true.)
      allocate(current_state%w_thetal_advection(column_size))
    end if
    if (w_thetal_viscosity_diffusion_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "w_thetal_viscosity_diffusion_local", .true.)
      allocate(current_state%w_thetal_viscosity_diffusion(column_size))
    end if
    if (w_thetal_buoyancy_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "w_thetal_buoyancy_local", .true.)
      allocate(current_state%w_thetal_buoyancy(column_size))
    end if
    if (ww_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "ww_thetal_local", .true.)
      allocate(current_state%ww_thetal(column_size))
    end if

    thetal_thetal_enabled = current_state%th%active
    sthetal_thetal_enabled = thetal_thetal_enabled
    thetal_thetal_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
                                      .and. thetal_thetal_enabled
    thetal_thetal_diffusion_enabled = current_state%diffusion_enabled .and. thetal_thetal_enabled
    wthetal_thetal_enabled = current_state%w%active .and. thetal_thetal_enabled
    some_thetal_diagnostics_enabled=some_thetal_diagnostics_enabled .or. thetal_thetal_enabled

    if (thetal_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "thetal_thetal_local", .true.)
      allocate(current_state%thetal_thetal(column_size))
    end if
    if (sthetal_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "sthetal_thetal_local", .true.)
      allocate(current_state%sthetal_thetal(column_size))
    end if
    if (thetal_thetal_advection_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "thetal_thetal_advection_local", .true.)
      allocate(current_state%thetal_thetal_advection(column_size))
    end if
    if (thetal_thetal_diffusion_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "thetal_thetal_diffusion_local", .true.)
      allocate(current_state%thetal_thetal_diffusion(column_size))
    end if
    if (wthetal_thetal_enabled) then
      !call set_published_field_enabled_state(thetal_fields, "wthetal_thetal_local", .true.)
      allocate(current_state%wthetal_thetal(column_size))
    end if
  end subroutine initialise_thetal_diagnostics


  !> Initialises the scalar diagnostics
  !! @param current_state The current model state
  subroutine initialise_scalar_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(DEFAULT_PRECISION) :: realnum
    integer :: iter

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "wmfcrit") then
          read(current_state%options_database_string(iter,2),*) realnum
          current_state%wmfcrit = realnum
      end if
    end do
  end subroutine initialise_scalar_diagnostics

  !> Initialises the TKE diagnostics
  !! @param current_state The current model state
  subroutine initialise_tke_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: column_size
    logical :: shprd_enabled, w_p_enabled, w_ke_enabled, buoy_enabled, tend_enabled

    column_size=current_state%local_grid%size(Z_INDEX)

    shprd_enabled = current_state%u%active .and. current_state%v%active
    w_ke_enabled = current_state%w%active
    buoy_enabled = current_state%w%active
    w_p_enabled = current_state%w%active
    tend_enabled = current_state%w%active
    some_tke_diagnostics_enabled = w_p_enabled .or. shprd_enabled .or. w_ke_enabled .or. buoy_enabled .or. tend_enabled

    if (w_p_enabled) then
      !call set_published_field_enabled_state(tke_fields, "resolved_pressure_transport_local", .true.)
      allocate(current_state%w_p(column_size))
    end if
    if (tend_enabled) then
      !call set_published_field_enabled_state(tke_fields, "tke_tendency_local", .true.)
      allocate(current_state%tke_tend(column_size))
    end if
    if (shprd_enabled) then
      !call set_published_field_enabled_state(tke_fields, "resolved_shear_production_local", .true.)
      allocate(current_state%shprd(column_size))
    end if
    if (w_ke_enabled) then
      !call set_published_field_enabled_state(tke_fields, "resolved_turbulent_transport_local", .true.)
      allocate(current_state%w_ke(column_size))
    end if
    if (buoy_enabled) then
      !call set_published_field_enabled_state(tke_fields, "resolved_buoyant_production_local", .true.)
      allocate(current_state%sub_buoy(column_size))
    end if
  end subroutine initialise_tke_diagnostics


  !> Initialises the mse diagnostics. For now we are assuming mse is the same as theta, which needs updating with moisture
  !! information
  !! @param current_state The current model state
!   subroutine initialise_mse_diagnostics(current_state)
!     type(model_state_type), target, intent(inout) :: current_state
!
!     integer :: column_size
!     logical :: u_mse_enabled, us_mse_enabled, u_mse_advection_enabled, u_mse_viscosity_diffusion_enabled, &
!          wu_mse_enabled, v_mse_enabled, vs_mse_enabled, v_mse_advection_enabled, &
!          v_mse_viscosity_diffusion_enabled, wv_mse_enabled, w_mse_enabled, ws_mse_enabled, &
!          w_mse_advection_enabled, w_mse_viscosity_diffusion_enabled, w_mse_buoyancy_enabled, ww_mse_enabled, &
!          mse_mse_enabled, smse_mse_enabled, mse_mse_advection_enabled, &
!          mse_mse_diffusion_enabled, wmse_mse_enabled
!
!     column_size=current_state%local_grid%size(Z_INDEX)
!
!     u_mse_enabled = current_state%u%active .and. current_state%th%active
!     us_mse_enabled = u_mse_enabled
!     u_mse_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
!                               .and. u_mse_enabled
!     u_mse_viscosity_diffusion_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
!                                       .and. u_mse_enabled
!     wu_mse_enabled = current_state%w%active .and. u_mse_enabled
!     some_mse_diagnostics_enabled=u_mse_enabled
!
!     if (u_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "u_mse_local", .true.)
!       allocate(current_state%u_mse(column_size))
!     end if
!     if (us_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "us_mse_local", .true.)
!       allocate(current_state%us_mse(column_size))
!     end if
!     if (u_mse_advection_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "u_mse_advection_local", .true.)
!       allocate(current_state%u_mse_advection(column_size))
!     end if
!     if (u_mse_viscosity_diffusion_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "u_mse_viscosity_diffusion_local", .true.)
!       allocate(current_state%u_mse_viscosity_diffusion(column_size))
!     end if
!     if (wu_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "wu_mse_local", .true.)
!       allocate(current_state%wu_mse(column_size))
!     end if
!
!     v_mse_enabled = current_state%v%active .and. current_state%th%active
!     vs_mse_enabled = v_mse_enabled
!     v_mse_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
!                               .and. v_mse_enabled
!     v_mse_viscosity_diffusion_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
!                                       .and. v_mse_enabled
!     wv_mse_enabled=current_state%w%active .and. v_mse_enabled
!     some_mse_diagnostics_enabled=some_mse_diagnostics_enabled .or. v_mse_enabled
!
!     if (v_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "v_mse_local", .true.)
!       allocate(current_state%v_mse(column_size))
!     end if
!     if (vs_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "vs_mse_local", .true.)
!       allocate(current_state%vs_mse(column_size))
!     end if
!     if (v_mse_advection_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "v_mse_advection_local", .true.)
!       allocate(current_state%v_mse_advection(column_size))
!     end if
!     if (v_mse_viscosity_diffusion_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "v_mse_viscosity_diffusion_local", .true.)
!       allocate(current_state%v_mse_viscosity_diffusion(column_size))
!     end if
!     if (wv_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "wv_mse_local", .true.)
!       allocate(current_state%wv_mse(column_size))
!     end if
!
!     w_mse_enabled = current_state%w%active .and. current_state%th%active
!     ws_mse_enabled = w_mse_enabled
!     w_mse_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
!                               .and. w_mse_enabled
!     w_mse_viscosity_diffusion_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
!                                       .and. w_mse_enabled
!     w_mse_buoyancy_enabled = current_state%th%active .and. current_state%buoyancy_enabled
!     ww_mse_enabled = w_mse_enabled
!     some_mse_diagnostics_enabled=some_mse_diagnostics_enabled .or. w_mse_enabled
!
!     if (w_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "w_mse_local", .true.)
!       allocate(current_state%w_mse(column_size))
!     end if
!     if (ws_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "ws_mse_local", .true.)
!       allocate(current_state%ws_mse(column_size))
!     end if
!     if (w_mse_advection_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "w_mse_advection_local", .true.)
!       allocate(current_state%w_mse_advection(column_size))
!     end if
!     if (w_mse_viscosity_diffusion_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "w_mse_viscosity_diffusion_local", .true.)
!       allocate(current_state%w_mse_viscosity_diffusion(column_size))
!     end if
!     if (w_mse_buoyancy_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "w_mse_buoyancy_local", .true.)
!       allocate(current_state%w_mse_buoyancy(column_size))
!     end if
!     if (ww_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "ww_mse_local", .true.)
!       allocate(current_state%ww_mse(column_size))
!     end if
!
!     mse_mse_enabled = current_state%th%active
!     smse_mse_enabled = mse_mse_enabled
!     mse_mse_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
!                                 .and. mse_mse_enabled
!     mse_mse_diffusion_enabled = current_state%diffusion_enabled .and. mse_mse_enabled
!     wmse_mse_enabled = current_state%w%active .and. mse_mse_enabled
!     some_mse_diagnostics_enabled=some_mse_diagnostics_enabled .or. mse_mse_enabled
!
!     if (mse_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "mse_mse_local", .true.)
!       allocate(current_state%mse_mse(column_size))
!     end if
!     if (smse_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "smse_mse_local", .true.)
!       allocate(current_state%smse_mse(column_size))
!     end if
!     if (mse_mse_advection_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "mse_mse_advection_local", .true.)
!       allocate(current_state%mse_mse_advection(column_size))
!     end if
!     if (mse_mse_diffusion_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "mse_mse_diffusion_local", .true.)
!       allocate(current_state%mse_mse_diffusion(column_size))
!     end if
!     if (wmse_mse_enabled) then
!       !call set_published_field_enabled_state(mse_fields, "wmse_mse_local", .true.)
!       allocate(current_state%wmse_mse(column_size))
!     end if
!   end subroutine initialise_mse_diagnostics

  !> Initialises the qt diagnostics. For now we are assuming qt is the same as theta, which needs updating with moisture
  !! information
  !! @param current_state The current model state
!   subroutine initialise_qt_diagnostics(current_state)
!     type(model_state_type), target, intent(inout) :: current_state
!
!     integer :: column_size
!     logical :: us_qt_enabled, u_qt_advection_enabled, u_qt_viscosity_diffusion_enabled, &
!          wu_qt_enabled, vs_qt_enabled, v_qt_advection_enabled, &
!          v_qt_viscosity_diffusion_enabled, wv_qt_enabled, w_qt_enabled, ws_qt_enabled, &
!          w_qt_advection_enabled, w_qt_viscosity_diffusion_enabled, w_qt_buoyancy_enabled, ww_qt_enabled, &
!          qt_qt_enabled, sqt_qt_enabled, qt_qt_advection_enabled, &
!          qt_qt_diffusion_enabled, wqt_qt_enabled
!
!     column_size=current_state%local_grid%size(Z_INDEX)
!
!     us_qt_enabled = current_state%u%active .and. current_state%th%active
!     u_qt_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
!                               .and. us_qt_enabled
!     u_qt_viscosity_diffusion_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
!                                       .and. us_qt_enabled
!     wu_qt_enabled = current_state%w%active .and. us_qt_enabled
!     some_qt_diagnostics_enabled = us_qt_enabled
!
!     if (us_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "us_qt_local", .true.)
!       allocate(current_state%us_qt(column_size))
!     end if
!     if (u_qt_advection_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "u_qt_advection_local", .true.)
!       allocate(current_state%u_qt_advection(column_size))
!     end if
!     if (u_qt_viscosity_diffusion_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "u_qt_viscosity_diffusion_local", .true.)
!       allocate(current_state%u_qt_viscosity_diffusion(column_size))
!     end if
!     if (wu_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "wu_qt_local", .true.)
!       allocate(current_state%wu_qt(column_size))
!     end if
!
!     vs_qt_enabled = current_state%v%active .and. current_state%th%active
!     v_qt_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
!                               .and. vs_qt_enabled
!     v_qt_viscosity_diffusion_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
!                                       .and. vs_qt_enabled
!     wv_qt_enabled = current_state%w%active .and. vs_qt_enabled
!     some_qt_diagnostics_enabled = some_qt_diagnostics_enabled .or. vs_qt_enabled
!
!     if (vs_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "vs_qt_local", .true.)
!       allocate(current_state%vs_qt(column_size))
!     end if
!     if (v_qt_advection_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "v_qt_advection_local", .true.)
!       allocate(current_state%v_qt_advection(column_size))
!     end if
!     if (v_qt_viscosity_diffusion_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "v_qt_viscosity_diffusion_local", .true.)
!       allocate(current_state%v_qt_viscosity_diffusion(column_size))
!     end if
!     if (wv_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "wv_qt_local", .true.)
!       allocate(current_state%wv_qt(column_size))
!     end if
!
!     w_qt_enabled = current_state%w%active .and. current_state%th%active
!     ws_qt_enabled = w_qt_enabled
!     w_qt_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
!                               .and. w_qt_enabled
!     w_qt_viscosity_diffusion_enabled = current_state%diffusion_enabled .and. current_state%viscosity_enabled &
!                                       .and. w_qt_enabled
!     w_qt_buoyancy_enabled = current_state%th%active .and. current_state%buoyancy_enabled
!     ww_qt_enabled = w_qt_enabled
!     some_qt_diagnostics_enabled = some_qt_diagnostics_enabled .or. w_qt_enabled
!
!     if (w_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "w_qt_local", .true.)
!       allocate(current_state%w_qt(column_size))
!     end if
!     if (ws_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "ws_qt_local", .true.)
!       allocate(current_state%ws_qt(column_size))
!     end if
!     if (w_qt_advection_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "w_qt_advection_local", .true.)
!       allocate(current_state%w_qt_advection(column_size))
!     end if
!     if (w_qt_viscosity_diffusion_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "w_qt_viscosity_diffusion_local", .true.)
!       allocate(current_state%w_qt_viscosity_diffusion(column_size))
!     end if
!     if (w_qt_buoyancy_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "w_qt_buoyancy_local", .true.)
!       allocate(current_state%w_qt_buoyancy(column_size))
!     end if
!     if (ww_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "ww_qt_local", .true.)
!       allocate(current_state%ww_qt(column_size))
!     end if
!
!     qt_qt_enabled = current_state%th%active
!     sqt_qt_enabled = qt_qt_enabled
!     qt_qt_advection_enabled = (current_state%pw_advection_enabled .or. current_state%th_advection_enabled) &
!                               .and. qt_qt_enabled
!     qt_qt_diffusion_enabled = current_state%diffusion_enabled .and. qt_qt_enabled
!     wqt_qt_enabled = current_state%w%active .and. qt_qt_enabled
!     some_qt_diagnostics_enabled=some_qt_diagnostics_enabled .or. qt_qt_enabled
!
!     if (qt_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "qt_qt_local", .true.)
!       allocate(current_state%qt_qt(column_size))
!     end if
!     if (sqt_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "sqt_qt_local", .true.)
!       allocate(current_state%sqt_qt(column_size))
!     end if
!     if (qt_qt_advection_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "qt_qt_advection_local", .true.)
!       allocate(current_state%qt_qt_advection(column_size))
!     end if
!     if (qt_qt_diffusion_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "qt_qt_diffusion_local", .true.)
!       allocate(current_state%qt_qt_diffusion(column_size))
!     end if
!     if (wqt_qt_enabled) then
!       !call set_published_field_enabled_state(qt_fields, "wqt_qt_local", .true.)
!       allocate(current_state%wqt_qt(column_size))
!     end if
!   end subroutine initialise_qt_diagnostics




    !> Populates the published field names in the appropriate map
!   subroutine populate_field_names()
!     call set_published_field_enabled_state(heat_flux_fields, "heat_flux_transport_local", .false.)
!     call set_published_field_enabled_state(heat_flux_fields, "heat_flux_gradient_local", .false.)
!     call set_published_field_enabled_state(heat_flux_fields, "heat_flux_dissipation_local", .false.)
!     call set_published_field_enabled_state(heat_flux_fields, "heat_flux_buoyancy_local", .false.)
!     call set_published_field_enabled_state(heat_flux_fields, "heat_flux_tendency_local", .false.)
!
!     call set_published_field_enabled_state(q_flux_fields, "q_flux_transport_local", .false.)
!     call set_published_field_enabled_state(q_flux_fields, "q_flux_gradient_local", .false.)
!     call set_published_field_enabled_state(q_flux_fields, "q_flux_dissipation_local", .false.)
!     call set_published_field_enabled_state(q_flux_fields, "q_flux_buoyancy_local", .false.)
!     call set_published_field_enabled_state(q_flux_fields, "q_flux_tendency_local", .false.)
!
!     call set_published_field_enabled_state(tke_fields, "resolved_shear_production_local", .false.)
!     call set_published_field_enabled_state(tke_fields, "resolved_buoyant_production_local", .false.)
!     call set_published_field_enabled_state(tke_fields, "resolved_turbulent_transport_local", .false.)
!     call set_published_field_enabled_state(tke_fields, "resolved_pressure_transport_local", .false.)
!     call set_published_field_enabled_state(tke_fields, "tke_tendency_local", .false.)
!
!     call set_published_field_enabled_state(uw_vw_fields, "uw_advection_local", .false.)
!     call set_published_field_enabled_state(uw_vw_fields, "vw_advection_local", .false.)
!     call set_published_field_enabled_state(uw_vw_fields, "uw_viscosity_local", .false.)
!     call set_published_field_enabled_state(uw_vw_fields, "vw_viscosity_local", .false.)
!     call set_published_field_enabled_state(uw_vw_fields, "uw_buoyancy_local", .false.)
!     call set_published_field_enabled_state(uw_vw_fields, "vw_buoyancy_local", .false.)
!     call set_published_field_enabled_state(uw_vw_fields, "uw_tendency_local", .false.)
!     call set_published_field_enabled_state(uw_vw_fields, "vw_tendency_local", .false.)
!     call set_published_field_enabled_state(uw_vw_fields, "uw_w_local", .false.)
!     call set_published_field_enabled_state(uw_vw_fields, "vw_w_local", .false.)
!
!     call set_published_field_enabled_state(prognostic_budget_fields, "tu_su_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "uu_advection_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "uu_viscosity_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "wu_u_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "tv_sv_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "vv_advection_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "vv_viscosity_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "wv_v_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "tw_sw_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "ww_advection_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "ww_viscosity_local", .false.)
!     call set_published_field_enabled_state(prognostic_budget_fields, "ww_buoyancy_local", .false.)
!
!     call set_published_field_enabled_state(thetal_fields, "u_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "us_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "u_thetal_advection_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "u_thetal_viscosity_diffusion_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "wu_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "v_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "vs_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "v_thetal_advection_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "v_thetal_viscosity_diffusion_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "wv_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "w_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "ws_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "w_thetal_advection_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "w_thetal_viscosity_diffusion_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "ww_thetal_buoyancy_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "ww_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "thetal_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "sthetal_thetal_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "thetal_thetal_advection_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "thetal_thetal_diffusion_local", .false.)
!     call set_published_field_enabled_state(thetal_fields, "wthetal_thetal_local", .false.)
!
!     call set_published_field_enabled_state(mse_fields, "u_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "us_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "u_mse_advection_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "u_mse_viscosity_diffusion_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "wu_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "v_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "vs_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "v_mse_advection_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "v_mse_viscosity_diffusion_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "wv_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "w_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "ws_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "w_mse_advection_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "w_mse_viscosity_diffusion_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "ww_mse_buoyancy_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "ww_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "mse_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "smse_mse_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "mse_mse_advection_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "mse_mse_diffusion_local", .false.)
!     call set_published_field_enabled_state(mse_fields, "wmse_mse_local", .false.)
!
!     call set_published_field_enabled_state(qt_fields, "us_qt_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "u_qt_advection_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "u_qt_viscosity_diffusion_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "wu_qt_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "vs_qt_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "v_qt_advection_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "v_qt_viscosity_diffusion_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "wv_qt_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "w_qt_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "ws_qt_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "w_qt_advection_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "w_qt_viscosity_diffusion_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "ww_qt_buoyancy_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "ww_qt_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "qt_qt_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "sqt_qt_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "qt_qt_advection_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "qt_qt_diffusion_local", .false.)
!     call set_published_field_enabled_state(qt_fields, "wqt_qt_local", .false.)
!
!     call set_published_field_enabled_state(scalar_fields, "mflux_local", .true.)
!   end subroutine populate_field_names

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
!   subroutine field_information_retrieval_callback(current_state, name, field_information)
!     type(model_state_type), target, intent(inout) :: current_state
!     character(len=*), intent(in) :: name
!     type(component_field_information_type), intent(out) :: field_information
!
!     ! Field description is the same regardless of the specific field being retrieved
!     field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
!     field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
!
!     if (is_field_heat_flux(name) .or. is_field_uw_vw(name) .or. is_field_prognostic_budget(name) &
!         .or. is_field_thetal(name) .or. is_field_mse(name) .or. is_field_qt(name) .or. is_field_scalar(name) &
!         .or. is_field_tke_flux(name)) then
!       field_information%number_dimensions=1
!       field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
!       if (is_field_heat_flux(name)) then
!         field_information%enabled=get_published_field_enabled_state(heat_flux_fields, name)
!       else if (is_field_uw_vw(name)) then
!         field_information%enabled=get_published_field_enabled_state(uw_vw_fields, name)
!       else if (is_field_prognostic_budget(name)) then
!         field_information%enabled=get_published_field_enabled_state(prognostic_budget_fields, name)
!       else if (is_field_thetal(name)) then
!         field_information%enabled=get_published_field_enabled_state(thetal_fields, name)
!       else if (is_field_tke_flux(name)) then
!         field_information%enabled=get_published_field_enabled_state(tke_fields, name)
!       else if (is_field_mse(name)) then
!         field_information%enabled=get_published_field_enabled_state(mse_fields, name)
!       else if (is_field_qt(name)) then
!         field_information%enabled=get_published_field_enabled_state(qt_fields, name)
!       else if (is_field_scalar(name)) then
!         field_information%enabled=get_published_field_enabled_state(scalar_fields, name)
!       end if
!     else if (is_field_q_flux(name)) then
!       field_information%number_dimensions=2
!       field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
!       field_information%dimension_sizes(2)=current_state%number_q_fields
!       field_information%enabled=get_published_field_enabled_state(q_flux_fields, name)
!     end if
!   end subroutine field_information_retrieval_callback

  !> Timestep call back, this will deduce the diagnostics for the current (non halo) column
  !! @param current_state Current model state
  subroutine timestep_callback_flux_budget(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep

    if (calculate_diagnostics) then
      if (current_state%first_timestep_column) then
        if (current_state%modulo_number_1d .ne. 0) then
          if (some_theta_flux_diagnostics_enabled) call clear_theta_fluxes(current_state)
          if (some_q_flux_diagnostics_enabled) call clear_q_fluxes(current_state)
          if (some_uw_vw_diagnostics_enabled) call clear_uw_vw(current_state)
          if (some_prognostic_budget_diagnostics_enabled) call clear_prognostic_budgets(current_state)
          if (some_thetal_diagnostics_enabled) call clear_thetal(current_state)
          if (some_tke_diagnostics_enabled) call clear_tke(current_state)
        end if
        if (current_state%modulo_number_0d .ne. 0) then
          call clear_scalars(current_state)
        end if
        !if (some_mse_diagnostics_enabled) call clear_mse()
        !if (some_qt_diagnostics_enabled) call clear_qt()
      end if
      if (.not. current_state%halo_column) then
        if (current_state%modulo_number_1d .ne. 0) then
          if (some_theta_flux_diagnostics_enabled) call compute_theta_flux_for_column(current_state)
          if (some_q_flux_diagnostics_enabled) call compute_q_flux_for_column(current_state)
          if (some_uw_vw_diagnostics_enabled) call compute_uw_vw_for_column(current_state)
          if (some_prognostic_budget_diagnostics_enabled) call compute_prognostic_budgets_for_column(current_state)
          if (some_thetal_diagnostics_enabled) call compute_thetal_for_column(current_state)
          if (some_tke_diagnostics_enabled) call compute_tke_for_column(current_state)
        end if
        if (current_state%modulo_number_0d .ne. 0) then
          call compute_scalars_for_column(current_state)
        end if
        !if (some_mse_diagnostics_enabled) call compute_mse_for_column(current_state)
        !if (some_qt_diagnostics_enabled) call compute_qt_for_column(current_state)
      end if
    end if
  end subroutine timestep_callback_flux_budget






  !> Clears the heat flux diagnostics at the start of a timestep
  subroutine clear_theta_fluxes(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    if (allocated(current_state%th_flux_values)) current_state%th_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%th_tendency)) current_state%th_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%th_gradient)) current_state%th_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%th_diff)) current_state%th_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%th_buoyancy)) current_state%th_buoyancy = 0.0_DEFAULT_PRECISION
  end subroutine clear_theta_fluxes

  !> Clears the Q flux diagnostics, called at the start of a timestep
  subroutine clear_q_fluxes(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    if (allocated(current_state%qv_flux_values)) current_state%qv_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ql_flux_values)) current_state%ql_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qr_flux_values)) current_state%qr_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qi_flux_values)) current_state%qi_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qs_flux_values)) current_state%qs_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qg_flux_values)) current_state%qg_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAitkenSolMass_flux_values)) current_state%qAitkenSolMass_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumSolMass_flux_values)) current_state%qAccumSolMass_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumInsolMass_flux_values)) current_state%qAccumInsolMass_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseSolMass_flux_values)) current_state%qCoarseSolMass_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseDustMass_flux_values)) current_state%qCoarseDustMass_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nl_flux_values)) current_state%nl_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nr_flux_values)) current_state%nr_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ni_flux_values)) current_state%ni_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ns_flux_values)) current_state%ns_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ng_flux_values)) current_state%ng_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAitkenSolNumber_flux_values)) current_state%nAitkenSolNumber_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumSolNumber_flux_values)) current_state%nAccumSolNumber_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumInsolNumber_flux_values)) current_state%nAccumInsolNumber_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseSolNumber_flux_values)) current_state%nCoarseSolNumber_flux_values = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseDustnumber_flux_values)) current_state%nCoarseDustnumber_flux_values = 0.0_DEFAULT_PRECISION

    if (allocated(current_state%qv_gradient)) current_state%qv_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ql_gradient)) current_state%ql_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qr_gradient)) current_state%qr_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qi_gradient)) current_state%qi_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qs_gradient)) current_state%qs_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qg_gradient)) current_state%qg_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAitkenSolMass_gradient)) current_state%qAitkenSolMass_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumSolMass_gradient)) current_state%qAccumSolMass_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumInsolMass_gradient)) current_state%qAccumInsolMass_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseSolMass_gradient)) current_state%qCoarseSolMass_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseDustMass_gradient)) current_state%qCoarseDustMass_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nl_gradient)) current_state%nl_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nr_gradient)) current_state%nr_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ni_gradient)) current_state%ni_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ns_gradient)) current_state%ns_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ng_gradient)) current_state%ng_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAitkenSolNumber_gradient)) current_state%nAitkenSolNumber_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumSolNumber_gradient)) current_state%nAccumSolNumber_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumInsolNumber_gradient)) current_state%nAccumInsolNumber_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseSolNumber_gradient)) current_state%nCoarseSolNumber_gradient = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseDustnumber_gradient)) current_state%nCoarseDustnumber_gradient = 0.0_DEFAULT_PRECISION

    if (allocated(current_state%qv_diff)) current_state%qv_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ql_diff)) current_state%ql_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qr_diff)) current_state%qr_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qi_diff)) current_state%qi_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qs_diff)) current_state%qs_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qg_diff)) current_state%qg_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAitkenSolMass_diff)) current_state%qAitkenSolMass_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumSolMass_diff)) current_state%qAccumSolMass_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumInsolMass_diff)) current_state%qAccumInsolMass_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseSolMass_diff)) current_state%qCoarseSolMass_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseDustMass_diff)) current_state%qCoarseDustMass_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nl_diff)) current_state%nl_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nr_diff)) current_state%nr_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ni_diff)) current_state%ni_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ns_diff)) current_state%ns_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ng_diff)) current_state%ng_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAitkenSolNumber_diff)) current_state%nAitkenSolNumber_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumSolNumber_diff)) current_state%nAccumSolNumber_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumInsolNumber_diff)) current_state%nAccumInsolNumber_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseSolNumber_diff)) current_state%nCoarseSolNumber_diff = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseDustnumber_diff)) current_state%nCoarseDustnumber_diff = 0.0_DEFAULT_PRECISION

    if (allocated(current_state%qv_buoyancy)) current_state%qv_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ql_buoyancy)) current_state%ql_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qr_buoyancy)) current_state%qr_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qi_buoyancy)) current_state%qi_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qs_buoyancy)) current_state%qs_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qg_buoyancy)) current_state%qg_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAitkenSolMass_buoyancy)) current_state%qAitkenSolMass_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumSolMass_buoyancy)) current_state%qAccumSolMass_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumInsolMass_buoyancy)) current_state%qAccumInsolMass_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseSolMass_buoyancy)) current_state%qCoarseSolMass_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseDustMass_buoyancy)) current_state%qCoarseDustMass_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nl_buoyancy)) current_state%nl_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nr_buoyancy)) current_state%nr_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ni_buoyancy)) current_state%ni_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ns_buoyancy)) current_state%ns_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ng_buoyancy)) current_state%ng_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAitkenSolNumber_buoyancy)) current_state%nAitkenSolNumber_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumSolNumber_buoyancy)) current_state%nAccumSolNumber_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumInsolNumber_buoyancy)) current_state%nAccumInsolNumber_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseSolNumber_buoyancy)) current_state%nCoarseSolNumber_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseDustnumber_buoyancy)) current_state%nCoarseDustnumber_buoyancy = 0.0_DEFAULT_PRECISION

    if (allocated(current_state%qv_tendency)) current_state%qv_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ql_tendency)) current_state%ql_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qr_tendency)) current_state%qr_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qi_tendency)) current_state%qi_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qs_tendency)) current_state%qs_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qg_tendency)) current_state%qg_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAitkenSolMass_tendency)) current_state%qAitkenSolMass_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumSolMass_tendency)) current_state%qAccumSolMass_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qAccumInsolMass_tendency)) current_state%qAccumInsolMass_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseSolMass_tendency)) current_state%qCoarseSolMass_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%qCoarseDustMass_tendency)) current_state%qCoarseDustMass_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nl_tendency)) current_state%nl_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nr_tendency)) current_state%nr_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ni_tendency)) current_state%ni_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ns_tendency)) current_state%ns_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ng_tendency)) current_state%ng_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAitkenSolNumber_tendency)) current_state%nAitkenSolNumber_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumSolNumber_tendency)) current_state%nAccumSolNumber_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nAccumInsolNumber_tendency)) current_state%nAccumInsolNumber_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseSolNumber_tendency)) current_state%nCoarseSolNumber_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%nCoarseDustnumber_tendency)) current_state%nCoarseDustnumber_tendency = 0.0_DEFAULT_PRECISION
  end subroutine clear_q_fluxes

  !> Clears the uw uv diagnostics
  subroutine clear_uw_vw(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    if (allocated(current_state%uw_advection)) current_state%uw_advection = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%vw_advection)) current_state%vw_advection = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%uw_viscosity)) current_state%uw_viscosity = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%vw_viscosity)) current_state%vw_viscosity = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%uw_buoyancy)) current_state%uw_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%vw_buoyancy)) current_state%vw_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%uw_tendency)) current_state%uw_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%vw_tendency)) current_state%vw_tendency = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%uw_w)) current_state%uw_w = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%vw_w)) current_state%vw_w = 0.0_DEFAULT_PRECISION
  end subroutine clear_uw_vw

  !> Clears the prognostic (uu, vv, ww) budgets
  subroutine clear_prognostic_budgets(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    if (allocated(current_state%tu_su)) current_state%tu_su = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%uu_advection)) current_state%uu_advection = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%uu_viscosity)) current_state%uu_viscosity = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%wu_u)) current_state%wu_u = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%tv_sv)) current_state%tv_sv = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%vv_advection)) current_state%vv_advection = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%vv_viscosity)) current_state%vv_viscosity = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%wv_v)) current_state%wv_v = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%tw_sw)) current_state%tw_sw = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ww_advection)) current_state%ww_advection = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ww_viscosity)) current_state%ww_viscosity = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ww_buoyancy)) current_state%ww_buoyancy = 0.0_DEFAULT_PRECISION
  end subroutine clear_prognostic_budgets

  !> Clears the thetal diagnostics
  subroutine clear_thetal(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    if (allocated(current_state%u_thetal)) current_state%u_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%us_thetal)) current_state%us_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%u_thetal_advection)) current_state%u_thetal_advection = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%u_thetal_viscosity_diffusion)) current_state%u_thetal_viscosity_diffusion = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%wu_thetal)) current_state%wu_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%v_thetal)) current_state%v_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%vs_thetal)) current_state%vs_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%v_thetal_advection)) current_state%v_thetal_advection = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%v_thetal_viscosity_diffusion)) current_state%v_thetal_viscosity_diffusion = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%wv_thetal)) current_state%wv_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%w_thetal)) current_state%w_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ws_thetal)) current_state%ws_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%w_thetal_advection)) current_state%w_thetal_advection=0.0_DEFAULT_PRECISION
    if (allocated(current_state%w_thetal_viscosity_diffusion)) current_state%w_thetal_viscosity_diffusion = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%w_thetal_buoyancy)) current_state%w_thetal_buoyancy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%ww_thetal)) current_state%ww_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%thetal_thetal)) current_state%thetal_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%sthetal_thetal)) current_state%sthetal_thetal = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%thetal_thetal_advection)) current_state%thetal_thetal_advection = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%thetal_thetal_diffusion)) current_state%thetal_thetal_diffusion = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%wthetal_thetal)) current_state%wthetal_thetal = 0.0_DEFAULT_PRECISION
  end subroutine clear_thetal


  !> Clears the TKE diagnostics
  subroutine clear_tke(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    if (allocated(current_state%shprd)) current_state%shprd = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%w_ke)) current_state%w_ke = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%w_p)) current_state%w_p = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%sub_buoy)) current_state%sub_buoy = 0.0_DEFAULT_PRECISION
    if (allocated(current_state%tke_tend)) current_state%tke_tend = 0.0_DEFAULT_PRECISION
  end subroutine clear_tke

  !> Clears the scalar diagnostics
  subroutine clear_scalars(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    current_state%mflux = 0.0_DEFAULT_PRECISION
  end subroutine clear_scalars


!     !> Clears the mse diagnostics
!   subroutine clear_mse()
!     if (allocated(current_state%u_mse)) current_state%u_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%us_mse)) current_state%us_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%u_mse_advection)) current_state%u_mse_advection = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%u_mse_viscosity_diffusion)) current_state%u_mse_viscosity_diffusion = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%wu_mse)) current_state%wu_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%v_mse)) current_state%v_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%vs_mse)) current_state%vs_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%v_mse_advection)) current_state%v_mse_advection = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%v_mse_viscosity_diffusion)) current_state%v_mse_viscosity_diffusion = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%wv_mse)) current_state%wv_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%w_mse)) current_state%w_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%ws_mse)) current_state%ws_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%w_mse_advection)) current_state%w_mse_advection = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%w_mse_viscosity_diffusion)) current_state%w_mse_viscosity_diffusion = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%w_mse_buoyancy))current_state%w_mse_buoyancy = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%ww_mse)) current_state%ww_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%mse_mse)) current_state%mse_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%smse_mse)) current_state%smse_mse = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%mse_mse_advection)) current_state%mse_mse_advection = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%mse_mse_diffusion)) current_state%mse_mse_diffusion = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%wmse_mse)) current_state%wmse_mse = 0.0_DEFAULT_PRECISION
!   end subroutine clear_mse
!
!   !> Clears the qt diagnostics
!   subroutine clear_qt()
!     if (allocated(current_state%us_qt)) current_state%us_qt = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%u_qt_advection)) current_state%u_qt_advection = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%u_qt_viscosity_diffusion)) current_state%u_qt_viscosity_diffusion = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%wu_qt)) current_state%wu_qt = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%vs_qt)) current_state%vs_qt = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%v_qt_advection)) current_state%v_qt_advection = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%v_qt_viscosity_diffusion)) current_state%v_qt_viscosity_diffusion = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%wv_qt)) current_state%wv_qt = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%w_qt)) current_state%w_qt = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%ws_qt)) current_state%ws_qt= 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%w_qt_advection)) current_state%w_qt_advection = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%w_qt_viscosity_diffusion)) current_state%w_qt_viscosity_diffusion = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%w_qt_buoyancy)) current_state%w_qt_buoyancy = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%ww_qt)) current_state%ww_qt = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%qt_qt)) current_state%qt_qt = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%sqt_qt)) current_state%sqt_qt = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%qt_qt_advection)) current_state%qt_qt_advection = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%qt_qt_diffusion)) current_state%qt_qt_diffusion = 0.0_DEFAULT_PRECISION
!     if (allocated(current_state%wqt_qt)) current_state%wqt_qt = 0.0_DEFAULT_PRECISION
!   end subroutine clear_qt






  !> Computes the heat flux diagnostics for a specific column
  !! @param current_state Current model state
  subroutine compute_theta_flux_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: thpr

    do k=1, current_state%local_grid%size(Z_INDEX)
      thpr(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
        thpr(k)=thpr(k)-current_state%global_grid%configuration%vertical%olthbar(k)
      end if
    end do
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(current_state%th_flux_values)) current_state%th_flux_values(k) = current_state%th_flux_values(k)+&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
           current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(thpr(k)+thpr(k+1))
      if (allocated(current_state%th_tendency)) current_state%th_tendency(k) = current_state%th_tendency(k)+0.5*&
           (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
           current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
           thpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
           current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      if (allocated(current_state%th_gradient)) then
        current_state%th_gradient(k) = current_state%th_gradient(k)+thpr(k)*0.5*(current_state%w_advection(k)+&
             current_state%w_advection(k-1))+0.5*&
             (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
             current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
             current_state%th_advection(k)
      end if
      if (allocated(current_state%th_diff)) then
        current_state%th_diff(k) = current_state%th_diff(k)+thpr(k)*0.5*(current_state%w_viscosity(k)+&
             current_state%w_viscosity(k-1))+&
             0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
             current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
             current_state%th_diffusion(k)
      end if
      if (allocated(current_state%th_buoyancy)) then
        current_state%th_buoyancy(k) = current_state%th_buoyancy(k)+thpr(k)*0.5*(current_state%w_buoyancy(k)+&
             current_state%w_buoyancy(k-1))
      end if
    end do
  end subroutine compute_theta_flux_for_column

  !> Computes the Q flux diagnostics for a specific column
  !! @param current_state Current model state
  subroutine compute_q_flux_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, n
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: qvpr, qlpr, &
      qrpr, qipr, qspr, qgpr, qAiktpr, qAccSolpr, qAccInsolpr, qCoaSolpr, qCoaDuspr, nlpr, nrpr, &
      nipr, nspr, ngpr, nAiktpr, nAccSolpr, nAccInsolpr, nCoaSolpr, nCoaDuspr

    do k=1, current_state%local_grid%size(Z_INDEX)
      qvpr(k)=current_state%qv%data(k,current_state%column_local_y,current_state%column_local_x)
      qlpr(k)=current_state%ql%data(k,current_state%column_local_y,current_state%column_local_x)
      qrpr(k)=current_state%qr%data(k,current_state%column_local_y,current_state%column_local_x)
      qipr(k)=current_state%qi%data(k,current_state%column_local_y,current_state%column_local_x)
      qspr(k)=current_state%qs%data(k,current_state%column_local_y,current_state%column_local_x)
      qgpr(k)=current_state%qg%data(k,current_state%column_local_y,current_state%column_local_x)
      qAiktpr(k)=current_state%qAitkenSolMass%data(k,current_state%column_local_y,current_state%column_local_x)
      qAccSolpr(k)=current_state%qAccumSolMass%data(k,current_state%column_local_y,current_state%column_local_x)
      qAccInsolpr(k)=current_state%qAccumInsolMass%data(k,current_state%column_local_y,current_state%column_local_x)
      qCoaSolpr(k)=current_state%qCoarseSolMass%data(k,current_state%column_local_y,current_state%column_local_x)
      qCoaDuspr(k)=current_state%qCoarseDustMass%data(k,current_state%column_local_y,current_state%column_local_x)
      nlpr(k)=current_state%nl%data(k,current_state%column_local_y,current_state%column_local_x)
      nrpr(k)=current_state%nr%data(k,current_state%column_local_y,current_state%column_local_x)
      nipr(k)=current_state%ni%data(k,current_state%column_local_y,current_state%column_local_x)
      nspr(k)=current_state%ns%data(k,current_state%column_local_y,current_state%column_local_x)
      ngpr(k)=current_state%ng%data(k,current_state%column_local_y,current_state%column_local_x)
      nAiktpr(k)=current_state%nAitkenSolNumber%data(k,current_state%column_local_y,current_state%column_local_x)
      nAccSolpr(k)=current_state%nAccumSolNumber%data(k,current_state%column_local_y,current_state%column_local_x)
      nAccInsolpr(k)=current_state%nAccumInsolNumber%data(k,current_state%column_local_y,current_state%column_local_x)
      nCoaSolpr(k)=current_state%nCoarseSolNumber%data(k,current_state%column_local_y,current_state%column_local_x)
      nCoaDuspr(k)=current_state%nCoarseDustnumber%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olqvbar)) &
        qvpr(k)=qvpr(k)-current_state%global_grid%configuration%vertical%olqvbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqlbar)) &
        qlpr(k)=qlpr(k)-current_state%global_grid%configuration%vertical%olqlbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqrbar)) &
        qrpr(k)=qrpr(k)-current_state%global_grid%configuration%vertical%olqrbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqibar)) &
        qipr(k)=qipr(k)-current_state%global_grid%configuration%vertical%olqibar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqsbar)) &
        qspr(k)=qspr(k)-current_state%global_grid%configuration%vertical%olqsbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqgbar)) &
        qgpr(k)=qgpr(k)-current_state%global_grid%configuration%vertical%olqgbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqAitkenSolMassbar)) &
        qAiktpr(k)=qAiktpr(k)-current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqAccumSolMassbar)) &
        qAccSolpr(k)=qAccSolpr(k)-current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqAccumInsolMassbar)) &
        qAccInsolpr(k)=qAccInsolpr(k)-current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqCoarseSolMassbar)) &
        qCoaSolpr(k)=qCoaSolpr(k)-current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olqCoarseDustMassbar)) &
        qCoaDuspr(k)=qCoaDuspr(k)-current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olnlbar)) &
        nlpr(k)=nlpr(k)-current_state%global_grid%configuration%vertical%olnlbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olnrbar)) &
        nrpr(k)=nrpr(k)-current_state%global_grid%configuration%vertical%olnrbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olnibar)) &
        nipr(k)=nipr(k)-current_state%global_grid%configuration%vertical%olnibar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olnsbar)) &
        nspr(k)=nspr(k)-current_state%global_grid%configuration%vertical%olnsbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olngbar)) &
        ngpr(k)=ngpr(k)-current_state%global_grid%configuration%vertical%olngbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar)) &
        nAiktpr(k)=nAiktpr(k)-current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olnAccumSolNumberbar)) &
        nAccSolpr(k)=nAccSolpr(k)-current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar)) &
        nAccInsolpr(k)=nAccInsolpr(k)-current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar)) &
        nCoaSolpr(k)=nCoaSolpr(k)-current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k)
      if (allocated(current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar)) &
        nCoaDuspr(k)=nCoaDuspr(k)-current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k)
    end do

    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(current_state%qv_flux_values)) current_state%qv_flux_values(k) = current_state%qv_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qvpr(k)+qvpr(k+1))
      if (allocated(current_state%ql_flux_values)) current_state%ql_flux_values(k) = current_state%ql_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qlpr(k)+qlpr(k+1))
      if (allocated(current_state%qr_flux_values)) current_state%qr_flux_values(k) = current_state%qr_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qrpr(k)+qrpr(k+1))
      if (allocated(current_state%qi_flux_values)) current_state%qi_flux_values(k) = current_state%qi_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qipr(k)+qipr(k+1))
      if (allocated(current_state%qs_flux_values)) current_state%qs_flux_values(k) = current_state%qs_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qspr(k)+qspr(k+1))
      if (allocated(current_state%qg_flux_values)) current_state%qg_flux_values(k) = current_state%qg_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qgpr(k)+qgpr(k+1))
      if (allocated(current_state%qAitkenSolMass_flux_values)) current_state%qAitkenSolMass_flux_values(k) = &
        current_state%qAitkenSolMass_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qAiktpr(k)+qAiktpr(k+1))
      if (allocated(current_state%qAitkenSolMass_flux_values)) current_state%qAccumSolMass_flux_values(k) = &
        current_state%qAccumSolMass_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qAccSolpr(k)+qAccSolpr(k+1))
      if (allocated(current_state%qAitkenSolMass_flux_values)) current_state%qAccumInsolMass_flux_values(k) = &
        current_state%qAccumInsolMass_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qAccInsolpr(k)+qAccInsolpr(k+1))
      if (allocated(current_state%qAitkenSolMass_flux_values)) current_state%qCoarseSolMass_flux_values(k) = &
        current_state%qCoarseSolMass_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qCoaSolpr(k)+qCoaSolpr(k+1))
      if (allocated(current_state%qAitkenSolMass_flux_values)) current_state%qCoarseDustMass_flux_values(k) = &
        current_state%qCoarseDustMass_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qCoaDuspr(k)+qCoaDuspr(k+1))
      if (allocated(current_state%nl_flux_values)) current_state%nl_flux_values(k) = current_state%nl_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(nlpr(k)+nlpr(k+1))
      if (allocated(current_state%nr_flux_values)) current_state%nr_flux_values(k) = current_state%nr_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(nrpr(k)+nrpr(k+1))
      if (allocated(current_state%ni_flux_values)) current_state%ni_flux_values(k) = current_state%ni_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(nipr(k)+nipr(k+1))
      if (allocated(current_state%ns_flux_values)) current_state%ns_flux_values(k) = current_state%ns_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(nspr(k)+nspr(k+1))
      if (allocated(current_state%ng_flux_values)) current_state%ng_flux_values(k) = current_state%ng_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(ngpr(k)+ngpr(k+1))
      if (allocated(current_state%nAitkenSolNumber_flux_values)) current_state%nAitkenSolNumber_flux_values(k) = &
        current_state%nAitkenSolNumber_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(nAiktpr(k)+nAiktpr(k+1))
      if (allocated(current_state%nAccumSolNumber_flux_values)) current_state%nAccumSolNumber_flux_values(k) = &
        current_state%nAccumSolNumber_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(nAccSolpr(k)+nAccSolpr(k+1))
      if (allocated(current_state%nAccumInsolNumber_flux_values)) current_state%nAccumInsolNumber_flux_values(k) = &
        current_state%nAccumInsolNumber_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(nAccInsolpr(k)+nAccInsolpr(k+1))
      if (allocated(current_state%nCoarseSolNumber_flux_values)) current_state%nCoarseSolNumber_flux_values(k) = &
        current_state%nCoarseSolNumber_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(nCoaSolpr(k)+nCoaSolpr(k+1))
      if (allocated(current_state%nCoarseSolNumber_flux_values)) current_state%nCoarseDustnumber_flux_values(k) = &
        current_state%nCoarseSolNumber_flux_values(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(nCoaDuspr(k)+nCoaDuspr(k+1))

      if (allocated(current_state%qv_tendency)) then
        current_state%qv_tendency(k) = current_state%qv_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqv%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qvpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%ql_tendency)) then
        current_state%ql_tendency(k) = current_state%ql_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sql%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qlpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%qr_tendency)) then
        current_state%qr_tendency(k) = current_state%qr_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqr%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qrpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%qi_tendency)) then
        current_state%qi_tendency(k) = current_state%qi_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqi%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qipr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%qs_tendency)) then
        current_state%qs_tendency(k) = current_state%qs_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqs%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qspr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%qg_tendency)) then
        current_state%qg_tendency(k) = current_state%qg_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqg%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qgpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%qAitkenSolMass_tendency)) then
        current_state%qAitkenSolMass_tendency(k) = current_state%qAitkenSolMass_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqAitkenSolMass%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qAiktpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%qAccumSolMass_tendency)) then
        current_state%qAccumSolMass_tendency(k) = current_state%qAccumSolMass_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqAccumSolMass%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qAccSolpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%qAccumInsolMass_tendency)) then
        current_state%qAccumInsolMass_tendency(k) = current_state%qAccumInsolMass_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqAccumInsolMass%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qAccInsolpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%qCoarseSolMass_tendency)) then
        current_state%qCoarseSolMass_tendency(k) = current_state%qCoarseSolMass_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqCoarseSolMass%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qCoaSolpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%qCoarseDustMass_tendency)) then
        current_state%qCoarseDustMass_tendency(k) = current_state%qCoarseDustMass_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sqCoarseDustMass%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        qCoaDuspr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%nl_tendency)) then
        current_state%nl_tendency(k) = current_state%nl_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%snl%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        nlpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%nr_tendency)) then
        current_state%nr_tendency(k) = current_state%nr_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%snr%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        nrpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%ni_tendency)) then
        current_state%ni_tendency(k) = current_state%ni_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sni%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        nipr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%ng_tendency)) then
        current_state%ng_tendency(k) = current_state%ng_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sng%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        ngpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%ns_tendency)) then
        current_state%ns_tendency(k) = current_state%ns_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%sns%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        nspr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%nAitkenSolNumber_tendency)) then
        current_state%nAitkenSolNumber_tendency(k) = current_state%nAitkenSolNumber_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%snAitkenSolNumber%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        nAiktpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%nAccumSolNumber_tendency)) then
        current_state%nAccumSolNumber_tendency(k) = current_state%nAccumSolNumber_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%snAccumSolNumber%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        nAccSolpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%nAccumInsolNumber_tendency)) then
        current_state%nAccumInsolNumber_tendency(k) = current_state%nAccumInsolNumber_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%snAccumInsolNumber%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        nAccInsolpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%nCoarseSolNumber_tendency)) then
        current_state%nCoarseSolNumber_tendency(k) = current_state%nCoarseSolNumber_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%snCoarseSolNumber%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        nCoaSolpr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if
      if (allocated(current_state%nCoarseDustnumber_tendency)) then
        current_state%nCoarseDustnumber_tendency(k) = current_state%nCoarseDustnumber_tendency(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
        current_state%snCoarseDustnumber%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*&
        nCoaDuspr(k)*(current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))
      end if

      if (allocated(current_state%qv_gradient)) then
          current_state%qv_gradient(k) = current_state%qv_gradient(k)+qvpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qv_advection(k)
      end if
      if (allocated(current_state%ql_gradient)) then
          current_state%ql_gradient(k) = current_state%ql_gradient(k)+qlpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%ql_advection(k)
      end if
      if (allocated(current_state%qr_gradient)) then
          current_state%qr_gradient(k) = current_state%qr_gradient(k)+qrpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qr_advection(k)
      end if
      if (allocated(current_state%qi_gradient)) then
          current_state%qi_gradient(k) = current_state%qi_gradient(k)+qipr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qi_advection(k)
      end if
      if (allocated(current_state%qs_gradient)) then
          current_state%qs_gradient(k) = current_state%qs_gradient(k)+qspr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qs_advection(k)
      end if
      if (allocated(current_state%qg_gradient)) then
          current_state%qg_gradient(k) = current_state%qg_gradient(k)+qgpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qg_advection(k)
      end if
      if (allocated(current_state%qAitkenSolMass_gradient)) then
          current_state%qAitkenSolMass_gradient(k) = current_state%qAitkenSolMass_gradient(k) + &
          qAiktpr(k)*0.5*(current_state%w_advection(k) + current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qAitkenSolMass_advection(k)
      end if
      if (allocated(current_state%qAccumSolMass_gradient)) then
          current_state%qAccumSolMass_gradient(k) = current_state%qAccumSolMass_gradient(k) + &
          qAccSolpr(k)*0.5*(current_state%w_advection(k) + current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qAccumSolMass_advection(k)
      end if
      if (allocated(current_state%qAccumInsolMass_gradient)) then
          current_state%qAccumInsolMass_gradient(k) = current_state%qAccumInsolMass_gradient(k) + &
          qAccInsolpr(k)*0.5*(current_state%w_advection(k) + current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qAccumInsolMass_advection(k)
      end if
      if (allocated(current_state%qCoarseSolMass_gradient)) then
          current_state%qCoarseSolMass_gradient(k) = current_state%qCoarseSolMass_gradient(k) + &
          qCoaSolpr(k)*0.5*(current_state%w_advection(k) + current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qCoarseSolMass_advection(k)
      end if
      if (allocated(current_state%qCoarseDustMass_gradient)) then
          current_state%qCoarseDustMass_gradient(k) = current_state%qCoarseDustMass_gradient(k) + &
          qCoaDuspr(k)*0.5*(current_state%w_advection(k) + current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qCoarseDustMass_advection(k)
      end if
      if (allocated(current_state%nl_gradient)) then
          current_state%nl_gradient(k) = current_state%nl_gradient(k)+nlpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nl_advection(k)
      end if
      if (allocated(current_state%nr_gradient)) then
          current_state%nr_gradient(k) = current_state%nr_gradient(k)+nrpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nr_advection(k)
      end if
      if (allocated(current_state%ni_gradient)) then
          current_state%ni_gradient(k) = current_state%ni_gradient(k)+nipr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%ni_advection(k)
      end if
      if (allocated(current_state%ns_gradient)) then
          current_state%ns_gradient(k) = current_state%ns_gradient(k)+nspr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%ns_advection(k)
      end if
      if (allocated(current_state%ng_gradient)) then
          current_state%ng_gradient(k) = current_state%ng_gradient(k)+ngpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%ng_advection(k)
      end if
      if (allocated(current_state%nAitkenSolNumber_gradient)) then
          current_state%nAitkenSolNumber_gradient(k) = current_state%nAitkenSolNumber_gradient(k) + &
          nAiktpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nAitkenSolNumber_advection(k)
      end if
      if (allocated(current_state%nAccumSolNumber_gradient)) then
          current_state%nAccumSolNumber_gradient(k) = current_state%nAccumSolNumber_gradient(k) + &
          nAccSolpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nAccumSolNumber_advection(k)
      end if
      if (allocated(current_state%nAccumInsolNumber_gradient)) then
          current_state%nAccumInsolNumber_gradient(k) = current_state%nAccumInsolNumber_gradient(k) + &
          nAccInsolpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nAccumInsolNumber_advection(k)
      end if
      if (allocated(current_state%nCoarseSolNumber_gradient)) then
          current_state%nCoarseSolNumber_gradient(k) = current_state%nCoarseSolNumber_gradient(k) + &
          nCoaSolpr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nCoarseSolNumber_advection(k)
      end if
      if (allocated(current_state%nCoarseDustnumber_gradient)) then
          current_state%nCoarseDustnumber_gradient(k) = current_state%nCoarseDustnumber_gradient(k) + &
          nCoaDuspr(k)*0.5*(current_state%w_advection(k)+&
          current_state%w_advection(k-1))+0.5*&
          (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nCoarseDustnumber_advection(k)
      end if

      if (allocated(current_state%qv_diff)) then
          current_state%qv_diff(k) = current_state%qv_diff(k)+qvpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qv_diffusion(k)
      end if
      if (allocated(current_state%ql_diff)) then
          current_state%ql_diff(k) = current_state%ql_diff(k)+qlpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%ql_diffusion(k)
      end if
      if (allocated(current_state%qr_diff)) then
          current_state%qr_diff(k) = current_state%qr_diff(k)+qrpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qr_diffusion(k)
      end if
      if (allocated(current_state%qi_diff)) then
          current_state%qi_diff(k) = current_state%qi_diff(k)+qipr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qi_diffusion(k)
      end if
      if (allocated(current_state%qg_diff)) then
          current_state%qg_diff(k) = current_state%qg_diff(k)+qgpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qg_diffusion(k)
      end if
      if (allocated(current_state%qs_diff)) then
          current_state%qs_diff(k) = current_state%qs_diff(k)+qspr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qs_diffusion(k)
      end if
      if (allocated(current_state%qAitkenSolMass_diff)) then
          current_state%qAitkenSolMass_diff(k) = current_state%qAitkenSolMass_diff(k) + &
          qAiktpr(k)*0.5*(current_state%w_viscosity(k) + current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qAitkenSolMass_diffusion(k)
      end if
      if (allocated(current_state%qAccumSolMass_diff)) then
          current_state%qAccumSolMass_diff(k) = current_state%qAccumSolMass_diff(k) + &
          qAccSolpr(k)*0.5*(current_state%w_viscosity(k) + current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qAccumSolMass_diffusion(k)
      end if
      if (allocated(current_state%qAccumInsolMass_diff)) then
          current_state%qAccumInsolMass_diff(k) = current_state%qAccumInsolMass_diff(k) + &
          qAccInsolpr(k)*0.5*(current_state%w_viscosity(k) + current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qAccumInsolMass_diffusion(k)
      end if
      if (allocated(current_state%qCoarseSolMass_diff)) then
          current_state%qCoarseSolMass_diff(k) = current_state%qCoarseSolMass_diff(k) + &
          qCoaSolpr(k)*0.5*(current_state%w_viscosity(k) + current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qCoarseSolMass_diffusion(k)
      end if
      if (allocated(current_state%qCoarseDustMass_diff)) then
          current_state%qCoarseDustMass_diff(k) = current_state%qCoarseDustMass_diff(k) + &
          qCoaDuspr(k)*0.5*(current_state%w_viscosity(k) + current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%qCoarseDustMass_diffusion(k)
      end if
      if (allocated(current_state%nl_diff)) then
          current_state%nl_diff(k) = current_state%nl_diff(k)+nlpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nl_diffusion(k)
      end if
      if (allocated(current_state%nr_diff)) then
          current_state%nr_diff(k) = current_state%nr_diff(k)+nrpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nr_diffusion(k)
      end if
      if (allocated(current_state%ni_diff)) then
          current_state%ni_diff(k) = current_state%ni_diff(k)+nipr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%ni_diffusion(k)
      end if
      if (allocated(current_state%ng_diff)) then
          current_state%ng_diff(k) = current_state%ng_diff(k)+ngpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%ng_diffusion(k)
      end if
      if (allocated(current_state%ns_diff)) then
          current_state%ns_diff(k) = current_state%ns_diff(k)+nspr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%ns_diffusion(k)
      end if
      if (allocated(current_state%nAitkenSolNumber_diff)) then
          current_state%nAitkenSolNumber_diff(k) = current_state%nAitkenSolNumber_diff(k) + &
          nAiktpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nAitkenSolNumber_diffusion(k)
      end if
      if (allocated(current_state%nAccumSolNumber_diff)) then
          current_state%nAccumSolNumber_diff(k) = current_state%nAccumSolNumber_diff(k) + &
          nAccSolpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nAccumSolNumber_diffusion(k)
      end if
      if (allocated(current_state%nAccumInsolNumber_diff)) then
          current_state%nAccumInsolNumber_diff(k) = current_state%nAccumInsolNumber_diff(k) + &
          nAccInsolpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nAccumInsolNumber_diffusion(k)
      end if
      if (allocated(current_state%nCoarseSolNumber_diff)) then
          current_state%nCoarseSolNumber_diff(k) = current_state%nCoarseSolNumber_diff(k) + &
          nCoaSolpr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nCoarseSolNumber_diffusion(k)
      end if
      if (allocated(current_state%nCoarseDustnumber_diff)) then
          current_state%nCoarseDustnumber_diff(k) = current_state%nCoarseDustnumber_diff(k) + &
          nCoaDuspr(k)*0.5*(current_state%w_viscosity(k)+&
          current_state%w_viscosity(k-1))+&
          0.5*(current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*&
          current_state%nCoarseDustnumber_diffusion(k)
      end if

      if (allocated(current_state%qv_buoyancy)) then
          current_state%qv_buoyancy(k) = current_state%qv_buoyancy(k)+qvpr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%ql_buoyancy)) then
          current_state%ql_buoyancy(k) = current_state%ql_buoyancy(k)+qlpr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%qr_buoyancy)) then
          current_state%qr_buoyancy(k) = current_state%qr_buoyancy(k)+qrpr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%qi_buoyancy)) then
          current_state%qi_buoyancy(k) = current_state%qi_buoyancy(k)+qipr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%qs_buoyancy)) then
          current_state%qs_buoyancy(k) = current_state%qs_buoyancy(k)+qspr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%qg_buoyancy)) then
          current_state%qg_buoyancy(k) = current_state%qg_buoyancy(k)+qgpr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%qAitkenSolMass_buoyancy)) then
          current_state%qAitkenSolMass_buoyancy(k) = current_state%qAitkenSolMass_buoyancy(k) + &
          qAiktpr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%qAccumSolMass_buoyancy)) then
          current_state%qAccumSolMass_buoyancy(k) = current_state%qAccumSolMass_buoyancy(k) + &
          qAccSolpr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%qAccumInsolMass_buoyancy)) then
          current_state%qAccumInsolMass_buoyancy(k) = current_state%qAccumInsolMass_buoyancy(k) + &
          qAccInsolpr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%qCoarseSolMass_buoyancy)) then
          current_state%qCoarseSolMass_buoyancy(k) = current_state%qCoarseSolMass_buoyancy(k) + &
          qCoaSolpr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%qCoarseDustMass_buoyancy)) then
          current_state%qCoarseDustMass_buoyancy(k) = current_state%qCoarseDustMass_buoyancy(k) + &
          qCoaDuspr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%nl_buoyancy)) then
          current_state%nl_buoyancy(k) = current_state%nl_buoyancy(k)+nlpr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%nr_buoyancy)) then
          current_state%nr_buoyancy(k) = current_state%nr_buoyancy(k)+nrpr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%ni_buoyancy)) then
          current_state%ni_buoyancy(k) = current_state%ni_buoyancy(k)+nipr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%ns_buoyancy)) then
          current_state%ns_buoyancy(k) = current_state%ns_buoyancy(k)+nspr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%ng_buoyancy)) then
          current_state%ng_buoyancy(k) = current_state%ng_buoyancy(k)+ngpr(k)*0.5*(current_state%w_buoyancy(k)+&
          current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%nAitkenSolNumber_buoyancy)) then
          current_state%nAitkenSolNumber_buoyancy(k) = current_state%nAitkenSolNumber_buoyancy(k) + &
          nAiktpr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%nAccumSolNumber_buoyancy)) then
          current_state%nAccumSolNumber_buoyancy(k) = current_state%nAccumSolNumber_buoyancy(k) + &
          nAccSolpr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%nAccumInsolNumber_buoyancy)) then
          current_state%nAccumInsolNumber_buoyancy(k) = current_state%nAccumInsolNumber_buoyancy(k) + &
          nAccInsolpr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%nCoarseSolNumber_buoyancy)) then
          current_state%nCoarseSolNumber_buoyancy(k) = current_state%nCoarseSolNumber_buoyancy(k) + &
          nCoaSolpr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
      if (allocated(current_state%nCoarseDustnumber_buoyancy)) then
          current_state%nCoarseDustnumber_buoyancy(k) = current_state%nCoarseDustnumber_buoyancy(k) + &
          nCoaDuspr(k)*0.5*(current_state%w_buoyancy(k) + current_state%w_buoyancy(k-1))
      end if
    end do
  end subroutine compute_q_flux_for_column


  !> Computes the uw uv diagnostics for a specific column
  !! @param current_state Current model state
  subroutine compute_uw_vw_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1
    integer :: k

    do k=1, current_state%local_grid%size(Z_INDEX)
      upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
      uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
      if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
        upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
        uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
      end if
      vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
      vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
        vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
        vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
      end if
    end do

    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(current_state%uw_advection)) current_state%uw_advection(k) = &
        current_state%uw_advection(k)+0.5*(upr(k)+uprm1(k))*0.5*&
        (current_state%w_advection(k) + current_state%w_advection(k-1))+0.25*(&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*current_state%u_advection(k)
      if (allocated(current_state%vw_advection)) current_state%vw_advection(k) = &
        current_state%vw_advection(k)+0.5*(vpr(k)+vprm1(k))*0.5*&
        (current_state%w_advection(k) + current_state%w_advection(k-1))+0.25*(&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*current_state%v_advection(k)
      if (allocated(current_state%uw_viscosity)) current_state%uw_viscosity(k) = &
        current_state%uw_viscosity(k)+0.5*(upr(k)+uprm1(k))*0.5*&
        (current_state%w_viscosity(k)+current_state%w_viscosity(k-1))+0.25*(&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*current_state%u_viscosity(k)
      if (allocated(current_state%vw_viscosity)) current_state%vw_viscosity(k) = &
        current_state%vw_viscosity(k)+0.5*(vpr(k)+vprm1(k))*0.5*&
        (current_state%w_viscosity(k)+current_state%w_viscosity(k-1))+0.25*(&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*current_state%v_viscosity(k)
      if (allocated(current_state%uw_buoyancy)) current_state%uw_buoyancy(k) = &
        current_state%uw_buoyancy(k)+0.5*(upr(k)+uprm1(k))*0.5*(&
        current_state%w_buoyancy(k)+current_state%w_buoyancy(k-1))
      if (allocated(current_state%vw_buoyancy)) current_state%vw_buoyancy(k) = &
        current_state%vw_buoyancy(k)+0.5*(vpr(k)+vprm1(k))*0.5*(&
        current_state%w_buoyancy(k)+current_state%w_buoyancy(k-1))
      if (allocated(current_state%uw_tendency)) current_state%uw_tendency(k) = &
        current_state%uw_tendency(k)+0.25*(&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*&
        current_state%su%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(&
        current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(upr(k)+uprm1(k))
      if (allocated(current_state%vw_tendency)) current_state%vw_tendency(k) = current_state%vw_tendency(k)+0.25*(&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x+1)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x+1))*&
        current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(&
        current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sw%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(vpr(k)+vprm1(k))
      if (allocated(current_state%uw_w)) current_state%uw_w(k) = current_state%uw_w(k)+0.25*(upr(k)+upr(k+1)+uprm1(k)+uprm1(k+1))*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(current_state%vw_w)) current_state%vw_w(k) = current_state%vw_w(k)+0.25*(vpr(k)+vpr(k+1)+vprm1(k)+vprm1(k+1))*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)
    end do

  end subroutine compute_uw_vw_for_column

  !> Computes the prognostic (uu, vv, ww) budgets for a specific column
  !! @param current_state The current model state
  subroutine compute_prognostic_budgets_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1
    integer :: k

    do k=1, current_state%local_grid%size(Z_INDEX)
      upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
      uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
      if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
        upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
        uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
      end if
      vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
      vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
        vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
        vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
      end if
    end do

    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(current_state%tu_su)) current_state%tu_su(k) = current_state%tu_su(k) + &
        2.0*upr(k)*current_state%su%data(k,current_state%column_local_y, current_state%column_local_x)
      if (allocated(current_state%uu_advection)) current_state%uu_advection(k) = &
        current_state%uu_advection(k) + 2.0*upr(k)*current_state%u_advection(k)
      if (allocated(current_state%uu_viscosity)) current_state%uu_viscosity(k) = &
        current_state%uu_viscosity(k) + 2.0*upr(k)*current_state%u_viscosity(k)
      if (allocated(current_state%wu_u)) current_state%wu_u(k) = &
        current_state%wu_u(k)+0.25*(upr(k)+upr(k+1)+uprm1(k)+uprm1(k+1))*0.25*&
        (upr(k)+upr(k+1)+uprm1(k)+uprm1(k+1))*current_state%w%data(k,current_state%column_local_y,&
        current_state%column_local_x)

      if (allocated(current_state%tv_sv)) current_state%tv_sv(k) = current_state%tv_sv(k) + &
        2.0*vpr(k)*current_state%sv%data(k,current_state%column_local_y, current_state%column_local_x)
      if (allocated(current_state%vv_advection)) current_state%vv_advection(k) = &
        current_state%vv_advection(k) + 2.0*vpr(k)*current_state%v_advection(k)
      if (allocated(current_state%vv_viscosity)) current_state%uu_viscosity(k) = &
        current_state%vv_viscosity(k) + 2.0*vpr(k)*current_state%v_viscosity(k)
      if (allocated(current_state%wv_v)) current_state%wv_v(k) = &
        current_state%wv_v(k)+0.25*(vpr(k)+vpr(k+1)+vprm1(k)+vprm1(k+1))*0.25*&
        (vpr(k)+vpr(k+1)+vprm1(k)+vprm1(k+1))*current_state%w%data(k,current_state%column_local_y,&
        current_state%column_local_x)

      if (allocated(current_state%tw_sw)) current_state%tw_sw(k) = current_state%tw_sw(k) + &
        2.0*current_state%w%data(k,current_state%column_local_y, current_state%column_local_x)*&
        current_state%sw%data(k,current_state%column_local_y, current_state%column_local_x)
      if (allocated(current_state%ww_advection)) current_state%ww_advection(k) = current_state%ww_advection(k) + &
        2.0*current_state%w%data(k,current_state%column_local_y, current_state%column_local_x)*current_state%w_advection(k)
      if (allocated(current_state%ww_viscosity)) current_state%ww_viscosity(k) = current_state%ww_viscosity(k) + &
        2.0*current_state%w%data(k,current_state%column_local_y, current_state%column_local_x)*current_state%w_viscosity(k)
      if (allocated(current_state%ww_buoyancy)) current_state%ww_buoyancy(k) = current_state%ww_buoyancy(k) + &
        2.0*current_state%w%data(k,current_state%column_local_y, current_state%column_local_x)*current_state%w_buoyancy(k)
    end do

  end subroutine compute_prognostic_budgets_for_column

  !> Computes the thetal diagnostics for a specific column
  !! @param current_state The current model state
  subroutine compute_thetal_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1, thlpr, thlprp1
    integer :: k

    do k=1, current_state%local_grid%size(Z_INDEX)
      upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
      uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
      if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
        upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
        uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
      end if
      vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
      vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
      if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
        vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
        vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
      end if
      thlpr(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x)
      thlprp1(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x+1)
      if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
        thlpr(k)=thlpr(k)-current_state%global_grid%configuration%vertical%olthbar(k)
        thlprp1(k)=thlprp1(k)-current_state%global_grid%configuration%vertical%olthbar(k)
      end if
    end do

    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(current_state%u_thetal)) current_state%u_thetal(k) = current_state%u_thetal(k)+0.5*(upr(k)+uprm1(k))*thlpr(k)
      if (allocated(current_state%us_thetal)) current_state%us_thetal(k) = current_state%us_thetal(k)+0.5*(upr(k)+uprm1(k))*&
        current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(thlpr(k)+thlprp1(k))*&
        current_state%su%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(current_state%u_thetal_advection)) current_state%u_thetal_advection(k) = &
        current_state%u_thetal_advection(k)+0.5*(upr(k)+uprm1(k))*&
        current_state%th_advection(k)+0.5*(thlpr(k)+thlprp1(k))*current_state%u_advection(k)
      if (allocated(current_state%u_thetal_viscosity_diffusion)) current_state%u_thetal_viscosity_diffusion(k) = &
        current_state%u_thetal_viscosity_diffusion(k)+0.5*&
        (thlpr(k)+thlprp1(k))*current_state%u_viscosity(k)+0.5*(upr(k)+uprm1(k))*current_state%th_diffusion(k)
      if (allocated(current_state%wu_thetal)) current_state%wu_thetal(k) = current_state%wu_thetal(k)+&
         current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(upr(k+1)+upr(k)+&
         uprm1(k+1)+uprm1(k))*0.5*(thlpr(k+1)+thlpr(k))

      if (allocated(current_state%v_thetal)) current_state%v_thetal(k) = current_state%v_thetal(k)+0.5*(vpr(k)+vprm1(k))*thlpr(k)
      if (allocated(current_state%vs_thetal)) current_state%vs_thetal(k) = current_state%vs_thetal(k)+0.5*(vpr(k)+vprm1(k))*&
        current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(thlpr(k)+thlprp1(k))*&
        current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(current_state%v_thetal_advection)) current_state%v_thetal_advection(k) = &
        current_state%v_thetal_advection(k)+0.5*(vpr(k)+vprm1(k))*&
        current_state%th_advection(k)+0.5*(thlpr(k)+thlprp1(k))*current_state%v_advection(k)
      if (allocated(current_state%v_thetal_viscosity_diffusion)) current_state%v_thetal_viscosity_diffusion(k) = &
        current_state%v_thetal_viscosity_diffusion(k)+0.5*&
        (thlpr(k)+thlprp1(k))*current_state%v_viscosity(k)+0.5*(vpr(k)+vprm1(k))*current_state%th_diffusion(k)
      if (allocated(current_state%wv_thetal)) current_state%wv_thetal(k) = current_state%wv_thetal(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(vpr(k+1)+vpr(k)+&
        vprm1(k+1)+vprm1(k))*0.5*(thlpr(k+1)+thlpr(k))

      if (allocated(current_state%w_thetal)) current_state%w_thetal(k) = current_state%w_thetal(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(thlpr(k)+thlpr(k+1))
      if (allocated(current_state%ws_thetal)) current_state%ws_thetal(k) = current_state%ws_thetal(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
        (current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%sth%data(k+1,current_state%column_local_y,current_state%column_local_x))+0.5*(thlpr(k)+thlpr(k+1))*&
        current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(current_state%w_thetal_advection)) current_state%w_thetal_advection(k) = current_state%w_thetal_advection(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(current_state%th_advection(k)+&
        current_state%th_advection(k+1))+0.5*(thlpr(k)+thlpr(k+1))*current_state%w_advection(k)
      if (allocated(current_state%w_thetal_viscosity_diffusion)) current_state%w_thetal_viscosity_diffusion(k) = &
        current_state%w_thetal_viscosity_diffusion(k)+0.5*&
        (thlpr(k)+thlpr(k+1))*current_state%w_viscosity(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
        (current_state%th_diffusion(k)+current_state%th_diffusion(k+1))
      if (allocated(current_state%w_thetal_buoyancy)) current_state%w_thetal_buoyancy(k) = &
        current_state%w_thetal_buoyancy(k)+0.5*(thlpr(k+1)+thlpr(k))*current_state%w_buoyancy(k)
      if (allocated(current_state%ww_thetal)) current_state%ww_thetal(k) = current_state%ww_thetal(k)+0.5*&
        (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
        current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*thlpr(k)

      if (allocated(current_state%thetal_thetal)) current_state%thetal_thetal(k) = current_state%thetal_thetal(k)+thlpr(k)*thlpr(k)
      if (allocated(current_state%sthetal_thetal)) current_state%sthetal_thetal(k) = current_state%sthetal_thetal(k)+2.0*thlpr(k)*&
        current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)
      if (allocated(current_state%thetal_thetal_advection)) current_state%thetal_thetal_advection(k) = &
        current_state%thetal_thetal_advection(k)+2.0*thlpr(k)*current_state%th_advection(k)
      if (allocated(current_state%thetal_thetal_diffusion)) current_state%thetal_thetal_diffusion(k) = &
        current_state%thetal_thetal_diffusion(k)+2.0*thlpr(k)*current_state%th_diffusion(k)
      if (allocated(current_state%wthetal_thetal)) current_state%wthetal_thetal(k) = current_state%wthetal_thetal(k)+&
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(thlpr(k+1)+thlpr(k))*0.5*(&
        thlpr(k+1)+thlpr(k))
    end do

  end subroutine compute_thetal_for_column

  !> Computes the tke diagnostics for a specific column.
  !! @param current_state The current model state
  subroutine compute_tke_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1, &
      uu_tendency,vv_tendency,ww_tendency
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: umean, wu_umean, vmean, wv_vmean, &
      w_pprime_at_p, rke1, w_qvprime_at_w, w_qclprime_at_w, w_thprime_at_w, wq, rho, rec_rho, rhon, rec_rhon, &
      uw_tot, vw_tot,w_upr_at_w,w_vpr_at_w, sg_w_buoyancy
    real(kind=DEFAULT_PRECISION) :: u_at_p, v_at_p, w_at_p
    real(kind=DEFAULT_PRECISION) ::  C_virtual

    integer :: k, n

    C_virtual = (ratio_mol_wts-1.0_DEFAULT_PRECISION)

    !Resolved diagnostics
! ***********************Buoyant production 1/2 ***************************

    do k=1, current_state%local_grid%size(Z_INDEX)
      rho(k)=current_state%global_grid%configuration%vertical%rho(k)
      rhon(k)=current_state%global_grid%configuration%vertical%rhon(k)
      rec_rho(k)=1.0_DEFAULT_PRECISION/rho(k)
      rec_rhon(k)=1.0_DEFAULT_PRECISION/rhon(k)
    end do

! ***********************TKE Tendency ***************************

    do k=2, current_state%local_grid%size(Z_INDEX)

      uu_tendency(k) = ((current_state%u%data(k,current_state%column_local_y,current_state%column_local_x) - &
                     current_state%global_grid%configuration%vertical%olubar(k) )**2 - &
                     (current_state%zu%data(k,current_state%column_local_y,current_state%column_local_x) - &
                     current_state%global_grid%configuration%vertical%olzubar(k) )**2 ) / &
                     current_state%dtm

      vv_tendency(k) = ((current_state%v%data(k,current_state%column_local_y,current_state%column_local_x) - &
                       current_state%global_grid%configuration%vertical%olvbar(k) )**2 - &
                       (current_state%zv%data(k,current_state%column_local_y,current_state%column_local_x) - &
                       current_state%global_grid%configuration%vertical%olzvbar(k) )**2 ) / &
                       current_state%dtm

      ww_tendency(k) = (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) * &
                       current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) - &
                       current_state%zw%data(k,current_state%column_local_y,current_state%column_local_x) * &
                       current_state%zw%data(k,current_state%column_local_y,current_state%column_local_x)) / &
                       current_state%dtm

    end do

    uu_tendency(1) = -uu_tendency(2)
    vv_tendency(1) = -vv_tendency(2)

    if (allocated(current_state%tke_tend)) then
      do k=2, current_state%local_grid%size(Z_INDEX)-1
        current_state%tke_tend(k) = current_state%tke_tend(k) + 0.5_DEFAULT_PRECISION * (&
          0.5_DEFAULT_PRECISION * (uu_tendency(k)+uu_tendency(k+1)) + &
          0.5_DEFAULT_PRECISION * (vv_tendency(k)+vv_tendency(k+1)) + &
          ww_tendency(k) )
      end do
    end if

! ***********************Shear production ***************************

    do k=2, current_state%local_grid%size(Z_INDEX)-1

      umean(k)=(current_state%global_grid%configuration%vertical%olubar(k+1) -&
               current_state%global_grid%configuration%vertical%olubar(k))* &
               current_state%global_grid%configuration%vertical%rdzn(k+1)

      w_upr_at_w(k) =current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) * &
            (0.25_DEFAULT_PRECISION * ( &
            (current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)  - &
            current_state%global_grid%configuration%vertical%olubar(k)) + &
            (current_state%u%data(k+1,current_state%column_local_y,current_state%column_local_x)  - &
            current_state%global_grid%configuration%vertical%olubar(k+1)) + &
            (current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)    - &
            current_state%global_grid%configuration%vertical%olubar(k)) + &
            (current_state%u%data(k+1,current_state%column_local_y,current_state%column_local_x-1)- &
            current_state%global_grid%configuration%vertical%olubar(k+1)) ) )

      vmean(k)=(current_state%global_grid%configuration%vertical%olvbar(k+1) - &
               current_state%global_grid%configuration%vertical%olvbar(k)) * &
               current_state%global_grid%configuration%vertical%rdzn(k+1)

      w_vpr_at_w(k) =current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) * &
            (0.25_DEFAULT_PRECISION * ( &
            (current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)  - &
            current_state%global_grid%configuration%vertical%olvbar(k)) + &
            (current_state%v%data(k+1,current_state%column_local_y,current_state%column_local_x)  - &
            current_state%global_grid%configuration%vertical%olvbar(k+1)) + &
            (current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)    - &
            current_state%global_grid%configuration%vertical%olvbar(k)) + &
            (current_state%v%data(k+1,current_state%column_local_y-1,current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olvbar(k+1)) ) )

      wu_umean(k)=(w_upr_at_w(k)*umean(k))
      wv_vmean(k)= (w_vpr_at_w(k)*vmean(k))

      if (allocated(current_state%shprd)) then
        current_state%shprd(k) = current_state%shprd(k) - (wv_vmean(k) + wu_umean(k))
      end if

    end do

    current_state%shprd(1) = 0.0_DEFAULT_PRECISION
    current_state%shprd(current_state%local_grid%size(Z_INDEX)) = 0.0_DEFAULT_PRECISION


! *********************** Pressure transport ***************************
    do k=2, current_state%local_grid%size(Z_INDEX)

      !In current state - p=p/rho (rho here is on same levels as p)
      ! Note - calculating on z levels (i.e. w)
      ! So need w'p' on p levels

      w_pprime_at_p(k) = 0.5_DEFAULT_PRECISION * &
          (current_state%w%data(k,  current_state%column_local_y,current_state%column_local_x) + &
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x)) * &
          (current_state%global_grid%configuration%vertical%rhon(k) * &
          current_state%p%data(k,current_state%column_local_y,current_state%column_local_x))

      u_at_p = 0.5_DEFAULT_PRECISION * &
          ((current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)- &
          current_state%global_grid%configuration%vertical%olubar(k)) + &
          (current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)  - &
          current_state%global_grid%configuration%vertical%olubar(k)))

      v_at_p = 0.5_DEFAULT_PRECISION * &
          ((current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)- &
          current_state%global_grid%configuration%vertical%olvbar(k)) + &
          (current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)  - &
          current_state%global_grid%configuration%vertical%olvbar(k)))

      w_at_p = 0.5_DEFAULT_PRECISION  * &
          (current_state%w%data(k,  current_state%column_local_y,current_state%column_local_x) + &
          current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))

      rke1(k)= 0.5_DEFAULT_PRECISION * w_at_p * &
          ( u_at_p*u_at_p + v_at_p*v_at_p + w_at_p*w_at_p) * rec_rhon(k)

    end do

    w_pprime_at_p(current_state%local_grid%size(Z_INDEX)) = 0.0_DEFAULT_PRECISION
    ! Zero gradient at surface
    w_pprime_at_p(1)=w_pprime_at_p(2)
    if (allocated(current_state%w_p)) then
      do k=1, current_state%local_grid%size(Z_INDEX)-1
        current_state%w_p(k)= current_state%w_p(k) - ((w_pprime_at_p(k+1) - w_pprime_at_p(k)) * &
          current_state%global_grid%configuration%vertical%rdzn(k+1) * rec_rho(k))
      end do
    end if

! ********************** Resolved turbulent transport ************************

    rke1(current_state%local_grid%size(Z_INDEX)) = 0.0_DEFAULT_PRECISION
    ! Zero gradient at surface
    rke1(1)=rke1(2)

    if (allocated(current_state%w_ke)) then
      do k=1, current_state%local_grid%size(Z_INDEX)-1
        current_state%w_ke(k) = current_state%w_ke(k) -(rho(k) * (rke1(k+1) - rke1(k) ) * &
          current_state%global_grid%configuration%vertical%rdzn(k+1))
      end do
    end if

! *********************** Subgrid buoyant production***************************
! !!!Using buoyancy.F90

#ifdef W_ACTIVE
    if (.not. current_state%passive_th .and. current_state%th%active) then
      do k=2,current_state%local_grid%size(Z_INDEX)-1
        sg_w_buoyancy(k)=(0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
             ((current_state%th%data(k, current_state%column_local_y, current_state%column_local_x) -       &
               current_state%global_grid%configuration%vertical%olthbar(k)) +                               &
              (current_state%th%data(k+1, current_state%column_local_y, current_state%column_local_x) -     &
               current_state%global_grid%configuration%vertical%olthbar(k+1)))
       end do
    end if

    if (.not. current_state%passive_q) then
      if (current_state%use_anelastic_equations) then
        do k=2,current_state%local_grid%size(Z_INDEX)-1
          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(1) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qv%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqvbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qv%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqvbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(2) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%ql%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqlbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%ql%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqlbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(3) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qr%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqrbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qr%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqrbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(4) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qi%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqibar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qi%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqibar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(5) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qs%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqsbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qs%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqsbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(6) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qg%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqgbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qg%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqgbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(7) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qAitkenSolMass%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qAitkenSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(8) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qAccumSolMass%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qAccumSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(9) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qAccumInsolMass%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qAccumInsolMass%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(10) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qCoarseSolMass%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qCoarseSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(11) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%qCoarseDustMass%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%qCoarseDustMass%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(12) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%nl%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnlbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%nl%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnlbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(13) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%nr%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnrbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%nr%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnrbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(14) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%ni%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnibar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%ni%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnibar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(15) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%ns%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnsbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%ns%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnsbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(16) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%ng%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olngbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%ng%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olngbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(17) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%nAitkenSolNumber%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%nAitkenSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(18) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%nAccumSolNumber%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%nAccumSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(19) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%nAccumInsolNumber%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%nAccumInsolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(20) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%nCoarseSolNumber%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%nCoarseSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k+1)))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            (0.5_DEFAULT_PRECISION*current_state%global_grid%configuration%vertical%buoy_co(k)) * &
            current_state%cq(21) * &
            (current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%nCoarseDustnumber%data(k, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k)) + &
            current_state%global_grid%configuration%vertical%thref(k+1) * &
            (current_state%nCoarseDustnumber%data(k+1, current_state%column_local_y, current_state%column_local_x) - &
            current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k+1)))
        end do
      else
        do k=2,current_state%local_grid%size(Z_INDEX)-1
          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(1)*&
            (current_state%qv%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqvbar(k) + &
            current_state%qv%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqvbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(2)*&
            (current_state%ql%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqlbar(k) + &
            current_state%ql%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqlbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(3)*&
            (current_state%qr%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqrbar(k) + &
            current_state%qr%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqrbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(4)*&
            (current_state%qi%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqibar(k) + &
            current_state%qi%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqibar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(5)*&
            (current_state%qs%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqsbar(k) + &
            current_state%qs%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqsbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(6)*&
            (current_state%qg%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqgbar(k) + &
            current_state%qg%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqgbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(7)*&
            (current_state%qAitkenSolMass%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k) + &
            current_state%qAitkenSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(8)*&
            (current_state%qAccumSolMass%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k) + &
            current_state%qAccumSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqAccumSolMassbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(9)*&
            (current_state%qAccumInsolMass%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k) + &
            current_state%qAccumInsolMass%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(10)*&
            (current_state%qCoarseSolMass%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k) + &
            current_state%qCoarseSolMass%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(11)*&
            (current_state%qCoarseDustMass%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k) + &
            current_state%qCoarseDustMass%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(12)*&
            (current_state%nl%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olnlbar(k) + &
            current_state%nl%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olnlbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(13)*&
            (current_state%nr%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olnrbar(k) + &
            current_state%nr%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olnrbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(14)*&
            (current_state%ni%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olnibar(k) + &
            current_state%ni%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olnibar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(15)*&
            (current_state%ns%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olnsbar(k) + &
            current_state%ns%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olnsbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(16)*&
            (current_state%ng%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olngbar(k) + &
            current_state%ng%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olngbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(17)*&
            (current_state%nAitkenSolNumber%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k) + &
            current_state%nAitkenSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(18)*&
            (current_state%nAccumSolNumber%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k) + &
            current_state%nAccumSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olnAccumSolNumberbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(19)*&
            (current_state%nAccumInsolNumber%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k) + &
            current_state%nAccumInsolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olnAccumInsolNumberbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(20)*&
            (current_state%nCoarseSolNumber%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k) + &
            current_state%nCoarseSolNumber%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olnCoarseSolNumberbar(k))

          sg_w_buoyancy(k) = sg_w_buoyancy(k) + &
            G*0.5_DEFAULT_PRECISION*current_state%cq(21)*&
            (current_state%nCoarseDustnumber%data(k, current_state%column_local_y, current_state%column_local_x) -&
            current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k) + &
            current_state%nCoarseDustnumber%data(k+1, current_state%column_local_y, current_state%column_local_x)- &
            current_state%global_grid%configuration%vertical%olnCoarseDustnumberbar(k))
        end do
      end if
    end if

    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (allocated(current_state%sub_buoy)) then
        current_state%sub_buoy(k) = current_state%sub_buoy(k) + &
        current_state%w%data(k,current_state%column_local_y,current_state%column_local_x) * &
        sg_w_buoyancy(k)
      end if
    end do

#endif

    current_state%sub_buoy(1) = 0.0_DEFAULT_PRECISION
    current_state%sub_buoy(current_state%local_grid%size(Z_INDEX)) = 0.0_DEFAULT_PRECISION

  end subroutine compute_tke_for_column


  !> Computes the scalar diagnostics for a specific column.
  !! @param current_state The current model state
  subroutine compute_scalars_for_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: k

    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (current_state%w%data(k, current_state%column_local_y,current_state%column_local_x) .gt. current_state%wmfcrit) then
        current_state%mflux = current_state%mflux + current_state%global_grid%configuration%vertical%rho(k)*&
             current_state%global_grid%configuration%vertical%dzn(k)*&
             current_state%w%data(k, current_state%column_local_y,current_state%column_local_x)
      end if
    end do
  end subroutine compute_scalars_for_column



!   !> Computes the qt diagnostics for a specific column. For now we are assuming qt is the same as theta,
!   !! which needs updating with moisture information as per resdgs
!   !! @param current_state The current model state
!   subroutine compute_qt_for_column(current_state)
!     type(model_state_type), target, intent(inout) :: current_state
!
!     real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1, qtpr, qtprp1
!     type(component_field_value_type) :: u_advection, u_viscosity, th_advection, th_diffusion, v_advection, v_viscosity, &
!          w_advection, w_viscosity, w_buoyancy
!     integer :: k
!
!     if (is_component_field_available("u_advection")) u_advection=get_component_field_value(current_state, "u_advection")
!     if (is_component_field_available("u_viscosity")) u_viscosity=get_component_field_value(current_state, "u_viscosity")
!     if (is_component_field_available("th_advection")) th_advection=get_component_field_value(current_state, "th_advection")
!     if (is_component_field_available("th_diffusion")) th_diffusion=get_component_field_value(current_state, "th_diffusion")
!     if (is_component_field_available("v_advection")) v_advection=get_component_field_value(current_state, "v_advection")
!     if (is_component_field_available("v_viscosity")) v_viscosity=get_component_field_value(current_state, "v_viscosity")
!     if (is_component_field_available("w_advection")) w_advection=get_component_field_value(current_state, "w_advection")
!     if (is_component_field_available("w_viscosity")) w_viscosity=get_component_field_value(current_state, "w_viscosity")
!     if (is_component_field_available("w_buoyancy")) w_buoyancy=get_component_field_value(current_state, "w_buoyancy")
!
!     do k=1, current_state%local_grid%size(Z_INDEX)
!       upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
!       uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
!       if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
!         upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
!         uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
!       end if
!       vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
!       vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
!       if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
!         vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
!         vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
!       end if
!       qtpr(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x)
!       qtprp1(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x+1)
!       if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
!         qtpr(k)=qtpr(k)-current_state%global_grid%configuration%vertical%olthbar(k)
!         qtprp1(k)=qtprp1(k)-current_state%global_grid%configuration%vertical%olthbar(k)
!       end if
!     end do
!     do k=2, current_state%local_grid%size(Z_INDEX)-1
!       if (allocated(current_state%us_qt)) current_state%us_qt(k) = current_state%us_qt(k)+0.5*(upr(k)+uprm1(k))*&
!            current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(qtpr(k)+qtprp1(k))*&
!            current_state%su%data(k,current_state%column_local_y,current_state%column_local_x)
!       if (allocated(current_state%u_qt_advection)) current_state%u_qt_advection(k) = current_state%u_qt_advection(k)+0.5*(upr(k)+uprm1(k))*&
!            th_advection%real_1d_array(k)+0.5*(qtpr(k)+qtprp1(k))*u_advection%real_1d_array(k)
!       if (allocated(current_state%u_qt_viscosity_diffusion)) current_state%u_qt_viscosity_diffusion(k) = current_state%u_qt_viscosity_diffusion(k)+0.5*&
!            (qtpr(k)+qtprp1(k))*u_viscosity%real_1d_array(k)+0.5*(upr(k)+uprm1(k))*th_diffusion%real_1d_array(k)
!       if (allocated(current_state%wu_qt)) current_state%wu_qt(k) = current_state%wu_qt(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(upr(k+1)+upr(k)+&
!            uprm1(k+1)+uprm1(k))*0.5*(qtpr(k+1)+qtpr(k))
!
!       if (allocated(current_state%vs_qt)) current_state%vs_qt(k) = current_state%vs_qt(k)+0.5*(vpr(k)+vprm1(k))*&
!            current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(qtpr(k)+qtprp1(k))*&
!            current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x)
!       if (allocated(current_state%v_qt_advection)) current_state%v_qt_advection(k) = current_state%v_qt_advection(k)+0.5*(vpr(k)+vprm1(k))*&
!            th_advection%real_1d_array(k)+0.5*(qtpr(k)+qtprp1(k))*v_advection%real_1d_array(k)
!       if (allocated(current_state%v_qt_viscosity_diffusion)) current_state%v_qt_viscosity_diffusion(k) = current_state%v_qt_viscosity_diffusion(k)+0.5*&
!            (qtpr(k)+qtprp1(k))*v_viscosity%real_1d_array(k)+0.5*(vpr(k)+vprm1(k))*th_diffusion%real_1d_array(k)
!       if (allocated(current_state%wv_qt)) current_state%wv_qt(k) = current_state%wv_qt(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(vpr(k+1)+vpr(k)+&
!            vprm1(k+1)+vprm1(k))*0.5*(qtpr(k+1)+qtpr(k))
!
!       if (allocated(current_state%w_qt)) current_state%w_qt(k) = current_state%w_qt(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qtpr(k)+qtpr(k+1))
!       if (allocated(current_state%ws_qt)) current_state%ws_qt(k) = current_state%ws_qt(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
!            (current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+&
!            current_state%sth%data(k+1,current_state%column_local_y,current_state%column_local_x))+0.5*(qtpr(k)+qtpr(k+1))*&
!            current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)
!       if (allocated(current_state%w_qt_advection)) current_state%w_qt_advection(k) = current_state%w_qt_advection(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(th_advection%real_1d_array(k)+&
!            th_advection%real_1d_array(k+1))+0.5*(qtpr(k)+qtpr(k+1))*w_advection%real_1d_array(k)
!       if (allocated(current_state%w_qt_viscosity_diffusion)) current_state%w_qt_viscosity_diffusion(k) = current_state%w_qt_viscosity_diffusion(k)+0.5*&
!            (qtpr(k)+qtpr(k+1))*w_viscosity%real_1d_array(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
!            (th_diffusion%real_1d_array(k)+th_diffusion%real_1d_array(k+1))
!       if (allocated(current_state%w_qt_buoyancy)) current_state%w_qt_buoyancy(k) = current_state%w_qt_buoyancy(k)+0.5*(qtpr(k+1)+qtpr(k))*&
!            w_buoyancy%real_1d_array(k)
!       if (allocated(current_state%ww_qt)) current_state%ww_qt(k) = current_state%ww_qt(k)+0.5*&
!            (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
!            current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
!            current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*qtpr(k)
!
!       if (allocated(current_state%qt_qt)) current_state%qt_qt(k) = current_state%qt_qt(k)+qtpr(k)*qtpr(k)
!       if (allocated(current_state%sqt_qt)) current_state%sqt_qt(k) = current_state%sqt_qt(k)+2.0*qtpr(k)*&
!            current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)
!       if (allocated(current_state%qt_qt_advection)) current_state%qt_qt_advection(k) = current_state%qt_qt_advection(k)+2.0*qtpr(k)*&
!            th_advection%real_1d_array(k)
!       if (allocated(current_state%qt_qt_diffusion)) current_state%qt_qt_diffusion(k) = current_state%qt_qt_diffusion(k)+2.0*qtpr(k)*&
!            th_diffusion%real_1d_array(k)
!       if (allocated(current_state%wqt_qt)) current_state%wqt_qt(k) = current_state%wqt_qt(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(qtpr(k+1)+qtpr(k))*0.5*(&
!            qtpr(k+1)+qtpr(k))
!     end do
!
!     if (allocated(u_advection%real_1d_array)) deallocate(u_advection%real_1d_array)
!     if (allocated(u_viscosity%real_1d_array)) deallocate(u_viscosity%real_1d_array)
!     if (allocated(th_advection%real_1d_array)) deallocate(th_advection%real_1d_array)
!     if (allocated(th_diffusion%real_1d_array)) deallocate(th_diffusion%real_1d_array)
!     if (allocated(v_advection%real_1d_array)) deallocate(v_advection%real_1d_array)
!     if (allocated(v_viscosity%real_1d_array)) deallocate(v_viscosity%real_1d_array)
!     if (allocated(w_advection%real_1d_array)) deallocate(w_advection%real_1d_array)
!     if (allocated(w_viscosity%real_1d_array)) deallocate(w_viscosity%real_1d_array)
!     if (allocated(w_buoyancy%real_1d_array)) deallocate(w_buoyancy%real_1d_array)
!   end subroutine compute_qt_for_column
!
!
!
!
!   !> Computes the mse diagnostics for a specific column. For now we are assuming mse is the same as theta,
!   !! which needs updating with moisture information
!   !! @param current_state The current model state
!   subroutine compute_mse_for_column(current_state)
!     type(model_state_type), target, intent(inout) :: current_state
!
!     real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: upr, vpr, uprm1, vprm1, msepr, mseprp1
!     type(component_field_value_type) :: u_advection, u_viscosity, th_advection, th_diffusion, v_advection, v_viscosity, &
!          w_advection, w_viscosity, w_buoyancy
!     integer :: k
!
!     if (is_component_field_available("u_advection")) u_advection=get_component_field_value(current_state, "u_advection")
!     if (is_component_field_available("u_viscosity")) u_viscosity=get_component_field_value(current_state, "u_viscosity")
!     if (is_component_field_available("th_advection")) th_advection=get_component_field_value(current_state, "th_advection")
!     if (is_component_field_available("th_diffusion")) th_diffusion=get_component_field_value(current_state, "th_diffusion")
!     if (is_component_field_available("v_advection")) v_advection=get_component_field_value(current_state, "v_advection")
!     if (is_component_field_available("v_viscosity")) v_viscosity=get_component_field_value(current_state, "v_viscosity")
!     if (is_component_field_available("w_advection")) w_advection=get_component_field_value(current_state, "w_advection")
!     if (is_component_field_available("w_viscosity")) w_viscosity=get_component_field_value(current_state, "w_viscosity")
!     if (is_component_field_available("w_buoyancy")) w_buoyancy=get_component_field_value(current_state, "w_buoyancy")
!
!     do k=1, current_state%local_grid%size(Z_INDEX)
!       upr(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x)
!       uprm1(k)=current_state%u%data(k,current_state%column_local_y,current_state%column_local_x-1)
!       if (allocated(current_state%global_grid%configuration%vertical%olubar)) then
!         upr(k)=upr(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
!         uprm1(k)=uprm1(k)-(current_state%global_grid%configuration%vertical%olubar(k)-current_state%ugal)
!       end if
!       vpr(k)=current_state%v%data(k,current_state%column_local_y,current_state%column_local_x)
!       vprm1(k)=current_state%v%data(k,current_state%column_local_y-1,current_state%column_local_x)
!       if (allocated(current_state%global_grid%configuration%vertical%olvbar)) then
!         vpr(k)=vpr(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
!         vprm1(k)=vprm1(k)-(current_state%global_grid%configuration%vertical%olvbar(k)-current_state%vgal)
!       end if
!       msepr(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x)
!       mseprp1(k)=current_state%th%data(k,current_state%column_local_y,current_state%column_local_x+1)
!       if (allocated(current_state%global_grid%configuration%vertical%olthbar)) then
!         msepr(k)=msepr(k)-current_state%global_grid%configuration%vertical%olthbar(k)
!         mseprp1(k)=mseprp1(k)-current_state%global_grid%configuration%vertical%olthbar(k)
!       end if
!     end do
!     do k=2, current_state%local_grid%size(Z_INDEX)-1
!       if (allocated(current_state%u_mse)) current_state%u_mse(k) = current_state%u_mse(k)+0.5*(upr(k)+uprm1(k))*msepr(k)
!       if (allocated(current_state%us_mse)) current_state%us_mse(k) = current_state%us_mse(k)+0.5*(upr(k)+uprm1(k))*&
!            current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(msepr(k)+mseprp1(k))*&
!            current_state%su%data(k,current_state%column_local_y,current_state%column_local_x)
!       if (allocated(current_state%u_mse_advection)) current_state%u_mse_advection(k) = current_state%u_mse_advection(k)+0.5*(upr(k)+uprm1(k))*&
!            th_advection%real_1d_array(k)+0.5*(msepr(k)+mseprp1(k))*u_advection%real_1d_array(k)
!       if (allocated(current_state%u_mse_viscosity_diffusion)) current_state%u_mse_viscosity_diffusion(k) = current_state%u_mse_viscosity_diffusion(k)+0.5*&
!            (msepr(k)+mseprp1(k))*u_viscosity%real_1d_array(k)+0.5*(upr(k)+uprm1(k))*th_diffusion%real_1d_array(k)
!       if (allocated(current_state%wu_mse)) current_state%wu_mse(k) = current_state%wu_mse(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(upr(k+1)+upr(k)+&
!            uprm1(k+1)+uprm1(k))*0.5*(msepr(k+1)+msepr(k))
!
!       if (allocated(current_state%v_mse)) current_state%v_mse(k) = current_state%v_mse(k)+0.5*(vpr(k)+vprm1(k))*msepr(k)
!       if (allocated(current_state%vs_mse)) current_state%vs_mse(k) = current_state%vs_mse(k)+0.5*(vpr(k)+vprm1(k))*&
!            current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+0.5*(msepr(k)+mseprp1(k))*&
!            current_state%sv%data(k,current_state%column_local_y,current_state%column_local_x)
!       if (allocated(current_state%v_mse_advection)) current_state%v_mse_advection(k) = current_state%v_mse_advection(k)+0.5*(vpr(k)+vprm1(k))*&
!            th_advection%real_1d_array(k)+0.5*(msepr(k)+mseprp1(k))*v_advection%real_1d_array(k)
!       if (allocated(current_state%v_mse_viscosity_diffusion)) current_state%v_mse_viscosity_diffusion(k) = current_state%v_mse_viscosity_diffusion(k)+0.5*&
!            (msepr(k)+mseprp1(k))*v_viscosity%real_1d_array(k)+0.5*(vpr(k)+vprm1(k))*th_diffusion%real_1d_array(k)
!       if (allocated(current_state%wv_mse)) current_state%wv_mse(k) = current_state%wv_mse(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.25*(vpr(k+1)+vpr(k)+&
!            vprm1(k+1)+vprm1(k))*0.5*(msepr(k+1)+msepr(k))
!
!       if (allocated(current_state%w_mse)) current_state%w_mse(k) = current_state%w_mse(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(msepr(k)+msepr(k+1))
!       if (allocated(current_state%ws_mse)) current_state%ws_mse(k) = current_state%ws_mse(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
!            (current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)+&
!            current_state%sth%data(k+1,current_state%column_local_y,current_state%column_local_x))+0.5*(msepr(k)+msepr(k+1))*&
!            current_state%sw%data(k,current_state%column_local_y,current_state%column_local_x)
!       if (allocated(current_state%w_mse_advection)) current_state%w_mse_advection(k) = current_state%w_mse_advection(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(th_advection%real_1d_array(k)+&
!            th_advection%real_1d_array(k+1))+0.5*(msepr(k)+msepr(k+1))*w_advection%real_1d_array(k)
!       if (allocated(current_state%w_mse_viscosity_diffusion)) current_state%w_mse_viscosity_diffusion(k) = current_state%w_mse_viscosity_diffusion(k)+0.5*&
!            (msepr(k)+msepr(k+1))*w_viscosity%real_1d_array(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*&
!            (th_diffusion%real_1d_array(k)+th_diffusion%real_1d_array(k+1))
!       if (allocated(current_state%w_mse_buoyancy)) current_state%w_mse_buoyancy(k) = current_state%w_mse_buoyancy(k)+0.5*(msepr(k+1)+msepr(k))*&
!            w_buoyancy%real_1d_array(k)
!       if (allocated(current_state%ww_mse)) current_state%ww_mse(k) = current_state%ww_mse(k)+0.5*&
!            (current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
!            current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*0.5*(&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)+&
!            current_state%w%data(k-1,current_state%column_local_y,current_state%column_local_x))*msepr(k)
!
!       if (allocated(current_state%mse_mse)) current_state%mse_mse(k) = current_state%mse_mse(k)+msepr(k)*msepr(k)
!       if (allocated(current_state%smse_mse)) current_state%smse_mse(k) = current_state%smse_mse(k)+2.0*msepr(k)*&
!            current_state%sth%data(k,current_state%column_local_y,current_state%column_local_x)
!       if (allocated(current_state%mse_mse_advection)) current_state%mse_mse_advection(k) = current_state%mse_mse_advection(k)+2.0*msepr(k)*&
!            th_advection%real_1d_array(k)
!       if (allocated(current_state%mse_mse_diffusion)) current_state%mse_mse_diffusion(k) = current_state%mse_mse_diffusion(k)+2.0*msepr(k)*&
!            th_diffusion%real_1d_array(k)
!       if (allocated(current_state%wmse_mse)) current_state%wmse_mse(k) = current_state%wmse_mse(k)+&
!            current_state%w%data(k,current_state%column_local_y,current_state%column_local_x)*0.5*(msepr(k+1)+msepr(k))*0.5*(&
!            msepr(k+1)+msepr(k))
!     end do
!
!     if (allocated(u_advection%real_1d_array)) deallocate(u_advection%real_1d_array)
!     if (allocated(u_viscosity%real_1d_array)) deallocate(u_viscosity%real_1d_array)
!     if (allocated(th_advection%real_1d_array)) deallocate(th_advection%real_1d_array)
!     if (allocated(th_diffusion%real_1d_array)) deallocate(th_diffusion%real_1d_array)
!     if (allocated(v_advection%real_1d_array)) deallocate(v_advection%real_1d_array)
!     if (allocated(v_viscosity%real_1d_array)) deallocate(v_viscosity%real_1d_array)
!     if (allocated(w_advection%real_1d_array)) deallocate(w_advection%real_1d_array)
!     if (allocated(w_viscosity%real_1d_array)) deallocate(w_viscosity%real_1d_array)
!     if (allocated(w_buoyancy%real_1d_array)) deallocate(w_buoyancy%real_1d_array)
!   end subroutine compute_mse_for_column



   !> Finalisation call back
  !! @param current_state Current model state
  subroutine finalisation_callback_flux_budget(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(current_state%th_flux_values)) deallocate(current_state%th_flux_values)
    if (allocated(current_state%th_gradient)) deallocate(current_state%th_gradient)
    if (allocated(current_state%th_diff)) deallocate(current_state%th_diff)
    if (allocated(current_state%th_buoyancy)) deallocate(current_state%th_buoyancy)
    if (allocated(current_state%th_tendency)) deallocate(current_state%th_tendency)


    if (allocated(current_state%qv_flux_values)) deallocate(current_state%qv_flux_values)
    if (allocated(current_state%ql_flux_values)) deallocate(current_state%ql_flux_values)
    if (allocated(current_state%qr_flux_values)) deallocate(current_state%qr_flux_values)
    if (allocated(current_state%qi_flux_values)) deallocate(current_state%qi_flux_values)
    if (allocated(current_state%qs_flux_values)) deallocate(current_state%qs_flux_values)
    if (allocated(current_state%qg_flux_values)) deallocate(current_state%qg_flux_values)
    if (allocated(current_state%qAitkenSolMass_flux_values)) deallocate(current_state%qAitkenSolMass_flux_values)
    if (allocated(current_state%qAccumSolMass_flux_values)) deallocate(current_state%qAccumSolMass_flux_values)
    if (allocated(current_state%qAccumInsolMass_flux_values)) deallocate(current_state%qAccumInsolMass_flux_values)
    if (allocated(current_state%qCoarseSolMass_flux_values)) deallocate(current_state%qCoarseSolMass_flux_values)
    if (allocated(current_state%qCoarseDustMass_flux_values)) deallocate(current_state%qCoarseDustMass_flux_values)
    if (allocated(current_state%nl_flux_values)) deallocate(current_state%nl_flux_values)
    if (allocated(current_state%nr_flux_values)) deallocate(current_state%nr_flux_values)
    if (allocated(current_state%ni_flux_values)) deallocate(current_state%ni_flux_values)
    if (allocated(current_state%ns_flux_values)) deallocate(current_state%ns_flux_values)
    if (allocated(current_state%ng_flux_values)) deallocate(current_state%ng_flux_values)
    if (allocated(current_state%nAitkenSolNumber_flux_values)) deallocate(current_state%nAitkenSolNumber_flux_values)
    if (allocated(current_state%nAccumSolNumber_flux_values)) deallocate(current_state%nAccumSolNumber_flux_values)
    if (allocated(current_state%nAccumInsolNumber_flux_values)) deallocate(current_state%nAccumInsolNumber_flux_values)
    if (allocated(current_state%nCoarseSolNumber_flux_values)) deallocate(current_state%nCoarseSolNumber_flux_values)
    if (allocated(current_state%nCoarseDustnumber_flux_values)) deallocate(current_state%nCoarseDustnumber_flux_values)

    if (allocated(current_state%qv_gradient)) deallocate(current_state%qv_gradient)
    if (allocated(current_state%ql_gradient)) deallocate(current_state%ql_gradient)
    if (allocated(current_state%qr_gradient)) deallocate(current_state%qr_gradient)
    if (allocated(current_state%qi_gradient)) deallocate(current_state%qi_gradient)
    if (allocated(current_state%qs_gradient)) deallocate(current_state%qs_gradient)
    if (allocated(current_state%qg_gradient)) deallocate(current_state%qg_gradient)
    if (allocated(current_state%qAitkenSolMass_gradient)) deallocate(current_state%qAitkenSolMass_gradient)
    if (allocated(current_state%qAccumSolMass_gradient)) deallocate(current_state%qAccumSolMass_gradient)
    if (allocated(current_state%qAccumInsolMass_gradient)) deallocate(current_state%qAccumInsolMass_gradient)
    if (allocated(current_state%qCoarseSolMass_gradient)) deallocate(current_state%qCoarseSolMass_gradient)
    if (allocated(current_state%qCoarseDustMass_gradient)) deallocate(current_state%qCoarseDustMass_gradient)
    if (allocated(current_state%nl_gradient)) deallocate(current_state%nl_gradient)
    if (allocated(current_state%nr_gradient)) deallocate(current_state%nr_gradient)
    if (allocated(current_state%ni_gradient)) deallocate(current_state%ni_gradient)
    if (allocated(current_state%ns_gradient)) deallocate(current_state%ns_gradient)
    if (allocated(current_state%ng_gradient)) deallocate(current_state%ng_gradient)
    if (allocated(current_state%nAitkenSolNumber_gradient)) deallocate(current_state%nAitkenSolNumber_gradient)
    if (allocated(current_state%nAccumSolNumber_gradient)) deallocate(current_state%nAccumSolNumber_gradient)
    if (allocated(current_state%nAccumInsolNumber_gradient)) deallocate(current_state%nAccumInsolNumber_gradient)
    if (allocated(current_state%nCoarseSolNumber_gradient)) deallocate(current_state%nCoarseSolNumber_gradient)
    if (allocated(current_state%nCoarseDustnumber_gradient)) deallocate(current_state%nCoarseDustnumber_gradient)

    if (allocated(current_state%qv_diff)) deallocate(current_state%qv_diff)
    if (allocated(current_state%ql_diff)) deallocate(current_state%ql_diff)
    if (allocated(current_state%qr_diff)) deallocate(current_state%qr_diff)
    if (allocated(current_state%qi_diff)) deallocate(current_state%qi_diff)
    if (allocated(current_state%qs_diff)) deallocate(current_state%qs_diff)
    if (allocated(current_state%qg_diff)) deallocate(current_state%qg_diff)
    if (allocated(current_state%qAitkenSolMass_diff)) deallocate(current_state%qAitkenSolMass_diff)
    if (allocated(current_state%qAccumSolMass_diff)) deallocate(current_state%qAccumSolMass_diff)
    if (allocated(current_state%qAccumInsolMass_diff)) deallocate(current_state%qAccumInsolMass_diff)
    if (allocated(current_state%qCoarseSolMass_diff)) deallocate(current_state%qCoarseSolMass_diff)
    if (allocated(current_state%qCoarseDustMass_diff)) deallocate(current_state%qCoarseDustMass_diff)
    if (allocated(current_state%nl_diff)) deallocate(current_state%nl_diff)
    if (allocated(current_state%nr_diff)) deallocate(current_state%nr_diff)
    if (allocated(current_state%ni_diff)) deallocate(current_state%ni_diff)
    if (allocated(current_state%ns_diff)) deallocate(current_state%ns_diff)
    if (allocated(current_state%ng_diff)) deallocate(current_state%ng_diff)
    if (allocated(current_state%nAitkenSolNumber_diff)) deallocate(current_state%nAitkenSolNumber_diff)
    if (allocated(current_state%nAccumSolNumber_diff)) deallocate(current_state%nAccumSolNumber_diff)
    if (allocated(current_state%nAccumInsolNumber_diff)) deallocate(current_state%nAccumInsolNumber_diff)
    if (allocated(current_state%nCoarseSolNumber_diff)) deallocate(current_state%nCoarseSolNumber_diff)
    if (allocated(current_state%nCoarseDustnumber_diff)) deallocate(current_state%nCoarseDustnumber_diff)

    if (allocated(current_state%qv_buoyancy)) deallocate(current_state%qv_buoyancy)
    if (allocated(current_state%ql_buoyancy)) deallocate(current_state%ql_buoyancy)
    if (allocated(current_state%qr_buoyancy)) deallocate(current_state%qr_buoyancy)
    if (allocated(current_state%qi_buoyancy)) deallocate(current_state%qi_buoyancy)
    if (allocated(current_state%qs_buoyancy)) deallocate(current_state%qs_buoyancy)
    if (allocated(current_state%qg_buoyancy)) deallocate(current_state%qg_buoyancy)
    if (allocated(current_state%qAitkenSolMass_buoyancy)) deallocate(current_state%qAitkenSolMass_buoyancy)
    if (allocated(current_state%qAccumSolMass_buoyancy)) deallocate(current_state%qAccumSolMass_buoyancy)
    if (allocated(current_state%qAccumInsolMass_buoyancy)) deallocate(current_state%qAccumInsolMass_buoyancy)
    if (allocated(current_state%qCoarseSolMass_buoyancy)) deallocate(current_state%qCoarseSolMass_buoyancy)
    if (allocated(current_state%qCoarseDustMass_buoyancy)) deallocate(current_state%qCoarseDustMass_buoyancy)
    if (allocated(current_state%nl_buoyancy)) deallocate(current_state%nl_buoyancy)
    if (allocated(current_state%nr_buoyancy)) deallocate(current_state%nr_buoyancy)
    if (allocated(current_state%ni_buoyancy)) deallocate(current_state%ni_buoyancy)
    if (allocated(current_state%ns_buoyancy)) deallocate(current_state%ns_buoyancy)
    if (allocated(current_state%ng_buoyancy)) deallocate(current_state%ng_buoyancy)
    if (allocated(current_state%nAitkenSolNumber_buoyancy)) deallocate(current_state%nAitkenSolNumber_buoyancy)
    if (allocated(current_state%nAccumSolNumber_buoyancy)) deallocate(current_state%nAccumSolNumber_buoyancy)
    if (allocated(current_state%nAccumInsolNumber_buoyancy)) deallocate(current_state%nAccumInsolNumber_buoyancy)
    if (allocated(current_state%nCoarseSolNumber_buoyancy)) deallocate(current_state%nCoarseSolNumber_buoyancy)
    if (allocated(current_state%nCoarseDustnumber_buoyancy)) deallocate(current_state%nCoarseDustnumber_buoyancy)

    if (allocated(current_state%qv_tendency)) deallocate(current_state%qv_tendency)
    if (allocated(current_state%ql_tendency)) deallocate(current_state%ql_tendency)
    if (allocated(current_state%qr_tendency)) deallocate(current_state%qr_tendency)
    if (allocated(current_state%qi_tendency)) deallocate(current_state%qi_tendency)
    if (allocated(current_state%qs_tendency)) deallocate(current_state%qs_tendency)
    if (allocated(current_state%qg_tendency)) deallocate(current_state%qg_tendency)
    if (allocated(current_state%qAitkenSolMass_tendency)) deallocate(current_state%qAitkenSolMass_tendency)
    if (allocated(current_state%qAccumSolMass_tendency)) deallocate(current_state%qAccumSolMass_tendency)
    if (allocated(current_state%qAccumInsolMass_tendency)) deallocate(current_state%qAccumInsolMass_tendency)
    if (allocated(current_state%qCoarseSolMass_tendency)) deallocate(current_state%qCoarseSolMass_tendency)
    if (allocated(current_state%qCoarseDustMass_tendency)) deallocate(current_state%qCoarseDustMass_tendency)
    if (allocated(current_state%nl_tendency)) deallocate(current_state%nl_tendency)
    if (allocated(current_state%nr_tendency)) deallocate(current_state%nr_tendency)
    if (allocated(current_state%ni_tendency)) deallocate(current_state%ni_tendency)
    if (allocated(current_state%ns_tendency)) deallocate(current_state%ns_tendency)
    if (allocated(current_state%ng_tendency)) deallocate(current_state%ng_tendency)
    if (allocated(current_state%nAitkenSolNumber_tendency)) deallocate(current_state%nAitkenSolNumber_tendency)
    if (allocated(current_state%nAccumSolNumber_tendency)) deallocate(current_state%nAccumSolNumber_tendency)
    if (allocated(current_state%nAccumInsolNumber_tendency)) deallocate(current_state%nAccumInsolNumber_tendency)
    if (allocated(current_state%nCoarseSolNumber_tendency)) deallocate(current_state%nCoarseSolNumber_tendency)
    if (allocated(current_state%nCoarseDustnumber_tendency)) deallocate(current_state%nCoarseDustnumber_tendency)


    if (allocated(current_state%uw_advection)) deallocate(current_state%uw_advection)
    if (allocated(current_state%vw_advection)) deallocate(current_state%vw_advection)
    if (allocated(current_state%uw_viscosity)) deallocate(current_state%uw_viscosity)
    if (allocated(current_state%vw_viscosity)) deallocate(current_state%vw_viscosity)
    if (allocated(current_state%uw_buoyancy)) deallocate(current_state%uw_buoyancy)
    if (allocated(current_state%vw_buoyancy)) deallocate(current_state%vw_buoyancy)
    if (allocated(current_state%uw_tendency)) deallocate(current_state%uw_tendency)
    if (allocated(current_state%vw_tendency)) deallocate(current_state%vw_tendency)
    if (allocated(current_state%uw_w)) deallocate(current_state%uw_w)
    if (allocated(current_state%vw_w)) deallocate(current_state%vw_w)

    if (allocated(current_state%tu_su)) deallocate(current_state%tu_su)
    if (allocated(current_state%uu_advection)) deallocate(current_state%uu_advection)
    if (allocated(current_state%uu_viscosity)) deallocate(current_state%uu_viscosity)
    if (allocated(current_state%wu_u)) deallocate(current_state%wu_u)
    if (allocated(current_state%tv_sv)) deallocate(current_state%tv_sv)
    if (allocated(current_state%vv_advection)) deallocate(current_state%vv_advection)
    if (allocated(current_state%vv_viscosity)) deallocate(current_state%vv_viscosity)
    if (allocated(current_state%wv_v)) deallocate(current_state%wv_v)
    if (allocated(current_state%tw_sw)) deallocate(current_state%tw_sw)
    if (allocated(current_state%ww_advection)) deallocate(current_state%ww_advection)
    if (allocated(current_state%ww_viscosity)) deallocate(current_state%ww_viscosity)
    if (allocated(current_state%ww_buoyancy)) deallocate(current_state%ww_buoyancy)

    if (allocated(current_state%u_thetal)) deallocate(current_state%u_thetal)
    if (allocated(current_state%us_thetal)) deallocate(current_state%us_thetal)
    if (allocated(current_state%u_thetal_advection)) deallocate(current_state%u_thetal_advection)
    if (allocated(current_state%u_thetal_viscosity_diffusion)) deallocate(current_state%u_thetal_viscosity_diffusion)
    if (allocated(current_state%wu_thetal)) deallocate(current_state%wu_thetal)
    if (allocated(current_state%v_thetal)) deallocate(current_state%v_thetal)
    if (allocated(current_state%vs_thetal)) deallocate(current_state%vs_thetal)
    if (allocated(current_state%v_thetal_advection)) deallocate(current_state%v_thetal_advection)
    if (allocated(current_state%v_thetal_viscosity_diffusion)) deallocate(current_state%v_thetal_viscosity_diffusion)
    if (allocated(current_state%wv_thetal)) deallocate(current_state%wv_thetal)
    if (allocated(current_state%w_thetal)) deallocate(current_state%w_thetal)
    if (allocated(current_state%ws_thetal)) deallocate(current_state%ws_thetal)
    if (allocated(current_state%w_thetal_advection)) deallocate(current_state%w_thetal_advection)
    if (allocated(current_state%w_thetal_viscosity_diffusion)) deallocate(current_state%w_thetal_viscosity_diffusion)
    if (allocated(current_state%w_thetal_buoyancy)) deallocate(current_state%w_thetal_buoyancy)
    if (allocated(current_state%ww_thetal)) deallocate(current_state%ww_thetal)
    if (allocated(current_state%thetal_thetal)) deallocate(current_state%thetal_thetal)
    if (allocated(current_state%sthetal_thetal)) deallocate(current_state%sthetal_thetal)
    if (allocated(current_state%thetal_thetal_advection)) deallocate(current_state%thetal_thetal_advection)
    if (allocated(current_state%thetal_thetal_diffusion)) deallocate(current_state%thetal_thetal_diffusion)
    if (allocated(current_state%wthetal_thetal)) deallocate(current_state%wthetal_thetal)

    if (allocated(current_state%shprd)) deallocate(current_state%shprd)
    if (allocated(current_state%sub_buoy)) deallocate(current_state%sub_buoy)
    if (allocated(current_state%w_ke)) deallocate(current_state%w_ke)
    if (allocated(current_state%w_p)) deallocate(current_state%w_p)
    if (allocated(current_state%tke_tend)) deallocate(current_state%tke_tend)

!     if (allocated(current_state%u_mse)) deallocate(current_state%u_mse)
!     if (allocated(current_state%us_mse)) deallocate(current_state%us_mse)
!     if (allocated(current_state%u_mse_advection)) deallocate(current_state%u_mse_advection)
!     if (allocated(current_state%u_mse_viscosity_diffusion)) deallocate(current_state%u_mse_viscosity_diffusion)
!     if (allocated(current_state%wu_mse)) deallocate(current_state%wu_mse)
!     if (allocated(current_state%v_mse)) deallocate(current_state%v_mse)
!     if (allocated(current_state%vs_mse)) deallocate(current_state%vs_mse)
!     if (allocated(current_state%v_mse_advection)) deallocate(current_state%v_mse_advection)
!     if (allocated(current_state%v_mse_viscosity_diffusion)) deallocate(current_state%v_mse_viscosity_diffusion)
!     if (allocated(current_state%wv_mse)) deallocate(current_state%wv_mse)
!     if (allocated(current_state%w_mse)) deallocate(current_state%w_mse)
!     if (allocated(current_state%ws_mse)) deallocate(current_state%ws_mse)
!     if (allocated(current_state%w_mse_advection)) deallocate(current_state%w_mse_advection)
!     if (allocated(current_state%w_mse_viscosity_diffusion)) deallocate(current_state%w_mse_viscosity_diffusion)
!     if (allocated(current_state%w_mse_buoyancy)) deallocate(current_state%w_mse_buoyancy)
!     if (allocated(current_state%ww_mse)) deallocate(current_state%ww_mse)
!     if (allocated(current_state%mse_mse)) deallocate(current_state%mse_mse)
!     if (allocated(current_state%smse_mse)) deallocate(current_state%smse_mse)
!     if (allocated(current_state%mse_mse_advection)) deallocate(current_state%mse_mse_advection)
!     if (allocated(current_state%mse_mse_diffusion)) deallocate(current_state%mse_mse_diffusion)
!     if (allocated(current_state%wmse_mse)) deallocate(current_state%wmse_mse)
!
!     if (allocated(current_state%us_qt)) deallocate(current_state%us_qt)
!     if (allocated(current_state%u_qt_advection)) deallocate(current_state%u_qt_advection)
!     if (allocated(current_state%u_qt_viscosity_diffusion)) deallocate(current_state%u_qt_viscosity_diffusion)
!     if (allocated(current_state%wu_qt)) deallocate(current_state%wu_qt)
!     if (allocated(current_state%vs_qt)) deallocate(current_state%vs_qt)
!     if (allocated(current_state%v_qt_advection)) deallocate(current_state%v_qt_advection)
!     if (allocated(current_state%v_qt_viscosity_diffusion)) deallocate(current_state%v_qt_viscosity_diffusion)
!     if (allocated(current_state%wv_qt)) deallocate(current_state%wv_qt)
!     if (allocated(current_state%w_qt)) deallocate(current_state%w_qt)
!     if (allocated(current_state%ws_qt)) deallocate(current_state%ws_qt)
!     if (allocated(current_state%w_qt_advection)) deallocate(current_state%w_qt_advection)
!     if (allocated(current_state%w_qt_viscosity_diffusion)) deallocate(current_state%w_qt_viscosity_diffusion)
!     if (allocated(current_state%w_qt_buoyancy)) deallocate(current_state%w_qt_buoyancy)
!     if (allocated(current_state%ww_qt)) deallocate(current_state%ww_qt)
!     if (allocated(current_state%qt_qt)) deallocate(current_state%qt_qt)
!     if (allocated(current_state%sqt_qt)) deallocate(current_state%sqt_qt)
!     if (allocated(current_state%qt_qt_advection)) deallocate(current_state%qt_qt_advection)
!     if (allocated(current_state%qt_qt_diffusion)) deallocate(current_state%qt_qt_diffusion)
!     if (allocated(current_state%wqt_qt)) deallocate(current_state%wqt_qt)
  end subroutine finalisation_callback_flux_budget

!   !> Determines whether a specific published field is a heat flux field
!   !! @param name The name of the field to check
!   !! @returns Whether the field name is a heat flux field
!   logical function is_field_heat_flux(name)
!     character(len=*), intent(in) :: name
!
!     is_field_heat_flux=c_contains(heat_flux_fields, name)
!   end function is_field_heat_flux
!
!   !> Determines whether a specific published field is a TKE budget field
!   !! @param name The name of the field to check
!   !! @returns Whether the field name is a TKE budget field
!   logical function is_field_tke_flux(name)
!     character(len=*), intent(in) :: name
!
!     is_field_tke_flux=c_contains(tke_fields, name)
!   end function is_field_tke_flux
!
!   !> Determines whether a specific published field is a q flux field
!   !! @param name The name of the field to check
!   !! @returns Whether the field name is a q flux field
!   logical function is_field_q_flux(name)
!     character(len=*), intent(in) :: name
!
!     is_field_q_flux=c_contains(q_flux_fields, name)
!   end function is_field_q_flux
!
!   !> Determines whether a specific published field is a uw or uv field
!   !! @param name The name of the field to check
!   !! @returns Whether the field name is a uw or uv field
!   logical function is_field_uw_vw(name)
!     character(len=*), intent(in) :: name
!
!     is_field_uw_vw=c_contains(uw_vw_fields, name)
!   end function is_field_uw_vw
!
!   !> Determines whether a specific published field is a uu, vv or ww field
!   !! @param name The name of the field to check
!   !! @returns Whether the field name is a uu, vv or ww field
!   logical function is_field_prognostic_budget(name)
!     character(len=*), intent(in) :: name
!
!     is_field_prognostic_budget=c_contains(prognostic_budget_fields, name)
!   end function is_field_prognostic_budget
!
!   !> Determines whether a specific published field is a thetal field
!   !! @param name The name of the field to check
!   !! @returns Whether the field name is a thetal field
!   logical function is_field_thetal(name)
!     character(len=*), intent(in) :: name
!
!     is_field_thetal=c_contains(thetal_fields, name)
!   end function is_field_thetal
!
!   !> Determines whether a specific published field is a mse field
!   !! @param name The name of the field to check
!   !! @returns Whether the field name is a mse field
!   logical function is_field_mse(name)
!     character(len=*), intent(in) :: name
!
!     is_field_mse=c_contains(mse_fields, name)
!   end function is_field_mse
!
!   !> Determines whether a specific published field is a mse field
!   !! @param name The name of the field to check
!   !! @returns Whether the field name is a mse field
!   logical function is_field_qt(name)
!     character(len=*), intent(in) :: name
!
!     is_field_qt=c_contains(qt_fields, name)
!   end function is_field_qt
!
!   !> Determines whether a specific published field is a scalar field
!   !! @param name The name of the field to check
!   !! @returns Whether the field name is a scalar field
!   logical function is_field_scalar(name)
!     character(len=*), intent(in) :: name
!
!     is_field_scalar=c_contains(scalar_fields, name)
!   end function is_field_scalar

  
  !> Field value retrieval callback, this returns the value of a specific published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve the value for
  !! @param field_value Populated with the value of the field
!   subroutine field_value_retrieval_callback(current_state, name, field_value)
!     type(model_state_type), target, intent(inout) :: current_state
!     character(len=*), intent(in) :: name
!     type(component_field_value_type), intent(out) :: field_value
!
!     if (name .eq. "heat_flux_transport_local" .and. allocated(th_flux_values)) then
!       call set_published_field_value(field_value, real_1d_field=th_flux_values)
!     else if (name .eq. "heat_flux_gradient_local" .and. allocated(th_gradient)) then
!       call set_published_field_value(field_value, real_1d_field=th_gradient)
!     else if (name .eq. "heat_flux_dissipation_local" .and. allocated(th_diff)) then
!       call set_published_field_value(field_value, real_1d_field=th_diff)
!     else if (name .eq. "heat_flux_buoyancy_local" .and. allocated(th_buoyancy)) then
!       call set_published_field_value(field_value, real_1d_field=th_buoyancy)
!     else if (name .eq. "heat_flux_tendency_local" .and. allocated(th_tendency)) then
!       call set_published_field_value(field_value, real_1d_field=th_tendency)
!     else if (name .eq. "q_flux_transport_local" .and. allocated(q_flux_values)) then
!       call set_published_field_value(field_value, real_2d_field=q_flux_values)
!     else if (name .eq. "q_flux_gradient_local" .and. allocated(q_gradient)) then
!       call set_published_field_value(field_value, real_2d_field=q_gradient)
!     else if (name .eq. "q_flux_dissipation_local" .and. allocated(q_diff)) then
!       call set_published_field_value(field_value, real_2d_field=q_diff)
!     else if (name .eq. "q_flux_buoyancy_local" .and. allocated(q_buoyancy)) then
!       call set_published_field_value(field_value, real_2d_field=q_buoyancy)
!     else if (name .eq. "q_flux_tendency_local" .and. allocated(q_tendency)) then
!       call set_published_field_value(field_value, real_2d_field=q_tendency)
!     else if (name .eq. "uw_advection_local" .and. allocated(uw_advection)) then
!       call set_published_field_value(field_value, real_1d_field=uw_advection)
!     else if (name .eq. "vw_advection_local" .and. allocated(vw_advection)) then
!       call set_published_field_value(field_value, real_1d_field=vw_advection)
!     else if (name .eq. "uw_viscosity_local" .and. allocated(uw_viscosity)) then
!       call set_published_field_value(field_value, real_1d_field=uw_viscosity)
!     else if (name .eq. "vw_viscosity_local" .and. allocated(vw_viscosity)) then
!       call set_published_field_value(field_value, real_1d_field=vw_viscosity)
!     else if (name .eq. "uw_buoyancy_local" .and. allocated(uw_buoyancy)) then
!       call set_published_field_value(field_value, real_1d_field=uw_buoyancy)
!     else if (name .eq. "vw_buoyancy_local" .and. allocated(vw_buoyancy)) then
!       call set_published_field_value(field_value, real_1d_field=vw_buoyancy)
!     else if (name .eq. "uw_tendency_local" .and. allocated(uw_tendency)) then
!       call set_published_field_value(field_value, real_1d_field=uw_tendency)
!     else if (name .eq. "vw_tendency_local" .and. allocated(vw_tendency)) then
!       call set_published_field_value(field_value, real_1d_field=vw_tendency)
!     else if (name .eq. "uw_w_local" .and. allocated(uw_w)) then
!       call set_published_field_value(field_value, real_1d_field=uw_w)
!     else if (name .eq. "vw_w_local" .and. allocated(vw_w)) then
!       call set_published_field_value(field_value, real_1d_field=vw_w)
!     else if (name .eq. "resolved_pressure_transport_local" .and. allocated(w_p)) then
!       call set_published_field_value(field_value, real_1d_field=w_p)
!     else if (name .eq. "tke_tendency_local" .and. allocated(tend)) then
!       call set_published_field_value(field_value, real_1d_field=tend)
!     else if (name .eq. "resolved_shear_production_local" .and. allocated(shprd)) then
!       call set_published_field_value(field_value, real_1d_field=shprd)
!     else if (name .eq. "resolved_turbulent_transport_local" .and. allocated(w_ke)) then
!       call set_published_field_value(field_value, real_1d_field=w_ke)
!     else if (name .eq. "resolved_buoyant_production_local" .and. allocated(buoy)) then
!       call set_published_field_value(field_value, real_1d_field=buoy)
!     else if (name .eq. "tu_su_local" .and. allocated(tu_su)) then
!       call set_published_field_value(field_value, real_1d_field=tu_su)
!     else if (name .eq. "uu_advection_local" .and. allocated(uu_advection)) then
!       call set_published_field_value(field_value, real_1d_field=uu_advection)
!     else if (name .eq. "uu_viscosity_local" .and. allocated(uu_viscosity)) then
!       call set_published_field_value(field_value, real_1d_field=uu_viscosity)
!     else if (name .eq. "wu_u_local" .and. allocated(wu_u)) then
!       call set_published_field_value(field_value, real_1d_field=wu_u)
!     else if (name .eq. "tv_sv_local" .and. allocated(tv_sv)) then
!       call set_published_field_value(field_value, real_1d_field=tv_sv)
!     else if (name .eq. "vv_advection_local" .and. allocated(vv_advection)) then
!       call set_published_field_value(field_value, real_1d_field=vv_advection)
!     else if (name .eq. "vv_viscosity_local" .and. allocated(vv_viscosity)) then
!       call set_published_field_value(field_value, real_1d_field=vv_viscosity)
!     else if (name .eq. "wv_v_local" .and. allocated(wv_v)) then
!       call set_published_field_value(field_value, real_1d_field=wv_v)
!     else if (name .eq. "tw_sw_local" .and. allocated(tw_sw)) then
!       call set_published_field_value(field_value, real_1d_field=tw_sw)
!     else if (name .eq. "ww_advection_local" .and. allocated(ww_advection)) then
!       call set_published_field_value(field_value, real_1d_field=ww_advection)
!     else if (name .eq. "ww_viscosity_local" .and. allocated(ww_viscosity)) then
!       call set_published_field_value(field_value, real_1d_field=ww_viscosity)
!     else if (name .eq. "ww_buoyancy_local" .and. allocated(ww_buoyancy)) then
!       call set_published_field_value(field_value, real_1d_field=ww_buoyancy)
!     else if (name .eq. "u_thetal_local" .and. allocated(u_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=u_thetal)
!     else if (name .eq. "us_thetal_local" .and. allocated(us_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=us_thetal)
!     else if (name .eq. "u_thetal_advection_local" .and. allocated(u_thetal_advection)) then
!       call set_published_field_value(field_value, real_1d_field=u_thetal_advection)
!     else if (name .eq. "u_thetal_viscosity_diffusion_local" .and. allocated(u_thetal_viscosity_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=u_thetal_viscosity_diffusion)
!     else if (name .eq. "wu_thetal_local" .and. allocated(wu_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=wu_thetal)
!     else if (name .eq. "v_thetal_local" .and. allocated(v_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=v_thetal)
!     else if (name .eq. "vs_thetal_local" .and. allocated(vs_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=vs_thetal)
!     else if (name .eq. "v_thetal_advection_local" .and. allocated(v_thetal_advection)) then
!       call set_published_field_value(field_value, real_1d_field=v_thetal_advection)
!     else if (name .eq. "v_thetal_viscosity_diffusion_local" .and. allocated(v_thetal_viscosity_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=v_thetal_viscosity_diffusion)
!     else if (name .eq. "wv_thetal_local" .and. allocated(wv_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=wv_thetal)
!     else if (name .eq. "w_thetal_local" .and. allocated(w_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=w_thetal)
!     else if (name .eq. "ws_thetal_local" .and. allocated(ws_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=ws_thetal)
!     else if (name .eq. "w_thetal_advection_local" .and. allocated(w_thetal_advection)) then
!       call set_published_field_value(field_value, real_1d_field=w_thetal_advection)
!     else if (name .eq. "w_thetal_viscosity_diffusion_local" .and. allocated(w_thetal_viscosity_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=w_thetal_viscosity_diffusion)
!     else if (name .eq. "w_thetal_buoyancy_local" .and. allocated(w_thetal_buoyancy)) then
!       call set_published_field_value(field_value, real_1d_field=w_thetal_buoyancy)
!     else if (name .eq. "ww_thetal_local" .and. allocated(ww_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=ww_thetal)
!     else if (name .eq. "thetal_thetal_local" .and. allocated(thetal_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=thetal_thetal)
!     else if (name .eq. "sthetal_thetal_local" .and. allocated(sthetal_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=sthetal_thetal)
!     else if (name .eq. "thetal_thetal_advection_local" .and. allocated(thetal_thetal_advection)) then
!       call set_published_field_value(field_value, real_1d_field=thetal_thetal_advection)
!     else if (name .eq. "thetal_thetal_diffusion_local" .and. allocated(thetal_thetal_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=thetal_thetal_diffusion)
!     else if (name .eq. "wthetal_thetal_local" .and. allocated(wthetal_thetal)) then
!       call set_published_field_value(field_value, real_1d_field=wthetal_thetal)
!     else if (name .eq. "u_mse_local" .and. allocated(u_mse)) then
!       call set_published_field_value(field_value, real_1d_field=u_mse)
!     else if (name .eq. "us_mse_local" .and. allocated(us_mse)) then
!       call set_published_field_value(field_value, real_1d_field=us_mse)
!     else if (name .eq. "u_mse_advection_local" .and. allocated(u_mse_advection)) then
!       call set_published_field_value(field_value, real_1d_field=u_mse_advection)
!     else if (name .eq. "u_mse_viscosity_diffusion_local" .and. allocated(u_mse_viscosity_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=u_mse_viscosity_diffusion)
!     else if (name .eq. "wu_mse_local" .and. allocated(wu_mse)) then
!       call set_published_field_value(field_value, real_1d_field=wu_mse)
!     else if (name .eq. "v_mse_local" .and. allocated(v_mse)) then
!       call set_published_field_value(field_value, real_1d_field=v_mse)
!     else if (name .eq. "vs_mse_local" .and. allocated(vs_mse)) then
!       call set_published_field_value(field_value, real_1d_field=vs_mse)
!     else if (name .eq. "v_mse_advection_local" .and. allocated(v_mse_advection)) then
!       call set_published_field_value(field_value, real_1d_field=v_mse_advection)
!     else if (name .eq. "v_mse_viscosity_diffusion_local" .and. allocated(v_mse_viscosity_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=v_mse_viscosity_diffusion)
!     else if (name .eq. "wv_mse_local" .and. allocated(wv_mse)) then
!       call set_published_field_value(field_value, real_1d_field=wv_mse)
!     else if (name .eq. "w_mse_local" .and. allocated(w_mse)) then
!       call set_published_field_value(field_value, real_1d_field=w_mse)
!     else if (name .eq. "ws_mse_local" .and. allocated(ws_mse)) then
!       call set_published_field_value(field_value, real_1d_field=ws_mse)
!     else if (name .eq. "w_mse_advection_local" .and. allocated(w_mse_advection)) then
!       call set_published_field_value(field_value, real_1d_field=w_mse_advection)
!     else if (name .eq. "w_mse_viscosity_diffusion_local" .and. allocated(w_mse_viscosity_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=w_mse_viscosity_diffusion)
!     else if (name .eq. "w_mse_buoyancy_local" .and. allocated(w_mse_buoyancy)) then
!       call set_published_field_value(field_value, real_1d_field=w_mse_buoyancy)
!     else if (name .eq. "ww_mse_local" .and. allocated(ww_mse)) then
!       call set_published_field_value(field_value, real_1d_field=ww_mse)
!     else if (name .eq. "mse_mse_local" .and. allocated(mse_mse)) then
!       call set_published_field_value(field_value, real_1d_field=mse_mse)
!     else if (name .eq. "smse_mse_local" .and. allocated(smse_mse)) then
!       call set_published_field_value(field_value, real_1d_field=smse_mse)
!     else if (name .eq. "mse_mse_advection_local" .and. allocated(mse_mse_advection)) then
!       call set_published_field_value(field_value, real_1d_field=mse_mse_advection)
!     else if (name .eq. "mse_mse_diffusion_local" .and. allocated(mse_mse_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=mse_mse_diffusion)
!     else if (name .eq. "wmse_mse_local" .and. allocated(wmse_mse)) then
!       call set_published_field_value(field_value, real_1d_field=wmse_mse)
!     else if (name .eq. "us_qt_local" .and. allocated(us_qt)) then
!       call set_published_field_value(field_value, real_1d_field=us_qt)
!     else if (name .eq. "u_qt_advection_local" .and. allocated(u_qt_advection)) then
!       call set_published_field_value(field_value, real_1d_field=u_qt_advection)
!     else if (name .eq. "u_qt_viscosity_diffusion_local" .and. allocated(u_qt_viscosity_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=u_qt_viscosity_diffusion)
!     else if (name .eq. "wu_qt_local" .and. allocated(wu_qt)) then
!       call set_published_field_value(field_value, real_1d_field=wu_qt)
!     else if (name .eq. "vs_qt_local" .and. allocated(vs_qt)) then
!       call set_published_field_value(field_value, real_1d_field=vs_qt)
!     else if (name .eq. "v_qt_advection_local" .and. allocated(v_qt_advection)) then
!       call set_published_field_value(field_value, real_1d_field=v_qt_advection)
!     else if (name .eq. "v_qt_viscosity_diffusion_local" .and. allocated(v_qt_viscosity_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=v_qt_viscosity_diffusion)
!     else if (name .eq. "wv_qt_local" .and. allocated(wv_qt)) then
!       call set_published_field_value(field_value, real_1d_field=wv_qt)
!     else if (name .eq. "w_qt_local" .and. allocated(w_qt)) then
!       call set_published_field_value(field_value, real_1d_field=w_qt)
!     else if (name .eq. "ws_qt_local" .and. allocated(ws_qt)) then
!       call set_published_field_value(field_value, real_1d_field=ws_qt)
!     else if (name .eq. "w_qt_advection_local" .and. allocated(w_qt_advection)) then
!       call set_published_field_value(field_value, real_1d_field=w_qt_advection)
!     else if (name .eq. "w_qt_viscosity_diffusion_local" .and. allocated(w_qt_viscosity_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=w_qt_viscosity_diffusion)
!     else if (name .eq. "w_qt_buoyancy_local" .and. allocated(w_qt_buoyancy)) then
!       call set_published_field_value(field_value, real_1d_field=w_qt_buoyancy)
!     else if (name .eq. "ww_qt_local" .and. allocated(ww_qt)) then
!       call set_published_field_value(field_value, real_1d_field=ww_qt)
!     else if (name .eq. "qt_qt_local" .and. allocated(qt_qt)) then
!       call set_published_field_value(field_value, real_1d_field=qt_qt)
!     else if (name .eq. "sqt_qt_local" .and. allocated(sqt_qt)) then
!       call set_published_field_value(field_value, real_1d_field=sqt_qt)
!     else if (name .eq. "qt_qt_advection_local" .and. allocated(qt_qt_advection)) then
!       call set_published_field_value(field_value, real_1d_field=qt_qt_advection)
!     else if (name .eq. "qt_qt_diffusion_local" .and. allocated(qt_qt_diffusion)) then
!       call set_published_field_value(field_value, real_1d_field=qt_qt_diffusion)
!     else if (name .eq. "wqt_qt_local" .and. allocated(wqt_qt)) then
!       call set_published_field_value(field_value, real_1d_field=wqt_qt)
!     else if (name .eq. "mflux_local") then
!       field_value%scalar_real=mflux
!     end if
!   end subroutine field_value_retrieval_callback

  !> Sets the published field value from the temporary diagnostic values held by this component.
  !! @param field_value Populated with the value of the field
  !! @param real_1d_field Optional one dimensional real of values to publish
  !! @param real_2d_field Optional two dimensional real of values to publish
!   subroutine set_published_field_value(field_value, real_1d_field, real_2d_field)
!     type(component_field_value_type), intent(inout) :: field_value
!     real(kind=DEFAULT_PRECISION), dimension(:), optional :: real_1d_field
!     real(kind=DEFAULT_PRECISION), dimension(:,:), optional :: real_2d_field
!
!     if (present(real_1d_field)) then
!       allocate(field_value%real_1d_array(size(real_1d_field)), source=real_1d_field)
!     else if (present(real_2d_field)) then
!       allocate(field_value%real_2d_array(size(real_2d_field, 1), size(real_2d_field, 2)), source=real_2d_field)
!     end if
!   end subroutine set_published_field_value
!
!   !> Sets the published value enabled state in the provided collection map
!   !! @param collection The map to set this in
!   !! @param field_name The name of the published field
!   !! @param enabled_state Whether the field is enabled which is stored
!   subroutine set_published_field_enabled_state(collection, field_name, enabled_state)
!     type(hashmap_type), intent(inout) :: collection
!     character(len=*), intent(in) :: field_name
!     logical, intent(in) :: enabled_state
!
!     call c_put_logical(collection, field_name, enabled_state)
!   end subroutine set_published_field_enabled_state
!
!   !> Retrieves whether a published field is enabled or not
!   !! @param collection The map to look up the published field in
!   !! @param field_name The name of the field to look up
!   !! @returns Whether this published field is enabled or not
!   logical function get_published_field_enabled_state(collection, field_name)
!     type(hashmap_type), intent(inout) :: collection
!     character(len=*), intent(in) :: field_name
!
!     get_published_field_enabled_state=c_get_logical(collection, field_name)
!   end function get_published_field_enabled_state

end module flux_budget_mod
