module casim_profile_dgs_mod
  use monc_component_mod, only :  COMPONENT_DOUBLE_DATA_TYPE, COMPONENT_ARRAY_FIELD_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type, &
       component_descriptor_type_v1
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use logging_mod, only : LOG_ERROR, log_master_log
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION

  ! needs casim modules
  use mphys_switches, only: l_warm
  ! need casdiags for the diag switches
  use generic_diagnostic_variables, ONLY: casdiags
  ! and casim component structure, which contains the process rate data
  use casim_monc_dgs_space, only: casim_monc_dgs

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: total_points, iqv, iql
!   real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
!        ! local process rate totals
!        phomc_tot, pinuc_tot, pidep_tot, psdep_tot, piacw_tot, psacw_tot, psacr_tot, pisub_tot,   &
!        pssub_tot, pimlt_tot, psmlt_tot, psaut_tot, psaci_tot, praut_tot, pracw_tot, prevp_tot,   &
!        pgacw_tot, pgacs_tot, pgmlt_tot, pgsub_tot, psedi_tot, pseds_tot, psedr_tot, psedg_tot,   &
!        psedl_tot, pcond_tot, &
!        ! local total tendencies
!        dth_mphys_tot, dth_cond_evap_tot, dqv_mphys_tot, dqv_cond_evap_tot, &
!        dqc_mphys_tot, dqr_mphys_tot, dqi_mphys_tot, dqs_mphys_tot, dqg_mphys_tot

  public casim_profile_dgs_get_descriptor, initialisation_callback_casim_profile_dgs, timestep_callback_casim_profile_dgs

contains

  !> Provides the component descriptor for the core to register
  !! @returns The descriptor describing this component
  type(component_descriptor_type_v1) function casim_profile_dgs_get_descriptor()
    casim_profile_dgs_get_descriptor%name="casim_profile_dgs"
    casim_profile_dgs_get_descriptor%version=0.1

    !casim_profile_dgs_get_descriptor%initialisation=>initialisation_callback
    !casim_profile_dgs_get_descriptor%timestep=>timestep_callback

    !casim_profile_dgs_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    !casim_profile_dgs_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    casim_profile_dgs_get_descriptor%published_fields_on_off = .true.
    allocate(casim_profile_dgs_get_descriptor%published_fields(35))

    casim_profile_dgs_get_descriptor%published_fields(1)="phomc_total"
    casim_profile_dgs_get_descriptor%published_fields(2)="pinuc_total"
    casim_profile_dgs_get_descriptor%published_fields(3)="pidep_total"
    casim_profile_dgs_get_descriptor%published_fields(4)="psdep_total"
    casim_profile_dgs_get_descriptor%published_fields(5)="piacw_total"
    casim_profile_dgs_get_descriptor%published_fields(6)="psacw_total"
    casim_profile_dgs_get_descriptor%published_fields(7)="psacr_total"
    casim_profile_dgs_get_descriptor%published_fields(8)="pisub_total"
    casim_profile_dgs_get_descriptor%published_fields(9)="pssub_total"
    casim_profile_dgs_get_descriptor%published_fields(10)="pimlt_total"
    casim_profile_dgs_get_descriptor%published_fields(11)="psmlt_total"
    casim_profile_dgs_get_descriptor%published_fields(12)="psaut_total"
    casim_profile_dgs_get_descriptor%published_fields(13)="psaci_total"
    casim_profile_dgs_get_descriptor%published_fields(14)="praut_total"
    casim_profile_dgs_get_descriptor%published_fields(15)="pracw_total"
    casim_profile_dgs_get_descriptor%published_fields(16)="prevp_total"
    casim_profile_dgs_get_descriptor%published_fields(17)="pgacw_total"
    casim_profile_dgs_get_descriptor%published_fields(18)="pgacs_total"
    casim_profile_dgs_get_descriptor%published_fields(19)="pgmlt_total"
    casim_profile_dgs_get_descriptor%published_fields(20)="pgsub_total"
    casim_profile_dgs_get_descriptor%published_fields(21)="psedi_total"
    casim_profile_dgs_get_descriptor%published_fields(22)="pseds_total"
    casim_profile_dgs_get_descriptor%published_fields(23)="psedr_total"
    casim_profile_dgs_get_descriptor%published_fields(24)="psedg_total"
    casim_profile_dgs_get_descriptor%published_fields(25)="psedl_total"
    casim_profile_dgs_get_descriptor%published_fields(26)="pcond_total"
    casim_profile_dgs_get_descriptor%published_fields(27)="dth_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(28)="dth_cond_evap_total"
    casim_profile_dgs_get_descriptor%published_fields(29)="dqv_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(30)="dqv_cond_evap_total"
    casim_profile_dgs_get_descriptor%published_fields(31)="dqc_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(32)="dqr_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(33)="dqi_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(34)="dqs_mphys_total"
    casim_profile_dgs_get_descriptor%published_fields(35)="dqg_mphys_total"

  end function casim_profile_dgs_get_descriptor


  subroutine initialisation_callback_casim_profile_dgs(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, iter
    logical :: logicnum, casim_enabled

    !if (.not. is_component_enabled(current_state%options_database, "casim")) then
    !  call log_master_log(LOG_ERROR, "Casim profile diags requested but casim not enabled - check config, STOP")
    !end if
    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "casim_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        casim_enabled = logicnum
        if(.not. casim_enabled) then
          call log_master_log(LOG_ERROR, "Casim profile diags requested but casim not enabled - check config, STOP")
        end if
      end if
    end do

    ! allocate local arrays for the horizontal wind averages
    allocate(current_state%psedl_tot(current_state%local_grid%size(Z_INDEX)), &
         current_state%pcond_tot(current_state%local_grid%size(Z_INDEX)), &
         current_state%praut_tot(current_state%local_grid%size(Z_INDEX)), &
         current_state%pracw_tot(current_state%local_grid%size(Z_INDEX)), &
         current_state%prevp_tot(current_state%local_grid%size(Z_INDEX)), &
         current_state%psedr_tot(current_state%local_grid%size(Z_INDEX)), &
         current_state%dth_mphys_tot(current_state%local_grid%size(Z_INDEX)),     &
         current_state%dth_cond_evap_tot(current_state%local_grid%size(Z_INDEX)), &
         current_state%dqv_mphys_tot(current_state%local_grid%size(Z_INDEX)),     &
         current_state%dqv_cond_evap_tot(current_state%local_grid%size(Z_INDEX)), &
         current_state%dqc_mphys_tot(current_state%local_grid%size(Z_INDEX)),     &
         current_state%dqr_mphys_tot(current_state%local_grid%size(Z_INDEX)))
    !if (.not. l_warm) then 
       allocate(current_state%phomc_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pinuc_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pidep_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%piacw_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pisub_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pimlt_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%psedi_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%psacw_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%psacr_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pssub_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%psmlt_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%psaut_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%psaci_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%psdep_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pseds_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pgacw_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pgacs_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pgmlt_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%pgsub_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%psedg_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%dqi_mphys_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%dqs_mphys_tot(current_state%local_grid%size(Z_INDEX)), &
            current_state%dqg_mphys_tot(current_state%local_grid%size(Z_INDEX)))
    !endif
       
  end subroutine initialisation_callback_casim_profile_dgs

  subroutine timestep_callback_casim_profile_dgs(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i
    integer :: icol, jcol, target_x_index, target_y_index
    
    icol=current_state%column_local_x
    jcol=current_state%column_local_y
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)

    if (current_state%first_timestep_column) then
       current_state%psedl_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%pcond_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%praut_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%pracw_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%prevp_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%psedr_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%dth_mphys_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%dth_cond_evap_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%dqv_mphys_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%dqv_cond_evap_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%dqc_mphys_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%dqr_mphys_tot(:)= 0.0_DEFAULT_PRECISION
       !if (.not. l_warm) then
          current_state%phomc_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pinuc_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pidep_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%piacw_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pisub_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pimlt_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%psedi_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%psacw_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%psacr_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pssub_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%psmlt_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%psaut_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%psaci_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%psdep_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pseds_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pgacw_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pgacs_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pgmlt_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%pgsub_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%psedg_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%dqi_mphys_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%dqs_mphys_tot(:)= 0.0_DEFAULT_PRECISION
          current_state%dqg_mphys_tot(:)= 0.0_DEFAULT_PRECISION
       !endif
    endif

    if (.not. current_state%halo_column) then
       if ( casdiags % l_psedl ) &
            current_state%psedl_tot(:)= current_state%psedl_tot(:) + &
            casim_monc_dgs % psedl(:,target_y_index,target_x_index)
       if ( casdiags % l_pcond ) &
            current_state%pcond_tot(:)= current_state%pcond_tot(:) + &
            casim_monc_dgs % pcond(:,target_y_index,target_x_index)
       if (  casdiags % l_praut ) & 
            current_state%praut_tot(:)= current_state%praut_tot(:) + &
            casim_monc_dgs % praut(:,target_y_index,target_x_index)
       if ( casdiags % l_pracw ) & 
            current_state%pracw_tot(:)= current_state%pracw_tot(:) + &
            casim_monc_dgs % pracw(:,target_y_index,target_x_index)
       if ( casdiags % l_prevp ) & 
            current_state%prevp_tot(:)= current_state%prevp_tot(:) + &
            casim_monc_dgs % prevp(:,target_y_index,target_x_index)
       if ( casdiags % l_psedr ) & 
            current_state%psedr_tot(:)= current_state%psedr_tot(:) + &
            casim_monc_dgs % psedr(:,target_y_index,target_x_index)
       if ( casdiags % l_dth ) then
          current_state%dth_mphys_tot(:)= current_state%dth_mphys_tot(:) + &
               casim_monc_dgs %  dth_total(:,target_y_index,target_x_index)
          
          current_state%dth_cond_evap_tot(:)= current_state%dth_cond_evap_tot(:) + &
               casim_monc_dgs % dth_cond_evap(:,target_y_index,target_x_index)
       endif

       if ( casdiags % l_dqv ) then
          current_state%dqv_mphys_tot(:)= current_state%dqv_mphys_tot(:) + &
               casim_monc_dgs % dqv_total(:,target_y_index,target_x_index)
          current_state%dqv_cond_evap_tot(:)= current_state%dqv_cond_evap_tot(:) + &
               casim_monc_dgs % dqv_cond_evap(:,target_y_index,target_x_index)
       endif

       if ( casdiags % l_dqc ) &
            current_state%dqc_mphys_tot(:)= current_state%dqc_mphys_tot(:) + &
            casim_monc_dgs % dqc(:,target_y_index,target_x_index)
       if ( casdiags % l_dqr ) & 
            current_state%dqr_mphys_tot(:)= current_state%dqr_mphys_tot(:) + &
            casim_monc_dgs % dqr(:,target_y_index,target_x_index)
       
       if (.not. l_warm) then
          if ( casdiags % l_phomc ) & 
               current_state%phomc_tot(:)= current_state%phomc_tot(:) + &
               casim_monc_dgs % phomc(:,target_y_index,target_x_index)
          if ( casdiags % l_pinuc ) & 
               current_state%pinuc_tot(:)= current_state%pinuc_tot(:) + &
               casim_monc_dgs % pinuc(:,target_y_index,target_x_index)
          if ( casdiags % l_pidep ) & 
               current_state%pidep_tot(:)= current_state%pidep_tot(:) + &
               casim_monc_dgs % pidep(:,target_y_index,target_x_index)
          if ( casdiags % l_piacw ) & 
               current_state%piacw_tot(:)= current_state%piacw_tot(:) + &
               casim_monc_dgs % piacw(:,target_y_index,target_x_index)
          if ( casdiags % l_pisub ) &
               current_state%pisub_tot(:)= current_state%pisub_tot(:) + &
               casim_monc_dgs % pisub(:,target_y_index,target_x_index)
          if ( casdiags % l_pimlt ) &
               current_state%pimlt_tot(:)= current_state%pimlt_tot(:) + &
               casim_monc_dgs % pimlt(:,target_y_index,target_x_index)
          if ( casdiags % l_psedi ) & 
               current_state%psedi_tot(:)= current_state%psedi_tot(:) + &
               casim_monc_dgs % psedi(:,target_y_index,target_x_index)
          if ( casdiags % l_psacw ) & 
               current_state%psacw_tot(:)= current_state%psacw_tot(:) + &
               casim_monc_dgs % psacw(:,target_y_index,target_x_index)
          if ( casdiags % l_psacr ) & 
               current_state%psacr_tot(:)= current_state%psacr_tot(:) + &
               casim_monc_dgs % psacr(:,target_y_index,target_x_index)
          if ( casdiags % l_pssub ) & 
               current_state%pssub_tot(:)= current_state%pssub_tot(:) + &
               casim_monc_dgs % pssub(:,target_y_index,target_x_index)
          if ( casdiags % l_psmlt ) & 
               current_state%psmlt_tot(:)= current_state%psmlt_tot(:) + &
               casim_monc_dgs % psmlt(:,target_y_index,target_x_index)
          if ( casdiags % l_psaut ) &
               current_state%psaut_tot(:)= current_state%psaut_tot(:) + &
               casim_monc_dgs % psaut(:,target_y_index,target_x_index)
          if ( casdiags % l_psaci ) &
               current_state%psaci_tot(:)= current_state%psaci_tot(:) + &
               casim_monc_dgs % psaci(:,target_y_index,target_x_index)
          if ( casdiags % l_psdep ) & 
               current_state%psdep_tot(:)= current_state%psdep_tot(:) + &
               casim_monc_dgs % psdep(:,target_y_index,target_x_index)
          if ( casdiags % l_pseds ) & 
               current_state%pseds_tot(:)= current_state%pseds_tot(:) + &
               casim_monc_dgs % pseds(:,target_y_index,target_x_index)
          if ( casdiags % l_pgacw ) & 
               current_state%pgacw_tot(:)= current_state%pgacw_tot(:) + &
               casim_monc_dgs % pgacw(:,target_y_index,target_x_index)
          if ( casdiags % l_pgacs ) &
               current_state%pgacs_tot(:)= current_state%pgacs_tot(:) + &
               casim_monc_dgs % pgacs(:,target_y_index,target_x_index)
          if ( casdiags % l_pgmlt ) &
               current_state%pgmlt_tot(:)= current_state%pgmlt_tot(:) + &
               casim_monc_dgs % pgmlt(:,target_y_index,target_x_index)
          if ( casdiags % l_pgsub ) &
               current_state%pgsub_tot(:)= current_state%pgsub_tot(:) + &
               casim_monc_dgs % pgsub(:,target_y_index,target_x_index)
          if ( casdiags % l_psedg ) & 
               current_state%psedg_tot(:)= current_state%psedg_tot(:) + &
               casim_monc_dgs % psedg(:,target_y_index,target_x_index)
          if ( casdiags % l_dqi ) & 
               current_state%dqi_mphys_tot(:)= current_state%dqi_mphys_tot(:) + &
               casim_monc_dgs % dqi(:,target_y_index,target_x_index)
          if ( casdiags % l_dqs ) & 
               current_state%dqs_mphys_tot(:)= current_state%dqs_mphys_tot(:) + &
               casim_monc_dgs % dqs(:,target_y_index,target_x_index)
          if ( casdiags % l_dqg ) &
               current_state%dqg_mphys_tot(:)= current_state%dqg_mphys_tot(:) + &
               casim_monc_dgs % dqg(:,target_y_index,target_x_index)
       endif
       
    endif
  end subroutine timestep_callback_casim_profile_dgs

  !> Field information retrieval callback, this returns information for a specific components published field
  !! @param current_state Current model state
  !! @param name The name of the field to retrieve information for
  !! @param field_information Populated with information about the field
!   subroutine field_information_retrieval_callback(current_state, name, field_information)
!     type(model_state_type), target, intent(inout) :: current_state
!     character(len=*), intent(in) :: name
!     type(component_field_information_type), intent(out) :: field_information
!
!     field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
!     field_information%number_dimensions=1
!     field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
!     field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE
!     if (l_warm) then
!        if (name .eq. "pcond_total" .or. name .eq. "praut_total" &
!             .or. name .eq. "pracw_total" .or. name .eq. "prevp_total" &
!             .or. name .eq. "psedl_total" .or. name .eq. "psedr_total" &
!             .or. name .eq. "dth_mphys_total" .or. name .eq. "dth_cond_evap_total" &
!             .or. name .eq. "dqv_mphys_total" .or. name .eq. "dqv_cond_evap_total" &
!             .or. name .eq. "dqc_mphys_total" .or. name .eq. "dqr_mphys_total") then
!           field_information%enabled=.true.
!        endif
!     else
!        field_information%enabled=.true.
!     endif
!
!   end subroutine field_information_retrieval_callback
!
!   !> Field value retrieval callback, this returns the value of a specific published field
!   !! @param current_state Current model state
!   !! @param name The name of the field to retrieve the value for
!   !! @param field_value Populated with the value of the field
!   subroutine field_value_retrieval_callback(current_state, name, field_value)
!     type(model_state_type), target, intent(inout) :: current_state
!     character(len=*), intent(in) :: name
!     type(component_field_value_type), intent(out) :: field_value
!
!     integer :: k
!
!     if (name .eq. "phomc_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=phomc_tot(:)
!     else if (name .eq. "pinuc_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pinuc_tot(:)
!     else if (name .eq. "pidep_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pidep_tot(:)
!     else if (name .eq. "psdep_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psdep_tot(:)
!     else if (name .eq. "piacw_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=piacw_tot(:)
!     else if (name .eq. "psacw_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psacw_tot(:)
!     else if (name .eq. "psacr_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psacr_tot(:)
!     else if (name .eq. "pisub_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pisub_tot(:)
!     else if (name .eq. "pssub_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pssub_tot(:)
!     else if (name .eq. "pimlt_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pimlt_tot(:)
!     else if (name .eq. "psmlt_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psmlt_tot(:)
!     else if (name .eq. "psaut_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psaut_tot(:)
!     else if (name .eq. "psaci_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psaci_tot(:)
!     else if (name .eq. "praut_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=praut_tot(:)
!     else if (name .eq. "pracw_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pracw_tot(:)
!     else if (name .eq. "prevp_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=prevp_tot(:)
!     else if (name .eq. "pgacw_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pgacw_tot(:)
!     else if (name .eq. "pgacs_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pgacs_tot(:)
!     else if (name .eq. "pgmlt_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pgmlt_tot(:)
!     else if (name .eq. "pgsub_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pgsub_tot(:)
!     else if (name .eq. "psedi_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psedi_tot(:)
!     else if (name .eq. "pseds_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pseds_tot(:)
!     else if (name .eq. "psedr_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psedr_tot(:)
!     else if (name .eq. "psedg_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psedg_tot(:)
!     else if (name .eq. "psedl_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=psedl_tot(:)
!     else if (name .eq. "pcond_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=pcond_tot(:)
!     else if (name .eq. "dth_mphys_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=dth_mphys_tot(:)
!     else if (name .eq. "dth_cond_evap_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=dth_cond_evap_tot(:)
!     else if (name .eq. "dqv_mphys_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=dqv_mphys_tot(:)
!     else if (name .eq. "dqv_cond_evap_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=dqv_cond_evap_tot(:)
!     else if (name .eq. "dqc_mphys_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=dqc_mphys_tot(:)
!     else if (name .eq. "dqr_mphys_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=dqr_mphys_tot(:)
!     else if (name .eq. "dqi_mphys_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=dqi_mphys_tot(:)
!     else if (name .eq. "dqs_mphys_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=dqs_mphys_tot(:)
!     else if (name .eq. "dqg_mphys_total") then
!        allocate(field_value%real_1d_array(current_state%local_grid%size(Z_INDEX)))
!        field_value%real_1d_array(:)=dqg_mphys_tot(:)
!     endif
!   end subroutine field_value_retrieval_callback
end module casim_profile_dgs_mod
