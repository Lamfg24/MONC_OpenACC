










!> This component sets the source term for the lowest level (Level 1) so 
!> that, depending on surface consitionm, there is consistent lower boundary 
!> condition.
module set_consistent_lowbc_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : PRESCRIBED_SURFACE_FLUXES, PRESCRIBED_SURFACE_VALUES, &
                        model_state_type
  use grids_mod, only : Z_INDEX
  use q_indices_mod, only: get_q_index, standard_q_names

  logical :: advect_flow, advect_th, advect_q

  integer :: iqv  ! index for vapour
  
  public initialisation_callback_set_consistent_lowbc, timestep_callback_set_consistent_lowbc

contains

  subroutine initialisation_callback_set_consistent_lowbc(current_state)
    ! copy of the initialisation callback from pwadvection, used to 
    ! identify whether acting on flow, theta or q fields
    type(model_state_type), target, intent(inout) :: current_state

    ! Determine vapour index
    if (.not. current_state%passive_q) then 
       iqv = 1!get_q_index(standard_q_names%VAPOUR, 'lowerbc')
    endif
    
  end subroutine initialisation_callback_set_consistent_lowbc

  subroutine timestep_callback_set_consistent_lowbc(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: current_x_index, current_y_index

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y

    if (current_state%halo_column) return
    if (current_state%immersed%ib_enabled)then
      if(current_state%immersed%ib_col(current_y_index, current_x_index)) return
    end if

    call set_flow_lowbc(current_state, current_x_index, current_y_index)
    if (current_state%th%active) then 
       call set_th_lowbc(current_state, current_x_index, current_y_index)
    endif
    if (current_state%number_q_fields .gt. 0) then
       call set_q_lowbc(current_state, current_x_index, current_y_index)
    end if
    if (current_state%n_tracers .gt. 0) then
       call set_tracer_lowbc(current_state, current_x_index, current_y_index)
    end if
  end subroutine timestep_callback_set_consistent_lowbc

  subroutine set_flow_lowbc(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index
    
    if (current_state%use_viscosity_and_diffusion .and. &
        current_state%use_surface_boundary_conditions) then 
        current_state%su%data(1, current_y_index, current_x_index)= &
               -current_state%su%data(2, current_y_index, current_x_index)
        current_state%sv%data(1, current_y_index, current_x_index)= &
               -current_state%sv%data(2, current_y_index, current_x_index)   
    else 
       current_state%su%data(1, current_y_index, current_x_index)= &               
                current_state%su%data(2, current_y_index, current_x_index)
       current_state%sv%data(1, current_y_index, current_x_index)= &
                current_state%sv%data(2, current_y_index, current_x_index)
    endif

  end subroutine set_flow_lowbc

  subroutine set_th_lowbc(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    if (current_state%use_surface_boundary_conditions) then 
       if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then
          current_state%sth%data(1, current_y_index, current_x_index)= &
               current_state%sth%data(2, current_y_index, current_x_index)
       else
          current_state%sth%data(1, current_y_index, current_x_index)= &
               -current_state%sth%data(2, current_y_index, current_x_index)        
       endif
    endif

  end subroutine set_th_lowbc

  subroutine set_q_lowbc(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: i  
    
    if (current_state%number_q_fields .gt. 0) then
       !do n=1,current_state%number_q_fields
      current_state%sqv%data(1, current_y_index, current_x_index)= &
            current_state%sqv%data(2,current_y_index, current_x_index)
      current_state%sql%data(1, current_y_index, current_x_index)= &
            current_state%sql%data(2,current_y_index, current_x_index)
      current_state%sqr%data(1, current_y_index, current_x_index)= &
            current_state%sqr%data(2,current_y_index, current_x_index)
      current_state%sqi%data(1, current_y_index, current_x_index)= &
            current_state%sqi%data(2,current_y_index, current_x_index)
      current_state%sqs%data(1, current_y_index, current_x_index)= &
            current_state%sqs%data(2,current_y_index, current_x_index)
      current_state%sqg%data(1, current_y_index, current_x_index)= &
            current_state%sqg%data(2,current_y_index, current_x_index)
      current_state%sqAitkenSolMass%data(1, current_y_index, current_x_index)= &
            current_state%sqAitkenSolMass%data(2,current_y_index, current_x_index)
      current_state%sqAccumSolMass%data(1, current_y_index, current_x_index)= &
            current_state%sqAccumSolMass%data(2,current_y_index, current_x_index)
      current_state%sqAccumInsolMass%data(1, current_y_index, current_x_index)= &
            current_state%sqAccumInsolMass%data(2,current_y_index, current_x_index)
      current_state%sqCoarseSolMass%data(1, current_y_index, current_x_index)= &
            current_state%sqCoarseSolMass%data(2,current_y_index, current_x_index)
      current_state%sqCoarseDustMass%data(1, current_y_index, current_x_index)= &
            current_state%sqCoarseDustMass%data(2,current_y_index, current_x_index)

      current_state%snl%data(1, current_y_index, current_x_index)= &
            current_state%snl%data(2,current_y_index, current_x_index)
      current_state%snr%data(1, current_y_index, current_x_index)= &
            current_state%snr%data(2,current_y_index, current_x_index)
      current_state%sni%data(1, current_y_index, current_x_index)= &
            current_state%sni%data(2,current_y_index, current_x_index)
      current_state%sns%data(1, current_y_index, current_x_index)= &
            current_state%sns%data(2,current_y_index, current_x_index)
      current_state%sng%data(1, current_y_index, current_x_index)= &
            current_state%sng%data(2,current_y_index, current_x_index)
      current_state%snAitkenSolNumber%data(1, current_y_index, current_x_index)= &
            current_state%snAitkenSolNumber%data(2,current_y_index, current_x_index)
      current_state%snAccumSolNumber%data(1, current_y_index, current_x_index)= &
            current_state%snAccumSolNumber%data(2,current_y_index, current_x_index)
      current_state%snAccumInsolNumber%data(1, current_y_index, current_x_index)= &
            current_state%snAccumInsolNumber%data(2,current_y_index, current_x_index)
      current_state%snCoarseSolNumber%data(1, current_y_index, current_x_index)= &
            current_state%snCoarseSolNumber%data(2,current_y_index, current_x_index)
      current_state%snCoarseDustnumber%data(1, current_y_index, current_x_index)= &
            current_state%snCoarseDustnumber%data(2,current_y_index, current_x_index)
       !end do
    end if
    if (current_state%use_surface_boundary_conditions .and. &
         current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then
       current_state%sqv%data(1, current_y_index, current_x_index)= &
            -(current_state%sqv%data(2,current_y_index, current_x_index))
    end if

  end subroutine set_q_lowbc

  subroutine set_tracer_lowbc(current_state, current_x_index, current_y_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  current_x_index, current_y_index

    integer :: i  
    
    if (current_state%n_tracers .gt. 0) then
       do n=1,current_state%n_tracers
          current_state%stracer(n)%data(1, current_y_index, current_x_index)= &
               current_state%stracer(n)%data(2,current_y_index, current_x_index)
       enddo
    endif

  end subroutine set_tracer_lowbc

end module set_consistent_lowbc_mod
