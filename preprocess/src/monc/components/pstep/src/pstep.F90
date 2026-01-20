










!> Stepping of the pressure field. Completes the time-stepping of the velocity fields
!! by adding the pressure term (dp/dx_i). In addition, ensures that l_zu and l_zv satisfy the
!! Galilean-transformed boundary condition. This does not do the flow field _p terms which are only
!! needed for diagnostics, nore does it do field halo swapping which is again only needed for diagnostics
module pstep_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type, &
      component_descriptor_type_v1
  use optionsdatabase_mod, only : options_get_integer
  use state_mod, only : model_state_type, CENTRED_STEPPING
  use grids_mod, only : Z_INDEX, X_INDEX, Y_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  use logging_mod, only : LOG_ERROR, log_master_log

  implicit none

  private

  ! Local tendency diagnostic variables for this component
  !
  ! This one is a little different.  tendp_ terms are after the tendency due to the pressure term,
  ! and the tend_ terms will collect the total tendency for the time step.
  !
  ! 3D tendency fields and logicals for their use
  logical :: l_tendp_3d_u, l_tendp_3d_v, l_tendp_3d_w, l_tend_3d_u, l_tend_3d_v, l_tend_3d_w
  logical :: l_tendp_pr_tot_u, l_tendp_pr_tot_v, l_tendp_pr_tot_w, l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w

  public initialisation_callback_pstep, timestep_callback_pstep, finalisation_callback_pstep

contains

  !> Initialisation callback hook which will check the diverr component is enabled (as this allocates p)
  !! @param current_state The current model state
  subroutine initialisation_callback_pstep(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: iter
    integer :: alloc_z, alloc_y, alloc_x
    logical :: logicnum

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "diverr_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (.not. logicnum) then
          call log_master_log(LOG_ERROR, "The pstep component requires the diverr component to be enabled")
        end if
      end if
    end do

    alloc_z=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    alloc_y=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    alloc_x=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2

    ! Tendency Logicals
    l_tendp_pr_tot_u   = current_state%u%active
    l_tendp_pr_tot_v   = current_state%v%active
    l_tendp_pr_tot_w   = current_state%w%active
    l_tend_pr_tot_u   = current_state%u%active
    l_tend_pr_tot_v   = current_state%v%active
    l_tend_pr_tot_w   = current_state%w%active

    l_tendp_3d_u   = current_state%u%active .or. l_tendp_pr_tot_u
    l_tendp_3d_v   = current_state%v%active .or. l_tendp_pr_tot_v
    l_tendp_3d_w   = current_state%w%active .or. l_tendp_pr_tot_w
    l_tend_3d_u   = current_state%u%active .or. l_tend_pr_tot_u
    l_tend_3d_v   = current_state%v%active .or. l_tend_pr_tot_v
    l_tend_3d_w   = current_state%w%active .or. l_tend_pr_tot_w

    ! Allocate 3d tendency fields upon availability
    if (l_tendp_3d_u) then
      !allocate( current_state%tendp_3d_u_pt(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tendp_3d_u_pt(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tendp_3d_v) then
      !allocate( current_state%tendp_3d_v_pt(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tendp_3d_v_pt(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tendp_3d_w) then
      !allocate( current_state%tendp_3d_w_pt(current_state%local_grid%size(Z_INDEX),  &
      !                     current_state%local_grid%size(Y_INDEX),  &
      !                     current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tendp_3d_w_pt(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_u) then
      !allocate( current_state%tend_3d_u_pt(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_u_pt(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_v) then
      !allocate( current_state%tend_3d_v_pt(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_v_pt(alloc_z, alloc_y, alloc_x))
    endif
    if (l_tend_3d_w) then
      !allocate( current_state%tend_3d_w_pt(current_state%local_grid%size(Z_INDEX),  &
      !                    current_state%local_grid%size(Y_INDEX),  &
      !                    current_state%local_grid%size(X_INDEX)   )    )
      allocate(current_state%tend_3d_w_pt(alloc_z, alloc_y, alloc_x))
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tendp_pr_tot_u) then
      !allocate( current_state%tendp_pr_tot_u_pt(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tendp_pr_tot_u_pt(alloc_z))
    endif
    if (l_tendp_pr_tot_v) then
      !allocate( current_state%tendp_pr_tot_v_pt(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tendp_pr_tot_v_pt(alloc_z))
    endif
    if (l_tendp_pr_tot_w) then
      !allocate( current_state%tendp_pr_tot_w_pt(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tendp_pr_tot_w_pt(alloc_z))
    endif
    if (l_tend_pr_tot_u) then
      !allocate( current_state%tend_pr_tot_u_pt(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_u_pt(alloc_z))
    endif
    if (l_tend_pr_tot_v) then
      !allocate( current_state%tend_pr_tot_v_pt(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_v_pt(alloc_z))
    endif
    if (l_tend_pr_tot_w) then
      !allocate( current_state%tend_pr_tot_w_pt(current_state%local_grid%size(Z_INDEX)) )
      allocate(current_state%tend_pr_tot_w_pt(alloc_z))
    endif

  end subroutine initialisation_callback_pstep


  subroutine finalisation_callback_pstep(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(current_state%tendp_3d_u_pt)) deallocate(current_state%tendp_3d_u_pt)
    if (allocated(current_state%tendp_3d_v_pt)) deallocate(current_state%tendp_3d_v_pt)
    if (allocated(current_state%tendp_3d_w_pt)) deallocate(current_state%tendp_3d_w_pt)
    if (allocated(current_state%tend_3d_u_pt)) deallocate(current_state%tend_3d_u_pt)
    if (allocated(current_state%tend_3d_v_pt)) deallocate(current_state%tend_3d_v_pt)
    if (allocated(current_state%tend_3d_w_pt)) deallocate(current_state%tend_3d_w_pt)

    if (allocated(current_state%tendp_pr_tot_u_pt)) deallocate(current_state%tendp_pr_tot_u_pt)
    if (allocated(current_state%tendp_pr_tot_v_pt)) deallocate(current_state%tendp_pr_tot_v_pt)
    if (allocated(current_state%tendp_pr_tot_w_pt)) deallocate(current_state%tendp_pr_tot_w_pt)
    if (allocated(current_state%tend_pr_tot_u_pt)) deallocate(current_state%tend_pr_tot_u_pt)
    if (allocated(current_state%tend_pr_tot_v_pt)) deallocate(current_state%tend_pr_tot_v_pt)
    if (allocated(current_state%tend_pr_tot_w_pt)) deallocate(current_state%tend_pr_tot_w_pt)

  end subroutine finalisation_callback_pstep

  !> Called each timestep, this will step the pressure field for the non halo columns
  !! @param current_state The current model state
  subroutine timestep_callback_pstep(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: local_y, local_x, target_x_index, target_y_index
    logical :: calculate_diagnostics
    real :: modulo_number

    !calculate_diagnostics = current_state%diagnostic_sample_timestep &
     !                       .and. .not. current_state%halo_column

    local_y=current_state%column_local_y
    local_x=current_state%column_local_x
    target_y_index=local_y-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=local_x-current_state%local_grid%halo_size(X_INDEX)

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tendp_pr_tot_u) then
        current_state%tendp_pr_tot_u_pt(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tendp_pr_tot_v) then
        current_state%tendp_pr_tot_v_pt(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tendp_pr_tot_w) then
        current_state%tendp_pr_tot_w_pt(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_u) then
        current_state%tend_pr_tot_u_pt(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_v) then
        current_state%tend_pr_tot_v_pt(:)= 0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_w) then
        current_state%tend_pr_tot_w_pt(:)= 0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals

    if (current_state%galilean_transformation) call perform_galilean_transformation(current_state, &
         current_state%column_local_y, current_state%column_local_x)

    if ((current_state%modulo_number_3d .eq. 0) &
        .and. (.not. current_state%halo_column)) then
      call save_precomponent_tendencies(current_state, local_x, local_y)!, target_x_index, target_y_index)
    end if

    if (.not. current_state%halo_column) call step_pressure_field(current_state)

    if ((current_state%modulo_number_3d .eq. 0) &
        .and. (.not. current_state%halo_column)) then
      call compute_component_tendencies(current_state, local_x, local_y)!, target_x_index, target_y_index)
    end if

  end subroutine timestep_callback_pstep


  !> Does the actual stepping of the pressure field
  !! @param current_state The current model state
  subroutine step_pressure_field(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, x_index, y_index
    real(kind=DEFAULT_PRECISION) :: dtmtmp

    x_index=current_state%column_local_x
    y_index=current_state%column_local_y

    dtmtmp=merge(current_state%dtm, 0.5_DEFAULT_PRECISION*current_state%dtm, current_state%field_stepping == CENTRED_STEPPING)

    if(current_state%immersed%ib_enabled)then
      if(current_state%immersed%ib_col(y_index, x_index))then
      do k=2,current_state%local_grid%size(Z_INDEX)
        if (current_state%immersed%indic_u(k,y_index,x_index).eq.0.and.&
            current_state%immersed%ghost_u(k,y_index,x_index).lt.1)then 
        current_state%zu%data(k, y_index, x_index)= current_state%zu%data(k, y_index, x_index)+ 2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%horizontal%cx*dtmtmp*(current_state%p%data(k, y_index, x_index)-&
             current_state%p%data(k, y_index, x_index+1))
        end if
        if (current_state%immersed%indic_v(k,y_index,x_index).eq.0.and.&
            current_state%immersed%ghost_v(k,y_index,x_index).lt.1)then 
          current_state%zv%data(k, y_index, x_index)=&
             current_state%zv%data(k, y_index, x_index)+2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%horizontal%cy*dtmtmp*&
             (current_state%p%data(k, y_index, x_index) - current_state%p%data(k, y_index+1, x_index))
        end if
       if (current_state%immersed%indic_w(k,y_index,x_index).eq.0.and.&
           current_state%immersed%ghost_w(k,y_index,x_index).lt.1)then 
          if (k .lt. current_state%local_grid%size(Z_INDEX)) then
          current_state%zw%data(k, y_index, x_index)=current_state%zw%data(k, y_index, x_index)+2.0_DEFAULT_PRECISION*&
               current_state%global_grid%configuration%vertical%rdzn(k+1)*dtmtmp*(current_state%p%data(k, y_index, x_index)-&
               current_state%p%data(k+1, y_index, x_index))
          end if
        end if
      end do
      end if
  
    else 
      ! IB not active
      do k=2,current_state%local_grid%size(Z_INDEX)
        current_state%zu%data(k, y_index, x_index)= current_state%zu%data(k, y_index, x_index)+ 2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%horizontal%cx*dtmtmp*(current_state%p%data(k, y_index, x_index)-&
             current_state%p%data(k, y_index, x_index+1))
        current_state%zv%data(k, y_index, x_index)=&
             current_state%zv%data(k, y_index, x_index)+2.0_DEFAULT_PRECISION*&
             current_state%global_grid%configuration%horizontal%cy*dtmtmp*&
             (current_state%p%data(k, y_index, x_index) - current_state%p%data(k, y_index+1, x_index))
        if (k .lt. current_state%local_grid%size(Z_INDEX)) then
          current_state%zw%data(k, y_index, x_index)=current_state%zw%data(k, y_index, x_index)+2.0_DEFAULT_PRECISION*&
               current_state%global_grid%configuration%vertical%rdzn(k+1)*dtmtmp*(current_state%p%data(k, y_index, x_index)-&
               current_state%p%data(k+1, y_index, x_index))
        end if
      end do
      if (current_state%use_viscosity_and_diffusion .and. current_state%use_surface_boundary_conditions) then
        current_state%zu%data(1, y_index, x_index)=-current_state%zu%data(2, y_index, x_index)-&
             2.0_DEFAULT_PRECISION*current_state%ugal
        current_state%zv%data(1, y_index, x_index)=-current_state%zv%data(2, y_index, x_index)-&
             2.0_DEFAULT_PRECISION*current_state%vgal
      else
        current_state%zu%data(1, y_index, x_index)=current_state%zu%data(2, y_index, x_index)
        current_state%zv%data(1, y_index, x_index)=current_state%zv%data(2, y_index, x_index)
      end if

    end if ! check for IB / no IB
  end subroutine step_pressure_field  

  !> Performs Galilean transformation of flow current and z fields.
  !! @param current_state The current model state
  !! @param y_index Local y index which we are working with
  !! @param x_index Local x index which we are working with
  subroutine perform_galilean_transformation(current_state, y_index, x_index)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: y_index, x_index

    integer :: k

    do k=1,current_state%local_grid%size(Z_INDEX)
      current_state%u%data(k, y_index, x_index)= current_state%u%data(k, y_index, x_index)-current_state%ugal
      current_state%zu%data(k, y_index, x_index)= current_state%zu%data(k, y_index, x_index)-current_state%ugal
      current_state%v%data(k, y_index, x_index)= current_state%v%data(k, y_index, x_index)-current_state%vgal
      current_state%zv%data(k, y_index, x_index)= current_state%zv%data(k, y_index, x_index)-current_state%vgal
    end do
  end subroutine perform_galilean_transformation


   !> Save the 3d tendencies coming into the component.
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine save_precomponent_tendencies(current_state, cxn, cyn)!, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  cxn, cyn!, txn, tyn

    ! Save 3d tendency fields upon request (of 3d or profiles) and availability
    if (l_tendp_3d_u) then
      current_state%tendp_3d_u_pt(:,cyn,cxn)=current_state%zu%data(:,cyn,cxn)
    endif
    if (l_tendp_3d_v) then
      current_state%tendp_3d_v_pt(:,cyn,cxn)=current_state%zv%data(:,cyn,cxn)
    endif
    if (l_tendp_3d_w) then
      current_state%tendp_3d_w_pt(:,cyn,cxn)=current_state%zw%data(:,cyn,cxn)
    endif
    if (l_tend_3d_u) then
      current_state%tend_3d_u_pt(:,cyn,cxn)=current_state%su%data(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      current_state%tend_3d_v_pt(:,cyn,cxn)=current_state%sv%data(:,cyn,cxn)
    endif
    if (l_tend_3d_w) then
      current_state%tend_3d_w_pt(:,cyn,cxn)=current_state%sw%data(:,cyn,cxn)
    endif

  end subroutine save_precomponent_tendencies


   !> Computation of component tendencies
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine compute_component_tendencies(current_state, cxn, cyn)!, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  cxn, cyn!, txn, tyn
    real(kind=DEFAULT_PRECISION) :: dtmtmp

    dtmtmp=2.0_DEFAULT_PRECISION* &
         merge(current_state%dtm, 0.5_DEFAULT_PRECISION*current_state%dtm, current_state%field_stepping == CENTRED_STEPPING)

    ! Calculate change in tendency due to component
    if (l_tendp_3d_u) then
      current_state%tendp_3d_u_pt(:,cyn,cxn)=(current_state%zu%data(:,cyn,cxn) - &
      current_state%tendp_3d_u_pt(:,cyn,cxn))/dtmtmp
    endif
    if (l_tendp_3d_v) then
      current_state%tendp_3d_v_pt(:,cyn,cxn)=(current_state%zv%data(:,cyn,cxn) - &
      current_state%tendp_3d_v_pt(:,cyn,cxn))/dtmtmp
    endif
    if (l_tendp_3d_w) then
      current_state%tendp_3d_w_pt(:,cyn,cxn)=(current_state%zw%data(:,cyn,cxn) - &
      current_state%tendp_3d_w_pt(:,cyn,cxn))/dtmtmp
    endif
    if (l_tend_3d_u) then
      current_state%tend_3d_u_pt(:,cyn,cxn)=current_state%tend_3d_u_pt(:,cyn,cxn) + &
      current_state%tendp_3d_u_pt(:,cyn,cxn)
    endif
    if (l_tend_3d_v) then
      current_state%tend_3d_v_pt(:,cyn,cxn)=current_state%tend_3d_v_pt(:,cyn,cxn) + &
      current_state%tendp_3d_v_pt(:,cyn,cxn)
    endif
    if (l_tend_3d_w) then
      current_state%tend_3d_w_pt(:,cyn,cxn)=current_state%tend_3d_w_pt(:,cyn,cxn) + &
      current_state%tendp_3d_w_pt(:,cyn,cxn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tendp_pr_tot_u) then
      current_state%tendp_pr_tot_u_pt(:)=current_state%tendp_pr_tot_u_pt(:) + &
      current_state%tendp_3d_u_pt(:,cyn,cxn)
    endif
    if (l_tendp_pr_tot_v) then
      current_state%tendp_pr_tot_v_pt(:)=current_state%tendp_pr_tot_v_pt(:) + &
      current_state%tendp_3d_v_pt(:,cyn,cxn)
    endif
    if (l_tendp_pr_tot_w) then
      current_state%tendp_pr_tot_w_pt(:)=current_state%tendp_pr_tot_w_pt(:) + &
      current_state%tendp_3d_w_pt(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_u) then
      current_state%tend_pr_tot_u_pt(:)=current_state%tend_pr_tot_u_pt(:) + &
      current_state%tend_3d_u_pt(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_v) then
      current_state%tend_pr_tot_v_pt(:)=current_state%tend_pr_tot_v_pt(:) + &
      current_state%tend_3d_v_pt(:,cyn,cxn)
    endif
    if (l_tend_pr_tot_w) then
      current_state%tend_pr_tot_w_pt(:)=current_state%tend_pr_tot_w_pt(:) + &
      current_state%tend_3d_w_pt(:,cyn,cxn)
    endif

  end subroutine compute_component_tendencies

end module pstep_mod
