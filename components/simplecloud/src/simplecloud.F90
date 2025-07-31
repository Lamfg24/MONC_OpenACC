!> A very simple saturation adjustment scheme without any microphysics
module simplecloud_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type, &
      component_descriptor_type_v1
  use optionsdatabase_mod, only : options_get_real, options_get_integer
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, X_INDEX, Y_INDEX
  use science_constants_mod, only : r_over_cp, rlvap_over_cp
  use saturation_mod, only: qsaturation, dqwsatdt
  use q_indices_mod, only: get_q_index, standard_q_names
  use logging_mod, only : LOG_ERROR, log_master_log

implicit none

#ifndef TEST_MODE
  private
#endif
    
  ! Indices for vapour and cloud
  integer :: iqv, iql

  integer :: k_cloudmax ! max k index for height
  real(kind=DEFAULT_PRECISION) :: max_height_cloud

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
                   tend_3d_th,      tend_3d_qv,      tend_3d_ql,      tend_3d_tabs
  logical ::     l_tend_3d_th,    l_tend_3d_qv,    l_tend_3d_ql,    l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
               tend_pr_tot_th,  tend_pr_tot_qv,  tend_pr_tot_ql,  tend_pr_tot_tabs
  logical :: l_tend_pr_tot_th,l_tend_pr_tot_qv,l_tend_pr_tot_ql,l_tend_pr_tot_tabs

  public initialisation_callback_simplecloud, timestep_callback_simplecloud, finalisation_callback_simplecloud

contains

  subroutine initialisation_callback_simplecloud(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, iter ! look counter
    logical :: l_qdiag, logicnum
    real(kind=DEFAULT_PRECISION) :: realnum

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "casim_enabled") then
          read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum) then
          call log_master_log(LOG_ERROR, "Casim and Simplecloud are enabled, this does not work yet. Please disable one")
        end if
      else if (current_state%options_database_string(iter,1) .eq. "max_height_cloud") then
        read(current_state%options_database_string(iter,2),*) realnum
        max_height_cloud = realnum
      end if
    end do
    !if (is_component_enabled(current_state%options_database, "casim")) then
    !  call log_master_log(LOG_ERROR, "Casim and Simplecloud are enabled, this does not work yet. Please disable one")
    !end if

    !iqv=get_q_index(standard_q_names%VAPOUR, 'simplecloud')
    !iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'simplecloud')

    ! set buoyancy coefficient (value for vapour should be set
    ! elsewhere for a moist model
    if (.not. allocated(current_state%cq))then
      allocate(current_state%cq(current_state%number_q_fields))
      current_state%cq=0.0_DEFAULT_PRECISION
    end if
    current_state%cq(1) = -1.0

    !max_height_cloud=options_get_real(current_state%options_database, "max_height_cloud")
    do k=2, current_state%local_grid%size(Z_INDEX)-1
      if (current_state%global_grid%configuration%vertical%zn(k) > max_height_cloud) exit
    end do
    k_cloudmax=k-1

    ! Set tendency diagnostic logicals based on availability
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated
    !      in the case where profiles are available
    l_qdiag =  (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0)

    l_tend_pr_tot_th  = current_state%th%active
    l_tend_pr_tot_qv  = l_qdiag
    l_tend_pr_tot_ql  = l_qdiag
    l_tend_pr_tot_tabs  = l_tend_pr_tot_th

    l_tend_3d_th  = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_qv  = l_qdiag .or. l_tend_pr_tot_qv
    l_tend_3d_ql  = l_qdiag .or. l_tend_pr_tot_ql
    l_tend_3d_tabs = l_tend_3d_th

    ! Allocate 3d tendency fields upon availability
    if (l_tend_3d_th) then
      allocate( tend_3d_th(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_qv) then
      allocate( tend_3d_qv(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_ql) then
      allocate( tend_3d_ql(current_state%local_grid%size(Z_INDEX),  &
                           current_state%local_grid%size(Y_INDEX),  &
                           current_state%local_grid%size(X_INDEX)   )    )
    endif
    if (l_tend_3d_tabs) then
      allocate( tend_3d_tabs(current_state%local_grid%size(Z_INDEX),  &
                             current_state%local_grid%size(Y_INDEX),  &
                             current_state%local_grid%size(X_INDEX)   )    )
    endif

    ! Allocate profile tendency fields upon availability
    if (l_tend_pr_tot_th) then
      allocate( tend_pr_tot_th(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_qv) then
      allocate( tend_pr_tot_qv(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_ql) then
      allocate( tend_pr_tot_ql(current_state%local_grid%size(Z_INDEX)) )
    endif
    if (l_tend_pr_tot_tabs) then
      allocate( tend_pr_tot_tabs(current_state%local_grid%size(Z_INDEX)) )
    endif

  end subroutine initialisation_callback_simplecloud


  subroutine finalisation_callback_simplecloud(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(tend_3d_th)) deallocate(tend_3d_th)
    if (allocated(tend_3d_qv)) deallocate(tend_3d_qv)
    if (allocated(tend_3d_ql)) deallocate(tend_3d_ql)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_th)) deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_qv)) deallocate(tend_pr_tot_qv)
    if (allocated(tend_pr_tot_ql)) deallocate(tend_pr_tot_ql)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)

  end subroutine finalisation_callback_simplecloud


  !> Called for each column per timestep this will apply a forcing term 
  !> to the aerosol fields
  !! @param current_state The current model state
  subroutine timestep_callback_simplecloud(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(DEFAULT_PRECISION) :: TdegK   ! Temperature in Kelvin
    real(DEFAULT_PRECISION) :: Pmb     ! Pressure in mb
    real(DEFAULT_PRECISION) :: exner   ! Exner pressure
    real(DEFAULT_PRECISION) :: one_over_exner ! Reciprocal of Exner pressure
    real(DEFAULT_PRECISION) :: qv,qc   ! Shorthand for vapour and cloud mass mixing ratio
    real(DEFAULT_PRECISION) :: qs      ! Saturation mixing ratio
    real(DEFAULT_PRECISION) :: dqsdT   ! Rate of change of qs with temperature
    real(DEFAULT_PRECISION) :: qsatfac ! Multiplicative factor
    real(DEFAULT_PRECISION) :: dmass   ! Mass transfer mixing ratio
    
    integer :: k          ! Loop counter
    integer :: icol, jcol ! Shorthand column indices

    real(DEFAULT_PRECISION) :: dtm  ! Local timestep variable
    integer :: target_x_index, target_y_index
    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep

    ! Zero profile tendency totals on first instance in the sum
    if (current_state%first_timestep_column) then
      if (l_tend_pr_tot_th) then
        tend_pr_tot_th(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_qv) then
        tend_pr_tot_qv(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_ql) then
        tend_pr_tot_ql(:)=0.0_DEFAULT_PRECISION
      endif
      if (l_tend_pr_tot_tabs) then
        tend_pr_tot_tabs(:)=0.0_DEFAULT_PRECISION
      endif
    endif  ! zero totals


    if (current_state%halo_column) return


    dtm = current_state%dtm*2.0
    if (current_state%field_stepping == FORWARD_STEPPING) dtm=current_state%dtm! Should this be revised to scalar_stepping

    icol=current_state%column_local_x
    jcol=current_state%column_local_y
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)

    if (calculate_diagnostics) call save_precomponent_tendencies(current_state, icol, jcol, target_x_index, target_y_index)

    do k=2,k_cloudmax

      exner = current_state%global_grid%configuration%vertical%rprefrcp(k)
      one_over_exner = current_state%global_grid%configuration%vertical%prefrcp(k)
      Pmb   = (current_state%global_grid%configuration%vertical%prefn(k)/100.) 

      if (current_state%field_stepping == FORWARD_STEPPING) then  ! Should this be revised to scalar_stepping
        qv    = current_state%qv%data(k, jcol, icol) + current_state%sqv%data(k, jcol, icol)*dtm
        qc    = current_state%ql%data(k, jcol, icol) + current_state%sql%data(k, jcol, icol)*dtm
        TdegK = (current_state%th%data(k, jcol, icol) + current_state%sth%data(k, jcol, icol)*dtm &
           + current_state%global_grid%configuration%vertical%thref(k))*exner
      else
        qv    = current_state%zqv%data(k, jcol, icol) + current_state%sqv%data(k, jcol, icol)*dtm
        qc    = current_state%zql%data(k, jcol, icol) + current_state%sql%data(k, jcol, icol)*dtm
        TdegK = (current_state%zth%data(k, jcol, icol) + current_state%sth%data(k, jcol, icol)*dtm &
           + current_state%global_grid%configuration%vertical%thref(k))*exner
      end if
      ! Calculate the cloud/vapour increments
      !print*,"k = ",k,"qv = ",qv
      !print*,"k = ",k,"Pmb = ",Pmb
      qs = qsaturation(TdegK, Pmb)
      !print*,"k = ",k,"qs = ",qs
      if (qv > qs .or. qc >0.0)then
        dqsdT = dqwsatdt(qs, TdegK)
        
        qsatfac = 1.0/(1.0 + rlvap_over_cp*dqsdT)
      
        dmass = MAX (-qc,(qv-qs)*qsatfac)/dtm
      
        current_state%sqv%data(k, jcol, icol) = current_state%sqv%data(k, jcol, icol) - dmass
        current_state%sql%data(k, jcol, icol) = current_state%sql%data(k, jcol, icol) + dmass

        current_state%sth%data(k, jcol, icol) = current_state%sth%data(k, jcol, icol) &
           + rlvap_over_cp*dmass*one_over_exner

      end if

    end do

    ! If there's any cloud above then evaporate it
    do k=k_cloudmax+1, current_state%local_grid%size(Z_INDEX)
      if (current_state%scalar_stepping == FORWARD_STEPPING) then
        qv    = current_state%qv%data(k, jcol, icol) + current_state%sq(iqv)%data(k, jcol, icol)*dtm
        qc    = current_state%ql%data(k, jcol, icol) + current_state%sq(iql)%data(k, jcol, icol)*dtm
        TdegK = (current_state%th%data(k, jcol, icol) + current_state%sth%data(k, jcol, icol)*dtm &
           + current_state%global_grid%configuration%vertical%thref(k))*exner
      else
        qv    = current_state%zqv%data(k, jcol, icol) + current_state%sqv%data(k, jcol, icol)*dtm
        qc    = current_state%zql%data(k, jcol, icol) + current_state%sql%data(k, jcol, icol)*dtm
        TdegK = (current_state%zth%data(k, jcol, icol) + current_state%sth%data(k, jcol, icol)*dtm &
           + current_state%global_grid%configuration%vertical%thref(k))*exner
      end if
      if (qc >0.0)then
        dmass = -qc/dtm
      
        current_state%sqv%data(k, jcol, icol) = current_state%sqv%data(k, jcol, icol) - dmass
        current_state%sql%data(k, jcol, icol) = current_state%sql%data(k, jcol, icol) + dmass

        current_state%sth%data(k, jcol, icol) = current_state%sth%data(k, jcol, icol) &
           + rlvap_over_cp*dmass*one_over_exner

      end if
    end do

    if (calculate_diagnostics) call compute_component_tendencies(current_state, icol, jcol, target_x_index, target_y_index)

  end subroutine timestep_callback_simplecloud


   !> Save the 3d tendencies coming into the component.
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine save_precomponent_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(in) :: current_state
    integer, intent(in) ::  cxn, cyn, txn, tyn

    ! Save 3d tendency fields upon request (of 3d or profiles) and availability
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)
    endif
    if (l_tend_3d_qv) then
      tend_3d_qv(:,tyn,txn)=current_state%sq(iqv)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_ql) then
      tend_3d_ql(:,tyn,txn)=current_state%sq(iql)%data(:,cyn,cxn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)=current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)
    endif

  end subroutine save_precomponent_tendencies


   !> Computation of component tendencies
  !! @param current_state Current model state
  !! @param cxn The current slice, x, index
  !! @param cyn The current column, y, index.
  !! @param txn target_x_index
  !! @param tyn target_y_index
  subroutine compute_component_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) ::  cxn, cyn, txn, tyn

    ! Calculate change in tendency due to component
    if (l_tend_3d_th) then
      tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)     - tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_3d_qv) then
      tend_3d_qv(:,tyn,txn)=current_state%sq(iqv)%data(:,cyn,cxn) - tend_3d_qv(:,tyn,txn)
    endif
    if (l_tend_3d_ql) then
      tend_3d_ql(:,tyn,txn)=current_state%sq(iql)%data(:,cyn,cxn) - tend_3d_ql(:,tyn,txn)
    endif
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)=   &
       current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)   &
        - tend_3d_tabs(:,tyn,txn)
    endif

   ! Add local tendency fields to the profile total
    if (l_tend_pr_tot_th) then
      tend_pr_tot_th(:)=tend_pr_tot_th(:) + tend_3d_th(:,tyn,txn)
    endif
    if (l_tend_pr_tot_qv) then
      tend_pr_tot_qv(:)=tend_pr_tot_qv(:) + tend_3d_qv(:,tyn,txn)
    endif
    if (l_tend_pr_tot_ql) then
      tend_pr_tot_ql(:)=tend_pr_tot_ql(:) + tend_3d_ql(:,tyn,txn)
    endif
    if (l_tend_pr_tot_tabs) then
      tend_pr_tot_tabs(:)=tend_pr_tot_tabs(:) + tend_3d_tabs(:,tyn,txn)
    endif

  end subroutine compute_component_tendencies

end module simplecloud_mod
