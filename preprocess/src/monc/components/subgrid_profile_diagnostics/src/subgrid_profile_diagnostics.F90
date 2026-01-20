










module subgrid_profile_diagnostics_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type, &
       component_descriptor_type_v1
  use grids_mod, only : vertical_grid_configuration_type, Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real, options_get_integer, &
       options_get_logical
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use q_indices_mod, only: get_q_index, standard_q_names
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_master_log
  use smagorinsky_mod, only: calculate_half_squared_strain_rate, &
                             calculate_richardson_number,        &
                             calculate_thermal_dissipation_rate
  use science_constants_mod, only: G, ratio_mol_wts
  

  implicit none

  private

  integer :: total_points, iqv=0, iql=0, iqr=0, iqi=0, iqs=0,    &
       iqg=0

! These arrays should in due course be allocated conditionally,
! but let's just get it working first.
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       epsth, ssq, elamr_sq, richardson_number,     &
       subgrid_tke
! Constants here provisionally
  real(kind=DEFAULT_PRECISION) :: a2_n, ath2_n, pr_n, ri_crit
  real(kind=DEFAULT_PRECISION) :: qlcrit

  type(vertical_grid_configuration_type) :: vertical_grid

  logical :: l_lem_dissipation_rate = .true.
  
  public initialisation_callback_subgrid_profile_diagnostics, &
         timestep_callback_subgrid_profile_diagnostics

contains

  subroutine initialisation_callback_subgrid_profile_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, iter
    logical :: logicnum
    real(kind=DEFAULT_PRECISION) :: realnum

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "smagorinsky_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (.not. logicnum) then
          call log_master_log(LOG_ERROR, "Subgrid model diags requested but subgrid model not enabled - check config")
        end if
      else if (current_state%options_database_string(iter,1) .eq. "l_lem_dissipation_rate") then
        read(current_state%options_database_string(iter,2),*) logicnum
        l_lem_dissipation_rate = logicnum
      end if
    end do

    vertical_grid=current_state%global_grid%configuration%vertical

    total_points=(current_state%local_grid%size(Y_INDEX) * current_state%local_grid%size(X_INDEX))
    
    allocate(current_state%uwsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   current_state%vwsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   current_state%uusg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   current_state%vvsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   current_state%wwsg_tot(current_state%local_grid%size(Z_INDEX))   &
         ,   current_state%tkesg_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   current_state%wkesg_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   current_state%vis_coef_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   current_state%diff_coef_tot(current_state%local_grid%size(Z_INDEX)) &
         ,   current_state%richardson_number_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   current_state%richardson_squared_tot(current_state%local_grid%size(Z_INDEX)) &
         ,   current_state%dissipation_tot(current_state%local_grid%size(Z_INDEX) ) &
         ,   current_state%ssub_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   current_state%sed_tot(current_state%local_grid%size(Z_INDEX))  &
         ,   current_state%buoysg_tot(current_state%local_grid%size(Z_INDEX)) )

!   Check consistency of any changes with the deallocations later.
    allocate(ssq(current_state%local_grid%size(Z_INDEX)))
    allocate(elamr_sq(current_state%local_grid%size(Z_INDEX)))
    allocate(richardson_number(current_state%local_grid%size(Z_INDEX)))
    allocate(subgrid_tke(current_state%local_grid%size(Z_INDEX)))
    ! subgrid tke 2d scalar field, included in this component to prevent copying
    ! of code. Need to revisit this
    allocate(current_state%subke_2d(current_state%local_grid%size(Y_INDEX),current_state%local_grid%size(X_INDEX)))
    if (current_state%th%active) then
       allocate(epsth(current_state%local_grid%size(Z_INDEX))         &
            ,   current_state%theta_dis_tot(current_state%local_grid%size(Z_INDEX)) &
            ,   current_state%wtsg_tot(current_state%local_grid%size(Z_INDEX))      &
            ,   current_state%th2sg_tot(current_state%local_grid%size(Z_INDEX)))
    endif
        
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then 
       allocate(current_state%wqsg_tot(current_state%local_grid%size(Z_INDEX)))
       iqv= 1 !get_q_index(standard_q_names%VAPOUR, 'subgrid_profile_diags')
       iql= 2 !get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'subgrid_profile_diags')
       do iter = 1,current_state%config_args
          if (current_state%options_database_string(iter,1) .eq. "qlcrit") then
            read(current_state%options_database_string(iter,2),*) realnum
            qlcrit = realnum
          end if
       end do
       allocate(current_state%wqv_sg_tot(current_state%local_grid%size(Z_INDEX)),   &
            current_state%wql_sg_tot(current_state%local_grid%size(Z_INDEX)))
    endif
    if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
      if (current_state%rain_water_mixing_ratio_index > 0) then
        iqr = current_state%rain_water_mixing_ratio_index
        allocate(current_state%wqr_sg_tot(current_state%local_grid%size(Z_INDEX)))
      endif
      if (current_state%ice_water_mixing_ratio_index > 0) then
        iqi = current_state%ice_water_mixing_ratio_index
        allocate(current_state%wqi_sg_tot(current_state%local_grid%size(Z_INDEX)))
      endif
      if (current_state%snow_water_mixing_ratio_index > 0) then
        iqs = current_state%snow_water_mixing_ratio_index
        allocate(current_state%wqs_sg_tot(current_state%local_grid%size(Z_INDEX)))
      endif
      if (current_state%graupel_water_mixing_ratio_index > 0) then
        iqg = current_state%graupel_water_mixing_ratio_index
        allocate(current_state%wqg_sg_tot(current_state%local_grid%size(Z_INDEX)))
      endif
    end if

    current_state%subke_2d = 0.0_DEFAULT_PRECISION ! to avoid NaN issues
   
  end subroutine initialisation_callback_subgrid_profile_diagnostics

  subroutine timestep_callback_subgrid_profile_diagnostics(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, i, j
    integer :: jcol, icol, target_x_index, target_y_index
    integer :: top_index
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: &
      S11, S22, S33, S12, S23, S13, S33_on_p
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: &
      tau11,tau22, tau33, tau12, tau23, tau13, tau33_on_p
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: u_i_prime_tau_i
    real(kind=DEFAULT_PRECISION) :: vistmp, vis12, vis12a, visonp2, visonp2a, vis13, vis23
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX)) :: umean, wu_umean, vmean, wv_vmean, &
         w_pprime_at_w, rke1, w_qvprime_at_w, w_qclprime_at_w, w_thprime_at_w, wq, rho, rec_rho, buoy_cof, &
         uw_tot, vw_tot,qvprime_at_w,qclprime_at_w,thprime_at_w, sed_eq13, sed_eq23, sed33,uwsg,vwsg,sed_swap
    real(kind=DEFAULT_PRECISION) :: C_virtual
    real(kind=DEFAULT_PRECISION) :: upr_at_w,vpr_at_w, vprime_w_local, uprime_w_local, U_1, U_1_bar
    real(kind=DEFAULT_PRECISION) :: buoy_prod, sg_shear_prod, dissipation_rate
    
    logical :: use_Ri_for_buoyant_prod=.TRUE.

    if ((current_state%modulo_number_0d .ne. 0) .or. (current_state%modulo_number_1d .ne. 0)) return
    !if (.not. current_state%diagnostic_sample_timestep) return

    C_virtual = (ratio_mol_wts-1.0_DEFAULT_PRECISION)
    jcol=current_state%column_local_y
    icol=current_state%column_local_x
    target_y_index=jcol-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=icol-current_state%local_grid%halo_size(X_INDEX)

      if (current_state%first_timestep_column) then
       current_state%sed_tot(:)      = 0.0_DEFAULT_PRECISION
       current_state%ssub_tot(:)     = 0.0_DEFAULT_PRECISION
       current_state%buoysg_tot(:)   = 0.0_DEFAULT_PRECISION
       sed_eq23(:)     = 0.0_DEFAULT_PRECISION
       sed_eq13(:)     = 0.0_DEFAULT_PRECISION
       sed_swap(:)     = 0.0_DEFAULT_PRECISION
       current_state%uwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       current_state%vwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       current_state%uusg_tot(:)     = 0.0_DEFAULT_PRECISION
       current_state%vvsg_tot(:)     = 0.0_DEFAULT_PRECISION
       current_state%wwsg_tot(:)     = 0.0_DEFAULT_PRECISION
       current_state%tkesg_tot(:)    = 0.0_DEFAULT_PRECISION
       current_state%vis_coef_tot(:) = 0.0_DEFAULT_PRECISION
       current_state%diff_coef_tot(:)= 0.0_DEFAULT_PRECISION
       current_state%dissipation_tot(:)      = 0.0_DEFAULT_PRECISION
       current_state%wkesg_tot(:)    = 0.0_DEFAULT_PRECISION
       current_state%richardson_number_tot(:) = 0.0_DEFAULT_PRECISION
       current_state%richardson_squared_tot(:) = 0.0_DEFAULT_PRECISION
       current_state%subke_2d(:,:)   = 0.0_DEFAULT_PRECISION
       if (current_state%th%active) then
         current_state%wtsg_tot(:)     = 0.0_DEFAULT_PRECISION
         current_state%th2sg_tot(:)    = 0.0_DEFAULT_PRECISION
         current_state%theta_dis_tot(:)= 0.0_DEFAULT_PRECISION
       endif
       if (.not. current_state%passive_q .and. &
         current_state%number_q_fields .gt. 0) then
         current_state%wqsg_tot(:)     = 0.0_DEFAULT_PRECISION
         current_state%wqv_sg_tot(:)=0.0_DEFAULT_PRECISION
         current_state%wql_sg_tot(:)=0.0_DEFAULT_PRECISION
         if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
          if (iqr > 0) current_state%wqr_sg_tot(:) = 0.0_DEFAULT_PRECISION
          if (iqi > 0) current_state%wqi_sg_tot(:) = 0.0_DEFAULT_PRECISION
          if (iqs > 0) current_state%wqs_sg_tot(:) = 0.0_DEFAULT_PRECISION
          if (iqg > 0) current_state%wqg_sg_tot(:) = 0.0_DEFAULT_PRECISION
         end if
       endif
    end if

      if (.not. current_state%halo_column) then
    
       do k=2,current_state%local_grid%size(Z_INDEX)-1

! Consistent with ssq calculation from calculate_half_squared_strain_rate
! It could usefully be moved into the Smagorinsky component as a subroutine.
       
      ! Average over 2 points W but 0.5 cancels 2
         S11(k)=current_state%global_grid%configuration%horizontal%cx*(&
           (current_state%u%data(k+1, jcol, icol)-&
           current_state%u%data(k+1, jcol, icol-1))+&
           (current_state%u%data(k, jcol, icol)-&
           current_state%u%data(k, jcol, icol-1)))
      ! Average over 2 points W but 0.5 cancels 2
         S22(k)=current_state%global_grid%configuration%horizontal%cy*(&
           (current_state%v%data(k+1, jcol, icol)-&
           current_state%v%data(k+1, jcol-1, icol))+&
           (current_state%v%data(k, jcol, icol)-&
           current_state%v%data(k, jcol-1, icol)))
      ! Average over 2 points W but 0.5 cancels 2
         S33(k)=((current_state%w%data(k, jcol, icol)-&
           current_state%w%data(k-1, jcol, icol))*&
           current_state%global_grid%configuration%vertical%rdz(k)) +&
           ((current_state%w%data(k+1, jcol, icol)-&
           current_state%w%data(k, jcol, icol))*&
           current_state%global_grid%configuration%vertical%rdz(k+1))
           
      ! No average needed
         S33_on_p(k) = 2.0_DEFAULT_PRECISION * &
         ((current_state%w%data(k,   jcol, icol)-&
           current_state%w%data(k-1, jcol, icol))*&
           current_state%global_grid%configuration%vertical%rdz(k))
      ! Average over 2 points U and W
         S13(k)=(((current_state%u%data(k+1, jcol, icol)-&
                   current_state%u%data(k, jcol, icol))*&
                  current_state%global_grid%configuration%vertical%rdzn(k+1)+&
                  (current_state%w%data(k, jcol, icol+1)-&
                   current_state%w%data(k, jcol, icol))*&
                  current_state%global_grid%configuration%horizontal%cx)+&
                 ((current_state%u%data(k+1, jcol, icol-1)-&
                   current_state%u%data(k, jcol, icol-1))*&
                  current_state%global_grid%configuration%vertical%rdzn(k+1)+&
                  (current_state%w%data(k, jcol, icol)-&
                   current_state%w%data(k, jcol, icol-1))*&
                  current_state%global_grid%configuration%horizontal%cx))*0.5_DEFAULT_PRECISION
      ! Average over 2 points W and V
         S23(k)=(((current_state%w%data(k, jcol, icol) - &
                   current_state%w%data(k, jcol-1, icol)) * &
                  current_state%global_grid%configuration%horizontal%cy + &
                  (current_state%v%data(k+1, jcol-1, icol) - &
                   current_state%v%data(k, jcol-1, icol)) * &
                  current_state%global_grid%configuration%vertical%rdzn(k+1)) + &
                 ((current_state%w%data(k, jcol+1, icol) - &
                   current_state%w%data(k, jcol, icol)) * &
                  current_state%global_grid%configuration%horizontal%cy + &
                  (current_state%v%data(k+1, jcol, icol)-&
                   current_state%v%data(k, jcol, icol))*&
                  current_state%global_grid%configuration%vertical%rdzn(k+1)))*0.5_DEFAULT_PRECISION
      ! Average over 8 points from U and V
         S12(k)=(((((current_state%u%data(k, jcol, icol-1)-&
           current_state%u%data(k, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k, jcol-1, icol)-&
           current_state%v%data(k, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx) +&
           ((current_state%u%data(k, jcol+1, icol-1)-&
           current_state%u%data(k, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k, jcol, icol)-&
           current_state%v%data(k, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx)) +(&
           ((current_state%u%data(k, jcol, icol)-&
           current_state%u%data(k, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k, jcol-1, icol+1)-&
           current_state%v%data(k, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cx) +&
           ((current_state%u%data(k, jcol+1, icol)-&
           current_state%u%data(k, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k, jcol, icol+1)-&
           current_state%v%data(k, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cx)))+((&
           ((current_state%u%data(k+1, jcol, icol-1)-&
           current_state%u%data(k+1, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k+1, jcol-1, icol)-&
           current_state%v%data(k+1, jcol-1, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx)+&
           ((current_state%u%data(k+1, jcol+1, icol-1)-&
           current_state%u%data(k+1, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k+1, jcol, icol)-&
           current_state%v%data(k+1, jcol, icol-1))*&
           current_state%global_grid%configuration%horizontal%cx))+(&
           ((current_state%u%data(k+1, jcol, icol)-&
           current_state%u%data(k+1, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k+1, jcol-1, icol+1)-&
           current_state%v%data(k+1, jcol-1, icol))*&
           current_state%global_grid%configuration%horizontal%cx)+&
           ((current_state%u%data(k+1, jcol+1, icol)-&
           current_state%u%data(k+1, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cy+&
           (current_state%v%data(k+1, jcol, icol+1)-&
           current_state%v%data(k+1, jcol, icol))*&
           current_state%global_grid%configuration%horizontal%cx))))*0.125_DEFAULT_PRECISION

       enddo   
       
      ! Average over 2 points U and W
       S13(1)=(current_state%u%data(2, jcol, icol) + &
               current_state%u%data(2, jcol, icol-1)) * &
               0.5_DEFAULT_PRECISION / current_state%global_grid%configuration%vertical%zn(2)
      ! Average over 2 points W and V
       S23(1)=(current_state%v%data(2, jcol, icol) + &
               current_state%v%data(2, jcol-1, icol)) * &
               0.5_DEFAULT_PRECISION / current_state%global_grid%configuration%vertical%zn(2)
       S12(1)=0.0_DEFAULT_PRECISION
    
       do k=1,current_state%local_grid%size(Z_INDEX)-1
         tau11(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S11(k)
         tau22(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S22(k)
         tau33(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S33(k)
         tau12(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S12(k)
         tau13(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S13(k)
         tau23(k) = current_state%global_grid%configuration%vertical%rho(k) * &
           current_state%vis_coefficient%data(k,jcol,icol) * &
           S23(k)
       enddo 
         
       do k=2,current_state%local_grid%size(Z_INDEX)
         tau33_on_p(k) = current_state%global_grid%configuration%vertical%rhon(k) * 0.5 *&
           (current_state%vis_coefficient%data(k-1,jcol,icol) + &
            current_state%vis_coefficient%data(k,  jcol,icol)) * &
           S33_on_p(k)
       enddo   

!      Subgrid diagnostics
       do k=1, current_state%local_grid%size(Z_INDEX)-1
!
!         Field on kth rho-level, so gradient evaluated from k+1st and kth rhon-levels. Horizontally aligned with
!         w-points.
          current_state%uwsg_tot(k) = current_state%uwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,jcol,icol-1) * &
            ( current_state%u%data(k+1,jcol,icol-1) - &
              current_state%u%data(k,jcol,icol-1) ) * &
            vertical_grid%rdzn(k+1)
          current_state%uwsg_tot(k) = current_state%uwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,jcol,icol) * &
            ( current_state%u%data(k+1,jcol,icol) - &
              current_state%u%data(k,jcol,icol) ) * &
            vertical_grid%rdzn(k+1)
!
!         Field on kth rho-level, so gradient evaluated from k+1st and kth rhon-levels. Horizontally aligned with
!         w-points.
          current_state%vwsg_tot(k) = current_state%vwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,jcol-1,icol) * &
            ( current_state%v%data(k+1,jcol-1,icol) - &
              current_state%v%data(k,jcol-1,icol) ) * &
            vertical_grid%rdzn(k+1)
          current_state%vwsg_tot(k) = current_state%vwsg_tot(k) - 0.5 * &
            current_state%vis_coefficient%data(k,jcol,icol) * &
            ( current_state%v%data(k+1,jcol,icol) - &
              current_state%v%data(k,jcol,icol) ) * &
            vertical_grid%rdzn(k+1)
          !
          ! viscosity and diffusion coefficients
          current_state%vis_coef_tot(k) = current_state%vis_coef_tot(k) + &
               current_state%vis_coefficient%data(k,jcol,icol)
          current_state%diff_coef_tot(k) = current_state%diff_coef_tot(k) + &
               current_state%diff_coefficient%data(k,jcol,icol)
       enddo
       !         Field on kth rho-level, so gradient evaluated from k+1st and kth rhon-levels. Horizontally aligned with
       !         w-points.
       if (current_state%th%active) then
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             current_state%wtsg_tot(k) = current_state%wtsg_tot(k) -  &
                  current_state%diff_coefficient%data(k,jcol,icol) * &
                  ( current_state%th%data(k+1,jcol,icol) + &
                  vertical_grid%thref(k+1) - &
                  current_state%th%data(k,jcol,icol) - &
                  vertical_grid%thref(k) ) * &
                  vertical_grid%rdzn(k+1)
             !
          enddo
       endif

       if (current_state%th%active .and. .not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
          
          do k=1, current_state%local_grid%size(Z_INDEX)-1
             current_state%wqv_sg_tot(k) = current_state%wqv_sg_tot(k) + &
                  current_state%diff_coefficient%data(k, jcol, icol)* &
                  vertical_grid%rdzn(k+1) * &
                  (current_state%qv%data(k,jcol,icol) - current_state%qv%data(k+1,jcol,icol))
             current_state%wql_sg_tot(k) = current_state%wql_sg_tot(k) + &
                  current_state%diff_coefficient%data(k, jcol, icol)* &
                  vertical_grid%rdzn(k+1) * &
                  (current_state%ql%data(k,jcol,icol) - current_state%ql%data(k+1,jcol,icol))
          enddo
          if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
            if (iqr > 0) then
              do k=1, current_state%local_grid%size(Z_INDEX)-1
                  current_state%wqr_sg_tot(k) = current_state%wqr_sg_tot(k) + &
                      current_state%diff_coefficient%data(k, jcol, icol)* &
                      vertical_grid%rdzn(k+1) * &
                      (current_state%qr%data(k,jcol,icol) - current_state%qr%data(k+1,jcol,icol))
              enddo
            endif
            if (iqi > 0) then
              do k=1, current_state%local_grid%size(Z_INDEX)-1
                  current_state%wqi_sg_tot(k) = current_state%wqi_sg_tot(k) + &
                      current_state%diff_coefficient%data(k, jcol, icol)* &
                      vertical_grid%rdzn(k+1) * &
                      (current_state%qi%data(k,jcol,icol) - current_state%qi%data(k+1,jcol,icol))
              enddo
            endif
            if (iqs > 0) then
              do k=1, current_state%local_grid%size(Z_INDEX)-1
                  current_state%wqs_sg_tot(k) = current_state%wqs_sg_tot(k) + &
                      current_state%diff_coefficient%data(k, jcol, icol)* &
                      vertical_grid%rdzn(k+1) * &
                      (current_state%qs%data(k,jcol,icol) - current_state%qs%data(k+1,jcol,icol))
              enddo
            endif
            if (iqg > 0) then
              do k=1, current_state%local_grid%size(Z_INDEX)-1
                  current_state%wqg_sg_tot(k) = current_state%wqg_sg_tot(k) + &
                      current_state%diff_coefficient%data(k, jcol, icol)* &
                      vertical_grid%rdzn(k+1) * &
                      (current_state%qg%data(k,jcol,icol) - current_state%qg%data(k+1,jcol,icol))
              enddo
            endif
          endif
       endif

! Level 1 wind speed - needed in various places

       U_1=SQRT(current_state%u%data(2,jcol,icol) * &
                current_state%u%data(2,jcol,icol) + &
                current_state%v%data(2,jcol,icol) * &
                current_state%v%data(2,jcol,icol))
       U_1_bar=SQRT(current_state%global_grid%configuration%vertical%olubar(2)*&
                    current_state%global_grid%configuration%vertical%olubar(2)+&
                    current_state%global_grid%configuration%vertical%olvbar(2)*&
                    current_state%global_grid%configuration%vertical%olvbar(2))

      ! Do p levels - calculate subgrid TKE diagnostics based on tau, interpolate to w levels
       do k=1, current_state%local_grid%size(Z_INDEX)-1
         rho(k) = current_state%global_grid%configuration%vertical%rho(k)
       end do

! Buoyancy coefficient on u (zn) levels
! note: this is not official MONC elsewhere
       do k=1, current_state%local_grid%size(Z_INDEX)
         buoy_cof(k)=G/current_state%global_grid%configuration%vertical%thref(k) 
       end do


!      =======================================================
!      Calculation of subgrid diagnostics dependent on the dissipation.
!
!      Rationalize this with smagorinsky.
       ri_crit = 0.25_DEFAULT_PRECISION
!
       ssq=calculate_half_squared_strain_rate(current_state, current_state%u, current_state%v, current_state%w)
       richardson_number=calculate_richardson_number(current_state, ssq, current_state%th, current_state%qv, &
       current_state%ql)
!
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          !        Mimic LEM: I think this is the square of the mixing length (Brown et al. 94, Eq 7)
          elamr_sq(k) = 0.0

          if ( ( &
               current_state%vis_coefficient%data(k, jcol, icol) &
               > 0.0) .and. (richardson_number(k) < ri_crit) ) then
             elamr_sq(k) = &
                  current_state%vis_coefficient%data(k, jcol, icol) / &
                  sqrt( 1.0 - richardson_number(k) * &
                  current_state%diff_coefficient%data(k, jcol, icol) / &
                  current_state%vis_coefficient%data(k, jcol, icol) )
          end if
          if (l_lem_dissipation_rate) then
             dissipation_rate = ssq(k) * ( &
                  current_state%vis_coefficient%data(k, jcol, icol) - &
                  richardson_number(k) * &
                  current_state%diff_coefficient%data(k, jcol, icol) )
          else ! CIRCLE-A code 
             if (use_Ri_for_buoyant_prod) then
      
                buoy_prod = -ssq(k) * richardson_number(k) *  &
                     current_state%diff_coefficient%data(k, jcol, icol)
                
             else
                if (current_state%th%active .and. &
                     .not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
                   buoy_prod =  &
                        -current_state%diff_coefficient%data(k,jcol,icol) * &
                        !d/dz - so inside needs to be on p points
                        current_state%global_grid%configuration%vertical%rdzn(k+1) * &
                        !G/th_ref on p
                        (buoy_cof(k+1) * &  ! 1
                        !(th')
                        (current_state%th%data(k+1,jcol,icol) + & ! 2
                        current_state%global_grid%configuration%vertical%thref(k+1) + & ! 2
                        !th_ref * (C_virtual* qv'(k)-qcl'(k))
                        current_state%global_grid%configuration%vertical%thref(k+1) * & ! 2
                        (C_virtual * current_state%qv%data(k+1,jcol,icol) -  & ! 3
                        current_state%ql%data(k+1,jcol,icol))) - & ! 1
                        !diff for d/dz
                        buoy_cof(k) * & ! 1
                        (current_state%th%data(k,jcol,icol) + & ! 2
                        current_state%global_grid%configuration%vertical%thref(k) + & ! 2
                        current_state%global_grid%configuration%vertical%thref(k) * & ! 2
                        (C_virtual * current_state%qv%data(k,jcol,icol) - & ! 3
                        current_state%ql%data(k,jcol,icol))) ) ! 0
                else
                   call log_master_log(LOG_ERROR, &
                        "Subgrid diags - buoy_prod calc needs theta and q active, STOP")
                endif
             endif
             ! Circle-A dissipation rate includes buoy_prod term, which is different
             ! to the original LEM version
             dissipation_rate = ssq(k) * &
                  current_state%vis_coefficient%data(k, jcol, icol) + &
                  buoy_prod

             current_state%buoysg_tot(k)=current_state%buoysg_tot(k) + buoy_prod

          endif
          
          current_state%dissipation_tot(k) =  current_state%dissipation_tot(k) + dissipation_rate
          
          ! PAC comment: Not clear the code below is correct (where doea a2_n fit in).
          !
          !        This constant needs to go somewhere accessible and where it can be set. The code in the LEM
          !        seems to go round and round in circles here.
          a2_n = 0.23
          subgrid_tke(k) = ( dissipation_rate * dissipation_rate * elamr_sq(k) ) ** (1.0/3.0) / a2_n
          
          !        Assume isotropy.
          current_state%tkesg_tot(k) = current_state%tkesg_tot(k) + subgrid_tke(k)
          current_state%uusg_tot(k) = current_state%uusg_tot(k) + (2.0/3.0) * subgrid_tke(k)
          current_state%vvsg_tot(k) = current_state%vvsg_tot(k) + (2.0/3.0) * subgrid_tke(k)
          current_state%wwsg_tot(k) = current_state%wwsg_tot(k) + (2.0/3.0) * subgrid_tke(k)
          current_state%wkesg_tot(k) = current_state%wkesg_tot(k) + subgrid_tke(k) * &
               current_state%w%data(k,jcol,icol)
          current_state%richardson_number_tot(k) = current_state%richardson_number_tot(k) +  richardson_number(k)
          current_state%richardson_squared_tot(k) = current_state%richardson_squared_tot(k) + &
          richardson_number(k)**2.0_DEFAULT_PRECISION
          
          ! Calculate the integrated subgrid TKE for the scalar.
          ! BAD place to put this but convenient
          current_state%subke_2d(target_y_index, target_x_index) = current_state%subke_2d(target_y_index, target_x_index) + &
               (subgrid_tke(k)*0.5_DEFAULT_PRECISION*                                           &
               vertical_grid%dzn(k+1)*                       &
               vertical_grid%rho(k))
       enddo
       ! *********************** Subgrid buoyant production***************************    
       ! Note - calculating on z levels (i.e. w) 
       ! So need dzn, rdzn.
       ! zn(k)-zn(k-1) is dzn(k), so nothing is stored in k=1
       !------
       if (.not. l_lem_dissipation_rate) then
          if (current_state%th%active .and. &
               .not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
             k = 1
             buoy_prod =  &
                  -current_state%diff_coefficient%data(k,jcol,icol) * &
                  !d/dz - so inside needs to be on p points
                  current_state%global_grid%configuration%vertical%rdzn(k+1) * &
                  !G/th_ref on p
                  (buoy_cof(k+1) * &  ! 1
                  !(th')
                  (current_state%th%data(k+1,jcol,icol) + & ! 2
                  current_state%global_grid%configuration%vertical%thref(k+1) + & ! 2
                  !th_ref * (C_virtual* qv'(k)-qcl'(k))
                  current_state%global_grid%configuration%vertical%thref(k+1) * & ! 2
                  (C_virtual * current_state%qv%data(k+1,jcol,icol)   & ! 3
                  - current_state%ql%data(k+1,jcol,icol))) - & ! 1
                  !diff for d/dz
                  buoy_cof(k) * & ! 1
                  (current_state%th%data(k,jcol,icol) + & ! 2
                  current_state%global_grid%configuration%vertical%thref(k) + & ! 2
                  current_state%global_grid%configuration%vertical%thref(k) * & ! 2
                  (C_virtual * current_state%qv%data(k,jcol,icol) & ! 3
                  - current_state%ql%data(k,jcol,icol))) ) ! 0
             
             current_state%buoysg_tot(k)=current_state%buoysg_tot(k) + buoy_prod

          
             current_state%dissipation_tot(1) = current_state%dissipation_tot(1) + &
                  current_state%vis_coefficient%data(1, jcol, icol) &
                  * U_1 * U_1 / &
                  (current_state%global_grid%configuration%vertical%zn(2) * &
                  current_state%global_grid%configuration%vertical%zn(2)) + buoy_prod
          else
             call log_master_log(LOG_ERROR, &
                "Subgrid diags - buoy_prod calc needs theta and q active, STOP")
          endif
       endif
          
       ! Needs buoyancy part added later
       
       ! divide subke by altitude to make it column mean (as in the LEM)
       current_state%subke_2d(target_y_index, target_x_index) =  current_state%subke_2d(target_y_index, target_x_index)/ &
            vertical_grid%z(current_state%local_grid%size(Z_INDEX))
       !
       !        Thermal equivalents:
       !        Again, this needs to go somewhere better. pr_n is used in Smagorinsky anyway (as a local variable).
       
       if (current_state%th%active) then
          epsth=calculate_thermal_dissipation_rate(current_state, current_state%th)         
          ath2_n = 0.3
          pr_n   = 0.7
          do k=2, current_state%local_grid%size(Z_INDEX)-1
             ! Theta dissipation rate
             current_state%theta_dis_tot(k) = current_state%theta_dis_tot(k) + epsth(k)
             if (subgrid_tke(k) > 0.0) &
                  current_state%th2sg_tot(k) = current_state%th2sg_tot(k) + sqrt( a2_n * elamr_sq(k) / subgrid_tke(k) ) * &
                  epsth(k) / ( ath2_n**2 * pr_n)
          enddo
       endif
       
       
       ! *********************** Subgrid shear production***************************    
       ! Note - calculating on z levels (i.e. w) 
       !
       do k=2,current_state%local_grid%size(Z_INDEX)-2
          
          
          !Subgrid shear-------
          umean(k)=(current_state%global_grid%configuration%vertical%olubar(k+1) - &
               current_state%global_grid%configuration%vertical%olubar(k)) * &
               current_state%global_grid%configuration%vertical%rdzn(k+1)
          vmean(k)=(current_state%global_grid%configuration%vertical%olvbar(k+1) -  &
               current_state%global_grid%configuration%vertical%olvbar(k)) * &
               current_state%global_grid%configuration%vertical%rdzn(k+1)
          
          sg_shear_prod = current_state%vis_coefficient%data(k,jcol,icol)* &
               (S13(k)*umean(k)+S23(k)*vmean(k))
          
          current_state%ssub_tot(k)=current_state%ssub_tot(k) + sg_shear_prod
          
          !tau is at p points - so convert to w points
          uwsg(k)=0.5*(tau13(k)+tau13(k+1))/rho(k)
          vwsg(k)=0.5*(tau23(k)+tau23(k+1))/rho(k)
          
       end do
       
       k=1
             
       sg_shear_prod = current_state%vis_coefficient%data(1, jcol, icol) * &
            ( 0.5 * & 
            (current_state%u%data(2,jcol,icol-1)  + & 
            current_state%u%data(2,jcol,icol) ) * &
            current_state%global_grid%configuration%vertical%olubar(2) + &
            0.5 * & 
            (current_state%v%data(2,jcol-1,icol)  + & 
            current_state%v%data(2,jcol,icol) ) * &
            current_state%global_grid%configuration%vertical%olvbar(2) ) / &
            (current_state%global_grid%configuration%vertical%zn(2) * &
            current_state%global_grid%configuration%vertical%zn(2)) 
       
       current_state%ssub_tot(1)=current_state%ssub_tot(1) + sg_shear_prod
       
       ! *********************** Subgrid turbulence transport***************************    
       ! Note - calculating on z levels (i.e. w) 
       ! so need  u_i_prime_tau_i on p levels
       
       do k=2, current_state%local_grid%size(Z_INDEX)-1
          
          u_i_prime_tau_i(k) = ( 0.5_DEFAULT_PRECISION * &
               (current_state%u%data(k,jcol,icol-1)  + &
               current_state%u%data(k,jcol,icol)  ) - &
               current_state%global_grid%configuration%vertical%olubar(k)) * &
               0.5_DEFAULT_PRECISION * ( tau13(k-1)+tau13(k))
          
          u_i_prime_tau_i(k) = u_i_prime_tau_i(k) + ( 0.5_DEFAULT_PRECISION * &
               (current_state%v%data(k,jcol-1,icol)  + &
               current_state%v%data(k,jcol,icol)  ) - &
               current_state%global_grid%configuration%vertical%olvbar(k)) * &
               0.5_DEFAULT_PRECISION * ( tau23(k-1)+tau23(k))
          
          u_i_prime_tau_i(k) = u_i_prime_tau_i(k) + ( 0.5_DEFAULT_PRECISION * &
               (current_state%w%data(k-1,jcol,icol)+&
               current_state%w%data(k,jcol,icol))) *&
               tau33_on_p(k)
          
       end do
       
       u_i_prime_tau_i(1) = &
            ( 0.5_DEFAULT_PRECISION  * & 
            (current_state%u%data(2,jcol,icol-1)  + & 
            current_state%u%data(2,jcol,icol) ) - &
           current_state%global_grid%configuration%vertical%olubar(2) ) * &
           ( 0.5_DEFAULT_PRECISION  * &
           (current_state%u%data(2,jcol,icol-1)  + & 
           current_state%u%data(2,jcol,icol) ) )
       
       u_i_prime_tau_i(1) = u_i_prime_tau_i(1) + &
            ( 0.5_DEFAULT_PRECISION  * & 
            (current_state%v%data(2,jcol-1,icol)  + & 
            current_state%v%data(2,jcol,icol) ) - &
            current_state%global_grid%configuration%vertical%olvbar(2) ) * &
         ( 0.5_DEFAULT_PRECISION  * &
         (current_state%v%data(2,jcol-1,icol)  + & 
         current_state%v%data(2,jcol,icol) ) )
       
       ! Note - here we use the surface vis_coefficient to ensure exact cancellation in
       ! numerical budget
       u_i_prime_tau_i(1) = u_i_prime_tau_i(1) * &        
            current_state%vis_coefficient%data(1,jcol,icol) / & 
            current_state%global_grid%configuration%vertical%zn(2)  
       
       
       do k=2, current_state%local_grid%size(Z_INDEX)-2
          current_state%sed_tot(k)=current_state%sed_tot(k) + (u_i_prime_tau_i(k+1)-u_i_prime_tau_i(k)) * &
               current_state%global_grid%configuration%vertical%rdzn(k+1) / &
               current_state%global_grid%configuration%vertical%rho(k)
       end do
       
       current_state%sed_tot(1)=current_state%sed_tot(1) + u_i_prime_tau_i(1) / &
            current_state%global_grid%configuration%vertical%zn(2)
       
       !      =======================================================
    endif ! (.not. current_state%halo_column)
  end subroutine timestep_callback_subgrid_profile_diagnostics
  
end module subgrid_profile_diagnostics_mod
