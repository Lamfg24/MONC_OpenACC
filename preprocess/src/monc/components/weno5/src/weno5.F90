










!> WENO5 advection scheme
module weno5_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type, &
      component_descriptor_type_v1
  use optionsdatabase_mod, only : options_get_string, options_get_integer
  use collections_mod, only : map_type
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use q_indices_mod, only: get_q_index, standard_q_names
  use prognostics_mod, only : prognostic_field_type

implicit none

  private

  logical :: advect_flow, advect_th, advect_q

  logical :: l_toplevel=.true.

  ! Local tendency diagnostic variables for this component
  ! 3D tendency fields and logicals for their use
  logical :: l_tend_3d_u, l_tend_3d_v, l_tend_3d_w, l_tend_3d_th,l_tend_3d_qv,       &
             l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,       &
             l_tend_3d_tabs
  ! Local mean tendency profile fields and logicals for their use
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w, l_tend_pr_tot_th,l_tend_pr_tot_qv,       &
             l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,l_tend_pr_tot_qs,l_tend_pr_tot_qg,       &
             l_tend_pr_tot_tabs
  ! q indices
  integer :: iqv=1, iql=2, iqr=3, iqi=4, iqs=5, iqg=6

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(DEFAULT_PRECISION) :: gamma1_weno3 = 1.0/3.0
  real(DEFAULT_PRECISION) :: gamma2_weno3 = 2.0/3.0


  real(DEFAULT_PRECISION) :: gamma1 = 1.0/10.0
  real(DEFAULT_PRECISION) :: gamma2 = 2.0/5.0
  real(DEFAULT_PRECISION) :: gamma3 = 3.0/10.0

  real(DEFAULT_PRECISION) :: sign_term = 0.5

  real(DEFAULT_PRECISION) :: epsilon_weno = 1.0E-6


 public initialisation_callback_weno5advection, timestep_callback_weno5advection, &
        finalisation_callback_weno5advection
contains


  !> Initialisation callback, will set up the configuration of this advection scheme
  !! @param current_state The current model state
  subroutine initialisation_callback_weno5advection(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    logical :: l_qdiag
    integer :: iter
    integer :: alloc_z, alloc_y, alloc_x

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "advection_flow_fields") then
        advect_flow=determine_if_advection_here(current_state%options_database_string(iter,2))
      else if (current_state%options_database_string(iter,1) .eq. "advection_theta_field") then
        advect_th=determine_if_advection_here(current_state%options_database_string(iter,2))
      else if (current_state%options_database_string(iter,1) .eq. "advection_q_fields") then
        advect_q=determine_if_advection_here(current_state%options_database_string(iter,2))
      end if
    end do

    ! Set tendency diagnostic logicals based on availability 
    ! Need to use 3d tendencies to compute the profiles, so they will be allocated 
    !      in the case where profiles are available
    !alloc_z=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    !alloc_y=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    !alloc_x=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2


  end subroutine initialisation_callback_weno5advection


  subroutine finalisation_callback_weno5advection(current_state)
    type(model_state_type), target, intent(inout) :: current_state


  end subroutine finalisation_callback_weno5advection


  !> Called per column of data, this will perform Piacsek-Williams advection on the applicable fields for non halo data
  !! @param current_state The current model state
  !! @param target_(x/y)_index This is the index with the halos subtracted. This is needed so that diagnostic does 
  !!                           not include halos and to prevent array out-of-bounds
  subroutine timestep_callback_weno5advection(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: current_x_index, current_y_index, target_x_index, target_y_index, k
    logical :: calculate_diagnostics

    real(DEFAULT_PRECISION) :: mean_u_cell_Lx, mean_u_cell_Rx, mean_v_cell_Ly, mean_v_cell_Ry, &
                               mean_w_cell_Bz, mean_w_cell_Tz, mean_rho_cell_Bz, mean_rho_cell_Tz
    real(DEFAULT_PRECISION) :: xu_flux_left_face, xu_flux_right_face, yu_flux_left_face, yu_flux_right_face, &
                               zu_flux_bottom_face, zu_flux_top_face
    real(DEFAULT_PRECISION) :: u_flux_total
    real(DEFAULT_PRECISION) :: xv_flux_left_face, xv_flux_right_face, yv_flux_left_face, yv_flux_right_face, &
                               zv_flux_bottom_face, zv_flux_top_face
    real(DEFAULT_PRECISION) :: v_flux_total
    real(DEFAULT_PRECISION) :: xw_flux_left_face, xw_flux_right_face, yw_flux_left_face, yw_flux_right_face, &
                               zw_flux_bottom_face, zw_flux_top_face
    real(DEFAULT_PRECISION) :: w_flux_total
    real(DEFAULT_PRECISION) :: xth_flux_left_face, xth_flux_right_face, yth_flux_left_face, yth_flux_right_face, &
                               zth_flux_bottom_face, zth_flux_top_face
    real(DEFAULT_PRECISION) :: th_flux_total
    real(DEFAULT_PRECISION) :: xqv_flux_left_face, xqv_flux_right_face, yqv_flux_left_face, yqv_flux_right_face, &
                               zqv_flux_bottom_face, zqv_flux_top_face
    real(DEFAULT_PRECISION) :: qv_flux_total
    real(DEFAULT_PRECISION) :: xql_flux_left_face, xql_flux_right_face, yql_flux_left_face, yql_flux_right_face, &
                               zql_flux_bottom_face, zql_flux_top_face
    real(DEFAULT_PRECISION) :: ql_flux_total
    real(DEFAULT_PRECISION) :: xqr_flux_left_face, xqr_flux_right_face, yqr_flux_left_face, yqr_flux_right_face, &
                               zqr_flux_bottom_face, zqr_flux_top_face
    real(DEFAULT_PRECISION) :: qr_flux_total
    real(DEFAULT_PRECISION) :: xqi_flux_left_face, xqi_flux_right_face, yqi_flux_left_face, yqi_flux_right_face, &
                               zqi_flux_bottom_face, zqi_flux_top_face
    real(DEFAULT_PRECISION) :: qi_flux_total
    real(DEFAULT_PRECISION) :: xqs_flux_left_face, xqs_flux_right_face, yqs_flux_left_face, yqs_flux_right_face, &
                               zqs_flux_bottom_face, zqs_flux_top_face
    real(DEFAULT_PRECISION) :: qs_flux_total
    real(DEFAULT_PRECISION) :: xqg_flux_left_face, xqg_flux_right_face, yqg_flux_left_face, yqg_flux_right_face, &
                               zqg_flux_bottom_face, zqg_flux_top_face
    real(DEFAULT_PRECISION) :: qg_flux_total
    real(DEFAULT_PRECISION) :: xqAitkenSolMass_flux_left_face, xqAitkenSolMass_flux_right_face, &
                               yqAitkenSolMass_flux_left_face, yqAitkenSolMass_flux_right_face, &
                               zqAitkenSolMass_flux_bottom_face, zqAitkenSolMass_flux_top_face
    real(DEFAULT_PRECISION) :: qAitkenSolMass_flux_total
    real(DEFAULT_PRECISION) :: xqAccumSolMass_flux_left_face, xqAccumSolMass_flux_right_face, &
                               yqAccumSolMass_flux_left_face, yqAccumSolMass_flux_right_face, &
                               zqAccumSolMass_flux_bottom_face, zqAccumSolMass_flux_top_face
    real(DEFAULT_PRECISION) :: qAccumSolMass_flux_total
    real(DEFAULT_PRECISION) :: xqAccumInsolMass_flux_left_face, xqAccumInsolMass_flux_right_face, &
                               yqAccumInsolMass_flux_left_face, yqAccumInsolMass_flux_right_face, &
                               zqAccumInsolMass_flux_bottom_face, zqAccumInsolMass_flux_top_face
    real(DEFAULT_PRECISION) :: qAccumInsolMass_flux_total
    real(DEFAULT_PRECISION) :: xqCoarseSolMass_flux_left_face, xqCoarseSolMass_flux_right_face, &
                               yqCoarseSolMass_flux_left_face, yqCoarseSolMass_flux_right_face, &
                               zqCoarseSolMass_flux_bottom_face, zqCoarseSolMass_flux_top_face
    real(DEFAULT_PRECISION) :: qCoarseSolMass_flux_total
    real(DEFAULT_PRECISION) :: xqCoarseDustMass_flux_left_face, xqCoarseDustMass_flux_right_face, &
                               yqCoarseDustMass_flux_left_face, yqCoarseDustMass_flux_right_face, &
                               zqCoarseDustMass_flux_bottom_face, zqCoarseDustMass_flux_top_face
    real(DEFAULT_PRECISION) :: qCoarseDustMass_flux_total
    real(DEFAULT_PRECISION) :: xnl_flux_left_face, xnl_flux_right_face, ynl_flux_left_face, ynl_flux_right_face, &
                               znl_flux_bottom_face, znl_flux_top_face
    real(DEFAULT_PRECISION) :: nl_flux_total
    real(DEFAULT_PRECISION) :: xnr_flux_left_face, xnr_flux_right_face, ynr_flux_left_face, ynr_flux_right_face, &
                               znr_flux_bottom_face, znr_flux_top_face
    real(DEFAULT_PRECISION) :: nr_flux_total
    real(DEFAULT_PRECISION) :: xni_flux_left_face, xni_flux_right_face, yni_flux_left_face, yni_flux_right_face, &
                               zni_flux_bottom_face, zni_flux_top_face
    real(DEFAULT_PRECISION) :: ni_flux_total
    real(DEFAULT_PRECISION) :: xns_flux_left_face, xns_flux_right_face, yns_flux_left_face, yns_flux_right_face, &
                               zns_flux_bottom_face, zns_flux_top_face
    real(DEFAULT_PRECISION) :: ns_flux_total
    real(DEFAULT_PRECISION) :: xng_flux_left_face, xng_flux_right_face, yng_flux_left_face, yng_flux_right_face, &
                               zng_flux_bottom_face, zng_flux_top_face
    real(DEFAULT_PRECISION) :: ng_flux_total
    real(DEFAULT_PRECISION) :: xnAitkenSolNumber_flux_left_face, xnAitkenSolNumber_flux_right_face, &
                               ynAitkenSolNumber_flux_left_face, ynAitkenSolNumber_flux_right_face, &
                               znAitkenSolNumber_flux_bottom_face, znAitkenSolNumber_flux_top_face
    real(DEFAULT_PRECISION) :: nAitkenSolNumber_flux_total
    real(DEFAULT_PRECISION) :: xnAccumSolNumber_flux_left_face, xnAccumSolNumber_flux_right_face, &
                               ynAccumSolNumber_flux_left_face, ynAccumSolNumber_flux_right_face, &
                               znAccumSolNumber_flux_bottom_face, znAccumSolNumber_flux_top_face
    real(DEFAULT_PRECISION) :: nAccumSolNumber_flux_total
    real(DEFAULT_PRECISION) :: xnAccumInsolNumber_flux_left_face, xnAccumInsolNumber_flux_right_face, &
                               ynAccumInsolNumber_flux_left_face, ynAccumInsolNumber_flux_right_face, &
                               znAccumInsolNumber_flux_bottom_face, znAccumInsolNumber_flux_top_face
    real(DEFAULT_PRECISION) :: nAccumInsolNumber_flux_total
    real(DEFAULT_PRECISION) :: xnCoarseSolNumber_flux_left_face, xnCoarseSolNumber_flux_right_face, &
                               ynCoarseSolNumber_flux_left_face, ynCoarseSolNumber_flux_right_face, &
                               znCoarseSolNumber_flux_bottom_face, znCoarseSolNumber_flux_top_face
    real(DEFAULT_PRECISION) :: nCoarseSolNumber_flux_total
    real(DEFAULT_PRECISION) :: xnCoarseDustnumber_flux_left_face, xnCoarseDustnumber_flux_right_face, &
                               ynCoarseDustnumber_flux_left_face, ynCoarseDustnumber_flux_right_face, &
                               znCoarseDustnumber_flux_bottom_face, znCoarseDustnumber_flux_top_face
    real(DEFAULT_PRECISION) :: nCoarseDustnumber_flux_total
    real(DEFAULT_PRECISION) :: somme, sum_column_sqv = 0.0, sum_column_qv = 0.0
    real(DEFAULT_PRECISION) :: zqv_flux_k2, zqv_flux_k3, zqv_flux_kn
    real(DEFAULT_PRECISION) :: surface, volume

    logical :: top_u, top_v, top_w, top_th, top_qv, top_ql, top_qr, top_qi, top_qs, top_qg, &
               top_qAitkenSolMass, top_qAccumSolMass, top_qAccumInsolMass, top_qCoarseSolMass, top_qCoarseDustMass, &
               top_nl, top_nr, top_ni, top_ns, top_ng, &
               top_nAitkenSolNumber, top_nAccumSolNumber, top_nAccumInsolNumber, top_nCoarseSolNumber, top_nCoarseDustnumber

    !calculate_diagnostics = current_state%diagnostic_sample_timestep

    surface = 60.0*60.0

    current_x_index=current_state%column_local_x
    current_y_index=current_state%column_local_y
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)


    mean_u_cell_Lx = 0.0_DEFAULT_PRECISION
    mean_u_cell_Rx = 0.0_DEFAULT_PRECISION
    mean_v_cell_Ly = 0.0_DEFAULT_PRECISION
    mean_v_cell_Ry = 0.0_DEFAULT_PRECISION
    mean_w_cell_Bz = 0.0_DEFAULT_PRECISION
    mean_w_cell_Tz = 0.0_DEFAULT_PRECISION
    mean_rho_cell_Bz = 0.0_DEFAULT_PRECISION
    mean_rho_cell_Tz = 0.0_DEFAULT_PRECISION

    u_flux_total = 0.0_DEFAULT_PRECISION
    v_flux_total = 0.0_DEFAULT_PRECISION
    w_flux_total = 0.0_DEFAULT_PRECISION
    th_flux_total = 0.0_DEFAULT_PRECISION
    qv_flux_total = 0.0_DEFAULT_PRECISION
    ql_flux_total = 0.0_DEFAULT_PRECISION
    qr_flux_total = 0.0_DEFAULT_PRECISION
    qi_flux_total = 0.0_DEFAULT_PRECISION
    qs_flux_total = 0.0_DEFAULT_PRECISION
    qg_flux_total = 0.0_DEFAULT_PRECISION
    qAitkenSolMass_flux_total = 0.0_DEFAULT_PRECISION
    qAccumSolMass_flux_total = 0.0_DEFAULT_PRECISION
    qAccumInsolMass_flux_total = 0.0_DEFAULT_PRECISION
    qCoarseSolMass_flux_total = 0.0_DEFAULT_PRECISION
    qCoarseDustMass_flux_total = 0.0_DEFAULT_PRECISION
    nl_flux_total = 0.0_DEFAULT_PRECISION
    nr_flux_total = 0.0_DEFAULT_PRECISION
    ni_flux_total = 0.0_DEFAULT_PRECISION
    ns_flux_total = 0.0_DEFAULT_PRECISION
    ng_flux_total = 0.0_DEFAULT_PRECISION
    nAitkenSolNumber_flux_total = 0.0_DEFAULT_PRECISION
    nAccumSolNumber_flux_total = 0.0_DEFAULT_PRECISION
    nAccumInsolNumber_flux_total = 0.0_DEFAULT_PRECISION
    nCoarseSolNumber_flux_total = 0.0_DEFAULT_PRECISION
    nCoarseDustnumber_flux_total = 0.0_DEFAULT_PRECISION
    somme = 0.0_DEFAULT_PRECISION

    top_u = .false.
    top_v = .false.
    top_w = .false.
    top_th = .false.
    top_qv = .false.
    top_ql = .false.
    top_qr = .false.
    top_qi = .false.
    top_qs = .false.
    top_qg = .false.
    top_qAitkenSolMass = .false.
    top_qAccumSolMass = .false.
    top_qAccumInsolMass = .false.
    top_qCoarseSolMass = .false.
    top_qCoarseDustMass = .false.
    top_nl = .false.
    top_nr = .false.
    top_ni = .false.
    top_ns = .false.
    top_ng = .false.
    top_nAitkenSolNumber = .false.
    top_nAccumSolNumber = .false.
    top_nAccumInsolNumber = .false.
    top_nCoarseSolNumber = .false.
    top_nCoarseDustnumber = .false.

    if (current_state%halo_column) return

!     do k = 2, current_state%local_grid%size(Z_INDEX)-1
!
!       mean_u_cell_Lx = 0.5 * (current_state%u%data(k, current_y_index, current_x_index) + &
!                               current_state%u%data(k, current_y_index, current_x_index-1))
!       mean_u_cell_Rx = 0.5 * (current_state%u%data(k, current_y_index, current_x_index) + &
!                               current_state%u%data(k, current_y_index, current_x_index+1))
!
!       mean_v_cell_Ly = 0.5 * (current_state%v%data(k, current_y_index, current_x_index) + &
!                               current_state%v%data(k, current_y_index-1, current_x_index))
!       mean_v_cell_Ry = 0.5 * (current_state%v%data(k, current_y_index, current_x_index) + &
!                               current_state%v%data(k, current_y_index+1, current_x_index))
!
!       mean_w_cell_Bz = 0.5 * (current_state%w%data(k, current_y_index, current_x_index) + &
!                               current_state%w%data(k-1, current_y_index, current_x_index))
!       mean_w_cell_Tz = 0.5 * (current_state%w%data(k, current_y_index, current_x_index) + &
!                               current_state%w%data(k+1, current_y_index, current_x_index))
!
!       mean_rho_cell_Bz = 0.5 * (current_state%global_grid%configuration%vertical%rho(k) + &
!                                current_state%global_grid%configuration%vertical%rho(k-1))
!       mean_rho_cell_Tz = 0.5 * (current_state%global_grid%configuration%vertical%rho(k) + &
!                               current_state%global_grid%configuration%vertical%rho(k-1))
!
!     end do

    advect_flow = .true.
    if (advect_flow) then
    !print*,"advect_flow_weno5"

    ! U component

    do k = 2, current_state%local_grid%size(Z_INDEX)

    u_flux_total = 0.0_DEFAULT_PRECISION
    v_flux_total = 0.0_DEFAULT_PRECISION
    w_flux_total = 0.0_DEFAULT_PRECISION

    top_u = .false.
    top_v = .false.
    top_w = .false.

    somme = 0.0_DEFAULT_PRECISION

    mean_u_cell_Lx = 0.5 * (current_state%u%data(k, current_y_index, current_x_index) + &
                              current_state%u%data(k, current_y_index, current_x_index-1))
    mean_u_cell_Rx = 0.5 * (current_state%u%data(k, current_y_index, current_x_index) + &
                            current_state%u%data(k, current_y_index, current_x_index+1))

    mean_v_cell_Ly = 0.5 * (current_state%v%data(k, current_y_index, current_x_index) + &
                            current_state%v%data(k, current_y_index-1, current_x_index))
    mean_v_cell_Ry = 0.5 * (current_state%v%data(k, current_y_index, current_x_index) + &
                            current_state%v%data(k, current_y_index+1, current_x_index))

    if (k .ge. 2 .and. k .le. (current_state%local_grid%size(Z_INDEX)-1)) then
      mean_w_cell_Bz = 0.5 * (current_state%w%data(k, current_y_index, current_x_index) + &
                              current_state%w%data(k-1, current_y_index, current_x_index))
      mean_w_cell_Tz = 0.5 * (current_state%w%data(k, current_y_index, current_x_index) + &
                              current_state%w%data(k+1, current_y_index, current_x_index))

      mean_rho_cell_Bz = 0.5 * (current_state%global_grid%configuration%vertical%rho(k) + &
                                current_state%global_grid%configuration%vertical%rho(k-1))
      mean_rho_cell_Tz = 0.5 * (current_state%global_grid%configuration%vertical%rho(k) + &
                              current_state%global_grid%configuration%vertical%rho(k+1))
    end if

    volume = surface * current_state%global_grid%configuration%vertical%dz(k)

    call x_left_face_flux_calc(current_state%u, mean_u_cell_Lx, current_x_index, current_y_index, k, xu_flux_left_face)
    call x_right_face_flux_calc(current_state%u, mean_u_cell_Rx, current_x_index, current_y_index, k, xu_flux_right_face)
    somme = (xu_flux_right_face - xu_flux_left_face)*current_state%global_grid%configuration%horizontal%cx ! (du_i+1/2*u - du_i-1/2*u)/delta_x
    if (xu_flux_right_face .eq. xu_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    u_flux_total = u_flux_total + somme

    call y_left_face_flux_calc(current_state%u, mean_v_cell_Ly, current_x_index, current_y_index, k, yu_flux_left_face)
    call y_right_face_flux_calc(current_state%u, mean_v_cell_Ry, current_x_index, current_y_index, k, yu_flux_right_face)
    somme = (yu_flux_right_face - yu_flux_left_face)*current_state%global_grid%configuration%horizontal%cy ! (du_j+1/2*v - du_j-1/2*v)/delta_y
    if (yu_flux_right_face .eq. yu_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    u_flux_total = u_flux_total + somme

    if (k .ge. 1 .and. k .le. 2) then

     call z_bottom_face_flux_calc_weno1(current_state, current_state%u, &
          current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zu_flux_bottom_face)
     call z_top_face_flux_calc_weno1(current_state, current_state%u, current_state%w%data(k, current_y_index, current_x_index), &
                                        current_x_index, current_y_index, k, zu_flux_top_face, top_u)
     somme = (zu_flux_top_face - zu_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

     if (zu_flux_top_face .eq. zu_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
     u_flux_total = u_flux_total + somme

    else if (k .eq. 3) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%u, mean_w_cell_Bz, &
                                         current_x_index, current_y_index, k, zu_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%u, mean_w_cell_Tz,&
                                      current_x_index, current_y_index, k, zu_flux_top_face)
      somme = (zu_flux_top_face - zu_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zu_flux_top_face .eq. zu_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      u_flux_total = u_flux_total + somme

    else if (k .ge. 4 .and. k .le. (current_state%local_grid%size(Z_INDEX)-3)) then

      call z_bottom_face_flux_calc_weno5(current_state, current_state%u, mean_w_cell_Bz, &
                                          current_x_index, current_y_index, k, zu_flux_bottom_face)
      call z_top_face_flux_calc_weno5(current_state, current_state%u, mean_w_cell_Tz, &
                                          current_x_index, current_y_index, k, zu_flux_top_face)
      somme = (zu_flux_top_face - zu_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zu_flux_top_face .eq. zu_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      u_flux_total = u_flux_total + somme

    else if (k .eq. (current_state%local_grid%size(Z_INDEX)-2)) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%u, mean_w_cell_Bz, &
                                         current_x_index, current_y_index, k, zu_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%u, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zu_flux_top_face)
      somme = (zu_flux_top_face - zu_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zu_flux_top_face .eq. zu_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      u_flux_total = u_flux_total + somme

    else if (k .ge. (current_state%local_grid%size(Z_INDEX)-1) .and. k .le. current_state%local_grid%size(Z_INDEX)) then

      if (k .eq. current_state%local_grid%size(Z_INDEX)) top_u = .true.
      call z_bottom_face_flux_calc_weno1(current_state, current_state%u, &
          current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zu_flux_bottom_face)
      call z_top_face_flux_calc_weno1(current_state, current_state%u, current_state%w%data(k, current_y_index, current_x_index), &
                                        current_x_index, current_y_index, k, zu_flux_top_face, top_u)
      somme = (zu_flux_top_face - zu_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zu_flux_top_face .eq. zu_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      u_flux_total = u_flux_total + somme

    end if

    current_state%su%data(k, current_y_index, current_x_index) = &
                current_state%su%data(k, current_y_index, current_x_index) - u_flux_total

    ! V component

    call x_left_face_flux_calc(current_state%v, mean_u_cell_Lx, current_x_index, current_y_index, k, xv_flux_left_face)
    call x_right_face_flux_calc(current_state%v, mean_u_cell_Rx, current_x_index, current_y_index, k, xv_flux_right_face)
    somme = (xv_flux_right_face - xv_flux_left_face)*current_state%global_grid%configuration%horizontal%cx ! (dv_i+1/2*u - dv_i-1/2*u)/delta_x
    if (xv_flux_right_face .eq. xv_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    v_flux_total = v_flux_total + somme

    call y_left_face_flux_calc(current_state%v, mean_v_cell_Ly, current_x_index, current_y_index, k, yv_flux_left_face)
    call y_right_face_flux_calc(current_state%v, mean_v_cell_Ry, current_x_index, current_y_index, k, yv_flux_right_face)
    somme = (yv_flux_right_face - yv_flux_left_face)*current_state%global_grid%configuration%horizontal%cy ! (dv_j+1/2*v - dv_j-1/2*v)/delta_y
    if (yv_flux_right_face .eq. yv_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    v_flux_total = v_flux_total + somme

    if (k .ge. 1 .and. k .le. 2) then

      call z_bottom_face_flux_calc_weno1(current_state, current_state%v, &
          current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zv_flux_bottom_face)
      call z_top_face_flux_calc_weno1(current_state, current_state%v, current_state%w%data(k, current_y_index, current_x_index), &
                                        current_x_index, current_y_index, k, zv_flux_top_face, top_v)
      somme = (zv_flux_top_face - zv_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zv_flux_top_face .eq. zv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      v_flux_total = v_flux_total + somme

    else if (k .eq. 3) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%v, mean_w_cell_Bz, &
                                        current_x_index, current_y_index, k, zv_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%v, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zv_flux_top_face)
      somme = (zv_flux_top_face - zv_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zv_flux_top_face .eq. zv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      v_flux_total = v_flux_total + somme

    else if (k .ge. 4 .and. k .le. (current_state%local_grid%size(Z_INDEX)-3)) then

      call z_bottom_face_flux_calc_weno5(current_state, current_state%v, mean_w_cell_Bz, &
                                      current_x_index, current_y_index, k, zv_flux_bottom_face)
      call z_top_face_flux_calc_weno5(current_state, current_state%v, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zv_flux_top_face)
      somme = (zv_flux_top_face - zv_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)
      if (zv_flux_top_face .eq. zv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      v_flux_total = v_flux_total + somme

    else if (k .eq. (current_state%local_grid%size(Z_INDEX)-2)) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%v, mean_w_cell_Bz, &
                                         current_x_index, current_y_index, k, zv_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%v, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zv_flux_top_face)
      somme = (zv_flux_top_face - zv_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zv_flux_top_face .eq. zv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      v_flux_total = v_flux_total + somme

    else if (k .ge. (current_state%local_grid%size(Z_INDEX)-1) .and. k .le. current_state%local_grid%size(Z_INDEX)) then

      if (k .eq. current_state%local_grid%size(Z_INDEX)) top_v = .true.
      call z_bottom_face_flux_calc_weno1(current_state, current_state%v, &
          current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zv_flux_bottom_face)
      call z_top_face_flux_calc_weno1(current_state, current_state%v, current_state%w%data(k, current_y_index, current_x_index), &
                                         current_x_index, current_y_index, k, zv_flux_top_face, top_v)
      somme = (zv_flux_top_face - zv_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zv_flux_top_face .eq. zv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      v_flux_total = v_flux_total + somme

    end if

    current_state%sv%data(k, current_y_index, current_x_index) = &
                   current_state%sv%data(k, current_y_index, current_x_index) - v_flux_total

    ! W component

    call x_left_face_flux_calc(current_state%w, mean_u_cell_Lx, current_x_index, current_y_index, k, xw_flux_left_face)
    call x_right_face_flux_calc(current_state%w, mean_u_cell_Rx, current_x_index, current_y_index, k, xw_flux_right_face)
    somme = (xw_flux_right_face - xw_flux_left_face)*current_state%global_grid%configuration%horizontal%cx ! (dw_i+1/2*u - dw_i-1/2*u)/delta_x
    if (xw_flux_right_face .eq. xw_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    w_flux_total = w_flux_total + somme

    call y_left_face_flux_calc(current_state%w, mean_v_cell_Ly, current_x_index, current_y_index, k, yw_flux_left_face)
    call y_right_face_flux_calc(current_state%w, mean_v_cell_Ry, current_x_index, current_y_index, k, yw_flux_right_face)
    somme = (yw_flux_right_face - yw_flux_left_face)*current_state%global_grid%configuration%horizontal%cy ! (dw_j+1/2*v - dw_j-1/2*v)/delta_y
    if (yw_flux_right_face .eq. yw_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    w_flux_total = w_flux_total + somme

    if (k .ge. 1 .and. k .le. 2) then

      call z_bottom_face_flux_calc_weno1(current_state, current_state%w, &
          current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zw_flux_bottom_face)
      call z_top_face_flux_calc_weno1(current_state, current_state%w, current_state%w%data(k, current_y_index, current_x_index), &
                                         current_x_index, current_y_index, k, zw_flux_top_face, top_w)
      somme = (zw_flux_top_face - zw_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zw_flux_top_face .eq. zw_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      w_flux_total = w_flux_total + somme

    else if (k .eq. 3) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%w, mean_w_cell_Bz, &
                                         current_x_index, current_y_index, k, zw_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%w, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zw_flux_top_face)
      somme = (zw_flux_top_face - zw_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zw_flux_top_face .eq. zw_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      w_flux_total = w_flux_total + somme

    else if (k .ge. 4 .and. k .le. (current_state%local_grid%size(Z_INDEX)-3)) then

      call z_bottom_face_flux_calc_weno5(current_state, current_state%w, mean_w_cell_Bz, &
                                      current_x_index, current_y_index, k, zw_flux_bottom_face)
      call z_top_face_flux_calc_weno5(current_state, current_state%w, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zw_flux_top_face)
      somme = (zw_flux_top_face - zw_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zw_flux_top_face .eq. zw_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      w_flux_total = w_flux_total + somme

    else if (k .eq. (current_state%local_grid%size(Z_INDEX)-2)) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%w, mean_w_cell_Bz, &
                                      current_x_index, current_y_index, k, zw_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%w, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zw_flux_top_face)
      somme = (zw_flux_top_face - zw_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zw_flux_top_face .eq. zw_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      w_flux_total = w_flux_total + somme

    else if (k .ge. (current_state%local_grid%size(Z_INDEX)-1) .and. k .le. current_state%local_grid%size(Z_INDEX)) then

      if (k .eq. current_state%local_grid%size(Z_INDEX)) top_w = .true.
      call z_bottom_face_flux_calc_weno1(current_state, current_state%w, &
          current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zw_flux_bottom_face)
      call z_top_face_flux_calc_weno1(current_state, current_state%w, current_state%w%data(k, current_y_index, current_x_index), &
                                         current_x_index, current_y_index, k, zw_flux_top_face, top_w)
      somme = (zw_flux_top_face - zw_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)
      if (zw_flux_top_face .eq. zw_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      w_flux_total = w_flux_total + somme

    end if

    current_state%sw%data(k, current_y_index, current_x_index) = &
                  current_state%sw%data(k, current_y_index, current_x_index) - w_flux_total

    end do
    end if


    advect_th = .true.
    if (advect_th) then
    !print*,"advect_th_weno5"

    do k = 2, current_state%local_grid%size(Z_INDEX)

    th_flux_total = 0.0_DEFAULT_PRECISION
    somme = 0.0_DEFAULT_PRECISION

    top_th = .false.

    mean_u_cell_Lx = 0.5 * (current_state%u%data(k, current_y_index, current_x_index) + &
                              current_state%u%data(k, current_y_index, current_x_index-1))
    mean_u_cell_Rx = 0.5 * (current_state%u%data(k, current_y_index, current_x_index) + &
                            current_state%u%data(k, current_y_index, current_x_index+1))

    mean_v_cell_Ly = 0.5 * (current_state%v%data(k, current_y_index, current_x_index) + &
                            current_state%v%data(k, current_y_index-1, current_x_index))
    mean_v_cell_Ry = 0.5 * (current_state%v%data(k, current_y_index, current_x_index) + &
                            current_state%v%data(k, current_y_index+1, current_x_index))

    if (k .ge. 2 .and. k .le. (current_state%local_grid%size(Z_INDEX)-1)) then
      mean_w_cell_Bz = 0.5 * (current_state%w%data(k, current_y_index, current_x_index) + &
                              current_state%w%data(k-1, current_y_index, current_x_index))
      mean_w_cell_Tz = 0.5 * (current_state%w%data(k, current_y_index, current_x_index) + &
                              current_state%w%data(k+1, current_y_index, current_x_index))

      mean_rho_cell_Bz = 0.5 * (current_state%global_grid%configuration%vertical%rho(k) + &
                                current_state%global_grid%configuration%vertical%rho(k-1))
      mean_rho_cell_Tz = 0.5 * (current_state%global_grid%configuration%vertical%rho(k) + &
                              current_state%global_grid%configuration%vertical%rho(k+1))
    end if

    volume = surface * current_state%global_grid%configuration%vertical%dzn(k)

    call x_left_face_flux_calc(current_state%th, mean_u_cell_Lx, current_x_index, current_y_index, k, xth_flux_left_face)
    call x_right_face_flux_calc(current_state%th, mean_u_cell_Rx, current_x_index, current_y_index, k, xth_flux_right_face)
    somme = (xth_flux_right_face - xth_flux_left_face)*current_state%global_grid%configuration%horizontal%cx ! (dth_i+1/2*u - dth_i-1/2*u)/delta_x
    if (xth_flux_right_face .eq. xth_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    th_flux_total = th_flux_total + somme


    call y_left_face_flux_calc(current_state%th, mean_v_cell_Ly, current_x_index, current_y_index, k, yth_flux_left_face)
    call y_right_face_flux_calc(current_state%th, mean_v_cell_Ry, current_x_index, current_y_index, k, yth_flux_right_face)
    somme = (yth_flux_right_face - yth_flux_left_face)*current_state%global_grid%configuration%horizontal%cy ! (dth_j+1/2*v - dth_j-1/2*v)/delta_y
    if (yth_flux_right_face .eq. yth_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    th_flux_total = th_flux_total + somme

    if (k .ge. 1 .and. k .le. 2) then

      call z_bottom_face_flux_calc_weno1(current_state, current_state%th, &
            current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zth_flux_bottom_face)
      call z_top_face_flux_calc_weno1(current_state, current_state%th, current_state%w%data(k, current_y_index, current_x_index), &
                                        current_x_index, current_y_index, k, zth_flux_top_face, top_th)
      somme = (zth_flux_top_face - zth_flux_bottom_face)* current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zth_flux_top_face .eq. zth_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      th_flux_total = th_flux_total + somme

    else if (k .eq. 3) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%th, mean_w_cell_Bz, &
                                      current_x_index, current_y_index, k, zth_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%th, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zth_flux_top_face)
      somme = (zth_flux_top_face - zth_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zth_flux_top_face .eq. zth_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      th_flux_total = th_flux_total + somme

    else if (k .ge. 4 .and. k .le. (current_state%local_grid%size(Z_INDEX)-3)) then

      call z_bottom_face_flux_calc_weno5(current_state, current_state%th, mean_w_cell_Bz, &
                                      current_x_index, current_y_index, k, zth_flux_bottom_face)
      call z_top_face_flux_calc_weno5(current_state, current_state%th, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zth_flux_top_face)
      somme = (zth_flux_top_face - zth_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zth_flux_top_face .eq. zth_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      th_flux_total = th_flux_total + somme

    else if (k .eq. (current_state%local_grid%size(Z_INDEX)-2)) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%th, mean_w_cell_Bz, &
                                      current_x_index, current_y_index, k, zth_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%th, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zth_flux_top_face)
      somme = (zth_flux_top_face - zth_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zth_flux_top_face .eq. zth_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      th_flux_total = th_flux_total + somme

    else if (k .ge. (current_state%local_grid%size(Z_INDEX)-1) .and. k .le. current_state%local_grid%size(Z_INDEX)) then

      if (k .eq. current_state%local_grid%size(Z_INDEX)) top_th = .true.
      call z_bottom_face_flux_calc_weno1(current_state, current_state%th, &
          current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zth_flux_bottom_face)
      call z_top_face_flux_calc_weno1(current_state, current_state%th, current_state%w%data(k, current_y_index, current_x_index), &
                                      current_x_index, current_y_index, k, zth_flux_top_face, top_th)
      somme = (zth_flux_top_face - zth_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zth_flux_top_face .eq. zth_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      th_flux_total = th_flux_total + somme

    end if

    current_state%sth%data(k, current_y_index, current_x_index) = &
                  current_state%sth%data(k, current_y_index, current_x_index) - th_flux_total

    end do
    end if

    advect_q = .false.
    if (advect_q) then
    !print*,"advect_q_weno5"

    do k = 2, current_state%local_grid%size(Z_INDEX)-1

    qv_flux_total = 0.0_DEFAULT_PRECISION
    ql_flux_total = 0.0_DEFAULT_PRECISION
    qr_flux_total = 0.0_DEFAULT_PRECISION
    qi_flux_total = 0.0_DEFAULT_PRECISION
    qs_flux_total = 0.0_DEFAULT_PRECISION
    qg_flux_total = 0.0_DEFAULT_PRECISION
    qAitkenSolMass_flux_total = 0.0_DEFAULT_PRECISION
    qAccumSolMass_flux_total = 0.0_DEFAULT_PRECISION
    qAccumInsolMass_flux_total = 0.0_DEFAULT_PRECISION
    qCoarseSolMass_flux_total = 0.0_DEFAULT_PRECISION
    qCoarseDustMass_flux_total = 0.0_DEFAULT_PRECISION
    nl_flux_total = 0.0_DEFAULT_PRECISION
    nr_flux_total = 0.0_DEFAULT_PRECISION
    ni_flux_total = 0.0_DEFAULT_PRECISION
    ns_flux_total = 0.0_DEFAULT_PRECISION
    ng_flux_total = 0.0_DEFAULT_PRECISION
    nAitkenSolNumber_flux_total = 0.0_DEFAULT_PRECISION
    nAccumSolNumber_flux_total = 0.0_DEFAULT_PRECISION
    nAccumInsolNumber_flux_total = 0.0_DEFAULT_PRECISION
    nCoarseSolNumber_flux_total = 0.0_DEFAULT_PRECISION
    nCoarseDustnumber_flux_total = 0.0_DEFAULT_PRECISION
    somme = 0.0_DEFAULT_PRECISION


    top_qv = .false.
    top_ql = .false.
    top_qr = .false.
    top_qi = .false.
    top_qs = .false.
    top_qg = .false.
    top_qAitkenSolMass = .false.
    top_qAccumSolMass = .false.
    top_qAccumInsolMass = .false.
    top_qCoarseSolMass = .false.
    top_qCoarseDustMass = .false.
    top_nl = .false.
    top_nr = .false.
    top_ni = .false.
    top_ns = .false.
    top_ng = .false.
    top_nAitkenSolNumber = .false.
    top_nAccumSolNumber = .false.
    top_nAccumInsolNumber = .false.
    top_nCoarseSolNumber = .false.
    top_nCoarseDustnumber = .false.

    mean_u_cell_Lx = 0.5 * (current_state%u%data(k, current_y_index, current_x_index) + &
                              current_state%u%data(k, current_y_index, current_x_index-1))
    mean_u_cell_Rx = 0.5 * (current_state%u%data(k, current_y_index, current_x_index) + &
                            current_state%u%data(k, current_y_index, current_x_index+1))

    mean_v_cell_Ly = 0.5 * (current_state%v%data(k, current_y_index, current_x_index) + &
                            current_state%v%data(k, current_y_index-1, current_x_index))
    mean_v_cell_Ry = 0.5 * (current_state%v%data(k, current_y_index, current_x_index) + &
                            current_state%v%data(k, current_y_index+1, current_x_index))

    if (k .ge. 2 .and. k .le. (current_state%local_grid%size(Z_INDEX)-1)) then
      mean_w_cell_Bz = 0.5 * (current_state%w%data(k, current_y_index, current_x_index) + &
                              current_state%w%data(k-1, current_y_index, current_x_index))
      mean_w_cell_Tz = 0.5 * (current_state%w%data(k, current_y_index, current_x_index) + &
                              current_state%w%data(k+1, current_y_index, current_x_index))

      mean_rho_cell_Bz = 0.5 * (current_state%global_grid%configuration%vertical%rho(k) + &
                                current_state%global_grid%configuration%vertical%rho(k-1))
      mean_rho_cell_Tz = 0.5 * (current_state%global_grid%configuration%vertical%rho(k) + &
                              current_state%global_grid%configuration%vertical%rho(k+1))
    end if

    volume = surface * current_state%global_grid%configuration%vertical%dzn(k)


    ! qv

    call x_left_face_flux_calc(current_state%qv, mean_u_cell_Lx, current_x_index, current_y_index, k, xqv_flux_left_face)
    call x_right_face_flux_calc(current_state%qv, mean_u_cell_Rx, current_x_index, current_y_index, k, xqv_flux_right_face)
    somme = (xqv_flux_right_face - xqv_flux_left_face)*current_state%global_grid%configuration%horizontal%cx ! (dth_i+1/2*u - dth_i-1/2*u)/delta_x
    if (xqv_flux_right_face .eq. xqv_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    qv_flux_total = qv_flux_total + somme


    call y_left_face_flux_calc(current_state%qv, mean_v_cell_Ly, current_x_index, current_y_index, k, yqv_flux_left_face)
    call y_right_face_flux_calc(current_state%qv, mean_v_cell_Ry, current_x_index, current_y_index, k, yqv_flux_right_face)
    somme = (yqv_flux_right_face - yqv_flux_left_face)*current_state%global_grid%configuration%horizontal%cy ! (dth_j+1/2*v - dth_j-1/2*v)/delta_y
    if (yqv_flux_right_face .eq. yqv_flux_left_face) somme = 0.0_DEFAULT_PRECISION
    qv_flux_total = qv_flux_total + somme


    if (k .ge. 1 .and. k .le. 2) then

      call z_bottom_face_flux_calc_weno1(current_state, current_state%th, &
            current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zqv_flux_bottom_face)
      call z_top_face_flux_calc_weno1(current_state, current_state%th, current_state%w%data(k, current_y_index, current_x_index), &
                                        current_x_index, current_y_index, k, zqv_flux_top_face, top_th)
      somme = (zqv_flux_top_face - zqv_flux_bottom_face)* current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zqv_flux_top_face .eq. zqv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      qv_flux_total = qv_flux_total + somme

    else if (k .eq. 3) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%th, mean_w_cell_Bz, &
                                      current_x_index, current_y_index, k, zqv_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%th, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zqv_flux_top_face)
      somme = (zqv_flux_top_face - zqv_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zqv_flux_top_face .eq. zqv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      qv_flux_total = qv_flux_total + somme

    else if (k .ge. 4 .and. k .le. (current_state%local_grid%size(Z_INDEX)-3)) then

      call z_bottom_face_flux_calc_weno5(current_state, current_state%th, mean_w_cell_Bz, &
                                      current_x_index, current_y_index, k, zqv_flux_bottom_face)
      call z_top_face_flux_calc_weno5(current_state, current_state%th, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zqv_flux_top_face)
      somme = (zqv_flux_top_face - zqv_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zqv_flux_top_face .eq. zqv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      qv_flux_total = qv_flux_total + somme

    else if (k .eq. (current_state%local_grid%size(Z_INDEX)-2)) then

      call z_bottom_face_flux_calc_weno3(current_state, current_state%th, mean_w_cell_Bz, &
                                      current_x_index, current_y_index, k, zqv_flux_bottom_face)
      call z_top_face_flux_calc_weno3(current_state, current_state%th, mean_w_cell_Tz, &
                                      current_x_index, current_y_index, k, zqv_flux_top_face)
      somme = (zqv_flux_top_face - zqv_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zqv_flux_top_face .eq. zqv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      qv_flux_total = qv_flux_total + somme

    else if (k .ge. (current_state%local_grid%size(Z_INDEX)-1) .and. k .le. current_state%local_grid%size(Z_INDEX)) then

      if (k .eq. current_state%local_grid%size(Z_INDEX)) top_th = .true.
      call z_bottom_face_flux_calc_weno1(current_state, current_state%th, &
          current_state%w%data(k, current_y_index, current_x_index), current_x_index, current_y_index, k, zqv_flux_bottom_face)
      call z_top_face_flux_calc_weno1(current_state, current_state%th, current_state%w%data(k, current_y_index, current_x_index), &
                                      current_x_index, current_y_index, k, zqv_flux_top_face, top_th)
      somme = (zqv_flux_top_face - zqv_flux_bottom_face) * current_state%global_grid%configuration%vertical%rdzn(k)/&
               current_state%global_grid%configuration%vertical%rho(k)

      if (zqv_flux_top_face .eq. zqv_flux_bottom_face) somme = 0.0_DEFAULT_PRECISION
      qv_flux_total = qv_flux_total + somme

    end if

    current_state%sqv%data(k, current_y_index, current_x_index) = &
                  current_state%sqv%data(k, current_y_index, current_x_index) - qv_flux_total



    end do
    end if



    
    !if (current_state%modulo_number_3d .eq. 0) &
    !    call save_precomponent_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)
    !if (advect_flow) call advect_flow_fields(current_state, current_x_index, current_y_index)
    !if (advect_th) call advect_th_field(current_state, current_x_index, current_y_index)
    !if (advect_q) call advect_q_field(current_state, current_x_index, current_y_index)

    !if (current_state%modulo_number_3d .eq. 0) &
    !    call compute_component_tendencies(current_state, current_x_index, current_y_index)!, target_x_index, target_y_index)

  end subroutine timestep_callback_weno5advection

  !Computation at the mass point on grid u(i+1/2,j,k)
  subroutine x_right_face_flux_calc(var, mean_var, current_x_index, current_y_index, k, flux_right_face)
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_right_face

    real(DEFAULT_PRECISION) :: alpha1, alpha2, alpha3, sum_alpha
    real(DEFAULT_PRECISION) :: beta1, beta2, beta3
    real(DEFAULT_PRECISION) :: w1, w2, w3

    real(DEFAULT_PRECISION) :: varL_1, varL_2, varL_3, var_L ! at i+1/2 and L=+
    real(DEFAULT_PRECISION) :: varR_1, varR_2, varR_3, var_R ! at i+1/2 and R=-
    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at i+1/2, terms results of the applying of Lax-Friedrich flux splitting method

    ! Positive fluxes at i+1/2

    !                             i+1/2
    !       i-2 ----- i-1 ----- i --|-- i+1 ----- i+2 ----> x direction (u wind)
    !                              L
    !        |------------------|
    !                   S1
    !                   |----------------|
    !                           S2
    !                           |------------------|
    !                                    S3

    ! Stencils S1 = [i-2, i-1, i], S2 = [i-1, i, i+1], S3 = [i, i+1, i+2]

    varL_1 = (1.0/3.0)*var%data(k, current_y_index, current_x_index-2) - &
             (7.0/6.0)*var%data(k, current_y_index, current_x_index-1) + &
             (11.0/6.0)*var%data(k, current_y_index, current_x_index)

    varL_2 = - (1.0/6.0)*var%data(k, current_y_index, current_x_index-1) + &
               (5.0/6.0)*var%data(k, current_y_index, current_x_index) + &
               (1.0/3.0)*var%data(k, current_y_index, current_x_index+1)

    varL_3 = (1.0/3.0)*var%data(k, current_y_index, current_x_index) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index+1) - &
             (1.0/6.0)*var%data(k, current_y_index, current_x_index+2)

    beta1 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index-2) - &
                        2*var%data(k, current_y_index, current_x_index-1) + &
                        var%data(k, current_y_index, current_x_index))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index-2) - &
                        4*var%data(k, current_y_index, current_x_index-1) + &
                        3*var%data(k, current_y_index, current_x_index))**2

    beta2 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index-1) - &
                        2*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index, current_x_index+1))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index-1) - &
                      var%data(k, current_y_index, current_x_index+1))**2

    beta3 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index) - &
                        2*var%data(k, current_y_index, current_x_index+1) + &
                        var%data(k, current_y_index, current_x_index+2))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index, current_x_index) - &
                        4*var%data(k, current_y_index, current_x_index+1) + &
                        var%data(k, current_y_index, current_x_index+2))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_L = w1*varL_1 + w2*varL_2 + w3*varL_3
    flux_L = var_L * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
              (0.5+sign(sign_term,mean_var))

    ! Negative fluxes at i+1/2 (take the symmetry of positive fluxes)

    !                             i+1/2
    !       i-2 ----- i-1 ----- i --|-- i+1 ----- i+2 ----- i+3 ----> x direction (u wind)
    !                                R
    !                                    |-------------------|
    !                                              S1
    !                           |------------------|
    !                                    S2
    !                  |-----------------|
    !                           S3

    ! Stencils S1 = [i+1, i+2, i+3], S2 = [i, i+1, i+2], S3 = [i-1, i, i+1]

    varR_1 = (11.0/6.0)*var%data(k, current_y_index, current_x_index+1) - &
             (7.0/6.0)*var%data(k, current_y_index, current_x_index+2) + &
             (1.0/3.0)*var%data(k, current_y_index, current_x_index+3)

    varR_2 = (1.0/3.0)*var%data(k, current_y_index, current_x_index) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index+1) - &
             (1.0/6.0)*var%data(k, current_y_index, current_x_index+2)

    varR_3 = - (1.0/6.0)*var%data(k, current_y_index, current_x_index-1) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index) + &
             (1.0/3.0)*var%data(k, current_y_index, current_x_index+1)

    beta1 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index+1) - &
                        2*var%data(k, current_y_index, current_x_index+2) + &
                        var%data(k, current_y_index, current_x_index+3))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index, current_x_index+1) - &
                        4*var%data(k, current_y_index, current_x_index+2) + &
                        var%data(k, current_y_index, current_x_index+3))**2

    beta2 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index) - &
                        2*var%data(k, current_y_index, current_x_index+1) + &
                        var%data(k, current_y_index, current_x_index+2))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index) - &
                      var%data(k, current_y_index, current_x_index+2))**2

    beta3 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index-1) - &
                        2*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index, current_x_index+1))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index-1) - &
                        4*var%data(k, current_y_index, current_x_index) + &
                        3*var%data(k, current_y_index, current_x_index+1))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_R = w1*varR_1 + w2*varR_2 + w3*varR_3
    flux_R = var_R * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
             (0.5-sign(sign_term,mean_var))

    flux_right_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to zero

  end subroutine x_right_face_flux_calc

  !Computation at the mass point on grid u(i-1/2,j,k)
  subroutine x_left_face_flux_calc(var, mean_var, current_x_index, current_y_index, k, flux_left_face)
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_left_face

    real(DEFAULT_PRECISION) :: alpha1, alpha2, alpha3, sum_alpha
    real(DEFAULT_PRECISION) :: beta1, beta2, beta3
    real(DEFAULT_PRECISION) :: w1, w2, w3

    real(DEFAULT_PRECISION) :: varL_1, varL_2, varL_3, var_L ! at i-1/2 and L=+
    real(DEFAULT_PRECISION) :: varR_1, varR_2, varR_3, var_R ! at i-1/2 and R=-
    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at i-1/2, terms results of the applying of Lax-Friedrich flux splitting method

    ! Positive fluxes at i-1/2

    !                             i-1/2
    !     i-3 ----- i-2 ----- i-1 --|-- i ----- i+1 ----- i+2 ----> x direction (u wind)
    !                              L
    !      |-------------------|
    !                S1
    !                |------------------|
    !                          S2
    !                          |-----------------|
    !                                   S3

    ! Stencils S1 = [i-3, i-2, i-1], S2 = [i-2, i-1, i], S3 = [i-1, i, i+1]

    varL_1 = (1.0/3.0)*var%data(k, current_y_index, current_x_index-3) - &
             (7.0/6.0)*var%data(k, current_y_index, current_x_index-2) + &
             (11.0/6.0)*var%data(k, current_y_index, current_x_index-1)

    varL_2 = - (1.0/6.0)*var%data(k, current_y_index, current_x_index-2) + &
               (5.0/6.0)*var%data(k, current_y_index, current_x_index-1) + &
               (1.0/3.0)*var%data(k, current_y_index, current_x_index)

    varL_3 = (1.0/3.0)*var%data(k, current_y_index, current_x_index-1) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index) - &
             (1.0/6.0)*var%data(k, current_y_index, current_x_index+1)

    beta1 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index-3) - &
                        2*var%data(k, current_y_index, current_x_index-2) + &
                        var%data(k, current_y_index, current_x_index-1))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index-3) - &
                        4*var%data(k, current_y_index, current_x_index-2) + &
                        3*var%data(k, current_y_index, current_x_index-1))**2

    beta2 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index-2) - &
                        2*var%data(k, current_y_index, current_x_index-1) + &
                        var%data(k, current_y_index, current_x_index))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index-2) - &
                      var%data(k, current_y_index, current_x_index))**2

    beta3 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index-1) - &
                        2*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index, current_x_index+1))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index, current_x_index-1) - &
                        4*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index, current_x_index+1))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_L = w1*varL_1 + w2*varL_2 + w3*varL_3
    flux_L = var_L * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
              (0.5+sign(sign_term,mean_var))

    ! Negative fluxes at i-1/2 (take the symmetry of positive fluxes)

    !                     i-1/2
    !       i-2 ----- i-1 --|-- i ----- i+1 ----- i+2 ----> x direction (u wind)
    !                        R
    !                           |------------------|
    !                                    S1
    !                  |-----------------|
    !                           S2
    !        |------------------|
    !                  S3

    ! Stencils S1 = [i, i+1, i+2], S2 = [i-1, i, i+1], S3 = [i-2, i-1, i]

    varR_1 = (11.0/6.0)*var%data(k, current_y_index, current_x_index) - &
             (7.0/6.0)*var%data(k, current_y_index, current_x_index+1) + &
             (1.0/3.0)*var%data(k, current_y_index, current_x_index+2)

    varR_2 = (1.0/3.0)*var%data(k, current_y_index, current_x_index-1) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index) - &
             (1.0/6.0)*var%data(k, current_y_index, current_x_index+1)

    varR_3 = - (1.0/6.0)*var%data(k, current_y_index, current_x_index-2) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index-1) + &
             (1.0/3.0)*var%data(k, current_y_index, current_x_index)

    beta1 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index) - &
                        2*var%data(k, current_y_index, current_x_index+1) + &
                        var%data(k, current_y_index, current_x_index+2))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index, current_x_index) - &
                        4*var%data(k, current_y_index, current_x_index+1) + &
                        var%data(k, current_y_index, current_x_index+2))**2

    beta2 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index-1) - &
                        2*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index, current_x_index+1))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index-1) - &
                      var%data(k, current_y_index, current_x_index+1))**2

    beta3 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index-2) - &
                        2*var%data(k, current_y_index, current_x_index-1) + &
                        var%data(k, current_y_index, current_x_index))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index-2) - &
                        4*var%data(k, current_y_index, current_x_index-1) + &
                        3*var%data(k, current_y_index, current_x_index))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_R = w1*varR_1 + w2*varR_2 + w3*varR_3
    flux_R = var_R * mean_var * & ! & var%data(k, current_y_index, current_x_index) * &
             (0.5-sign(sign_term,mean_var))

    flux_left_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to zero

  end subroutine x_left_face_flux_calc

  !Computation at the mass point on grid v(i,j+1/2,k)
  subroutine y_right_face_flux_calc(var, mean_var, current_x_index, current_y_index, k, flux_right_face)
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_right_face ! at i+1/2

    real(DEFAULT_PRECISION) :: alpha1, alpha2, alpha3, sum_alpha
    real(DEFAULT_PRECISION) :: beta1, beta2, beta3
    real(DEFAULT_PRECISION) :: w1, w2, w3

    real(DEFAULT_PRECISION) :: varL_1, varL_2, varL_3, var_L ! at j+1/2 and L=+
    real(DEFAULT_PRECISION) :: varR_1, varR_2, varR_3, var_R ! at j+1/2 and R=-
    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at j11/2, terms results of the applying of Lax-Friedrich flux splitting method

    ! Positive fluxes at j+1/2

    !                             j+1/2
    !       j-2 ----- j-1 ----- j --|-- j+1 ----- j+2 ----> y direction (v wind)
    !                              L
    !        |------------------|
    !                   S1
    !                   |----------------|
    !                           S2
    !                           |------------------|
    !                                    S3

    ! Stencils S1 = [j-2, j-1, j], S2 = [j-1, j, j+1], S3 = [j, j+1, j+2]

    varL_1 = (1.0/3.0)*var%data(k, current_y_index-2, current_x_index) - &
             (7.0/6.0)*var%data(k, current_y_index-1, current_x_index) + &
             (11.0/6.0)*var%data(k, current_y_index, current_x_index)

    varL_2 = - (1.0/6.0)*var%data(k, current_y_index-1, current_x_index) + &
               (5.0/6.0)*var%data(k, current_y_index, current_x_index) + &
               (1.0/3.0)*var%data(k, current_y_index+1, current_x_index)

    varL_3 = (1.0/3.0)*var%data(k, current_y_index, current_x_index) + &
             (5.0/6.0)*var%data(k, current_y_index+1, current_x_index) - &
             (1.0/6.0)*var%data(k, current_y_index+2, current_x_index)

    beta1 = (13.0/12.0)*(var%data(k, current_y_index-2, current_x_index) - &
                        2*var%data(k, current_y_index-1, current_x_index) + &
                        var%data(k, current_y_index, current_x_index))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index-2, current_x_index) - &
                        4*var%data(k, current_y_index-1, current_x_index) + &
                        3*var%data(k, current_y_index, current_x_index))**2

    beta2 = (13.0/12.0)*(var%data(k, current_y_index-1, current_x_index) - &
                        2*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index+1, current_x_index))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index-1, current_x_index) - &
                      var%data(k, current_y_index+1, current_x_index))**2

    beta3 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index) - &
                        2*var%data(k, current_y_index+1, current_x_index) + &
                        var%data(k, current_y_index+2, current_x_index))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index, current_x_index) - &
                        4*var%data(k, current_y_index+1, current_x_index) + &
                        var%data(k, current_y_index+2, current_x_index))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_L = w1*varL_1 + w2*varL_2 + w3*varL_3
    flux_L = var_L * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
              (0.5+sign(sign_term,mean_var))

    ! Negative fluxes at j+1/2 (take the symmetry of positive fluxes)

    !                             j+1/2
    !       j-2 ----- j-1 ----- j --|-- j+1 ----- j+2 ----- j+3 ----> y direction (v wind)
    !                                R
    !                                    |-------------------|
    !                                              S1
    !                           |------------------|
    !                                    S2
    !                  |-----------------|
    !                           S3

    ! Stencils S1 = [j+1, j+2, j+3], S2 = [j, j+1, j+2], S3 = [j-1, j, j+1]

    varR_1 = (11.0/6.0)*var%data(k, current_y_index+1, current_x_index) - &
             (7.0/6.0)*var%data(k, current_y_index+2, current_x_index) + &
             (1.0/3.0)*var%data(k, current_y_index+3, current_x_index)

    varR_2 = (1.0/3.0)*var%data(k, current_y_index, current_x_index) + &
             (5.0/6.0)*var%data(k, current_y_index+1, current_x_index) - &
             (1.0/6.0)*var%data(k, current_y_index+2, current_x_index)

    varR_3 = - (1.0/6.0)*var%data(k, current_y_index-1, current_x_index) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index) + &
             (1.0/3.0)*var%data(k, current_y_index+1, current_x_index)

    beta1 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index+1) - &
                        2*var%data(k, current_y_index, current_x_index+2) + &
                        var%data(k, current_y_index, current_x_index+3))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index, current_x_index+1) - &
                        4*var%data(k, current_y_index, current_x_index+2) + &
                        var%data(k, current_y_index, current_x_index+3))**2

    beta2 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index) - &
                        2*var%data(k, current_y_index, current_x_index+1) + &
                        var%data(k, current_y_index, current_x_index+2))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index) - &
                      var%data(k, current_y_index, current_x_index+2))**2

    beta3 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index-1) - &
                        2*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index, current_x_index+1))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index-1) - &
                        4*var%data(k, current_y_index, current_x_index) + &
                        3*var%data(k, current_y_index, current_x_index+1))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_R = w1*varR_1 + w2*varR_2 + w3*varR_3
    flux_R = var_R * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
             (0.5-sign(sign_term,mean_var))

    flux_right_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to zero

  end subroutine y_right_face_flux_calc


  !Computation at the mass point on grid v(i,j-1/2,k)
  subroutine y_left_face_flux_calc(var, mean_var, current_x_index, current_y_index, k, flux_left_face)
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_left_face ! at j-1/2

    real(DEFAULT_PRECISION) :: alpha1, alpha2, alpha3, sum_alpha
    real(DEFAULT_PRECISION) :: beta1, beta2, beta3
    real(DEFAULT_PRECISION) :: w1, w2, w3

    real(DEFAULT_PRECISION) :: varL_1, varL_2, varL_3, var_L ! at j-1/2 and L=+
    real(DEFAULT_PRECISION) :: varR_1, varR_2, varR_3, var_R ! at j-1/2 and R=-
    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at j-1/2, terms results of the applying of Lax-Friedrich flux splitting method


    ! Positive fluxes at j-1/2

    !                             j-1/2
    !     j-3 ----- j-2 ----- j-1 --|-- j ----- j+1 ----- j+2 ----> y direction (v wind)
    !                              L
    !      |-------------------|
    !                S1
    !                |------------------|
    !                          S2
    !                          |-----------------|
    !                                   S3

    ! Stencils S1 = [j-3, j-2, j-1], S2 = [j-2, j-1, j], S3 = [j-1, j, j+1]

    varL_1 = (1.0/3.0)*var%data(k, current_y_index-3, current_x_index) - &
             (7.0/6.0)*var%data(k, current_y_index-2, current_x_index) + &
             (11.0/6.0)*var%data(k, current_y_index-1, current_x_index)

    varL_2 = - (1.0/6.0)*var%data(k, current_y_index-2, current_x_index) + &
               (5.0/6.0)*var%data(k, current_y_index-1, current_x_index) + &
               (1.0/3.0)*var%data(k, current_y_index, current_x_index)

    varL_3 = (1.0/3.0)*var%data(k, current_y_index-1, current_x_index) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index) - &
             (1.0/6.0)*var%data(k, current_y_index+1, current_x_index)

    beta1 = (13.0/12.0)*(var%data(k, current_y_index-3, current_x_index) - &
                        2*var%data(k, current_y_index-2, current_x_index) + &
                        var%data(k, current_y_index-1, current_x_index))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index-3, current_x_index) - &
                        4*var%data(k, current_y_index-2, current_x_index) + &
                        3*var%data(k, current_y_index-1, current_x_index))**2

    beta2 = (13.0/12.0)*(var%data(k, current_y_index-2, current_x_index) - &
                        2*var%data(k, current_y_index-1, current_x_index) + &
                        var%data(k, current_y_index, current_x_index))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index-2, current_x_index) - &
                      var%data(k, current_y_index, current_x_index))**2

    beta3 = (13.0/12.0)*(var%data(k, current_y_index-1, current_x_index) - &
                        2*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index+1, current_x_index))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index-1, current_x_index) - &
                        4*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index+1, current_x_index))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_L = w1*varL_1 + w2*varL_2 + w3*varL_3
    flux_L = var_L * mean_var * & !var%data(k, current_y_index, current_x_index) * &
              (0.5+sign(sign_term,mean_var))

    ! Negative fluxes at j-1/2 (take the symmetry of positive fluxes)

    !                     j-1/2
    !       j-2 ----- j-1 --|-- j ----- j+1 ----- j+2 ----> y direction (v wind)
    !                        R
    !                           |------------------|
    !                                    S1
    !                  |-----------------|
    !                           S2
    !        |------------------|
    !                  S3

    ! Stencils S1 = [j, j+1, j+2], S2 = [j-1, j, j+1], S3 = [j-2, j-1, j]

    varR_1 = (11.0/6.0)*var%data(k, current_y_index, current_x_index) - &
             (7.0/6.0)*var%data(k, current_y_index+1, current_x_index) + &
             (1.0/3.0)*var%data(k, current_y_index+2, current_x_index)

    varR_2 = (1.0/3.0)*var%data(k, current_y_index-1, current_x_index) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index) - &
             (1.0/6.0)*var%data(k, current_y_index+1, current_x_index)

    varR_3 = - (1.0/6.0)*var%data(k, current_y_index-2, current_x_index) + &
             (5.0/6.0)*var%data(k, current_y_index-1, current_x_index) + &
             (1.0/3.0)*var%data(k, current_y_index, current_x_index)

    beta1 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index) - &
                        2*var%data(k, current_y_index+1, current_x_index) + &
                        var%data(k, current_y_index+2, current_x_index))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index, current_x_index) - &
                        4*var%data(k, current_y_index+1, current_x_index) + &
                        var%data(k, current_y_index+2, current_x_index))**2

    beta2 = (13.0/12.0)*(var%data(k, current_y_index-1, current_x_index) - &
                        2*var%data(k, current_y_index, current_x_index) + &
                        var%data(k, current_y_index+1, current_x_index))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index-1, current_x_index) - &
                      var%data(k, current_y_index+1, current_x_index))**2

    beta3 = (13.0/12.0)*(var%data(k, current_y_index-2, current_x_index) - &
                        2*var%data(k, current_y_index-1, current_x_index) + &
                        var%data(k, current_y_index, current_x_index))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index-2, current_x_index) - &
                        4*var%data(k, current_y_index-1, current_x_index) + &
                        3*var%data(k, current_y_index, current_x_index))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_R = w1*varR_1 + w2*varR_2 + w3*varR_3
    flux_R = var_R * mean_var * & !var%data(k, current_y_index, current_x_index) * &
             (0.5-sign(sign_term,mean_var))

    flux_left_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to zero

  end subroutine y_left_face_flux_calc

  subroutine z_top_face_flux_calc_weno1(current_state, var, mean_var, current_x_index, current_y_index, k, flux_top_face, top)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_top_face
    logical, intent(in) :: top

    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at k+1/2, terms results of the applying of Lax-Friedrich flux splitting method

    ! Positive fluxes at k+1/2

    flux_L = var%data(k, current_y_index, current_x_index) * current_state%global_grid%configuration%vertical%rho(k) * &
             mean_var * (0.5+sign(sign_term,mean_var))

    ! Negative fluxes at k+1/2 (take the symmetry of positive fluxes)

    if (top .eqv. .true.) then
      flux_R = var%data(k, current_y_index, current_x_index) * current_state%global_grid%configuration%vertical%rho(k) * &
               mean_var * (0.5-sign(sign_term,mean_var))
    else
      flux_R = var%data(k+1, current_y_index, current_x_index) * current_state%global_grid%configuration%vertical%rho(k+1) * &
               mean_var * (0.5-sign(sign_term,mean_var))
    end if

    flux_top_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to zero

  end subroutine z_top_face_flux_calc_weno1

  subroutine z_top_face_flux_calc_weno3(current_state, var, mean_var, current_x_index, current_y_index, k, flux_top_face)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_top_face

    real(DEFAULT_PRECISION) :: alpha1, alpha2, sum_alpha
    real(DEFAULT_PRECISION) :: beta1, beta2
    real(DEFAULT_PRECISION) :: w1, w2

    real(DEFAULT_PRECISION) :: varL_1, varL_2, var_L ! at k+1/2 and L=+
    real(DEFAULT_PRECISION) :: varR_1, varR_2, var_R ! at k+1/2 and R=-
    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at k+1/2, terms results of the applying of Lax-Friedrich flux splitting method

    ! Positive fluxes at k+1/2

    !                   k+1/2
    !       k-1 ----- k --|-- k+1  ----- k+2 ----> z direction (w wind)
    !                    L
    !        |--------|
    !            S1
    !                 |--------|
    !                     S2

    varL_1 = (3.0/2.0)*var%data(k, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k) - &
             (1.0/2.0)*var%data(k-1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k-1)
    varL_2 = (1.0/2.0)*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k) + &
             (1.0/2.0)*var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1)

    beta1 = (var%data(k, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k) - &
             var%data(k-1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k-1))**2
    beta2 = (var%data(k, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k) - &
             var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1))**2

    alpha1 = gamma1_weno3/(beta1 + epsilon_weno)**2
    alpha2 = gamma2_weno3/(beta2 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)

    var_L = w1*varL_1 + w2*varL_2
    flux_L = var_L * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
              (0.5+sign(sign_term,mean_var))

    ! Negative fluxes at k+1/2 (take the symmetry of positive fluxes)

    !                    k+1/2
    !       k-1 ----- k --|-- k+1 ----- k+2 ----> z direction (w wind)
    !                      R
    !                 |--------|
    !                     S1
    !                          |---------|
    !                               S2

    varR_1 = (3.0/2.0)*var%data(k+1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+1) - &
             (1.0/2.0)*var%data(k+2, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+2)

    varR_2 = (1.0/2.0)*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k) + &
             (1.0/2.0)*var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1)

    beta1 = (var%data(k+1, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k+1) - &
             var%data(k+2, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k+2))**2
    beta2 = (var%data(k, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k) - &
             var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1))**2

    alpha1 = gamma1_weno3/(beta1 + epsilon_weno)**2
    alpha2 = gamma2_weno3/(beta2 + epsilon_weno)**2

    var_R = w1*varR_1 + w2*varR_2
    flux_R = var_R * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
             (0.5-sign(sign_term,mean_var))

    flux_top_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to zero

  end subroutine z_top_face_flux_calc_weno3


  !Computation at the mass point on grid w(i,j,k+1/2)
  subroutine z_top_face_flux_calc_weno5(current_state, var, mean_var, current_x_index, current_y_index, k, flux_top_face)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_top_face

    real(DEFAULT_PRECISION) :: alpha1, alpha2, alpha3, sum_alpha
    real(DEFAULT_PRECISION) :: beta1, beta2, beta3
    real(DEFAULT_PRECISION) :: w1, w2, w3

    real(DEFAULT_PRECISION) :: varL_1, varL_2, varL_3, var_L ! at k+1/2 and L=+
    real(DEFAULT_PRECISION) :: varR_1, varR_2, varR_3, var_R ! at k+1/2 and R=-
    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at k+1/2, terms results of the applying of Lax-Friedrich flux splitting method

    ! Positive fluxes at k+1/2

    !                             k+1/2
    !       k-2 ----- k-1 ----- k --|-- k+1 ----- k+2 ----> z direction (w wind)
    !                              L
    !        |------------------|
    !                   S1
    !                   |----------------|
    !                           S2
    !                           |------------------|
    !                                    S3

    ! Stencils S1 = [k-2, k-1, k], S2 = [k-1, k, k+1], S3 = [k, k+1, k+2]

    varL_1 = (1.0/3.0)*var%data(k-2, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k-2) - &
             (7.0/6.0)*var%data(k-1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k-1) + &
             (11.0/6.0)*var%data(k, current_y_index, current_x_index)* &
             current_state%global_grid%configuration%vertical%rho(k)

    varL_2 = - (1.0/6.0)*var%data(k-1, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k-1) + &
               (5.0/6.0)*var%data(k, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k) + &
               (1.0/3.0)*var%data(k+1, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k+1)

    varL_3 = (1.0/3.0)*var%data(k, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k) + &
             (5.0/6.0)*var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1) - &
             (1.0/6.0)*var%data(k+2, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+2)

    beta1 = (13.0/12.0)*(var%data(k-2, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k-2) - &
                        2*var%data(k-1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-1) + &
                        var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k))**2 + &
            (1.0/4.0)*(var%data(k-2, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k-2) - &
                       4*var%data(k-1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-1)+ &
                        3*var%data(k, current_y_index, current_x_index)* &
                current_state%global_grid%configuration%vertical%rho(k))**2

    beta2 = (13.0/12.0)*(var%data(k-1, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k-1)- &
                        2*var%data(k, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k) + &
                        var%data(k+1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+1))**2 + &
            (1.0/4.0)*(var%data(k-1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-1)- &
                      var%data(k+1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+1))**2

    beta3 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k) - &
                        2*var%data(k+1, current_y_index, current_x_index) * &
                  current_state%global_grid%configuration%vertical%rho(k+1)+ &
                        var%data(k+2, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+2))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index, current_x_index) * &
                  current_state%global_grid%configuration%vertical%rho(k) - &
                        4*var%data(k+1, current_y_index, current_x_index) * &
                  current_state%global_grid%configuration%vertical%rho(k+1)+ &
                        var%data(k+2, current_y_index, current_x_index) * &
                  current_state%global_grid%configuration%vertical%rho(k+2))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_L = w1*varL_1 + w2*varL_2 + w3*varL_3
    !flux_L = var_L * var%data(k, current_y_index, current_x_index) * &
    flux_L = var_L * mean_var * &
              (0.5+sign(sign_term,mean_var))

    ! Negative fluxes at k+1/2 (take the symmetry of positive fluxes)

    !                             k+1/2
    !       k-2 ----- k-1 ----- k --|-- k+1 ----- k+2 ----- k+3 ----> z direction (w wind)
    !                                R
    !                                    |-------------------|
    !                                              S1
    !                           |------------------|
    !                                    S2
    !                  |-----------------|
    !                           S3

    ! Stencils S1 = [k+1, k+2, k+3], S2 = [k, k+1, k+2], S3 = [k-1, k, k+1]

    varR_1 = (11.0/6.0)*var%data(k+1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+1) - &
             (7.0/6.0)*var%data(k+2, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+2) + &
             (1.0/3.0)*var%data(k+3, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+3)

    varR_2 = (1.0/3.0)*var%data(k, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k) + &
             (5.0/6.0)*var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1) - &
             (1.0/6.0)*var%data(k+2, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+2)

    varR_3 = - (1.0/6.0)*var%data(k-1, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k-1) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k) + &
             (1.0/3.0)*var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1)

    beta1 = (13.0/12.0)*(var%data(k+1, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k+1) - &
                        2*var%data(k+2, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k+2) + &
                        var%data(k+3, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+3))**2 + &
            (1.0/4.0)*(3*var%data(k+1, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k+1) - &
                        4*var%data(k+2, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+2)+ &
                        var%data(k+3, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+3))**2

    beta2 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k) - &
                        2*var%data(k+1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+1) + &
                        var%data(k+2, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+2))**2 + &
            (1.0/4.0)*(var%data(k, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k) - &
                      var%data(k+2, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+2))**2

    beta3 = (13.0/12.0)*(var%data(k-1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-1) - &
                        2*var%data(k, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k) + &
                        var%data(k+1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+1))**2 + &
            (1.0/4.0)*(var%data(k-1, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k-1) - &
                        4*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k) + &
                        3*var%data(k+1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+1))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_R = w1*varR_1 + w2*varR_2 + w3*varR_3
    !flux_R = var_R * var%data(k, current_y_index, current_x_index) * &
    flux_R = var_R * mean_var * &
             (0.5-sign(sign_term,mean_var))

    flux_top_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to zero

  end subroutine z_top_face_flux_calc_weno5

  subroutine z_bottom_face_flux_calc_weno1(current_state, var, mean_var, current_x_index, current_y_index, k, flux_bottom_face)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_bottom_face

    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at k-1/2, terms results of the applying of Lax-Friedrich flux splitting method

    if (k .eq. 1) then
      flux_L = var%data(k, current_y_index, current_x_index) * current_state%global_grid%configuration%vertical%rho(k) * &
               mean_var * (0.5-sign(sign_term,mean_var))
    else
      flux_L = var%data(k-1, current_y_index, current_x_index) * current_state%global_grid%configuration%vertical%rho(k) * &
               mean_var * (0.5+sign(sign_term,mean_var))
    end if

    ! Negative fluxes at k+1/2 (take the symmetry of positive fluxes)

    flux_R = var%data(k, current_y_index, current_x_index) * current_state%global_grid%configuration%vertical%rho(k) * &
              mean_var * (0.5-sign(sign_term,mean_var))

    flux_bottom_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to ze

  end subroutine z_bottom_face_flux_calc_weno1

  subroutine z_bottom_face_flux_calc_weno3(current_state, var, mean_var, current_x_index, current_y_index, k, flux_bottom_face)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_bottom_face

    real(DEFAULT_PRECISION) :: alpha1, alpha2, sum_alpha
    real(DEFAULT_PRECISION) :: beta1, beta2
    real(DEFAULT_PRECISION) :: w1, w2

    real(DEFAULT_PRECISION) :: varL_1, varL_2, var_L ! at k-1/2 and L=+
    real(DEFAULT_PRECISION) :: varR_1, varR_2, var_R ! at k-1/2 and R=-
    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at k-1/2, terms results of the applying of Lax-Friedrich flux splitting method

    ! Positive fluxes at k-1/2

    !                   k-1/2
    !     k-2 ----- k-1 --|-- k ----- k+1 ----> z direction (w wind)
    !                    L
    !
    !      |---------|
    !           S1
    !                |--------|
    !                    S2

    varL_1 = (3.0/2.0)*var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1) - &
             (1.0/2.0)*var%data(k-2, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k-2)
    varL_2 = (1.0/2.0)*var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1) + &
             (1.0/2.0)*var%data(k, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k)

    beta1 = (var%data(k-1, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k-1) - &
             var%data(k-2, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k-2))**2
    beta2 = (var%data(k-1, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k-1) - &
             var%data(k, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k))**2

    alpha1 = gamma1_weno3/(beta1 + epsilon_weno)**2
    alpha2 = gamma2_weno3/(beta2 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)

    var_L = w1*varL_1 + w2*varL_2
    flux_L = var_L * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
              (0.5+sign(sign_term,mean_var))

    ! Negative fluxes at k+1/2 (take the symmetry of positive fluxes)

    !                   k-1/2
    !     k-2 ----- k-1 --|-- k ----- k+1 ----> z direction (w wind)
    !                      R
    !
    !                |--------|
    !                    S1
    !                         |--------|
    !                              S2

    varR_1 = (3.0/2.0)*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k) - &
             (1.0/2.0)*var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1)

    varR_2 = (1.0/2.0)*var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1) + &
             (1.0/2.0)*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k)

    beta1 = (var%data(k, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k) - &
             var%data(k+1, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k+1))**2
    beta2 = (var%data(k-1, current_y_index, current_x_index) * &
            current_state%global_grid%configuration%vertical%rho(k-1) - &
             var%data(k, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k))**2

    alpha1 = gamma1_weno3/(beta1 + epsilon_weno)**2
    alpha2 = gamma2_weno3/(beta2 + epsilon_weno)**2

    var_R = w1*varR_1 + w2*varR_2
    flux_R = var_R * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
             (0.5-sign(sign_term,mean_var))

    flux_bottom_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to zero

  end subroutine z_bottom_face_flux_calc_weno3



  !Computation at the mass point on grid w(i,j,k-1/2)
  subroutine z_bottom_face_flux_calc_weno5(current_state, var, mean_var, current_x_index, current_y_index, k, flux_bottom_face)
    type(model_state_type), target, intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: var
    real(DEFAULT_PRECISION), intent(in) :: mean_var
    integer, intent(in) :: current_x_index, current_y_index, k
    real(DEFAULT_PRECISION), intent(out) :: flux_bottom_face

    real(DEFAULT_PRECISION) :: alpha1, alpha2, alpha3, sum_alpha
    real(DEFAULT_PRECISION) :: beta1, beta2, beta3
    real(DEFAULT_PRECISION) :: w1, w2, w3

    real(DEFAULT_PRECISION) :: varL_1, varL_2, varL_3, var_L ! at k-1/2 and L=+
    real(DEFAULT_PRECISION) :: varR_1, varR_2, varR_3, var_R ! at k-1/2 and R=-
    real(DEFAULT_PRECISION) :: flux_L, flux_R ! at k-1/2, terms results of the applying of Lax-Friedrich flux splitting method

    ! Positive fluxes at k-1/2

    !                             k-1/2
    !     k-3 ----- k-2 ----- k-1 --|-- k ----- k+1 ----- k+2 ----> z direction (w wind)
    !                              L
    !      |-------------------|
    !                S1
    !                |------------------|
    !                          S2
    !                          |-----------------|
    !                                   S3

    ! Stencils S1 = [k-3, k-2, k-1], S2 = [k-2, k-1, k], S3 = [k-1, k, k+1]

    varL_1 = (1.0/3.0)*var%data(k-3, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k-3) - &
             (7.0/6.0)*var%data(k-2, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k-2) + &
             (11.0/6.0)*var%data(k-1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k-1)

    varL_2 = - (1.0/6.0)*var%data(k-2, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-2) + &
               (5.0/6.0)*var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1) + &
               (1.0/3.0)*var%data(k, current_y_index, current_x_index) * &
               current_state%global_grid%configuration%vertical%rho(k)

    varL_3 = (1.0/3.0)*var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1) + &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k) - &
             (1.0/6.0)*var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1)

    beta1 = (13.0/12.0)*(var%data(k-3, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-3)- &
                        2*var%data(k-2, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-2) + &
                        var%data(k-1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-1))**2 + &
            (1.0/4.0)*(var%data(k-3, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-3) - &
                        4*var%data(k-2, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-2) + &
                        3*var%data(k-1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-1))**2

    beta2 = (13.0/12.0)*(var%data(k-2, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-2) - &
                        2*var%data(k-1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-1) + &
                        var%data(k, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k))**2 + &
            (1.0/4.0)*(var%data(k-2, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-2) - &
                      var%data(k, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k))**2

    beta3 = (13.0/12.0)*(var%data(k-1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-1) - &
                        2*var%data(k, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k) + &
                        var%data(k+1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+1))**2 + &
            (1.0/4.0)*(3*var%data(k-1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k-1) - &
                        4*var%data(k, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k) + &
                        var%data(k+1, current_y_index, current_x_index) * &
                current_state%global_grid%configuration%vertical%rho(k+1))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_L = w1*varL_1 + w2*varL_2 + w3*varL_3
    flux_L = var_L * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
              (0.5+sign(sign_term,mean_var))

    ! Negative fluxes at k-1/2 (take the symmetry of positive fluxes)

    !                     k-1/2
    !       k-2 ----- k-1 --|-- k ----- k+1 ----- k+2 ----> z direction (w wind)
    !                        R
    !                           |------------------|
    !                                    S1
    !                  |-----------------|
    !                           S2
    !        |------------------|
    !                  S3

    ! Stencils S1 = [k, k+1, k+2], S2 = [k-1, k, k+1], S3 = [k-2, k-1, i]

    varR_1 = (11.0/6.0)*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k) - &
             (7.0/6.0)*var%data(k+1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+1) + &
             (1.0/3.0)*var%data(k+2, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+2)

    varR_2 = (1.0/3.0)*var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1)+ &
             (5.0/6.0)*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k) - &
             (1.0/6.0)*var%data(k+1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k+1)

    varR_3 = - (1.0/6.0)*var%data(k-2, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-2) + &
             (5.0/6.0)*var%data(k-1, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k-1) + &
             (1.0/3.0)*var%data(k, current_y_index, current_x_index) * &
             current_state%global_grid%configuration%vertical%rho(k)

    beta1 = (13.0/12.0)*(var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k) - &
                        2*var%data(k+1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+1) + &
                        var%data(k+2, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+2))**2 + &
            (1.0/4.0)*(3*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k) - &
                        4*var%data(k+1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+1)+ &
                        var%data(k+2, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+2))**2

    beta2 = (13.0/12.0)*(var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1) - &
                        2*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k) + &
                        var%data(k+1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+1))**2 + &
            (1.0/4.0)*(var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1)- &
                      var%data(k+1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k+1))**2

    beta3 = (13.0/12.0)*(var%data(k-2, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-2) - &
                        2*var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1) + &
                        var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k))**2 + &
            (1.0/4.0)*(var%data(k-2, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-2)- &
                        4*var%data(k-1, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k-1)+ &
                        3*var%data(k, current_y_index, current_x_index) * &
              current_state%global_grid%configuration%vertical%rho(k))**2

    alpha1 = gamma1/(beta1 + epsilon_weno)**2
    alpha2 = gamma2/(beta2 + epsilon_weno)**2
    alpha3 = gamma3/(beta3 + epsilon_weno)**2

    sum_alpha = alpha1 + alpha2 + alpha3

    w1 = (alpha1/sum_alpha)
    w2 = (alpha2/sum_alpha)
    w3 = (alpha3/sum_alpha)

    var_R = w1*varR_1 + w2*varR_2 + w3*varR_3
    flux_R = var_R * mean_var * & ! var%data(k, current_y_index, current_x_index) * &
             (0.5-sign(sign_term,mean_var))

    flux_bottom_face = 0.5*(flux_L + flux_R) ! depending on the velocity sign one will of the term will be equal to zero

  end subroutine z_bottom_face_flux_calc_weno5



  !> Parses a field string (read in from the configuration file) and determines whether this algorithm should be used
  !! for advecting that field
  !! @param field The string configuration of field advection
  !! @returns Whether or not the field is advected here
  logical function determine_if_advection_here(field)
    character(len=*), intent(in) :: field

    if (len_trim(field) .ne. 0) then
      if (trim(field) .eq. "weno5" .or. trim(field) .eq. "any") then
        determine_if_advection_here=.true.
      else
        determine_if_advection_here=.false.
      end if
    else
      determine_if_advection_here=.true.
    end if
  end function determine_if_advection_here

end module weno5_mod

