










!> The model state which represents the current state of a run
module state_mod
  use collections_mod, only : hashmap_type
  use grids_mod, only : global_grid_type, local_grid_type, immersed_boundary_type
  use prognostics_mod, only : prognostic_field_type
  use communication_types_mod, only : halo_communication_type
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  ! casim modules...
  use variable_precision, ONLY: wp
  implicit none

  private

  !> Stepping parameter values which determine centred or forward stepping
  integer, parameter, public :: CENTRED_STEPPING=0, FORWARD_STEPPING=1, PRESCRIBED_SURFACE_FLUXES=0, PRESCRIBED_SURFACE_VALUES=1
  !> The constants defining the reason why the model has terminated
  integer, parameter, public :: TIME_TERMINATION_REASON=0, TIMESTEP_TERMINATION_REASON=1, MESSAGE_TERMINATION_REASON=2, &
       WALLTIME_TERMINATION_REASON=3

  !> Information about the parallel aspects of the system
  type, public :: parallel_state_type
     integer :: processes, & !> Total number of processes
          my_rank,&           !> My process rank in the MONC system
          my_global_rank,&    !> My process rank in the global system
          neighbour_comm,&     !> Neighbour communicator
          monc_communicator=-1, io_communicator=-1, corresponding_io_server_process
     integer, dimension(3) :: &
          my_coords,&         !> My process coordinates in each dimension
          dim_sizes           !> Number of processes in each dimension
     logical, dimension(3,2) :: wrapped_around
     procedure(), nopass, pointer :: decomposition_procedure => null() !> The decomposition procedure to use
  end type parallel_state_type


  !> Information about the non-zero sampling intervals
  !!   Also used to track radiation timings when time_basis=.true.
  type, public :: sampling_interval_type
    integer :: interval = 0        ! sampling interval [ts, s if time_basis]
    integer, dimension(:), allocatable :: output   ! output intervals associated with %interval [s]
                                                   ! nint(output_frequency)
    integer :: next_time = 0       ! if time_basis, the next sample time for this %interval [s]
                                   !   if force_output_on_interval, the next output time for 
                                   !   this %interval [s]
    integer :: next_step = 0       ! the next sample timestep for this %interval [ts]
    logical :: active = .false.    ! .true.: sampling for this %interval on the current timestep
    logical :: radiation = .false. ! .true.: this %interval is used to track radiation calculations
  end type sampling_interval_type

  !> The ModelState which represents the current state of a run
  !!
  !! This state is provided to each callback and may be used and modified as required by
  !! the callbacks. Apart from this state, there should be no other state (global) variables
  !! declared. This allows us to simply persist and retrieve the ModelState when suspending
  !! and reactivating MONC.
  type, public :: model_state_type
    logical :: continue_timestep=.true., initialised=.false., continuation_run=.false.
    
    logical :: reconfig_run=.false.    ! whether this is the first cycle of a reconfigured 
                                       ! continuation run
    logical :: retain_model_time=.false.    ! by default, reconfigurations have model time 
                                            ! reset to zero 
    logical :: only_compute_on_sample_timestep=.false.    ! by default, diagnostics are available
                                  ! on every timestep.  When .true., certain diagnostics are only
                                  ! computed on specified diagnostic_sample_timesteps
    logical :: diagnostic_sample_timestep=.false.    ! diagnostics should be computed on the 
                                                     ! current timestep
    logical :: normal_step=.true.    ! the current timestep is a typical timestep, not a special, 
                                     ! shortened timestep due to proximity to a
                                     ! diagnostic_sample_timestep
    logical :: force_output_on_interval=.false.   ! allows the model to adjust the dtm to 
                                                  ! ensure that samples are sent to the IO
                                                  ! server on the output_frequency
                                                  ! time_basis=.true. does this automatically
    logical :: radiation_timestep=.false.  ! The current timestep is used for radiation
                                           ! calculations (determination is time_basis-sensitive)
    logical :: print_debug_data=.false.    ! Prints data for specific variables/points for
                                           ! debugging.  See registry.F90:execute_callbacks
    logical :: use_viscosity_and_diffusion=.true., &
       use_surface_boundary_conditions=.false., backscatter=.true.

    type(hashmap_type) :: options_database
    real(kind=DEFAULT_PRECISION), dimension(1500,1500) :: options_database_real
    character(len=STRING_LENGTH), dimension(1500,2) :: options_database_string
    character(len=200) :: checkpoint_path, diagnostic_path
    logical, dimension(100) :: options_database_logical
    type(global_grid_type) :: global_grid
    type(local_grid_type) :: local_grid
    type(parallel_state_type) :: parallel
    type(immersed_boundary_type) :: immersed
    type(prognostic_field_type) :: u, w, v, th, p, zu, zw, zv, zth, su, sw, sv, sth, savu, savv, & 
         savw, vis_coefficient, &
         diff_coefficient, dis, dis_th, &
         ! q species
         qv, zqv, sqv, dis_qv, &
         ql, zql, sql, dis_ql, &
         qr, zqr, sqr, dis_qr, &
         qi, zqi, sqi, dis_qi, &
         qs, zqs, sqs, dis_qs, &
         qg, zqg, sqg, dis_qg, &
         nl, znl, snl, dis_nl, &
         nr, znr, snr, dis_nr, &
         ni, zni, sni, dis_ni, &
         ns, zns, sns, dis_ns, &
         ng, zng, sng, dis_ng, &
         qAitkenSolMass, zqAitkenSolMass, sqAitkenSolMass, dis_qAitkenSolMass, &
         qAccumSolMass, zqAccumSolMass, sqAccumSolMass, dis_qAccumSolMass, &
         qAccumInsolMass, zqAccumInsolMass, sqAccumInsolMass, dis_qAccumInsolMass, &
         qCoarseSolMass, zqCoarseSolMass, sqCoarseSolMass, dis_qCoarseSolMass, &
         qCoarseDustMass, zqCoarseDustMass, sqCoarseDustMass, dis_qCoarseDustMass, &
         nAitkenSolNumber, znAitkenSolNumber, snAitkenSolNumber, dis_nAitkenSolNumber, &
         nAccumSolNumber, znAccumSolNumber, snAccumSolNumber, dis_nAccumSolNumber, &
         nAccumInsolNumber, znAccumInsolNumber, snAccumInsolNumber, dis_nAccumInsolNumber, &
         nCoarseSolNumber, znCoarseSolNumber, snCoarseSolNumber, dis_nCoarseSolNumber, &
         nCoarseDustnumber, znCoarseDustnumber, snCoarseDustnumber, dis_nCoarseDustnumber, &
         rdAitkenSol, rdAccumSol, rdCoarseSol, Tk, RH, RI, qv_saturation, qi_saturation, &
         D0_cloud, D0_rain, D0_ice, D0_snow, D0_graupel, &
         ! Heating rates from socrates contribution to sth
         sth_lw, sth_sw, cloud_reff, lwrad_hr, swrad_hr
    type(prognostic_field_type), dimension(:), allocatable :: q, zq, sq, disq
    type(prognostic_field_type), dimension(:), allocatable :: tracer, ztracer, stracer
    ! longwave and shortwave downwelling flux at the surface
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: sw_down_surf, lw_down_surf
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: cloud_reff_tot, longwave_hr_tot, shortwave_hr_tot

    ! advection contribution variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: u_advection, v_advection, &
       w_advection, th_advection, qv_advection, ql_advection, qr_advection, qi_advection, qs_advection, &
       qg_advection, qAitkenSolMass_advection, qAccumSolMass_advection, qAccumInsolMass_advection, &
       qCoarseSolMass_advection, qCoarseDustMass_advection, &
       nl_advection, nr_advection, ni_advection, ns_advection, ng_advection, &
       nAitkenSolNumber_advection, nAccumSolNumber_advection, nAccumInsolNumber_advection, &
       nCoarseSolNumber_advection, nCoarseDustnumber_advection

    ! viscosity contribution variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: u_viscosity, v_viscosity, w_viscosity

    ! diffusion contribution variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: th_diffusion, qv_diffusion, ql_diffusion, qr_diffusion, &
                               qi_diffusion, qs_diffusion, qg_diffusion, qAitkenSolMass_diffusion, qAccumSolMass_diffusion, &
                               qAccumInsolMass_diffusion, qCoarseSolMass_diffusion, qCoarseDustMass_diffusion, &
                               nl_diffusion, nr_diffusion, ni_diffusion, ns_diffusion, ng_diffusion, &
                               nAitkenSolNumber_diffusion, nAccumSolNumber_diffusion, nAccumInsolNumber_diffusion, &
                               nCoarseSolNumber_diffusion, nCoarseDustnumber_diffusion

    ! buoyancy contribution variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: w_buoyancy

    !!!!!!! diagnostics and checkpoint variables !!!!!!!!!!
    ! buoyancy
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: tend_3d_w_buoy
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tend_pr_tot_w_buoy
    ! casim.F90
    REAL(wp), allocatable :: surface_precip(:,:), surface_cloudsed(:,:), surface_rainsed(:,:)
    ! casim_profile_dgs
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       ! local process rate totals
       phomc_tot, & ! mass homogeneous freezing of cloud droplet rate
       pinuc_tot, & ! mass ice nucleation rate
       pidep_tot, & ! mass ice deposition rate
       psdep_tot, & ! mass snow deposition rate
       piacw_tot, & ! mass ice -> cloud -> ice accretion rate
       psacw_tot, & ! mass snow -> cloud -> snow accretion rate
       psacr_tot, & ! mass snow -> rain -> graupel and snow -> rain -> snow accretion rate
       pisub_tot, & ! mass ice sublimation rate
       pssub_tot, & ! mass snow sublimation rate
       pimlt_tot, & ! mass ice melting rate
       psmlt_tot, & ! mass snow melting rate
       psaut_tot, & ! mass autoconversion to snow rate
       psaci_tot, & ! mass snow -> ice -> snow accretion rate
       praut_tot, & ! mass autoconversion to rain rate
       pracw_tot, & ! mass rain accreting cloud rate
       prevp_tot, & ! mass evaporation of rain rate
       pgacw_tot, & ! mass graupel -> cloud -> graupel accretion rate
       pgacs_tot, & ! mass graupel -> snow accretion rate
       pgmlt_tot, & ! mass graupel melting rate
       pgsub_tot, & ! mass graupel sublimation rate
       psedi_tot, & ! mass ice sedimentation rate
       pseds_tot, & ! mass snow sedimentation rate
       psedr_tot, & ! mass rain sedimentation rate
       psedg_tot, & ! mass graupel sedimentation rate
       psedl_tot, & ! mass liquid/cloud sedimentation rate
       pcond_tot, & ! mass condensation rate
       ! local total tendencies
       dth_mphys_tot, dth_cond_evap_tot, dqv_mphys_tot, dqv_cond_evap_tot, &
       dqc_mphys_tot, dqr_mphys_tot, dqi_mphys_tot, dqs_mphys_tot, dqg_mphys_tot
    ! cond_diag_column
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: w_zn_cln, w_zn2_cln, tmp_th_cln, &
       wth_cln, th_pr_cln, wth_pr_cln, thv_pr_cln, wthv_pr_cln, th_pr2_cln, wthsg_cln, w_zn3_cln, &
       relhum_cln, tmp_u_cln, tmp_v_cln, wu_cln, wv_cln, wusg_cln, wvsg_cln, TdegK_cln, th_h_cln, &
       th_h_pr1_cln, th_h_pr2_cln,  qvli_cln, qvli_pr_cln, qvli_pr2_cln, qppt_cln, qppt_pr_cln, &
       qppt_pr2_cln, wqvli_pr_cln, wqppt_pr_cln
    ! coriolis
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: tend_3d_u_corio, tend_3d_v_corio
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tend_pr_tot_u_corio, tend_pr_tot_v_corio
    ! diagnostics_3d
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       TdegK,               & ! absolute temperature in kelvin
       theta,               & ! potential temperature in kelvin (th + thref)
       liquid_ice_theta       ! liquid-ice potential temperature in kelvin
    ! diffusion
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_th_diff,tend_3d_qv_diff,tend_3d_ql_diff,tend_3d_qi_diff,tend_3d_qr_diff,tend_3d_qs_diff, &
       tend_3d_qg_diff,tend_3d_tabs_diff
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
       tend_pr_tot_th_diff,tend_pr_tot_qv_diff,tend_pr_tot_ql_diff,tend_pr_tot_qi_diff,tend_pr_tot_qr_diff, &
       tend_pr_tot_qs_diff,tend_pr_tot_qg_diff, tend_pr_tot_tabs_diff
    ! forcing
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_u_forc, tend_3d_v_forc, tend_3d_th_forc,tend_3d_qv_forc, &
       tend_3d_ql_forc,tend_3d_qi_forc,tend_3d_qr_forc,tend_3d_qs_forc,tend_3d_qg_forc, &
       tend_3d_tabs_forc
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::                         &
       tend_pr_tot_u_forc, tend_pr_tot_v_forc, tend_pr_tot_th_forc,tend_pr_tot_qv_forc, &
       tend_pr_tot_ql_forc,tend_pr_tot_qi_forc,tend_pr_tot_qr_forc,tend_pr_tot_qs_forc,tend_pr_tot_qg_forc, &
       tend_pr_tot_tabs_forc
    ! profile_diag arrays used to store the change in field due to forcing
    real(kind=DEFAULT_PRECISION), allocatable :: du_profile_diag(:), dv_profile_diag(:), dtheta_profile_diag(:), &
       dqv_profile_diag(:), dql_profile_diag(:), dqr_profile_diag(:), dqi_profile_diag(:), dqs_profile_diag(:), &
       dqg_profile_diag(:), dqAitkenSolMass_profile_diag(:), dqAccumSolMass_profile_diag(:), dqAccumInsolMass_profile_diag(:), &
       dqCoarseSolMass_profile_diag(:), dqCoarseDustMass_profile_diag(:)
    ! subs_profile_diag arrays used to store the change in field due to subsidence
    real(kind=DEFAULT_PRECISION), allocatable :: du_subs_profile_diag(:), dv_subs_profile_diag(:), &
       dtheta_subs_profile_diag(:), dqv_subs_profile_diag(:), dql_subs_profile_diag(:), dqr_subs_profile_diag(:), &
       dqi_subs_profile_diag(:), dqs_subs_profile_diag(:), dqg_subs_profile_diag(:), dqAitkenSolMass_subs_profile_diag(:), &
       dqAccumSolMass_subs_profile_diag(:) , dqAccumInsolMass_subs_profile_diag(:), dqCoarseSolMass_subs_profile_diag(:), &
       dqCoarseDustMass_subs_profile_diag(:)
    ! profile_diag.F90
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       u_wind_tot, uprime_tot, v_wind_tot, vprime_tot,  &
       uprime, vprime, wke_tot, wwww_tot, www_tot, ww_tot,       &
       theta_tot, w_wind_tot, rh_tot, wtheta_ad_tot,             &
       wtheta_cn_tot, uw_tot, vw_tot, uv_tot, th2_tot,           &
       thinit, uinit, vinit,            &
       ! mositure means
       q_temp, qv_tot, ql_tot, qr_tot, qi_tot, qs_tot, qg_tot,   &
       ! number concentrations
       nl_tot, nr_tot, ni_tot, ns_tot, ng_tot,   &
       ! moisture flux terms
       wqv_cn_tot, wql_cn_tot, wqr_cn_tot, wqi_cn_tot,           &
       wqs_cn_tot, wqg_cn_tot,                                   &
       wqv_ad_tot, wql_ad_tot, wqr_ad_tot, wqi_ad_tot,           &
       wqs_ad_tot, wqg_ad_tot
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: cloud_mask
    real(kind=DEFAULT_PRECISION), dimension(:)    , allocatable :: cloud_mask_tot
    real(kind=DEFAULT_PRECISION), dimension(:)    , allocatable :: cloud_liq_mask_tot
    real(kind=DEFAULT_PRECISION), dimension(:)    , allocatable :: cloud_ice_mask_tot
    ! pstep diag terms
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
        tendp_3d_u_pt, tendp_3d_v_pt, tendp_3d_w_pt, tend_3d_u_pt, tend_3d_v_pt, tend_3d_w_pt
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
        tendp_pr_tot_u_pt, tendp_pr_tot_v_pt, tendp_pr_tot_w_pt, tend_pr_tot_u_pt, tend_pr_tot_v_pt, tend_pr_tot_w_pt
    ! pwad diag terms
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_u_pwad, tend_3d_v_pwad, tend_3d_w_pwad, tend_3d_th_pwad,tend_3d_qv_pwad,       &
       tend_3d_ql_pwad,tend_3d_qi_pwad,tend_3d_qr_pwad,tend_3d_qs_pwad,tend_3d_qg_pwad,       &
       tend_3d_tabs_pwad
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::                             &
       tend_pr_tot_u_pwad, tend_pr_tot_v_pwad, tend_pr_tot_w_pwad, tend_pr_tot_th_pwad,tend_pr_tot_qv_pwad,       &
       tend_pr_tot_ql_pwad,tend_pr_tot_qi_pwad,tend_pr_tot_qr_pwad,tend_pr_tot_qs_pwad,tend_pr_tot_qg_pwad,       &
       tend_pr_tot_tabs_pwad
    ! scalar diag terms
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: vwp, lwp, wmax, wmin, &
       qlmax, hqlmax, cltop, clbas,  senhf, lathf, rwp, iwp, swp, gwp, tot_iwp,      &
       reske, &
       vwp_usr, lwp_usr, rwp_usr, iwp_usr, swp_usr, gwp_usr, tot_iwp_usr
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: subke_2d
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: trsfflux
    ! stepfield diag terms
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_th_sf,tend_3d_qv_sf,       &
       tend_3d_ql_sf,tend_3d_qi_sf,tend_3d_qr_sf,tend_3d_qs_sf,tend_3d_qg_sf,       &
       tend_3d_tabs_sf
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
       tend_pr_tot_th_sf,tend_pr_tot_qv_sf,       &
       tend_pr_tot_ql_sf,tend_pr_tot_qi_sf,tend_pr_tot_qr_sf,tend_pr_tot_qs_sf,tend_pr_tot_qg_sf,       &
       tend_pr_tot_tabs_sf
    ! subgrid diag terms
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::     &
       uwsg_tot, vwsg_tot, uusg_tot, vvsg_tot, wwsg_tot,         &
       tkesg_tot, wtsg_tot, th2sg_tot, wqsg_tot,                 &
       ! subgrid tke fluxes
       sed_tot,ssub_tot, dissipation_tot,buoysg_tot, wkesg_tot,  &
       theta_dis_tot, vis_coef_tot, diff_coef_tot,               &
       richardson_number_tot, richardson_squared_tot,            &
       ! subgrid moisture fluxes
       wqv_sg_tot, wql_sg_tot, wqr_sg_tot, wqi_sg_tot,           &
       wqs_sg_tot, wqg_sg_tot
    ! thad diag terms
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: tend_3d_th_thad, tend_3d_tabs_thad
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tend_pr_tot_th_thad, tend_pr_tot_tabs_thad
    ! tvad diag terms
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_u_tvad, tend_3d_v_tvad, tend_3d_w_tvad, tend_3d_th_tvad,tend_3d_qv_tvad, &
       tend_3d_ql_tvad,tend_3d_qi_tvad,tend_3d_qr_tvad,tend_3d_qs_tvad,tend_3d_qg_tvad, &
       tend_3d_tabs_tvad
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::                             &
       tend_pr_tot_u_tvad, tend_pr_tot_v_tvad, tend_pr_tot_w_tvad,tend_pr_tot_th_tvad,tend_pr_tot_qv_tvad, &
       tend_pr_tot_ql_tvad,tend_pr_tot_qi_tvad,tend_pr_tot_qr_tvad,tend_pr_tot_qs_tvad,tend_pr_tot_qg_tvad, &
       tend_pr_tot_tabs_tvad
    ! visc diag terms
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_u_visc, tend_3d_v_visc, tend_3d_w_visc
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
       tend_pr_tot_u_visc, tend_pr_tot_v_visc, tend_pr_tot_w_visc
    ! SOCRATES diag terms
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
      tend_3d_tabs_lw, tend_3d_tabs_sw, tend_3d_tabs_total
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
      tend_pr_tot_th_lw,    tend_pr_tot_tabs_lw, tend_pr_tot_th_sw,    tend_pr_tot_tabs_sw,   &
      tend_pr_tot_th_total, tend_pr_tot_tabs_total
    ! flux budget terms
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: th_flux_values, th_gradient, th_diff, th_buoyancy, th_tendency, &
       uw_advection, vw_advection, uw_viscosity, vw_viscosity, uw_buoyancy, vw_buoyancy, uw_tendency, vw_tendency, uw_w, vw_w, &
       tu_su, uu_advection, uu_viscosity, wu_u, tv_sv, vv_advection, vv_viscosity, wv_v, tw_sw, ww_advection, &
       ww_viscosity, ww_buoyancy, &
       u_thetal, us_thetal, u_thetal_advection, u_thetal_viscosity_diffusion, wu_thetal, v_thetal, vs_thetal, &
       v_thetal_advection, v_thetal_viscosity_diffusion, wv_thetal, w_thetal, ws_thetal, &
       w_thetal_advection, w_thetal_viscosity_diffusion, w_thetal_buoyancy, ww_thetal, thetal_thetal, sthetal_thetal, &
       thetal_thetal_advection, thetal_thetal_diffusion, wthetal_thetal, &
       u_mse, us_mse, u_mse_advection, u_mse_viscosity_diffusion, wu_mse, v_mse, vs_mse, v_mse_advection, &
       v_mse_viscosity_diffusion, wv_mse, w_mse, ws_mse, w_mse_advection, w_mse_viscosity_diffusion, w_mse_buoyancy, &
       ww_mse, mse_mse, smse_mse, mse_mse_advection, mse_mse_diffusion, wmse_mse, &
       us_qt, u_qt_advection, u_qt_viscosity_diffusion, wu_qt, vs_qt, v_qt_advection, &
       v_qt_viscosity_diffusion, wv_qt, w_qt, ws_qt, w_qt_advection, w_qt_viscosity_diffusion, w_qt_buoyancy, &
       ww_qt, qt_qt, sqt_qt, qt_qt_advection, qt_qt_diffusion, wqt_qt, &
       shprd, w_ke, sub_buoy, w_p, tke_tend
    real(kind=DEFAULT_PRECISION) :: mflux, wmfcrit
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: qv_flux_values, qv_gradient, qv_diff, qv_buoyancy, qv_tendency, &
      ql_flux_values, ql_gradient, ql_diff, ql_buoyancy, ql_tendency, &
      qr_flux_values, qr_gradient, qr_diff, qr_buoyancy, qr_tendency, &
      qi_flux_values, qi_gradient, qi_diff, qi_buoyancy, qi_tendency, &
      qs_flux_values, qs_gradient, qs_diff, qs_buoyancy, qs_tendency, &
      qg_flux_values, qg_gradient, qg_diff, qg_buoyancy, qg_tendency, &
      qAitkenSolMass_flux_values, qAitkenSolMass_gradient, qAitkenSolMass_diff, qAitkenSolMass_buoyancy, qAitkenSolMass_tendency, &
      qAccumSolMass_flux_values, qAccumSolMass_gradient, qAccumSolMass_diff, qAccumSolMass_buoyancy, qAccumSolMass_tendency, &
      qAccumInsolMass_flux_values, qAccumInsolMass_gradient, qAccumInsolMass_diff, &
      qAccumInsolMass_buoyancy, qAccumInsolMass_tendency, &
      qCoarseSolMass_flux_values, qCoarseSolMass_gradient, qCoarseSolMass_diff, qCoarseSolMass_buoyancy, qCoarseSolMass_tendency, &
      qCoarseDustMass_flux_values, qCoarseDustMass_gradient, qCoarseDustMass_diff, &
      qCoarseDustMass_buoyancy, qCoarseDustMass_tendency, &
      nl_flux_values, nl_gradient, nl_diff, nl_buoyancy, nl_tendency, &
      nr_flux_values, nr_gradient, nr_diff, nr_buoyancy, nr_tendency, &
      ni_flux_values, ni_gradient, ni_diff, ni_buoyancy, ni_tendency, &
      ns_flux_values, ns_gradient, ns_diff, ns_buoyancy, ns_tendency, &
      ng_flux_values, ng_gradient, ng_diff, ng_buoyancy, ng_tendency, &
      nAitkenSolNumber_flux_values, nAitkenSolNumber_gradient, nAitkenSolNumber_diff, &
      nAitkenSolNumber_buoyancy, nAitkenSolNumber_tendency, &
      nAccumSolNumber_flux_values, nAccumSolNumber_gradient, nAccumSolNumber_diff, &
      nAccumSolNumber_buoyancy, nAccumSolNumber_tendency, &
      nAccumInsolNumber_flux_values, nAccumInsolNumber_gradient, nAccumInsolNumber_diff, &
      nAccumInsolNumber_buoyancy, nAccumInsolNumber_tendency, &
      nCoarseSolNumber_flux_values, nCoarseSolNumber_gradient, nCoarseSolNumber_diff, &
      nCoarseSolNumber_buoyancy, nCoarseSolNumber_tendency, &
      nCoarseDustnumber_flux_values, nCoarseDustnumber_gradient, nCoarseDustnumber_diff, &
      nCoarseDustnumber_buoyancy, nCoarseDustnumber_tendency
    ! pdf analysis terminated
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: w_histogram_profile_local
    ! averaged data
    real(kind=DEFAULT_PRECISION) :: mean_mo_l, mean_friction_vel, mean_z0, mean_z0th
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable  :: mean_th, mean_u, mean_v, &
                                  mean_qv, mean_ql, mean_qr, mean_qi, mean_qs, mean_qg, &
                                  mean_qAitkenSolMass, mean_qAccumSolMass, mean_qAccumInsolMass, &
                                  mean_qCoarseSolMass, mean_qCoarseDustMass, &
                                  mean_nl, mean_nr, mean_ni, mean_ns, mean_ng, &
                                  mean_nAitkenSolNumber, mean_nAccumSolNumber, mean_nAccumInsolNumber, &
                                  mean_nCoarseSolNumber, mean_nCoarseDustnumber
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! state of components (active of not)
    logical ::                                                         &
       buoyancy_enabled = .false., forcing_enabled = .false., profile_diagnostics_enabled = .false., &
       casim_profile_dgs_enabled = .false., scalar_diagnostics_enabled = .false., stepfields_enabled = .false., &
       mean_profiles_enabled = .false., pstep_enabled = .false., conditional_diagnostics_column_enabled = .false., &
       diagnostics_3d_enabled = .false., conditional_diagnostics_whole_enabled = .false., &
       subgrid_profile_diagnostics_enabled = .false., casim_enabled = .false., coriolis_enabled = .false., &
       diffusion_enabled = .false., pw_advection_enabled = .false., th_advection_enabled = .false., &
       tvd_advection_enabled = .false., viscosity_enabled = .false., socrates_enabled = .false., &
       simplecloud_enabled = .false., flux_budget_enabled= .false., pdf_analysis_enabled = .false., &
       diagnostics_averaged_enabled = .false.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type(halo_communication_type) :: viscosity_halo_swap_state, diffusion_halo_swap_state
    type(sampling_interval_type), dimension(:), allocatable :: sampling
    real(kind=DEFAULT_PRECISION) :: time=.0_DEFAULT_PRECISION,& ! Model time in seconds
            dtm,& ! Modeltimestep (s)
            absolute_new_dtm, &
            thref0,&
            rhobous,&
            tsmth=1e-2_DEFAULT_PRECISION,&
            timestep_runtime,&
            local_divmax, global_divmax, cvis=0.0_DEFAULT_PRECISION, surface_temperature_flux, &
            surface_vapour_flux, theta_surf, surface_vapour_mixing_ratio, fbuoy, &
            fbuoynew, theta_virtual_surf, cmbc, rcmbc, ellmocon, velmax, velmin, aloginv, cneut, cfc, mo_l = 0.0, &
            friction_vel = 0.0, &
            surface_pressure=100000.0_DEFAULT_PRECISION, surface_reference_pressure = 100000.0_DEFAULT_PRECISION, &
            cvel, cvel_x, cvel_y, cvel_z, dtm_new, rmlmax, geostrophic_wind_rate_of_change_in_x, &
            geostrophic_wind_rate_of_change_in_y, surface_geostrophic_wind_x, surface_geostrophic_wind_y, &
            local_zumin, local_zumax, local_zvmin, local_zvmax, local_cvel_z, termination_time = 0.0
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: lookup_table_velocity, &
         lookup_table_ustr, cq, abswmax
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: psrce_recv_buffer_x, psrce_recv_buffer_y
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tracer_decay_rate, tracer_surf_bc
    integer :: timestep=1, column_global_x, column_global_y, column_local_x, column_local_y, field_stepping, scalar_stepping, &
         momentum_stepping, number_q_fields=0, start_timestep=1, type_of_surface_boundary_conditions = 0, lookup_table_entries, &
         cfl_frequency, termination_reason, last_cfl_timestep=0, reconfig_timestep_offset=0, checkpoint_frequency = 0, &
         diagnostic_file_1d_write_frequency = 0, diagnostic_file_2d_write_frequency = 0, diagnostic_file_3d_write_frequency = 0, &
         diagnostic_file_0d_write_frequency = 0, diagnostic_file_averaged_sample_frequency = 0, &
         diagnostic_file_averaged_write_frequency = 0, &
         modulo_number_0d = -1, modulo_number_1d = -1, modulo_number_2d = -1, modulo_number_3d = -1, &
         modulo_number_check = -1, modulo_number_averaged_sample = -1, modulo_number_averaged_write = -1, &
         dim_var_to_mean = -1, &
         average_sample_iter = 1
    integer :: water_vapour_mixing_ratio_index=0, liquid_water_mixing_ratio_index=0, &
         rain_water_mixing_ratio_index=0, ice_water_mixing_ratio_index=0, &
         snow_water_mixing_ratio_index=0, graupel_water_mixing_ratio_index=0, & 
         ! integer index for number concentrations
         liquid_water_nc_index=0, &
         rain_water_nc_index=0, ice_water_nc_index=0, &
         snow_water_nc_index=0, graupel_water_nc_index=0, &
         psrce_x_hs_send_request, psrce_y_hs_send_request, psrce_x_hs_recv_request, psrce_y_hs_recv_request
    integer :: n_tracers=0, n_radioactive_tracers=0
    integer :: traj_tracer_index=0, radioactive_tracer_index=0
    integer, dimension(:), allocatable:: tracer_surf_bc_opt
    logical :: first_timestep_column, last_timestep_column, halo_column, first_nonhalo_timestep_column, &
         passive_q=.false., passive_th=.false., &
         use_time_varying_surface_values, use_anelastic_equations, & ! use_anelastic_equations or use Boussinesq
         saturated_surface, update_dtm=.false., calculate_th_and_q_init, origional_vertical_grid_setup, &
         new_vertical_grid_setup, special_vertical_grid_setup, &
         io_server_enabled, reinit_tracer=.false., time_basis=.false., time_frequency_enabled=.false., &
         diag_0d_done=.false., diag_1d_done=.false., diag_2d_done=.false., diag_3d_done=.false., checkpt_done=.false.
    logical, allocatable :: l_forceq(:)
    double precision :: model_start_wtime

    logical :: galilean_transformation=.true., fix_ugal=.false., fix_vgal=.false., init_to_apply_random_noise = .true.
    logical :: tracking_variables_enabled = .false., constant_w_profile_enabled = .false.
    real(kind=DEFAULT_PRECISION) :: ugal=0.,vgal=0., time_to_apply_random_noise=0.0, constant_w_profile_timer = 0.0, &
                                    duration_random_noise = 0.0
    ! SOCRATES time variables are included in state since they need to be dumped
    real(kind=DEFAULT_PRECISION) :: rad_last_time=0.0
    ! Global grid location for print_debug_data
    integer :: pdd_z=-999, pdd_y=-999, pdd_x=-999
    integer :: config_args
    integer :: n_w_bins

  end type model_state_type
end module state_mod
