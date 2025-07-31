!> Conditionally averaged diagnostics, Part 1 of 2.
!  On vertical levels, identifies diagnostic values under specified conditions and 
!  records the values and their count associated with that condition and its complement (.not. condition).
module conditional_diagnostics_column_mod
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type, &
       component_descriptor_type_v1
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use q_indices_mod, only: get_q_index, standard_q_names
  use saturation_mod, only: qsaturation
  use optionsdatabase_mod, only : options_has_key, options_get_logical, options_get_integer, &
                                  options_get_string, options_get_real, options_get_array_size
  use science_constants_mod, only : cp, rlvap
  use logging_mod, only : LOG_ERROR, log_master_log

  implicit none

#ifndef TEST_MODE
  private
#endif


  integer ::  iqv, iql, iqr, iqi, iqs, iqg                      

  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: CondDiags_tot

  real(kind=DEFAULT_PRECISION) :: qlcrit
  real(kind=DEFAULT_PRECISION) :: qicrit
  real(kind=DEFAULT_PRECISION) :: qpptcrit
  real(kind=DEFAULT_PRECISION) :: vpptcrit
  real(kind=DEFAULT_PRECISION) :: thvprcrit
  real(kind=DEFAULT_PRECISION) :: wSdwncrit
  real(kind=DEFAULT_PRECISION) :: wSupcrit
  real(kind=DEFAULT_PRECISION) :: wupcrit
  real(kind=DEFAULT_PRECISION) :: wdwncrit

  integer :: ncond 
  integer :: ndiag  
  integer :: x_size
  integer :: y_size

  real(kind=DEFAULT_PRECISION) :: gpts_total

  logical :: thv_from_th_with_liqice
  logical :: l_qi_qr_qs_qg
  logical :: l_do_near

  character(len=STRING_LENGTH), dimension(23) :: master_conditions_list
  character(len=STRING_LENGTH), dimension(31) :: master_diagnostics_list

  character(len=STRING_LENGTH), dimension(:), allocatable :: cond_request, cond_long
  character(len=STRING_LENGTH), dimension(:), allocatable :: diag_request, diag_long
  integer, dimension(:), allocatable :: diag_locations, cond_locations
  integer :: requested_area

  public CondDiags_tot, ncond, ndiag, gpts_total, requested_area, &
         cond_request, diag_request, cond_long, diag_long, initialisation_callback_conditional_diagnostics_column, &
         timestep_callback_diagnostics_column, finalisation_callback_conditional_diagnostics_column


contains

  subroutine initialisation_callback_conditional_diagnostics_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: total_points, inc, iter, intnum
    real :: realnum
    logical :: logicnum

    !----------------------------------------------------------------------
    !>> This is section contains the full list of available conditions and 
    !   diagnostics.  
    !
    !   If you add more, you will need to update some of these:
    !         1. master_conditions_list (and its declared size)
    !         2. master_diagnostics_list (and its declared size)
    !         3. cond_locations CASE checks
    !         4. diag_locations CASE checks
    !         5. local_diag value
    !         6. l_cond logical test
    !         7. add short name to global_config default list
    !
    !   TAKE CARE TO ENSURE that indices correspond directly between:
    !         1. master_conditions_list / cond_locations
    !         2. master_diagnostics_list / diag_locations / local_diag
    !
    !----------------------------------------------------------------------

    master_conditions_list(1) = "all points"
    master_conditions_list(2) = "buoyant updraft"
    master_conditions_list(3) = "buoyant cloudy updraft"
    master_conditions_list(4) = "near buoyant cloudy updraft"
    master_conditions_list(5) = "cloudy"
    master_conditions_list(6) = "cloudy updraft"
    master_conditions_list(7) = "cloudy downdraft"
    master_conditions_list(8) = "strong updraft"
    master_conditions_list(9) = "strong downdraft"
    master_conditions_list(10) = "updraft"
    master_conditions_list(11) = "downdraft"
    master_conditions_list(12) = "liquid cloudy updraft"
    master_conditions_list(13) = "liquid cloudy downdraft"
    master_conditions_list(14) = "hydrometeors"
    master_conditions_list(15) = "cloud liquid water"
    master_conditions_list(16) = "cloud ice"
    master_conditions_list(17) = "precipitating downdraft"
    master_conditions_list(18) = "large precipitating downdraft"
    master_conditions_list(19) = "precipitating strong downdraft"
    master_conditions_list(20) = "saturated updraft"
    master_conditions_list(21) = "buoyant saturated updraft"
    master_conditions_list(22) = "cloud-free precipitating downdraft"
    master_conditions_list(23) = "cloud-free large precipitating downdraft"

    master_diagnostics_list(1) = "Fractional area"
    master_diagnostics_list(2) = "Vertical velocity"
    master_diagnostics_list(3) = "Vertical velocity variance"
    master_diagnostics_list(4) = "Potential temperature"
    master_diagnostics_list(5) = "w * th"
    master_diagnostics_list(6) = "Potential temperature anomalies"
    master_diagnostics_list(7) = "Flux of potential temperature"
    master_diagnostics_list(8) = "Virtual potential temperature anomalies"
    master_diagnostics_list(9) = "Flux of virtual potential temperature"
    master_diagnostics_list(10) = "Variance of potential temperature"
    master_diagnostics_list(11) = "Diff_coeff * dth/dz"
    master_diagnostics_list(12) = "w**3"
    master_diagnostics_list(13) = "Relative humidity"
    master_diagnostics_list(14) = "Zonal velocity"
    master_diagnostics_list(15) = "Meridional velovity"
    master_diagnostics_list(16) = "w * u"
    master_diagnostics_list(17) = "w * v"
    master_diagnostics_list(18) = "Vis_coeff * du/dz"
    master_diagnostics_list(19) = "Vis_coeff * dv/dz"
    master_diagnostics_list(20) = "Temperature"
    master_diagnostics_list(21) = "th+lvap*qv/cp"
    master_diagnostics_list(22) = "Anomalies of THL"
    master_diagnostics_list(23) = "Variance of THL"
    master_diagnostics_list(24) = "qv+ql+qi"
    master_diagnostics_list(25) = "Anomalies of QVLI"
    master_diagnostics_list(26) = "Variance of QVLI"
    master_diagnostics_list(27) = "qr+qs+qg"
    master_diagnostics_list(28) = "Anomalies of QRSG"
    master_diagnostics_list(29) = "Variance of QRSG"
    master_diagnostics_list(30) = "Flux of QVLI"
    master_diagnostics_list(31) = "Flux of QRSG"


    total_points=(current_state%local_grid%size(Y_INDEX) * current_state%local_grid%size(X_INDEX))
 
    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then                        
      iqv       = 1  !get_q_index(standard_q_names%VAPOUR, 'profile_diags')
      iql       = 2 !get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'profile_diags')
      iqi       = 3 !get_q_index(standard_q_names%ICE_MASS, 'profile_diags')
      iqr       = 4 !get_q_index(standard_q_names%RAIN_MASS, 'profile_diags')
      iqs       = 5 !get_q_index(standard_q_names%SNOW_MASS, 'profile_diags')
      iqg       = 6 !get_q_index(standard_q_names%GRAUPEL_MASS, 'profile_diags')

      do iter = 1,current_state%config_args
        if (current_state%options_database_string(iter,1) .eq. "qlcrit") then
          read(current_state%options_database_string(iter,2),*) realnum
          qlcrit = realnum
        else if (current_state%options_database_string(iter,1) .eq. "qicrit") then
          read(current_state%options_database_string(iter,2),*) realnum
          qicrit = realnum
        else if (current_state%options_database_string(iter,1) .eq. "qpptcrit") then
          read(current_state%options_database_string(iter,2),*) realnum
          qpptcrit = realnum
        else if (current_state%options_database_string(iter,1) .eq. "vpptcrit") then
          read(current_state%options_database_string(iter,2),*) realnum
          vpptcrit = realnum
        end if
      end do
      !qlcrit    =options_get_real(current_state%options_database, "qlcrit")
      !qicrit    =options_get_real(current_state%options_database, "qicrit")
      !qpptcrit  =options_get_real(current_state%options_database, "qpptcrit")
      !vpptcrit  =options_get_real(current_state%options_database, "vpptcrit")
    endif

    !calculate qi, qr, qs and qg if and only if number_q_fields=11
    l_qi_qr_qs_qg =  (.not. current_state%passive_q .and. current_state%number_q_fields .ge. 11)

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "x_size") then
        read(current_state%options_database_string(iter,2),*) intnum
        x_size = intnum
      else if (current_state%options_database_string(iter,1) .eq. "y_size") then
        read(current_state%options_database_string(iter,2),*) intnum
        y_size = intnum
      else if (current_state%options_database_string(iter,1) .eq. "thvprcrit") then
        read(current_state%options_database_string(iter,2),*) realnum
        thvprcrit = realnum
      else if (current_state%options_database_string(iter,1) .eq. "wSdwncrit") then
        read(current_state%options_database_string(iter,2),*) realnum
        wSdwncrit = realnum
      else if (current_state%options_database_string(iter,1) .eq. "wSupcrit") then
        read(current_state%options_database_string(iter,2),*) realnum
        wSupcrit = realnum
      else if (current_state%options_database_string(iter,1) .eq. "wupcrit") then
        read(current_state%options_database_string(iter,2),*) realnum
        wupcrit = realnum
      else if (current_state%options_database_string(iter,1) .eq. "wdwncrit") then
        read(current_state%options_database_string(iter,2),*) realnum
        wdwncrit = realnum
      else if (current_state%options_database_string(iter,1) .eq. "cond_request") then
        read(current_state%options_database_string(iter,2),*) intnum
        ncond = intnum
      else if (current_state%options_database_string(iter,1) .eq. "diag_request") then
        read(current_state%options_database_string(iter,2),*) intnum
        ndiag = intnum
      else if (current_state%options_database_string(iter,1) .eq. "thv_from_th_with_liqice") then
        read(current_state%options_database_string(iter,2),*) logicnum
        thv_from_th_with_liqice = logicnum
      end if
    end do
    !x_size       =options_get_integer(current_state%options_database, "x_size")
    !y_size       =options_get_integer(current_state%options_database, "y_size")
    !thvprcrit    =options_get_real(current_state%options_database, "thvprcrit")
    !wSdwncrit    =options_get_real(current_state%options_database, "wSdwncrit")
    !wSupcrit     =options_get_real(current_state%options_database, "wSupcrit")
    !wupcrit      =options_get_real(current_state%options_database, "wupcrit")
    !wdwncrit     =options_get_real(current_state%options_database, "wdwncrit")
    !ncond        =options_get_array_size(current_state%options_database, "cond_request")
    !ndiag        =options_get_array_size(current_state%options_database, "diag_request")
    !thv_from_th_with_liqice =options_get_logical(current_state%options_database, "thv_from_th_with_liqice")

    gpts_total   = 1.0_DEFAULT_PRECISION * x_size * y_size

    !>Check for requested diagnostics.
    if (ncond .lt. 1 .or. ndiag .lt. 1) call log_master_log(LOG_ERROR, &
         "When conditional_diagnostics_column are enabled, condition and diagnostic request numbers must be provided.")

    if (current_state%th%active) then
      allocate(current_state%w_zn_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%w_zn2_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%tmp_th_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wth_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%th_pr_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wth_pr_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%thv_pr_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wthv_pr_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%th_pr2_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wthsg_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%w_zn3_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%relhum_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%tmp_u_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%tmp_v_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wu_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wv_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wusg_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wvsg_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%TdegK_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%th_h_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%th_h_pr1_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%th_h_pr2_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%qvli_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%qvli_pr_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%qvli_pr2_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%qppt_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%qppt_pr_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%qppt_pr2_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wqvli_pr_cln(current_state%local_grid%size(Z_INDEX)), &
               current_state%wqppt_pr_cln(current_state%local_grid%size(Z_INDEX)))
    end if

    !if (current_state%th%active) then
    !  allocate(CondDiags_tot(current_state%local_grid%size(Z_INDEX),  &
    !                         ncond * 2 ,  &
    !                         ndiag))
    !endif


    l_do_near = .false.
    requested_area = 0
    !> Record which conditions and diagnostics have been requested
    allocate(cond_request(ncond))
    allocate(cond_long(ncond))
    allocate(cond_locations(ncond))
    allocate(diag_request(ndiag))
    allocate(diag_long(ndiag))
    allocate(diag_locations(ndiag))

    !do inc = 1, ncond
      !cond_request(inc) = trim(options_get_string(current_state%options_database, "cond_request", inc))
      !cond_locations(inc) = inc
      !SELECT CASE ( trim(cond_request(inc)) )
      !  CASE DEFAULT
      !    call log_master_log(LOG_ERROR, &
      !        "Condition '"//trim(cond_request(inc))//"' has not been set up.")
      !  CASE ("ALL")
      !    cond_locations(inc) = 1
      !  CASE ("BYu")
      !    cond_locations(inc) = 2
      !  CASE ("BCu")
      !    cond_locations(inc) = 3
      !  CASE ("NrBCu")
      !    cond_locations(inc) = 4
      !    l_do_near = .true.
      !  CASE ("AC")
      !    cond_locations(inc) = 5
      !  CASE ("ACu")
      !    cond_locations(inc) = 6
      !  CASE ("ACd")
      !    cond_locations(inc) = 7
      !  CASE ("WG1")
      !    cond_locations(inc) = 8
      !  CASE ("WL1")
      !    cond_locations(inc) = 9
      !  CASE ("ALu")
      !    cond_locations(inc) = 10
      !  CASE ("ALd")
      !    cond_locations(inc) = 11
      !  CASE ("CLu")
      !    cond_locations(inc) = 12
      !  CASE ("CLd")
      !    cond_locations(inc) = 13
      !  CASE ("AH")
      !    cond_locations(inc) = 14
      !  CASE ("AL")
      !    cond_locations(inc) = 15
      !  CASE ("AI")
      !    cond_locations(inc) = 16
      !  CASE ("PPd")
      !    cond_locations(inc) = 17
      !  CASE ("VPd")
      !    cond_locations(inc) = 18
      !  CASE ("PVd")
      !    cond_locations(inc) = 19
      !  CASE ("MO")
      !    cond_locations(inc) = 20
      !  CASE ("BM")
      !    cond_locations(inc) = 21
      !  CASE ("AA")
      !    cond_locations(inc) = 22
      !  CASE ("AV")
      !    cond_locations(inc) = 23
      !END SELECT
    !end do ! loop over requested conditions


    !do inc = 1, ndiag
      !diag_request(inc) = trim(options_get_string(current_state%options_database, "diag_request", inc))
      !if ( trim(diag_request(inc)) .eq. "area" ) requested_area = inc
      !diag_locations(inc) = inc
      !SELECT CASE ( trim(diag_request(inc)) )
      !  CASE DEFAULT
      !    call log_master_log(LOG_ERROR, &
      !               "Diagnostic '"//trim(diag_request(inc))//"' has not been set up.")
      !  CASE ("area")
      !    diag_locations(inc) = 1
      !  CASE ("W")
      !    diag_locations(inc) = 2
      !  CASE ("W2")
      !    diag_locations(inc) = 3
      !  CASE ("TH")
      !    diag_locations(inc) = 4
      !  CASE ("WTH")
      !    diag_locations(inc) = 5
      !  CASE ("THP")
      !    diag_locations(inc) = 6
      !  CASE ("WTHP")
      !    diag_locations(inc) = 7
      !  CASE ("THVP")
      !    diag_locations(inc) = 8
      !  CASE ("WTHVP")
      !    diag_locations(inc) = 9
      !  CASE ("THP2")
      !    diag_locations(inc) = 10
      !  CASE ("WTHSG")
      !    diag_locations(inc) = 11
      !  CASE ("W3")
      !    diag_locations(inc) = 12
      !  CASE ("RH")
      !    diag_locations(inc) = 13
      !  CASE ("U")
      !    diag_locations(inc) = 14
      !  CASE ("V")
      !    diag_locations(inc) = 15
      !  CASE ("WU")
      !    diag_locations(inc) = 16
      !  CASE ("WV")
      !    diag_locations(inc) = 17
      !  CASE ("WUSG")
      !    diag_locations(inc) = 18
      !  CASE ("WVSG")
      !    diag_locations(inc) = 19
      !  CASE ("TEMP")
      !    diag_locations(inc) = 20
      !  CASE ("THL")
      !    diag_locations(inc) = 21
      !  CASE ("THLP")
      !    diag_locations(inc) = 22
      !  CASE ("THLP2")
      !    diag_locations(inc) = 23
      !  CASE ("QVLI")
      !    diag_locations(inc) = 24
      !  CASE ("QVLIP")
      !    diag_locations(inc) = 25
      !  CASE ("QVLIP2")
      !    diag_locations(inc) = 26
      !  CASE ("QRSG")
      !    diag_locations(inc) = 27
      !  CASE ("QRSGP")
      !    diag_locations(inc) = 28
      !  CASE ("QRSGP2")
      !    diag_locations(inc) = 29
      !  CASE ("WQVLIP")
      !    diag_locations(inc) = 30
      !  CASE ("WQRSGP")
      !    diag_locations(inc) = 31
      !END SELECT
    !end do ! loop over requested diagnostics

    !if (requested_area == 0) call log_master_log(LOG_ERROR, &
    !      "The diagnostic 'area' must be provided to complete the conditional diagnostics process.")

    !cond_long=master_conditions_list(cond_locations)
    !diag_long=master_diagnostics_list(diag_locations)
  end subroutine initialisation_callback_conditional_diagnostics_column


  subroutine finalisation_callback_conditional_diagnostics_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    print*,"finalisation_callback_conditional_diagnostics_column"
    if (allocated(CondDiags_tot)) deallocate(CondDiags_tot)
    if (allocated(cond_request)) deallocate(cond_request)
    if (allocated(cond_long)) deallocate(cond_long)
    if (allocated(cond_locations)) deallocate(cond_locations)
    if (allocated(diag_request)) deallocate(diag_request)
    if (allocated(diag_long)) deallocate(diag_long)
    if (allocated(diag_locations)) deallocate(diag_locations)

  end subroutine finalisation_callback_conditional_diagnostics_column


  subroutine timestep_callback_diagnostics_column(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    logical ::  l_cond                                     ! condition at local position (k,j,i)
    real(kind=DEFAULT_PRECISION), dimension(size(master_diagnostics_list)) ::  local_diag ! diagnostics at local position (k,j,i)

    real(kind=DEFAULT_PRECISION) ::  qi, ql,qv,qli,qsat,qppt,qvl
    real(kind=DEFAULT_PRECISION) ::  qi_pr, ql_pr,qv_pr,qli_pr,qvli, qvli_pr, qppt_pr
    real(kind=DEFAULT_PRECISION) ::  thv_pr,thv,exner,TdegK,th_pr,th_pr2
    real(kind=DEFAULT_PRECISION) ::  tmp_th !=thref+th
    real(kind=DEFAULT_PRECISION) ::  pottemp !=thref+olzthbar
    real(kind=DEFAULT_PRECISION) ::  relhum,Pmb,wth,wthsg, wthv_pr, wth_pr
    real(kind=DEFAULT_PRECISION) ::  w_zn ! vertical velocity at theta levels
    real(kind=DEFAULT_PRECISION) ::  tmp_u,tmp_v ! horizontal veleocities
    real(kind=DEFAULT_PRECISION) ::  w_zn2, w_zn3, wu, wv, wusg, wvsg, wqvli_pr,wqppt_pr
    real(kind=DEFAULT_PRECISION) ::  th_h,th_h_pr1,th_h_pr2 
    real(kind=DEFAULT_PRECISION) ::  qvli_pr2, qppt_pr2

    real(kind=DEFAULT_PRECISION) ::  qv_jip1,qv_jim1,qv_jp1i,qv_jm1i,qv_pr_jip1,qv_pr_jim1,qv_pr_jp1i,qv_pr_jm1i
    real(kind=DEFAULT_PRECISION) ::  qi_jip1,qi_jim1,qi_jp1i,qi_jm1i
    real(kind=DEFAULT_PRECISION) ::  ql_jip1,ql_jim1,ql_jp1i,ql_jm1i
    real(kind=DEFAULT_PRECISION) ::  qli_jip1,qli_jim1,qli_jp1i,qli_jm1i,qli_pr_jip1,qli_pr_jim1,qli_pr_jp1i,qli_pr_jm1i
    real(kind=DEFAULT_PRECISION) ::  w_zn_jip1,w_zn_jim1,w_zn_jp1i,w_zn_jm1i
    real(kind=DEFAULT_PRECISION) ::  thv_pr_jip1,thv_pr_jim1,thv_pr_jp1i,thv_pr_jm1i
    real(kind=DEFAULT_PRECISION) ::  tmp_th_jip1,tmp_th_jim1,tmp_th_jp1i,tmp_th_jm1i
    real(kind=DEFAULT_PRECISION) ::  th_pr_jip1,th_pr_jim1,th_pr_jp1i,th_pr_jm1i

    integer :: k, j, i
    integer :: inc     ! loop increment variable
    logical :: calculate_diagnostics

    !calculate_diagnostics = current_state%diagnostic_sample_timestep

    j=current_state%column_local_y
    i=current_state%column_local_x

    !> Reset sums on this process
    !if (current_state%first_timestep_column) then
    !  CondDiags_tot(:,:,:)    = 0.0_DEFAULT_PRECISION
    !end if

    if (current_state%first_timestep_column) then
      current_state%w_zn_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%w_zn2_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%tmp_th_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wth_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%th_pr_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wth_pr_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%thv_pr_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wthv_pr_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%th_pr2_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wthsg_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%w_zn3_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%relhum_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%tmp_u_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%tmp_v_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wu_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wv_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wusg_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wvsg_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%TdegK_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%th_h_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%th_h_pr1_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%th_h_pr2_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%qvli_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%qvli_pr_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%qvli_pr2_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%qppt_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%qppt_pr_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%qppt_pr2_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wqvli_pr_cln(:) = 0.0_DEFAULT_PRECISION
      current_state%wqppt_pr_cln(:) = 0.0_DEFAULT_PRECISION
    end if

    !> Decide if conditions are appropriate to proceed with calculations
    if (current_state%halo_column) return
    if ( .not. (current_state%th%active .and.       &
                .not. current_state%passive_q .and. &
                current_state%number_q_fields .gt. 0) ) return

    if (current_state%modulo_number_1d .ne. 0) return

    !> Begin the calculations
    !> Loop over levels
    !  Due to above/below level solution dependency, no data will exist at first/last levels,
    !  but k-dimension size kept equal to size(Z_INDEX) for convenience
    do k = 2, current_state%local_grid%size(Z_INDEX)-1

      !> work out all potential terms at this location
      qv      = current_state%zqv%data(k,j,i)
      ql      = current_state%zql%data(k,j,i)
      qvl     = qv + ql
      qv_pr   = qv - current_state%global_grid%configuration%vertical%olzqvbar(k)
      ql_pr   = ql - current_state%global_grid%configuration%vertical%olzqlbar(k)

      if (.not. l_qi_qr_qs_qg) then
        qi      = 0.0
        qli     = ql
        qvli    = qvl
        qppt    = 0.0
        qi_pr   = 0.0
        qli_pr  = ql_pr
        qvli_pr = qv_pr + ql_pr
        qppt_pr = 0.0
      end if

      if (l_qi_qr_qs_qg) then
        qi      = current_state%zqi%data(k,j,i)
        qli     = ql + qi
        qvli    = qv + qli
        qppt    = current_state%zqr%data(k,j,i) + current_state%zqs%data(k,j,i)+current_state%zqg%data(k,j,i)
        qi_pr   = qli - current_state%global_grid%configuration%vertical%olzqibar(k)
        qli_pr  = ql_pr + qi_pr
        qvli_pr = qv_pr + qli_pr
        qppt_pr = current_state%zqr%data(k,j,i) - current_state%global_grid%configuration%vertical%olzqrbar(k)+   &
                  current_state%zqs%data(k,j,i) - current_state%global_grid%configuration%vertical%olzqsbar(k)+   &
                  current_state%zqg%data(k,j,i) - current_state%global_grid%configuration%vertical%olzqgbar(k)
      end if

      w_zn     = 0.5 * (current_state%zw%data(k,j,i) + current_state%zw%data(k-1,j,i)) !w at theta levels
      w_zn2    = w_zn * w_zn
      w_zn3    = w_zn * w_zn * w_zn

      qvli_pr2 = qvli_pr * qvli_pr
      qppt_pr2 = qppt_pr * qppt_pr
      wqvli_pr = w_zn * qvli_pr
      wqppt_pr = w_zn * qppt_pr

      tmp_u    = 0.5 * (current_state%zu%data(k,j,i-1) + current_state%zu%data(k,j,i)) + current_state%ugal !u at theta points
      tmp_v    = 0.5 * (current_state%zv%data(k,j-1,i) + current_state%zv%data(k,j,i)) + current_state%vgal !v at theta points
      wu       = w_zn * tmp_u
      wv       = w_zn * tmp_v
      tmp_th   = current_state%zth%data(k,j,i) ! + current_state%global_grid%configuration%vertical%thref(k) +
      wth      = w_zn * tmp_th
      wthsg    = -0.5 * current_state%diff_coefficient%data(k,j,i) * (current_state%zth%data(k+1,j,i) +          &
                 current_state%global_grid%configuration%vertical%thref(k+1) - current_state%zth%data(k,j,i) -   &
                 current_state%global_grid%configuration%vertical%thref(k)) *                                    &
                 current_state%global_grid%configuration%vertical%rdzn(k+1)                                      &
                 -0.5 * current_state%diff_coefficient%data(k-1,j,i) * (current_state%zth%data(k,j,i) +          &
                 current_state%global_grid%configuration%vertical%thref(k) - current_state%zth%data(k-1,j,i) -   &
                 current_state%global_grid%configuration%vertical%thref(k-1)) *                                  &
                 current_state%global_grid%configuration%vertical%rdzn(k)

      wusg     = -0.25 * current_state%vis_coefficient%data(k  ,j,i-1) *                      &
                 (current_state%zu%data(k+1,j,i-1) - current_state%zu%data(k  ,j,i-1)) *      &
                 current_state%global_grid%configuration%vertical%rdzn(k+1)                   &
                 -0.25 * current_state%vis_coefficient%data(k  ,j,i  ) *                      &
                 (current_state%zu%data(k+1,j,i  ) - current_state%zu%data(k  ,j,i  )) *      &
                 current_state%global_grid%configuration%vertical%rdzn(k+1)                   &
                 -0.25 * current_state%vis_coefficient%data(k-1,j,i-1) *                      &
                 (current_state%zu%data(k  ,j,i-1) - current_state%zu%data(k-1,j,i-1)) *      &
                 current_state%global_grid%configuration%vertical%rdzn(k)                     &
                 -0.25 * current_state%vis_coefficient%data(k-1,j,i  ) *                      &
                 (current_state%zu%data(k  ,j,i  ) - current_state%zu%data(k-1,j,i  )) *      &
                 current_state%global_grid%configuration%vertical%rdzn(k)  

      wvsg     = -0.25 * current_state%vis_coefficient%data(k  ,j-1,i) *                      &
                 (current_state%zv%data(k+1,j-1,i) - current_state%zv%data(k  ,j-1,i)) *      &
                 current_state%global_grid%configuration%vertical%rdzn(k+1)                   &
                 -0.25 * current_state%vis_coefficient%data(k  ,j  ,i) *                      &
                 (current_state%zv%data(k+1,j  ,i) - current_state%zv%data(k  ,j  ,i)) *      &
                 current_state%global_grid%configuration%vertical%rdzn(k+1)                   &
                 -0.25 * current_state%vis_coefficient%data(k-1,j-1,i) *                      &
                 (current_state%zv%data(k  ,j-1,i) - current_state%zv%data(k-1,j-1,i)) *      &
                 current_state%global_grid%configuration%vertical%rdzn(k)                     &
                 -0.25 * current_state%vis_coefficient%data(k-1,j  ,i) *                      &
                 (current_state%zv%data(k  ,j  ,i) - current_state%zv%data(k-1,j  ,i)) *      &
                 current_state%global_grid%configuration%vertical%rdzn(k)

      pottemp  = current_state%global_grid%configuration%vertical%thref(k) +    &
                 current_state%global_grid%configuration%vertical%olzthbar(k)

      th_pr    = current_state%zth%data(k,j,i) - current_state%global_grid%configuration%vertical%olzthbar(k)
      wth_pr   = w_zn * th_pr
      th_pr2   = th_pr * th_pr
      exner    = current_state%global_grid%configuration%vertical%rprefrcp(k)
      TdegK    = (current_state%global_grid%configuration%vertical%thref(k) + current_state%zth%data(k,j,i)) * exner
      Pmb      = current_state%global_grid%configuration%vertical%prefn(k) / 100.
      qsat     = qsaturation(TdegK,Pmb)
      relhum   = qv / qsat
      th_h     = tmp_th + rlvap * qv / cp
      th_h_pr1 = th_pr + rlvap * qv_pr / cp
      th_h_pr2 = th_h_pr1 * th_h_pr1

      !calculate values at nearby points
      if (l_do_near) then
        qv_jip1      = current_state%zqv%data(k,j,i+1) !qv at the position j and i+1
        qv_jim1      = current_state%zqv%data(k,j,i-1) !qv at the position j and i-1
        qv_jp1i      = current_state%zqv%data(k,j+1,i) !qv at the position j+1 and i
        qv_jm1i      = current_state%zqv%data(k,j-1,i) !qv at the position j-1 and i
        ql_jip1      = current_state%zql%data(k,j,i+1)
        ql_jim1      = current_state%zql%data(k,j,i-1)
        ql_jp1i      = current_state%zql%data(k,j+1,i)
        ql_jm1i      = current_state%zql%data(k,j-1,i)
        qv_pr_jip1   = current_state%zqv%data(k,j,i+1)-current_state%global_grid%configuration%vertical%olzqvbar(k)
        qv_pr_jim1   = current_state%zqv%data(k,j,i-1)-current_state%global_grid%configuration%vertical%olzqvbar(k)
        qv_pr_jp1i   = current_state%zqv%data(k,j+1,i)-current_state%global_grid%configuration%vertical%olzqvbar(k)
        qv_pr_jm1i   = current_state%zqv%data(k,j-1,i)-current_state%global_grid%configuration%vertical%olzqvbar(k)

        if (.not. l_qi_qr_qs_qg) then
          qi_jip1      = 0.0
          qi_jim1      = 0.0
          qi_jp1i      = 0.0
          qi_jm1i      = 0.0
          qli_jip1     = current_state%zql%data(k,j,i+1)+0.0
          qli_jim1     = current_state%zql%data(k,j,i-1)+0.0
          qli_jp1i     = current_state%zql%data(k,j+1,i)+0.0
          qli_jm1i     = current_state%zql%data(k,j-1,i)+0.0
          qli_pr_jip1  = current_state%zql%data(k,j,i+1)-current_state%global_grid%configuration%vertical%olzqlbar(k)+0.0
          qli_pr_jim1  = current_state%zql%data(k,j,i-1)-current_state%global_grid%configuration%vertical%olzqlbar(k)+0.0
          qli_pr_jp1i  = current_state%zql%data(k,j+1,i)-current_state%global_grid%configuration%vertical%olzqlbar(k)+0.0
          qli_pr_jm1i  = current_state%zql%data(k,j-1,i)-current_state%global_grid%configuration%vertical%olzqlbar(k)+0.0
        end if

        if (l_qi_qr_qs_qg) then
          qi_jip1      = current_state%zqi%data(k,j,i+1)
          qi_jim1      = current_state%zqi%data(k,j,i-1)
          qi_jp1i      = current_state%zqi%data(k,j+1,i)
          qi_jm1i      = current_state%zqi%data(k,j-1,i)
          qli_jip1     = current_state%zql%data(k,j,i+1)+current_state%zqi%data(k,j,i+1)
          qli_jim1     = current_state%zql%data(k,j,i-1)+current_state%zqi%data(k,j,i-1)
          qli_jp1i     = current_state%zql%data(k,j+1,i)+current_state%zqi%data(k,j+1,i)
          qli_jm1i     = current_state%zql%data(k,j-1,i)+current_state%zqi%data(k,j-1,i)
          qli_pr_jip1  = current_state%zql%data(k,j,i+1)-current_state%global_grid%configuration%vertical%olzqlbar(k)+ &
                         current_state%zqi%data(k,j,i+1)-current_state%global_grid%configuration%vertical%olzqibar(k)
          qli_pr_jim1  = current_state%zql%data(k,j,i-1)-current_state%global_grid%configuration%vertical%olzqlbar(k)+ &
                         current_state%zqi%data(k,j,i-1)-current_state%global_grid%configuration%vertical%olzqibar(k)
          qli_pr_jp1i  = current_state%zql%data(k,j+1,i)-current_state%global_grid%configuration%vertical%olzqlbar(k)+ &
                         current_state%zqi%data(k,j+1,i)-current_state%global_grid%configuration%vertical%olzqibar(k)
          qli_pr_jm1i  = current_state%zql%data(k,j-1,i)-current_state%global_grid%configuration%vertical%olzqlbar(k)+ &
                         current_state%zqi%data(k,j-1,i)-current_state%global_grid%configuration%vertical%olzqibar(k)
        end if

        tmp_th_jip1  = current_state%global_grid%configuration%vertical%thref(k) + current_state%zth%data(k,j,i+1) 
        tmp_th_jim1  = current_state%global_grid%configuration%vertical%thref(k) + current_state%zth%data(k,j,i-1) 
        tmp_th_jp1i  = current_state%global_grid%configuration%vertical%thref(k) + current_state%zth%data(k,j+1,i) 
        tmp_th_jm1i  = current_state%global_grid%configuration%vertical%thref(k) + current_state%zth%data(k,j-1,i) 
        th_pr_jip1   = current_state%zth%data(k,j,i+1)-current_state%global_grid%configuration%vertical%olzthbar(k)
        th_pr_jim1   = current_state%zth%data(k,j,i-1)-current_state%global_grid%configuration%vertical%olzthbar(k)
        th_pr_jp1i   = current_state%zth%data(k,j-1,i)-current_state%global_grid%configuration%vertical%olzthbar(k)
        th_pr_jm1i   = current_state%zth%data(k,j+1,i)-current_state%global_grid%configuration%vertical%olzthbar(k)
        w_zn_jip1    = 0.5*(current_state%zw%data(k,j,i+1)+current_state%zw%data(k-1,j,i+1)) 
        w_zn_jim1    = 0.5*(current_state%zw%data(k,j,i-1)+current_state%zw%data(k-1,j,i-1)) 
        w_zn_jp1i    = 0.5*(current_state%zw%data(k,j+1,i)+current_state%zw%data(k-1,j+1,i)) 
        w_zn_jm1i    = 0.5*(current_state%zw%data(k,j-1,i)+current_state%zw%data(k-1,j-1,i)) 

        if (thv_from_th_with_liqice) then  
          thv       = tmp_th*(1.0+0.61*qv-qli)
          !the term th_pr*(0.61*qv-qli) is very small and can be neglected
          thv_pr    = th_pr+tmp_th*(0.61*qv_pr-qli_pr)+th_pr*(0.61*qv-qli)
          wthv_pr   = w_zn*thv_pr

          thv_pr_jip1   =  th_pr_jip1+tmp_th_jip1*(0.61*qv_pr_jip1-qli_pr_jip1)+th_pr_jip1*(0.61*qv_jip1-qli_jip1)
          thv_pr_jim1   =  th_pr_jim1+tmp_th_jim1*(0.61*qv_pr_jim1-qli_pr_jim1)+th_pr_jim1*(0.61*qv_jim1-qli_jim1)
          thv_pr_jp1i   =  th_pr_jp1i+tmp_th_jp1i*(0.61*qv_pr_jp1i-qli_pr_jp1i)+th_pr_jp1i*(0.61*qv_jp1i-qli_jp1i)
          thv_pr_jm1i   =  th_pr_jm1i+tmp_th_jm1i*(0.61*qv_pr_jm1i-qli_pr_jm1i)+th_pr_jm1i*(0.61*qv_jm1i-qli_jm1i)
        end if

        if (.not. thv_from_th_with_liqice) then  
          thv      = tmp_th*(1.0+0.61*qv)
          ! the 0.61*qv*th_pr is very small and can be neglected  
          thv_pr   = th_pr + 0.61 * tmp_th * qv_pr + 0.61 * qv * th_pr 
          wthv_pr  = w_zn * thv_pr
          thv_pr_jip1   = th_pr_jip1+0.61*tmp_th_jip1*qv_pr_jip1+0.61*qv_jip1*th_pr_jip1 
          thv_pr_jim1   = th_pr_jim1+0.61*tmp_th_jim1*qv_pr_jim1+0.61*qv_jim1*th_pr_jim1 
          thv_pr_jp1i   = th_pr_jp1i+0.61*tmp_th_jp1i*qv_pr_jp1i+0.61*qv_jp1i*th_pr_jp1i 
          thv_pr_jm1i   = th_pr_jm1i+0.61*tmp_th_jm1i*qv_pr_jm1i+0.61*qv_jm1i*th_pr_jm1i 
        end if 

      end if ! l_do_near
      !> end of potential term calculations


      !> Store the diagnostics
      !local_diag(1) = 1.0_DEFAULT_PRECISION
      !local_diag(2) = w_zn
      !local_diag(3) = w_zn2
      !local_diag(4) = tmp_th
      !local_diag(5) = wth
      !local_diag(6) = th_pr
      !local_diag(7) = wth_pr
      !local_diag(8) = thv_pr
      !local_diag(9) = wthv_pr
      !local_diag(10) = th_pr2
      !local_diag(11) = wthsg
      !local_diag(12) = w_zn3
      !local_diag(13) = relhum
      !local_diag(14) = tmp_u
      !local_diag(15) = tmp_v
      !local_diag(16) = wu
      !local_diag(17) = wv
      !local_diag(18) = wusg
      !local_diag(19) = wvsg
      !local_diag(20) = TdegK
      !local_diag(21) = th_h
      !local_diag(22) = th_h_pr1
      !local_diag(23) = th_h_pr2
      !local_diag(24) = qvli
      !local_diag(25) = qvli_pr
      !local_diag(26) = qvli_pr2
      !local_diag(27) = qppt
      !local_diag(28) = qppt_pr
      !local_diag(29) = qppt_pr2
      !local_diag(30) = wqvli_pr
      !local_diag(31) = wqppt_pr

      !> Determine status of requested conditions
      !do inc = 1, ncond

      !  l_cond = .false.  ! initialization

        !Condition processing on requested condition
        !SELECT CASE ( trim(cond_request(inc)) )

          !CASE DEFAULT
          !  call log_master_log(LOG_ERROR, &
          !        "Condition '"//trim(cond_request(inc))//"' has not been set up.")

          !CASE ("ALL")
!             l_cond = .true.
!
!           !CASE ("BYu")
!             l_cond = ( (w_zn .gt. wupcrit) .and. (thv_pr .gt. thvprcrit) )
!
!           !CASE ("BCu")
!             l_cond = ( (w_zn .gt. wupcrit .and. thv_pr .gt. thvprcrit ) .and.      &
!                        ( ql .gt. qlcrit .or. qi .gt. qicrit ) )
!
!           !CASE ("NrBCu")
!             !> Near it, but not it.
!             if ( .not. ( (w_zn .gt. wupcrit .and. thv_pr .gt. thvprcrit ) .and.      &
!                        ( ql .gt. qlcrit .or. qi .gt. qicrit ) )    ) then
!
!               l_cond = (  ( (w_zn_jip1 .gt. wupcrit .and. thv_pr_jip1 .gt. thvprcrit) .and.   &
!                             (ql_jip1 .gt. qlcrit .or. qi_jip1 .gt. qicrit) )                  &
!                           .or.                                                                &
!                           ( (w_zn_jim1 .gt. wupcrit .and. thv_pr_jim1 .gt. thvprcrit) .and.   &
!                             (ql_jim1 .gt. qlcrit .or. qi_jim1 .gt. qicrit) )                  &
!                           .or.                                                                &
!                           ( (w_zn_jp1i .gt. wupcrit .and. thv_pr_jp1i .gt. thvprcrit) .and.   &
!                             (ql_jp1i .gt. qlcrit .or. qi_jp1i .gt. qicrit) )                  &
!                           .or.                                                                &
!                           ( (w_zn_jm1i .gt. wupcrit .and. thv_pr_jm1i .gt. thvprcrit) .and.   &
!                             (ql_jm1i .gt. qlcrit .or. qi_jm1i .gt. qicrit) )  )
!
!             end if
!
!           !CASE ("AC")
!             l_cond = ( (ql .gt. qlcrit) .or. (qi .gt. qicrit) )
!
!           !CASE ("ACu")
!             l_cond = ( (w_zn .gt. wupcrit) .and.          &
!                        ( (ql .gt. qlcrit) .or. (qi .gt. qicrit) ) )
!
!           !CASE ("ACd")
!             l_cond = ( (w_zn .lt. wdwncrit) .and.         &
!                        ( (ql .gt. qlcrit) .or. (qi .gt. qicrit) ) )
!
!           !CASE ("WG1")
!             l_cond = ( w_zn .gt. wSupcrit )
!
!           !CASE ("WL1")
!             l_cond = ( w_zn .lt. wSdwncrit )
!
!           !CASE ("ALu")
!             l_cond = ( w_zn .gt. wupcrit )
!
!           !CASE ("ALd")
!             l_cond = ( w_zn .lt. wdwncrit )
!
!           !CASE ("CLu")
!             l_cond = ( (w_zn .gt. wupcrit) .and. (ql .gt. qlcrit) )
!
!           !CASE ("CLd")
!             l_cond = ( (w_zn .lt. wdwncrit) .and. (ql .gt. qlcrit) )
!
!           !CASE ("AH")
!             l_cond = ( (ql .gt. qlcrit) .or. (qi .gt. qicrit) .or. (qppt .gt. qpptcrit) )
!
!           !CASE ("AL")
!             l_cond = ( ql .gt. qlcrit )
!
!           !CASE ("AI")
!             l_cond = ( qi .gt. qicrit )
!
!           !CASE ("PPd")
!             l_cond = ( (w_zn .lt. wdwncrit) .and. (qppt .gt. qpptcrit) )
!
!           !CASE ("VPd")
!             l_cond = ( (w_zn .lt. wdwncrit) .and. (qppt .gt. vpptcrit) )
!
!           !CASE ("PVd")
!             l_cond = ( (w_zn .lt. wSdwncrit) .and. (qppt .gt. qpptcrit) )
!
!           !CASE ("MO")
!             l_cond = ( (w_zn .gt. wupcrit) .and. (qvl .gt. qsat) )
!
!           !CASE ("BM")
!             l_cond = ( (w_zn .gt. wupcrit) .and. (thv_pr .gt. thvprcrit) .and. (qvl .gt. qsat) )
!
!           !CASE ("AA")
!             l_cond = ( (ql .lt. qlcrit) .and. (qi .lt. qicrit) .and.   &
!                        (w_zn .lt. wdwncrit) .and. (qppt .gt. qpptcrit)  )
!
!           !CASE ("AV")
!             l_cond = ( (ql .lt. qlcrit) .and. (qi .lt. qicrit) .and.   &
!                        (w_zn .lt. wdwncrit) .and. (qppt .gt. vpptcrit)  )

        !END SELECT

            
        !> Record the value for this diagnostic meeting this condition and its complement by adding to cumulative sum
        !if (l_cond) then
        !  CondDiags_tot(k,inc,:) = CondDiags_tot(k,inc,:) + local_diag(diag_locations)
        !else
        !  CondDiags_tot(k,ncond+inc,:) = CondDiags_tot(k,ncond+inc,:) + local_diag(diag_locations)
        !end if

      !end do ! inc loop over requested conditions
      current_state%w_zn_cln(k) = current_state%w_zn_cln(k) + w_zn
      current_state%w_zn2_cln(k) = current_state%w_zn2_cln(k) + w_zn2
      current_state%tmp_th_cln(k) = current_state%tmp_th_cln(k) + tmp_th
      current_state%wth_cln(k) = current_state%wth_cln(k) + wth
      current_state%th_pr_cln(k) = current_state%th_pr_cln(k) + th_pr
      current_state%wth_pr_cln(k) = current_state%wth_pr_cln(k) + wth_pr
      current_state%thv_pr_cln(k) = current_state%thv_pr_cln(k) + thv_pr
      current_state%wthv_pr_cln(k) = current_state%wthv_pr_cln(k) + wthv_pr
      current_state%th_pr2_cln(k) = current_state%th_pr2_cln(k) + th_pr2
      current_state%wthsg_cln(k) = current_state%wthsg_cln(k) + wthsg
      current_state%w_zn3_cln(k) = current_state%w_zn3_cln(k) + w_zn3
      current_state%relhum_cln(k) = current_state%relhum_cln(k) + relhum
      current_state%tmp_u_cln(k) = current_state%tmp_u_cln(k) + tmp_u
      current_state%tmp_v_cln(k) = current_state%tmp_v_cln(k) + tmp_v
      current_state%wu_cln(k) = current_state%wu_cln(k) + wu
      current_state%wv_cln(k) = current_state%wv_cln(k) + wv
      current_state%wusg_cln(k) = current_state%wusg_cln(k) + wusg
      current_state%wvsg_cln(k) = current_state%wvsg_cln(k) + wvsg
      current_state%TdegK_cln(k) = current_state%TdegK_cln(k) + TdegK
      current_state%th_h_cln(k) = current_state%th_h_cln(k) + th_h
      current_state%th_h_pr1_cln(k) = current_state%th_h_pr1_cln(k) + th_h_pr1
      current_state%th_h_pr2_cln(k) = current_state%th_h_pr2_cln(k) + th_h_pr2
      current_state%qvli_cln(k) = current_state%qvli_cln(k) + qvli
      current_state%qvli_pr_cln(k) = current_state%qvli_pr_cln(k) + qvli_pr
      current_state%qvli_pr2_cln(k) = current_state%qvli_pr2_cln(k) + qvli_pr2
      current_state%qppt_cln(k) = current_state%qppt_cln(k) + qppt
      current_state%qppt_pr_cln(k) = current_state%qppt_pr_cln(k) + qppt_pr
      current_state%qppt_pr2_cln(k) = current_state%qppt_pr2_cln(k) + qppt_pr2
      current_state%wqvli_pr_cln(k) = current_state%wqvli_pr_cln(k) + wqvli_pr
      current_state%wqppt_pr_cln(k) = current_state%wqppt_pr_cln(k) + wqppt_pr
    end do ! k loop over vertical levels

  end subroutine timestep_callback_diagnostics_column



end module conditional_diagnostics_column_mod
