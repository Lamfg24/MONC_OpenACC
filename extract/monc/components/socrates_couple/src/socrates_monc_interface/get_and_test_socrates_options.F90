module get_and_test_socrates_options_mod

 use datadefn_mod, only : DEFAULT_PRECISION
 use state_mod, only : model_state_type
 use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, &
       LOG_DEBUG, log_master_log, log_log, log_get_logging_level, &
       log_master_log
 use conversions_mod, only : conv_to_string
 use q_indices_mod, only: get_q_index, standard_q_names
 use optionsdatabase_mod, only : options_get_string, options_get_integer, &
       options_get_real, options_get_logical, options_has_key
 
 use def_socrates_options, only: str_socrates_options
 

contains
  
  subroutine set_and_test_socrates_monc_options(current_state, socrates_opt)

    type(model_state_type), target, intent(inout) :: current_state

    type (str_socrates_options), intent(inout) :: socrates_opt

    integer ::  intnum, iter
    real(kind=DEFAULT_PRECISION) :: realnum
    logical :: logicnum, casim_enabled, simplecloud_enabled, rcemip_gases_enabled
    
    ! First check the MONC switches are sensible to run socrates
    !
    ! check that lwrad_exp is not on, as this will double the cloud top
    ! longwave cooling
    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "lwrad_exponential_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (logicnum .eqv. .true.) then
          call log_master_log &
            (LOG_ERROR, "Socrates and lwrad_exponential both enabled, switch off on in config - STOP")
        end if
    ! determine if the radiation is expecting cloud
      else if (current_state%options_database_string(iter,1) .eq. "i_cloud_representation") then
        read(current_state%options_database_string(iter,2),*) intnum
        socrates_opt%cloud_representation = intnum
      else if (current_state%options_database_string(iter,1) .eq. "casim_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        casim_enabled = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "casim_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        simplecloud_enabled = logicnum
      endif
    end do
    
    ! check whether potential temperature is active, if not stop run
    if (.not. current_state%th%active) then
       call log_master_log &
            (LOG_ERROR, "Socrates called with theta inactive, this is not supported - STOP")
    endif

    if (socrates_opt%cloud_representation == 5) then
       ! only clear sky radiation calculation, initialise vapour index
       if (current_state%number_q_fields < 1) then 
          call log_master_log(LOG_ERROR, "Socrates called for clear sky calc but no vapour field - STOP")
       !else
       !  socrates_opt%iqv=get_q_index(standard_q_names%VAPOUR, 'socrates_couple')
       endif
    else if (socrates_opt%cloud_representation == 2 .or. socrates_opt%cloud_representation == 1) then
       ! Do some tests to make sure config is sensible
       ! (these tests are a bit overkill, better safe than sorry)
       if (.not. casim_enabled .and. .not. simplecloud_enabled) then
          call log_master_log &
               (LOG_ERROR, "Socrates called for cloudy sky but no microphysics scheme enabled - STOP")
       endif
       if (current_state%passive_q ) then
          call log_master_log &
               (LOG_ERROR, "Socrates called for cloudy sky but q is passive so not cloud or vapour - STOP")
       endif
       if (current_state%number_q_fields < 2) then 
          call log_master_log &
               (LOG_ERROR, "Socrates called for clear and cloud sky calc but no vapour or cloud field - STOP")
       endif
       ! Now monc is happy the set-up seems OK, set-up the microphysics logicals and allocate
       ! arrays
       !
       ! set vapour by default
       !socrates_opt%iqv=get_q_index(standard_q_names%VAPOUR, 'socrates_couple')
       ! read from options database to see which hydrometeors are available
       ! for radiation calculation. NOTE: This can differ from the number of hydrometeors
       ! active in the microphysics scheme
       do iter = 1,current_state%config_args
        if (current_state%options_database_string(iter,1) .eq. "mphys_nq_l") then
                read(current_state%options_database_string(iter,2),*) intnum
                socrates_opt%mphys_nq_l = intnum
        else if (current_state%options_database_string(iter,1) .eq. "mphys_nd_l") then
                read(current_state%options_database_string(iter,2),*) intnum
                socrates_opt%mphys_nd_l = intnum
        else if (current_state%options_database_string(iter,1) .eq. "mphys_nq_r") then
                read(current_state%options_database_string(iter,2),*) intnum
                socrates_opt%mphys_nq_r = intnum
        else if (current_state%options_database_string(iter,1) .eq. "mphys_nq_i") then
                read(current_state%options_database_string(iter,2),*) intnum
                socrates_opt%mphys_nq_i = intnum
        else if (current_state%options_database_string(iter,1) .eq. "mphys_nq_s") then
                read(current_state%options_database_string(iter,2),*) intnum
                socrates_opt%mphys_nq_s = intnum
        else if (current_state%options_database_string(iter,1) .eq. "mphys_nq_g") then
                read(current_state%options_database_string(iter,2),*) intnum
                socrates_opt%mphys_nq_g = intnum
        else if (current_state%options_database_string(iter,1) .eq. "l_fix_re") then
                read(current_state%options_database_string(iter,2),*) logicnum
                socrates_opt%l_fix_re = logicnum
        else if (current_state%options_database_string(iter,1) .eq. "l_use_ndrop") then
                read(current_state%options_database_string(iter,2),*) logicnum
                socrates_opt%l_use_ndrop = logicnum
        else if (current_state%options_database_string(iter,1) .eq. "fixed_cloud_re") then
                read(current_state%options_database_string(iter,2),*) realnum
                socrates_opt%fixed_cloud_re = realnum
        else if (current_state%options_database_string(iter,1) .eq. "fixed_ice_re") then
                read(current_state%options_database_string(iter,2),*) realnum
                socrates_opt%fixed_ice_re = realnum
        else if (current_state%options_database_string(iter,1) .eq. "rho_water") then
                read(current_state%options_database_string(iter,2),*) realnum
                socrates_opt%rho_water = realnum
        else if (current_state%options_database_string(iter,1) .eq. "kparam") then
                read(current_state%options_database_string(iter,2),*) realnum
                socrates_opt%kparam = realnum
        else if (current_state%options_database_string(iter,1) .eq. "l_use_liu_spec") then
                read(current_state%options_database_string(iter,2),*) logicnum
                socrates_opt%l_use_liu_spec = logicnum
        else if (current_state%options_database_string(iter,1) .eq. "fixed_cloud_number") then
                read(current_state%options_database_string(iter,2),*) realnum
                socrates_opt%fixed_cloud_number = realnum
        end if
       end do
       !if (socrates_opt%mphys_nq_l > 0) &
       !     socrates_opt%iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'socrates_couple')
       !if (socrates_opt%mphys_nq_r > 0) &
       !     socrates_opt%iqr=get_q_index(standard_q_names%RAIN_MASS, 'socrates_couple')
       !if (socrates_opt%mphys_nq_i > 0) &
       !     socrates_opt%iqi=get_q_index(standard_q_names%ICE_MASS, 'socrates_couple')
       !if (socrates_opt%mphys_nq_s > 0) &
       !     socrates_opt%iqs=get_q_index(standard_q_names%SNOW_MASS, 'socrates_couple')
       !if (socrates_opt%mphys_nq_g > 0) &
       !     socrates_opt%iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'socrates_couple')
       ! test for fixed effective radius in config
       if ((socrates_opt%l_fix_re .and. socrates_opt%l_use_ndrop) .or. &
            (.not. socrates_opt%l_fix_re .and. .not. socrates_opt%l_use_ndrop)) then
          call log_master_log &
                     (LOG_ERROR, "Socrates - l_fix_re and l_use_ndrop both true or both false, please pick one - STOP")
       endif
       if (socrates_opt%l_fix_re) then
          call log_master_log &
               (LOG_INFO, "Socrates - using fixed cloud effective radius"//trim(conv_to_string(socrates_opt%fixed_cloud_re)))
          call log_master_log &
               (LOG_INFO, "Socrates - using fixed ice effective radius"//trim(conv_to_string(socrates_opt%fixed_ice_re)))
       endif
       if (socrates_opt%l_use_ndrop) then
          if (socrates_opt%mphys_nd_l > 0) then
             if (casim_enabled .eqv. .true.) then
                !socrates_opt%inl = get_q_index(standard_q_names%CLOUD_LIQUID_NUMBER, 'socrates_couple')
                call log_master_log &
                     (LOG_INFO, "Socrates using prognostic cloud number from CASIM, fixed ice re="&
                     //trim(conv_to_string(socrates_opt%fixed_ice_re))//" microns")
             else
                call log_master_log &
                     (LOG_ERROR, "Socrates - casim not enabled so no prognostic nd - STOP")
             endif
          else
             !socrates_opt%inl = 0
             call log_master_log &
                  (LOG_INFO, "Socrates using prescribed fixed_cloud_number="&
                  //trim(conv_to_string(socrates_opt%fixed_cloud_number))//" /m**3")
          endif
       endif
    else
       call log_master_log &
            (LOG_ERROR, "Socrates config using unrecognised i_cloud_representation, check config - STOP")    
    endif
    
    ! Get options for time and location variables. These are set in the MONC
    ! configuration. By default they are set to -999.0, which will fail. Hence
    ! user must set them appropriately
    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "l_360") then
        read(current_state%options_database_string(iter,2),*) logicnum
        socrates_opt%l_360 = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "l_solar_fixed") then
        read(current_state%options_database_string(iter,2),*) logicnum
        socrates_opt%l_solar_fixed = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "l_no_solar") then
        read(current_state%options_database_string(iter,2),*) logicnum
        socrates_opt%l_no_solar = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "solar_fixed") then
        read(current_state%options_database_string(iter,2),*) realnum
        socrates_opt%solar_fixed = realnum
      else if (current_state%options_database_string(iter,1) .eq. "sec_fixed") then
        read(current_state%options_database_string(iter,2),*) realnum
        socrates_opt%sec_fixed = realnum
      else if (current_state%options_database_string(iter,1) .eq. "l_variable_srf_albedo") then
        read(current_state%options_database_string(iter,2),*) logicnum
        socrates_opt%l_variable_srf_albedo = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "surface_albedo") then
        read(current_state%options_database_string(iter,2),*) realnum
        socrates_opt%surface_albedo = realnum
      else if (current_state%options_database_string(iter,1) .eq. "rad_start_year") then
        read(current_state%options_database_string(iter,2),*) realnum
        socrates_opt%rad_year = realnum
      else if (current_state%options_database_string(iter,1) .eq. "rad_start_day") then
        read(current_state%options_database_string(iter,2),*) realnum
        socrates_opt%rad_start_day = realnum
      else if (current_state%options_database_string(iter,1) .eq. "rad_start_time") then
        read(current_state%options_database_string(iter,2),*) realnum
        socrates_opt%rad_start_time = realnum
      else if (current_state%options_database_string(iter,1) .eq. "latitude") then
        read(current_state%options_database_string(iter,2),*) realnum
        socrates_opt%latitude = realnum
      else if (current_state%options_database_string(iter,1) .eq. "longitude") then
        read(current_state%options_database_string(iter,2),*) realnum
        socrates_opt%longitude = realnum
      ! Read the radiation call interval
      else if (current_state%options_database_string(iter,1) .eq. "rad_interval") then
        read(current_state%options_database_string(iter,2),*) intnum
        socrates_opt%rad_interval = intnum
      else if (current_state%options_database_string(iter,1) .eq. "default_solar_constant") then
        read(current_state%options_database_string(iter,2),*) realnum
        socrates_opt%default_solar_constant = realnum
      else if (current_state%options_database_string(iter,1) .eq. "l_rcemip_gases") then
        read(current_state%options_database_string(iter,2),*) logicnum
        rcemip_gases_enabled = logicnum
      end if
    end do

    socrates_opt%l_rad_calc = .false.

    if (socrates_opt%l_solar_fixed) then
       if (socrates_opt%solar_fixed .lt. -1.0 .or. socrates_opt%sec_fixed .lt. -1.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - l_solar_fixed but solar_fixed and/or sec_fixed not set, check config - STOP")
       endif
       if (socrates_opt%l_variable_srf_albedo) then
          call log_master_log &
               (LOG_INFO, "Socrates - solar_fixed but variable srf albedo, default to fixed srf albedo")
       endif
    else ! solar angle will vary so get all the time variables
       ! first get the radiation initial year, day and time variables and add to
       if (socrates_opt%rad_year < 0.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - start year is negative, which is wrong, check config - STOP")
       endif
       if (socrates_opt%rad_start_day < 0.0 .or. socrates_opt%rad_start_day > 360.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - start day is outside sensible range, check config - STOP")
       endif
       if (socrates_opt%rad_time_hours < 0.0 .or. socrates_opt%rad_time_hours > 24.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - start time is outside sensible range, check config - STOP")
       endif
       ! Now get the surface albedo variables
       if (socrates_opt%l_variable_srf_albedo) then
          call log_master_log &
            (LOG_ERROR, "Socrates config using variable surface albedo, but this has not been developed. Set to false - STOP")
       endif
       ! Now get the longitude and latitude variables
       if (socrates_opt%latitude < -90.0 .or. socrates_opt%latitude > 90.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - latitude is outside sensible range, check config - STOP")
       endif
       if (socrates_opt%longitude < -180.0 .or. socrates_opt%longitude > 180.0) then
          call log_master_log &
               (LOG_ERROR, "Socrates - longitude is outside sensible range, check config - STOP")
       endif
    endif ! end l_solar_fixed

    if (socrates_opt%rad_interval <= 0) then
       call log_master_log &
            (LOG_WARN, "Socrates - rad_interval <= 0 ; SOCRATES will be called every timestep")
    endif
    !if (options_has_key(current_state%options_database, "rad_int_time")) then
    !   call log_master_log &
    !        (LOG_ERROR, "Socrates - option key 'rad_int_time' is deprecated and no longer functions.  "//&
    !                    "Please remove this from your configuration, and use 'rad_interval', which "//&
    !                    "has functionality dependent upon 'time_basis'.")
    !end if

    if (socrates_opt%surface_albedo < 0.0 .or. socrates_opt%surface_albedo > 1.0) then
       call log_master_log &
            (LOG_ERROR, "Socrates - surface albedo outside sensible range, check config - STOP")
    endif

     ! set the well mixed gases. Values are based on UM GA settings
     ! This should be moved to configuration once reading well mixed
     ! gases from configuration can be shown to bit compare (see #306)
     !
    socrates_opt%co2_mmr = 5.94100e-4
    socrates_opt%n2o_mmr = 4.925e-7
    socrates_opt%ch4_mmr = 9.994e-7
    socrates_opt%o2_mmr = 0.2314
    socrates_opt%cfc12_mmr = 1.129e-9
    socrates_opt%cfc11_mmr = 2.225e-9
    socrates_opt%cfc113_mmr = 0.0
    socrates_opt%cfc114_mmr = 0.0
    socrates_opt%hcfc22_mmr = 0.0
    socrates_opt%hfc125_mmr = 0.0
    socrates_opt%hfc134a_mmr = 0.0

     ! Change the well-mixed gas concentrations in the RCEMIP case
     ! As above, see #306
    if (rcemip_gases_enabled .eqv. .true.) then
      socrates_opt%co2_mmr = 5.288e-4
      socrates_opt%n2o_mmr = 4.651e-7
      socrates_opt%ch4_mmr = 9.139e-7
      socrates_opt%o2_mmr = 0.2314
      socrates_opt%cfc12_mmr = 0.0
      socrates_opt%cfc11_mmr = 0.0
      socrates_opt%cfc113_mmr = 0.0
      socrates_opt%cfc114_mmr = 0.0
      socrates_opt%hcfc22_mmr = 0.0
      socrates_opt%hfc125_mmr = 0.0
      socrates_opt%hfc134a_mmr = 0.0
    end if

  end subroutine set_and_test_socrates_monc_options


  end module get_and_test_socrates_options_mod

