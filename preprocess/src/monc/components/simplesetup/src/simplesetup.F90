










module simplesetup_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use datadefn_mod, only : DEFAULT_PRECISION, LONG_STRING_LENGTH, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use state_mod, only : PRESCRIBED_SURFACE_FLUXES, model_state_type
  use conversions_mod, only : conv_to_string
  use logging_mod, only :  LOG_INFO, LOG_ERROR, log_log, log_master_log
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX, PRIMAL_GRID, DUAL_GRID
  use prognostics_mod, only : prognostic_field_type
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
       options_get_integer_array, options_get_real_array
  use tracers_mod, only : get_tracer_options
  use q_indices_mod, only: get_q_index, standard_q_names

  implicit none

  private

  integer :: x_size, y_size, z_size
  real(kind=DEFAULT_PRECISION) :: zztop, dxx, dyy
  logical :: enable_theta=.false.

  public initialisation_callback_simplesetup
contains


  subroutine initialisation_callback_simplesetup(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer :: iter
    real(kind=DEFAULT_PRECISION) :: realnum

    call read_configuration(current_state)
    if (.not. current_state%initialised) then
      do iter = 1,current_state%config_args
        if (current_state%options_database_string(iter,1) .eq. "dtm") then
          read(current_state%options_database_string(iter,2),*) realnum
          current_state%dtm = realnum
        end if
      end do
    !  current_state%dtm=options_get_real(current_state%options_database, "dtm")
      current_state%dtm_new=current_state%dtm
      call create_grid(current_state, current_state%global_grid)
      call decompose_grid(current_state)
      call allocate_prognostics(current_state)
      current_state%initialised=.true.
    end if

    ! Isolate MONC process that has the requested print_debug_data location by changing the value of the flag to be true
    !   only on that process.
    if (current_state%print_debug_data) then
      current_state%print_debug_data =                                        &
        current_state%pdd_x .ge. current_state%local_grid%start(X_INDEX).and. &
        current_state%pdd_x .le. current_state%local_grid%end(X_INDEX)  .and. &
        current_state%pdd_y .ge. current_state%local_grid%start(Y_INDEX).and. &
        current_state%pdd_y .le. current_state%local_grid%end(Y_INDEX)
    end if
  end subroutine initialisation_callback_simplesetup

  
  subroutine allocate_prognostics(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: alloc_z, alloc_y, alloc_x, i

    alloc_z=current_state%local_grid%size(Z_INDEX) + current_state%local_grid%halo_size(Z_INDEX) * 2
    alloc_y=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    alloc_x=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    
    call allocate_prognostic(current_state%u, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, PRIMAL_GRID)
    call allocate_prognostic(current_state%zu, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, PRIMAL_GRID)
    call allocate_prognostic(current_state%su, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, PRIMAL_GRID)
    call allocate_prognostic(current_state%savu, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, PRIMAL_GRID)
    call allocate_prognostic(current_state%v, alloc_z, alloc_y, alloc_x, DUAL_GRID, PRIMAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zv, alloc_z, alloc_y, alloc_x, DUAL_GRID, PRIMAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sv, alloc_z, alloc_y, alloc_x, DUAL_GRID, PRIMAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%savv, alloc_z, alloc_y, alloc_x, DUAL_GRID, PRIMAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%w, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)    
    call allocate_prognostic(current_state%zw, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)        
    call allocate_prognostic(current_state%sw, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)        
    call allocate_prognostic(current_state%savw, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    if (enable_theta) then
      call allocate_prognostic(current_state%th, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
      call allocate_prognostic(current_state%zth, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
      call allocate_prognostic(current_state%sth, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    end if

    ! q species
    call allocate_prognostic(current_state%qv, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqv, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqv, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%ql, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zql, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sql, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%qr, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqr, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqr, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%qi, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqi, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqi, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%qs, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqs, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqs, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%qg, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqg, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqg, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%qAitkenSolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqAitkenSolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqAitkenSolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%qAccumSolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqAccumSolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqAccumSolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%qAccumInsolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqAccumInsolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqAccumInsolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%qCoarseSolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqCoarseSolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqCoarseSolMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%qCoarseDustMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zqCoarseDustMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sqCoarseDustMass, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%nl, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%znl, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%snl, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%nr, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%znr, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%snr, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%ni, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zni, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sni, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%ns, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zns, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sns, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%ng, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%zng, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sng, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%nAitkenSolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%znAitkenSolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%snAitkenSolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%nAccumSolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%znAccumSolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%snAccumSolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%nAccumInsolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%znAccumInsolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%snAccumInsolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%nCoarseSolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%znCoarseSolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%snCoarseSolNumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%nCoarseDustnumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%znCoarseDustnumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%snCoarseDustnumber, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%rdAitkenSol, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%rdAccumSol, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%rdCoarseSol, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%D0_cloud, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%D0_rain, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%D0_ice, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%D0_snow, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%D0_graupel, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%Tk, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%RH, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%RI, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%qv_saturation, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%qi_saturation, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)

    call allocate_prognostic(current_state%cloud_reff, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%lwrad_hr, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%swrad_hr, alloc_z, alloc_y, alloc_x, PRIMAL_GRID, DUAL_GRID, DUAL_GRID)


      !! Set standard q indices
      current_state%water_vapour_mixing_ratio_index=1!get_q_index(standard_q_names%VAPOUR, 'simplesetup')
      current_state%liquid_water_mixing_ratio_index=2!get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'simplesetup')
    !end if0
    
    if (current_state%n_tracers .gt. 0) then
      allocate( current_state%tracer(current_state%n_tracers),  &
                current_state%ztracer(current_state%n_tracers), &
                current_state%stracer(current_state%n_tracers))
      do i=1, current_state%n_tracers
        call allocate_prognostic(current_state%tracer(i), alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
        call allocate_prognostic(current_state%ztracer(i), alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
        call allocate_prognostic(current_state%stracer(i), alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
      end do
      
    endif ! allocate tracers

    ! Set arrays for radiative heating rates - Note: this should be protected by a switch
    call allocate_prognostic(current_state%sth_lw, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID)
    call allocate_prognostic(current_state%sth_sw, alloc_z, alloc_y, alloc_x, DUAL_GRID, DUAL_GRID, DUAL_GRID) 

  end subroutine allocate_prognostics

  subroutine allocate_prognostic(field, alloc_z, alloc_y, alloc_x, z_grid, y_grid, x_grid)
    type(prognostic_field_type), intent(inout) :: field
    integer, intent(in) :: alloc_z, alloc_y, alloc_x, z_grid, y_grid, x_grid

    field%active=.true.
    field%grid(Z_INDEX) = z_grid
    field%grid(Y_INDEX) = y_grid
    field%grid(X_INDEX) = x_grid
    allocate(field%data(alloc_z, alloc_y, alloc_x))
    field%data=0.0_DEFAULT_PRECISION
  end subroutine allocate_prognostic  

  subroutine decompose_grid(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (associated(current_state%parallel%decomposition_procedure)) then
      call current_state%parallel%decomposition_procedure(current_state)
    else
      call log_log(LOG_ERROR, "No decomposition specified")
    end if
  end subroutine decompose_grid

  subroutine create_grid(current_state, specific_grid) !Lambert
    type(model_state_type), intent(inout) :: current_state
    type(global_grid_type), intent(inout) :: specific_grid

    integer :: iter1, iter2, size_array
    real(kind=DEFAULT_PRECISION) :: realnum

    do iter1 = 1,current_state%config_args
      size_array = 0
      if (current_state%options_database_string(iter1,1) .eq. "origional_vertical_grid_interp_levels") then
        if (isnan(current_state%options_database_real(iter1,1))) then
          allocate(current_state%global_grid%configuration%vertical%origional_vertical_grid_interp_levels(1))
          read(current_state%options_database_string(iter1,2),*) realnum
          current_state%global_grid%configuration%vertical%origional_vertical_grid_interp_levels(1) = realnum
        else
          do iter2 = 1,size(current_state%options_database_real(iter1,:))
            if (isnan(current_state%options_database_real(iter1,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          allocate(current_state%global_grid%configuration%vertical%origional_vertical_grid_interp_levels(size_array))
          current_state%global_grid%configuration%vertical%origional_vertical_grid_interp_levels = &
              current_state%options_database_real(iter1,1:size_array)
        end if
      else if (current_state%options_database_string(iter1,1) .eq. "origional_vertical_grid_interp_altitudes") then
        if (isnan(current_state%options_database_real(iter1,1))) then
          allocate(current_state%global_grid%configuration%vertical%origional_vertical_grid_interp_altitudes(1))
          read(current_state%options_database_string(iter1,2),*) realnum
          current_state%global_grid%configuration%vertical%origional_vertical_grid_interp_altitudes = realnum
        else
          do iter2 = 1,size(current_state%options_database_real(iter1,:))
            if (isnan(current_state%options_database_real(iter1,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          allocate(current_state%global_grid%configuration%vertical%origional_vertical_grid_interp_altitudes(size_array))
          current_state%global_grid%configuration%vertical%origional_vertical_grid_interp_altitudes = &
              current_state%options_database_real(iter1,1:size_array)
        end if
      end if
    end do

    specific_grid%bottom(Z_INDEX) = 0
    specific_grid%bottom(Y_INDEX) = 0
    specific_grid%bottom(X_INDEX) = 0

    specific_grid%top(Z_INDEX) = zztop
    specific_grid%top(Y_INDEX) = dyy * y_size
    specific_grid%top(X_INDEX) = dxx * x_size

    !specific_grid%resolution(Z_INDEX) = zztop / z_size
    specific_grid%resolution(Y_INDEX) = dyy
    specific_grid%resolution(X_INDEX) = dxx

    specific_grid%size(Z_INDEX) = z_size
    specific_grid%size(Y_INDEX) = y_size
    specific_grid%size(X_INDEX) = x_size

    specific_grid%active(Z_INDEX) = .true.
    specific_grid%active(Y_INDEX) = .true.
    specific_grid%active(X_INDEX) = .true.

    specific_grid%dimensions = 3
  end subroutine create_grid

  subroutine read_configuration(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer ::  intnum, iter
    real(kind=DEFAULT_PRECISION) :: realnum
    logical :: logicnum

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "rhobous") then
        read(current_state%options_database_string(iter,2),*) realnum
        current_state%rhobous = realnum
      else if (current_state%options_database_string(iter,1) .eq. "thref0") then
        read(current_state%options_database_string(iter,2),*) realnum
        current_state%thref0 = realnum
      else if (current_state%options_database_string(iter,1) .eq. "number_q_fields") then
        read(current_state%options_database_string(iter,2),*) intnum
        current_state%number_q_fields = intnum
      else if (current_state%options_database_string(iter,1) .eq. "surface_pressure") then
        read(current_state%options_database_string(iter,2),*) realnum
        current_state%surface_pressure = realnum
      else if (current_state%options_database_string(iter,1) .eq. "surface_reference_pressure") then
        read(current_state%options_database_string(iter,2),*) realnum
        current_state%surface_reference_pressure = realnum
      else if (current_state%options_database_string(iter,1) .eq. "use_anelastic_equations") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%use_anelastic_equations = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "origional_vertical_grid_setup") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%origional_vertical_grid_setup = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "new_vertical_grid_setup") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%new_vertical_grid_setup = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "special_vertical_grid_setup") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%special_vertical_grid_setup = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "passive_q") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%passive_q = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "passive_th") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%passive_th = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "rmlmax") then
        read(current_state%options_database_string(iter,2),*) realnum
        current_state%rmlmax = realnum
      else if (current_state%options_database_string(iter,1) .eq. "calculate_th_and_q_init") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%calculate_th_and_q_init = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "use_viscosity_and_diffusion") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%use_viscosity_and_diffusion = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "backscatter") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%backscatter = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "x_size") then
        read(current_state%options_database_string(iter,2),*) intnum
        x_size = intnum
      else if (current_state%options_database_string(iter,1) .eq. "y_size") then
        read(current_state%options_database_string(iter,2),*) intnum
        y_size = intnum
      else if (current_state%options_database_string(iter,1) .eq. "z_size") then
        read(current_state%options_database_string(iter,2),*) intnum
        z_size = intnum
        current_state%global_grid%configuration%vertical%z_size = z_size
      else if (current_state%options_database_string(iter,1) .eq. "dxx") then
        read(current_state%options_database_string(iter,2),*) realnum
        dxx = realnum
      else if (current_state%options_database_string(iter,1) .eq. "dyy") then
        read(current_state%options_database_string(iter,2),*) realnum
        dyy = realnum
      else if (current_state%options_database_string(iter,1) .eq. "zztop") then
        read(current_state%options_database_string(iter,2),*) realnum
        zztop = realnum
      else if (current_state%options_database_string(iter,1) .eq. "enable_theta") then
        read(current_state%options_database_string(iter,2),*) logicnum
        enable_theta = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "galilean_transformation") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%galilean_transformation = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "print_debug_data") then
        read(current_state%options_database_string(iter,2),*) logicnum
        current_state%print_debug_data = logicnum
     end if
    end do

    if (current_state%rmlmax<=0.0)current_state%rmlmax=0.23 * max(dxx, dyy)

    if (.not. enable_theta) current_state%passive_th=.true.
    if ( current_state%number_q_fields == 0) current_state%passive_q=.true.


    if (current_state%galilean_transformation)then
      do iter = 1,current_state%config_args
        if (current_state%options_database_string(iter,1) .eq. "fix_ugal") then
          read(current_state%options_database_string(iter,2),*) logicnum
          current_state%fix_ugal = logicnum
        else if (current_state%options_database_string(iter,1) .eq. "fix_vgal") then
          read(current_state%options_database_string(iter,2),*) logicnum
          current_state%fix_vgal = logicnum
        end if
      end do

      do iter = 1,current_state%config_args
        if (current_state%fix_ugal .eqv. .true.) then
          if (current_state%options_database_string(iter,1) .eq. "ugal") then
            read(current_state%options_database_string(iter,2),*) realnum
            current_state%ugal = realnum
          else if (current_state%options_database_string(iter,1) .eq. "vgal") then
            read(current_state%options_database_string(iter,2),*) realnum
            current_state%vgal = realnum
          end if
        end if
      end do
    end if

    ! Parameters for print_debug_data

    if (current_state%print_debug_data) then
      do iter = 1,current_state%config_args
        if (current_state%options_database_string(iter,1) .eq. "pdd_z") then
          read(current_state%options_database_string(iter,2),*) intnum
          current_state%pdd_z = intnum
          if (current_state%pdd_z .lt. 0) current_state%pdd_z = z_size/2
        else if (current_state%options_database_string(iter,1) .eq. "pdd_y") then
          read(current_state%options_database_string(iter,2),*) intnum
          current_state%pdd_y = intnum
          if (current_state%pdd_y .lt. 0) current_state%pdd_y = y_size/2
        else if (current_state%options_database_string(iter,1) .eq. "pdd_x") then
          read(current_state%options_database_string(iter,2),*) intnum
          current_state%pdd_x = intnum
          if (current_state%pdd_x .lt. 0) current_state%pdd_x = x_size/2
        end if
      end do

      current_state%column_global_x = current_state%pdd_x
      current_state%column_global_y = current_state%pdd_y
      current_state%halo_column = .false.
    end if

    !! a voir pour le restart
    !if (.not. current_state%reconfig_run) then
    !  call get_tracer_options(current_state)
    !end if ! not reconfig

  end subroutine read_configuration
end module simplesetup_mod
