!> Derive 3-D diagnostic fields which are available in current_state 
module diagnostics_3d_mod
  !!
  !! module to derive 3-D fields which are not already stored
  !! in current_state and output as a 3-D diagnostic
  !!
  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, &
       component_descriptor_type, component_field_value_type, component_field_information_type,&
       component_descriptor_type_v1
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_real
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION
  use q_indices_mod, only: get_q_index, standard_q_names
  use saturation_mod, only: qsaturation
  use science_constants_mod, only : rlvap_over_cp

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: iqv, iql, iqr
  !real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
  !     TdegK,               & ! absolute temperature in kelvin
  !     theta,               & ! potential temperature in kelvin (th + thref)
  !     liquid_ice_theta       ! liquid-ice potential temperature in kelvin
  real, parameter :: p0 = 101300 ! Pa pression at sel level
  real, parameter :: T0 = 373.0 ! K ebulition temperature of water at one atmosphere (p0)
  real, parameter :: Lv = 2.26E6 ! J.kg-1 water vaporization laten heat
  real, parameter :: R = 8.314 ! J.K-1.mol-1 water vaporization laten heat
  real, parameter :: M = 0.018 ! kg.mol-1 water vaporization laten heat
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &    
       total_condensate
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: e ! Pa water vapor pressure
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: e_saturation ! Pa saturation water vapor pressure
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: Tkelv ! Temperature (K)
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: qv_saturation ! saturated water vapour mixing ratio (K)
  real(kind=DEFAULT_PRECISION) :: qlcrit

  public initialisation_callback_diagnostics_3d, timestep_callback_diagnostics_3d

contains
  
  subroutine initialisation_callback_diagnostics_3d(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (current_state%th%active) then
       allocate(Tkelv(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),            &
            current_state%local_grid%size(X_INDEX)))
       allocate(qv_saturation(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),            &
            current_state%local_grid%size(X_INDEX)))
       allocate(current_state%TdegK(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),            &
            current_state%local_grid%size(X_INDEX)))
       allocate(current_state%theta(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),            &
            current_state%local_grid%size(X_INDEX)))
       allocate(current_state%liquid_ice_theta(current_state%local_grid%size(Z_INDEX),  &
            current_state%local_grid%size(Y_INDEX),            &
            current_state%local_grid%size(X_INDEX)))
    endif

    if (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) then
       allocate(total_condensate(current_state%local_grid%size(Z_INDEX)))
       iqv= 1 !get_q_index(standard_q_names%VAPOUR, 'diagnostics_3d')
       iql= 2 !get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'diagnostics_3d')
       if (current_state%rain_water_mixing_ratio_index > 0) & 
            iqr = current_state%rain_water_mixing_ratio_index
    endif
    
  end subroutine initialisation_callback_diagnostics_3d
  
  subroutine timestep_callback_diagnostics_3d(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    real(kind=DEFAULT_PRECISION) :: exner, pmb, qv, qs, T_K

    integer :: k
    integer :: current_y_index, current_x_index, target_x_index, target_y_index

    if (current_state%halo_column) return
    !if (.not. current_state%diagnostic_sample_timestep) return
    if (current_state%modulo_number_3d .ne. 0) return
   
    current_y_index=current_state%column_local_y
    current_x_index=current_state%column_local_x
    target_y_index=current_y_index-current_state%local_grid%halo_size(Y_INDEX)
    target_x_index=current_x_index-current_state%local_grid%halo_size(X_INDEX)

    !

    if (current_state%th%active) then
       current_state%theta(:,target_y_index, target_x_index) = &
            (current_state%th%data(:,current_y_index,current_x_index)       &
            + current_state%global_grid%configuration%vertical%thref(:))
       current_state%liquid_ice_theta(:,target_y_index, target_x_index) =   &
            (current_state%th%data(:,current_y_index,current_x_index)       &
            + current_state%global_grid%configuration%vertical%thref(:))  
            ! test for the qfields
       if (.not. current_state%passive_q .and. &
            current_state%number_q_fields .gt. 0) then
          Tkelv(:,target_y_index, target_x_index) = current_state%theta(:,target_y_index, target_x_index) * &
                                              current_state%global_grid%configuration%vertical%rprefrcp(:)
          !do k = 1,current_state%local_grid%size(Z_INDEX)
          !  qv_saturation(k,target_y_index, target_x_index) = qsaturation(Tkelv(k,target_y_index, target_x_index), &
          !                    current_state%global_grid%configuration%vertical%prefn(k)/100.0)
          !end do

          total_condensate(:) =                                               &
               current_state%ql%data(:,current_y_index,current_x_index) + &
               current_state%qr%data(:,current_y_index,current_x_index)
          current_state%liquid_ice_theta(:,target_y_index, target_x_index) =     &
               current_state%liquid_ice_theta(:,target_y_index, target_x_index) - &
               ( total_condensate(:) * (rlvap_over_cp) )
       endif
       current_state%TdegK(:,target_y_index, target_x_index) =              &
            (current_state%th%data(:,current_y_index,current_x_index)       &
            + current_state%global_grid%configuration%vertical%thref(:))    &
            * current_state%global_grid%configuration%vertical%rprefrcp(:)       
    endif

  end subroutine timestep_callback_diagnostics_3d

end module diagnostics_3d_mod
