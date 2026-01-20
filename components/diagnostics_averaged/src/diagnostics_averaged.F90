!> Average dignostics
module diagnostics_averaged_mod
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
  use science_constants_mod, only : rlvap_over_cp, z0, z0th

  implicit none

#ifndef TEST_MODE
  private
#endif
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: mo_l_to_ope, friction_vel_to_ope, z0_to_ope, z0th_to_ope
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: th_to_ope, u_to_ope, v_to_ope, &
                                qv_to_ope, ql_to_ope, qr_to_ope, qi_to_ope, qs_to_ope, qg_to_ope, &
                                qAitkenSolMass_to_ope, qAccumSolMass_to_ope, qAccumInsolMass_to_ope, &
                                qCoarseSolMass_to_ope, qCoarseDustMass_to_ope, &
                                nl_to_ope, nr_to_ope, ni_to_ope, ns_to_ope, ng_to_ope, &
                                nAitkenSolNumber_to_ope, nAccumSolNumber_to_ope, nAccumInsolNumber_to_ope, &
                                nCoarseSolNumber_to_ope, nCoarseDustnumber_to_ope

  public initialisation_callback_diagnostics_averaged, timestep_callback_diagnostics_averaged

contains

  subroutine initialisation_callback_diagnostics_averaged(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    allocate(th_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(u_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(v_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qv_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(ql_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qr_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qi_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qs_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qg_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qAitkenSolMass_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qAccumSolMass_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qAccumInsolMass_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qCoarseSolMass_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(qCoarseDustMass_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(nl_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(nr_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(ni_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(ns_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(ng_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(nAitkenSolNumber_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(nAccumSolNumber_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(nAccumInsolNumber_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(nCoarseSolNumber_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(nCoarseDustnumber_to_ope(current_state%dim_var_to_mean + 1, current_state%local_grid%size(Z_INDEX)))
    allocate(mo_l_to_ope(current_state%dim_var_to_mean + 1))
    allocate(friction_vel_to_ope(current_state%dim_var_to_mean + 1))
    allocate(z0_to_ope(current_state%dim_var_to_mean + 1))
    allocate(z0th_to_ope(current_state%dim_var_to_mean + 1))

    allocate(current_state%mean_th(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_u(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_v(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qv(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_ql(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qr(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qi(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qs(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qg(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qAitkenSolMass(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qAccumSolMass(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qAccumInsolMass(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qCoarseSolMass(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_qCoarseDustMass(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_nl(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_nr(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_ni(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_ns(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_ng(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_nAitkenSolNumber(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_nAccumSolNumber(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_nAccumInsolNumber(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_nCoarseSolNumber(current_state%local_grid%size(Z_INDEX)))
    allocate(current_state%mean_nCoarseDustnumber(current_state%local_grid%size(Z_INDEX)))

    th_to_ope = 0.0
    u_to_ope = 0.0
    v_to_ope = 0.0
    qv_to_ope = 0.0
    ql_to_ope = 0.0
    qr_to_ope = 0.0
    qi_to_ope = 0.0
    qs_to_ope = 0.0
    qg_to_ope = 0.0
    qAitkenSolMass_to_ope = 0.0
    qAccumSolMass_to_ope = 0.0
    qAccumInsolMass_to_ope = 0.0
    qCoarseSolMass_to_ope = 0.0
    qCoarseDustMass_to_ope = 0.0
    nl_to_ope = 0.0
    nr_to_ope = 0.0
    ni_to_ope = 0.0
    ns_to_ope = 0.0
    ng_to_ope = 0.0
    nAitkenSolNumber_to_ope = 0.0
    nAccumSolNumber_to_ope = 0.0
    nAccumInsolNumber_to_ope = 0.0
    nCoarseSolNumber_to_ope = 0.0
    nCoarseDustnumber_to_ope = 0.0
    mo_l_to_ope = 0.0
    friction_vel_to_ope = 0.0
    z0_to_ope = 0.0
    z0th_to_ope =0.0


    current_state%mean_th = 0.0
    current_state%mean_u = 0.0
    current_state%mean_v = 0.0
    current_state%mean_qv = 0.0
    current_state%mean_ql = 0.0
    current_state%mean_qr = 0.0
    current_state%mean_qi = 0.0
    current_state%mean_qs = 0.0
    current_state%mean_qg = 0.0
    current_state%mean_qAitkenSolMass = 0.0
    current_state%mean_qAccumSolMass = 0.0
    current_state%mean_qAccumInsolMass = 0.0
    current_state%mean_qCoarseSolMass = 0.0
    current_state%mean_qCoarseDustMass = 0.0
    current_state%mean_nl = 0.0
    current_state%mean_nr = 0.0
    current_state%mean_ni = 0.0
    current_state%mean_ns = 0.0
    current_state%mean_ng = 0.0
    current_state%mean_nAitkenSolNumber = 0.0
    current_state%mean_nAccumSolNumber = 0.0
    current_state%mean_nAccumInsolNumber = 0.0
    current_state%mean_nCoarseSolNumber = 0.0
    current_state%mean_nCoarseDustnumber = 0.0
    current_state%mean_mo_l = 0.0
    current_state%mean_friction_vel = 0.0
    current_state%mean_z0 = 0.0
    current_state%mean_z0th = 0.0

  end subroutine initialisation_callback_diagnostics_averaged

  subroutine timestep_callback_diagnostics_averaged(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer ::  iter
    real(kind=DEFAULT_PRECISION) :: sum_mo_l_mean, sum_friction_vel_mean, sum_z0_mean, sum_z0th_mean
    real(kind=DEFAULT_PRECISION), dimension(current_state%local_grid%size(Z_INDEX))  :: sum_th_mean, sum_u_mean, sum_v_mean, &
                                  sum_qv_mean, sum_ql_mean, sum_qr_mean, sum_qi_mean, sum_qs_mean, sum_qg_mean, &
                                  sum_qAitkenSolMass_mean, sum_qAccumSolMass_mean, sum_qAccumInsolMass_mean, &
                                  sum_qCoarseSolMass_mean, sum_qCoarseDustMass_mean, &
                                  sum_nl_mean, sum_nr_mean, sum_ni_mean, sum_ns_mean, sum_ng_mean, &
                                  sum_nAitkenSolNumber_mean, sum_nAccumSolNumber_mean, sum_nAccumInsolNumber_mean, &
                                  sum_nCoarseSolNumber_mean, sum_nCoarseDustnumber_mean

    th_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olthbar(:)
    u_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olubar(:)
    v_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olvbar(:)
    qv_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olqvbar(:)
    ql_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olqlbar(:)
    qr_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olqrbar(:)
    qi_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olqibar(:)
    qs_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olqsbar(:)
    qg_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olqgbar(:)
    qAitkenSolMass_to_ope(current_state%average_sample_iter, :) = &
                          current_state%global_grid%configuration%vertical%olqAitkenSolMassbar(:)
    qAccumSolMass_to_ope(current_state%average_sample_iter, :) = &
                          current_state%global_grid%configuration%vertical%olqAccumSolMassbar(:)
    qAccumInsolMass_to_ope(current_state%average_sample_iter, :) = &
                          current_state%global_grid%configuration%vertical%olqAccumInsolMassbar(:)
    qCoarseSolMass_to_ope(current_state%average_sample_iter, :) = &
                          current_state%global_grid%configuration%vertical%olqCoarseSolMassbar(:)
    qCoarseDustMass_to_ope(current_state%average_sample_iter, :) = &
                          current_state%global_grid%configuration%vertical%olqCoarseDustMassbar(:)
    nl_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olnlbar(:)
    nr_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olnrbar(:)
    ni_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olnibar(:)
    ns_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olnsbar(:)
    ng_to_ope(current_state%average_sample_iter, :) = current_state%global_grid%configuration%vertical%olngbar(:)
    nAitkenSolNumber_to_ope(current_state%average_sample_iter, :) = &
                            current_state%global_grid%configuration%vertical%olnAitkenSolNumberbar(:)
    nAccumSolNumber_to_ope(current_state%average_sample_iter, :) = &
                            current_state%global_grid%configuration%vertical%olznAccumSolNumberbar(:)
    nAccumInsolNumber_to_ope(current_state%average_sample_iter, :) = &
                            current_state%global_grid%configuration%vertical%olznAccumInsolNumberbar(:)
    nCoarseSolNumber_to_ope(current_state%average_sample_iter, :) = &
                            current_state%global_grid%configuration%vertical%olznCoarseSolNumberbar(:)
    nCoarseDustnumber_to_ope(current_state%average_sample_iter, :) = &
                            current_state%global_grid%configuration%vertical%olznCoarseDustnumberbar(:)
    mo_l_to_ope(current_state%average_sample_iter) = current_state%mo_l
    friction_vel_to_ope(current_state%average_sample_iter) = current_state%friction_vel
    z0_to_ope(current_state%average_sample_iter) = z0
    z0th_to_ope(current_state%average_sample_iter) = z0th

    current_state%average_sample_iter = current_state%average_sample_iter + 1

    if ((current_state%dim_var_to_mean + 1) .eq. current_state%average_sample_iter) then
      do iter = 1, current_state%dim_var_to_mean
        sum_th_mean(:) = sum_th_mean(:) + th_to_ope(iter,:)
        sum_u_mean(:) = sum_u_mean(:) + u_to_ope(iter,:)
        sum_v_mean(:) = sum_v_mean(:) + v_to_ope(iter,:)
        sum_qv_mean(:) = sum_qv_mean(:) + qv_to_ope(iter,:)
        sum_ql_mean(:) = sum_ql_mean(:) + ql_to_ope(iter,:)
        sum_qr_mean(:) = sum_qr_mean(:) + qr_to_ope(iter,:)
        sum_qi_mean(:) = sum_qi_mean(:) + qi_to_ope(iter,:)
        sum_qs_mean(:) = sum_qs_mean(:) + qs_to_ope(iter,:)
        sum_qg_mean(:) = sum_qg_mean(:) + qg_to_ope(iter,:)
        sum_qAitkenSolMass_mean(:) = sum_qAitkenSolMass_mean(:) + qAitkenSolMass_to_ope(iter,:)
        sum_qAccumSolMass_mean(:) = sum_qAccumSolMass_mean(:) + qAccumSolMass_to_ope(iter,:)
        sum_qAccumInsolMass_mean(:) = sum_qAccumInsolMass_mean(:) + qAccumInsolMass_to_ope(iter,:)
        sum_qCoarseSolMass_mean(:) = sum_qCoarseSolMass_mean(:) + qCoarseSolMass_to_ope(iter,:)
        sum_qCoarseDustMass_mean(:) = sum_qCoarseDustMass_mean(:) + qCoarseDustMass_to_ope(iter,:)
        sum_nl_mean(:) = sum_nl_mean(:) + nl_to_ope(iter,:)
        sum_nr_mean(:) = sum_nr_mean(:) + nr_to_ope(iter,:)
        sum_ni_mean(:) = sum_ni_mean(:) + ni_to_ope(iter,:)
        sum_ns_mean(:) = sum_ns_mean(:) + ns_to_ope(iter,:)
        sum_ng_mean(:) = sum_ng_mean(:) + ng_to_ope(iter,:)
        sum_nAitkenSolNumber_mean(:) = sum_nAitkenSolNumber_mean(:) + nAitkenSolNumber_to_ope(iter,:)
        sum_nAccumSolNumber_mean(:) = sum_nAccumSolNumber_mean(:) + nAccumSolNumber_to_ope(iter,:)
        sum_nAccumInsolNumber_mean(:) = sum_nAccumInsolNumber_mean(:) + nAccumInsolNumber_to_ope(iter,:)
        sum_nCoarseSolNumber_mean(:) = sum_nCoarseSolNumber_mean(:) + nCoarseSolNumber_to_ope(iter,:)
        sum_nCoarseDustnumber_mean(:) = sum_nCoarseDustnumber_mean(:) + nCoarseDustnumber_to_ope(iter,:)
        sum_mo_l_mean = sum_mo_l_mean + mo_l_to_ope(iter)
        sum_friction_vel_mean = sum_friction_vel_mean + friction_vel_to_ope(iter)
        sum_z0_mean = sum_z0_mean + z0_to_ope(iter)
        sum_z0th_mean = sum_z0th_mean + z0th_to_ope(iter)
      end do

      current_state%mean_th(:) = sum_th_mean(:) / current_state%dim_var_to_mean
      current_state%mean_u(:) = sum_u_mean(:) / current_state%dim_var_to_mean
      current_state%mean_v(:) = sum_v_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qv(:) = sum_qv_mean(:) / current_state%dim_var_to_mean
      current_state%mean_ql(:) = sum_ql_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qr(:) = sum_qr_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qi(:) = sum_qi_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qs(:) = sum_qs_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qg(:) = sum_qg_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qAitkenSolMass(:) = sum_qAitkenSolMass_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qAccumSolMass(:) = sum_qAccumSolMass_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qAccumInsolMass(:) = sum_qAccumInsolMass_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qCoarseSolMass(:) = sum_qCoarseSolMass_mean(:) / current_state%dim_var_to_mean
      current_state%mean_qCoarseDustMass(:) = sum_qCoarseDustMass_mean(:) / current_state%dim_var_to_mean
      current_state%mean_nl(:) = sum_nl_mean(:) / current_state%dim_var_to_mean
      current_state%mean_nr(:) = sum_nr_mean(:) / current_state%dim_var_to_mean
      current_state%mean_ni(:) = sum_ni_mean(:) / current_state%dim_var_to_mean
      current_state%mean_ns(:) = sum_ns_mean(:) / current_state%dim_var_to_mean
      current_state%mean_ng(:) = sum_ng_mean(:) / current_state%dim_var_to_mean
      current_state%mean_nAitkenSolNumber(:) = sum_nAitkenSolNumber_mean(:) / current_state%dim_var_to_mean
      current_state%mean_nAccumSolNumber(:) = sum_nAccumSolNumber_mean(:) / current_state%dim_var_to_mean
      current_state%mean_nAccumInsolNumber(:) = sum_nAccumInsolNumber_mean(:) / current_state%dim_var_to_mean
      current_state%mean_nCoarseSolNumber(:) = sum_nCoarseSolNumber_mean(:) / current_state%dim_var_to_mean
      current_state%mean_nCoarseDustnumber(:) = sum_nCoarseDustnumber_mean(:) / current_state%dim_var_to_mean
      current_state%mean_mo_l = sum_mo_l_mean / current_state%dim_var_to_mean
      current_state%mean_friction_vel = sum_friction_vel_mean / current_state%dim_var_to_mean
      current_state%mean_z0 = sum_z0_mean / current_state%dim_var_to_mean
      current_state%mean_z0th = sum_z0th_mean / current_state%dim_var_to_mean

      current_state%average_sample_iter = 1
    end if

  end subroutine timestep_callback_diagnostics_averaged

end module diagnostics_averaged_mod
