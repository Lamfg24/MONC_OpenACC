!> The main IO server Netcdf writer

module io_server_netcdf_writer_mod
    use netcdf
    use state_mod, only : model_state_type
    use conversions_mod, only : conv_to_string
    use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH, LONG_STRING_LENGTH
    use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX, local_grid_type, vertical_grid_configuration_type
    use mpi, only : MPI_COMM_WORLD, MPI_STATUSES_IGNORE, MPI_BYTE, MPI_INT, MPI_STATUS_IGNORE, MPI_REAL, &
                    MPI_INFO_NULL, MPI_DOUBLE, MPI_LOGICAL, MPI_CHAR

#ifndef TEST_MODE
  private
#endif

integer :: time_id
integer :: global_grid_z_size_netcdf, global_grid_y_size_netcdf, global_grid_x_size_netcdf, size_bin, timeline

public diagnostic_file_0d_generation, diagnostic_file_1d_generation, diagnostic_file_2d_generation, &
       diagnostic_file_3d_generation, checkpoint_file_generation, diagnostic_file_average_generation

contains

  subroutine diagnostic_file_0d_generation(current_state, vertical_grid, io_communicator_arg, diagnostic_path)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: io_communicator_arg
    character(len=200), intent(in) :: diagnostic_path

    !integer :: scalar_size_id, time_id
    integer :: ncdf_id, ierr
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename
    real(kind=DEFAULT_PRECISION) :: gwp_min, gwp_mean, gwp_max, iwp_min, iwp_mean, iwp_max, lathf_min, lathf_mean, lathf_max, &
                                    lwp_min, lwp_mean, lwp_max, reske_min, reske_mean, reske_max, rwp_min, rwp_mean, rwp_max, &
                                    senhf_min, senhf_mean, senhf_max, subke_2d_min, subke_2d_mean, subke_2d_max, &
                                    surface_precip_min, surface_precip_mean, surface_precip_max, swp_min, swp_mean, swp_max, &
                                    tot_iwp_min, tot_iwp_mean, tot_iwp_max, vwp_min, vwp_mean, vwp_max, w_min, w_max

    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(rwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(rwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(rwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(iwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(iwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(iwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(swp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(swp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(swp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(gwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(gwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(gwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(tot_iwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(tot_iwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(tot_iwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(lathf_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lathf_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lathf_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(lwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(reske_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(reske_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(reske_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(senhf_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(senhf_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(senhf_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if (current_state%casim_enabled .eqv. .true.) then
        call mpi_recv(surface_precip_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(surface_precip_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(surface_precip_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(vwp_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vwp_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vwp_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(w_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(w_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(subke_2d_min, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(subke_2d_mean, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(subke_2d_max, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ls1 = len_trim(diagnostic_path)
    ls2 = 0
    do i = 1,ls1
        if(diagnostic_path(i:i).ne.' ') then
          ls2 = ls2 + 1
        end if
    end do

    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = diagnostic_path(:ls2)//"/full_diag_0d_instantaneous_time_"//trim(conv_to_string(&
                                                                                            int(current_state%time)))//".nc"
    else
      unique_filename = diagnostic_path(:ls2)//"/full_diag_0d_instantaneous_timestep_"//trim(conv_to_string(&
                                                                                              current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "time", 1, time_id))
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))

    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)
    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "rwp_min", "minimum rain water path", rwp_min, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "rwp_mean", "mean rain water path", rwp_mean, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "rwp_max", "maximum rain water path", rwp_max, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "iwp_min", "minimum ice water path", iwp_min, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "iwp_mean", "mean ice water path", iwp_mean, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "iwp_max", "maximum ice water path", iwp_max, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "swp_min", "minimum snow water path", swp_min, "g.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "swp_mean", "mean snow water path", swp_mean, "g.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "swp_max", "maximum snow water path", swp_max, "g.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "gwp_min", "minimum graupel water path", gwp_min, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "gwp_mean", "mean graupel water path", gwp_mean, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "gwp_max", "maximum graupel water path", gwp_max, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "tot_iwp_min", "minimum total ice water path", tot_iwp_min, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "tot_iwp_mean", "mean total ice water path", tot_iwp_mean, "kg.m-2", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "tot_iwp_max", "maximum total ice water path", tot_iwp_max, "kg.m-2", current_state%time)
      end if
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lathf_min", "minimum surface latent heat flux", lathf_min, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lathf_mean", "mean surface latent heat flux", lathf_mean, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lathf_max", "maximum surface latent heat flux", lathf_max, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lwp_min", "minimum liquid water path", lwp_min, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lwp_mean", "mean liquid water path", lwp_mean, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "lwp_max", "maximum liquid water path", lwp_max, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "reske_min", "minimum resolved ke", reske_min, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "reske_mean", "mean resolved ke", reske_mean, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "reske_max", "maximum resolved ke", reske_max, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "senhf_min", "minimum surface sensible heat flux", senhf_min, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "senhf_mean", "mean surface sensible heat flux", senhf_mean, "W.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "senhf_max", "maximum surface sensible heat flux", senhf_max, "W.m-2", current_state%time)
      if (current_state%casim_enabled .eqv. .true.) then
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "surface_precip_min", "minimum surface precipitation", surface_precip_min, "???", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "surface_precip_mean", "mean surface precipitation", surface_precip_mean, "???", current_state%time)
        call define_and_write_variable_real_scalar(ncdf_id, &
                        "surface_precip_max", "maximum surface precipitation", surface_precip_max, "???", current_state%time)
      end if
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "vwp_min", "minimum water vapour path", vwp_min, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "vwp_mean", "mean water vapour path", vwp_mean, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "vwp_max", "maximum water vapour path", vwp_max, "kg.m-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "w_min", "minimum vertical wind speed", w_min, "m.s-1", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "w_max", "maximum vertical wind speed", w_max, "m.s-1", current_state%time)
    end if
    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "subke_2d_min", "minimum subgrid ke", subke_2d_min, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "subke_2d_mean", "mean subgrid ke", subke_2d_mean, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                      "subke_2d_max", "maximum subgrid ke", subke_2d_max, "J or kg.m2.s-2", current_state%time)
    end if

    call check(nf90_close(ncdf_id))

    current_state%diag_0d_done = .true.
    call mpi_send(current_state%diag_0d_done, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, ierr)
    current_state%diag_0d_done = .false.

  end subroutine diagnostic_file_0d_generation





  subroutine diagnostic_file_1d_generation(current_state, vertical_grid, io_communicator_arg, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size, diagnostic_path)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           io_communicator_arg
    character(len=200), intent(in) :: diagnostic_path

    integer :: z_size_id, y_size_id, x_size_id, bin_size_id!, scalar_size_id, time_id
    integer :: ncdf_id, ierr, surface
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename

    global_grid_z_size_netcdf = global_grid_z_size
    surface = global_grid_x_size * global_grid_y_size
    size_bin = current_state%n_w_bins

    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      allocate(current_state%uwsg_tot(global_grid_z_size))
      allocate(current_state%vwsg_tot(global_grid_z_size))
      allocate(current_state%uusg_tot(global_grid_z_size))
      allocate(current_state%vvsg_tot(global_grid_z_size))
      allocate(current_state%wwsg_tot(global_grid_z_size))
      allocate(current_state%tkesg_tot(global_grid_z_size))
      allocate(current_state%wtsg_tot(global_grid_z_size))
      allocate(current_state%th2sg_tot(global_grid_z_size))
      allocate(current_state%sed_tot(global_grid_z_size))
      allocate(current_state%ssub_tot(global_grid_z_size))
      allocate(current_state%dissipation_tot(global_grid_z_size))
      allocate(current_state%buoysg_tot(global_grid_z_size))
      allocate(current_state%wkesg_tot(global_grid_z_size))
      allocate(current_state%theta_dis_tot(global_grid_z_size))
      allocate(current_state%vis_coef_tot(global_grid_z_size))
      allocate(current_state%diff_coef_tot(global_grid_z_size))
      allocate(current_state%richardson_number_tot(global_grid_z_size))
      allocate(current_state%richardson_squared_tot(global_grid_z_size))
      allocate(current_state%wqv_sg_tot(global_grid_z_size))
      allocate(current_state%wql_sg_tot(global_grid_z_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%wqr_sg_tot(global_grid_z_size))
        allocate(current_state%wqi_sg_tot(global_grid_z_size))
        allocate(current_state%wqs_sg_tot(global_grid_z_size))
        allocate(current_state%wqg_sg_tot(global_grid_z_size))
      end if
    end if

    if (current_state%casim_profile_dgs_enabled .eqv. .true.) then
      allocate(current_state%dqc_mphys_tot(global_grid_z_size))
      allocate(current_state%dqg_mphys_tot(global_grid_z_size))
      allocate(current_state%dqi_mphys_tot(global_grid_z_size))
      allocate(current_state%dqr_mphys_tot(global_grid_z_size))
      allocate(current_state%dqs_mphys_tot(global_grid_z_size))
      allocate(current_state%dqv_mphys_tot(global_grid_z_size))
      allocate(current_state%dth_mphys_tot(global_grid_z_size))
      allocate(current_state%dth_cond_evap_tot(global_grid_z_size))
      allocate(current_state%dqv_cond_evap_tot(global_grid_z_size))
      allocate(current_state%phomc_tot(global_grid_z_size))
      allocate(current_state%pinuc_tot(global_grid_z_size))
      allocate(current_state%pidep_tot(global_grid_z_size))
      allocate(current_state%psdep_tot(global_grid_z_size))
      allocate(current_state%piacw_tot(global_grid_z_size))
      allocate(current_state%psacw_tot(global_grid_z_size))
      allocate(current_state%psacr_tot(global_grid_z_size))
      allocate(current_state%pisub_tot(global_grid_z_size))
      allocate(current_state%pssub_tot(global_grid_z_size))
      allocate(current_state%pimlt_tot(global_grid_z_size))
      allocate(current_state%psmlt_tot(global_grid_z_size))
      allocate(current_state%psaut_tot(global_grid_z_size))
      allocate(current_state%psaci_tot(global_grid_z_size))
      allocate(current_state%praut_tot(global_grid_z_size))
      allocate(current_state%pracw_tot(global_grid_z_size))
      allocate(current_state%prevp_tot(global_grid_z_size))
      allocate(current_state%pgacw_tot(global_grid_z_size))
      allocate(current_state%pgacs_tot(global_grid_z_size))
      allocate(current_state%pgmlt_tot(global_grid_z_size))
      allocate(current_state%pgsub_tot(global_grid_z_size))
      allocate(current_state%psedi_tot(global_grid_z_size))
      allocate(current_state%pseds_tot(global_grid_z_size))
      allocate(current_state%psedr_tot(global_grid_z_size))
      allocate(current_state%psedg_tot(global_grid_z_size))
      allocate(current_state%psedl_tot(global_grid_z_size))
      allocate(current_state%pcond_tot(global_grid_z_size))
    end if

    if (current_state%profile_diagnostics_enabled .eqv. .true.) then
      allocate(current_state%u_wind_tot(global_grid_z_size))
      allocate(current_state%uprime_tot(global_grid_z_size))
      allocate(current_state%v_wind_tot(global_grid_z_size))
      allocate(current_state%vprime_tot(global_grid_z_size))
      allocate(current_state%wke_tot(global_grid_z_size))
      allocate(current_state%ww_tot(global_grid_z_size))
      allocate(current_state%www_tot(global_grid_z_size))
      allocate(current_state%wwww_tot(global_grid_z_size))
      allocate(current_state%theta_tot(global_grid_z_size))
      allocate(current_state%w_wind_tot(global_grid_z_size))
      allocate(current_state%rh_tot(global_grid_z_size))
      allocate(current_state%wtheta_ad_tot(global_grid_z_size))
      allocate(current_state%wtheta_cn_tot(global_grid_z_size))
      allocate(current_state%uw_tot(global_grid_z_size))
      allocate(current_state%vw_tot(global_grid_z_size))
      allocate(current_state%uv_tot(global_grid_z_size))
      allocate(current_state%th2_tot(global_grid_z_size))
      allocate(current_state%thinit(global_grid_z_size))
      allocate(current_state%qv_tot(global_grid_z_size))
      allocate(current_state%ql_tot(global_grid_z_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%qr_tot(global_grid_z_size))
        allocate(current_state%qi_tot(global_grid_z_size))
        allocate(current_state%qs_tot(global_grid_z_size))
        allocate(current_state%qg_tot(global_grid_z_size))
      end if
      allocate(current_state%wqv_cn_tot(global_grid_z_size))
      allocate(current_state%wql_cn_tot(global_grid_z_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%wqr_cn_tot(global_grid_z_size))
        allocate(current_state%wqi_cn_tot(global_grid_z_size))
        allocate(current_state%wqs_cn_tot(global_grid_z_size))
        allocate(current_state%wqg_cn_tot(global_grid_z_size))
      end if
      allocate(current_state%wqv_ad_tot(global_grid_z_size))
      allocate(current_state%wql_ad_tot(global_grid_z_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%wqr_ad_tot(global_grid_z_size))
        allocate(current_state%wqi_ad_tot(global_grid_z_size))
        allocate(current_state%wqs_ad_tot(global_grid_z_size))
        allocate(current_state%wqg_ad_tot(global_grid_z_size))
        allocate(current_state%cloud_mask_tot(global_grid_z_size))
        allocate(current_state%cloud_liq_mask_tot(global_grid_z_size))
        allocate(current_state%cloud_ice_mask_tot(global_grid_z_size))
      end if
    end if

    if (current_state%forcing_enabled .eqv. .true.) then
      allocate(current_state%du_subs_profile_diag(global_grid_z_size))
      allocate(current_state%dv_subs_profile_diag(global_grid_z_size))
      allocate(current_state%dtheta_subs_profile_diag(global_grid_z_size))
      allocate(current_state%dqv_subs_profile_diag(global_grid_z_size))
      allocate(current_state%dql_subs_profile_diag(global_grid_z_size))
      allocate(current_state%dqr_subs_profile_diag(global_grid_z_size))
      allocate(current_state%dqi_subs_profile_diag(global_grid_z_size))
      allocate(current_state%dqs_subs_profile_diag(global_grid_z_size))
      allocate(current_state%dqg_subs_profile_diag(global_grid_z_size))
      allocate(current_state%tend_pr_tot_u_forc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_v_forc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_th_forc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qv_forc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_ql_forc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qi_forc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qr_forc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qs_forc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qg_forc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_tabs_forc(global_grid_z_size))
    end if

    if (current_state%diffusion_enabled .eqv. .true.) then
      allocate(current_state%tend_pr_tot_th_diff(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qv_diff(global_grid_z_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%tend_pr_tot_ql_diff(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qi_diff(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qr_diff(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qs_diff(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qg_diff(global_grid_z_size))
      end if
      allocate(current_state%tend_pr_tot_tabs_diff(global_grid_z_size))
    end if

    if (current_state%buoyancy_enabled .eqv. .true.) then
      allocate(current_state%tend_pr_tot_w_buoy(global_grid_z_size))
    end if

    if (current_state%coriolis_enabled .eqv. .true.) then
      allocate(current_state%tend_pr_tot_u_corio(global_grid_z_size))
      allocate(current_state%tend_pr_tot_v_corio(global_grid_z_size))
    end if

    if (current_state%pstep_enabled .eqv. .true.) then
      allocate(current_state%tendp_pr_tot_u_pt(global_grid_z_size))
      allocate(current_state%tendp_pr_tot_v_pt(global_grid_z_size))
      allocate(current_state%tendp_pr_tot_w_pt(global_grid_z_size))
    end if

    if (current_state%pw_advection_enabled .eqv. .true.) then
      allocate(current_state%tend_pr_tot_u_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_v_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_w_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_th_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qv_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_ql_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qi_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qr_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qs_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qg_pwad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_tabs_pwad(global_grid_z_size))
    end if

    if (current_state%stepfields_enabled .eqv. .true.) then
      allocate(current_state%tend_pr_tot_th_sf(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qv_sf(global_grid_z_size))
      allocate(current_state%tend_pr_tot_ql_sf(global_grid_z_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%tend_pr_tot_qi_sf(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qr_sf(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qs_sf(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qg_sf(global_grid_z_size))
      end if
      allocate(current_state%tend_pr_tot_tabs_sf(global_grid_z_size))
    end if

    if (current_state%th_advection_enabled .eqv. .true.) then
      allocate(current_state%tend_pr_tot_th_thad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_tabs_thad(global_grid_z_size))
    end if

    if (current_state%tvd_advection_enabled .eqv. .true.) then
      allocate(current_state%tend_pr_tot_u_tvad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_v_tvad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_w_tvad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_th_tvad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_qv_tvad(global_grid_z_size))
      allocate(current_state%tend_pr_tot_ql_tvad(global_grid_z_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%tend_pr_tot_qi_tvad(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qr_tvad(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qs_tvad(global_grid_z_size))
        allocate(current_state%tend_pr_tot_qg_tvad(global_grid_z_size))
      end if
      allocate(current_state%tend_pr_tot_tabs_tvad(global_grid_z_size))
    end if

    if (current_state%viscosity_enabled .eqv. .true.) then
      allocate(current_state%tend_pr_tot_u_visc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_v_visc(global_grid_z_size))
      allocate(current_state%tend_pr_tot_w_visc(global_grid_z_size))
    end if

    if (current_state%socrates_enabled .eqv. .true.) then
      allocate(current_state%cloud_reff_tot(global_grid_z_size))
      allocate(current_state%longwave_hr_tot(global_grid_z_size))
      allocate(current_state%shortwave_hr_tot(global_grid_z_size))
      allocate(current_state%tend_pr_tot_th_lw(global_grid_z_size))
      allocate(current_state%tend_pr_tot_tabs_lw(global_grid_z_size))
      allocate(current_state%tend_pr_tot_th_sw(global_grid_z_size))
      allocate(current_state%tend_pr_tot_tabs_sw(global_grid_z_size))
      allocate(current_state%tend_pr_tot_th_total(global_grid_z_size))
      allocate(current_state%tend_pr_tot_tabs_total(global_grid_z_size))
    end if

    if (current_state%flux_budget_enabled .eqv. .true.) then
      allocate(current_state%th_flux_values(global_grid_z_size))
      allocate(current_state%th_gradient(global_grid_z_size))
      allocate(current_state%th_diff(global_grid_z_size))
      allocate(current_state%th_buoyancy(global_grid_z_size))
      allocate(current_state%th_tendency(global_grid_z_size))

      allocate(current_state%qv_flux_values(global_grid_z_size))
      allocate(current_state%ql_flux_values(global_grid_z_size))
      allocate(current_state%qr_flux_values(global_grid_z_size))
      allocate(current_state%qi_flux_values(global_grid_z_size))
      allocate(current_state%qs_flux_values(global_grid_z_size))
      allocate(current_state%qg_flux_values(global_grid_z_size))
      allocate(current_state%qAitkenSolMass_flux_values(global_grid_z_size))
      allocate(current_state%qAccumSolMass_flux_values(global_grid_z_size))
      allocate(current_state%qAccumInsolMass_flux_values(global_grid_z_size))
      allocate(current_state%qCoarseSolMass_flux_values(global_grid_z_size))
      allocate(current_state%qCoarseDustMass_flux_values(global_grid_z_size))
      allocate(current_state%nl_flux_values(global_grid_z_size))
      allocate(current_state%nr_flux_values(global_grid_z_size))
      allocate(current_state%ni_flux_values(global_grid_z_size))
      allocate(current_state%ns_flux_values(global_grid_z_size))
      allocate(current_state%ng_flux_values(global_grid_z_size))
      allocate(current_state%nAitkenSolNumber_flux_values(global_grid_z_size))
      allocate(current_state%nAccumSolNumber_flux_values(global_grid_z_size))
      allocate(current_state%nAccumInsolNumber_flux_values(global_grid_z_size))
      allocate(current_state%nCoarseSolNumber_flux_values(global_grid_z_size))
      allocate(current_state%nCoarseDustnumber_flux_values(global_grid_z_size))

      allocate(current_state%qv_gradient(global_grid_z_size))
      allocate(current_state%ql_gradient(global_grid_z_size))
      allocate(current_state%qr_gradient(global_grid_z_size))
      allocate(current_state%qi_gradient(global_grid_z_size))
      allocate(current_state%qs_gradient(global_grid_z_size))
      allocate(current_state%qg_gradient(global_grid_z_size))
      allocate(current_state%qAitkenSolMass_gradient(global_grid_z_size))
      allocate(current_state%qAccumSolMass_gradient(global_grid_z_size))
      allocate(current_state%qAccumInsolMass_gradient(global_grid_z_size))
      allocate(current_state%qCoarseSolMass_gradient(global_grid_z_size))
      allocate(current_state%qCoarseDustMass_gradient(global_grid_z_size))
      allocate(current_state%nl_gradient(global_grid_z_size))
      allocate(current_state%nr_gradient(global_grid_z_size))
      allocate(current_state%ni_gradient(global_grid_z_size))
      allocate(current_state%ns_gradient(global_grid_z_size))
      allocate(current_state%ng_gradient(global_grid_z_size))
      allocate(current_state%nAitkenSolNumber_gradient(global_grid_z_size))
      allocate(current_state%nAccumSolNumber_gradient(global_grid_z_size))
      allocate(current_state%nAccumInsolNumber_gradient(global_grid_z_size))
      allocate(current_state%nCoarseSolNumber_gradient(global_grid_z_size))
      allocate(current_state%nCoarseDustnumber_gradient(global_grid_z_size))

      allocate(current_state%qv_diff(global_grid_z_size))
      allocate(current_state%ql_diff(global_grid_z_size))
      allocate(current_state%qr_diff(global_grid_z_size))
      allocate(current_state%qi_diff(global_grid_z_size))
      allocate(current_state%qs_diff(global_grid_z_size))
      allocate(current_state%qg_diff(global_grid_z_size))
      allocate(current_state%qAitkenSolMass_diff(global_grid_z_size))
      allocate(current_state%qAccumSolMass_diff(global_grid_z_size))
      allocate(current_state%qAccumInsolMass_diff(global_grid_z_size))
      allocate(current_state%qCoarseSolMass_diff(global_grid_z_size))
      allocate(current_state%qCoarseDustMass_diff(global_grid_z_size))
      allocate(current_state%nl_diff(global_grid_z_size))
      allocate(current_state%nr_diff(global_grid_z_size))
      allocate(current_state%ni_diff(global_grid_z_size))
      allocate(current_state%ns_diff(global_grid_z_size))
      allocate(current_state%ng_diff(global_grid_z_size))
      allocate(current_state%nAitkenSolNumber_diff(global_grid_z_size))
      allocate(current_state%nAccumSolNumber_diff(global_grid_z_size))
      allocate(current_state%nAccumInsolNumber_diff(global_grid_z_size))
      allocate(current_state%nCoarseSolNumber_diff(global_grid_z_size))
      allocate(current_state%nCoarseDustnumber_diff(global_grid_z_size))

      allocate(current_state%qv_buoyancy(global_grid_z_size))
      allocate(current_state%ql_buoyancy(global_grid_z_size))
      allocate(current_state%qr_buoyancy(global_grid_z_size))
      allocate(current_state%qi_buoyancy(global_grid_z_size))
      allocate(current_state%qs_buoyancy(global_grid_z_size))
      allocate(current_state%qg_buoyancy(global_grid_z_size))
      allocate(current_state%qAitkenSolMass_buoyancy(global_grid_z_size))
      allocate(current_state%qAccumSolMass_buoyancy(global_grid_z_size))
      allocate(current_state%qAccumInsolMass_buoyancy(global_grid_z_size))
      allocate(current_state%qCoarseSolMass_buoyancy(global_grid_z_size))
      allocate(current_state%qCoarseDustMass_buoyancy(global_grid_z_size))
      allocate(current_state%nl_buoyancy(global_grid_z_size))
      allocate(current_state%nr_buoyancy(global_grid_z_size))
      allocate(current_state%ni_buoyancy(global_grid_z_size))
      allocate(current_state%ns_buoyancy(global_grid_z_size))
      allocate(current_state%ng_buoyancy(global_grid_z_size))
      allocate(current_state%nAitkenSolNumber_buoyancy(global_grid_z_size))
      allocate(current_state%nAccumSolNumber_buoyancy(global_grid_z_size))
      allocate(current_state%nAccumInsolNumber_buoyancy(global_grid_z_size))
      allocate(current_state%nCoarseSolNumber_buoyancy(global_grid_z_size))
      allocate(current_state%nCoarseDustnumber_buoyancy(global_grid_z_size))

      allocate(current_state%qv_tendency(global_grid_z_size))
      allocate(current_state%ql_tendency(global_grid_z_size))
      allocate(current_state%qr_tendency(global_grid_z_size))
      allocate(current_state%qi_tendency(global_grid_z_size))
      allocate(current_state%qs_tendency(global_grid_z_size))
      allocate(current_state%qg_tendency(global_grid_z_size))
      allocate(current_state%qAitkenSolMass_tendency(global_grid_z_size))
      allocate(current_state%qAccumSolMass_tendency(global_grid_z_size))
      allocate(current_state%qAccumInsolMass_tendency(global_grid_z_size))
      allocate(current_state%qCoarseSolMass_tendency(global_grid_z_size))
      allocate(current_state%qCoarseDustMass_tendency(global_grid_z_size))
      allocate(current_state%nl_tendency(global_grid_z_size))
      allocate(current_state%nr_tendency(global_grid_z_size))
      allocate(current_state%ni_tendency(global_grid_z_size))
      allocate(current_state%ns_tendency(global_grid_z_size))
      allocate(current_state%ng_tendency(global_grid_z_size))
      allocate(current_state%nAitkenSolNumber_tendency(global_grid_z_size))
      allocate(current_state%nAccumSolNumber_tendency(global_grid_z_size))
      allocate(current_state%nAccumInsolNumber_tendency(global_grid_z_size))
      allocate(current_state%nCoarseSolNumber_tendency(global_grid_z_size))
      allocate(current_state%nCoarseDustnumber_tendency(global_grid_z_size))

      allocate(current_state%uw_advection(global_grid_z_size))
      allocate(current_state%vw_advection(global_grid_z_size))
      allocate(current_state%uw_viscosity(global_grid_z_size))
      allocate(current_state%vw_viscosity(global_grid_z_size))
      allocate(current_state%uw_buoyancy(global_grid_z_size))
      allocate(current_state%vw_buoyancy(global_grid_z_size))
      allocate(current_state%uw_tendency(global_grid_z_size))
      allocate(current_state%vw_tendency(global_grid_z_size))
      allocate(current_state%uw_w(global_grid_z_size))
      allocate(current_state%vw_w(global_grid_z_size))

      allocate(current_state%tu_su(global_grid_z_size))
      allocate(current_state%uu_advection(global_grid_z_size))
      allocate(current_state%uu_viscosity(global_grid_z_size))
      allocate(current_state%wu_u(global_grid_z_size))
      allocate(current_state%tv_sv(global_grid_z_size))
      allocate(current_state%vv_advection(global_grid_z_size))
      allocate(current_state%vv_viscosity(global_grid_z_size))
      allocate(current_state%wv_v(global_grid_z_size))
      allocate(current_state%tw_sw(global_grid_z_size))
      allocate(current_state%ww_advection(global_grid_z_size))
      allocate(current_state%ww_viscosity(global_grid_z_size))
      allocate(current_state%ww_buoyancy(global_grid_z_size))

      allocate(current_state%u_thetal(global_grid_z_size))
      allocate(current_state%us_thetal(global_grid_z_size))
      allocate(current_state%u_thetal_advection(global_grid_z_size))
      allocate(current_state%u_thetal_viscosity_diffusion(global_grid_z_size))
      allocate(current_state%wu_thetal(global_grid_z_size))
      allocate(current_state%v_thetal(global_grid_z_size))
      allocate(current_state%vs_thetal(global_grid_z_size))
      allocate(current_state%v_thetal_advection(global_grid_z_size))
      allocate(current_state%v_thetal_viscosity_diffusion(global_grid_z_size))
      allocate(current_state%wv_thetal(global_grid_z_size))
      allocate(current_state%w_thetal(global_grid_z_size))
      allocate(current_state%ws_thetal(global_grid_z_size))
      allocate(current_state%w_thetal_advection(global_grid_z_size))
      allocate(current_state%w_thetal_viscosity_diffusion(global_grid_z_size))
      allocate(current_state%w_thetal_buoyancy(global_grid_z_size))
      allocate(current_state%ww_thetal(global_grid_z_size))
      allocate(current_state%thetal_thetal(global_grid_z_size))
      allocate(current_state%sthetal_thetal(global_grid_z_size))
      allocate(current_state%thetal_thetal_advection(global_grid_z_size))
      allocate(current_state%thetal_thetal_diffusion(global_grid_z_size))
      allocate(current_state%wthetal_thetal(global_grid_z_size))

      allocate(current_state%shprd(global_grid_z_size))
      allocate(current_state%w_ke(global_grid_z_size))
      allocate(current_state%w_p(global_grid_z_size))
      allocate(current_state%sub_buoy(global_grid_z_size))
      allocate(current_state%tke_tend(global_grid_z_size))
    end if

    if (current_state%conditional_diagnostics_column_enabled .eqv. .true.) then
      allocate(current_state%w_zn_cln(global_grid_z_size))
      allocate(current_state%w_zn2_cln(global_grid_z_size))
      allocate(current_state%tmp_th_cln(global_grid_z_size))
      allocate(current_state%wth_cln(global_grid_z_size))
      allocate(current_state%th_pr_cln(global_grid_z_size))
      allocate(current_state%wth_pr_cln(global_grid_z_size))
      allocate(current_state%thv_pr_cln(global_grid_z_size))
      allocate(current_state%wthv_pr_cln(global_grid_z_size))
      allocate(current_state%th_pr2_cln(global_grid_z_size))
      allocate(current_state%wthsg_cln(global_grid_z_size))
      allocate(current_state%w_zn3_cln(global_grid_z_size))
      allocate(current_state%relhum_cln(global_grid_z_size))
      allocate(current_state%tmp_u_cln(global_grid_z_size))
      allocate(current_state%tmp_v_cln(global_grid_z_size))
      allocate(current_state%wu_cln(global_grid_z_size))
      allocate(current_state%wv_cln(global_grid_z_size))
      allocate(current_state%wusg_cln(global_grid_z_size))
      allocate(current_state%wvsg_cln(global_grid_z_size))
      allocate(current_state%TdegK_cln(global_grid_z_size))
      allocate(current_state%th_h_cln(global_grid_z_size))
      allocate(current_state%th_h_pr1_cln(global_grid_z_size))
      allocate(current_state%th_h_pr2_cln(global_grid_z_size))
      allocate(current_state%qvli_cln(global_grid_z_size))
      allocate(current_state%qvli_pr_cln(global_grid_z_size))
      allocate(current_state%qvli_pr2_cln(global_grid_z_size))
      allocate(current_state%qppt_cln(global_grid_z_size))
      allocate(current_state%qppt_pr_cln(global_grid_z_size))
      allocate(current_state%qppt_pr2_cln(global_grid_z_size))
      allocate(current_state%wqvli_pr_cln(global_grid_z_size))
      allocate(current_state%wqppt_pr_cln(global_grid_z_size))
    end if

    if (current_state%pdf_analysis_enabled .eqv. .true.) then
      allocate(current_state%global_grid%configuration%vertical%w_dwn(global_grid_z_size))
      allocate(current_state%global_grid%configuration%vertical%w_up(global_grid_z_size))
      allocate(current_state%w_histogram_profile_local(global_grid_z_size, current_state%n_w_bins))
    end if


    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      current_state%uwsg_tot = 0.0_DEFAULT_PRECISION
      current_state%vwsg_tot = 0.0_DEFAULT_PRECISION
      current_state%uusg_tot = 0.0_DEFAULT_PRECISION
      current_state%vvsg_tot = 0.0_DEFAULT_PRECISION
      current_state%wwsg_tot = 0.0_DEFAULT_PRECISION
      current_state%tkesg_tot = 0.0_DEFAULT_PRECISION
      current_state%wtsg_tot = 0.0_DEFAULT_PRECISION
      current_state%th2sg_tot = 0.0_DEFAULT_PRECISION
      current_state%sed_tot = 0.0_DEFAULT_PRECISION
      current_state%ssub_tot = 0.0_DEFAULT_PRECISION
      current_state%dissipation_tot = 0.0_DEFAULT_PRECISION
      current_state%buoysg_tot = 0.0_DEFAULT_PRECISION
      current_state%wkesg_tot = 0.0_DEFAULT_PRECISION
      current_state%theta_dis_tot = 0.0_DEFAULT_PRECISION
      current_state%vis_coef_tot = 0.0_DEFAULT_PRECISION
      current_state%diff_coef_tot = 0.0_DEFAULT_PRECISION
      current_state%richardson_number_tot = 0.0_DEFAULT_PRECISION
      current_state%richardson_squared_tot = 0.0_DEFAULT_PRECISION
      current_state%wqv_sg_tot = 0.0_DEFAULT_PRECISION
      current_state%wql_sg_tot = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%wqr_sg_tot = 0.0_DEFAULT_PRECISION
        current_state%wqi_sg_tot = 0.0_DEFAULT_PRECISION
        current_state%wqs_sg_tot = 0.0_DEFAULT_PRECISION
        current_state%wqg_sg_tot = 0.0_DEFAULT_PRECISION
      end if
    end if

    if (current_state%casim_profile_dgs_enabled .eqv. .true.) then
      current_state%dqc_mphys_tot = 0.0_DEFAULT_PRECISION
      current_state%dqg_mphys_tot = 0.0_DEFAULT_PRECISION
      current_state%dqi_mphys_tot = 0.0_DEFAULT_PRECISION
      current_state%dqr_mphys_tot = 0.0_DEFAULT_PRECISION
      current_state%dqs_mphys_tot = 0.0_DEFAULT_PRECISION
      current_state%dqv_mphys_tot = 0.0_DEFAULT_PRECISION
      current_state%dth_mphys_tot = 0.0_DEFAULT_PRECISION
      current_state%dth_cond_evap_tot = 0.0_DEFAULT_PRECISION
      current_state%dqv_cond_evap_tot = 0.0_DEFAULT_PRECISION
      current_state%phomc_tot = 0.0_DEFAULT_PRECISION
      current_state%pinuc_tot = 0.0_DEFAULT_PRECISION
      current_state%pidep_tot = 0.0_DEFAULT_PRECISION
      current_state%psdep_tot = 0.0_DEFAULT_PRECISION
      current_state%piacw_tot = 0.0_DEFAULT_PRECISION
      current_state%psacw_tot = 0.0_DEFAULT_PRECISION
      current_state%psacr_tot = 0.0_DEFAULT_PRECISION
      current_state%pisub_tot = 0.0_DEFAULT_PRECISION
      current_state%pssub_tot = 0.0_DEFAULT_PRECISION
      current_state%pimlt_tot = 0.0_DEFAULT_PRECISION
      current_state%psmlt_tot = 0.0_DEFAULT_PRECISION
      current_state%psaut_tot = 0.0_DEFAULT_PRECISION
      current_state%psaci_tot = 0.0_DEFAULT_PRECISION
      current_state%praut_tot = 0.0_DEFAULT_PRECISION
      current_state%pracw_tot = 0.0_DEFAULT_PRECISION
      current_state%prevp_tot = 0.0_DEFAULT_PRECISION
      current_state%pgacw_tot = 0.0_DEFAULT_PRECISION
      current_state%pgacs_tot = 0.0_DEFAULT_PRECISION
      current_state%pgmlt_tot = 0.0_DEFAULT_PRECISION
      current_state%pgsub_tot = 0.0_DEFAULT_PRECISION
      current_state%psedi_tot = 0.0_DEFAULT_PRECISION
      current_state%pseds_tot = 0.0_DEFAULT_PRECISION
      current_state%psedr_tot = 0.0_DEFAULT_PRECISION
      current_state%psedg_tot = 0.0_DEFAULT_PRECISION
      current_state%psedl_tot = 0.0_DEFAULT_PRECISION
      current_state%pcond_tot = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%profile_diagnostics_enabled .eqv. .true.) then
      current_state%u_wind_tot = 0.0_DEFAULT_PRECISION
      current_state%uprime_tot = 0.0_DEFAULT_PRECISION
      current_state%v_wind_tot = 0.0_DEFAULT_PRECISION
      current_state%vprime_tot = 0.0_DEFAULT_PRECISION
      current_state%wke_tot = 0.0_DEFAULT_PRECISION
      current_state%ww_tot = 0.0_DEFAULT_PRECISION
      current_state%www_tot = 0.0_DEFAULT_PRECISION
      current_state%wwww_tot = 0.0_DEFAULT_PRECISION
      current_state%theta_tot = 0.0_DEFAULT_PRECISION
      current_state%w_wind_tot = 0.0_DEFAULT_PRECISION
      current_state%rh_tot = 0.0_DEFAULT_PRECISION
      current_state%wtheta_ad_tot = 0.0_DEFAULT_PRECISION
      current_state%wtheta_cn_tot = 0.0_DEFAULT_PRECISION
      current_state%uw_tot = 0.0_DEFAULT_PRECISION
      current_state%vw_tot = 0.0_DEFAULT_PRECISION
      current_state%uv_tot = 0.0_DEFAULT_PRECISION
      current_state%th2_tot = 0.0_DEFAULT_PRECISION
      current_state%thinit = 0.0_DEFAULT_PRECISION
      current_state%qv_tot = 0.0_DEFAULT_PRECISION
      current_state%ql_tot = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%qr_tot = 0.0_DEFAULT_PRECISION
        current_state%qi_tot = 0.0_DEFAULT_PRECISION
        current_state%qs_tot = 0.0_DEFAULT_PRECISION
        current_state%qg_tot = 0.0_DEFAULT_PRECISION
      end if
      current_state%wqv_cn_tot = 0.0_DEFAULT_PRECISION
      current_state%wql_cn_tot = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%wqr_cn_tot = 0.0_DEFAULT_PRECISION
        current_state%wqi_cn_tot = 0.0_DEFAULT_PRECISION
        current_state%wqs_cn_tot = 0.0_DEFAULT_PRECISION
        current_state%wqg_cn_tot = 0.0_DEFAULT_PRECISION
      end if
      current_state%wqv_ad_tot = 0.0_DEFAULT_PRECISION
      current_state%wql_ad_tot = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%wqr_ad_tot = 0.0_DEFAULT_PRECISION
        current_state%wqi_ad_tot = 0.0_DEFAULT_PRECISION
        current_state%wqs_ad_tot = 0.0_DEFAULT_PRECISION
        current_state%wqg_ad_tot = 0.0_DEFAULT_PRECISION
        current_state%cloud_mask_tot = 0.0_DEFAULT_PRECISION
        current_state%cloud_liq_mask_tot = 0.0_DEFAULT_PRECISION
        current_state%cloud_ice_mask_tot = 0.0_DEFAULT_PRECISION
      end if
    end if

    if (current_state%forcing_enabled .eqv. .true.) then
      current_state%du_subs_profile_diag = 0.0_DEFAULT_PRECISION
      current_state%dv_subs_profile_diag = 0.0_DEFAULT_PRECISION
      current_state%dtheta_subs_profile_diag = 0.0_DEFAULT_PRECISION
      current_state%dqv_subs_profile_diag = 0.0_DEFAULT_PRECISION
      current_state%dql_subs_profile_diag = 0.0_DEFAULT_PRECISION
      current_state%dqr_subs_profile_diag = 0.0_DEFAULT_PRECISION
      current_state%dqi_subs_profile_diag = 0.0_DEFAULT_PRECISION
      current_state%dqs_subs_profile_diag = 0.0_DEFAULT_PRECISION
      current_state%dqg_subs_profile_diag = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_u_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_v_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_th_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qv_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_ql_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qi_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qr_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qs_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qg_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_tabs_forc = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%diffusion_enabled .eqv. .true.) then
      current_state%tend_pr_tot_th_diff = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qv_diff = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_ql_diff = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%tend_pr_tot_qi_diff = 0.0_DEFAULT_PRECISION
        current_state%tend_pr_tot_qr_diff = 0.0_DEFAULT_PRECISION
        current_state%tend_pr_tot_qs_diff = 0.0_DEFAULT_PRECISION
        current_state%tend_pr_tot_qg_diff = 0.0_DEFAULT_PRECISION
      end if
      current_state%tend_pr_tot_tabs_diff = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%buoyancy_enabled .eqv. .true.) then
      current_state%tend_pr_tot_w_buoy = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%coriolis_enabled .eqv. .true.) then
      current_state%tend_pr_tot_u_corio = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_v_corio = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%pstep_enabled .eqv. .true.) then
      current_state%tendp_pr_tot_u_pt = 0.0_DEFAULT_PRECISION
      current_state%tendp_pr_tot_v_pt = 0.0_DEFAULT_PRECISION
      current_state%tendp_pr_tot_w_pt = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%pw_advection_enabled .eqv. .true.) then
      current_state%tend_pr_tot_u_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_v_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_w_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_th_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qv_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_ql_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qi_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qr_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qs_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qg_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_tabs_pwad = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%stepfields_enabled .eqv. .true.) then
      current_state%tend_pr_tot_th_sf = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qv_sf = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_ql_sf = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%tend_pr_tot_qi_sf = 0.0_DEFAULT_PRECISION
        current_state%tend_pr_tot_qr_sf = 0.0_DEFAULT_PRECISION
        current_state%tend_pr_tot_qs_sf = 0.0_DEFAULT_PRECISION
        current_state%tend_pr_tot_qg_sf = 0.0_DEFAULT_PRECISION
      end if
      current_state%tend_pr_tot_tabs_sf = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%th_advection_enabled .eqv. .true.) then
      current_state%tend_pr_tot_th_thad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_tabs_thad = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%tvd_advection_enabled .eqv. .true.) then
      current_state%tend_pr_tot_u_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_v_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_w_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_th_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_qv_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_ql_tvad = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%tend_pr_tot_qi_tvad = 0.0_DEFAULT_PRECISION
        current_state%tend_pr_tot_qr_tvad = 0.0_DEFAULT_PRECISION
        current_state%tend_pr_tot_qs_tvad = 0.0_DEFAULT_PRECISION
        current_state%tend_pr_tot_qg_tvad = 0.0_DEFAULT_PRECISION
      end if
      current_state%tend_pr_tot_tabs_tvad = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%viscosity_enabled .eqv. .true.) then
      current_state%tend_pr_tot_u_visc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_v_visc = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_w_visc = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%socrates_enabled .eqv. .true.) then
      current_state%cloud_reff_tot = 0.0_DEFAULT_PRECISION
      current_state%longwave_hr_tot = 0.0_DEFAULT_PRECISION
      current_state%shortwave_hr_tot = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_th_lw = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_tabs_lw = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_th_sw = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_tabs_sw = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_th_total = 0.0_DEFAULT_PRECISION
      current_state%tend_pr_tot_tabs_total = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%flux_budget_enabled .eqv. .true.) then
      current_state%th_flux_values = 0.0_DEFAULT_PRECISION
      current_state%th_gradient = 0.0_DEFAULT_PRECISION
      current_state%th_diff = 0.0_DEFAULT_PRECISION
      current_state%th_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%th_tendency = 0.0_DEFAULT_PRECISION

      current_state%qv_flux_values = 0.0_DEFAULT_PRECISION
      current_state%ql_flux_values = 0.0_DEFAULT_PRECISION
      current_state%qr_flux_values = 0.0_DEFAULT_PRECISION
      current_state%qi_flux_values = 0.0_DEFAULT_PRECISION
      current_state%qs_flux_values = 0.0_DEFAULT_PRECISION
      current_state%qg_flux_values = 0.0_DEFAULT_PRECISION
      current_state%qAitkenSolMass_flux_values = 0.0_DEFAULT_PRECISION
      current_state%qAccumSolMass_flux_values = 0.0_DEFAULT_PRECISION
      current_state%qAccumInsolMass_flux_values = 0.0_DEFAULT_PRECISION
      current_state%qCoarseSolMass_flux_values = 0.0_DEFAULT_PRECISION
      current_state%qCoarseDustMass_flux_values = 0.0_DEFAULT_PRECISION
      current_state%nl_flux_values = 0.0_DEFAULT_PRECISION
      current_state%nr_flux_values = 0.0_DEFAULT_PRECISION
      current_state%ni_flux_values = 0.0_DEFAULT_PRECISION
      current_state%ns_flux_values = 0.0_DEFAULT_PRECISION
      current_state%ng_flux_values = 0.0_DEFAULT_PRECISION
      current_state%nAitkenSolNumber_flux_values = 0.0_DEFAULT_PRECISION
      current_state%nAccumSolNumber_flux_values = 0.0_DEFAULT_PRECISION
      current_state%nAccumInsolNumber_flux_values = 0.0_DEFAULT_PRECISION
      current_state%nCoarseSolNumber_flux_values = 0.0_DEFAULT_PRECISION
      current_state%nCoarseDustnumber_flux_values = 0.0_DEFAULT_PRECISION

      current_state%qv_gradient = 0.0_DEFAULT_PRECISION
      current_state%ql_gradient = 0.0_DEFAULT_PRECISION
      current_state%qr_gradient = 0.0_DEFAULT_PRECISION
      current_state%qi_gradient = 0.0_DEFAULT_PRECISION
      current_state%qs_gradient = 0.0_DEFAULT_PRECISION
      current_state%qg_gradient = 0.0_DEFAULT_PRECISION
      current_state%qAitkenSolMass_gradient = 0.0_DEFAULT_PRECISION
      current_state%qAccumSolMass_gradient = 0.0_DEFAULT_PRECISION
      current_state%qAccumInsolMass_gradient = 0.0_DEFAULT_PRECISION
      current_state%qCoarseSolMass_gradient = 0.0_DEFAULT_PRECISION
      current_state%qCoarseDustMass_gradient = 0.0_DEFAULT_PRECISION
      current_state%nl_gradient = 0.0_DEFAULT_PRECISION
      current_state%nr_gradient = 0.0_DEFAULT_PRECISION
      current_state%ni_gradient = 0.0_DEFAULT_PRECISION
      current_state%ns_gradient = 0.0_DEFAULT_PRECISION
      current_state%ng_gradient = 0.0_DEFAULT_PRECISION
      current_state%nAitkenSolNumber_gradient = 0.0_DEFAULT_PRECISION
      current_state%nAccumSolNumber_gradient = 0.0_DEFAULT_PRECISION
      current_state%nAccumInsolNumber_gradient = 0.0_DEFAULT_PRECISION
      current_state%nCoarseSolNumber_gradient = 0.0_DEFAULT_PRECISION
      current_state%nCoarseDustnumber_gradient = 0.0_DEFAULT_PRECISION

      current_state%qv_diff = 0.0_DEFAULT_PRECISION
      current_state%ql_diff = 0.0_DEFAULT_PRECISION
      current_state%qr_diff = 0.0_DEFAULT_PRECISION
      current_state%qi_diff = 0.0_DEFAULT_PRECISION
      current_state%qs_diff = 0.0_DEFAULT_PRECISION
      current_state%qg_diff = 0.0_DEFAULT_PRECISION
      current_state%qAitkenSolMass_diff = 0.0_DEFAULT_PRECISION
      current_state%qAccumSolMass_diff = 0.0_DEFAULT_PRECISION
      current_state%qAccumInsolMass_diff = 0.0_DEFAULT_PRECISION
      current_state%qCoarseSolMass_diff = 0.0_DEFAULT_PRECISION
      current_state%qCoarseDustMass_diff = 0.0_DEFAULT_PRECISION
      current_state%nl_diff = 0.0_DEFAULT_PRECISION
      current_state%nr_diff = 0.0_DEFAULT_PRECISION
      current_state%ni_diff = 0.0_DEFAULT_PRECISION
      current_state%ns_diff = 0.0_DEFAULT_PRECISION
      current_state%ng_diff = 0.0_DEFAULT_PRECISION
      current_state%nAitkenSolNumber_diff = 0.0_DEFAULT_PRECISION
      current_state%nAccumSolNumber_diff = 0.0_DEFAULT_PRECISION
      current_state%nAccumInsolNumber_diff = 0.0_DEFAULT_PRECISION
      current_state%nCoarseSolNumber_diff = 0.0_DEFAULT_PRECISION
      current_state%nCoarseDustnumber_diff = 0.0_DEFAULT_PRECISION

      current_state%qv_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%ql_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%qr_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%qi_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%qs_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%qg_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%qAitkenSolMass_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%qAccumSolMass_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%qAccumInsolMass_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%qCoarseSolMass_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%qCoarseDustMass_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%nl_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%nr_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%ni_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%ns_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%ng_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%nAitkenSolNumber_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%nAccumSolNumber_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%nAccumInsolNumber_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%nCoarseSolNumber_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%nCoarseDustnumber_buoyancy = 0.0_DEFAULT_PRECISION

      current_state%qv_tendency = 0.0_DEFAULT_PRECISION
      current_state%ql_tendency = 0.0_DEFAULT_PRECISION
      current_state%qr_tendency = 0.0_DEFAULT_PRECISION
      current_state%qi_tendency = 0.0_DEFAULT_PRECISION
      current_state%qs_tendency = 0.0_DEFAULT_PRECISION
      current_state%qg_tendency = 0.0_DEFAULT_PRECISION
      current_state%qAitkenSolMass_tendency = 0.0_DEFAULT_PRECISION
      current_state%qAccumSolMass_tendency = 0.0_DEFAULT_PRECISION
      current_state%qAccumInsolMass_tendency = 0.0_DEFAULT_PRECISION
      current_state%qCoarseSolMass_tendency = 0.0_DEFAULT_PRECISION
      current_state%qCoarseDustMass_tendency = 0.0_DEFAULT_PRECISION
      current_state%nl_tendency = 0.0_DEFAULT_PRECISION
      current_state%nr_tendency = 0.0_DEFAULT_PRECISION
      current_state%ni_tendency = 0.0_DEFAULT_PRECISION
      current_state%ns_tendency = 0.0_DEFAULT_PRECISION
      current_state%ng_tendency = 0.0_DEFAULT_PRECISION
      current_state%nAitkenSolNumber_tendency = 0.0_DEFAULT_PRECISION
      current_state%nAccumSolNumber_tendency = 0.0_DEFAULT_PRECISION
      current_state%nAccumInsolNumber_tendency = 0.0_DEFAULT_PRECISION
      current_state%nCoarseSolNumber_tendency = 0.0_DEFAULT_PRECISION
      current_state%nCoarseDustnumber_tendency = 0.0_DEFAULT_PRECISION

      current_state%uw_advection = 0.0_DEFAULT_PRECISION
      current_state%vw_advection = 0.0_DEFAULT_PRECISION
      current_state%uw_viscosity = 0.0_DEFAULT_PRECISION
      current_state%vw_viscosity = 0.0_DEFAULT_PRECISION
      current_state%uw_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%vw_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%uw_tendency = 0.0_DEFAULT_PRECISION
      current_state%vw_tendency = 0.0_DEFAULT_PRECISION
      current_state%uw_w = 0.0_DEFAULT_PRECISION
      current_state%vw_w = 0.0_DEFAULT_PRECISION

      current_state%tu_su = 0.0_DEFAULT_PRECISION
      current_state%uu_advection = 0.0_DEFAULT_PRECISION
      current_state%uu_viscosity = 0.0_DEFAULT_PRECISION
      current_state%wu_u = 0.0_DEFAULT_PRECISION
      current_state%tv_sv = 0.0_DEFAULT_PRECISION
      current_state%vv_advection = 0.0_DEFAULT_PRECISION
      current_state%vv_viscosity = 0.0_DEFAULT_PRECISION
      current_state%wv_v = 0.0_DEFAULT_PRECISION
      current_state%tw_sw = 0.0_DEFAULT_PRECISION
      current_state%ww_advection = 0.0_DEFAULT_PRECISION
      current_state%ww_viscosity = 0.0_DEFAULT_PRECISION
      current_state%ww_buoyancy = 0.0_DEFAULT_PRECISION

      current_state%u_thetal = 0.0_DEFAULT_PRECISION
      current_state%us_thetal = 0.0_DEFAULT_PRECISION
      current_state%u_thetal_advection = 0.0_DEFAULT_PRECISION
      current_state%u_thetal_viscosity_diffusion = 0.0_DEFAULT_PRECISION
      current_state%wu_thetal = 0.0_DEFAULT_PRECISION
      current_state%v_thetal = 0.0_DEFAULT_PRECISION
      current_state%vs_thetal = 0.0_DEFAULT_PRECISION
      current_state%v_thetal_advection = 0.0_DEFAULT_PRECISION
      current_state%v_thetal_viscosity_diffusion = 0.0_DEFAULT_PRECISION
      current_state%wv_thetal = 0.0_DEFAULT_PRECISION
      current_state%w_thetal = 0.0_DEFAULT_PRECISION
      current_state%ws_thetal = 0.0_DEFAULT_PRECISION
      current_state%w_thetal_advection = 0.0_DEFAULT_PRECISION
      current_state%w_thetal_viscosity_diffusion = 0.0_DEFAULT_PRECISION
      current_state%w_thetal_buoyancy = 0.0_DEFAULT_PRECISION
      current_state%ww_thetal = 0.0_DEFAULT_PRECISION
      current_state%thetal_thetal = 0.0_DEFAULT_PRECISION
      current_state%sthetal_thetal = 0.0_DEFAULT_PRECISION
      current_state%thetal_thetal_advection = 0.0_DEFAULT_PRECISION
      current_state%thetal_thetal_diffusion = 0.0_DEFAULT_PRECISION
      current_state%wthetal_thetal = 0.0_DEFAULT_PRECISION

      current_state%shprd = 0.0_DEFAULT_PRECISION
      current_state%w_ke = 0.0_DEFAULT_PRECISION
      current_state%w_p = 0.0_DEFAULT_PRECISION
      current_state%sub_buoy = 0.0_DEFAULT_PRECISION
      current_state%tke_tend = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%conditional_diagnostics_column_enabled .eqv. .true.) then
      current_state%w_zn_cln = 0.0_DEFAULT_PRECISION
      current_state%w_zn2_cln = 0.0_DEFAULT_PRECISION
      current_state%tmp_th_cln = 0.0_DEFAULT_PRECISION
      current_state%wth_cln = 0.0_DEFAULT_PRECISION
      current_state%th_pr_cln = 0.0_DEFAULT_PRECISION
      current_state%wth_pr_cln = 0.0_DEFAULT_PRECISION
      current_state%thv_pr_cln = 0.0_DEFAULT_PRECISION
      current_state%wthv_pr_cln = 0.0_DEFAULT_PRECISION
      current_state%th_pr2_cln = 0.0_DEFAULT_PRECISION
      current_state%wthsg_cln = 0.0_DEFAULT_PRECISION
      current_state%w_zn3_cln = 0.0_DEFAULT_PRECISION
      current_state%relhum_cln = 0.0_DEFAULT_PRECISION
      current_state%tmp_u_cln = 0.0_DEFAULT_PRECISION
      current_state%tmp_v_cln = 0.0_DEFAULT_PRECISION
      current_state%wu_cln = 0.0_DEFAULT_PRECISION
      current_state%wv_cln = 0.0_DEFAULT_PRECISION
      current_state%wusg_cln = 0.0_DEFAULT_PRECISION
      current_state%wvsg_cln = 0.0_DEFAULT_PRECISION
      current_state%TdegK_cln = 0.0_DEFAULT_PRECISION
      current_state%th_h_cln = 0.0_DEFAULT_PRECISION
      current_state%th_h_pr1_cln = 0.0_DEFAULT_PRECISION
      current_state%th_h_pr2_cln = 0.0_DEFAULT_PRECISION
      current_state%qvli_cln = 0.0_DEFAULT_PRECISION
      current_state%qvli_pr_cln = 0.0_DEFAULT_PRECISION
      current_state%qvli_pr2_cln = 0.0_DEFAULT_PRECISION
      current_state%qppt_cln = 0.0_DEFAULT_PRECISION
      current_state%qppt_pr_cln = 0.0_DEFAULT_PRECISION
      current_state%qppt_pr2_cln = 0.0_DEFAULT_PRECISION
      current_state%wqvli_pr_cln = 0.0_DEFAULT_PRECISION
      current_state%wqppt_pr_cln = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%pdf_analysis_enabled .eqv. .true.) then
      current_state%global_grid%configuration%vertical%w_dwn = 0.0_DEFAULT_PRECISION
      current_state%global_grid%configuration%vertical%w_up = 0.0_DEFAULT_PRECISION
      current_state%w_histogram_profile_local = 0.0_DEFAULT_PRECISION
    end if


    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(current_state%uwsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vwsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uusg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vvsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wwsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tkesg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wtsg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th2sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%sed_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ssub_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dissipation_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%buoysg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wkesg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%theta_dis_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vis_coef_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%diff_coef_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%richardson_number_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%richardson_squared_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wqv_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wql_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%wqr_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqi_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqs_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqg_sg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
    end if

    if (current_state%casim_profile_dgs_enabled .eqv. .true.) then
      call mpi_recv(current_state%dqc_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqg_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqi_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqr_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqs_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqv_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dth_mphys_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dth_cond_evap_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqv_cond_evap_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%phomc_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pinuc_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pidep_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psdep_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%piacw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psacw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psacr_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pisub_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pssub_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pimlt_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psmlt_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psaut_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psaci_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%praut_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pracw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%prevp_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pgacw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pgacs_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pgmlt_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pgsub_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psedi_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pseds_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psedr_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psedg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%psedl_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%pcond_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%profile_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(current_state%u_wind_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uprime_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%v_wind_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vprime_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wke_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ww_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%www_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wwww_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%theta_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_wind_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%rh_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wtheta_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wtheta_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vw_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uv_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th2_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%thinit, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qv_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ql_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%qr_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%qi_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%qs_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%qg_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%wqv_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wql_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%wqr_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqi_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqs_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqg_cn_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%wqv_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wql_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%wqr_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqi_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqs_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%wqg_ad_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%cloud_mask_tot, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%cloud_liq_mask_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%cloud_ice_mask_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
    end if

    if (current_state%forcing_enabled .eqv. .true.) then
      call mpi_recv(current_state%du_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dv_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dtheta_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqv_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dql_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqr_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqi_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqs_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%dqg_subs_profile_diag, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_u_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qi_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qr_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qs_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qg_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%diffusion_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_th_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_pr_tot_qi_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qr_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qs_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qg_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_pr_tot_tabs_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%buoyancy_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_w_buoy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%coriolis_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_u_corio, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_corio, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%pstep_enabled .eqv. .true.) then
      call mpi_recv(current_state%tendp_pr_tot_u_pt, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tendp_pr_tot_v_pt, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tendp_pr_tot_w_pt, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%pw_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_u_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_w_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qi_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qr_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qs_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qg_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%stepfields_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_th_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_pr_tot_qi_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qr_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qs_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qg_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_pr_tot_tabs_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%th_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_th_thad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_thad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_u_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_w_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_qv_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_ql_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_pr_tot_qi_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qr_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qs_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_pr_tot_qg_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_pr_tot_tabs_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%viscosity_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_pr_tot_u_visc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_v_visc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_w_visc, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%socrates_enabled .eqv. .true.) then
      call mpi_recv(current_state%cloud_reff_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%longwave_hr_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%shortwave_hr_tot, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_lw, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_lw, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_sw, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_sw, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_th_total, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_pr_tot_tabs_total, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%flux_budget_enabled .eqv. .true.) then
      call mpi_recv(current_state%th_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%qv_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ql_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qr_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qi_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qs_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qg_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAitkenSolMass_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumSolMass_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumInsolMass_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseSolMass_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseDustMass_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nl_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nr_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ni_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ns_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ng_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAitkenSolNumber_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumSolNumber_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumInsolNumber_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseSolNumber_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseDustnumber_flux_values, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%qv_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ql_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qr_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qi_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qs_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qg_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAitkenSolMass_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumSolMass_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumInsolMass_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseSolMass_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseDustMass_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nl_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nr_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ni_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ns_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ng_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAitkenSolNumber_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumSolNumber_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumInsolNumber_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseSolNumber_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseDustnumber_gradient, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%qv_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ql_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qr_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qi_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qs_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qg_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAitkenSolMass_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumSolMass_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumInsolMass_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseSolMass_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseDustMass_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nl_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nr_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ni_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ns_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ng_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAitkenSolNumber_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumSolNumber_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumInsolNumber_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseSolNumber_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseDustnumber_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%qv_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ql_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qr_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qi_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qs_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qg_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAitkenSolMass_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumSolMass_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumInsolMass_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseSolMass_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseDustMass_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nl_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nr_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ni_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ns_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ng_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAitkenSolNumber_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumSolNumber_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumInsolNumber_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseSolNumber_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseDustnumber_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%qv_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ql_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qr_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qi_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qs_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qg_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAitkenSolMass_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumSolMass_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumInsolMass_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseSolMass_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseDustMass_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nl_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nr_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ni_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ns_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ng_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAitkenSolNumber_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumSolNumber_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumInsolNumber_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseSolNumber_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseDustnumber_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%uw_advection, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vw_advection, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uw_viscosity, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vw_viscosity, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uw_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vw_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uw_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vw_tendency, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uw_w, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vw_w, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%tu_su, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uu_advection, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%uu_viscosity, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wu_u, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tv_sv, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vv_advection, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vv_viscosity, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wv_v, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tw_sw, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ww_advection, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ww_viscosity, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ww_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%u_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%us_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%u_thetal_advection, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%u_thetal_viscosity_diffusion, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wu_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%v_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vs_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%v_thetal_advection, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%v_thetal_viscosity_diffusion, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wv_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ws_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_thetal_advection, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_thetal_viscosity_diffusion, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_thetal_buoyancy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ww_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%thetal_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%sthetal_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%thetal_thetal_advection, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%thetal_thetal_diffusion, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wthetal_thetal, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%shprd, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_ke, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_p, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%sub_buoy, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tke_tend, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%conditional_diagnostics_column_enabled .eqv. .true.) then
      call mpi_recv(current_state%w_zn_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_zn2_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tmp_th_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wth_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_pr_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wth_pr_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%thv_pr_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wthv_pr_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_pr2_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wthsg_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_zn3_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%relhum_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tmp_u_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tmp_v_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wu_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wv_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wusg_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wvsg_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%TdegK_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_h_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_h_pr1_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%th_h_pr2_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qvli_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qvli_pr_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qvli_pr2_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qppt_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qppt_pr_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qppt_pr2_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wqvli_pr_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wqppt_pr_cln, global_array_size*2, MPI_REAL, 1, 1000, &
                                                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%pdf_analysis_enabled .eqv. .true.) then
      call mpi_recv(current_state%global_grid%configuration%vertical%w_dwn, &
                    global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%global_grid%configuration%vertical%w_up, &
                    global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%w_histogram_profile_local, global_array_size*current_state%n_w_bins*2, &
                    MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ls1 = len_trim(diagnostic_path)
      ls2 = 0
      do i = 1,ls1
          if(diagnostic_path(i:i).ne.' ') then
            ls2 = ls2 + 1
          endif
      enddo
    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = diagnostic_path(:ls2)//"/full_diag_1d_instantaneous_time_"//trim(&
                                                                          conv_to_string(int(current_state%time)))//".nc"
    else
      unique_filename = diagnostic_path(:ls2)//"/full_diag_1d_instantaneous_timestep_"//trim(&
                                                                            conv_to_string(current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "t", 1, time_id))
    call check(nf90_def_dim(ncdf_id, "z_size", global_grid_z_size, z_size_id))
    if (current_state%pdf_analysis_enabled .eqv. .true.) then
      call check(nf90_def_dim(ncdf_id, "bin_size", current_state%n_w_bins, bin_size_id))
    end if
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))


    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "zh", "heights at w levels", vertical_grid%z, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "zhn", "heights at pressure levels", vertical_grid%zn, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "dzh", "height of cells", vertical_grid%dz, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "rho", "cells density", vertical_grid%rho, "kg.m-3", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                          "rhon", "cells density (pressure levels)", vertical_grid%rho, "kg.m-3", current_state%time)

    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uwsg_mean", "uwsg_mean", &
                                          current_state%uwsg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vwsg_tot", "vwsg_mean", &
                                          current_state%vwsg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uusg_mean", "uusg_mean", &
                                          current_state%uusg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vvsg_tot", "vvsg_mean", &
                                          current_state%vvsg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wwsg_mean", "wwsg_mean", &
                                          current_state%wwsg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tkesg_mean", "tkesg_mean", &
                                          current_state%tkesg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wtsg_mean", "wtsg_mean", &
                                          current_state%wtsg_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th2sg_mean", "th2sg_mean", &
                                          current_state%th2sg_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "sed_mean", "mean subgrid turbulent transport", &
                                          current_state%sed_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                         "ssub_mean", "mean subgrid shear stress", &
                                          current_state%ssub_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dissipation_mean", "dissipation", &
                                          current_state%dissipation_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "buoysg_mean", "mean subgrid buoyant production", &
                                          current_state%buoysg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wkesg_mean", "wkesg_mean", &
                                          current_state%wkesg_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "theta_dis_mean", "theta_dis_mean", current_state%theta_dis_tot, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vis_coef_tot", "viscosity_coef_mean", &
                                          current_state%vis_coef_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "diff_coef_mean", "mean diffusion coefficient",&
                                          current_state%diff_coef_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                   "richardson_number_mean", "mean Richardson number", &
                                    current_state%richardson_number_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "richardson_squared_mean", "mean squared Richardson number", &
                                      current_state%richardson_squared_tot/surface, "K2", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqv_sg_mean", "wqv_sg_mean", &
                                            current_state%wqv_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wql_sg_mean", "wql_sg_mean", &
                                            current_state%wql_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqr_sg_mean", "wqr_sg_mean", &
                                            current_state%wqr_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqi_sg_mean", "wqi_sg_mean", &
                                            current_state%wqi_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqs_sg_mean", "wqs_sg_mean", &
                                            current_state%wqs_sg_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqg_sg_mean", "wqg_sg_mean", &
                                            current_state%wqg_sg_tot/surface, "K.m.s-1", current_state%time)
      end if
    end if

    if (current_state%casim_profile_dgs_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqc_mphys_mean", "mean qv rate", &
                                          current_state%dqc_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqg_mphys_mean", "mean qg rate", &
                                          current_state%dqg_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqi_mphys_mean", "mean qi rate", &
                                          current_state%dqi_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqr_mphys_mean", "mean qr rate", &
                                          current_state%dqr_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqs_mphys_mean", "mean qs rate", &
                                          current_state%dqs_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqv_mphys_mean", " mean qv rate", &
                                          current_state%dqv_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dth_mphys_mean", "mean th rate", &
                                          current_state%dth_mphys_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dth_cond_evap_mean", "mean condensation/evaporation rate due to th", &
                                          current_state%dth_cond_evap_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dqv_cond_evap_mean", "mean condensation/evaporation rate due to qv", &
                                          current_state%dqv_cond_evap_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "phomc_mean", "mean mass homogeneous freezing of cloud droplet rate", &
                                          current_state%phomc_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pinuc_mean", "mean mass ice nucleation rate", &
                                          current_state%pinuc_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pidep_mean", "mean mass ice deposition rate", &
                                          current_state%pidep_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psdep_mean", "mean mass snow deposition rate", &
                                          current_state%psdep_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "piacw_mean", "mean mass ice -> cloud -> ice accretion rate", &
                                          current_state%piacw_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psacw_mean", "mean mass snow -> cloud -> snow accretion rate", &
                                          current_state%psacw_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                            "psacr_mean", "mean mass snow -> rain -> graupel and snow -> rain -> snow accretion rate", &
                                          current_state%psacr_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pisub_mean", "mean mass ice sublimation rate", &
                                          current_state%pisub_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pssub_mean", "mean mass snow sublimation rate", &
                                          current_state%pssub_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pimlt_mean", "mean mass ice melting rate", &
                                          current_state%pimlt_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psmlt_mean", "mean mass snow melting rate", &
                                          current_state%psmlt_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psaut_mean", "mean mass autoconversion to snow rate", &
                                          current_state%psaut_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psaci_mean", "mean mass snow -> ice -> snow accretion rate", &
                                          current_state%psaci_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "praut_mean", "mean mass autoconversion to rain rate", &
                                          current_state%praut_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pracw_mean", "mean mass rain accreting cloud rate", &
                                          current_state%pracw_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "prevp_mean", "mean mass evaporation of rain rate", &
                                          current_state%prevp_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pgacw_mean", "mean mass graupel -> cloud -> graupel accretion rate", &
                                          current_state%pgacw_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pgacs_mean", "mean mass graupel -> snow accretion rate", &
                                          current_state%pgacs_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pgmlt_mean", "mean mass graupel melting rate", &
                                          current_state%pgmlt_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pgsub_mean", "mean mass graupel sublimation rate", &
                                          current_state%pgsub_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psedi_mean", "mean mass ice sedimentation rate", &
                                          current_state%psedi_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pseds_mean", "mean mass snow sedimentation rate", &
                                          current_state%pseds_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psedr_mean", "mean mass rain sedimentation rate", &
                                          current_state%psedr_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psedg_mean", "mean mass graupel sedimentation rate", &
                                          current_state%psedg_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "psedl_mean", "mean mass liquid/cloud sedimentation rate", &
                                          current_state%psedl_tot/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "pcond_mean", "mean mass condensation rate", &
                                          current_state%pcond_tot/surface, "kg.kg-1.s-1", current_state%time)
    end if

    if (current_state%profile_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "u_wind_mean", "u_wind_mean", &
                                          current_state%u_wind_tot/surface, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uu_mean", "uu_mean", &
                                          current_state%uprime_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "v_wind_tot", "v_wind_mean", &
                                          current_state%v_wind_tot/surface, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vprime_tot", "vv_mean", &
                                          current_state%vprime_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wke_mean", "wke_mean", &
                                          current_state%wke_tot/surface, "J.m-2.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ww_mean", "ww_mean", &
                                          current_state%ww_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "www_mean", "www_mean", &
                                          current_state%www_tot/surface, "m3.s-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wwww_mean", "wwww_mean", &
                                          current_state%wwww_tot/surface, "m4.s-4", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "theta_mean", "theta_mean", &
                                          current_state%theta_tot/surface, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "w_wind_mean", "w_wind_mean", &
                                          current_state%w_wind_tot/surface, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "rh_mean", "mean Relative Humidity", &
                                          current_state%rh_tot/surface, "%", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wtheta_ad_mean", "wtheta_ad_mean", &
                                          current_state%wtheta_ad_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wtheta_cn_mean", "wtheta_cn_mean", &
                                          current_state%wtheta_cn_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uw_mean", "uw_mean", &
                                          current_state%uw_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vw_tot", "vw_mean", &
                                          current_state%vw_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uv_mean", "uv_mean", &
                                          current_state%uv_tot/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th2_mean", "th2_mean", &
                                          current_state%th2_tot/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "thinit", "thinit", current_state%thinit, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qv_tot", "mean vapour mass", &
                                          current_state%qv_tot/surface, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ql_mean", "mean liquid mass", current_state%ql_tot/surface, &
                                          "kg.kg-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "qr_mean", "mean rain mass", current_state%qr_tot/surface, &
                                            "kg.kg-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "qi_mean", "mean ice mass", &
                                            current_state%qi_tot/surface, "kg.kg-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "qs_mean", "snow_mmr_mean", current_state%qs_tot/surface, &
                                            "kg.kg-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "qg_mean", "graupel_mmr_mean", current_state%qg_tot/surface, &
                                            "kg.kg-1", current_state%time)
      end if
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wqv_cn_mean", "wqv_cn_mean", &
                                          current_state%wqv_cn_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wql_cn_mean", "wql_cn_mean", &
                                          current_state%wql_cn_tot/surface, "K.m.s-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqr_cn_mean", "wqr_cn_mean", &
                                            current_state%wqr_cn_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqi_cn_mean", "wqi_cn_mean", &
                                            current_state%wqi_cn_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqs_cn_mean", "wqs_cn_mean", &
                                            current_state%wqs_cn_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqg_cn_mean", "wqg_cn_mean", &
                                            current_state%wqg_cn_tot/surface, "K.m.s-1", current_state%time)
      end if
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wqv_ad_mean", "wqv_ad_mean", &
                                          current_state%wqv_ad_tot/surface, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wql_ad_mean", "wql_ad_mean", &
                                          current_state%wql_ad_tot/surface, "K.m.s-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqr_ad_mean", "wqr_ad_mean", &
                                            current_state%wqr_ad_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqi_ad_mean", "wqi_ad_mean", &
                                            current_state%wqi_ad_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqs_ad_mean", "wqs_ad_mean", &
                                            current_state%wqs_ad_tot/surface, "K.m.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "wqg_ad_mean", "wqg_ad_mean", &
                                            current_state%wqg_ad_tot/surface, "K.m.s-1", current_state%time)

        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "cloud_mask_tot", "total_cloud_fraction", &
                                            current_state%cloud_mask_tot/surface, "no unit", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "cloud_liq_mask_mean", "liquid_cloud_fraction", &
                                            current_state%cloud_liq_mask_tot/surface, "no unit", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                            "cloud_ice_mask_mean", "ice_cloud_fraction", &
                                            current_state%cloud_ice_mask_tot/surface, "no unit", current_state%time)
      end if
    end if

    if (current_state%forcing_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "du_subs_profile_diag", "u_subsidence_mean", current_state%du_subs_profile_diag, &
                                          "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "dv_subs_profile_diag", "v_subsidence_mean", current_state%dv_subs_profile_diag, &
                                          "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                            "dtheta_subs_profile_diag", "th_subsidence_mean", current_state%dtheta_subs_profile_diag, &
                            "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqv_subs_profile_diag", "vapour_mmr_subsidence_mean", &
                                          current_state%dqv_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dql_subs_profile_diag", "liquid_mmr_subsidence_mean", &
                                          current_state%dql_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqr_subs_profile_diag", "mean qr rate (subsidence)", &
                                    current_state%dqr_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqi_subs_profile_diag", "mean qi rate (subsidence)",&
                                    current_state%dqi_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqs_subs_profile_diag", "mean qs rate (subsidence)", &
                                      current_state%dqs_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "dqg_subs_profile_diag", "dqg_subs_profile_diag", &
                                    current_state%dqg_subs_profile_diag, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_forc_mean", "mean tendancy u from forcing", &
                                      current_state%tend_pr_tot_u_forc/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_forc_mean", "mean tendancy v from forcing", &
                                      current_state%tend_pr_tot_v_forc/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_forc_mean", "mean tendancy th from forcing", &
                                      current_state%tend_pr_tot_th_forc/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_forc_mean", "mean tendancy qv from forcing", &
                                      current_state%tend_pr_tot_qv_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_forc_mean", "mean tendancy ql from forcing", &
                                      current_state%tend_pr_tot_ql_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qi_forc_mean", "mean tendancy qi from forcing", &
                                      current_state%tend_pr_tot_qi_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qr_forc_mean", "mean tendancy qr from forcing", &
                                      current_state%tend_pr_tot_qr_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qs_forc_mean", "mean tendancy qs from forcing", &
                                      current_state%tend_pr_tot_qs_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qg_forc_mean", "mean tendancy qg from forcing", &
                                      current_state%tend_pr_tot_qg_forc/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "tend_tabs_forc_mean", "mean tendancy tabs from forcing", &
                                    current_state%tend_pr_tot_tabs_forc/surface, "???", current_state%time)
    end if

    if (current_state%diffusion_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_th_diff_mean", "mean tendancy th from diffusion", &
                                          current_state%tend_pr_tot_th_diff/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_diff_mean", "mean tendancy qv from diffusion", &
                                      current_state%tend_pr_tot_qv_diff/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_diff_mean", "mean tendancy ql from diffusion", &
                                      current_state%tend_pr_tot_ql_diff/surface, "kg.kg-1.s-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qi_diff_mean", "mean tendancy qi from diffusion", &
                                        current_state%tend_pr_tot_qi_diff/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qr_diff_mean", "mean tendancy qr from diffusion", &
                                        current_state%tend_pr_tot_qr_diff/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qs_diff_mean", "mean tendancy qs from diffusion", &
                                        current_state%tend_pr_tot_qs_diff/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qg_diff_mean", "mean tendancy qg from diffusion", &
                                        current_state%tend_pr_tot_qg_diff/surface, "kg.kg-1.s-1", current_state%time)
      end if
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_diff_mean", "mean tendancy tabs from diffusion", &
                                      current_state%tend_pr_tot_tabs_diff/surface, "???", current_state%time)
    end if

    if (current_state%buoyancy_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_pr_tot_w_buoy", "mean tendancy w from bouyancy", &
                                      current_state%tend_pr_tot_w_buoy/surface, "m.s-2", current_state%time)
    end if

    if (current_state%coriolis_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_corio_mean", "mean tendancy u from coriolis", &
                                      current_state%tend_pr_tot_u_corio/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_corio_mean", "mean tendancy v from coriolis", &
                                      current_state%tend_pr_tot_v_corio/surface, "m.s-2", current_state%time)
    end if

    if (current_state%pstep_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_pt_mean", "mean tendancy u from pstep", &
                                      current_state%tendp_pr_tot_u_pt/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_pt_mean", "mean tendancy v from pstep", &
                                      current_state%tendp_pr_tot_v_pt/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_w_pt_mean", "mean tendancy w from pstep", &
                                      current_state%tendp_pr_tot_w_pt/surface, "m.s-2", current_state%time)
    end if

    if (current_state%pw_advection_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_pwad_mean", "mean tendancy u from pwadvection", &
                                      current_state%tend_pr_tot_u_pwad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_pwad_mean", "mean tendancy v from pwadvection", &
                                      current_state%tend_pr_tot_v_pwad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_w_pwad_mean", "mean tendancy w from pwadvection", &
                                      current_state%tend_pr_tot_w_pwad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_pwad_mean", "mean tendancy th from pwadvection", &
                                      current_state%tend_pr_tot_th_pwad/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_pwad_mean", "mean tendancy qv from pwadvection", &
                                      current_state%tend_pr_tot_qv_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_pwad_mean", "mean tendancy ql from pwadvection", &
                                      current_state%tend_pr_tot_ql_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qi_pwad_mean", "mean tendancy qi from pwadvection", &
                                      current_state%tend_pr_tot_qi_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qr_pwad_mean", "mean tendancy qr from pwadvection", &
                                      current_state%tend_pr_tot_qr_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qs_pwad_mean", "mean tendancy qs from pwadvection", &
                                      current_state%tend_pr_tot_qs_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qg_pwad_mean", "mean tendancy qg from pwadvection", &
                                      current_state%tend_pr_tot_qg_pwad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_pwad_mean", "mean tendancy tabs from pwadvection", &
                                      current_state%tend_pr_tot_tabs_pwad/surface, "???", current_state%time)
    end if

    if (current_state%stepfields_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_sf_mean", "mean tendancy th from stepfield", &
                                      current_state%tend_pr_tot_th_sf/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_sf_mean", "mean tendancy qv from stepfield", &
                                      current_state%tend_pr_tot_qv_sf/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_sf_mean", "mean tendancy ql from stepfield", &
                                      current_state%tend_pr_tot_ql_sf/surface, "kg.kg-1.s-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qi_sf_mean", "mean tendancy qi from stepfield", &
                                        current_state%tend_pr_tot_qi_sf/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qr_sf_mean", "mean tendancy qr from stepfield", &
                                        current_state%tend_pr_tot_qr_sf/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qs_sf_mean", "mean tendancy qs from stepfield", &
                                        current_state%tend_pr_tot_qs_sf/surface, "kg.kg-1.s-1", current_state%time)
        call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "tend_qg_sf_mean", "mean tendancy qg from stepfield", &
                                        current_state%tend_pr_tot_qg_sf/surface, "kg.kg-1.s-1", current_state%time)
      end if
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_sf_mean", "mean tendancy tabs from stepfield", &
                                      current_state%tend_pr_tot_tabs_sf/surface, "???", current_state%time)
    end if

    if (current_state%th_advection_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_thad_mean", "mean tendancy th from thadvection", &
                                      current_state%tend_pr_tot_th_thad/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_thad_mean", "mean tendancy tabs from thadvection", &
                                      current_state%tend_pr_tot_tabs_thad/surface, "???", current_state%time)
    end if

    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_u_tvad_mean", "mean tendancy u from tvdadvection", &
                                      current_state%tend_pr_tot_u_tvad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_v_tvad_mean", "mean tendancy v from tvdadvection", &
                                      current_state%tend_pr_tot_v_tvad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_w_tvad_mean", "mean tendancy w from tvdadvection", &
                                      current_state%tend_pr_tot_w_tvad/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_th_tvad_mean", "mean tendancy th from tvdadvection", &
                                      current_state%tend_pr_tot_th_tvad/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qv_tvad_mean", "mean tendancy qv from tvdadvection", &
                                      current_state%tend_pr_tot_qv_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_ql_tvad_mean", "mean tendancy ql from tvdadvection", &
                                      current_state%tend_pr_tot_ql_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qi_tvad_mean", "mean tendancy qi from tvdadvection", &
                                      current_state%tend_pr_tot_qi_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qr_tvad_mean", "mean tendancy qr from tvdadvection", &
                                      current_state%tend_pr_tot_qr_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qs_tvad_mean", "mean tendancy qs from tvdadvection", &
                                      current_state%tend_pr_tot_qs_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_qg_tvad_mean", "mean tendancy qg from tvdadvection", &
                                      current_state%tend_pr_tot_qs_tvad/surface, "kg.kg-1.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "tend_tabs_tvad_mean", "mean tendancy tabs from tvdadvection", &
                                      current_state%tend_pr_tot_tabs_tvad/surface, "???", current_state%time)
    end if

    if (current_state%viscosity_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_u_visc_mean", "mean tendancy u from viscosity", &
                                          current_state%tend_pr_tot_u_visc/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_v_visc_mean", "mean tendancy v from viscosity", &
                                          current_state%tend_pr_tot_v_visc/surface, "m.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_w_visc_mean", "mean tendancy w from viscosity", &
                                          current_state%tend_pr_tot_w_visc/surface, "m.s-2", current_state%time)
    end if

    if (current_state%socrates_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "cloud_reff_mean", "mean cloud reflectivity", &
                                          current_state%cloud_reff_tot/surface, " ", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "longwave_hr_mean", "mean lw heating rate", &
                                          current_state%longwave_hr_tot/surface, "K.d-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "shortwave_hr_mean", "mean sw heating rate", &
                                          current_state%shortwave_hr_tot/surface, "K.d-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_th_lw", "mean tendency th from lw heating rate", &
                                          current_state%tend_pr_tot_th_lw/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_tabs_lw", "mean tendency tabs from lw heating rate", &
                                          current_state%tend_pr_tot_tabs_lw/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_th_sw", "mean tendency th from sw heating rate", &
                                          current_state%tend_pr_tot_th_sw/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_tabs_sw", "mean tendency tabs from sw heating rate", &
                                          current_state%tend_pr_tot_tabs_sw/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_th_total", "mean tendency th from total heating rate", &
                                          current_state%tend_pr_tot_th_total/surface, "K.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tend_pr_mean_tabs_total", "mean tendency th from tabs heating rate", &
                                          current_state%tend_pr_tot_tabs_total/surface, "K.s-1", current_state%time)
    end if

    if (current_state%flux_budget_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_flux_values", "transport of kinematic heat flux", &
                                          current_state%th_flux_values/surface, "m2.s-2.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_gradient", "kinematic heat flux tendency", &
                                          current_state%th_gradient/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_diff", "dissipation kinematic heat flux", &
                                          current_state%th_diff/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_buoyancy", "advection kinematic heat flux", &
                                          current_state%th_buoyancy/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_tendency", "buoyancy kinematic heat flux", &
                                          current_state%th_tendency/surface, "m.s-1.K", current_state%time)

      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qv_flux_values", "transport of kinematic water vapor mass flux", &
                                          current_state%qv_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ql_flux_values", "transport of kinematic cloud water mass flux", &
                                          current_state%ql_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qr_flux_values", "transport of kinematic rain mass flux", &
                                          current_state%qr_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qi_flux_values", "transport of kinematic ice mass flux", &
                                          current_state%qi_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qs_flux_values", "transport of kinematic snow mass flux", &
                                          current_state%qs_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qg_flux_values", "transport of kinematic graupel mass flux", &
                                          current_state%qg_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAitkenSolMass_flux_values", "transport of kinematic Aitken mode mass flux", &
                                current_state%qAitkenSolMass_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumSolMass_flux_values", "transport of kinematic soluble Accumulation mode mass flux", &
                                current_state%qAccumSolMass_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumInsolMass_flux_values", "transport of kinematic insoluble Accumulation mode mass flux", &
                                current_state%qAccumInsolMass_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseSolMass_flux_values", "transport of kinematic soluble Coarse mode mass flux", &
                                current_state%qCoarseSolMass_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseDustMass_flux_values", "transport of kinematic dust Coarse mode mass flux", &
                                current_state%qCoarseDustMass_flux_values/surface, "m2.s-2.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nl_flux_values", "transport of kinematic cloud water number concentration flux", &
                                          current_state%nl_flux_values/surface, "m2.s-2.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nr_flux_values", "transport of kinematic rain number concentration flux", &
                                          current_state%nr_flux_values/surface, "m2.s-2.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ni_flux_values", "transport of kinematic ice number concentration flux", &
                                          current_state%ni_flux_values/surface, "m2.s-2.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ns_flux_values", "transport of kinematic snow number concentration flux", &
                                          current_state%ns_flux_values/surface, "m2.s-2.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ng_flux_values", "transport of kinematic graupel number concentration flux", &
                                          current_state%ng_flux_values/surface, "m2.s-2.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                              "nAitkenSolNumber_flux_values", "transport of kinematic Aitken mode number concentration flux", &
                              current_state%nAitkenSolNumber_flux_values/surface, "m2.s-2.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                      "nAccumSolNumber_flux_values", "transport of kinematic soluble Accumulation mode number concentration flux", &
                      current_state%nAccumSolNumber_flux_values/surface, "m2.s-2.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                  "nAccumInsolNumber_flux_values", "transport of kinematic insoluble Accumulation mode number concentration flux", &
                  current_state%nAccumInsolNumber_flux_values/surface, "m2.s-2.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseSolNumber_flux_values", "transport of kinematic soluble Coarse mode number concentration flux", &
                        current_state%nCoarseSolNumber_flux_values/surface, "m2.s-2.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseDustnumber_flux_values", "transport of kinematic dust Coarse mode number concentration flux", &
                        current_state%nCoarseDustnumber_flux_values/surface, "m2.s-2.kg-1", current_state%time)

      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qv_gradient", "advection kinematic water vapor mass flux", &
                                          current_state%qv_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ql_gradient", "advection kinematic cloud water mass flux", &
                                          current_state%ql_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qr_gradient", "advection kinematic rain mass flux", &
                                          current_state%qr_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qi_gradient", "advection kinematic ice mass flux", &
                                          current_state%qi_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qs_gradient", "advection kinematic snow mass flux", &
                                          current_state%qs_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qg_gradient", "advection kinematic graupel mass flux", &
                                          current_state%qg_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAitkenSolMass_gradient", "advection kinematic Aitken mode mass flux", &
                                current_state%qAitkenSolMass_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumSolMass_gradient", "advection kinematic soluble Accumulation mode mass flux", &
                                current_state%qAccumSolMass_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumInsolMass_gradient", "advection kinematic insoluble Accumulation mode mass flux", &
                                current_state%qAccumInsolMass_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseSolMass_gradient", "advection kinematic soluble Coarse mode mass flux", &
                                current_state%qCoarseSolMass_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseDustMass_gradient", "advection kinematic dust Coarse mode mass flux", &
                                current_state%qCoarseDustMass_gradient/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nl_gradient", "advection kinematic cloud water number concentration flux", &
                                          current_state%nl_gradient/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nr_gradient", "advection kinematic rain number concentration flux", &
                                          current_state%nr_gradient/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ni_gradient", "advection kinematic ice number concentration flux", &
                                          current_state%ni_gradient/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ns_gradient", "advection kinematic snow number concentration flux", &
                                          current_state%ns_gradient/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ng_gradient", "advection kinematic graupel number concentration flux", &
                                          current_state%ng_gradient/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                              "nAitkenSolNumber_gradient", "advection kinematic Aitken mode number concentration flux", &
                              current_state%nAitkenSolNumber_gradient/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                      "nAccumSolNumber_gradient", "advection kinematic soluble Accumulation mode number concentration flux", &
                      current_state%nAccumSolNumber_gradient/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                  "nAccumInsolNumber_gradient", "advection kinematic insoluble Accumulation mode number concentration flux", &
                  current_state%nAccumInsolNumber_gradient/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseSolNumber_gradient", "advection kinematic soluble Coarse mode number concentration flux", &
                        current_state%nCoarseSolNumber_gradient/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseDustnumber_gradient", "advection kinematic dust Coarse mode number concentration flux", &
                        current_state%nCoarseDustnumber_gradient/surface, "m.s-1.kg-1", current_state%time)

      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qv_diff", "dissipation kinematic water vapor mass flux", &
                                          current_state%qv_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ql_diff", "dissipation kinematic cloud water mass flux", &
                                          current_state%ql_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qr_diff", "dissipation kinematic rain mass flux", &
                                          current_state%qr_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qi_diff", "dissipation kinematic ice mass flux", &
                                          current_state%qi_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qs_diff", "dissipation kinematic snow mass flux", &
                                          current_state%qs_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qg_diff", "dissipation kinematic graupel mass flux", &
                                          current_state%qg_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAitkenSolMass_diff", "dissipation kinematic Aitken mode mass flux", &
                                current_state%qAitkenSolMass_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumSolMass_diff", "dissipation kinematic soluble Accumulation mode mass flux", &
                                current_state%qAccumSolMass_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumInsolMass_diff", "dissipation kinematic insoluble Accumulation mode mass flux", &
                                current_state%qAccumInsolMass_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseSolMass_diff", "dissipation kinematic soluble Coarse mode mass flux", &
                                current_state%qCoarseSolMass_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseDustMass_diff", "dissipation kinematic dust Coarse mode mass flux", &
                                current_state%qCoarseDustMass_diff/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nl_diff", "dissipation kinematic cloud water number concentration flux", &
                                          current_state%nl_diff/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nr_diff", "dissipation kinematic rain number concentration flux", &
                                          current_state%nr_diff/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ni_diff", "dissipation kinematic ice number concentration flux", &
                                          current_state%ni_diff/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ns_diff", "dissipation kinematic snow number concentration flux", &
                                          current_state%ns_diff/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ng_diff", "dissipation kinematic graupel number concentration flux", &
                                          current_state%ng_diff/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                              "nAitkenSolNumber_diff", "dissipation kinematic Aitken mode number concentration flux", &
                              current_state%nAitkenSolNumber_diff/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                      "nAccumSolNumber_diff", "dissipation kinematic soluble Accumulation mode number concentration flux", &
                      current_state%nAccumSolNumber_diff/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                  "nAccumInsolNumber_diff", "dissipation kinematic insoluble Accumulation mode number concentration flux", &
                  current_state%nAccumInsolNumber_diff/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseSolNumber_diff", "dissipation kinematic soluble Coarse mode number concentration flux", &
                        current_state%nCoarseSolNumber_diff/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseDustnumber_diff", "dissipation kinematic dust Coarse mode number concentration flux", &
                        current_state%nCoarseDustnumber_diff/surface, "m.s-1.kg-1", current_state%time)

      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qv_buoyancy", "buoyancy kinematic water vapor mass flux", &
                                          current_state%qv_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ql_buoyancy", "buoyancy kinematic cloud water mass flux", &
                                          current_state%ql_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qr_buoyancy", "buoyancy kinematic rain mass flux", &
                                          current_state%qr_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qi_buoyancy", "buoyancy kinematic ice mass flux", &
                                          current_state%qi_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qs_buoyancy", "buoyancy kinematic snow mass flux", &
                                          current_state%qs_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qg_buoyancy", "buoyancy kinematic graupel mass flux", &
                                          current_state%qg_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAitkenSolMass_buoyancy", "buoyancy kinematic Aitken mode mass flux", &
                                current_state%qAitkenSolMass_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumSolMass_buoyancy", "buoyancy kinematic soluble Accumulation mode mass flux", &
                                current_state%qAccumSolMass_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumInsolMass_buoyancy", "buoyancy kinematic insoluble Accumulation mode mass flux", &
                                current_state%qAccumInsolMass_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseSolMass_buoyancy", "buoyancy kinematic soluble Coarse mode mass flux", &
                                current_state%qCoarseSolMass_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseDustMass_buoyancy", "buoyancy kinematic dust Coarse mode mass flux", &
                                current_state%qCoarseDustMass_buoyancy/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nl_buoyancy", "buoyancy kinematic cloud water number concentration flux", &
                                          current_state%nl_buoyancy/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nr_buoyancy", "buoyancy kinematic rain number concentration flux", &
                                          current_state%nr_buoyancy/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ni_buoyancy", "buoyancy kinematic ice number concentration flux", &
                                          current_state%ni_buoyancy/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ns_buoyancy", "buoyancy kinematic snow number concentration flux", &
                                          current_state%ns_buoyancy/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ng_buoyancy", "buoyancy kinematic graupel number concentration flux", &
                                          current_state%ng_buoyancy/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                              "nAitkenSolNumber_buoyancy", "buoyancy kinematic Aitken mode number concentration flux", &
                              current_state%nAitkenSolNumber_buoyancy/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                      "nAccumSolNumber_buoyancy", "buoyancy kinematic soluble Accumulation mode number concentration flux", &
                      current_state%nAccumSolNumber_buoyancy/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                  "nAccumInsolNumber_buoyancy", "buoyancy kinematic insoluble Accumulation mode number concentration flux", &
                  current_state%nAccumInsolNumber_buoyancy/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseSolNumber_buoyancy", "buoyancy kinematic soluble Coarse mode number concentration flux", &
                        current_state%nCoarseSolNumber_buoyancy/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseDustnumber_buoyancy", "buoyancy kinematic dust Coarse mode number concentration flux", &
                        current_state%nCoarseDustnumber_buoyancy/surface, "m.s-1.kg-1", current_state%time)

      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qv_tendency", "kinematic water vapor mass flux tendency", &
                                          current_state%qv_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ql_tendency", "kinematic cloud water mass flux tendency", &
                                          current_state%ql_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qr_tendency", "kinematic rain mass flux tendency", &
                                          current_state%qr_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qi_tendency", "kinematic ice mass flux tendency", &
                                          current_state%qi_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qs_tendency", "kinematic snow mass flux tendency", &
                                          current_state%qs_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qg_tendency", "kinematic graupel mass flux tendency", &
                                          current_state%qg_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAitkenSolMass_tendency", "kinematic Aitken mode mass flux tendency", &
                                current_state%qAitkenSolMass_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumSolMass_tendency", "kinematic soluble Accumulation mode mass flux tendency", &
                                current_state%qAccumSolMass_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qAccumInsolMass_tendency", "kinematic insoluble Accumulation mode mass flux tendency", &
                                current_state%qAccumInsolMass_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseSolMass_tendency", "kinematic soluble Coarse mode mass flux tendency", &
                                current_state%qCoarseSolMass_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                "qCoarseDustMass_tendency", "kinematic dust Coarse mode mass flux tendency", &
                                current_state%qCoarseDustMass_tendency/surface, "m.s-1.kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nl_tendency", "kinematic cloud water number concentration flux tendency", &
                                          current_state%nl_tendency/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "nr_tendency", "kinematic rain number concentration flux tendency", &
                                          current_state%nr_tendency/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ni_tendency", "kinematic ice number concentration flux tendency", &
                                          current_state%ni_tendency/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ns_tendency", "kinematic snow number concentration flux tendency", &
                                          current_state%ns_tendency/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ng_tendency", "kinematic graupel number concentration flux tendency", &
                                          current_state%ng_tendency/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                              "nAitkenSolNumber_tendency", "kinematic Aitken mode number concentration flux tendency", &
                              current_state%nAitkenSolNumber_tendency/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                      "nAccumSolNumber_tendency", "kinematic soluble Accumulation mode number concentration flux tendency", &
                      current_state%nAccumSolNumber_tendency/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                  "nAccumInsolNumber_tendency", "kinematic insoluble Accumulation mode number concentration flux tendency", &
                  current_state%nAccumInsolNumber_tendency/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseSolNumber_tendency", "kinematic soluble Coarse mode number concentration flux tendency", &
                        current_state%nCoarseSolNumber_tendency/surface, "m.s-1.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                        "nCoarseDustnumber_tendency", "kinematic dust Coarse mode number concentration flux tendency", &
                        current_state%nCoarseDustnumber_tendency/surface, "m.s-1.kg-1", current_state%time)

      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uw_advection", "advection kinematic momentum flux uw", &
                                          current_state%uw_advection/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vw_advection", "advection kinematic momentum flux vw", &
                                          current_state%vw_advection/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uw_viscosity", "viscosity kinematic momentum flux uw", &
                                          current_state%uw_viscosity/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vw_viscosity", "viscosity kinematic momentum flux vw", &
                                          current_state%vw_viscosity/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uw_buoyancy", "buoyancy kinematic momentum flux uw", &
                                          current_state%uw_buoyancy/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vw_buoyancy", "buoyancy kinematic momentum flux vw", &
                                          current_state%vw_buoyancy/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uw_tendency", "kinematic momentum flux tendency uw", &
                                          current_state%uw_tendency/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vw_tendency", "kinematic momentum flux tendency vw", &
                                          current_state%vw_tendency/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "uw_w", "transport (along w) of kinematic momentum flux uw due to vertical speed", &
                                      current_state%uw_w/surface, "m3.s-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "vw_w", "transport (along w) of kinematic momentum flux vw due to vertical speed", &
                                      current_state%vw_w/surface, "m3.s-3", current_state%time)

      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tu_su", "kinematic momentum su flux tendency", &
                                          current_state%tu_su/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uu_advection", "advection kinematic momentum flux uu", &
                                          current_state%uu_advection/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "uu_viscosity", "viscosity kinematic momentum flux uu", &
                                          current_state%uu_viscosity/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wu_u", "transport (along u) of kinematic momentum flux wu due to horizontal speed", &
                                          current_state%wu_u/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tv_sv", "kinematic momentum sv flux tendency", &
                                          current_state%tv_sv/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vv_advection", "advection kinematic momentum flux vv", &
                                          current_state%vv_advection/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vv_viscosity", "viscosity kinematic momentum flux vv", &
                                          current_state%vv_viscosity/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wv_v", "transport (along v) of kinematic momentum flux wv due to horizontal speed", &
                                          current_state%wv_v/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tw_sw", "kinematic momentum sw flux tendency", &
                                          current_state%tw_sw/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ww_advection", "advection kinematic momentum flux ww", &
                                          current_state%ww_advection/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ww_viscosity", "viscosity kinematic momentum flux ww", &
                                          current_state%ww_viscosity/surface, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ww_buoyancy", "buoyancy kinematic momentum flux ww", &
                                          current_state%ww_buoyancy/surface, "m2.s-2", current_state%time)

      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "u_thetal", "kinematic heat flux", &
                                          current_state%u_thetal/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "us_thetal", "kinematic heat flux tendency", &
                                          current_state%us_thetal/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "u_thetal_advection", "advection kinematic heat flux", &
                                          current_state%u_thetal_advection/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "u_thetal_viscosity_diffusion", "viscosity_diffusion kinematic heat flux", &
                                        current_state%u_thetal_viscosity_diffusion/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wu_thetal", "transport of kinematic momentum flux due to temperature", &
                                          current_state%wu_thetal/surface, "m2.s-2.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "v_thetal", "kinematic heat flux", &
                                          current_state%v_thetal/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "vs_thetal", "kinematic heat flux tendency", &
                                          current_state%vs_thetal/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "v_thetal_advection", "advection kinematic heat flux", &
                                          current_state%v_thetal_advection/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "v_thetal_viscosity_diffusion", "viscosity_diffusion kinematic heat flux", &
                                        current_state%v_thetal_viscosity_diffusion/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wv_thetal", "transport of kinematic momentum flux due to temperature", &
                                          current_state%wv_thetal/surface, "m2.s-2.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "w_thetal", "kinematic heat flux", &
                                          current_state%w_thetal/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ws_thetal", "kinematic heat flux tendency", &
                                          current_state%ws_thetal/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "w_thetal_advection", "advection kinematic heat flux", &
                                          current_state%w_thetal_advection/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "w_thetal_viscosity_diffusion", "viscosity_diffusion kinematic heat flux", &
                                        current_state%w_thetal_viscosity_diffusion/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "w_thetal_buoyancy", "buoyancy kinematic heat flux", &
                                        current_state%w_thetal_buoyancy/surface, "m.s-1.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "ww_thetal", "transport of kinematic momentum flux due to temperature", &
                                          current_state%ww_thetal/surface, "m2.s-2.K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "thetal_thetal", "kinematic heat flux", &
                                          current_state%thetal_thetal/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "sthetal_thetal", "kinematic heat flux tendency", &
                                          current_state%sthetal_thetal/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "thetal_thetal_advection", "advection kinematic heat flux", &
                                          current_state%thetal_thetal_advection/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "thetal_thetal_diffusion", "diffusion kinematic heat flux", &
                                        current_state%thetal_thetal_diffusion/surface, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wthetal_thetal", "transport of kinematic momentum flux due to temperature", &
                                          current_state%wthetal_thetal/surface, "m.s-1.K2", current_state%time)

      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "shprd", "shear production", &
                                          current_state%shprd/surface, "m2.s-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "w_ke", "resolved turbulent transport", &
                                          current_state%w_ke/surface, "m3.s-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "sub_buoy", "subgrid buoyant production", &
                                          current_state%sub_buoy/surface, "m2.s-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "w_p", "pressure transport", &
                                          current_state%w_p/surface, "m2.s-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tke_tend", "TKE tendency", &
                                          current_state%tke_tend/surface, "m2.s-2", current_state%time)
    end if

    if (current_state%conditional_diagnostics_column_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "w_zn_cln", "vertical velocity", &
                                          current_state%w_zn_cln, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "w_zn2_cln", "vertical velocity variance", &
                                          current_state%w_zn2_cln, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tmp_th_cln", "potential temperature", &
                                          current_state%tmp_th_cln, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wth_cln", "vertical flux of potential temperature", &
                                          current_state%wth_cln, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_pr_cln", "potential temperature anomalies", &
                                          current_state%th_pr_cln, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wth_pr_cln", "vertical flux of potential temperature anomalies", &
                                          current_state%wth_pr_cln, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "thv_pr_cln", "virtual potential temperature anomalies", &
                                          current_state%thv_pr_cln, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wthv_pr_cln", "vertical flux of virtual potential temperature", &
                                          current_state%wthv_pr_cln, "K.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_pr2_cln", "variance of potential temperature", &
                                          current_state%th_pr2_cln, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wthsg_cln", "flux diffusion coefficient by temperature", &
                                          current_state%wthsg_cln, "???", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "w_zn3_cln", "w_zn***3", &
                                          current_state%w_zn3_cln, "m3.s-3", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "relhum_cln", "relative humidity", &
                                          current_state%relhum_cln, "no unit", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tmp_u_cln", "zonal velocity", &
                                          current_state%tmp_u_cln, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "tmp_v_cln", "meridional velovity", &
                                          current_state%tmp_v_cln, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wu_cln", "flux components wu", &
                                          current_state%wu_cln, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wv_cln", "flux components wv", &
                                          current_state%wv_cln, "m2.s-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wusg_cln", "flux viscosity coefficient follonwing u", &
                                          current_state%wusg_cln, "???", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wvsg_cln", "flux viscosity coefficient follonwing v", &
                                          current_state%wvsg_cln, "???", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "TdegK_cln", "Degree temperature", &
                                          current_state%TdegK_cln, "°C", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_h_cln", "air temperature variation", &
                                          current_state%th_h_cln, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_h_pr1_cln", "potential temperature anomalies variation", &
                                          current_state%th_h_pr1_cln, "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "th_h_pr2_cln", "variance ofpotential temperature anomalies variation", &
                                          current_state%th_h_pr2_cln, "K2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qvli_cln", "???", &
                                          current_state%qvli_cln, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qvli_pr_cln", "??? anomalies", &
                                          current_state%qvli_pr_cln, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qvli_pr2_cln", "??? variance anomalies", &
                                          current_state%qvli_pr2_cln, "kg2.kg-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qppt_cln", "???", &
                                          current_state%qppt_cln, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qppt_pr_cln", "??? anomalies", &
                                          current_state%qppt_pr_cln, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "qppt_pr2_cln", "??? variance anomalies", &
                                          current_state%qppt_pr2_cln, "kg2.kg-2", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wqvli_pr_cln", "??? flux anomalies by w component", &
                                          current_state%wqvli_pr_cln, "kg.kg-1.m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                          "wqppt_pr_cln", "??? flux anomalies by w component", &
                                          current_state%wqppt_pr_cln, "kg.kg-1.m.s-1", current_state%time)
    end if

    if (current_state%pdf_analysis_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                      "w_dwn", "downward vertical wind speed", &
                      current_state%global_grid%configuration%vertical%w_dwn, "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                      "w_up", "upward vertical wind speed", &
                      current_state%global_grid%configuration%vertical%w_up, "m.s-1", current_state%time)
      call define_and_write_variable_1D_bin(ncdf_id, z_size_id, bin_size_id, &
                                          "w_histogram_profile_local", "vertical wind speed histogram", &
                                          current_state%w_histogram_profile_local, "no unit", current_state%time)
    end if

    call check(nf90_close(ncdf_id))

    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      deallocate(current_state%uwsg_tot, current_state%vwsg_tot, current_state%uusg_tot, current_state%vvsg_tot, &
              current_state%wwsg_tot, current_state%tkesg_tot, current_state%wtsg_tot, current_state%th2sg_tot, &
              current_state%sed_tot, current_state%ssub_tot, current_state%dissipation_tot, current_state%buoysg_tot, &
              current_state%wkesg_tot, current_state%theta_dis_tot, current_state%vis_coef_tot, &
              current_state%diff_coef_tot, current_state%richardson_number_tot, current_state%richardson_squared_tot, &
              current_state%wqv_sg_tot, current_state%wql_sg_tot)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        deallocate(current_state%wqr_sg_tot, current_state%wqi_sg_tot, current_state%wqs_sg_tot, current_state%wqg_sg_tot)
      end if
    end if

    if (current_state%casim_profile_dgs_enabled .eqv. .true.) then
      deallocate(current_state%dqc_mphys_tot, current_state%dqg_mphys_tot, current_state%dqi_mphys_tot, &
                current_state%dqr_mphys_tot, current_state%dqs_mphys_tot, current_state%dqv_mphys_tot, &
                current_state%dth_mphys_tot, current_state%dth_cond_evap_tot, current_state%dqv_cond_evap_tot, &
                current_state%phomc_tot, current_state%pinuc_tot, current_state%pidep_tot, current_state%psdep_tot, &
                current_state%piacw_tot, current_state%psacw_tot, current_state%psacr_tot, current_state%pisub_tot, &
                current_state%pssub_tot, current_state%pimlt_tot, current_state%psmlt_tot, current_state%psaut_tot, &
                current_state%psaci_tot, current_state%praut_tot, current_state%pracw_tot, current_state%prevp_tot, &
                current_state%pgacw_tot, current_state%pgacs_tot, current_state%pgmlt_tot, current_state%pgsub_tot, &
                current_state%psedi_tot, current_state%pseds_tot, current_state%psedr_tot, current_state%psedg_tot, &
                current_state%psedl_tot, current_state%pcond_tot)
    end if

    if (current_state%profile_diagnostics_enabled .eqv. .true.) then
      deallocate(current_state%u_wind_tot, current_state%uprime_tot, current_state%v_wind_tot, current_state%vprime_tot, &
                current_state%wke_tot, current_state%ww_tot,  current_state%www_tot, current_state%wwww_tot, &
                current_state%theta_tot, current_state%w_wind_tot, current_state%rh_tot, current_state%wtheta_ad_tot, &
                current_state%wtheta_cn_tot, current_state%uw_tot, current_state%vw_tot, current_state%uv_tot, &
                current_state%th2_tot, current_state%thinit, current_state%qv_tot, &
                current_state%ql_tot, current_state%wqv_cn_tot, current_state%wql_cn_tot, current_state%wqv_ad_tot, &
                current_state%wql_ad_tot)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        deallocate(current_state%qr_tot, current_state%qi_tot, current_state%qs_tot, current_state%qg_tot, &
                  current_state%wqr_cn_tot, current_state%wqi_cn_tot, current_state%wqs_cn_tot, &
                  current_state%wqg_cn_tot, current_state%wqr_ad_tot, current_state%wqi_ad_tot, &
                  current_state%wqs_ad_tot, current_state%wqg_ad_tot, current_state%cloud_mask_tot, &
                  current_state%cloud_liq_mask_tot, current_state%cloud_ice_mask_tot)
      end if
    end if

    if (current_state%forcing_enabled .eqv. .true.) then
      deallocate(current_state%du_subs_profile_diag, current_state%dv_subs_profile_diag, current_state%dtheta_subs_profile_diag, &
                current_state%dqv_subs_profile_diag, current_state%dql_subs_profile_diag, current_state%dqr_subs_profile_diag, &
                current_state%dqi_subs_profile_diag, current_state%dqs_subs_profile_diag, current_state%dqg_subs_profile_diag, &
                current_state%tend_pr_tot_u_forc, current_state%tend_pr_tot_v_forc, current_state%tend_pr_tot_th_forc, &
                current_state%tend_pr_tot_qv_forc, current_state%tend_pr_tot_ql_forc, current_state%tend_pr_tot_qi_forc, &
                current_state%tend_pr_tot_qr_forc, current_state%tend_pr_tot_qs_forc, current_state%tend_pr_tot_qg_forc, &
                current_state%tend_pr_tot_tabs_forc)
    end if

    if (current_state%diffusion_enabled .eqv. .true.) then
      deallocate(current_state%tend_pr_tot_th_diff, current_state%tend_pr_tot_qv_diff, current_state%tend_pr_tot_ql_diff, &
                current_state%tend_pr_tot_tabs_diff)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        deallocate(current_state%tend_pr_tot_qi_diff, current_state%tend_pr_tot_qr_diff, current_state%tend_pr_tot_qs_diff, &
                current_state%tend_pr_tot_qg_diff)
      end if
    end if

    if (current_state%buoyancy_enabled .eqv. .true.) then
      deallocate(current_state%tend_pr_tot_w_buoy)
    end if

    if (current_state%coriolis_enabled .eqv. .true.) then
      deallocate(current_state%tend_pr_tot_u_corio,  current_state%tend_pr_tot_v_corio)
    end if

    if (current_state%pstep_enabled .eqv. .true.) then
      deallocate(current_state%tendp_pr_tot_u_pt, current_state%tendp_pr_tot_v_pt, current_state%tendp_pr_tot_w_pt)
    end if

    if (current_state%pw_advection_enabled .eqv. .true.) then
      deallocate(current_state%tend_pr_tot_u_pwad, current_state%tend_pr_tot_v_pwad, current_state%tend_pr_tot_w_pwad, &
                current_state%tend_pr_tot_th_pwad, current_state%tend_pr_tot_qv_pwad, current_state%tend_pr_tot_ql_pwad, &
                current_state%tend_pr_tot_qi_pwad, current_state%tend_pr_tot_qr_pwad, current_state%tend_pr_tot_qs_pwad, &
                current_state%tend_pr_tot_qg_pwad, current_state%tend_pr_tot_tabs_pwad)
    end if

    if (current_state%stepfields_enabled .eqv. .true.) then
      deallocate(current_state%tend_pr_tot_th_sf, current_state%tend_pr_tot_qv_sf, current_state%tend_pr_tot_ql_sf, &
                current_state%tend_pr_tot_tabs_sf)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        deallocate(current_state%tend_pr_tot_qi_sf, current_state%tend_pr_tot_qr_sf, &
                current_state%tend_pr_tot_qs_sf, current_state%tend_pr_tot_qg_sf)
      end if
    end if

    if (current_state%th_advection_enabled .eqv. .true.) then
      deallocate(current_state%tend_pr_tot_th_thad, current_state%tend_pr_tot_tabs_thad)
    end if

    if (current_state%tvd_advection_enabled .eqv. .true.) then
      deallocate(current_state%tend_pr_tot_u_tvad, current_state%tend_pr_tot_v_tvad, current_state%tend_pr_tot_w_tvad, &
                current_state%tend_pr_tot_th_tvad, current_state%tend_pr_tot_qv_tvad, current_state%tend_pr_tot_ql_tvad, &
                current_state%tend_pr_tot_tabs_tvad)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        deallocate(current_state%tend_pr_tot_qi_tvad, current_state%tend_pr_tot_qr_tvad, current_state%tend_pr_tot_qs_tvad, &
                  current_state%tend_pr_tot_qg_tvad)
      end if
    end if

    if (current_state%viscosity_enabled .eqv. .true.) then
      deallocate(current_state%tend_pr_tot_u_visc, current_state%tend_pr_tot_v_visc, current_state%tend_pr_tot_w_visc)
    end if

    if (current_state%socrates_enabled .eqv. .true.) then
      deallocate(current_state%cloud_reff_tot, current_state%longwave_hr_tot, current_state%shortwave_hr_tot, &
                current_state%tend_pr_tot_th_lw, current_state%tend_pr_tot_tabs_lw, current_state%tend_pr_tot_th_sw, &
                current_state%tend_pr_tot_tabs_sw, current_state%tend_pr_tot_th_total, current_state%tend_pr_tot_tabs_total)
    end if

    if (current_state%flux_budget_enabled .eqv. .true.) then
      deallocate(current_state%th_flux_values, current_state%th_gradient, current_state%th_diff, &
                current_state%th_buoyancy, current_state%th_tendency, &
                current_state%qv_flux_values, current_state%ql_flux_values, current_state%qr_flux_values, &
                current_state%qi_flux_values, current_state%qs_flux_values, current_state%qg_flux_values, &
                current_state%qAitkenSolMass_flux_values, current_state%qAccumSolMass_flux_values, &
                current_state%qAccumInsolMass_flux_values, current_state%qCoarseSolMass_flux_values, &
                current_state%qCoarseDustMass_flux_values, current_state%nl_flux_values, &
                current_state%nr_flux_values, current_state%ni_flux_values, current_state%ns_flux_values, &
                current_state%ng_flux_values, current_state%nAitkenSolNumber_flux_values, &
                current_state%nAccumSolNumber_flux_values, current_state%nAccumInsolNumber_flux_values, &
                current_state%nCoarseSolNumber_flux_values, current_state%nCoarseDustnumber_flux_values, &
                current_state%qv_gradient, current_state%ql_gradient, current_state%qr_gradient, &
                current_state%qi_gradient, current_state%qs_gradient, current_state%qg_gradient, &
                current_state%qAitkenSolMass_gradient, current_state%qAccumSolMass_gradient, &
                current_state%qAccumInsolMass_gradient, current_state%qCoarseSolMass_gradient, &
                current_state%qCoarseDustMass_gradient, current_state%nl_gradient, &
                current_state%nr_gradient, current_state%ni_gradient, current_state%ns_gradient, &
                current_state%ng_gradient, current_state%nAitkenSolNumber_gradient, &
                current_state%nAccumSolNumber_gradient, current_state%nAccumInsolNumber_gradient, &
                current_state%nCoarseSolNumber_gradient, current_state%nCoarseDustnumber_gradient, &
                current_state%qv_diff, current_state%ql_diff, current_state%qr_diff, &
                current_state%qi_diff, current_state%qs_diff, current_state%qg_diff, &
                current_state%qAitkenSolMass_diff, current_state%qAccumSolMass_diff, &
                current_state%qAccumInsolMass_diff, current_state%qCoarseSolMass_diff, &
                current_state%qCoarseDustMass_diff, current_state%nl_diff, &
                current_state%nr_diff, current_state%ni_diff, current_state%ns_diff, &
                current_state%ng_diff, current_state%nAitkenSolNumber_diff, &
                current_state%nAccumSolNumber_diff, current_state%nAccumInsolNumber_diff, &
                current_state%nCoarseSolNumber_diff, current_state%nCoarseDustnumber_diff, &
                current_state%qv_buoyancy, current_state%ql_buoyancy, current_state%qr_buoyancy, &
                current_state%qi_buoyancy, current_state%qs_buoyancy, current_state%qg_buoyancy, &
                current_state%qAitkenSolMass_buoyancy, current_state%qAccumSolMass_buoyancy, &
                current_state%qAccumInsolMass_buoyancy, current_state%qCoarseSolMass_buoyancy, &
                current_state%qCoarseDustMass_buoyancy, current_state%nl_buoyancy, &
                current_state%nr_buoyancy, current_state%ni_buoyancy, current_state%ns_buoyancy, &
                current_state%ng_buoyancy, current_state%nAitkenSolNumber_buoyancy, &
                current_state%nAccumSolNumber_buoyancy, current_state%nAccumInsolNumber_buoyancy, &
                current_state%nCoarseSolNumber_buoyancy, current_state%nCoarseDustnumber_buoyancy, &
                current_state%qv_tendency, current_state%ql_tendency, current_state%qr_tendency, &
                current_state%qi_tendency, current_state%qs_tendency, current_state%qg_tendency, &
                current_state%qAitkenSolMass_tendency, current_state%qAccumSolMass_tendency, &
                current_state%qAccumInsolMass_tendency, current_state%qCoarseSolMass_tendency, &
                current_state%qCoarseDustMass_tendency, current_state%nl_tendency, &
                current_state%nr_tendency, current_state%ni_tendency, current_state%ns_tendency, &
                current_state%ng_tendency, current_state%nAitkenSolNumber_tendency, &
                current_state%nAccumSolNumber_tendency, current_state%nAccumInsolNumber_tendency, &
                current_state%nCoarseSolNumber_tendency, current_state%nCoarseDustnumber_tendency, &
                current_state%uw_advection, current_state%vw_advection, current_state%uw_viscosity, &
                current_state%vw_viscosity, current_state%uw_buoyancy, current_state%vw_buoyancy, &
                current_state%uw_tendency, current_state%vw_tendency, current_state%uw_w, &
                current_state%vw_w, current_state%tu_su, current_state%uu_advection, current_state%uu_viscosity, &
                current_state%wu_u, current_state%tv_sv, current_state%vv_advection, current_state%vv_viscosity, &
                current_state%wv_v, current_state%tw_sw, current_state%ww_advection, current_state%ww_viscosity, &
                current_state%ww_buoyancy, current_state%u_thetal, current_state%us_thetal, &
                current_state%u_thetal_advection, current_state%u_thetal_viscosity_diffusion, &
                current_state%wu_thetal, current_state%v_thetal, current_state%vs_thetal, &
                current_state%v_thetal_advection, current_state%v_thetal_viscosity_diffusion, &
                current_state%wv_thetal, current_state%w_thetal, current_state%ws_thetal, &
                current_state%w_thetal_advection, current_state%w_thetal_viscosity_diffusion, &
                current_state%w_thetal_buoyancy, current_state%ww_thetal, current_state%thetal_thetal, &
                current_state%sthetal_thetal, current_state%thetal_thetal_advection, &
                current_state%thetal_thetal_diffusion, current_state%wthetal_thetal, &
                current_state%shprd, current_state%w_ke, current_state%sub_buoy, &
                current_state%w_p, current_state%tke_tend)
    end if
    if (current_state%conditional_diagnostics_column_enabled .eqv. .true.) then
      deallocate(current_state%w_zn_cln, current_state%w_zn2_cln, current_state%tmp_th_cln, &
                current_state%wth_cln, current_state%th_pr_cln, current_state%wth_pr_cln, &
                current_state%thv_pr_cln, current_state%wthv_pr_cln, current_state%th_pr2_cln, &
                current_state%wthsg_cln, current_state%w_zn3_cln, current_state%relhum_cln, &
                current_state%tmp_u_cln, current_state%tmp_v_cln, current_state%wu_cln, &
                current_state%wv_cln, current_state%wusg_cln, current_state%wvsg_cln, &
                current_state%TdegK_cln, current_state%th_h_cln, current_state%th_h_pr1_cln, &
                current_state%th_h_pr2_cln, current_state%qvli_cln, current_state%qvli_pr_cln, &
                current_state%qvli_pr2_cln, current_state%qppt_cln, current_state%qppt_pr_cln, &
                current_state%qppt_pr2_cln, current_state%wqvli_pr_cln, current_state%wqppt_pr_cln)
    end if

    if (current_state%pdf_analysis_enabled .eqv. .true.) then
      deallocate(current_state%global_grid%configuration%vertical%w_dwn, &
                current_state%global_grid%configuration%vertical%w_up, current_state%w_histogram_profile_local)
    end if

    current_state%diag_1d_done = .true.
    call mpi_send(current_state%diag_1d_done, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, ierr)
    current_state%diag_1d_done = .false.

  end subroutine diagnostic_file_1d_generation




  subroutine diagnostic_file_2d_generation(current_state, vertical_grid, io_communicator_arg, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size, diagnostic_path)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           io_communicator_arg
    character(len=200), intent(in) :: diagnostic_path

    integer :: z_size_id, y_size_id, x_size_id!, scalar_size_id, time_id
    integer :: ncdf_id, ierr
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename

    global_grid_z_size_netcdf = global_grid_z_size
    global_grid_y_size_netcdf = global_grid_y_size
    global_grid_x_size_netcdf = global_grid_x_size

    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      allocate(current_state%qlmax(global_grid_y_size, global_grid_x_size))
      allocate(current_state%hqlmax(global_grid_y_size, global_grid_x_size))
      allocate(current_state%cltop(global_grid_y_size, global_grid_x_size))
      allocate(current_state%clbas(global_grid_y_size, global_grid_x_size))
      allocate(current_state%vwp(global_grid_y_size, global_grid_x_size))
      allocate(current_state%vwp_usr(global_grid_y_size, global_grid_x_size))
      allocate(current_state%lwp(global_grid_y_size, global_grid_x_size))
      allocate(current_state%lwp_usr(global_grid_y_size, global_grid_x_size))
      allocate(current_state%wmax(global_grid_y_size, global_grid_x_size))
      allocate(current_state%wmin(global_grid_y_size, global_grid_x_size))
      allocate(current_state%reske(global_grid_y_size, global_grid_x_size))
      allocate(current_state%senhf(global_grid_y_size, global_grid_x_size))
      allocate(current_state%lathf(global_grid_y_size, global_grid_x_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%rwp(global_grid_y_size, global_grid_x_size))
        allocate(current_state%rwp_usr(global_grid_y_size, global_grid_x_size))
        allocate(current_state%iwp(global_grid_y_size, global_grid_x_size))
        allocate(current_state%iwp_usr(global_grid_y_size, global_grid_x_size))
        allocate(current_state%swp(global_grid_y_size, global_grid_x_size))
        allocate(current_state%swp_usr(global_grid_y_size, global_grid_x_size))
        allocate(current_state%gwp(global_grid_y_size, global_grid_x_size))
        allocate(current_state%gwp_usr(global_grid_y_size, global_grid_x_size))
        allocate(current_state%tot_iwp(global_grid_y_size, global_grid_x_size))
        allocate(current_state%tot_iwp_usr(global_grid_y_size, global_grid_x_size))
      end if
    end if

    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      allocate(current_state%subke_2d(global_grid_y_size, global_grid_x_size))
    end if

    if (current_state%casim_enabled .eqv. .true.) then
      allocate(current_state%surface_precip(global_grid_y_size, global_grid_x_size))
      allocate(current_state%surface_cloudsed(global_grid_y_size, global_grid_x_size))
      allocate(current_state%surface_rainsed(global_grid_y_size, global_grid_x_size))
    end if

    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      current_state%qlmax = 0.0_DEFAULT_PRECISION
      current_state%hqlmax = 0.0_DEFAULT_PRECISION
      current_state%cltop = 0.0_DEFAULT_PRECISION
      current_state%clbas = 0.0_DEFAULT_PRECISION
      current_state%vwp = 0.0_DEFAULT_PRECISION
      current_state%vwp_usr = 0.0_DEFAULT_PRECISION
      current_state%lwp = 0.0_DEFAULT_PRECISION
      current_state%lwp_usr = 0.0_DEFAULT_PRECISION
      current_state%wmax = 0.0_DEFAULT_PRECISION
      current_state%wmin = 0.0_DEFAULT_PRECISION
      current_state%reske = 0.0_DEFAULT_PRECISION
      current_state%senhf = 0.0_DEFAULT_PRECISION
      current_state%lathf = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%rwp = 0.0_DEFAULT_PRECISION
        current_state%rwp_usr = 0.0_DEFAULT_PRECISION
        current_state%iwp = 0.0_DEFAULT_PRECISION
        current_state%iwp_usr = 0.0_DEFAULT_PRECISION
        current_state%swp = 0.0_DEFAULT_PRECISION
        current_state%swp_usr = 0.0_DEFAULT_PRECISION
        current_state%gwp = 0.0_DEFAULT_PRECISION
        current_state%gwp_usr = 0.0_DEFAULT_PRECISION
        current_state%tot_iwp = 0.0_DEFAULT_PRECISION
        current_state%tot_iwp_usr = 0.0_DEFAULT_PRECISION
      end if
    end if

    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      current_state%surface_precip = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%casim_enabled .eqv. .true.) then
      current_state%surface_cloudsed = 0.0_DEFAULT_PRECISION
      current_state%surface_rainsed = 0.0_DEFAULT_PRECISION
      current_state%subke_2d = 0.0_DEFAULT_PRECISION
    end if

    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(current_state%qlmax, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%hqlmax, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%cltop, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%clbas, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%vwp_usr, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%lwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%lwp_usr, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wmax, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%wmin, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%reske, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%senhf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%lathf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%rwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%rwp_usr, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%iwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%iwp_usr, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%swp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%swp_usr, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%gwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%gwp_usr, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tot_iwp, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tot_iwp_usr, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
    end if

    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call mpi_recv(current_state%subke_2d, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    if (current_state%casim_enabled .eqv. .true.) then
      call mpi_recv(current_state%surface_precip, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%surface_cloudsed, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%surface_rainsed, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ls1 = len_trim(diagnostic_path)
    ls2 = 0
    do i = 1,ls1
        if(diagnostic_path(i:i).ne.' ') then
          ls2 = ls2 + 1
        endif
    enddo
    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = diagnostic_path(:ls2)//"/full_diag_2d_instantaneous_time_"//trim(&
                                                                                conv_to_string(int(current_state%time)))//".nc"
    else
      unique_filename = diagnostic_path(:ls2)//"/full_diag_2d_instantaneous_timestep_"//trim(&
                                                                                  conv_to_string(current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "t", 1, time_id))
    call check(nf90_def_dim(ncdf_id, "z", global_grid_z_size, z_size_id))
    call check(nf90_def_dim(ncdf_id, "y", global_grid_y_size, y_size_id))
    call check(nf90_def_dim(ncdf_id, "x", global_grid_x_size, x_size_id))
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))


    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)

    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "zh", "heights at w levels", vertical_grid%z, "m", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "zhn", "heights at pressure levels", vertical_grid%zn, "m", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "dzh", "height of cells", vertical_grid%dz, "m", current_state%time)

      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "qlmax", "maximum liquid water content in a column", current_state%qlmax, &
                                        "kg.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                      "hqlmax", "height of the maximum liquid water content in a column", current_state%hqlmax, &
                      "m", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "cltop", "cloud top height", current_state%cltop, "m.s-1", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "clbas", "cloud base height", current_state%clbas, "m.s-1", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "vwp", "water vapour path for each column", current_state%vwp, "kg.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "vwp_usr", "user water vapour path for each column", current_state%vwp_usr, &
                                        "kg.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "lwp", "liquid water path for each column", current_state%lwp, "kg.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "lwp_usr", "user liquid water path for each column", current_state%lwp_usr, &
                                        "kg.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "wmax", "maximum vertical velocity for each column", current_state%wmax, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "wmin", "minimum vertical velocity for each column", current_state%wmin, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "reske", "resolved ke", current_state%reske, "J or kg.m2.s-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "senhf", "surface sensible heat flux", current_state%senhf, "W.m-2", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "lathf", "surface latent heat flux", current_state%lathf, "W.m-2", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "rwp", "rain water path for each column", current_state%rwp, "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "rwp_usr", "user rain water path for each column", current_state%rwp_usr, &
                                          "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "iwp", "ice water path for each column", current_state%iwp, "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "iwp_usr", "user ice water path for each column", current_state%iwp_usr, &
                                          "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "swp", "snow water path for each column", current_state%swp, "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "swp_usr", "user snow water path for each column", current_state%swp_usr, &
                                          "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "gwp", "graupel water path for each column", &
                                          current_state%gwp, "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                          "gwp_usr", "user graupel water path for each column", &
                                          current_state%gwp_usr, "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                        "tot_iwp", "total ice water path (iwp + swp + gwp) for each column", current_state%tot_iwp, &
                        "kg.m-2", current_state%time)
        call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                        "usr_tot_iwp", "user total ice water path (iwp + swp + gwp) for each column", current_state%tot_iwp_usr, &
                        "kg.m-2", current_state%time)
      end if
    end if

    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "subke_2d", "subke_2d", current_state%subke_2d, "J or kg.m2.s-2", current_state%time)
    end if

    if (current_state%casim_enabled .eqv. .true.) then
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "surface_precip", "surface_precip", current_state%surface_precip, &
                                        "???", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "surface_cloudsed", "surface_cloudsed", current_state%surface_cloudsed, &
                                        "???", current_state%time)
      call define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, &
                                        "surface_rainsed", "surface_rainsed", current_state%surface_rainsed, &
                                        "???", current_state%time)
    end if

    call check(nf90_close(ncdf_id))

    if (current_state%scalar_diagnostics_enabled .eqv. .true.) then
      deallocate(current_state%qlmax, current_state%hqlmax, current_state%cltop, current_state%clbas, current_state%vwp, &
              current_state%vwp_usr, current_state%lwp, current_state%lwp_usr, current_state%wmax, current_state%wmin, &
              current_state%reske, current_state%senhf, current_state%lathf)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        deallocate(current_state%rwp, current_state%rwp_usr, current_state%iwp, current_state%iwp_usr, current_state%swp, &
                  current_state%swp_usr, current_state%gwp, current_state%gwp_usr, current_state%tot_iwp, &
                  current_state%tot_iwp_usr)
      end if
    end if

    if (current_state%subgrid_profile_diagnostics_enabled .eqv. .true.) then
      deallocate(current_state%subke_2d)
    end if

    if (current_state%casim_enabled .eqv. .true.) then
      deallocate(current_state%surface_precip, current_state%surface_cloudsed, current_state%surface_rainsed,)
    end if

    current_state%diag_2d_done = .true.
    call mpi_send(current_state%diag_2d_done, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, ierr)
    current_state%diag_2d_done = .false.

    !print*,"send_data_for_diag_2D MONC"
  end subroutine diagnostic_file_2d_generation




  subroutine diagnostic_file_3d_generation(current_state, vertical_grid, io_communicator_arg, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size, diagnostic_path)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           io_communicator_arg
    character(len=200), intent(in) :: diagnostic_path

    integer :: z_size_id, y_size_id, x_size_id!, scalar_size_id, time_id
    integer :: ncdf_id, ierr
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename
    real(DEFAULT_PRECISION), dimension(global_grid_z_size, global_grid_y_size, global_grid_x_size) :: theta, temp_K

    global_grid_z_size_netcdf = global_grid_z_size
    global_grid_y_size_netcdf = global_grid_y_size
    global_grid_x_size_netcdf = global_grid_x_size

    allocate(current_state%u%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%v%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%w%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%th%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%p%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qv%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ql%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qi%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qs%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qg%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAitkenSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAccumSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAccumInsolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qCoarseSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qCoarseDustMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nl%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ni%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ns%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ng%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAitkenSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAccumSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAccumInsolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nCoarseSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nCoarseDustnumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    ! diag bouyancy
    if (current_state%buoyancy_enabled .eqv. .true.) then
      allocate(current_state%tend_3d_w_buoy(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if
    ! diag coriolis
    if (current_state%coriolis_enabled .eqv. .true.) then
      allocate(current_state%tend_3d_u_corio(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_v_corio(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if
    ! diag diffusion
    if (current_state%diffusion_enabled .eqv. .true.) then
      allocate(current_state%tend_3d_th_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qv_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_ql_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%tend_3d_qr_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
        allocate(current_state%tend_3d_qi_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
        allocate(current_state%tend_3d_qs_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
        allocate(current_state%tend_3d_qg_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      end if
      allocate(current_state%tend_3d_tabs_diff(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if
    ! diag forcing
    if (current_state%forcing_enabled .eqv. .true.) then
      allocate(current_state%tend_3d_u_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_v_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_th_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qv_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_ql_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qr_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qi_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qs_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qg_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_tabs_forc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if
    ! diag pstep
    if (current_state%pstep_enabled .eqv. .true.) then
      allocate(current_state%tendp_3d_u_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tendp_3d_v_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tendp_3d_w_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_u_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_v_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_w_pt(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if
    ! diag pwadvection
    if (current_state%pw_advection_enabled .eqv. .true.) then
      allocate(current_state%tend_3d_u_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_v_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_w_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_th_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qv_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_ql_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qr_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qi_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qs_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qg_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_tabs_pwad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if
    ! diag stepfields
    if (current_state%stepfields_enabled .eqv. .true.) then
      allocate(current_state%tend_3d_th_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qv_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_ql_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%tend_3d_qr_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
        allocate(current_state%tend_3d_qi_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
        allocate(current_state%tend_3d_qs_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
        allocate(current_state%tend_3d_qg_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      end if
      allocate(current_state%tend_3d_tabs_sf(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if
    ! diag thadvection
    if (current_state%th_advection_enabled .eqv. .true.) then
      allocate(current_state%tend_3d_th_thad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_tabs_thad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if
    ! diag tvadvection
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      allocate(current_state%tend_3d_u_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_v_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_w_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_th_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_qv_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_ql_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        allocate(current_state%tend_3d_qr_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
        allocate(current_state%tend_3d_qi_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
        allocate(current_state%tend_3d_qs_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
        allocate(current_state%tend_3d_qg_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      end if
      allocate(current_state%tend_3d_tabs_tvad(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if
    ! diag viscosity
    if (current_state%viscosity_enabled .eqv. .true.) then
      allocate(current_state%tend_3d_u_visc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_v_visc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_w_visc(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if

    allocate(current_state%rdAitkenSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdAccumSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdCoarseSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_cloud%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_rain%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_ice%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_snow%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_graupel%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%RH%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%RI%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))

    ! SOCRATES
    if (current_state%socrates_enabled .eqv. .true.) then
      allocate(current_state%cloud_reff%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%lwrad_hr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%swrad_hr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_tabs_lw(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_tabs_sw(global_grid_z_size, global_grid_y_size, global_grid_x_size))
      allocate(current_state%tend_3d_tabs_total(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    end if

    current_state%u%data = 0.0_DEFAULT_PRECISION
    current_state%v%data = 0.0_DEFAULT_PRECISION
    current_state%w%data = 0.0_DEFAULT_PRECISION
    current_state%th%data = 0.0_DEFAULT_PRECISION
    current_state%p%data = 0.0_DEFAULT_PRECISION
    current_state%qv%data = 0.0_DEFAULT_PRECISION
    current_state%ql%data = 0.0_DEFAULT_PRECISION
    current_state%qr%data = 0.0_DEFAULT_PRECISION
    current_state%qi%data = 0.0_DEFAULT_PRECISION
    current_state%qs%data = 0.0_DEFAULT_PRECISION
    current_state%qg%data = 0.0_DEFAULT_PRECISION
    current_state%qAitkenSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qAccumSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qAccumInsolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qCoarseSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qCoarseDustMass%data = 0.0_DEFAULT_PRECISION
    current_state%nl%data = 0.0_DEFAULT_PRECISION
    current_state%nr%data = 0.0_DEFAULT_PRECISION
    current_state%ni%data = 0.0_DEFAULT_PRECISION
    current_state%ns%data = 0.0_DEFAULT_PRECISION
    current_state%ng%data = 0.0_DEFAULT_PRECISION
    current_state%nAitkenSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nAccumSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nAccumInsolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nCoarseSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nCoarseDustnumber%data = 0.0_DEFAULT_PRECISION
    ! diag bouyancy
    if (current_state%buoyancy_enabled .eqv. .true.) then
      current_state%tend_3d_w_buoy = 0.0_DEFAULT_PRECISION
    end if
    ! diag coriolis
    if (current_state%coriolis_enabled .eqv. .true.) then
      current_state%tend_3d_u_corio = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_v_corio = 0.0_DEFAULT_PRECISION
    end if
    ! diag diffusion
    if (current_state%diffusion_enabled .eqv. .true.) then
      current_state%tend_3d_th_diff = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qv_diff = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_ql_diff = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%tend_3d_qr_diff = 0.0_DEFAULT_PRECISION
        current_state%tend_3d_qi_diff = 0.0_DEFAULT_PRECISION
        current_state%tend_3d_qs_diff = 0.0_DEFAULT_PRECISION
        current_state%tend_3d_qg_diff = 0.0_DEFAULT_PRECISION
      end if
      current_state%tend_3d_tabs_diff = 0.0_DEFAULT_PRECISION
    end if
    ! diag forcing
    if (current_state%forcing_enabled .eqv. .true.) then
      current_state%tend_3d_u_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_v_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_th_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qv_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_ql_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qr_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qi_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qs_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qg_forc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_tabs_forc = 0.0_DEFAULT_PRECISION
    end if
    ! diag pstep
    if (current_state%pstep_enabled .eqv. .true.) then
      current_state%tendp_3d_u_pt = 0.0_DEFAULT_PRECISION
      current_state%tendp_3d_v_pt = 0.0_DEFAULT_PRECISION
      current_state%tendp_3d_w_pt = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_u_pt = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_v_pt = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_w_pt = 0.0_DEFAULT_PRECISION
    end if
    ! diag pwadvection
    if (current_state%pw_advection_enabled .eqv. .true.) then
      current_state%tend_3d_u_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_v_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_w_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_th_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qv_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_ql_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qr_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qi_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qs_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qg_pwad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_tabs_pwad = 0.0_DEFAULT_PRECISION
    end if
    ! diag stepfields
    if (current_state%stepfields_enabled .eqv. .true.) then
      current_state%tend_3d_th_sf = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qv_sf = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_ql_sf = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%tend_3d_qr_sf = 0.0_DEFAULT_PRECISION
        current_state%tend_3d_qi_sf = 0.0_DEFAULT_PRECISION
        current_state%tend_3d_qs_sf = 0.0_DEFAULT_PRECISION
        current_state%tend_3d_qg_sf = 0.0_DEFAULT_PRECISION
      end if
      current_state%tend_3d_tabs_sf = 0.0_DEFAULT_PRECISION
    end if
    ! diag thadvection
    if (current_state%th_advection_enabled .eqv. .true.) then
      current_state%tend_3d_th_thad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_tabs_thad = 0.0_DEFAULT_PRECISION
    end if
    ! diag tvadvection
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      current_state%tend_3d_u_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_v_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_w_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_th_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_qv_tvad = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_ql_tvad = 0.0_DEFAULT_PRECISION
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        current_state%tend_3d_qr_tvad = 0.0_DEFAULT_PRECISION
        current_state%tend_3d_qi_tvad = 0.0_DEFAULT_PRECISION
        current_state%tend_3d_qs_tvad = 0.0_DEFAULT_PRECISION
        current_state%tend_3d_qg_tvad = 0.0_DEFAULT_PRECISION
      end if
      current_state%tend_3d_tabs_tvad = 0.0_DEFAULT_PRECISION
    end if
    ! diag viscosity
    if (current_state%viscosity_enabled .eqv. .true.) then
      current_state%tend_3d_u_visc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_v_visc = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_w_visc = 0.0_DEFAULT_PRECISION
    end if

    current_state%rdAitkenSol%data = 0.0_DEFAULT_PRECISION
    current_state%rdAccumSol%data = 0.0_DEFAULT_PRECISION
    current_state%rdCoarseSol%data = 0.0_DEFAULT_PRECISION
    current_state%D0_cloud%data = 0.0_DEFAULT_PRECISION
    current_state%D0_rain%data = 0.0_DEFAULT_PRECISION
    current_state%D0_ice%data = 0.0_DEFAULT_PRECISION
    current_state%D0_snow%data = 0.0_DEFAULT_PRECISION
    current_state%D0_graupel%data = 0.0_DEFAULT_PRECISION
    current_state%RH%data = 0.0_DEFAULT_PRECISION
    current_state%RI%data = 0.0_DEFAULT_PRECISION

    ! SOCRATES
    if (current_state%socrates_enabled .eqv. .true.) then
      current_state%cloud_reff%data = 0.0_DEFAULT_PRECISION
      current_state%lwrad_hr%data = 0.0_DEFAULT_PRECISION
      current_state%swrad_hr%data = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_tabs_lw = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_tabs_sw = 0.0_DEFAULT_PRECISION
      current_state%tend_3d_tabs_total = 0.0_DEFAULT_PRECISION
    end if

    call mpi_recv(current_state%u%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%v%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%w%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%th%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%p%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qv%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ql%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qi%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qs%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qg%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qAitkenSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qAccumSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qAccumInsolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qCoarseSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%qCoarseDustMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nl%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ni%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ns%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ng%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nAitkenSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nAccumSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nAccumInsolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nCoarseSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%nCoarseDustnumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    ! diag bouyancy
    if (current_state%buoyancy_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_w_buoy, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag coriolis
    if (current_state%coriolis_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_corio, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_corio, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag diffusion
    if (current_state%diffusion_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_th_diff, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_diff, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_diff, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_3d_qr_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qi_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qs_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qg_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_3d_tabs_diff, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag forcing
    if (current_state%forcing_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_th_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qr_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qi_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qs_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qg_forc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_forc, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag pstep
    if (current_state%pstep_enabled .eqv. .true.) then
      call mpi_recv(current_state%tendp_3d_u_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tendp_3d_v_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tendp_3d_w_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_u_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_w_pt, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag pwadvection
    if (current_state%pw_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_w_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_th_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qr_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qi_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qs_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qg_pwad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_pwad, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag stepfields
    if (current_state%stepfields_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_th_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_3d_qr_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qi_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qs_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qg_sf, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_3d_tabs_sf, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag thadvection
    if (current_state%th_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_th_thad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_thad, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag tvadvection
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_w_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_th_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_qv_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_ql_tvad, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call mpi_recv(current_state%tend_3d_qr_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qi_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qs_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(current_state%tend_3d_qg_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      call mpi_recv(current_state%tend_3d_tabs_tvad, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    ! diag viscosity
    if (current_state%viscosity_enabled .eqv. .true.) then
      call mpi_recv(current_state%tend_3d_u_visc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_v_visc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_w_visc, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    call mpi_recv(current_state%rdAitkenSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%rdAccumSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%rdCoarseSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_cloud%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_rain%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_ice%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_snow%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_graupel%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%RH%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%RI%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

    ! SOCRATES
    if (current_state%socrates_enabled .eqv. .true.) then
      call mpi_recv(current_state%cloud_reff%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%lwrad_hr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%swrad_hr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_lw, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_sw, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%tend_3d_tabs_total, global_array_size*2, MPI_REAL, &
                    1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ls1 = len_trim(diagnostic_path)
    ls2 = 0
    do i = 1,ls1
        if(diagnostic_path(i:i).ne.' ') then
          ls2 = ls2 + 1
        endif
    enddo
    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = diagnostic_path(:ls2)//"/full_diag_3d_instantaneous_time_"//trim(&
                                                                            conv_to_string(int(current_state%time)))//".nc"
    else
      unique_filename = diagnostic_path(:ls2)//"/full_diag_3d_instantaneous_timestep_"//trim(&
                                                                            conv_to_string(current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "t", 1, time_id))
    call check(nf90_def_dim(ncdf_id, "z", global_grid_z_size, z_size_id))
    call check(nf90_def_dim(ncdf_id, "y", global_grid_y_size, y_size_id))
    call check(nf90_def_dim(ncdf_id, "x", global_grid_x_size, x_size_id))
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))


    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "zh", "heights at w levels", vertical_grid%z, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "zhn", "heights at pressure levels", vertical_grid%zn, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "dzh", "height of cells", vertical_grid%dz, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "rho", "cells density", vertical_grid%rho, "kg.m-3", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                          "rhon", "cells density (pressure levels)", vertical_grid%rho, "kg.m-3", current_state%time)

    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_bottom", "x bottom coordinate", current_state%global_grid%bottom(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_top", "x top coordinate", current_state%global_grid%top(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_resolution", "resolution along x axis", current_state%global_grid%resolution(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_bottom", "y bottom coordinate", current_state%global_grid%bottom(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_top", "y top coordinate", current_state%global_grid%top(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_resolution", "resolution along y axis", current_state%global_grid%resolution(2), "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "u", "u wind component", current_state%u%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                     "v", "v wind component", current_state%v%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "w", "w wind component", current_state%w%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "thp", "potentiel temperature perturbations", current_state%th%data, "K", current_state%time)
    do iter_x = 1, global_grid_x_size
      do iter_y = 1,global_grid_y_size
        theta(:,iter_y, iter_x) = current_state%th%data(:,iter_y, iter_x) + vertical_grid%thref(:)
        temp_K(:,iter_y, iter_x) = theta(:,iter_y, iter_x) * vertical_grid%rprefrcp(:)
      end do
    end do
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "theta", "potentiel temperature", theta, "K", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "temp_K", "air temperature", temp_K, "K", current_state%time)

    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "p", "pressure", current_state%p%data, "Pa", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qv", "water vapour mass mixing ratio", current_state%qv%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ql", "liquid water mass mixing ratio", current_state%ql%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qr", "rain water mass mixing ratio", current_state%qr%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qi", "ice water mass mixing ratio", current_state%qi%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qs", "snow water mass mixing ratio", current_state%qs%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qg", "graupel water mass mixing ratio", current_state%qg%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAitkenSolMass", "Aitken soluble aerosol mass mixing ratio", &
                                      current_state%qAitkenSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAccumSolMass", "accumulation soluble aerosol mass mixing ratio",&
                                      current_state%qAccumSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAccumInsolMass", "accumulation insoluble aerosol mass mixing ratio",&
                                      current_state%qAccumInsolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qCoarseSolMass", "coarse soluble aerosol mixing ratio",&
                                      current_state%qCoarseSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qCoarseDustMass", "dust coarse aerosol mixing ratio",&
                                      current_state%qCoarseDustMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nl", "liquid water number concentration", current_state%nl%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nr", "rain water number concentration", current_state%nr%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ni", "ice water number concentration", current_state%ni%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ns", "snow water number concentration", current_state%ns%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ng", "graupel water number concentration", current_state%ng%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAitkenSolNumber", "Aitken soluble aerosol number concentration", &
                                      current_state%nAitkenSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAccumSolNumber", "accumulation soluble aerosol number concentration",&
                                      current_state%nAccumSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAccumInsolNumber", "accumulation insoluble aerosol number concentration",&
                                      current_state%nAccumInsolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nCoarseSolNumber", "coarse soluble aerosol number concentration",&
                                      current_state%nCoarseSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nCoarseDustnumber", "dust coarse aerosol number concentration",&
                                      current_state%nCoarseDustnumber%data, "kg-1", current_state%time)
    ! diag bouyancy
    if (current_state%buoyancy_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_buoy", "tend_3d_w_buoy", current_state%tend_3d_w_buoy, &
                                        "m.s-1", current_state%time)
    end if
    ! diag coriolis
    if (current_state%coriolis_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_corio", "3d tendency u from coriolis", current_state%tend_3d_u_corio, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_corio", "3d tendency v from coriolis", current_state%tend_3d_v_corio, &
                                        "m.s-1", current_state%time)
    end if
    ! diag diffusion
    if (current_state%diffusion_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_diff", "3d tendency th from diffusion", current_state%tend_3d_th_diff, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_diff", "3d tendency qv from diffusion", current_state%tend_3d_qv_diff, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_diff", "3d tendency ql from diffusion", current_state%tend_3d_ql_diff, &
                                        "kg.kg-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qr_diff", "3d tendency qr from diffusion", current_state%tend_3d_qr_diff, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qi_diff", "3d tendency qi from diffusion", current_state%tend_3d_qi_diff, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qs_diff", "3d tendency qs from diffusion", current_state%tend_3d_qs_diff, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qg_diff", "3d tendency qg from diffusion", current_state%tend_3d_qg_diff, &
                                          "kg.kg-1", current_state%time)
      end if
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_diff", "3d tendency tabs from diffusion", &
                                        current_state%tend_3d_tabs_diff, "???", current_state%time)
    end if
    ! diag forcing
    if (current_state%forcing_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_forc", "3d tendency u from forcing", current_state%tend_3d_u_forc, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_forc", "3d tendency v from forcing", current_state%tend_3d_v_forc, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_forc", "3d tendency th from forcing", current_state%tend_3d_th_forc, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_forc", "3d tendency qv from forcing", current_state%tend_3d_qv_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_forc", "3d tendency ql from forcing", current_state%tend_3d_ql_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qr_forc", "3d tendency qr from forcing", current_state%tend_3d_qr_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qi_forc", "3d tendency qi from forcing", current_state%tend_3d_qi_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qs_forc", "3d tendency qs from forcing", current_state%tend_3d_qs_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qg_forc", "3d tendency qg from forcing", current_state%tend_3d_qg_forc, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                "tend_3d_tabs_forc", "3d tendency tabs from forcing", current_state%tend_3d_tabs_forc, &
                                "???", current_state%time)
    end if
    ! diag pstep
    if (current_state%pstep_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tendp_3d_u_pt", "3d previous tendency u from pstep", current_state%tendp_3d_u_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tendp_3d_v_pt", "3d previous tendency v from pstep", current_state%tendp_3d_v_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tendp_3d_w_pt", "3d previous tendency w from pstep", current_state%tendp_3d_w_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_pt", "3d tendency u from pstep", current_state%tend_3d_u_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_pt", "3d tendency v from pstep", current_state%tend_3d_v_pt, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_pt", "3d tendency w from pstep", current_state%tend_3d_w_pt, &
                                        "m.s-1", current_state%time)
    end if
    ! diag pwadvection
    if (current_state%pw_advection_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_pwad", "3d tendency u from pwadvection", current_state%tend_3d_u_pwad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_pwad", "3d tendency v from pwadvection", current_state%tend_3d_v_pwad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_pwad", "3d tendency w from pwadvection", current_state%tend_3d_w_pwad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_pwad", "3d tendency th from pwadvection", current_state%tend_3d_th_pwad, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_pwad", "3d tendency qv from pwadvection", &
                                        current_state%tend_3d_qv_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_pwad", "3d tendency ql from pwadvection", &
                                        current_state%tend_3d_ql_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qr_pwad", "3d tendency qr from pwadvection", &
                                        current_state%tend_3d_qr_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qi_pwad", "3d tendency qi from pwadvection", &
                                        current_state%tend_3d_qi_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qs_pwad", "3d tendency qs from pwadvection", &
                                        current_state%tend_3d_qs_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qg_pwad", "3d tendency qg from pwadvection", &
                                        current_state%tend_3d_qg_pwad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_pwad", "3d tendency tabs from pwadvection", &
                                        current_state%tend_3d_tabs_pwad, "???", current_state%time)
    end if
    ! diag stepfields
    if (current_state%stepfields_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_sf", "3d tendency th from stepfield", current_state%tend_3d_th_sf, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_sf", "3d tendency qv from stepfield", current_state%tend_3d_qv_sf, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_sf", "3d tendency ql from stepfield", current_state%tend_3d_ql_sf, &
                                        "kg.kg-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qr_sf", "3d tendency qr from stepfield", current_state%tend_3d_qr_sf, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qi_sf", "3d tendency qi from stepfield", current_state%tend_3d_qi_sf, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qs_sf", "3d tendency qs from stepfield", current_state%tend_3d_qs_sf, &
                                          "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qg_sf", "3d tendency qg from stepfield", current_state%tend_3d_qg_sf, &
                                          "kg.kg-1", current_state%time)
      end if
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_sf", "3d tendency tabs from stepfield", current_state%tend_3d_tabs_sf, &
                                        "???", current_state%time)
    end if
    ! diag thadvection
    if (current_state%th_advection_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_thad", "3d tendency th from thadvection", current_state%tend_3d_th_thad, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_thad", "3d tendency th tabs thadvection", &
                                        current_state%tend_3d_tabs_thad, "???", current_state%time)
    end if
    ! diag tvadvection
    if (current_state%tvd_advection_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_tvad", "3d tendency u from tvdadvection", current_state%tend_3d_u_tvad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_tvad", "3d tendency v from tvdadvection", current_state%tend_3d_v_tvad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_tvad", "3d tendency w from tvdadvection", current_state%tend_3d_w_tvad, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_th_tvad", "3d tendency th from tvdadvection", current_state%tend_3d_th_tvad, &
                                        "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_qv_tvad", "3d tendency qv from tvdadvection", &
                                        current_state%tend_3d_qv_tvad, "kg.kg-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_ql_tvad", "3d tendency ql from tvdadvection", &
                                        current_state%tend_3d_ql_tvad, "kg.kg-1", current_state%time)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qr_tvad", "3d tendency qr from tvdadvection", &
                                          current_state%tend_3d_qr_tvad, "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qi_tvad", "3d tendency qi from tvdadvection", &
                                          current_state%tend_3d_qi_tvad, "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qs_tvad", "3d tendency qs from tvdadvection", &
                                          current_state%tend_3d_qs_tvad, "kg.kg-1", current_state%time)
        call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                          "tend_3d_qg_tvad", "3d tendency qg from tvdadvection", &
                                          current_state%tend_3d_qg_tvad, "kg.kg-1", current_state%time)
      end if
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_tvad", "3d tendency tabs from tvdadvection", &
                                        current_state%tend_3d_tabs_tvad, "???", current_state%time)
    end if
    ! diag viscosity
    if (current_state%viscosity_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_u_visc", "3d tendency u from viscosity", current_state%tend_3d_u_visc, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_v_visc", "3d tendency v from viscosity", current_state%tend_3d_v_visc, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_w_visc", "3d tendency w from viscosity", current_state%tend_3d_w_visc, &
                                        "m.s-1", current_state%time)
    end if
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdAitkenSol", "Aitken soluble aerosol radius", &
                                      current_state%rdAitkenSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdAccumSol", "accumSol soluble aerosol radius", &
                                      current_state%rdAccumSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdCoarseSol", "coarseSol soluble aerosol radius", &
                                      current_state%rdCoarseSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_cloud", "cloud droplet diameter", &
                                      current_state%D0_cloud%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_rain", "rain droplet diameter", &
                                      current_state%D0_rain%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_ice", "ice particle diameter", &
                                      current_state%D0_ice%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_snow", "snow particle diameter", &
                                      current_state%D0_snow%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_graupel", "graupel particle diameter", &
                                      current_state%D0_graupel%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "RH", "relative Humidity", &
                                      current_state%RH%data, "%", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "RHI", "relative Himidity with respect to Ice", &
                                      current_state%RI%data, "%", current_state%time)

    ! SOCRATES cloud_reff, lwrad_hr, swrad_hr
    if (current_state%socrates_enabled .eqv. .true.) then
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "cloud_reff", "cloud reflectivity", &
                                        current_state%cloud_reff%data, " ", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "lwrad_hr", "longwave heating rate", &
                                        current_state%lwrad_hr%data, "K.d-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "swrad_hr", "shortwave heating rate", &
                                        current_state%swrad_hr%data, "K.d-1", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_lw", "3d tendency th from Socrates (lw)", &
                                        current_state%tend_3d_tabs_lw, "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_sw", "3d tendency th from Socrates (sw)", &
                                        current_state%tend_3d_tabs_sw, "K", current_state%time)
      call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                        "tend_3d_tabs_total", "3d tendency th from Socrates", &
                                        current_state%tend_3d_tabs_total, "K", current_state%time)
    end if

    call check(nf90_close(ncdf_id))

    deallocate(current_state%u%data, current_state%v%data, current_state%w%data, current_state%th%data, &
              current_state%p%data, current_state%qv%data, current_state%ql%data, current_state%qr%data, &
              current_state%qi%data, current_state%qs%data, current_state%qg%data, &
              current_state%qAitkenSolMass%data, current_state%qAccumSolMass%data, current_state%qAccumInsolMass%data, &
              current_state%qCoarseSolMass%data, current_state%qCoarseDustMass%data,&
              current_state%nl%data, current_state%nr%data, current_state%ni%data, &
              current_state%ns%data, current_state%ng%data, &
              current_state%nAitkenSolNumber%data, current_state%nAccumSolNumber%data, current_state%nAccumInsolNumber%data, &
              current_state%nCoarseSolNumber%data, current_state%nCoarseDustnumber%data)

    if (current_state%buoyancy_enabled .eqv. .true.) then
      deallocate(current_state%tend_3d_w_buoy)
    end if

    if (current_state%coriolis_enabled .eqv. .true.) then
      deallocate(current_state%tend_3d_u_corio, current_state%tend_3d_v_corio)
    end if

    if (current_state%diffusion_enabled .eqv. .true.) then
      deallocate(current_state%tend_3d_th_diff, current_state%tend_3d_qv_diff, current_state%tend_3d_ql_diff, &
                current_state%tend_3d_tabs_diff)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        deallocate(current_state%tend_3d_qr_diff, current_state%tend_3d_qi_diff, current_state%tend_3d_qs_diff, &
              current_state%tend_3d_qg_diff)
      end if
    end if

    if (current_state%forcing_enabled .eqv. .true.) then
      deallocate(current_state%tend_3d_u_forc, current_state%tend_3d_v_forc, current_state%tend_3d_th_forc, &
              current_state%tend_3d_qv_forc, current_state%tend_3d_ql_forc, current_state%tend_3d_qr_forc, &
              current_state%tend_3d_qi_forc, current_state%tend_3d_qs_forc, current_state%tend_3d_qg_forc, &
              current_state%tend_3d_tabs_forc)
    end if

    if (current_state%pstep_enabled .eqv. .true.) then
      deallocate(current_state%tendp_3d_u_pt, current_state%tendp_3d_v_pt, current_state%tendp_3d_w_pt, &
                current_state%tend_3d_u_pt, current_state%tend_3d_v_pt, current_state%tend_3d_w_pt)
    end if

    if (current_state%pw_advection_enabled .eqv. .true.) then
      deallocate(current_state%tend_3d_u_pwad, current_state%tend_3d_v_pwad, current_state%tend_3d_w_pwad, &
              current_state%tend_3d_th_pwad, current_state%tend_3d_qv_pwad, current_state%tend_3d_ql_pwad, &
              current_state%tend_3d_qr_pwad, current_state%tend_3d_qi_pwad, current_state%tend_3d_qs_pwad,&
              current_state%tend_3d_qg_pwad, current_state%tend_3d_tabs_pwad)
    end if

    if (current_state%stepfields_enabled .eqv. .true.) then
      deallocate(current_state%tend_3d_th_sf, current_state%tend_3d_qv_sf, current_state%tend_3d_ql_sf, &
                current_state%tend_3d_tabs_sf)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        deallocate(current_state%tend_3d_qr_sf, current_state%tend_3d_qi_sf, current_state%tend_3d_qs_sf, &
                  current_state%tend_3d_qg_sf)
      end if
    end if

    if (current_state%th_advection_enabled .eqv. .true.) then
      deallocate(current_state%tend_3d_th_thad, current_state%tend_3d_tabs_thad)
    end if

    if (current_state%tvd_advection_enabled .eqv. .true.) then
      deallocate(current_state%tend_3d_u_tvad, current_state%tend_3d_v_tvad, current_state%tend_3d_w_tvad, &
                current_state%tend_3d_th_tvad, current_state%tend_3d_qv_tvad, current_state%tend_3d_ql_tvad, &
                current_state%tend_3d_tabs_tvad)
      if ((current_state%casim_enabled .eqv. .true.) .or. (current_state%simplecloud_enabled .eqv. .true.)) then
        deallocate(current_state%tend_3d_qr_tvad, current_state%tend_3d_qi_tvad, current_state%tend_3d_qs_tvad, &
                  current_state%tend_3d_qg_tvad)
      end if
    end if

    if (current_state%viscosity_enabled .eqv. .true.) then
      deallocate(current_state%tend_3d_u_visc, current_state%tend_3d_v_visc, current_state%tend_3d_w_visc)
    end if

    if (current_state%socrates_enabled .eqv. .true.) then
      deallocate(current_state%cloud_reff%data, current_state%lwrad_hr%data, current_state%swrad_hr%data, &
              current_state%tend_3d_tabs_lw, current_state%tend_3d_tabs_sw, current_state%tend_3d_tabs_total)
    end if

    deallocate(current_state%rdAitkenSol%data, current_state%rdAccumSol%data, current_state%rdCoarseSol%data, &
              current_state%D0_cloud%data, current_state%D0_rain%data, current_state%D0_ice%data, &
              current_state%D0_snow%data, current_state%D0_graupel%data, &
              current_state%RH%data, current_state%RI%data)


    current_state%diag_3d_done = .true.
    call mpi_send(current_state%diag_3d_done, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, ierr)
    current_state%diag_3d_done = .false.

  end subroutine diagnostic_file_3d_generation



  subroutine diagnostic_file_average_generation(current_state, vertical_grid, io_communicator_arg, &
                                                                            global_grid_z_size, diagnostic_path)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_grid_z_size, io_communicator_arg
    character(len=200), intent(in) :: diagnostic_path

    integer :: z_size_id
    integer :: ncdf_id, ierr
    integer :: i,ls1,ls2
    character(len=250) :: time_char, timestep_char
    character(len=LONG_STRING_LENGTH) :: unique_filename

    global_grid_z_size_netcdf = global_grid_z_size

    allocate(current_state%mean_th(global_grid_z_size))
    allocate(current_state%mean_u(global_grid_z_size))
    allocate(current_state%mean_v(global_grid_z_size))
    allocate(current_state%mean_qv(global_grid_z_size))
    allocate(current_state%mean_ql(global_grid_z_size))
    allocate(current_state%mean_qr(global_grid_z_size))
    allocate(current_state%mean_qi(global_grid_z_size))
    allocate(current_state%mean_qs(global_grid_z_size))
    allocate(current_state%mean_qg(global_grid_z_size))
    allocate(current_state%mean_qAitkenSolMass(global_grid_z_size))
    allocate(current_state%mean_qAccumSolMass(global_grid_z_size))
    allocate(current_state%mean_qAccumInsolMass(global_grid_z_size))
    allocate(current_state%mean_qCoarseSolMass(global_grid_z_size))
    allocate(current_state%mean_qCoarseDustMass(global_grid_z_size))
    allocate(current_state%mean_nl(global_grid_z_size))
    allocate(current_state%mean_nr(global_grid_z_size))
    allocate(current_state%mean_ni(global_grid_z_size))
    allocate(current_state%mean_ns(global_grid_z_size))
    allocate(current_state%mean_ng(global_grid_z_size))
    allocate(current_state%mean_nAitkenSolNumber(global_grid_z_size))
    allocate(current_state%mean_nAccumSolNumber(global_grid_z_size))
    allocate(current_state%mean_nAccumInsolNumber(global_grid_z_size))
    allocate(current_state%mean_nCoarseSolNumber(global_grid_z_size))
    allocate(current_state%mean_nCoarseDustnumber(global_grid_z_size))

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

    if (current_state%diagnostics_averaged_enabled .eqv. .true.) then
      call mpi_recv(current_state%mean_th, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_u, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_v, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qv, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_ql, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qr, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qi, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qs, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qg, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qAitkenSolMass, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qAccumSolMass, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qAccumInsolMass, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qCoarseSolMass, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_qCoarseDustMass, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_nl, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_nr, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_ni, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_ns, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_ng, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_nAitkenSolNumber, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_nAccumSolNumber, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_nAccumInsolNumber, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_nCoarseSolNumber, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_nCoarseDustnumber, (global_grid_z_size), &
                                                      MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(current_state%mean_mo_l, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_friction_vel, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_z0, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%mean_z0th, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if

    ls1 = len_trim(diagnostic_path)
    ls2 = 0
    do i = 1,ls1
        if(diagnostic_path(i:i).ne.' ') then
          ls2 = ls2 + 1
        end if
    end do

    write (unit=time_char,fmt=*) int(current_state%time)
    write (unit=timestep_char,fmt=*) current_state%timestep

    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = diagnostic_path(:ls2)//"/full_diag_1d_average_time_"//trim(adjustl(time_char))//".nc"
    else
      unique_filename = diagnostic_path(:ls2)//"/full_diag_1d_average_timestep_"//trim(adjustl(timestep_char))//".nc"
    end if

    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "t", 1, time_id))
    call check(nf90_def_dim(ncdf_id, "z", global_grid_z_size, z_size_id))


    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)

    if (current_state%diagnostics_averaged_enabled .eqv. .true.) then
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                       "mean_th", "average potential temperature profile", current_state%mean_th, &
                                       "K", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                       "mean_u", "average u profile", current_state%mean_u, &
                                       "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_v", "average v profile", current_state%mean_v, &
                                        "m.s-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qv", "average water vapour mass profile", current_state%mean_qv, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_ql", "average cloud water mass profile", current_state%mean_ql, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qr", "average rain mass profile", current_state%mean_qr, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qi", "average ice mass profile", current_state%mean_qi, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qs", "average snow mass profile", current_state%mean_qs, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qg", "average graupel mass profile", current_state%mean_qg, &
                                        "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qAitkenSolMass", "average mass Aitken mode profile", &
                                        current_state%mean_qAitkenSolMass, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qAccumSolMass", "average mass soluble Accumulation mode profile", &
                                        current_state%mean_qAccumSolMass, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qAccumInsolMass", "average mass insoluble Accumulation mode profile", &
                                        current_state%mean_qAccumInsolMass, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qCoarseSolMass", "average mass soluble Coarse mode profile", &
                                        current_state%mean_qCoarseSolMass, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_qCoarseDustMass", "average mass Dust profile", &
                                        current_state%mean_qCoarseDustMass, "kg.kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_nl", "average cloud water number concentration profile", current_state%mean_nl, &
                                        "kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_nr", "average rain number concentration profile", current_state%mean_nr, &
                                        "kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_ni", "average ice number concentration profile", current_state%mean_ni, &
                                        "kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_ns", "average snow number concentration profile", current_state%mean_ns, &
                                        "kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_ng", "average graupel number concentration profile", current_state%mean_ng, &
                                        "kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_nAitkenSolNumber", "average number concentration Aitken mode profile", &
                                        current_state%mean_nAitkenSolNumber, "kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_nAccumSolNumber", "average number concentration soluble Accumulation mode profile", &
                                        current_state%mean_nAccumSolNumber, "kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                    "mean_nAccumInsolNumber", "average number concentration insoluble Accumulation mode profile", &
                                    current_state%mean_nAccumInsolNumber, "kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_nCoarseSolNumber", "average number concentration soluble Coarse mode profile", &
                                        current_state%mean_nCoarseSolNumber, "kg-1", current_state%time)
      call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                        "mean_nCoarseDustnumber", "average number concentration Dust profile", &
                                        current_state%mean_nCoarseDustnumber, "kg-1", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                        "MO_l", "Monin-Obukhov length", current_state%mean_mo_l, "m", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                        "friction_vel", "friction velocity", current_state%mean_friction_vel, "m.s-1", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                        "z0", "roughness length for velocity", current_state%mean_z0, "m", current_state%time)
      call define_and_write_variable_real_scalar(ncdf_id, &
                        "z0th", "roughness length for scalars", current_state%mean_z0th, "m", current_state%time)
    end if

    call check(nf90_close(ncdf_id))

    deallocate(current_state%mean_th, current_state%mean_u, current_state%mean_v, current_state%mean_qv, &
              current_state%mean_ql, current_state%mean_qr, current_state%mean_qi, current_state%mean_qs, &
              current_state%mean_qg, current_state%mean_qAitkenSolMass, current_state%mean_qAccumSolMass, &
              current_state%mean_qAccumInsolMass, current_state%mean_qCoarseSolMass, current_state%mean_qCoarseDustMass, &
              current_state%mean_nl, current_state%mean_nr, current_state%mean_ni, current_state%mean_ns, &
              current_state%mean_ng, current_state%mean_nAitkenSolNumber, current_state%mean_nAccumSolNumber, &
              current_state%mean_nAccumInsolNumber, current_state%mean_nCoarseSolNumber, current_state%mean_nCoarseDustnumber)

  end subroutine diagnostic_file_average_generation




  subroutine checkpoint_file_generation(current_state, vertical_grid, io_communicator_arg, global_array_size, &
                                global_grid_z_size, global_grid_y_size, global_grid_x_size, time_model, checkpoint_path)
    type(model_state_type), target, intent(inout) :: current_state
    type(vertical_grid_configuration_type), target, intent(inout) :: vertical_grid
    integer, intent(in) :: global_array_size, global_grid_z_size, global_grid_y_size, global_grid_x_size, &
                           io_communicator_arg
    real(kind=DEFAULT_PRECISION), intent(out):: time_model
    character(len=200) :: checkpoint_path

    integer :: z_size_id, y_size_id, x_size_id!, scalar_size_id, time_id
    integer :: ncdf_id, ierr
    integer :: i,ls1,ls2
    character(len=LONG_STRING_LENGTH) :: unique_filename

    global_grid_z_size_netcdf = global_grid_z_size
    global_grid_y_size_netcdf = global_grid_y_size
    global_grid_x_size_netcdf = global_grid_x_size

    call mpi_recv(current_state%dtm_new, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%absolute_new_dtm, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%ugal, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%vgal, 1, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    time_model = current_state%time! + current_state%dtm

    allocate(current_state%global_grid%configuration%vertical%thref(global_grid_z_size))
    allocate(vertical_grid%prefn(global_grid_z_size))
    allocate(vertical_grid%olubar(global_grid_z_size))
    allocate(vertical_grid%olvbar(global_grid_z_size))
    allocate(vertical_grid%olthbar(global_grid_z_size))
    allocate(vertical_grid%olqvbar(global_grid_z_size))
    allocate(vertical_grid%olqlbar(global_grid_z_size))
    allocate(vertical_grid%olqrbar(global_grid_z_size))
    allocate(vertical_grid%olqibar(global_grid_z_size))
    allocate(vertical_grid%olqsbar(global_grid_z_size))
    allocate(vertical_grid%olqgbar(global_grid_z_size))
    allocate(vertical_grid%olqAitkenSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olqAccumSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olqAccumInsolMassbar(global_grid_z_size))
    allocate(vertical_grid%olqCoarseSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olqCoarseDustMassbar(global_grid_z_size))
    allocate(vertical_grid%olnlbar(global_grid_z_size))
    allocate(vertical_grid%olnrbar(global_grid_z_size))
    allocate(vertical_grid%olnibar(global_grid_z_size))
    allocate(vertical_grid%olnsbar(global_grid_z_size))
    allocate(vertical_grid%olngbar(global_grid_z_size))
    allocate(vertical_grid%olnAitkenSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olnAccumSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olnAccumInsolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olnCoarseSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olnCoarseDustnumberbar(global_grid_z_size))
    allocate(vertical_grid%olzubar(global_grid_z_size))
    allocate(vertical_grid%olzvbar(global_grid_z_size))
    allocate(vertical_grid%olzthbar(global_grid_z_size))
    allocate(vertical_grid%olzqvbar(global_grid_z_size))
    allocate(vertical_grid%olzqlbar(global_grid_z_size))
    allocate(vertical_grid%olzqrbar(global_grid_z_size))
    allocate(vertical_grid%olzqibar(global_grid_z_size))
    allocate(vertical_grid%olzqsbar(global_grid_z_size))
    allocate(vertical_grid%olzqgbar(global_grid_z_size))
    allocate(vertical_grid%olzqAitkenSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olzqAccumSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olzqAccumInsolMassbar(global_grid_z_size))
    allocate(vertical_grid%olzqCoarseSolMassbar(global_grid_z_size))
    allocate(vertical_grid%olzqCoarseDustMassbar(global_grid_z_size))
    allocate(vertical_grid%olznlbar(global_grid_z_size))
    allocate(vertical_grid%olznrbar(global_grid_z_size))
    allocate(vertical_grid%olznibar(global_grid_z_size))
    allocate(vertical_grid%olznsbar(global_grid_z_size))
    allocate(vertical_grid%olzngbar(global_grid_z_size))
    allocate(vertical_grid%olznAitkenSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olznAccumSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olznAccumInsolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olznCoarseSolNumberbar(global_grid_z_size))
    allocate(vertical_grid%olznCoarseDustnumberbar(global_grid_z_size))
    current_state%global_grid%configuration%vertical%thref = 0.0_DEFAULT_PRECISION
    vertical_grid%prefn = 0.0_DEFAULT_PRECISION
    vertical_grid%olubar = 0.0_DEFAULT_PRECISION
    vertical_grid%olvbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olthbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqvbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqlbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqrbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqibar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqsbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqgbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqAitkenSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqAccumSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqAccumInsolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqCoarseSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olqCoarseDustMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnlbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnrbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnibar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnsbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olngbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnAitkenSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnAccumSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnAccumInsolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnCoarseSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olnCoarseDustnumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzubar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzvbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzthbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqvbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqlbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqrbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqibar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqsbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqgbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqAitkenSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqAccumSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqAccumInsolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqCoarseSolMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzqCoarseDustMassbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznlbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznrbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznibar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznsbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olzngbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznAitkenSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznAccumSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznAccumInsolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznCoarseSolNumberbar = 0.0_DEFAULT_PRECISION
    vertical_grid%olznCoarseDustnumberbar = 0.0_DEFAULT_PRECISION
    call mpi_recv(current_state%global_grid%configuration%vertical%thref, &
        global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(vertical_grid%prefn, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(vertical_grid%olubar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(vertical_grid%olvbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%th%active) then
      call mpi_recv(vertical_grid%olthbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%number_q_fields .ne. 0) then
      call mpi_recv(vertical_grid%olqvbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqlbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqrbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqibar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqsbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqgbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqAitkenSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqAccumSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqAccumInsolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqCoarseSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olqCoarseDustMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnlbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnrbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnibar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnsbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olngbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnAitkenSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnAccumSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnAccumInsolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnCoarseSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olnCoarseDustnumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    call mpi_recv(vertical_grid%olzubar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(vertical_grid%olzvbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%th%active) then
      call mpi_recv(vertical_grid%olzthbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%number_q_fields .ne. 0) then
      call mpi_recv(vertical_grid%olzqvbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqlbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqrbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqibar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqsbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqgbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqAitkenSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqAccumSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqAccumInsolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqCoarseSolMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzqCoarseDustMassbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznlbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznrbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznibar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznsbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olzngbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznAitkenSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznAccumSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznAccumInsolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznCoarseSolNumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(vertical_grid%olznCoarseDustnumberbar, global_grid_z_size, MPI_DOUBLE, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if


    allocate(current_state%u%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%v%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%w%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%th%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%p%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qv%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ql%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qi%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qs%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qg%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAitkenSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAccumSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qAccumInsolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qCoarseSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%qCoarseDustMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nl%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ni%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ns%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%ng%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAitkenSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAccumSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nAccumInsolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nCoarseSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%nCoarseDustnumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zu%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zv%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zw%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zth%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqv%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zql%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqi%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqs%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqg%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqAitkenSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqAccumSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqAccumInsolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqCoarseSolMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zqCoarseDustMass%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znl%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znr%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zni%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zns%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%zng%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znAitkenSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znAccumSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znAccumInsolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znCoarseSolNumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%znCoarseDustnumber%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdAitkenSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdAccumSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%rdCoarseSol%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_cloud%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_rain%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_ice%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_snow%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%D0_graupel%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%RH%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    allocate(current_state%RI%data(global_grid_z_size, global_grid_y_size, global_grid_x_size))
    current_state%u%data = 0.0_DEFAULT_PRECISION
    current_state%v%data = 0.0_DEFAULT_PRECISION
    current_state%w%data = 0.0_DEFAULT_PRECISION
    current_state%th%data = 0.0_DEFAULT_PRECISION
    current_state%p%data = 0.0_DEFAULT_PRECISION
    current_state%qv%data = 0.0_DEFAULT_PRECISION
    current_state%ql%data = 0.0_DEFAULT_PRECISION
    current_state%qr%data = 0.0_DEFAULT_PRECISION
    current_state%qi%data = 0.0_DEFAULT_PRECISION
    current_state%qs%data = 0.0_DEFAULT_PRECISION
    current_state%qg%data = 0.0_DEFAULT_PRECISION
    current_state%qAitkenSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qAccumSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qAccumInsolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qCoarseSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%qCoarseDustMass%data = 0.0_DEFAULT_PRECISION
    current_state%nl%data = 0.0_DEFAULT_PRECISION
    current_state%nr%data = 0.0_DEFAULT_PRECISION
    current_state%ni%data = 0.0_DEFAULT_PRECISION
    current_state%ns%data = 0.0_DEFAULT_PRECISION
    current_state%ng%data = 0.0_DEFAULT_PRECISION
    current_state%nAitkenSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nAccumSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nAccumInsolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nCoarseSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%nCoarseDustnumber%data = 0.0_DEFAULT_PRECISION
    current_state%zu%data = 0.0_DEFAULT_PRECISION
    current_state%zv%data = 0.0_DEFAULT_PRECISION
    current_state%zw%data = 0.0_DEFAULT_PRECISION
    current_state%zth%data = 0.0_DEFAULT_PRECISION
    current_state%zqv%data = 0.0_DEFAULT_PRECISION
    current_state%zql%data = 0.0_DEFAULT_PRECISION
    current_state%zqr%data = 0.0_DEFAULT_PRECISION
    current_state%zqi%data = 0.0_DEFAULT_PRECISION
    current_state%zqs%data = 0.0_DEFAULT_PRECISION
    current_state%zqg%data = 0.0_DEFAULT_PRECISION
    current_state%zqAitkenSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%zqAccumSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%zqAccumInsolMass%data = 0.0_DEFAULT_PRECISION
    current_state%zqCoarseSolMass%data = 0.0_DEFAULT_PRECISION
    current_state%zqCoarseDustMass%data = 0.0_DEFAULT_PRECISION
    current_state%znl%data = 0.0_DEFAULT_PRECISION
    current_state%znr%data = 0.0_DEFAULT_PRECISION
    current_state%zni%data = 0.0_DEFAULT_PRECISION
    current_state%zns%data = 0.0_DEFAULT_PRECISION
    current_state%zng%data = 0.0_DEFAULT_PRECISION
    current_state%znAitkenSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%znAccumSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%znAccumInsolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%znCoarseSolNumber%data = 0.0_DEFAULT_PRECISION
    current_state%znCoarseDustnumber%data = 0.0_DEFAULT_PRECISION
    current_state%rdAitkenSol%data = 0.0_DEFAULT_PRECISION
    current_state%rdAccumSol%data = 0.0_DEFAULT_PRECISION
    current_state%rdCoarseSol%data = 0.0_DEFAULT_PRECISION
    current_state%D0_cloud%data = 0.0_DEFAULT_PRECISION
    current_state%D0_rain%data = 0.0_DEFAULT_PRECISION
    current_state%D0_ice%data = 0.0_DEFAULT_PRECISION
    current_state%D0_snow%data = 0.0_DEFAULT_PRECISION
    current_state%D0_graupel%data = 0.0_DEFAULT_PRECISION
    current_state%RH%data = 0.0_DEFAULT_PRECISION
    current_state%RI%data = 0.0_DEFAULT_PRECISION
    call mpi_recv(current_state%u%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%v%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%w%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%th%active) then
      call mpi_recv(current_state%th%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    call mpi_recv(current_state%p%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%number_q_fields .ne. 0) then
      call mpi_recv(current_state%qv%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ql%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qi%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qs%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qg%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAitkenSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qAccumInsolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%qCoarseDustMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nl%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ni%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ns%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%ng%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAitkenSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nAccumInsolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%nCoarseDustnumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    call mpi_recv(current_state%zu%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%zv%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%zw%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    if (current_state%th%active) then
      call mpi_recv(current_state%zth%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    if (current_state%number_q_fields .ne. 0) then
      call mpi_recv(current_state%zqv%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zql%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqi%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqs%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqg%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqAitkenSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqAccumSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqAccumInsolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqCoarseSolMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zqCoarseDustMass%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znl%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znr%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zni%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zns%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%zng%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znAitkenSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znAccumSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znAccumInsolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znCoarseSolNumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      call mpi_recv(current_state%znCoarseDustnumber%data, global_array_size*2, MPI_REAL, 1, 1000, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    end if
    call mpi_recv(current_state%rdAitkenSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%rdAccumSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%rdCoarseSol%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_cloud%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_rain%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_ice%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_snow%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%D0_graupel%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%RH%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_recv(current_state%RI%data, global_array_size*2, MPI_REAL, 1, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

    ls1 = len_trim(checkpoint_path)
    ls2 = 0
    do i = 1,ls1
        if(checkpoint_path(i:i).ne.' ') then
          ls2 = ls2 + 1
        endif
    enddo
    if (current_state%time_frequency_enabled .eqv. .true.) then
      unique_filename = checkpoint_path(:ls2)//"/checkpt_run_time_"//trim(conv_to_string(int(current_state%time)))//".nc"
    else
      unique_filename = checkpoint_path(:ls2)//"/checkpt_run_timestep_"//trim(conv_to_string(current_state%timestep))//".nc"
    end if
    call check(nf90_create(unique_filename, ior(NF90_NETCDF4, NF90_MPIIO), ncdf_id, &
            comm = io_communicator_arg, info = MPI_INFO_NULL))

    call check(nf90_def_dim(ncdf_id, "t", 1, time_id))
    call check(nf90_def_dim(ncdf_id, "z", global_grid_z_size, z_size_id))
    call check(nf90_def_dim(ncdf_id, "y", global_grid_y_size, y_size_id))
    call check(nf90_def_dim(ncdf_id, "x", global_grid_x_size, x_size_id))
    !call check(nf90_def_dim(ncdf_id, "scalar_size", 1, scalar_size_id))


    call define_and_write_variable_real_scalar(ncdf_id, &
            "time", "time", current_state%time, "s", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "timestep", "timestep number", current_state%timestep, "no unit", current_state%time)
    call define_and_write_variable_integer_scalar(ncdf_id, &
            "last_cfl_timestep", "last cfl timestep", current_state%last_cfl_timestep, "no unit", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "zh", "heights at w levels", vertical_grid%z, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "zhn", "heights at pressure levels", vertical_grid%zn, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "dzh", "height of cells", vertical_grid%dz, "m", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "rho", "cells density", vertical_grid%rho, "kg.m-3", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                          "rhon", "cells density (pressure levels)", vertical_grid%rho, "kg.m-3", current_state%time)


    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_bottom", "x bottom coordinate", current_state%global_grid%bottom(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_top", "x top coordinate", current_state%global_grid%top(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "x_resolution", "resolution along x axis", current_state%global_grid%resolution(3), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_bottom", "y bottom coordinate", current_state%global_grid%bottom(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_top", "y top coordinate", current_state%global_grid%top(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "y_resolution", "resolution along y axis", current_state%global_grid%resolution(2), "m", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "dtm", "time timestep", current_state%dtm, "s", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "dtm_new", "new time timestep (cfltest activated)", current_state%absolute_new_dtm, "s", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "absolute_new_dtm", "absolute time timestep", current_state%absolute_new_dtm, "s", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "ugal", "galilean u wind component", current_state%ugal, "m.s-1", current_state%time)
    call define_and_write_variable_real_scalar(ncdf_id, &
                    "vgal", "galilean v wind component", current_state%vgal, "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "thref", "reference potential temperature profile", &
                                      current_state%global_grid%configuration%vertical%thref, "K", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "prefn", "reference pressure profil", vertical_grid%prefn, "Pa", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olubar", "u mean along z axis", vertical_grid%olubar, "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olvbar", "v mean along z axis", vertical_grid%olvbar, "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olthbar", "th mean along z axis", vertical_grid%olthbar, "K", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqvbar", "qv mean along z axis", vertical_grid%olqvbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqlbar", "ql mean along z axis", vertical_grid%olqlbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqrbar", "qr mean along z axis", vertical_grid%olqrbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqibar", "qi mean along z axis", vertical_grid%olqibar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqsbar", "qs mean along z axis", vertical_grid%olqsbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olqgbar", "qg mean along z axis", vertical_grid%olqgbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqAitkenSolMassbar", "qAitkenSolMass mean along z axis", vertical_grid%olqAitkenSolMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqAccumSolMassbar", "qAccumSolMass mean along z axis", vertical_grid%olqAccumSolMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqAccumInsolMassbar", "qAccumInsolMass mean along z axis", vertical_grid%olqAccumInsolMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqCoarseSolMassbar", "qCoarseSolMass mean along z axis", vertical_grid%olqCoarseSolMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olqCoarseDustMassbar", "qCoarseDustMass mean along z axis", vertical_grid%olqCoarseDustMassbar, &
                    "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olnlbar", "nl mean along z axis", vertical_grid%olnlbar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olnrbar", "nr mean along z axis", vertical_grid%olnrbar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olnibar", "ni mean along z axis", vertical_grid%olnibar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olnsbar", "ns mean along z axis", vertical_grid%olnsbar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olngbar", "ng mean along z axis", vertical_grid%olngbar, "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnAitkenSolNumberbar", "nAitkenSolNumber mean along z axis", vertical_grid%olnAitkenSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnAccumSolNumberbar", "nAccumSolNumber mean along z axis", vertical_grid%olnAccumSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnAccumInsolNumberbar", "nAccumInsolNumber mean along z axis", vertical_grid%olnAccumInsolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnCoarseSolNumberbar", "nCoarseSolNumber mean along z axis", vertical_grid%olnCoarseSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olnCoarseDustnumberbar", "nCoarseDustnumber mean along z axis", vertical_grid%olnCoarseDustnumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzubar", "previous timestep u mean along z axis", vertical_grid%olzubar, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzvbar", "previous timestep v mean along z axis", vertical_grid%olzvbar, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzthbar", "previous timestep th mean along z axis", vertical_grid%olzthbar, &
                                      "K", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqvbar", "previous timestep qv mean along z axis", vertical_grid%olzqvbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqlbar", "previous timestep ql mean along z axis", vertical_grid%olzqlbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqrbar", "previous timestep qr mean along z axis", vertical_grid%olzqrbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqibar", "previous timestep qi mean along z axis", vertical_grid%olzqibar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqsbar", "previous timestep qs mean along z axis", vertical_grid%olzqsbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzqgbar", "previous timestep qg mean along z axis", vertical_grid%olzqgbar, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqAitkenSolMassbar", "previous timestep qAitkenSolMass mean along z axis", &
                                                      vertical_grid%olzqAitkenSolMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqAccumSolMassbar", "previous timestep qAccumSolMass mean along z axis", &
                                                    vertical_grid%olzqAccumSolMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqAccumInsolMassbar", "previous timestep qAccumInsolMass mean along z axis", &
                                                      vertical_grid%olzqAccumInsolMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqCoarseSolMassbar", "previous timestep qCoarseSolMass mean along z axis", &
                                                     vertical_grid%olzqCoarseSolMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                    "olzqCoarseDustMassbar", "previous timestep qCoarseDustMass mean along z axis", &
                                                      vertical_grid%olzqCoarseDustMassbar, "kg.kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olznlbar", "previous timestep nl mean along z axis", vertical_grid%olznlbar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olznrbar", "previous timestep nr mean along z axis", vertical_grid%olznrbar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olznibar", "previous timestep ni mean along z axis", vertical_grid%olznibar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olznsbar", "previous timestep ns mean along z axis", vertical_grid%olznsbar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, &
                                      "olzngbar", "previous timestep ng mean along z axis", vertical_grid%olzngbar, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id,"olznAitkenSolNumberbar", &
                    "previous timestep nAitkenSolNumber mean along z axis", vertical_grid%olznAitkenSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, "olznAccumSolNumberbar", &
                    "previous timestep nAccumSolNumber mean along z axis", vertical_grid%olznAccumSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, "olznAccumInsolNumberbar", &
                    "previous timestep nAccumInsolNumber mean along z axis", vertical_grid%olznAccumInsolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, "olznCoarseSolNumberbar", &
                    "previous timestep nCoarseSolNumber mean along z axis", vertical_grid%olznCoarseSolNumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_1D(ncdf_id, z_size_id, "olznCoarseDustnumberbar", &
                    "previous timestep nCoarseDustnumber mean along z axis", vertical_grid%olznCoarseDustnumberbar, &
                    "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "u", "u wind component", current_state%u%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "v", "v wind component", current_state%v%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "w", "w wind component", current_state%w%data, "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "th", "potentiel temperature perturbations", current_state%th%data, "K", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "p", "pressure perturbations", current_state%p%data, "Pa", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qv", "water vapour mass mixing ratio", current_state%qv%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ql", "liquid water mass mixing ratio", current_state%ql%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qr", "rain water mass mixing ratio", current_state%qr%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qi", "ice water mass mixing ratio", current_state%qi%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qs", "snow water mass mixing ratio", current_state%qs%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qg", "graupel water mass mixing ratio", current_state%qg%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAitkenSolMass", "Aitken soluble aerosol mass mixing ratio", &
                                      current_state%qAitkenSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAccumSolMass", "accumulation soluble aerosol mass mixing ratio",&
                                      current_state%qAccumSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qAccumInsolMass", "accumulation insoluble aerosol mass mixing ratio",&
                                      current_state%qAccumInsolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qCoarseSolMass", "coarse soluble aerosol mixing ratio",&
                                      current_state%qCoarseSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "qCoarseDustMass", "dust coarse aerosol mixing ratio",&
                                      current_state%qCoarseDustMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nl", "liquid water number concentration", current_state%nl%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nr", "rain water number concentration", current_state%nr%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ni", "ice water number concentration", current_state%ni%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ns", "snow water number concentration", current_state%ns%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "ng", "graupel water number concentration", current_state%ng%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAitkenSolNumber", "Aitken soluble aerosol number concentration", &
                                      current_state%nAitkenSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAccumSolNumber", "accumulation soluble aerosol number concentration",&
                                      current_state%nAccumSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nAccumInsolNumber", "accumulation insoluble aerosol number concentration",&
                                      current_state%nAccumInsolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nCoarseSolNumber", "coarse soluble aerosol number concentration",&
                                      current_state%nCoarseSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "nCoarseDustnumber", "dust coarse aerosol number concentration",&
                                      current_state%nCoarseDustnumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zu", "previous timestep u wind component", current_state%zu%data, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zv", "previous timestep v wind component", current_state%zv%data, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zw", "previous timestep w wind component", current_state%zw%data, &
                                      "m.s-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zth", "previous timestep potentiel temperature", current_state%zth%data, &
                                      "K", current_state%time)

    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqv", "previous timestep water vapour mass mixing ratio", current_state%zqv%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zql", "previous timestep liquid water mass mixing ratio", current_state%zql%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqr", "previous timestep rain water mass mixing ratio", current_state%zqr%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqi", "previous timestep ice water mass mixing ratio", current_state%zqi%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqs", "previous timestep snow water mass mixing ratio", current_state%zqs%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqg", "previous timestep graupel water mass mixing ratio", current_state%zqg%data, &
                                      "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqAitkenSolMass", "previous timestep Aitken soluble aerosol mass mixing ratio", &
                                      current_state%zqAitkenSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqAccumSolMass", "previous timestep accumulation soluble aerosol mass mixing ratio",&
                                      current_state%zqAccumSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqAccumInsolMass", "previous timestep accumulation insoluble aerosol mass mixing ratio",&
                                      current_state%zqAccumInsolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqCoarseSolMass", "previous timestep coarse soluble aerosol mixing ratio",&
                                      current_state%zqCoarseSolMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zqCoarseDustMass", "previous timestep dust coarse aerosol mixing ratio",&
                                      current_state%zqCoarseDustMass%data, "kg.kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znl", "previous timestep liquid water number concentration", current_state%znl%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znr", "previous timestep rain water number concentration", current_state%znr%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zni", "previous timestep ice water number concentration", current_state%zni%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zns", "previous timestep snow water number concentration", current_state%zns%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "zng", "previous timestep graupel water number concentration", current_state%zng%data, &
                                      "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znAitkenSolNumber", "previous timestep Aitken soluble aerosol number concentration", &
                                      current_state%znAitkenSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znAccumSolNumber", "previous timestep accumulation soluble aerosol number concentration",&
                                      current_state%znAccumSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                  "znAccumInsolNumber", "previous timestep accumulation insoluble aerosol number concentration",&
                                  current_state%znAccumInsolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znCoarseSolNumber", "previous timestep coarse soluble aerosol number concentration",&
                                      current_state%znCoarseSolNumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "znCoarseDustnumber", "previous timestep dust coarse aerosol number concentration",&
                                      current_state%znCoarseDustnumber%data, "kg-1", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdAitkenSol", "Aitken soluble aerosol radius", &
                                      current_state%rdAitkenSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdAccumSol", "accumSol soluble aerosol radius", &
                                      current_state%rdAccumSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "rdCoarseSol", "coarseSol soluble aerosol radius", &
                                      current_state%rdCoarseSol%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_cloud", "cloud droplet diameter", &
                                      current_state%D0_cloud%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_rain", "rain droplet diameter", &
                                      current_state%D0_rain%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_ice", "ice particle diameter", &
                                      current_state%D0_ice%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_snow", "snow particle diameter", &
                                      current_state%D0_snow%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "D0_graupel", "graupel particle diameter", &
                                      current_state%D0_graupel%data, "m", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "RH", "relative Humidity", &
                                      current_state%RH%data, "%", current_state%time)
    call define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, &
                                      "RHI", "relative Himidity with respect to Ice", &
                                      current_state%RI%data, "%", current_state%time)

    call check(nf90_close(ncdf_id))

    deallocate( current_state%global_grid%configuration%vertical%thref, vertical_grid%prefn, &
                vertical_grid%olubar, vertical_grid%olvbar, vertical_grid%olthbar, &
                vertical_grid%olqvbar, vertical_grid%olqlbar, vertical_grid%olqrbar, &
                vertical_grid%olqibar, vertical_grid%olqsbar, vertical_grid%olqgbar, &
                vertical_grid%olqAitkenSolMassbar, vertical_grid%olqAccumSolMassbar, vertical_grid%olqAccumInsolMassbar, &
                vertical_grid%olqCoarseSolMassbar, vertical_grid%olqCoarseDustMassbar, &
                vertical_grid%olnlbar, vertical_grid%olnrbar, vertical_grid%olnibar,&
                vertical_grid%olnsbar, vertical_grid%olngbar, &
                vertical_grid%olnAitkenSolNumberbar, vertical_grid%olnAccumSolNumberbar, vertical_grid%olnAccumInsolNumberbar, &
                vertical_grid%olnCoarseSolNumberbar, vertical_grid%olnCoarseDustnumberbar, &
                vertical_grid%olzubar, vertical_grid%olzvbar, vertical_grid%olzthbar, &
                vertical_grid%olzqvbar, vertical_grid%olzqlbar, vertical_grid%olzqrbar, &
                vertical_grid%olzqibar, vertical_grid%olzqsbar, vertical_grid%olzqgbar, &
                vertical_grid%olzqAitkenSolMassbar, vertical_grid%olzqAccumSolMassbar, vertical_grid%olzqAccumInsolMassbar, &
                vertical_grid%olzqCoarseSolMassbar, vertical_grid%olzqCoarseDustMassbar, &
                vertical_grid%olznlbar, vertical_grid%olznrbar, vertical_grid%olznibar,&
                vertical_grid%olznsbar, vertical_grid%olzngbar, &
                vertical_grid%olznAitkenSolNumberbar, vertical_grid%olznAccumSolNumberbar, &
                vertical_grid%olznAccumInsolNumberbar, &
                vertical_grid%olznCoarseSolNumberbar, vertical_grid%olznCoarseDustnumberbar, &
              )
    deallocate(current_state%u%data, current_state%v%data, current_state%w%data, current_state%th%data, &
              current_state%p%data, current_state%qv%data, current_state%ql%data, current_state%qr%data, &
              current_state%qi%data, current_state%qs%data, current_state%qg%data, &
              current_state%qAitkenSolMass%data, current_state%qAccumSolMass%data, current_state%qAccumInsolMass%data, &
              current_state%qCoarseSolMass%data, current_state%qCoarseDustMass%data,&
              current_state%nl%data, current_state%nr%data, current_state%ni%data, &
              current_state%ns%data, current_state%ng%data, &
              current_state%nAitkenSolNumber%data, current_state%nAccumSolNumber%data, current_state%nAccumInsolNumber%data, &
              current_state%nCoarseSolNumber%data, current_state%nCoarseDustnumber%data,&
              current_state%zu%data, current_state%zv%data, current_state%zw%data, current_state%zth%data, &
              current_state%zqv%data, current_state%zql%data, current_state%zqr%data, &
              current_state%zqi%data, current_state%zqs%data, current_state%zqg%data, &
              current_state%zqAitkenSolMass%data, current_state%zqAccumSolMass%data, current_state%zqAccumInsolMass%data, &
              current_state%zqCoarseSolMass%data, current_state%zqCoarseDustMass%data,&
              current_state%znl%data, current_state%znr%data, current_state%zni%data, &
              current_state%zns%data, current_state%zng%data, &
              current_state%znAitkenSolNumber%data, current_state%znAccumSolNumber%data, current_state%znAccumInsolNumber%data, &
              current_state%znCoarseSolNumber%data, current_state%znCoarseDustnumber%data,&
              current_state%rdAitkenSol%data, current_state%rdAccumSol%data, current_state%rdCoarseSol%data, &
              current_state%D0_cloud%data, current_state%D0_rain%data, current_state%D0_ice%data, &
              current_state%D0_snow%data, current_state%D0_graupel%data, &
              current_state%RH%data, current_state%RI%data)

    current_state%checkpt_done = .true.
    call mpi_send(current_state%checkpt_done, 1, MPI_LOGICAL, 1, 1000, MPI_COMM_WORLD, ierr)
    current_state%checkpt_done = .false.

    !print*,"checkpoint_file_generation IO"
  end subroutine checkpoint_file_generation




  subroutine define_and_write_variable_integer_scalar(ncdf_id, name_var, standard_name, data_scalar, unit, time)
    integer, intent(in) :: ncdf_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    integer, intent(in) :: data_scalar

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,1) :: data_scalar_time
    data_scalar_time(1,1) = data_scalar
    call check(nf90_def_var(ncdf_id, name_var, NF90_INT, (/ time_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_scalar_time))
  end subroutine define_and_write_variable_integer_scalar




  subroutine define_and_write_variable_real_scalar(ncdf_id, name_var, standard_name, data_scalar, unit, time)
    integer, intent(in) :: ncdf_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    real(kind=DEFAULT_PRECISION), intent(in) :: data_scalar

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,1) :: data_scalar_time
    data_scalar_time(1,1) = data_scalar

    call check(nf90_def_var(ncdf_id, name_var, NF90_DOUBLE, (/ time_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_scalar_time))
  end subroutine define_and_write_variable_real_scalar




  subroutine define_and_write_variable_1D(ncdf_id, z_size_id, name_var, standard_name, data_1D, unit, time)
    integer, intent(in) :: ncdf_id, z_size_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: data_1D

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,global_grid_z_size_netcdf) :: data_1D_time
    data_1D_time(1,:) = data_1D

    call check(nf90_def_var(ncdf_id, name_var, NF90_DOUBLE, (/ time_id, z_size_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_1D_time))
  end subroutine define_and_write_variable_1D




  subroutine define_and_write_variable_1D_bin(ncdf_id, z_size_id, bin_size_id, name_var, standard_name, data_1D_bin, unit, time)
    integer, intent(in) :: ncdf_id, z_size_id, bin_size_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: data_1D_bin

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,global_grid_z_size_netcdf,size_bin) :: data_1D_time
    data_1D_time(1,:,:) = data_1D_bin

    call check(nf90_def_var(ncdf_id, name_var, NF90_DOUBLE, (/ time_id, z_size_id, bin_size_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_1D_time))
  end subroutine define_and_write_variable_1D_bin





  subroutine define_and_write_variable_2D(ncdf_id, y_size_id, x_size_id, name_var, standard_name, data_2D, unit, time)
    integer, intent(in) :: ncdf_id, y_size_id, x_size_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: data_2D

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,global_grid_y_size_netcdf,global_grid_x_size_netcdf) :: data_2D_time
    data_2D_time(1,:,:) = data_2D

    call check(nf90_def_var(ncdf_id, name_var, NF90_DOUBLE, (/ time_id, y_size_id, x_size_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_2D_time))
  end subroutine define_and_write_variable_2D




  subroutine define_and_write_variable_3D(ncdf_id, z_size_id, y_size_id, x_size_id, name_var, standard_name, data_3D, unit, time)
    integer, intent(in) :: ncdf_id, z_size_id, y_size_id, x_size_id
    character(len = *), intent( in) :: name_var
    character(len = *), intent( in) :: standard_name
    character(len = *), intent( in) :: unit
    real(kind=DEFAULT_PRECISION), intent(in) :: time
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: data_3D

    integer :: var_id
    real(kind=DEFAULT_PRECISION), dimension(1,global_grid_z_size_netcdf,global_grid_y_size_netcdf, global_grid_x_size_netcdf) :: &
                                  data_3D_time
    data_3D_time(1,:,:,:) = data_3D

    call check(nf90_def_var(ncdf_id, name_var, NF90_DOUBLE, (/ time_id, z_size_id, y_size_id, x_size_id /), var_id))
    call check(nf90_put_att(ncdf_id, var_id, "standard_name", standard_name))
    call check(nf90_put_att(ncdf_id, var_id, "unit", unit))
    call check(nf90_put_var(ncdf_id, var_id , data_3D_time))
  end subroutine define_and_write_variable_3D




  subroutine check(status)
    integer(4), intent ( in) :: status
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine

end module io_server_netcdf_writer_mod
