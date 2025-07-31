!> Conditionally averaged diagnostics, Part 2 of 2.
!  Using the totalled values from Part 1, averages the results over the horizontal domain,
!  creating profiles of the requested diagnostics under the requested conditions and .not. conditions.
!  The "area" diagnostic is a special case, recording the fraction of the domain meeting the associated
!  condition and .not. condition.
module conditional_diagnostics_whole_mod
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use state_mod, only : model_state_type

  use conditional_diagnostics_column_mod, only : CondDiags_tot, ncond, ndiag, gpts_total, requested_area
  use grids_mod, only : Z_INDEX
  use datadefn_mod, only : PRECISION_TYPE, DEFAULT_PRECISION
  use mpi, only : MPI_SUM, MPI_IN_PLACE, MPI_INT, MPI_REAL, MPI_DOUBLE
  use missing_data_mod, only: rmdi
  use optionsdatabase_mod, only : options_get_integer

  implicit none

#ifndef TEST_MODE
  private
#endif

  public initialisation_callback_conditional_diagnostics_whole, &
        timestep_callback_conditional_diagnostics_whole, finalisation_callback_conditional_diagnostics_whole

contains

  !> Initialisation hook: currently doesn't need to do anything
  !! @param current_state The current model state
  subroutine initialisation_callback_conditional_diagnostics_whole(current_state)
    type(model_state_type), target, intent(inout) :: current_state

  end subroutine initialisation_callback_conditional_diagnostics_whole


  !> The timestep hook will perform averaging of the conditional diagnostics
  !! @param current_state The current model state
  subroutine timestep_callback_conditional_diagnostics_whole(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: k,cnc,dnc,ierr
    real(kind=DEFAULT_PRECISION) :: temp
    logical :: calculate_diagnostics

    !calculate_diagnostics = current_state%diagnostic_sample_timestep

    !> Decide if conditions are appropriate to proceed with calculations
    if (current_state%modulo_number_1d .ne. 0) return

    !> Sum conditional diagnostics total array (horizontally), placing the result on process 0 of monc processes which is fact equal to 1
    !! Reduction call on process 0 requires special MPI_IN_PLACE handling
    if (current_state%parallel%my_rank == 0) then
      call mpi_reduce(MPI_IN_PLACE , current_state%w_zn_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%w_zn2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%tmp_th_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wth_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wth_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%thv_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wthv_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%th_pr2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wthsg_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%w_zn3_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%relhum_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%tmp_u_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%tmp_v_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wu_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wv_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wusg_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wvsg_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%TdegK_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%th_h_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%th_h_pr1_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%th_h_pr2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%qvli_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%qvli_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%qvli_pr2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%qppt_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%qppt_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%qppt_pr2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wqvli_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(MPI_IN_PLACE , current_state%wqppt_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
    else
      call mpi_reduce(current_state%w_zn_cln, current_state%w_zn_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%w_zn2_cln, current_state%w_zn2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%tmp_th_cln, current_state%tmp_th_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wth_cln, current_state%wth_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wth_pr_cln, current_state%wth_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%thv_pr_cln, current_state%thv_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wthv_pr_cln, current_state%wthv_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%th_pr2_cln, current_state%th_pr2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wthsg_cln, current_state%wthsg_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%w_zn3_cln, current_state%w_zn3_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%relhum_cln, current_state%relhum_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%tmp_u_cln, current_state%tmp_u_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%tmp_v_cln, current_state%tmp_v_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wu_cln, current_state%wu_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wv_cln, current_state%wv_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wusg_cln, current_state%wusg_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wvsg_cln, current_state%wvsg_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%TdegK_cln, current_state%TdegK_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%th_h_cln, current_state%th_h_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%th_h_pr1_cln, current_state%th_h_pr1_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%th_h_pr2_cln, current_state%th_h_pr2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%qvli_cln, current_state%qvli_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%qvli_pr_cln, current_state%qvli_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%qvli_pr2_cln, current_state%qvli_pr2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%qppt_cln, current_state%qppt_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%qppt_pr_cln, current_state%qppt_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%qppt_pr2_cln, current_state%qppt_pr2_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wqvli_pr_cln, current_state%wqvli_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)
      call mpi_reduce(current_state%wqppt_pr_cln, current_state%wqppt_pr_cln, current_state%local_grid%size(Z_INDEX), &
                      PRECISION_TYPE, MPI_SUM, 0, current_state%parallel%monc_communicator, ierr)

    end if

    !> Average summed diagnostics over the domain by dividing the total diagnostic for each condition 
    !! by the total number of points for the associated conditions.
    !! This is NOT done for the area diagnostic, identified by the requested_area position in the array.
    !! Note: CondDiags_tot(k, ncond, ndiag)
    !do dnc = 1, ndiag
    !  if (dnc /= requested_area) then
    !    do cnc = 1, ncond
    !      do k = 2, current_state%local_grid%size(Z_INDEX) - 1
    !        temp = CondDiags_tot(k,cnc,requested_area)
    !        if (temp .gt. 0) then
    !          CondDiags_tot(k,cnc,dnc) = CondDiags_tot(k,cnc,dnc) / temp
    !        else
    !          CondDiags_tot(k,cnc,dnc) = rmdi
    !        end if
     !     end do ! k
    !    end do ! cnc over ncond*2
    !  end if ! check requested_area
    !end do ! dnc over ndiag
    if (current_state%parallel%my_rank == 0) then
      do k = 1, current_state%local_grid%size(Z_INDEX)

        current_state%w_zn_cln(k) = current_state%w_zn_cln(k)/gpts_total
        current_state%w_zn2_cln(k) = current_state%w_zn2_cln(k)/gpts_total
        current_state%tmp_th_cln(k) = current_state%tmp_th_cln(k)/gpts_total
        current_state%wth_cln(k) = current_state%wth_cln(k)/gpts_total
        current_state%th_pr_cln(k) = current_state%th_pr_cln(k)/gpts_total
        current_state%wth_pr_cln(k) = current_state%wth_pr_cln(k)/gpts_total
        current_state%thv_pr_cln(k) = current_state%thv_pr_cln(k)/gpts_total
        current_state%wthv_pr_cln(k) = current_state%wthv_pr_cln(k)/gpts_total
        current_state%th_pr2_cln(k) = current_state%th_pr2_cln(k)/gpts_total
        current_state%wthsg_cln(k) = current_state%wthsg_cln(k)/gpts_total
        current_state%w_zn3_cln(k) = current_state%w_zn3_cln(k)/gpts_total
        current_state%relhum_cln(k) = current_state%relhum_cln(k)/gpts_total
        current_state%tmp_u_cln(k) = current_state%tmp_u_cln(k)/gpts_total
        current_state%tmp_v_cln(k) = current_state%tmp_v_cln(k)/gpts_total
        current_state%wu_cln(k) = current_state%wu_cln(k)/gpts_total
        current_state%wv_cln(k) = current_state%wv_cln(k)/gpts_total
        current_state%wusg_cln(k) = current_state%wusg_cln(k)/gpts_total
        current_state%wvsg_cln(k) = current_state%wvsg_cln(k)/gpts_total
        current_state%TdegK_cln(k) = current_state%TdegK_cln(k)/gpts_total
        current_state%th_h_cln(k) = current_state%th_h_cln(k)/gpts_total
        current_state%th_h_pr1_cln(k) = current_state%th_h_pr1_cln(k)/gpts_total
        current_state%th_h_pr2_cln(k) = current_state%th_h_pr2_cln(k)/gpts_total
        current_state%qvli_cln(k) = current_state%qvli_cln(k)/gpts_total
        current_state%qvli_pr_cln(k) = current_state%qvli_pr_cln(k)/gpts_total
        current_state%qvli_pr2_cln(k) = current_state%qvli_pr2_cln(k)/gpts_total
        current_state%qppt_cln(k) = current_state%qppt_cln(k)/gpts_total
        current_state%qppt_pr_cln(k) = current_state%qppt_pr_cln(k)/gpts_total
        current_state%qppt_pr2_cln(k) = current_state%qppt_pr2_cln(k)/gpts_total
        current_state%wqvli_pr_cln(k) = current_state%wqvli_pr_cln(k)/gpts_total
        current_state%wqppt_pr_cln(k) = current_state%wqppt_pr_cln(k)/gpts_total
      end do
      !print*,"current_state%tmp_th_cln = ", current_state%tmp_th_cln
    end if

    !> Convert the area total number of points for each condition to fraction of the horizontal domain.
    !CondDiags_tot(:,:,requested_area) = CondDiags_tot(:,:,requested_area) / gpts_total

    !> Apply missing data mask to top/bottom
    !CondDiags_tot(1,:,:) = rmdi
    !CondDiags_tot(current_state%local_grid%size(Z_INDEX),:,:) = rmdi

    !> Since the xml handling of CondDiags_tot will perform a sum over processes, divide by the number
    !! of proceses.
    !CondDiags_tot = CondDiags_tot / current_state%parallel%processes

    !> Broadcast the process-fractional solution to all processes.
    call mpi_bcast(current_state%w_zn_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%w_zn2_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%tmp_th_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wth_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%th_pr_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wth_pr_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%thv_pr_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wthv_pr_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%th_pr2_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wthsg_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%w_zn3_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%relhum_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%tmp_u_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%tmp_v_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wu_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wv_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wusg_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wvsg_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%TdegK_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%th_h_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%th_h_pr1_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%th_h_pr2_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%qvli_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%qvli_pr_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%qvli_pr2_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%qppt_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%qppt_pr_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%qppt_pr2_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wqvli_pr_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(current_state%wqppt_pr_cln,  current_state%local_grid%size(Z_INDEX), &
                   PRECISION_TYPE, 0, current_state%parallel%monc_communicator, ierr)

  end subroutine timestep_callback_conditional_diagnostics_whole

  
  !> Called on termination: currently doesn't need to do anything
  !! @param current_state The current model state
  subroutine finalisation_callback_conditional_diagnostics_whole(current_state)
    type(model_state_type), target, intent(inout) :: current_state

  end subroutine finalisation_callback_conditional_diagnostics_whole

end module conditional_diagnostics_whole_mod
