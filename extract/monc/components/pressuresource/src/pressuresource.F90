!> Calculates the gradient of the source flow fields (SU, SV, SW.) This is based upon the P field values already
!! set for the divergence error in the diverr component.
!! Note that some communication is required, as a process needs to know the x and y -1 values, which are held on the
!! -1 neighbour in that dimension. These are computed here and sent to neighbour +1, the recvs from a process to neighbour -1
!! are also registered here and the handles are waited upon the the results combined into P in the solver. If the neighbour
!! is local in any dimension then it is just a local memory copy
module pressuresource_mod
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use state_mod, only : model_state_type
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use logging_mod, only : LOG_ERROR, log_master_log
  use mpi, only : MPI_REQUEST_NULL
  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: send_buffer_x, send_buffer_y

  public initialisation_callback_pressuresource, timestep_callback_pressuresource
contains

  !> On initialisation this will allocate the buffer areas required and set communication handles to null
  !! @param current_state The current model state
  subroutine initialisation_callback_pressuresource(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: iter
    logical :: logicnum

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "diverr_enabled") then
        read(current_state%options_database_string(iter,2),*) logicnum
        if (.not. logicnum) then
          call log_master_log(LOG_ERROR, "The pressure source component requires the diverr component to be enabled")
        end if
      end if
    end do

#ifdef U_ACTIVE
    allocate(send_buffer_x(current_state%local_grid%size(Z_INDEX)-1, current_state%local_grid%size(Y_INDEX)), &
         current_state%psrce_recv_buffer_x(current_state%local_grid%size(Z_INDEX)-1, current_state%local_grid%size(Y_INDEX)))
#endif
#ifdef V_ACTIVE
    allocate(send_buffer_y(current_state%local_grid%size(Z_INDEX)-1, current_state%local_grid%size(X_INDEX)), &
         current_state%psrce_recv_buffer_y(current_state%local_grid%size(Z_INDEX)-1, current_state%local_grid%size(X_INDEX)))
#endif
    current_state%psrce_x_hs_send_request=MPI_REQUEST_NULL
    current_state%psrce_y_hs_send_request=MPI_REQUEST_NULL
    current_state%psrce_x_hs_recv_request=MPI_REQUEST_NULL
    current_state%psrce_y_hs_recv_request=MPI_REQUEST_NULL
  end subroutine initialisation_callback_pressuresource

  !> The timestep callback will update the values of P for each column
  !! @param current_state The current model state
  subroutine timestep_callback_pressuresource(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (current_state%first_timestep_column) call register_neighbouring_pressure_data_recv(current_state)
    if (.not. current_state%halo_column) call calculate_psrce(current_state)
    if (current_state%last_timestep_column) call send_neighbouring_pressure_data(current_state)
  end subroutine timestep_callback_pressuresource

  !> Frees up the allocated buffers (if such were allocated)
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(send_buffer_x)) deallocate(send_buffer_x)
    if (allocated(current_state%psrce_recv_buffer_x)) deallocate(current_state%psrce_recv_buffer_x)
    if (allocated(send_buffer_y)) deallocate(send_buffer_y)
    if (allocated(current_state%psrce_recv_buffer_y)) deallocate(current_state%psrce_recv_buffer_y)
  end subroutine finalisation_callback      

  !> Combines the source fields with the pressure values. For U and V, if this is on the low boundary then delay dealing with the
  !! -1 values as they will be communicated. Equally, for those fields if this is the high boundary then compute the values
  !! and send them on to the neighbouring process.
  !! @param current_state The current model state
  subroutine calculate_psrce(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: k, local_y, local_x, corrected_y, corrected_x
    logical :: last_x, last_y

    local_y=current_state%column_local_y
    local_x=current_state%column_local_x
    if (local_y .eq. current_state%local_grid%local_domain_end_index(Y_INDEX)) last_y=.true.
    if (local_x .eq. current_state%local_grid%local_domain_end_index(X_INDEX)) last_x=.true.
    !last_y = local_y == current_state%local_grid%local_domain_end_index(Y_INDEX)
    !last_x = local_x == current_state%local_grid%local_domain_end_index(X_INDEX)
    if (last_x .or. last_y) then
      corrected_x=local_x-current_state%local_grid%halo_size(X_INDEX)
      corrected_y=local_y-current_state%local_grid%halo_size(Y_INDEX)
    end if

    if(current_state%immersed%ib_enabled) then
      if (current_state%immersed%ib_col(local_y, local_x)) then
      ! IB enabled
      do k=2,current_state%local_grid%size(Z_INDEX)
#ifdef W_ACTIVE
        if (current_state%immersed%indic_s(k,local_y,local_x).eq.0)then
        current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)+&
             4.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc2(k)*&
             current_state%sw%data(k, local_y, local_x)-&
             current_state%global_grid%configuration%vertical%tzc1(k)*current_state%sw%data(k-1, local_y, local_x))
         end if
#endif
#ifdef U_ACTIVE
        if (current_state%immersed%indic_s(k,local_y,local_x).eq.0)then
        current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)+&
             current_state%global_grid%configuration%horizontal%cx * current_state%su%data(k, local_y, local_x)
         end if
#endif
#ifdef V_ACTIVE
        if (current_state%immersed%indic_s(k,local_y,local_x).eq.0)then
        current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)+&
             current_state%global_grid%configuration%horizontal%cy * current_state%sv%data(k, local_y, local_x)
         end if
#endif
#ifdef U_ACTIVE
        if (local_x .gt. 3) then
          if (current_state%immersed%indic_s(k,local_y,local_x).eq.0)then
            current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)-&
               current_state%global_grid%configuration%horizontal%cx * current_state%su%data(k, local_y, local_x-1)
          end if
        end if

        if (last_x) then
          send_buffer_x(k-1, corrected_y)=&
               current_state%global_grid%configuration%horizontal%cx * current_state%su%data(k, local_y, local_x)
        end if
#endif
#ifdef V_ACTIVE
        if (local_y .gt. 3 .and. local_x .gt. 3) then
          if (current_state%immersed%indic_s(k,local_y,local_x).eq.0)then
            current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)-&
               current_state%global_grid%configuration%horizontal%cy * current_state%sv%data(k, local_y-1, local_x)
          end if
        end if
        if (last_y) then
          send_buffer_y(k-1, corrected_x)=&
               current_state%global_grid%configuration%horizontal%cy * current_state%sv%data(k, local_y, local_x)
        end if
#endif
      end do
      end if
    else
      ! No IB
      do k=2,current_state%local_grid%size(Z_INDEX)
#ifdef W_ACTIVE
        current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)+&
             4.0_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%tzc2(k)*&
             current_state%sw%data(k, local_y, local_x)-&
             current_state%global_grid%configuration%vertical%tzc1(k)*current_state%sw%data(k-1, local_y, local_x))
#endif
#ifdef U_ACTIVE
        current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)+&
             current_state%global_grid%configuration%horizontal%cx * current_state%su%data(k, local_y, local_x)
#endif
#ifdef V_ACTIVE
        current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)+&
             current_state%global_grid%configuration%horizontal%cy * current_state%sv%data(k, local_y, local_x)
#endif
#ifdef U_ACTIVE
        if (local_x .gt. 3) then
          current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)-&
               current_state%global_grid%configuration%horizontal%cx * current_state%su%data(k, local_y, local_x-1)
        end if

        if (last_x) then
          send_buffer_x(k-1, corrected_y)=&
               current_state%global_grid%configuration%horizontal%cx * current_state%su%data(k, local_y, local_x)
        end if
#endif
#ifdef V_ACTIVE
        if (local_y .gt. 3 .and. local_x .gt. 3) then
          current_state%p%data(k, local_y, local_x)=current_state%p%data(k, local_y, local_x)-&
               current_state%global_grid%configuration%horizontal%cy * current_state%sv%data(k, local_y-1, local_x)
        end if
        if (last_y) then
          send_buffer_y(k-1, corrected_x)=&
               current_state%global_grid%configuration%horizontal%cy * current_state%sv%data(k, local_y, local_x)
        end if
#endif

      end do

    end if
  end subroutine calculate_psrce

  !> Sends the computed source pressure data terms to the p+1 process
  !! @param current_state The current model state
  subroutine send_neighbouring_pressure_data(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: ierr

#ifdef U_ACTIVE
    if (current_state%local_grid%neighbours(X_INDEX,3) .eq. current_state%parallel%my_rank) then
      current_state%psrce_recv_buffer_x=send_buffer_x
    else
      call mpi_isend(send_buffer_x, size(send_buffer_x), PRECISION_TYPE, current_state%local_grid%neighbours(X_INDEX,3), &
           10, current_state%parallel%neighbour_comm, current_state%psrce_x_hs_send_request, ierr)
    end if
#endif
#ifdef V_ACTIVE
    if (current_state%local_grid%neighbours(Y_INDEX,3) .eq. current_state%parallel%my_rank) then
      current_state%psrce_recv_buffer_y=send_buffer_y
    else
      call mpi_isend(send_buffer_y, size(send_buffer_y), PRECISION_TYPE, current_state%local_grid%neighbours(Y_INDEX,3), &
           10, current_state%parallel%neighbour_comm, current_state%psrce_y_hs_send_request, ierr)
    end if
#endif
  end subroutine send_neighbouring_pressure_data

  !> Registers the receive requests for each neighbouring process if that is not local, is recieves from p-1
  !! @param current_state The current model state
  subroutine register_neighbouring_pressure_data_recv(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: ierr

#ifdef U_ACTIVE
    if (current_state%local_grid%neighbours(X_INDEX,2) .ne. current_state%parallel%my_rank) then
      call mpi_irecv(current_state%psrce_recv_buffer_x, size(current_state%psrce_recv_buffer_x), PRECISION_TYPE, &
           current_state%local_grid%neighbours(X_INDEX,2), 10, current_state%parallel%neighbour_comm, &
           current_state%psrce_x_hs_recv_request, ierr)
    end if
#endif
#ifdef V_ACTIVE
    if (current_state%local_grid%neighbours(Y_INDEX,2) .ne. current_state%parallel%my_rank) then
      call mpi_irecv(current_state%psrce_recv_buffer_y, size(current_state%psrce_recv_buffer_y), PRECISION_TYPE, &
           current_state%local_grid%neighbours(Y_INDEX,2), 10, current_state%parallel%neighbour_comm, &
           current_state%psrce_y_hs_recv_request, ierr)
    end if
#endif
  end subroutine register_neighbouring_pressure_data_recv
end module pressuresource_mod
