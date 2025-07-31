










!> Add random noise into the fields
module randomnoise_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type, component_descriptor_type_v1
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
       options_get_logical_array, options_get_real_array, options_get_string_array, options_get_array_size, &
       options_compare_profile_arrays
  use interpolation_mod, only: piecewise_linear_1d
  use q_indices_mod, only: get_q_index, standard_q_names
  use logging_mod, only: LOG_ERROR, log_master_log
  use conversions_mod, only : conv_to_string

  implicit none

  private

  integer, parameter :: MAX_SIZE_SEED_ARRAY=256,  &  ! large value to work on multiple systems
                        I_SEED     = 7,           &  ! initial seed value,        non-reproducible case
                        THETA_SEED = -1731191804, &  ! initial seed for theta,    reproducible case
                        Q_SEED     =  2011234875, &  ! initial seed for q-fields, reproducible case
                        W_SEED     =  -163411914     ! initial seed for w,        reproducible case

  public initialisation_callback_randomnoise
contains

  !> The initialisation callback sets up the buoyancy coefficient
  !! @param current_state The current model state
  subroutine initialisation_callback_randomnoise(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer, dimension(MAX_SIZE_SEED_ARRAY) :: iranseed
    real(kind=DEFAULT_PRECISION) :: random_num

    integer :: nq_rand   ! The number of q fields to add noise to
    integer :: nzq       ! The number of input levels for noise
    integer :: i,j,k,n,s ! loop counters
    integer :: iq        ! temporary q varible index
    integer :: iloc, jloc
    integer :: iter, size_array = 0, iter2

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_rand_pl_q     ! Random Noise node height values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qv   ! Random Noise node amplitude for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_ql
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qr
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qi
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qs
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qg
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qAitkenSolMass
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qAccumSolMass
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qAccumInsolMass
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qCoarseSolMass
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_qCoarseDustMass
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_theta ! Random Noise node amplitude for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_rand_pl_theta ! Random Noise node height values for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_w     ! Random Noise node amplitude for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_rand_pl_w     ! Random Noise node height values for theta variable

    logical :: l_rand_pl_theta ! if .true. then random noise added to theta field
    logical :: l_rand_pl_q     ! if .true. then random noise added to q fields
    logical :: l_rand_pl_w     ! if .true. then random noise added to w field
    logical :: l_rand_bit_reproducible ! if .true. then is bit reproducible between runs
    logical :: logicnum

    character(len=STRING_LENGTH), dimension(:), allocatable :: names_rand_pl_q ! names of q variables to add random noise to

    real(kind=DEFAULT_PRECISION), allocatable :: f_rand_pl_q_tmp(:) !temporary 1D storage of random noise for q field
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid to use in interpolation

    if (current_state%continuation_run) return

    allocate(zgrid(current_state%local_grid%local_domain_end_index(Z_INDEX)))

    !l_rand_pl_theta=options_get_logical(current_state%options_database, "l_rand_pl_theta")
    !l_rand_pl_w=options_get_logical(current_state%options_database, "l_rand_pl_w")
    !l_rand_pl_q=options_get_logical(current_state%options_database, "l_rand_pl_q")
    !l_rand_bit_reproducible=options_get_logical(current_state%options_database, "l_rand_bit_reproducible")

    do iter = 1,current_state%config_args
      if (current_state%options_database_string(iter,1) .eq. "l_rand_pl_theta") then
        read(current_state%options_database_string(iter,2),*) logicnum
        l_rand_pl_theta = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "l_rand_pl_w") then
        read(current_state%options_database_string(iter,2),*) logicnum
        l_rand_pl_w = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "l_rand_pl_q") then
        read(current_state%options_database_string(iter,2),*) logicnum
        l_rand_pl_q = logicnum
      else if (current_state%options_database_string(iter,1) .eq. "l_rand_bit_reproducible") then
        read(current_state%options_database_string(iter,2),*) logicnum
        l_rand_bit_reproducible = logicnum
      end if
    end do

    ! Initialise random seed to be different in every MONC
    !  This is overwritten if l_rand_bit_reproducible=.true.
    iranseed = I_SEED + current_state%parallel%my_rank
    call random_seed(put = iranseed)

    !---------------------------------------------------
    ! Handle randomnoise for theta
    if (l_rand_pl_theta)then
      do iter = 1,current_state%config_args
        size_array = 0
        if (current_state%options_database_string(iter,1) .eq. "z_rand_pl_theta") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          z_rand_pl_theta = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_theta") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_theta = current_state%options_database_real(iter,1:size_array)
        end if
      end do


      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_rand_pl_theta(1:size(z_rand_pl_theta)), f_rand_pl_theta(1:size(f_rand_pl_theta)), zgrid, &
           current_state%global_grid%configuration%vertical%theta_rand)

      ! Cycle through all points in the global grid, record if on the local grid
      if (l_rand_bit_reproducible) then
        iranseed(:) = THETA_SEED
        call random_seed(put = iranseed)
        do i=1, current_state%global_grid%size(X_INDEX)
          iloc = i - current_state%local_grid%start(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) + 1
          do j=1, current_state%global_grid%size(Y_INDEX)
            jloc = j - current_state%local_grid%start(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) + 1
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
              call random_number(random_num)

              if (iloc .ge. current_state%local_grid%local_domain_start_index(X_INDEX) .and. &
                  iloc .le. current_state%local_grid%local_domain_end_index(X_INDEX)   .and. &
                  jloc .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) .and. &
                  jloc .le. current_state%local_grid%local_domain_end_index(Y_INDEX) ) then
                  current_state%th%data(k,jloc,iloc) = current_state%th%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%theta_rand(k) * 2.0 * (random_num-0.5)

              end if ! in local space?
            end do ! k
          end do ! j
        end do ! i

      else ! not l_rand_bit_reproducible
            
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
  
              ! Apply random number for this grid point
              !   if not l_rand_bit_reproducible, then random_num is determined by local rank
              call random_number(random_num)
              current_state%th%data(k,j,i) = current_state%th%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%theta_rand(k) * 2.0 * (random_num-0.5)
             end do
          end do
        end do
      end if !check l_rand_bit_reproducible
      deallocate(z_rand_pl_theta, f_rand_pl_theta)
    end if ! l_rand_pl_theta



    !---------------------------------------------------
    ! Handle randomnoise for q-fields
    if (l_rand_pl_q)then
      do iter = 1,current_state%config_args
        size_array = 0
        if (current_state%options_database_string(iter,1) .eq. "z_rand_pl_q") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          z_rand_pl_q = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qv") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qv = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_ql") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_ql = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qr") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qr = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qi") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qi = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qs") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qs = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qg") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qg = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qAitkenSolMass") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qAitkenSolMass = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qAccumSolMass") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qAccumSolMass = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qAccumInsolMass") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qAccumInsolMass = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qCoarseSolMass") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qCoarseSolMass = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_qCoarseDustMass") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_qCoarseDustMass = current_state%options_database_real(iter,1:size_array)
        end if
      end do
    !  nq_rand=size(names_rand_pl_q)
    !  allocate(z_rand_pl_q(options_get_array_size(current_state%options_database, "z_rand_pl_q")))
    !  call options_get_real_array(current_state%options_database, "z_rand_pl_q", z_rand_pl_q)
    !  nzq=size(z_rand_pl_q)
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
    !  allocate(f_rand_pl_q_tmp(options_get_array_size(current_state%options_database, "f_rand_pl_q")))
    !  if (nq_rand*nzq .ne. size(f_rand_pl_q_tmp)) then
    !    call log_master_log(LOG_ERROR, "There is a mismatch between the number of moisture perturbation heights, "//&
    !                                   "size(z_rand_pl_q)="//trim(conv_to_string(nzq))//                            &
    !                                   ", and the perturbation values, "//                                               &
    !                                   "size(f_rand_pl_q)="//trim(conv_to_string(size(f_rand_pl_q_tmp)))//          &
    !                                   ".  The length of f_rand_pl_q should equal the length of z_rand_pl_q "//     &
    !                                   "multiplied by the number of names_rand_pl_q.")
    !  end if
    !  call options_get_real_array(current_state%options_database, "f_rand_pl_q", f_rand_pl_q_tmp)
    !  allocate(f_rand_pl_q(nzq, nq_rand))
    !  f_rand_pl_q(1:nzq, 1:nq_rand)=reshape(f_rand_pl_q_tmp, (/nzq, nq_rand/))
    !  do n=1,nq_rand

    !    iq=get_q_index(trim(names_rand_pl_q(n)), 'random noise')
    !    zgrid=current_state%global_grid%configuration%vertical%zn(:)

      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_qv(1:size(f_rand_pl_qv)), zgrid, &
                               current_state%global_grid%configuration%vertical%qv_rand)
      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_ql(1:size(f_rand_pl_ql)), zgrid, &
                                current_state%global_grid%configuration%vertical%ql_rand)
      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_qr(1:size(f_rand_pl_qr)), zgrid, &
                                current_state%global_grid%configuration%vertical%qr_rand)
      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_qi(1:size(f_rand_pl_qi)), zgrid, &
                                current_state%global_grid%configuration%vertical%qi_rand)
      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_qs(1:size(f_rand_pl_qs)), zgrid, &
                                current_state%global_grid%configuration%vertical%qs_rand)
      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_qAitkenSolMass(1:size(f_rand_pl_qAitkenSolMass)), &
                                zgrid, current_state%global_grid%configuration%vertical%qAitkenSolMass_rand)
      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_qAccumSolMass(1:size(f_rand_pl_qAccumSolMass)), &
                                zgrid, current_state%global_grid%configuration%vertical%qAccumSolMass_rand)
      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_qAccumInsolMass(1:size(f_rand_pl_qAccumInsolMass)), &
                                zgrid, current_state%global_grid%configuration%vertical%qAccumInsolMass_rand)
      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_qCoarseSolMass(1:size(f_rand_pl_qCoarseSolMass)), &
                                zgrid, current_state%global_grid%configuration%vertical%qCoarseSolMass_rand)
      call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_qCoarseDustMass(1:size(f_rand_pl_qCoarseDustMass)), &
                                zgrid, current_state%global_grid%configuration%vertical%qCoarseDustMass_rand)


        ! Cycle through all points in the global grid, record if on the local grid
      if (l_rand_bit_reproducible) then
        iranseed(:) = Q_SEED
        call random_seed(put = iranseed)
        do i=1, current_state%global_grid%size(X_INDEX)
          iloc = i - current_state%local_grid%start(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) + 1
          do j=1, current_state%global_grid%size(Y_INDEX)
            jloc = j - current_state%local_grid%start(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) + 1
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
              call random_number(random_num)

              if (iloc .ge. current_state%local_grid%local_domain_start_index(X_INDEX) .and. &
                  iloc .le. current_state%local_grid%local_domain_end_index(X_INDEX)   .and. &
                  jloc .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) .and. &
                  jloc .le. current_state%local_grid%local_domain_end_index(Y_INDEX) ) then
                  current_state%qv%data(k,jloc,iloc) = current_state%qv%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qv_rand(k) * 2.0 * (random_num-0.5)
                  current_state%ql%data(k,jloc,iloc) = current_state%ql%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%ql_rand(k) * 2.0 * (random_num-0.5)
                  current_state%qr%data(k,jloc,iloc) = current_state%qr%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qr_rand(k) * 2.0 * (random_num-0.5)
                  current_state%qi%data(k,jloc,iloc) = current_state%qi%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qi_rand(k) * 2.0 * (random_num-0.5)
                  current_state%qs%data(k,jloc,iloc) = current_state%qs%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qs_rand(k) * 2.0 * (random_num-0.5)
                  current_state%qg%data(k,jloc,iloc) = current_state%qg%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qg_rand(k) * 2.0 * (random_num-0.5)

                  current_state%qAitkenSolMass%data(k,jloc,iloc) = current_state%qAitkenSolMass%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qAitkenSolMass_rand(k) * 2.0 * (random_num-0.5)
                  current_state%qAccumSolMass%data(k,jloc,iloc) = current_state%qAccumSolMass%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qAccumSolMass_rand(k) * 2.0 * (random_num-0.5)
                  current_state%qAccumInsolMass%data(k,jloc,iloc) = current_state%qAccumInsolMass%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qAccumInsolMass_rand(k) * 2.0 * (random_num-0.5)
                  current_state%qCoarseSolMass%data(k,jloc,iloc) = current_state%qCoarseSolMass%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qCoarseSolMass_rand(k) * 2.0 * (random_num-0.5)
                  current_state%qCoarseDustMass%data(k,jloc,iloc) = current_state%qCoarseDustMass%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%qCoarseDustMass_rand(k) * 2.0 * (random_num-0.5)

              end if ! in local space?
            end do ! k
          end do ! j
        end do ! i

      else ! not l_rand_bit_reproducible
  
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), &
               current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), &
                 current_state%local_grid%local_domain_end_index(Y_INDEX)
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
                ! Apply random number for this grid point
                !   if not l_rand_bit_reproducible, then random_num is determined by local rank
              call random_number(random_num)
              current_state%qv%data(k,j,i) = current_state%qv%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qv_rand(k) * 2.0 * (random_num-0.5)
              current_state%ql%data(k,j,i) = current_state%ql%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%ql_rand(k) * 2.0 * (random_num-0.5)
              current_state%qr%data(k,j,i) = current_state%qr%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qr_rand(k) * 2.0 * (random_num-0.5)
              current_state%qi%data(k,j,i) = current_state%qi%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qi_rand(k) * 2.0 * (random_num-0.5)
              current_state%qs%data(k,j,i) = current_state%qs%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qs_rand(k) * 2.0 * (random_num-0.5)
              current_state%qg%data(k,j,i) = current_state%qg%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qg_rand(k) * 2.0 * (random_num-0.5)

              current_state%qAitkenSolMass%data(k,j,i) = current_state%qAitkenSolMass%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qAitkenSolMass_rand(k) * 2.0 * (random_num-0.5)
              current_state%qAccumSolMass%data(k,j,i) = current_state%qAccumSolMass%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qAccumSolMass_rand(k) * 2.0 * (random_num-0.5)
              current_state%qAccumInsolMass%data(k,j,i) = current_state%qAccumInsolMass%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qAccumInsolMass_rand(k) * 2.0 * (random_num-0.5)
              current_state%qCoarseSolMass%data(k,j,i) = current_state%qCoarseSolMass%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qCoarseSolMass_rand(k) * 2.0 * (random_num-0.5)
              current_state%qCoarseDustMass%data(k,j,i) = current_state%qCoarseDustMass%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%qCoarseDustMass_rand(k) * 2.0 * (random_num-0.5)
            end do ! k
          end do ! j
        end do ! i
      end if ! l_rand_bit_reproducible
    !  end do ! n
      !deallocate(z_rand_pl_q, f_rand_pl_q_tmp, f_rand_pl_q, names_rand_pl_q)
      deallocate(z_rand_pl_q, f_rand_pl_qv, f_rand_pl_ql, f_rand_pl_qr, f_rand_pl_qi, f_rand_pl_qs, &
                 f_rand_pl_qg, f_rand_pl_qAitkenSolMass, f_rand_pl_qAccumSolMass, f_rand_pl_qAccumInsolMass, &
                 f_rand_pl_qCoarseSolMass, f_rand_pl_qCoarseDustMass)
    end if ! l_rand_pl_q


    !---------------------------------------------------
    ! Handle randomnoise for w
    if (l_rand_pl_w)then
      do iter = 1,current_state%config_args
        size_array = 0
        if (current_state%options_database_string(iter,1) .eq. "z_rand_pl_w") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          z_rand_pl_w = current_state%options_database_real(iter,1:size_array)
        else if (current_state%options_database_string(iter,1) .eq. "f_rand_pl_w") then
          do iter2 = 1,size(current_state%options_database_real(iter,:))
            if (isnan(current_state%options_database_real(iter,iter2))) then
              exit
            else
              size_array = size_array + 1
            end if
          end do
          f_rand_pl_w = current_state%options_database_real(iter,1:size_array)
        end if
      end do

      ! Get amplitude profiles
    !  allocate(z_rand_pl_w(options_get_array_size(current_state%options_database, "z_rand_pl_w")), &
    !       f_rand_pl_w(options_get_array_size(current_state%options_database, "f_rand_pl_w")))
    !  call options_compare_profile_arrays(current_state%options_database, &
    !                            "z_rand_pl_w", "f_rand_pl_w", "w perturbation")
    !  call options_get_real_array(current_state%options_database, "z_rand_pl_w", z_rand_pl_w)
    !  call options_get_real_array(current_state%options_database, "f_rand_pl_w", f_rand_pl_w)

      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_rand_pl_w(1:size(z_rand_pl_w)), f_rand_pl_w(1:size(f_rand_pl_w)), zgrid, &
           current_state%global_grid%configuration%vertical%w_rand)

      ! Cycle through all points in the global grid, record if on the local grid
      if (l_rand_bit_reproducible) then
        iranseed(:) = W_SEED
        call random_seed(put = iranseed)
        do i=1, current_state%global_grid%size(X_INDEX)
          iloc = i - current_state%local_grid%start(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) + 1
          do j=1, current_state%global_grid%size(Y_INDEX)
            jloc = j - current_state%local_grid%start(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) + 1
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
              call random_number(random_num)

              if (iloc .ge. current_state%local_grid%local_domain_start_index(X_INDEX) .and. &
                  iloc .le. current_state%local_grid%local_domain_end_index(X_INDEX)   .and. &
                  jloc .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) .and. &
                  jloc .le. current_state%local_grid%local_domain_end_index(Y_INDEX) ) then
                current_state%w%data(k,jloc,iloc) = current_state%w%data(k,jloc,iloc) + &
                   current_state%global_grid%configuration%vertical%w_rand(k) * (random_num-0.5)
              end if ! in local space?
            end do ! k
          end do ! j
        end do ! i
      else ! not l_rand_bit_reproducible

        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)

              ! Apply random number for this grid point
              !   if not l_rand_bit_reproducible, then random_num is determined by local rank
              call random_number(random_num)
              current_state%w%data(k,j,i) = current_state%w%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%w_rand(k) * (random_num-0.5)
            end do
            current_state%w%data(current_state%local_grid%local_domain_end_index(Z_INDEX),j,i)=0.0_DEFAULT_PRECISION
            current_state%w%data(1,j,i)=0.0_DEFAULT_PRECISION
            if (current_state%use_viscosity_and_diffusion) then
              current_state%u%data(1,j,i)=-current_state%u%data(2,j,i)
              current_state%v%data(1,j,i)=-current_state%v%data(2,j,i)
            else
              current_state%u%data(1,j,i)=current_state%u%data(2,j,i)
              current_state%v%data(1,j,i)=current_state%v%data(2,j,i)
            end if
          end do ! j
        end do ! i
      end if ! l_rand_bit_reproducible
      deallocate(z_rand_pl_w, f_rand_pl_w)
    end if ! l_rand_pl_w
    deallocate(zgrid)
  end subroutine initialisation_callback_randomnoise
end module randomnoise_mod
