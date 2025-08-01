










module forthread_types
  use iso_c_binding
  implicit none

  type, bind(c) :: sched_param
     integer(c_int) :: sched_priority
  end type sched_param

  ! wrapping the C timespec type - maybe we can
  ! skip this and only use primite fortran types
  type, bind(c) :: timespec
     integer(c_int)  :: tv_sec  ! seconds
     integer(c_long) :: tv_nsec ! nanoseconds
  end type timespec

  integer, parameter :: size_t = c_size_t
end module forthread_types
