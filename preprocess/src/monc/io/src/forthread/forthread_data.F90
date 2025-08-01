










module forthread_data
  use iso_c_binding
  use forthread_ciface_mod
  implicit none

  abstract interface 
     subroutine i_run(arg)
       integer :: arg
     end subroutine i_run
  end interface

  abstract interface 
     subroutine i_once() bind(c)
     end subroutine i_once
  end interface

  abstract interface 
     subroutine i_destructor(arg) bind(c)
       use iso_c_binding
       integer(c_int) :: arg
     end subroutine i_destructor
  end interface

  type t_run
     procedure(i_run), pointer, nopass :: run
     integer, pointer :: arg
  end type t_run

  type ptr_t_run
     type(t_run), pointer :: t => NULL()
  end type ptr_t_run

  integer, parameter                       :: init_size = 16
  type(ptr_t_run), dimension(:), pointer   :: routine_table => NULL()
  integer                                  :: routine_table_size
  integer                                  :: routine_table_mutex
contains

  ! the return value still needs work to be done
  ! currently the start_routine cannot return any value.
  ! This should be handled similarly to thread_exit
  type(c_ptr) function start_routine(arg) bind(c)

    type(c_ptr), value, intent(in)  :: arg
    !returns?

    type(t_run), pointer :: exec
    integer,pointer      :: ret

    call c_f_pointer(arg,exec)
    call exec%run(exec%arg)
    start_routine = c_null_ptr
  end function start_routine
end module forthread_data
