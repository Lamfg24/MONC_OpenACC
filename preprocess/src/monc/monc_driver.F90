










!> MONC program entry point which simply calls the main procedure in the MONC module
!  
program monc_driver
  use monc_mod, only : monc_core_bootstrap

  implicit none

  call monc_core_bootstrap()

end program monc_driver
