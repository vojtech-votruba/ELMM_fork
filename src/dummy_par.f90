module custom_par
  !Dummy version  of global synchronization routines for non-distributed version
  ! to avoid excessive need of preprocessing macros around common calls.

  integer, parameter :: nims = 1, npx = 1, npy = 1, npz = 1
  integer, parameter :: npxyz(3) = 1
  integer, parameter :: nxims = 1, nyims = 1, nzims = 1
  integer, parameter :: myim = 1, iim = 1, jim = 1, kim = 1  

contains

  subroutine par_init()
    
  end subroutine


  subroutine par_finalize()
    
  end subroutine


  subroutine par_sync_all()
    
  end subroutine

  subroutine par_sync_out(str)
    character(*), intent(in) :: str
    write(*,*) str
  end subroutine

end module
