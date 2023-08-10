module Stop_procedures
  use iso_fortran_env, only: error_unit
#ifdef MPI
  use mpi
#endif
  private
  
  public master, error_stop

#ifdef PAR
  logical :: master = .true.
#else
  logical, parameter :: master = .true.
#endif

  interface error_stop
    module procedure error_stop_null
    module procedure error_stop_int
    module procedure error_stop_char
    module procedure error_stop_char_int
  end interface
  
contains
    
  subroutine error_stop_null
    
#ifdef MPI
    call abort_mpi(1)
#endif
    stop
  end subroutine
  

  subroutine error_stop_int(n)
    integer,intent(in) :: n

    write(error_unit,'(*(g0))') "ELMM Abort: ", n
    
#ifdef MPI
    call abort_mpi(n)
#endif
    stop
  end subroutine
  

  subroutine error_stop_char(ch)
    character(*),intent(in) :: ch

    write(error_unit,'(*(g0))') "ELMM Abort: ", ch

#ifdef MPI
    call abort_mpi(1)
#endif
    stop
  end subroutine

  
  subroutine error_stop_char_int(ch,n)
    character(*),intent(in) :: ch
    integer,intent(in) :: n
    
    write(error_unit,'(2g0,1x,g0)') "ELMM Abort: ", ch, n

#ifdef MPI
    call abort_mpi(n)
#endif
    stop
  end subroutine

#ifdef MPI
  subroutine abort_mpi(n)
    integer, intent(in) :: n
    integer :: ie

    call MPI_abort(MPI_COMM_WORLD,n,ie)

  end subroutine
#endif

end module

