program joinvtkframes

  character(:), allocatable :: domain, filename, cmd
  character(80) :: arg
  integer :: i, stat
  
  if (command_argument_count()<1) then
    write(*,*) "Usage: joinvtkframes label npx,npy,npz"
    stop 1
  end if
 
  call get_command_argument(1, value=arg)
  
  domain = trim(arg)

  if (command_argument_count()==2) then
    call get_command_argument(2, value=arg)
  else
    call get_shape_from_directories(arg)
  end if

  i = 0
  do
    filename = 'frame-'//domain//'-'//itoa(i)//'.vtk'
    cmd = 'joinvtk '//filename//' '//arg
    
    write(*,*) filename

    call execute_command_line(cmd, exitstat=stat)

    if (stat/=0) exit
    i = i + 1
  end do
    
contains
  
  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function
  
  subroutine get_shape_from_directories(str)
    use iso_c_binding
    character(30) :: dirname
    character(*) :: str
    integer :: i, nxims, nyims, nzims
    INTERFACE
      SUBROUTINE file_info(filename,mode,exist,time) BIND(C,name="file_info")
        USE iso_c_binding
        CHARACTER(kind=C_CHAR),INTENT(in) :: filename(*)
        INTEGER(C_INT),INTENT(out) :: mode,exist,time
      END SUBROUTINE
    END INTERFACE
    integer(c_int) :: mode, ex, time

    
    i = 1
    do
      write(dirname,'(*(g0))') "im-",i,"-",1,"-",1,"/",achar(0)
      call file_info(dirname,mode,ex,time)      
      if (ex==0) exit
      i = i + 1
    end do
    nxims = i-1
    
    i = 1
    do
      write(dirname,'(*(g0))') "im-",1,"-",i,"-",1,"/",achar(0)
      call file_info(dirname,mode,ex,time)
      if (ex==0) exit
      i = i + 1
    end do
    nyims = i-1
    
    i = 1
    do
      write(dirname,'(*(g0))') "im-",1,"-",1,"-",i,"/",achar(0)
      call file_info(dirname,mode,ex,time)
      if (ex==0) exit
      i = i + 1
    end do
    nzims = i-1
  
    str = itoa(nxims)//itoa(nyims)//itoa(nzims)
  end subroutine
  
end program
