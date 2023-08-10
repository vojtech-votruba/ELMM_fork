module kinds
  use iso_fortran_env
  
  integer, parameter :: rp = real32
end module

program series2txt
  use kinds

  implicit none

  character(1024) :: arg, fname, path, valmsg, timmsg, resmsg
  
  integer :: uval, utim, ures, ioval, iotim, iores, n
  
  real(rp) :: tim
  real(rp), allocatable :: val(:)

  if (command_argument_count()<1) then
    write(*,*) "Supply the basename of the file on the command line."
  end if
    
  call get_command_argument(1, fname)

  if (command_argument_count()>1) then
    call get_command_argument(2, arg)
    read(arg,*) n
  else
    n = 1
  end if
    
    
  allocate(val(n))
  
  path = fname( : scan(fname,'/',back=.true.)-1 )
  
  open(newunit = uval, &
       file = trim(fname)//'.unf', &
       action = 'read', &
       access = 'stream', &
       form = 'unformatted', &
       status = 'old', &
       iostat = ioval, &
       iomsg = valmsg)
  open(newunit = utim, &
       file = trim(path)//'times.unf', &
       action = 'read', &
       access = 'stream', &
       form = 'unformatted', &
       status = 'old', &
       iostat = iotim, &
       iomsg = timmsg)
  open(newunit = ures, &
       file = trim(fname)//'.txt', &
       action = 'write', &
       status = 'replace', &
       iostat = iores, &
       iomsg = resmsg)
       
  if (ioval/=0) then
    write(*,*) "Error opening file '",trim(fname)//'.unf',"'"
    write(*,*) trim(valmsg)
    stop
  end if
  
  if (iotim/=0) then
    write(*,*) "Error opening file '",trim(path)//'times.unf',"'"
    write(*,*) trim(timmsg)
    stop
  end if
  
  if (iores/=0) then
    write(*,*) "Error opening file '",trim(fname)//'.txt',"'"
    write(*,*) trim(resmsg)
    stop
  end if
  
  do
    read(utim, iostat=iotim) tim
    read(uval, iostat=ioval,iomsg=valmsg) val   
    if (iotim/=0 .or. ioval/=0) exit
    write(ures, *) tim, val
  end do
  
end program

