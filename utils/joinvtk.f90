module Kinds
  use iso_fortran_env
  
  integer,parameter :: rp = real32
end module

module Endianness
  use iso_fortran_env
  
  implicit none

  private
  public :: GetEndianness, BigEnd

  logical,save :: littleendian = .true.

  interface BigEnd
    module procedure BigEnd32
    module procedure BigEnd64
  end interface

contains

  subroutine GetEndianness
    integer(int8),dimension(4):: bytes !may not work on some processors

    bytes=transfer(1_int32,bytes,4)
    if (bytes(4)==1) then
      littleendian=.false.
    else
      littleendian=.true.
    endif
  end subroutine GetEndianness

  elemental function BigEnd32(x) result(res)
    real(real32) :: res
    real(real32),intent(in)::x
    integer(int8),dimension(4):: bytes !may not work on some processors

    if (.not.littleendian) then
      res = x
    else
      bytes = transfer(x,bytes,4)
      res = transfer(bytes(4:1:-1),res)
    endif
  end function BigEnd32

  elemental function BigEnd64(x) result(res)
    real(real64) :: res
    real(real64),intent(in)::x
    integer(int8),dimension(8):: bytes !may not work on some processors

    if (.not.littleendian) then
      res = x
    else
      bytes = transfer(x,bytes,8)
      res = transfer(bytes(8:1:-1),res)
    endif
  end function BigEnd64

end module Endianness




module Types
  use Kinds
  
  implicit none

  type grid
    integer :: nx,ny,nz
    integer :: offx,offy,offz
    integer :: unit
    real(rp),allocatable :: x(:),y(:),z(:)
    character(100) :: fname
  contains
    procedure :: read_header
    procedure :: read_title
    procedure :: read_scalar
    procedure :: read_vector
  end type
  
  integer, parameter :: SCALAR = 1, VECTOR = 2
  character,parameter :: lf = achar(10)

contains
  

  subroutine read_header(g)
    use Endianness
    class(grid),intent(inout) :: g
    integer :: io
    character(70) :: str
    character :: ch

    open(newunit=g%unit,file=g%fname,access="stream",form="unformatted",status="old",action="read",iostat=io)

    if (io/=0) then
      g%nx = 0
      g%ny = 0
      g%nz = 0
      return
    end if

    call skip_to("X_COORDINATES", io)
    if (io/=0) then
      write(*,*) "Error, cannot find X_COORDINATES in '",g%fname,"'"
      stop
    end if
    call get_number(g%nx)
    call skip_line
    
    allocate(g%x(g%nx))
    read(g%unit,iostat=io) g%x
    g%x=BigEnd(g%x)

    call skip_to("Y_COORDINATES", io)
    if (io/=0) then
      write(*,*) "Error, cannot find Y_COORDINATES in '",g%fname,"'"
      stop
    end if
    call get_number(g%ny)
    call skip_line
    
    allocate(g%y(g%ny))
    read(g%unit,iostat=io) g%y
    g%y=BigEnd(g%y)

    call skip_to("Z_COORDINATES", io)
    if (io/=0) then
      write(*,*) "Error, cannot find Z_COORDINATES in '",g%fname,"'"
      stop
    end if
    call get_number(g%nz)
    call skip_line
    
    allocate(g%z(g%nz))
    read(g%unit,iostat=io) g%z
    g%z=BigEnd(g%z)
    
    call skip_line
    call skip_line
    
  contains
    subroutine skip_line
      character :: ch
      do
        read(g%unit) ch
        if (ch==new_line("a")) return
      end do
    end subroutine

    subroutine get_number(n)
      integer, intent(out) :: n
      character :: ch
      character(:), allocatable :: num_str
      
      do
        read(g%unit) ch
        if (ch/=' ') exit
      end do
      
      num_str=ch
      do
        read(g%unit) ch
        if (ch<'0' .or. ch>'9') exit
        num_str = num_str // ch
      end do

      read(num_str,*) n
    end subroutine

    subroutine skip_to(str, stat)
      character(*), intent(in) :: str
      integer, intent(out) :: stat
      character :: ch
      integer :: io

      do
        read(g%unit, iostat=io) ch

        if (io/=0) then
          stat = 1
          return
        end if

        if (ch==str(1:1)) then
          call check(str(2:), stat)
          if (stat == 0) return
        end if

      end do
    end subroutine
    
    subroutine check(str, stat)
      character(*), intent(in) :: str
      integer, intent(out) :: stat
      character :: ch
      integer :: i, io

      stat = 1
      i = 0

      do
        i = i + 1

        read(g%unit, iostat=io) ch

        if (io/=0) return

        if (ch/=str(i:i)) return

        if (i==len(str)) then
          stat = 0
          return
        end if
      end do
    end subroutine
   
  end subroutine
  
  subroutine read_title(g,status,title)
    use iso_fortran_env
    class(grid),intent(in) :: g
    integer :: status
    character(:),allocatable,intent(out) :: title
    character(7) :: vtype
    character :: ch
    integer :: io, n, lf_pos
    character(256) :: msg
  
    do
      vtype=""
      read(g%unit,iostat=io,iomsg=msg) vtype
      
      if (io==iostat_end) then
        status = 0
        return
      else if (io/=0) then
        write(*,'(*(g0))') io, trim(msg), "  ", g%fname, " ", "'",vtype,"'"
        stop "Error reading title."
      end if
      
      if (vtype=='VECTORS') then
      
        status = VECTOR
        title = vtype
        do
          read(g%unit) ch
          title = title // ch
          if (ch==lf) exit
        end do
        
        exit
        
      else if (vtype=='SCALARS') then
      
        status = SCALAR
        title = vtype
        n = 0
        do
          read(g%unit) ch
          title = title // ch
          if (ch==lf) n = n + 1
          if (n==2) exit
        end do
        
        exit
        
      else
        !if our vtype buffer contains the end-of-line LF character, move back just behind it
        !otherwise skip to the nearest one
        lf_pos = scan(vtype, achar(10),.true.)
        if (lf_pos>0) then
          call stream_back(len(vtype)-lf_pos+1)
        else
          call skip_line
        end if
        
      end if
    end do
    
  contains
  
    subroutine skip_line
      character :: ch
      do
        read(g%unit) ch
        if (ch==new_line("a")) return
      end do
    end subroutine
    
    subroutine stream_back(n)
      integer, intent(in) :: n
      integer :: pos
      character :: dummy
      inquire(unit=g%unit,pos=pos)
      read(g%unit,pos=pos-n) dummy
    end subroutine

  end subroutine
  
  subroutine read_scalar(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:)
    character :: ch
 
    read(g%unit) buf(g%offx:g%offx+g%nx-1,g%offy:g%offy+g%ny-1,g%offz:g%offz+g%nz-1)
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:,:)
    character :: ch
    
    read(g%unit) buf(:,g%offx:g%offx+g%nx-1,g%offy:g%offy+g%ny-1,g%offz:g%offz+g%nz-1)
    read(g%unit) ch
  end subroutine
  
end module Types





program joinvtk
  use Kinds
  use Types
  use Endianness
  
  implicit none
  
  type(grid), allocatable :: ims(:,:,:)
  integer :: unit
  integer :: nxims, nyims, nzims, nims
  integer :: offx, offy, offz
  integer :: nx, ny, nz
  real(rp), allocatable :: sc(:,:,:),vec(:,:,:,:)
 
  character(:), allocatable :: title
  integer :: status, io
  
  character(80) :: file_name, arg
  
  integer :: i, j, k
  
  nxims = 1
  
  nyims = 2
  
  nzims = 2
  
  call GetEndianness

  if (command_argument_count()<1) then
    write(*,*) "Usage: joinvtk file_name npx,npy,npz"
    stop 1
  end if
 
  call get_command_argument(1, value=file_name)

  call get_command_argument(2, value=arg, status=io)
  if (io==0) read(arg,*,iostat=io) nxims, nyims, nzims
  
  if (io/=0) then
    call get_shape_from_directories()
  end if
  
  nims = nxims * nyims * nzims
 
  allocate(ims(nxims,nyims,nzims))
 
  do k = 1, nzims
    do j = 1, nyims
      do i = 1, nxims
        write(ims(i,j,k)%fname,'(*(g0))') "im-",i,"-",j,"-",k,"/",file_name
        call ims(i,j,k)%read_header
      end do
    end do
  end do
  
  if (all(ims%nx==0)) then
    write(*,*) "No input file found."
    stop 3
  end if
    
  
  if (.not.all(ims%nx>0)) then
    call crop_images
  end if
  
  offx = 1
  do i = 1, nxims
    ims(i,:,:)%offx = offx
    offx = offx + ims(i,1,1)%nx
  end do
  offy = 1
  do j = 1, nyims
    ims(:,j,:)%offy = offy
    offy = offy + ims(1,j,1)%ny
  end do
  offz = 1
  do k = 1, nzims
    ims(:,:,k)%offz = offz
    offz = offz + ims(1,1,k)%nz
  end do
  
  nx = offx-1
  ny = offy-1
  nz = offz-1
  
  allocate(sc(1:nx,1:ny,1:nz))
  allocate(vec(1:3,1:nx,1:ny,1:nz))
  
  call save_header
  
  do
    call get_next(status,title)
    if (status==0) exit
    call get_buffer(status,sc,vec)
    call save_buffer(status,title,sc,vec)
  end do
 
  do k = 1, nzims
    do j = 1, nyims
      do i = 1, nxims
        close(ims(i,j,k)%unit)
      end do
    end do
  end do
  close(unit)

contains

  subroutine crop_images
    integer :: i, j, k
    integer :: mini, minj, mink
    integer :: maxi, maxj, maxk
    
    mini = 1
    do i = 1, nxims
      if (all(ims(i,:,:)%nx==0)) then
        mini = i + 1
      else
        exit
      end if
    end do
    maxi = nxims
    do i = nxims, 1, -1
      if (all(ims(i,:,:)%nx==0)) then
        maxi = i - 1
      else
        exit
      end if
    end do
    
    minj = 1
    do j = 1, nyims
      if (all(ims(:,j,:)%nx==0)) then
        minj = j + 1
      else
        exit
      end if
    end do
    maxj = nyims
    do j = nyims, 1, -1
      if (all(ims(:,j,:)%nx==0)) then
        maxj = j - 1
      else
        exit
      end if
    end do
    
    mink = 1
    do k = 1, nzims
      if (all(ims(:,:,k)%nx==0)) then
        mink = k + 1
      else
        exit
      end if
    end do
    maxk = nzims
    do k = nzims, 1, -1
      if (all(ims(:,:,k)%nx==0)) then
        maxk = k - 1
      else
        exit
      end if
    end do
    
    ims = ims(mini:maxi,minj:maxj,mink:maxk)

    nxims = size(ims, 1)
    nyims = size(ims, 2)
    nzims = size(ims, 3)
  end subroutine

  subroutine get_next(status,title)
    integer :: status, all_status
    character(:),allocatable,intent(out) :: title
    integer :: i,j,k
    
    all_status = -1
    do k = 1, nzims
      do j = 1, nyims
        do i = 1, nxims

          call ims(i,j,k)%read_title(status,title)

          if (all_status==-1) then
            all_status = status
          else if (status /= all_status) then
            stop "Error status /= all_status."
          end if
        end do
      end do
    end do
  end subroutine
  
  subroutine get_buffer(status,sc,vec)
    integer,intent(in) :: status
    real(rp),intent(out) :: sc(:,:,:),vec(:,:,:,:)
    integer :: i,j,k
    
    do k = 1, nzims
      do j = 1, nyims
        do i = 1, nxims
          if (status==SCALAR) then
            call ims(i,j,k)%read_scalar(sc)
          else
            call ims(i,j,k)%read_vector(vec)
          end if
        end do
      end do
    end do
  end subroutine
  
  subroutine save_header
    integer :: i
    character(70) :: str
    real(rp), allocatable :: x(:),y(:),z(:)
  
    allocate(x(0))
    do i=1,nxims
      if (.not.allocated(ims(i,1,1)%x)) then
        write(*,*) "Error at image",i,1,1," - x not allocated."
        stop
      end if
      x = [x,ims(i,1,1)%x]
    end do
    allocate(y(0))
    do i=1,nyims
      if (.not.allocated(ims(1,i,1)%y)) then
        write(*,*) "Error at image",1,i,1," - y not allocated."
        stop
      end if
      y = [y,ims(1,i,1)%y]
    end do
    allocate(z(0))
    do i=1,nzims
      if (.not.allocated(ims(1,1,i)%z)) then
        write(*,*) "Error at image",1,1,i," - z not allocated."
        stop
      end if
      z = [z,ims(1,1,i)%z]
    end do
    

    open(newunit=unit,file=file_name, &
      access='stream',status='replace',form="unformatted",action="write")

    write(unit) "# vtk DataFile Version 2.0",lf
    write(unit) "CLMM output file",lf
    write(unit) "BINARY",lf
    write(unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write(str(12:),*) nx,ny,nz
    write(unit) str,lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(x, real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(y, real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(z, real32)),lf
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit) str,lf 
  end subroutine
  
  subroutine save_buffer(status,title,sc,vec)
    integer,intent(in) :: status
    character(*),intent(in) :: title
    real(rp),intent(in) :: sc(:,:,:),vec(:,:,:,:)
    
    write(unit) title
    if (status==SCALAR) then
      write(unit) sc,lf
    else
      write(unit) vec,lf
    end if
  end subroutine
  
  subroutine get_shape_from_directories()
    use iso_c_binding
    character(30) :: dirname
    integer i
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
  
  end subroutine
  
end program
