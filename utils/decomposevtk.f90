module WorkKinds
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
  use WorkKinds
  
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
    integer :: io,nx,ny,nz
    character(5) :: ch5
    character(70) :: str
    character :: ch

    open(newunit=g%unit,file=g%fname,access="stream",form="unformatted",status="old",action="read",iostat=io)

    if (io/=0) then
      g%nx = 0
      g%ny = 0
      g%nz = 0
      return
    end if

    read(g%unit,pos=162,iostat=io) ch5
    if (io/=0) call error

    read(ch5,'(i5)') g%nx
    nx=g%nx
    allocate(g%x(g%nx))
    read(g%unit,pos=219,iostat=io) g%x
    g%x=BigEnd(g%x)

    read(g%unit,pos=234+nx*4,iostat=io) ch5
    read(ch5,'(i5)') g%ny
    ny=g%ny

    allocate(g%y(g%ny))
    read(g%unit,pos=291+nx*4,iostat=io) g%y
    g%y=BigEnd(g%y)

    read(g%unit,pos=306+nx*4+ny*4,iostat=io) ch5
    read(ch5,'(i5)') g%nz
    nz=g%nz
    allocate(g%z(g%nz))
    read(g%unit,pos=363+nx*4+ny*4,iostat=io) g%z
    g%z=BigEnd(g%z)
    read(g%unit) ch,str
    read(g%unit) ch
    return
    
  contains
    subroutine error()
      use custom_par
      if (master) write(*,*) "Error reading from ",trim(g%fname)
      call error_stop
    end subroutine    
  end subroutine
  
  subroutine read_title(g,status,title)
    use iso_fortran_env
    class(grid),intent(in) :: g
    integer :: status
    character(:),allocatable,intent(out) :: title
    character(7) :: vtype
    character :: ch
    integer :: io, n
    character(256) :: msg
    
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
    else
      write(*,'(*(g0))') "'",vtype,"'"
      stop "Error reading title."
    end if
  end subroutine
  
  subroutine read_scalar(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:)
    character :: ch
 
    read(g%unit) buf(g%offx+1:g%offx+g%nx, g%offy+1:g%offy+g%ny, g%offz+1:g%offz+g%nz)
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:,:)
    character :: ch
    
    read(g%unit) buf(:, g%offx+1:g%offx+g%nx, g%offy+1:g%offy+g%ny, g%offz+1:g%offz+g%nz)
    read(g%unit) ch
  end subroutine
  
end module Types





program decomposevtk
  use WorkKinds
  use Types
  use Endianness
  use custom_par
  use Parameters
  use Strings
  
  implicit none
  
  type(grid) :: glob
  integer :: unit
  real(rp), allocatable :: sc(:,:,:),vec(:,:,:,:)
 
  character(:), allocatable :: title
  integer :: status, io
  
  character(80) :: file_name, arg
  
  integer :: im

  integer :: nx, ny, nz

  call par_init

  
  call GetEndianness

  if (command_argument_count()<2) then
    call error_stop("Usage: joinvtk file_name npx,npy,npz")
  end if
 
  call get_command_argument(1, value=file_name)

  call get_command_argument(2, value=arg)
  read(arg,*,iostat=io) npxyz
  
  if (io/=0) then
    call error_stop("Error reading image grid shape, provide three integers 'npx,npy,npz'.")
  end if
  
  !The grid shape is necessary, if we do not want to parse the boundary conditions
  ! and still want to be able to decompose fields in U, V and W staggered grids.
  call get_command_argument(3, value=arg)
  read(arg,*,iostat=io) Prnx, Prny, Prnz
  
  if (io/=0) then
    call error_stop("Error reading grid shape, provide three integers 'Prnx,Prny,Prnz'.")
  end if

  call par_init_grid


  glob%fname = file_name

  call glob%read_header
  
  nx = glob%nx
  ny = glob%ny
  nz = glob%nz

  if (nx==0 .or. ny==0 .or. nz==0) then
    write(*,*) "File ",trim(glob%fname)," does not exist or cannot be read."
    call error_stop
  end if

  if (nx>gPrnx.or.nx<gPrnx-1) then
    write(*,*) "Inconsistent grid size, supplied Prnx:",gPrnx,"nx in ",trim(glob%fname),':',nx
    call error_stop
  end if
  if (ny>gPrny.or.ny<gPrny-1) then
    write(*,*) "Inconsistent grid size, supplied Prny:",gPrny,"ny in ",trim(glob%fname),':',ny
    call error_stop
  end if
  if (nz>gPrnz.or.nz<gPrnz-1) then
    write(*,*) "Inconsistent grid size, supplied Prnz:",gPrnz,"nz in ",trim(glob%fname),':',nz
    call error_stop
  end if

  if (Prnx + offset_to_global_x == glob%nx+1) then
    nx = Prnx - 1
  else
    nx = Prnx
  end if
  
  if (Prny + offset_to_global_y == glob%ny+1) then
    ny = Prny - 1
  else
    ny = Prny
  end if

  if (Prnz + offset_to_global_z == glob%nz+1) then
    nz = Prnz - 1
  else
    nz = Prnz
  end if
  
  !Serialized to avoid memory overflow. Hopefully, the deallocated memory is returned to the OS (mmap).
  do im = 1, nims

    call par_sync_all

    if (im==myim) then

      call execute_command_line("mkdir -p "//"im-" // itoa(iim) // "-" // itoa(jim) // "-" // itoa(kim))
      file_name = "im-" // itoa(iim) // "-" // itoa(jim) // "-" // itoa(kim) // "/" // file_name

      allocate(sc(1:glob%nx,1:glob%ny,1:glob%nz))
      allocate(vec(1:3,1:glob%nx,1:glob%ny,1:glob%nz))
      
      call save_header
      
      do
        call get_next(status,title)
        if (status==0) exit
        call get_buffer(status,sc,vec)
        call save_buffer(status,title,sc,vec)
      end do

      close(glob%unit)

      close(unit)

      deallocate(sc)
      deallocate(vec)

    end if

    call par_sync_all

  end do

  call par_finalize

contains

 
  subroutine get_next(status,title)
    integer :: status
    character(:), allocatable, intent(out) :: title

    call glob%read_title(status,title)

  end subroutine
  
  subroutine get_buffer(status, sc, vec)
    integer, intent(in) :: status
    real(rp), intent(out) :: sc(:,:,:), vec(:,:,:,:)

    if (status==SCALAR) then
      call glob%read_scalar(sc)
    else
      call glob%read_vector(vec)
    end if

  end subroutine
  
  subroutine save_header
    character(70) :: str
    real(rp), allocatable :: x(:),y(:),z(:)
  
    x = glob%x(offset_to_global_x+1:offset_to_global_x+nx)
    y = glob%y(offset_to_global_y+1:offset_to_global_y+ny)

    z = glob%z(offset_to_global_z+1:offset_to_global_z+nz)
    

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
      write(unit) sc(offset_to_global_x+1:offset_to_global_x+nx, &
                     offset_to_global_y+1:offset_to_global_y+ny, &
                     offset_to_global_z+1:offset_to_global_z+nz), &
                  lf
    else
      write(unit) vec(:, &
                      offset_to_global_x+1:offset_to_global_x+nx, &
                      offset_to_global_y+1:offset_to_global_y+ny, &
                      offset_to_global_z+1:offset_to_global_z+nz), &
                  lf
    end if
  end subroutine


  
end program
