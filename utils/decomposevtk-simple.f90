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
      write(*,*) "Error reading from ",trim(g%fname)
      stop 2
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
  use Parameters
  use Strings
  
  implicit none
  
  type(grid) :: glob

  real(rp), allocatable :: sc(:,:,:),vec(:,:,:,:)
 
  character(:), allocatable :: title, base_name
  integer :: status, io
  
  character(80) :: file_name = "", arg

  integer, allocatable, dimension(:,:,:) :: nx, ny, nz, unit, offx, offy, offz

  integer :: nxims, nyims, nzims

  integer :: iim, jim, kim

   
  call GetEndianness

  if (command_argument_count()<2) then
    write(*,*) "Usage: joinvtk file_name npx,npy,npz nx,ny,nz"
    stop 1
  end if
 
  call get_command_argument(1, value=file_name)

  call get_command_argument(2, value=arg)
  read(arg,*,iostat=io) nxims, nyims, nzims
  
  allocate(nx(nxims, nyims, nzims))
  allocate(ny(nxims, nyims, nzims))
  allocate(nz(nxims, nyims, nzims))
  allocate(offx(nxims, nyims, nzims))
  allocate(offy(nxims, nyims, nzims))
  allocate(offz(nxims, nyims, nzims))
  allocate(unit(nxims, nyims, nzims))

  if (io/=0) then
    write(*,*) "Error reading image grid shape, provide three integers 'npx,npy,npz'."
    stop 2
  end if
  
  !The grid shape is necessary, if we do not want to parse the boundary conditions
  ! and still want to be able to decompose fields in U, V and W staggered grids.
  call get_command_argument(3, value=arg)
  read(arg,*,iostat=io) Prnx, Prny, Prnz
  
  if (io/=0) then
    write(*,*) "Error reading grid shape, provide three integers 'Prnx,Prny,Prnz'."
    stop 3
  end if

  !Simplified decomposition. Usable only in evenly divisible grids and some other less predictable cases.
  gPrnx = Prnx
  gPrny = Prny
  gPrnz = Prnz

  glob%fname = file_name
  base_name = trim(file_name)

  call glob%read_header

  if (glob%nx==0 .or. glob%ny==0 .or. glob%nz==0) then
    write(*,*) "File ",trim(glob%fname)," does not exist or cannot be read."
    stop 4
  end if

  if (glob%nx>gPrnx.or.glob%nx<gPrnx-1) then
    write(*,*) "Inconsistent grid size, supplied Prnx:",gPrnx,"nx in ",trim(glob%fname),':',glob%nx
    stop 4
  end if
  if (glob%ny>gPrny.or.glob%ny<gPrny-1) then
    write(*,*) "Inconsistent grid size, supplied Prny:",gPrny,"ny in ",trim(glob%fname),':',glob%ny
    stop 4
  end if
  if (glob%nz>gPrnz.or.glob%nz<gPrnz-1) then
    write(*,*) "Inconsistent grid size, supplied Prnz:",gPrnz,"nz in ",trim(glob%fname),':',glob%nz
    stop 4
  end if

  allocate(sc(1:glob%nx,1:glob%ny,1:glob%nz))
  allocate(vec(1:3,1:glob%nx,1:glob%ny,1:glob%nz))
        

  do kim = 1, nzims
    do jim = 1, nyims
      do iim = 1, nxims

        Prnx = ceiling(1.d0 * gPrnx / nxims)
        Prny = ceiling(1.d0 * gPrny / nyims)
        Prnz = ceiling(1.d0 * gPrnz / nzims)

        offx(iim,jim,kim) = (iim-1) * Prnx
        offy(iim,jim,kim) = (jim-1) * Prny
        offz(iim,jim,kim) = (kim-1) * Prnz

        Prnx = min(Prnx, gPrnx - offx(iim,jim,kim))
        Prny = min(Prny, gPrny - offy(iim,jim,kim))
        Prnz = min(Prnz, gPrnz - offz(iim,jim,kim))

        if (Prnx + offx(iim,jim,kim) == glob%nx+1) then
          nx(iim,jim,kim) = Prnx - 1
        else
          nx(iim,jim,kim) = Prnx
        end if
        
        if (Prny + offy(iim,jim,kim) == glob%ny+1) then
          ny(iim,jim,kim) = Prny - 1
        else
          ny(iim,jim,kim) = Prny
        end if

        if (Prnz + offz(iim,jim,kim) == glob%nz+1) then
          nz(iim,jim,kim) = Prnz - 1
        else
          nz(iim,jim,kim) = Prnz
        end if
  
        call execute_command_line("mkdir -p "//"im-" // itoa(iim) // "-" // itoa(jim) // "-" // itoa(kim))
        file_name = "im-" // itoa(iim) // "-" // itoa(jim) // "-" // itoa(kim) // "/" // base_name

        call save_header
      end do
    end do
  end do

        
  do
    call get_next(status,title)
    if (status==0) exit

    call get_buffer(status,sc,vec)

    do kim = 1, nzims
      do jim = 1, nyims
        do iim = 1, nxims
          call save_buffer(status,title,sc,vec)
        end do
      end do
    end do

  end do


  do kim = 1, nzims
    do jim = 1, nyims
      do iim = 1, nxims
        close(unit(iim,jim,kim))
      end do
    end do
  end do

  close(glob%unit)

  deallocate(sc)
  deallocate(vec)

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
  
    x = glob%x(offx(iim,jim,kim)+1:offx(iim,jim,kim)+nx(iim,jim,kim))
    y = glob%y(offy(iim,jim,kim)+1:offy(iim,jim,kim)+ny(iim,jim,kim))

    z = glob%z(offz(iim,jim,kim)+1:offz(iim,jim,kim)+nz(iim,jim,kim))
    

    open(newunit=unit(iim,jim,kim),file=file_name, &
      access='stream',status='replace',form="unformatted",action="write")

    write(unit(iim,jim,kim)) "# vtk DataFile Version 2.0",lf
    write(unit(iim,jim,kim)) "CLMM output file",lf
    write(unit(iim,jim,kim)) "BINARY",lf
    write(unit(iim,jim,kim)) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write(str(12:),*) nx(iim,jim,kim),ny(iim,jim,kim),nz(iim,jim,kim)
    write(unit(iim,jim,kim)) str,lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx(iim,jim,kim),"float"
    write(unit(iim,jim,kim)) str,lf
    write(unit(iim,jim,kim)) BigEnd(real(x, real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny(iim,jim,kim),"float"
    write(unit(iim,jim,kim)) str,lf
    write(unit(iim,jim,kim)) BigEnd(real(y, real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz(iim,jim,kim),"float"
    write(unit(iim,jim,kim)) str,lf
    write(unit(iim,jim,kim)) BigEnd(real(z, real32)),lf
    str="POINT_DATA"
    write(str(12:),*) nx(iim,jim,kim)*ny(iim,jim,kim)*nz(iim,jim,kim)
    write(unit(iim,jim,kim)) str,lf 
  end subroutine
  
  subroutine save_buffer(status,title,sc,vec)
    integer,intent(in) :: status
    character(*),intent(in) :: title
    real(rp),intent(in) :: sc(:,:,:),vec(:,:,:,:)
    
    write(unit(iim,jim,kim)) title
    if (status==SCALAR) then
      write(unit(iim,jim,kim)) sc(offx(iim,jim,kim)+1:offx(iim,jim,kim)+nx(iim,jim,kim), &
                     offy(iim,jim,kim)+1:offy(iim,jim,kim)+ny(iim,jim,kim), &
                     offz(iim,jim,kim)+1:offz(iim,jim,kim)+nz(iim,jim,kim)), &
                  lf
    else
      write(unit(iim,jim,kim)) vec(:, &
                      offx(iim,jim,kim)+1:offx(iim,jim,kim)+nx(iim,jim,kim), &
                      offy(iim,jim,kim)+1:offy(iim,jim,kim)+ny(iim,jim,kim), &
                      offz(iim,jim,kim)+1:offz(iim,jim,kim)+nz(iim,jim,kim)), &
                  lf
    end if
  end subroutine


  
end program
