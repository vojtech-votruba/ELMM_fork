module kinds
  use, intrinsic :: iso_fortran_env
  integer, parameter :: rp = real32
end module

module vtk_kinds
  use iso_fortran_env
  
  integer,parameter :: vtk_rp = real32
end module

module Endianness
  use iso_fortran_env
  
  implicit none

  private
  public :: GetEndianness, BigEnd

  logical,save :: littleendian

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
  use Endianness
  use vtk_kinds
  
  implicit none

  type grid
    integer :: nx,ny,nz
    integer :: unit
    real(vtk_rp),allocatable :: x(:),y(:),z(:)
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
  

  subroutine read_header(g, fname)
    class(grid),intent(inout) :: g
    character(*), intent(in) :: fname
    integer :: io,nx,ny,nz
    character(5) :: ch5
    character(70) :: str
    character :: ch

    g%fname = fname
    open(newunit=g%unit,file=g%fname,access="stream",form="unformatted",status="old",action="read",iostat=io)

    if (io/=0) call error

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
    
    vtype = ""
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
    real(vtk_rp),intent(out) :: buf(:,:,:)
    character :: ch
 
    read(g%unit) buf
    buf = BigEnd(buf)
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
    class(grid),intent(in) :: g
    real(vtk_rp),intent(out) :: buf(:,:,:,:)
    character :: ch
    
    read(g%unit) buf
    buf = BigEnd(buf)
    read(g%unit) ch
  end subroutine
  
end module Types


module Periodic
  use Kinds
  use Types
  
  implicit none
  

contains


  subroutine compute_and_save(input_file, px, py, pz)
    type(grid) :: g, mean
    character(*), intent(in) :: input_file
    integer, intent(in) :: px, py, pz
    real(vtk_rp), allocatable :: sc(:,:,:), vec(:,:,:,:)
    real(vtk_rp), allocatable :: mean_sc(:,:,:), mean_vec(:,:,:,:)
    integer :: status
    character(:), allocatable :: title, dir
    integer :: unit
    
    call g%read_header(input_file)

    mean%nx = g%nx / px
    mean%ny = g%ny / py
    mean%nz = g%nz / pz

    mean%x = g%x(1:mean%nx)
    mean%y = g%y(1:mean%ny)
    mean%z = g%z(1:mean%nz)

    allocate(character(0)::dir)
    dir = dirname(input_file)
    if (len(dir)>0) dir = dir // '/'
    mean%fname = dir // "mean_" // basename(input_file)

    allocate(sc(g%nx,g%ny,g%nz))
    allocate(vec(3,g%nx,g%ny,g%nz))

    allocate(mean_sc(mean%nx,mean%ny,mean%nz))
    allocate(mean_vec(3,mean%nx,mean%ny,mean%nz))

    call save_header

    do
      call get_next(status,title)
      if (status==0) exit
      call get_buffer(status,sc,vec)

      select case (status)
        case (SCALAR)
          call mean_scalar(sc, mean_sc)
        case (VECTOR)
          call mean_vector(vec, mean_vec)
      end select
      call save_buffer(status, title, mean_sc, mean_vec)
    end do

    close(unit)

  contains

    subroutine get_next(status,title)
      integer :: status, all_status
      character(:),allocatable,intent(out) :: title
      
      all_status = -1


      call g%read_title(status,title)

      if (all_status==-1) then
        all_status = status
      else if (status /= all_status) then
        stop "Error status /= all_status."
      end if

    end subroutine
    
    subroutine get_buffer(status,sc,vec)
      integer,intent(in) :: status
      real(vtk_rp),intent(out) :: sc(:,:,:),vec(:,:,:,:)

      if (status==SCALAR) then
        call g%read_scalar(sc)
      else
        call g%read_vector(vec)
      end if

    end subroutine

    subroutine save_header
      character(70) :: str

      open(newunit=unit,file=mean%fname, &
        access='stream',status='replace',form="unformatted",action="write")

      write(unit) "# vtk DataFile Version 2.0",lf
      write(unit) "CLMM output file",lf
      write(unit) "BINARY",lf
      write(unit) "DATASET RECTILINEAR_GRID",lf
      str="DIMENSIONS"
      write(str(12:),*) mean%nx,mean%ny,mean%nz
      write(unit) str,lf
      str="X_COORDINATES"
      write(str(15:),'(i5,2x,a)') mean%nx,"float"
      write(unit) str,lf
      write(unit) BigEnd(real(mean%x, real32)),lf
      str="Y_COORDINATES"
      write(str(15:),'(i5,2x,a)') mean%ny,"float"
      write(unit) str,lf
      write(unit) BigEnd(real(mean%y, real32)),lf
      str="Z_COORDINATES"
      write(str(15:),'(i5,2x,a)') mean%nz,"float"
      write(unit) str,lf
      write(unit) BigEnd(real(mean%z, real32)),lf
      str="POINT_DATA"
      write(str(12:),*) mean%nx*mean%ny*mean%nz
      write(unit) str,lf 
    end subroutine
    
    subroutine save_buffer(status,title,sc,vec)
      integer,intent(in) :: status
      character(*),intent(in) :: title
      real(rp),intent(in) :: sc(:,:,:),vec(:,:,:,:)
      
      write(unit) title
      if (status==SCALAR) then
        write(unit) BigEnd(sc),lf
      else
        write(unit) BigEnd(vec),lf
      end if
    end subroutine
    
    subroutine mean_scalar(sc, mean_sc)
      real(vtk_rp),intent(in) :: sc(:,:,:)
      real(vtk_rp),intent(out) :: mean_sc(:,:,:)
      integer :: pi, pj, pk, ox, oy, oz

      mean_sc = 0
      do pk = 1, pz
        do pj = 1, py
          do pi = 1, px
            ox = (pi-1)*mean%nx
            oy = (pj-1)*mean%ny
            oz = (pk-1)*mean%nz
            mean_sc = mean_sc + sc(1+ox:mean%nx+ox, &
                                   1+oy:mean%ny+oy, &
                                   1+oz:mean%nz+oz)
          end do
        end do
      end do
      mean_sc = mean_sc / (px*py*pz)
    end subroutine

    subroutine mean_vector(vec, mean_vec)
      real(vtk_rp),intent(in) :: vec(:,:,:,:)
      real(vtk_rp),intent(out) :: mean_vec(:,:,:,:)
      integer :: pi, pj, pk, ox, oy, oz
      integer :: comp

      mean_sc = 0
      do pk = 1, pz
        do pj = 1, py
          do pi = 1, px
            ox = (pi-1)*mean%nx
            oy = (pj-1)*mean%ny
            oz = (pk-1)*mean%nz
            do comp = 1, size(vec, 1)
              mean_vec(comp,:,:,:) = mean_vec(comp,:,:,:) + vec(comp, &
                                                                1+ox:mean%nx+ox,&
                                                                1+oy:mean%ny+oy, &
                                                                1+oz:mean%nz+oz)
            end do
          end do
        end do
      end do
      mean_vec = mean_vec / (px*py*pz)
    end subroutine
  end subroutine compute_and_save
  

  function basename(path) result(res)
    character(:), allocatable :: res
    character(*), intent(in) :: path
    integer :: n, pos

    res = trim(path)

    n = len(res)
    if (res(n:n)=='/') res = res(1:n-1)

    pos = scan(res, "/", .true.)
    res = res(pos+1:)
  end function

  function dirname(path) result(res)
    character(:), allocatable :: res
    character(*), intent(in) :: path
    integer :: n, pos

    res = trim(path)

    n = len(res)
    if (res(n:n)=='/') res = res(1:n-1)

    pos = scan(res, "/", .true.)
    res = res(:pos-1)
  end function


end module Periodic


program vtk_periodic_mean
  use Periodic
  
  implicit none

  character(1024) :: input_file, arg

  integer :: periodsx, periodsy, periodsz

  call GetEndianness

  if (COMMAND_ARGUMENT_COUNT()>=2) then

    call get_command_argument(1, input_file)
  
    call get_command_argument(2, arg)

    read(arg,*) periodsx, periodsy, periodsz

  else

    stop "Usage: vtk_periodic_mean input_file perx,pery,perz"

  end if
    
  call compute_and_save(trim(input_file), periodsx, periodsy, periodsz)

  
end program
