module WorkKinds
  use iso_fortran_env
  
  integer,parameter :: rp = real32
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
  use WorkKinds
  use Endianness
  
  implicit none

  type grid
    integer nx,ny,nz
    integer unit
    real(rp),allocatable :: x(:),y(:),z(:)
    real(rp) :: x1, y1, z1
    real(rp) :: dx, dy, dz
    character(1024) :: fname
  contains
    procedure :: read_header
    procedure :: read_title
    procedure :: read_scalar
    procedure :: read_vector
  end type
  
  integer, parameter :: SCALAR = 1, VECTOR = 2
  character,parameter :: lf = achar(10)

  real(real32), parameter :: nan32 = transfer(-4194304, 1._real32)

contains
  

  subroutine read_header(g)
    use Endianness
    class(grid),intent(inout) :: g
    integer :: io, nx, ny, nz
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
    g%x = BigEnd(g%x)

    read(g%unit,pos=234+nx*4,iostat=io) ch5
    read(ch5,'(i5)') g%ny
    ny = g%ny

    allocate(g%y(g%ny))
    read(g%unit,pos=291+nx*4,iostat=io) g%y
    g%y = BigEnd(g%y)

    read(g%unit,pos=306+nx*4+ny*4,iostat=io) ch5
    read(ch5,'(i5)') g%nz
    nz = g%nz
    allocate(g%z(g%nz))
    read(g%unit,pos=363+nx*4+ny*4,iostat=io) g%z
    g%z = BigEnd(g%z)
    read(g%unit) ch,str
    read(g%unit) ch
    
    g%x1 = g%x(1)
    g%y1 = g%y(1)
    g%z1 = g%z(1)
    
    g%dx = (g%x(nx) - g%x(1)) / (nx-1)
    g%dy = (g%y(ny) - g%y(1)) / (ny-1)
    g%dz = (g%z(nz) - g%z(1)) / (nz-1)
   
    
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
    real(rp),intent(inout) :: buf(:,:,:)
    character :: ch
 
    read(g%unit) buf
    buf = BigEnd(buf)
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(inout) :: buf(:,:,:,:)
    character :: ch
    
    read(g%unit) buf
    buf = BigEnd(buf)
    read(g%unit) ch
  end subroutine
  
end module Types



program extend_periodic_x
  use WorkKinds
  use Types
  use Endianness
  
  implicit none
  
  type(grid) :: old, new

  real(rp), allocatable :: old_sc(:,:,:), old_vec(:,:,:,:)
  real(rp), allocatable :: new_sc(:,:,:), new_vec(:,:,:,:)

  character(1024) :: arg

  character(:), allocatable :: file_name, base_name, dir_name, title, bc_dir
  
  integer :: itx, ity, io, status, arg_len
  
  integer :: nx_times = 2, ny_times = 1

  call GetEndianness

  if (command_argument_count()<1) then
    write(*,*) "Usage: extend_periodic_x file_name nx_times,ny_times"
    stop 1
  end if
 
  call get_command_argument(1, length=arg_len)
  allocate(character(arg_len) :: file_name)
  
  call get_command_argument(1, value=file_name)
  base_name = trim(file_name(scan(file_name,'/',back=.true.) + 1 :  ))
  
  if (command_argument_count()>=2) then
    call get_command_argument(2, value=arg)
    read(arg,*,iostat=io) nx_times, ny_times
    if (io/=0) then
      write(*,*) "Error reading the number of repetitions."
      write(*,*) "Usage: extend_periodic_x file_name nx_times,ny_times"
      stop
    end if
  else
    nx_times = 2
    ny_times = 1
  end if

  dir_name = "extended"

  call execute_command_line("mkdir -p "//dir_name)
  
  old%fname = file_name

  call old%read_header

  new%nx = nx_times * old%nx
  new%ny = ny_times * old%ny
  new%nz = old%nz

  new%fname = dir_name // '/' // base_name

  allocate(new%x(new%nx))
  allocate(new%y(new%ny))
  new%x(1:old%nx) = old%x
  new%x(old%nx+1:) = [(new%x(old%nx)+itx*old%dx, itx = 1, new%nx-old%nx)]
  new%y(1:old%ny) = old%y
  new%y(old%ny+1:) = [(new%y(old%ny)+ity*old%dy, ity = 1, new%ny-old%ny)]
  new%z = old%z

  call save_header


  allocate(old_sc(old%nx,old%ny,old%nz), &
           old_vec(3,old%nx,old%ny,old%nz))
  allocate(new_sc(new%nx,new%ny,new%nz), &
           new_vec(3,new%nx,new%ny,new%nz))

  old_sc = 0
  old_vec = 0
  new_sc = 0
  new_vec = 0

  do
    call get_next(status,title)

    if (status==0) exit


    call get_buffer(status, old_sc, old_vec)
    
    if (status==SCALAR) then
      do ity = 1, ny_times
      do itx = 1, nx_times
        new_sc(1+(itx-1)*old%nx:itx*old%nx,1+(ity-1)*old%ny:ity*old%ny,:) = old_sc
      end do
      end do
    else if (status==VECTOR) then
      do ity = 1, ny_times
      do itx = 1, nx_times
        new_vec(:,1+(itx-1)*old%nx:itx*old%nx,1+(ity-1)*old%ny:ity*old%ny,:) = old_vec
      end do
      end do
    end if

    call save_buffer(status,title, new_sc, new_vec)
  end do
  
  close(old%unit)
  close(new%unit)


  deallocate(new_sc, new_vec, old_sc, old_vec)

  if (allocated(file_name)) deallocate(file_name)
  if (allocated(base_name)) deallocate(base_name)
  if (allocated(dir_name)) deallocate(dir_name)
  if (allocated(title)) deallocate(title)
  if (allocated(bc_dir)) deallocate(bc_dir)
  
contains

  subroutine get_next(status,title)
    integer :: status
    character(:),allocatable,intent(out) :: title

    call old%read_title(status,title)

  end subroutine
  
  subroutine get_buffer(status,sc,vec)
    integer,intent(in) :: status
    real(rp),intent(inout) :: sc(:,:,:),vec(:,:,:,:)
    
    if (status==SCALAR) then
      call old%read_scalar(sc)
    else
      call old%read_vector(vec)
    end if
  end subroutine
  
  subroutine save_header
    character(70) :: str
   
    open(newunit=new%unit,file=new%fname, &
      access='stream',status='replace',form="unformatted",action="write")

    write(new%unit) "# vtk DataFile Version 2.0",lf
    write(new%unit) "CLMM output file",lf
    write(new%unit) "BINARY",lf
    write(new%unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write(str(12:),*) new%nx, new%ny, new%nz
    write(new%unit) str,lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') new%nx,"float"
    write(new%unit) str,lf
    write(new%unit) BigEnd(real(new%x(1:new%nx), real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') new%ny,"float"
    write(new%unit) str,lf
    write(new%unit) BigEnd(real(new%y(1:new%ny), real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') new%nz,"float"
    write(new%unit) str,lf
    write(new%unit) BigEnd(real(new%z(1:new%nz), real32)),lf
    str="POINT_DATA"
    write(str(12:),*) new%nx * new%ny * new%nz
    write(new%unit) str,lf 
  end subroutine
  
  subroutine save_buffer(status,title,sc,vec)
    integer, intent(in) :: status
    character(*), intent(in) :: title
    real(rp), intent(inout) :: sc(:,:,:), vec(:,:,:,:)
    
    write(new%unit) title
    if (status==SCALAR) then
      sc = BigEnd(sc)
      write(new%unit) sc, lf
    else
      vec = BigEnd(vec)
      write(new%unit) vec, lf
    end if
  end subroutine

end program
