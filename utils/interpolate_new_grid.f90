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


module Interpolation
  use WorkKinds

  implicit none

contains

  pure real(rp) function TriLinInt(a, b, c, &
                                   vel000, vel100, vel010, vel001, vel110, vel101, vel011, vel111)
    real(rp), intent(in) :: a, b, c
    real(rp), intent(in) :: vel000, vel100, vel010, vel001, vel110, vel101, vel011, vel111

    TriLinInt =  (1-a) * (1-b) * (1-c) * vel000 + &
                 a     * (1-b) * (1-c) * vel100 + &
                 (1-a) * b     * (1-c) * vel010 + &
                 (1-a) * (1-b) * c     * vel001 + &
                 a     * b     * (1-c) * vel110 + &
                 a     * (1-b) * c     * vel101 + &
                 (1-a) * b     * c     * vel011 + &
                 a     * b     * c     * vel111

  end function TriLinInt

end module


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
    allocate(g%x(-2:g%nx+3))
    read(g%unit,pos=219,iostat=io) g%x(1:g%nx)
    g%x(1:g%nx) = BigEnd(g%x(1:g%nx))

    read(g%unit,pos=234+nx*4,iostat=io) ch5
    read(ch5,'(i5)') g%ny
    ny = g%ny

    allocate(g%y(-2:g%ny+3))
    read(g%unit,pos=291+nx*4,iostat=io) g%y(1:g%ny)
    g%y(1:g%ny) = BigEnd(g%y(1:g%ny))

    read(g%unit,pos=306+nx*4+ny*4,iostat=io) ch5
    read(ch5,'(i5)') g%nz
    nz = g%nz
    allocate(g%z(-2:g%nz+3))
    read(g%unit,pos=363+nx*4+ny*4,iostat=io) g%z(1:g%nz)
    g%z(1:g%nz) = BigEnd(g%z(1:g%nz))
    read(g%unit) ch,str
    read(g%unit) ch
    
    g%x1 = g%x(1)
    g%y1 = g%y(1)
    g%z1 = g%z(1)
    
    g%dx = (g%x(nx) - g%x(1)) / (nx-1)
    g%dy = (g%y(ny) - g%y(1)) / (ny-1)
    g%dz = (g%z(nz) - g%z(1)) / (nz-1)
 
    g%x(-2:0) = g%x(1:3) - 3*g%dx     
    g%y(-2:0) = g%y(1:3) - 3*g%dy     
    g%z(-2:0) = g%z(1:3) - 3*g%dz

    g%x(nx+1:nx+3) = g%x(nx-2:nx) + 3*g%dx     
    g%y(ny+1:ny+3) = g%y(ny-2:ny) + 3*g%dy     
    g%z(nz+1:nz+3) = g%z(nz-2:nz) + 3*g%dz     
  
    
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
    real(rp),intent(inout) :: buf(-2:,-2:,-2:)
    character :: ch
 
    read(g%unit) buf(1:g%nx, 1:g%ny, 1:g%nz)
    buf(1:g%nx, 1:g%ny, 1:g%nz) = BigEnd(buf(1:g%nx, 1:g%ny, 1:g%nz))
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(inout) :: buf(:,-2:,-2:,-2:)
    character :: ch
    
    read(g%unit) buf(:, 1:g%nx, 1:g%ny, 1:g%nz)
    buf(:, 1:g%nx, 1:g%ny, 1:g%nz) = BigEnd(buf(:, 1:g%nx, 1:g%ny, 1:g%nz))
    read(g%unit) ch
  end subroutine
  
end module Types


module InterpBounds

  use Parameters
  use Boundaries
  use ScalarBoundaries
  use Types

  implicit none

  logical :: enable_bc = .false.

  integer :: staggered = 0 !0..cell centres, 1..U, 2..V, 3..W

contains

  subroutine ReadBounds(dir, g)
    character(*), intent(in) :: dir
    type(grid), intent(in) :: g
    integer :: unit, io
    real(knd) :: Prandtl

    interface get
      procedure chget1
      procedure lget1, lget2, lget3
      procedure iget1, iget2, iget3
      procedure rget1, rget2, rget3
      procedure rgetv3
    end interface

    unit = 11

    open(unit,file=dir//"boundconds.conf",status="old",action="read")
    call get(Btype(We))
    call get(Btype(Ea))
    call get(Btype(So))
    call get(Btype(No))
    call get(Btype(Bo))
    call get(Btype(To))
    call get(sideU(1,So))
    call get(sideU(2,So))
    call get(sideU(3,So))
    call get(sideU(1,No))
    call get(sideU(2,No))
    call get(sideU(3,No))
    call get(sideU(1,Bo))
    call get(sideU(2,Bo))
    call get(sideU(3,Bo))
    call get(sideU(1,To))
    call get(sideU(2,To))
    call get(sideU(3,To))
    call get(z0W)
    call get(z0E)
    call get(z0S)
    call get(z0N)
    call get(z0B)
    call get(z0T)
    close(unit)

    open(unit,file=dir//"thermal.conf",status="old",action="read",iostat = io)
    if (io==0) then
      call get(enable_buoyancy)
      call get(Prandtl)
      call get(grav_acc)
      call get(temperature_ref)
      call get(TempBtype(We))
      call get(TempBtype(Ea))
      call get(TempBtype(So))
      call get(TempBtype(No))
      call get(TempBtype(Bo))
      call get(TempBtype(To))
      call get(sideTemp(We))
      call get(sideTemp(Ea))
      call get(sideTemp(So))
      call get(sideTemp(No))
      call get(sideTemp(Bo))
      call get(sideTemp(To))
      close(unit)
    else
      enable_buoyancy = .false.
    end if

    open(unit,file=dir//"moisture.conf",status="old",action="read",iostat = io)
    if (io==0) then
      call get(enable_moisture)
      call get(moisture_ref)
      call get(MoistBtype(We))
      call get(MoistBtype(Ea))
      call get(MoistBtype(So))
      call get(MoistBtype(No))
      call get(MoistBtype(Bo))
      call get(MoistBtype(To))
      call get(sideMoist(We))
      call get(sideMoist(Ea))
      call get(sideMoist(So))
      call get(sideMoist(No))
      call get(sideMoist(Bo))
      call get(sideMoist(To))
      close(unit)
    else
      enable_moisture = .false.
    end if

   open(unit,file=dir//"scalars.conf",status="old",action="read",iostat = io)
   if (io==0) then
     call get(num_of_scalars)
     call get(computedeposition)
     call get(computegravsettling)
     call get(partdistrib)
     call get(totalscalsource)
     call get(scalsourcetype)

     call get(ScalBtype(We))
     call get(ScalBtype(Ea))
     call get(ScalBtype(So))
     call get(ScalBtype(No))
     call get(ScalBtype(Bo))
     call get(ScalBtype(To))
     call get(sideScal(We))
     call get(sideScal(Ea))
     call get(sideScal(So))
     call get(sideScal(No))
     call get(sideScal(Bo))
     call get(sideScal(To))
     close(unit)
   else
     num_of_scalars = 0
   end if

   where (Btype==BC_TURBULENT_INLET.or.Btype==BC_INLET_FROM_FILE) Btype = BC_NEUMANN

   if (staggered == 1) then
     Prnx = g%nx
     Prny = g%ny
     Prnz = g%nz
     Unx = g%nx
     Uny = g%ny
     Unz = g%nz
   else if (staggered == 2) then
     Prnx = g%nx
     Prny = g%ny
     Prnz = g%nz
     Vnx = g%nx
     Vny = g%ny
     Vnz = g%nz
   else if (staggered == 2) then
     Prnx = g%nx
     Prny = g%ny
     Prnz = g%nz
     Vnx = g%nx
     Vny = g%ny
     Vnz = g%nz
   else
     Prnx = g%nx
     Prny = g%ny
     Prnz = g%nz

     Unx = Prnx
     Uny = Prny
     Unz = Prnz

     Vnx = Prnx
     Vny = Prny
     Vnz = Prnz

     Wnx = Prnx
     Wny = Prny
     Wnz = Prnz
   end if

   allocate(Uin(-2:Uny+3,-2:Unz+3),Vin(-2:Vny+3,-2:Vnz+3),Win(-2:Wny+3,-2:Wnz+3))
   Uin = 0; Vin = 0; Win = 0

   allocate(Viscosity(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
   Viscosity = 1e5
   allocate(TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
   TDiff = 1e5

   allocate(dxU(-2:Unx+3), dyV(-2:Vny+3), dzW(-2:Wnz+3))
   dxU = g%dx
   dyV = g%dy
   dzW = g%dz

  contains

     subroutine chget1(x)
       character(*),intent(out) :: x
       read(unit,fmt='(/)')
       read(unit,*) x
     end subroutine
     subroutine lget1(x)
       logical,intent(out) :: x
       character(120) :: line, fname
       integer :: ierr, tmp
       read(unit,fmt='(/)')
       read(unit,'(a)') line
       read(line,*,iostat=ierr) x
       if (ierr/=0) then
         read(line,*,iostat=ierr) tmp
         if (ierr/=0) then
           inquire(unit,name=fname)
           if (master) write(*,*) "Stop expected a boolean flag in file "//trim(fname)
           if (master) write(*,*) "Received '"//trim(line)//"' instead."
           call error_stop
         end if
         x = tmp /=0
       end if
     end subroutine
     subroutine lget2(x,y)
       logical,intent(out) :: x,y
       read(unit,fmt='(/)')
       read(unit,*) x,y
     end subroutine
     subroutine lget3(x,y,z)
       logical,intent(out) :: x,y,z
       read(unit,fmt='(/)')
       read(unit,*) x,y,z
     end subroutine
     subroutine iget1(x)
       integer,intent(out) :: x
       read(unit,fmt='(/)')
       read(unit,*) x
     end subroutine
     subroutine iget2(x,y)
       integer,intent(out) :: x,y
       read(unit,fmt='(/)')
       read(unit,*) x,y
     end subroutine
     subroutine iget3(x,y,z)
       integer,intent(out) :: x,y,z
       read(unit,fmt='(/)')
       read(unit,*) x,y,z
     end subroutine
     subroutine rget1(x)
       real(knd),intent(out) :: x
       read(unit,fmt='(/)')
       read(unit,*) x
     end subroutine
     subroutine rget2(x,y)
       real(knd),intent(out) :: x,y
       read(unit,fmt='(/)')
       read(unit,*) x,y
     end subroutine
     subroutine rget3(x,y,z)
       real(knd),intent(out) :: x,y,z
       read(unit,fmt='(/)')
       read(unit,*) x,y,z
     end subroutine
     subroutine rgetv3(v)
       real(knd),intent(out) :: v(3)
       read(unit,fmt='(/)')
       read(unit,*) v
     end subroutine
  end subroutine


end module



program interpolate_new_grid
  use WorkKinds
  use Types
  use Endianness
  use Interpolation
  use Strings
  use InterpBounds
  
  implicit none
  
  type(grid) :: old, new

  real(rp), allocatable :: old_sc(:,:,:), old_vec(:,:,:,:)
  real(rp), allocatable :: new_sc(:,:,:), new_vec(:,:,:,:)

  character(1024) :: arg

  character(:), allocatable :: file_name, base_name, dir_name, title, bc_dir
  
  integer :: it, io, status, arg_len

  real(rp) :: lo, up

  call GetEndianness

  if (command_argument_count()<2) then
    write(*,*) "Usage: interpolate_new_grid file_name nx,ny,nz"
    write(*,*) "nx,ny,nz are the numbers of cells in each dimension."
    stop 1
  end if
 
  call get_command_argument(1, length=arg_len)
  allocate(character(arg_len) :: file_name)
  call get_command_argument(1, value=file_name)
  base_name = trim(file_name(scan(file_name,'/',back=.true.) + 1 :  ))

  select case (base_name(1:1))
    case ('U')
      staggered = 1
    case ('V')
      staggered = 2
    case ('W')
      staggered = 3
    case default
      staggered = 0
  end select

  call get_command_argument(2, value=arg)
  read(arg,*,iostat=io) new%nx, new%ny, new%nz
  
  if (new%nx<=0 .or. new%ny<=0 .or. new%nz<=0) then
    write(*,*) "Error, supply three positive gride sizes as 'interpolate_new_grid file_name nx,ny,nz'."
  end if
  
  if (command_argument_count()>=3) then
    call get_command_argument(3, value=arg)
    enable_bc = arg == "-bc"
  end if
  
  if (enable_bc.and.command_argument_count()>=4) then
    call get_command_argument(4, value=arg)
    bc_dir = trim(arg)
  end if
  
  dir_name = "grid-"//itoa(new%nx)//"-"//itoa(new%ny)//"-"//itoa(new%nz)

  call execute_command_line("mkdir -p "//dir_name)
  
  old%fname = file_name

  call old%read_header



  if (enable_bc) call ReadBounds(bc_dir, old)



  new%fname = dir_name // '/' // base_name


  
  if (staggered==1) then
    if (enable_bc.and.Btype(Ea)==BC_PERIODIC) then
      lo = old%x(1) - old%dx
      up = old%x(old%nx)
      new%dx = (up - lo) / new%nx
    else
      new%nx = new%nx - 1

      lo = old%x(1) - old%dx
      up = old%x(old%nx) + old%dx
      new%dx = (up - lo) / (new%nx+1)
    end if

    allocate(new%x(-2: new%nx+3))
    new%x(:) = [ (lo + new%dx * it, it = -2, new%nx+3) ]
    new%x1 = new%x(1)
  else
    lo = old%x(1) - old%dx / 2
    up = old%x(old%nx) + old%dx / 2

    new%dx = (up - lo) / new%nx
    allocate(new%x(-2: new%nx+3))
    new%x(:) = [ (lo + new%dx * (it - 0.5_rp), it = -2, new%nx+3) ]
    new%x1 = new%x(1)
  end if
  
  if (staggered==2) then
    if (enable_bc.and.Btype(No)==BC_PERIODIC) then
      lo = old%y(1) - old%dy
      up = old%y(old%ny)
      new%dy = (up - lo) / new%ny
    else
      new%ny = new%ny - 1

      lo = old%y(1) - old%dy
      up = old%y(old%ny) + old%dy 
      new%dy = (up - lo) / (new%ny+1)
    end if

    allocate(new%y(-2: new%ny+3))
    new%y(:) = [ (lo + new%dy * it, it = -2, new%ny+3) ]
    new%y1 = new%y(1)
  else
    lo = old%y(1) - old%dy / 2
    up = old%y(old%ny) + old%dy / 2
    
    new%dy = (up - lo) / new%ny
    allocate(new%y(-2: new%ny+3))
    new%y(:) = [ (lo + new%dy * (it - 0.5_rp), it = -2, new%ny+3) ]
    new%y1 = new%y(1)
  end if

  if (staggered==3) then
    if (enable_bc.and.Btype(To)==BC_PERIODIC) then
      lo = old%z(1) - old%dz
      up = old%z(old%nz)
      new%dz = (up - lo) / new%nz
    else
      new%nz = new%nz - 1

      lo = old%z(1) - old%dz
      up = old%z(old%nz) + old%dz
      new%dz = (up - lo) / (new%nz+1)
    end if

    allocate(new%z(-2: new%nz+3))
    new%z(:) = [ (lo + new%dz * it, it = -2, new%nz+3) ]
    new%z1 = new%z(1)
  else  
    lo = old%z(1) - old%dz / 2
    up = old%z(old%nz) + old%dz / 2

    new%dz = (up - lo) / new%nz
    allocate(new%z(-2: new%nz+3))
    new%z(:) = [ (lo + new%dz * (it - 0.5_rp), it = -2, new%nz+3) ]
    new%z1 = new%z(1)
  end if
  
  call save_header


  allocate(old_sc(-2:old%nx+3,-2:old%ny+3,-2:old%nz+3), &
           old_vec(3,-2:old%nx+3,-2:old%ny+3,-2:old%nz+3))
  allocate(new_sc(-2:new%nx+3,-2:new%ny+3,-2:new%nz+3), &
           new_vec(3,-2:new%nx+3,-2:new%ny+3,-2:new%nz+3))

  old_sc = 0
  old_vec = 0
  new_sc = 0
  new_vec = 0

  do
    call get_next(status,title)

    if (status==0) exit


    call get_buffer(status, old_sc, old_vec)


    if (enable_bc) then
      if (staggered==1.and.status==SCALAR.and.starts_with_ins(title,"scalars u")) then
        Uin(:,:) = old_sc(1,:,:)
        call BoundU(1,old_sc,Uin)
      else if (staggered==2.and.status==SCALAR.and.starts_with_ins(title,"scalars v")) then
        Vin(:,:) = old_sc(1,:,:)
        call BoundU(2,old_sc,Vin)
      else if (staggered==3.and.status==SCALAR.and.starts_with_ins(title,"scalars w")) then
        Win(:,:) = old_sc(1,:,:)
        call BoundU(3,old_sc,Win)
      else if (staggered==0.and.status==SCALAR.and.starts_with_ins(title,"scalars temperature")) then
        call BoundTemperature(old_sc)
      else if (staggered==0.and.status==SCALAR.and.starts_with_ins(title,"scalars moisture")) then
        call BoundMoisture(old_sc)
      else if (staggered==0.and.status==SCALAR.and.starts_with_ins(title,"scalars scalar")) then
        call BoundScalar(old_sc)
      else if (staggered==0.and.status==SCALAR) then
        call bound_simple(old_sc)
      else if (staggered==0.and.status==VECTOR) then
        call bound_simple(old_vec(1,:,:,:))
        call bound_simple(old_vec(2,:,:,:))
        call bound_simple(old_vec(3,:,:,:))
      end if
    else
      if (status==SCALAR) then
        call bound_simple(old_sc)
      else if (status==VECTOR) then
        call bound_simple(old_vec(1,:,:,:))
        call bound_simple(old_vec(2,:,:,:))
        call bound_simple(old_vec(3,:,:,:))
      end if
    end if

    call interpolate_buffer(status, new_sc, new_vec, old_sc, old_vec)

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
    real(rp),intent(inout) :: sc(-2:,-2:,-2:),vec(:,-2:,-2:,-2:)
    
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
    real(rp), intent(inout) :: sc(-2:,-2:,-2:), vec(:,-2:,-2:,-2:)
    
    write(new%unit) title
    if (status==SCALAR) then
      sc = BigEnd(sc)
      write(new%unit) sc(1:new%nx, 1:new%ny, 1:new%nz),lf
    else
      vec = BigEnd(vec)
      write(new%unit) vec(:,1:new%nx, 1:new%ny, 1:new%nz),lf
    end if
  end subroutine

  subroutine old_index(x, y, z, xi, yj, zk)
    real(rp), intent(in) :: x, y, z
    integer, intent(out) :: xi, yj, zk

    xi = min(max(floor( (x - old%x(1))/old%dx )+1, 0), old%nx)
    yj = min(max(floor( (y - old%y(1))/old%dy )+1, 0), old%ny)
    zk = min(max(floor( (z - old%z(1))/old%dz )+1, 0), old%nz)
  end subroutine

  function interpolate_trilinear(arr, x, y, z) result(res)
    real(rp) :: res
    real(rp), intent(in) :: arr(-2:,-2:,-2:)
    real(rp), intent(in) :: x, y, z
    integer :: xi, yj, zk

    call old_index(x, y, z, xi, yj, zk)

    res = TriLinInt((x - old%x(xi)) / old%dx, &
                    (y - old%y(yj)) / old%dy, &
                    (z - old%z(zk)) / old%dz, &
                    arr(xi, yj, zk), &
                    arr(xi+1, yj  , zk  ), &
                    arr(xi  , yj+1, zk  ), &
                    arr(xi  , yj  , zk+1), &
                    arr(xi+1, yj+1, zk  ), &
                    arr(xi+1, yj  , zk+1), &
                    arr(xi  , yj+1, zk+1), &
                    arr(xi+1, yj+1, zk+1))
  end function

  subroutine interpolate_scalar(new_sc, old_sc)
    real(rp), intent(out) :: new_sc(-2:,-2:,-2:)
    real(rp), intent(in)  :: old_sc(-2:,-2:,-2:)
    integer :: i, j, k
    do k = 1, new%nz
      do j = 1, new%ny
        do i = 1, new%nx
          new_sc(i,j,k) = interpolate_trilinear(old_sc, new%x(i), new%y(j), new%z(k))
        end do
      end do
    end do
  end subroutine

  subroutine interpolate_vector(new_vec, old_vec)
    real(rp),intent(out) :: new_vec(:,-2:,-2:,-2:)
    real(rp),intent(in)  :: old_vec(:,-2:,-2:,-2:)
    integer :: i, j, k, comp
    do k = 1, new%nz
      do j = 1, new%ny
        do i = 1, new%nx
          do comp = 1, 3
            new_vec(comp,i,j,k) = interpolate_trilinear(old_vec(comp,:,:,:), new%x(i), new%y(j), new%z(k))
          end do
        end do
      end do
    end do
  end subroutine

  subroutine interpolate_scalar_spline(new_sc, old_sc)
    use bspline_oo_module
    real(rp), intent(out) :: new_sc(-2:,-2:,-2:)
    real(rp), intent(in)  :: old_sc(-2:,-2:,-2:)
    integer :: i, j, k
    type(bspline_3d) :: spl
    integer :: iflag, idx, idy, idz
    real(real64) :: val

    idx = 0; idy = 0; idz = 0
    
    call spl%initialize(real(old%x(-1:old%nx+2), real64), &
                        real(old%y(-1:old%ny+2), real64), &
                        real(old%z(-1:old%nz+2), real64), &
                        real(old_sc(-1:old%nx+2,-1:old%ny+2,-1:old%nz+2), real64), &
                        3, 3, 3, iflag)
    if (iflag/=1) then
      write(*,*) "Error in spline initialization, iflag:",iflag
      stop
    end if
    
    
    do k = 1, new%nz
      do j = 1, new%ny
        do i = 1, new%nx
          call spl%evaluate(real(new%x(i), real64), real(new%y(j), real64), real(new%z(k), real64), &
                            idx, idy, idz, val, iflag)

          new_sc(i,j,k) = real(val, rp)
        end do
      end do
    end do

    call spl%destroy()
  end subroutine

  subroutine interpolate_vector_spline(new_vec, old_vec)
    use bspline_oo_module
    real(rp),intent(out) :: new_vec(:,-2:,-2:,-2:)
    real(rp),intent(in)  :: old_vec(:,-2:,-2:,-2:)
    integer :: i, j, k, comp
    type(bspline_3d) :: spl(3)
    integer :: iflag, idx(3), idy(3), idz(3)
    real(real64) :: val(3)

    idx = 0; idy = 0; idz = 0

    do comp = 1, 3
      call spl(comp)%initialize(real(old%x(-1:old%nx+2), real64), &
                                real(old%y(-1:old%ny+2), real64), &
                                real(old%z(-1:old%nz+2), real64), &
                                real(old_vec(comp,-1:old%nx+2,-1:old%ny+2,-1:old%nz+2), real64), &
                                3, 3, 3, iflag)
      if (iflag/=1) then
        write(*,*) "Error in spline initialization, iflag:",iflag
        stop
      end if
    end do
    
    do k = 1, new%nz
      do j = 1, new%ny
        do i = 1, new%nx
          do comp = 1, 3
            call spl(comp)%evaluate(real(new%x(i), real64), real(new%y(j), real64), real(new%z(k), real64), &
                                    idx(comp), idy(comp), idz(comp), val(comp), iflag)

            new_vec(comp,i,j,k) = real(val(comp), rp)
          end do
        end do
      end do
    end do

    do comp = 1, 3
      call spl(comp)%destroy()
    end do
  end subroutine

  subroutine interpolate_buffer(status, new_sc, new_vec, old_sc, old_vec)
    integer,intent(in) :: status
    real(rp),intent(out) :: new_sc(:,:,:), new_vec(:,:,:,:)
    real(rp),intent(in) :: old_sc(:,:,:), old_vec(:,:,:,:)
    
    if (status==SCALAR) then
      call interpolate_scalar_spline(new_sc, old_sc)
    else
      call interpolate_vector_spline(new_vec, old_vec)
    end if
  end subroutine

  function starts_with_ins(ch1, ch2) result(res)
    logical :: res
    character(*), intent(in) :: ch1, ch2
    integer :: l1, l2
    l1 = len_trim(ch1)
    l2 = len_trim(ch2)
    res = downcase(ch1(1:min(l1,l2))) == downcase(ch2(1:l2))
  end function

  subroutine bound_simple(a)
    real(knd),intent(inout) :: a(-2:,-2:,-2:)
    integer i,j,k,nx,ny,nz

    nx = old%nx
    ny = old%ny
    nz = old%nz
    
    if (Btype(Ea)==BC_PERIODIC) then
      do k = 1,nz
       do j = 1,ny
         a(-2:0,j,k) = a(nx-2:nx,j,k)
         a(nx+1:nx+3,j,k) = a(1:3,j,k)
       end do
      end do
    else
      do k = 1,nz
       do j = 1,ny
         a(-2:0,j,k) = a(3:1:-1,j,k)
         a(nx+1:nx+3,j,k) = a(nx:nx-2:-1,j,k)
       end do
      end do
    end if

    if (Btype(No)==BC_PERIODIC) then
      do k = 1,nz
       do i = -2,nx+3
         a(i,-2:0,k) = a(i,ny-2:ny,k)
         a(i,ny+1:ny+3,k) = a(i,1:3,k)
       end do
      end do
    else
      do k = 1,nz
       do i = -2,nx+3
         a(i,-2:0,k) = a(i,3:1:-1,k)
         a(i,ny+1:ny+3,k) = a(i,ny:ny-2:-1,k)
       end do
      end do
    end if
    
    if (Btype(To)==BC_PERIODIC) then
      do j = -2,ny+3
       do i = -2,nx+3
         a(i,j,-2:0) = a(i,j,nz-2:nz)
         a(i,j,nz+1:nz+3) = a(i,j,1:3)
       end do
      end do
    else
      do j = -2,ny+3
       do i = -2,nx+3
         a(i,j,-2:0) = a(i,j,3:1:-1)
         a(i,j,nz+1:nz+3) = a(i,j,nz:nz-2:-1)
       end do
      end do
    end if
    
  end subroutine bound_simple


end program
