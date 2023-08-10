module Kinds
  use iso_fortran_env
  
  integer,parameter :: rp = real32
end module

module Endianness

  use iso_fortran_env
  
  implicit none

  private
  public :: GetEndianness, BigEnd, SwapB, littleendian

  logical, save :: littleendian = .true.

  interface BigEnd
    module procedure BigEnd16
    module procedure BigEnd32
    module procedure BigEnd32_a1
    module procedure BigEnd32_a2
    module procedure BigEnd32_a3
    module procedure BigEnd32_a4
    module procedure BigEnd64
    module procedure BigEnd64_a1
    module procedure BigEnd64_a2
    module procedure BigEnd64_a3
    module procedure BigEnd64_a4
  end interface

  interface SwapB
    module procedure SwapB32
    module procedure SwapB64
  end interface

  contains

    subroutine GetEndianness
      character(4) :: bytes !may not work on some processors

      bytes = transfer(1_int32,bytes)
      if (ichar(bytes(4:4))==1) then
        littleendian=.false.
      else
        littleendian=.true.
      endif
    end subroutine GetEndianness

   
    elemental function BigEnd16(x) result(res)
      integer(int16) :: res
      integer(int16),intent(in)::x
      character(2) :: bytes
      
      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes)
        res = ichar(bytes(2:2),int16)
        res = ior( ishft(ichar(bytes(1:1),int16),8), res )
      endif
    end function
    
    function BigEnd32(x) result(res)
      real(real32),intent(in) :: x
      real(real32) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a1(x) result(res)
      real(real32),intent(in) :: x(:)
      real(real32) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a2(x) result(res)
      real(real32),intent(in) :: x(:,:)
      real(real32) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a3(x) result(res)
      real(real32),intent(in) :: x(:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a4(x) result(res)
      real(real32),intent(in) :: x(:,:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64(x) result(res)
      real(real64),intent(in) :: x
      real(real64) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a1(x) result(res)
      real(real64),intent(in) :: x(:)
      real(real64) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a2(x) result(res)
      real(real64),intent(in) :: x(:,:)
      real(real64) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a3(x) result(res)
      real(real64),intent(in) :: x(:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a4(x) result(res)
      real(real64),intent(in) :: x(:,:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    elemental function SwapB32(x) result(res)
      real(real32) :: res
      real(real32),intent(in) :: x
      character(4) :: bytes
      integer(int32) :: t
      real(real32) :: rbytes, rt
      !nicer looking TRANSFER is problematic to optimize by ifort
      ! gfortran would be OK with that
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(4:4),int32)

      t = ior( ishftc(ichar(bytes(3:3),int32),8),  t )

      t = ior( ishftc(ichar(bytes(2:2),int32),16), t )

      t = ior( ishftc(ichar(bytes(1:1),int32),24), t )

      res = rt
        
    end function
    
    elemental function SwapB64(x) result(res)
      real(real64) :: res
      real(real64),intent(in) :: x
      character(8) :: bytes
      integer(int64) :: t
      real(real64) :: rbytes, rt
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(8:8),int64)

      t = ior( ishftc(ichar(bytes(7:7),int64),8),  t )

      t = ior( ishftc(ichar(bytes(6:6),int64),16), t )

      t = ior( ishftc(ichar(bytes(5:5),int64),24), t )

      t = ior( ishftc(ichar(bytes(4:4),int64),32), t )

      t = ior( ishftc(ichar(bytes(3:3),int64),40), t )

      t = ior( ishftc(ichar(bytes(2:2),int64),48), t )

      t = ior( ishftc(ichar(bytes(1:1),int64),56), t )

      res = rt

    end function

end module Endianness




module Types
  use Kinds
  use Endianness
  
  implicit none

  type grid
    type(grid), pointer :: g_out => null()
    integer nx,ny,nz
    integer offx, offy, offz
    integer hu, fu
    integer header_pos
    real(rp),allocatable :: x(:),y(:),z(:)
    character(1024) :: fname
    real(rp), pointer, contiguous :: vec(:,:,:,:)
    real(rp), pointer, contiguous :: sc(:,:,:) => null() ! to save space in memory
  contains
    procedure :: process_next
    procedure :: read_var_from_header
    procedure :: read_var_from_eaf, save_var_to_eaf
    procedure :: save_header
    procedure :: open_frame, create_frame
    procedure :: rewind_header
    procedure :: close_header, close_frame
    final :: finalize
  end type
  
  interface grid
    procedure grid_init
  end interface
  
  ! in the efh header
  integer, parameter :: scalar_flag = 1001, vector_flag = 1002
  ! local here
  integer, parameter :: SCALAR = 1, VECTOR = 2, EOF=-1
  
  character,parameter :: lf = achar(10)

  real(real32), parameter :: nan32 = transfer(-4194304, 1._real32)

contains
  

  function grid_init(fname) result (g)
    use Endianness
    type(grid), target :: g
    character(*) :: fname
    integer :: io
    integer(int32) :: itmp, nxs(3)
    character(256) :: msg

    g%fname = fname
    
    open(newunit=g%hu,file=g%fname,access="stream",form="unformatted",status="old",action="read",iostat=io, iomsg=msg)
    if (io/=0) then
      stop msg
    end if

    read(g%hu) itmp
    
    if (itmp/=1) then
      stop "Running vtk_subset on a computer with different endianness not yet supported."
    end if
    
    read(g%hu) itmp
    if (itmp /= storage_size(1_rp)/CHARACTER_STORAGE_SIZE) then
      write(*,*) "vtk_subset compiled for different real precision.Storage sizes:"
      write(*,*) "  compiled:", storage_size(1_rp)/CHARACTER_STORAGE_SIZE
      write(*,*) "  data file uses:", itmp
      stop
    end if
      
    read(g%hu) nxs
    
    allocate(g%x(nxs(1)))
    allocate(g%y(nxs(2)))
    allocate(g%z(nxs(3)))
    
    g%nx = nxs(1)
    g%ny = nxs(2)
    g%nz = nxs(3)
    
    read(g%hu) g%x, g%y, g%z
    
    allocate(g%vec(3,g%nx, g%ny, g%nz))
    
    g%sc(1:g%nx,1:g%ny,1:g%nz) => g%vec
    
    inquire(unit=g%hu, pos=g%header_pos)

  end function

  
  subroutine read_var_from_header(g, title, status)
    class(grid),intent(in) :: g
    integer(int32) :: flag
    character(len=:), allocatable, intent(out) :: title
    integer, intent(out) :: status
    character(1) :: ch
    integer :: ios
    
    title = ""
    
    read(g%hu, iostat=ios) flag
    
    if(ios==IOSTAT_END) then
      status = EOF
      return
    else if (ios/=0) then
      stop "Error reading the next item from the header."
    end if
    
    if (flag==scalar_flag) then
      status = SCALAR
    else if (flag==vector_flag) then
      status = VECTOR
    else  
      stop "Unknown flag for the variable type in the header file."
    end if
    
    do
      read(g%hu) ch
      if (iachar(ch)==10) exit
      title = title // ch
    end do
  end subroutine
  
  subroutine read_var_from_eaf(g, status)
    class(grid),intent(inout) :: g
    integer, intent(in) :: status

    if (status == SCALAR) then
      read(g%fu) g%sc
    else
      read(g%fu) g%vec
    end if
  end subroutine

  subroutine save_var_to_eaf(g, status)
    class(grid),intent(inout) :: g
    integer, intent(in) :: status

    associate(o => g%g_out)
      if (status == SCALAR) then
        write(o%fu) g%sc(1+o%offx:o%nx+o%offx,1+o%offy:o%ny+o%offy,1+o%offz:o%nz+o%offz)
      else
        write(o%fu) g%vec(:,1+o%offx:o%nx+o%offx,1+o%offy:o%ny+o%offy,1+o%offz:o%nz+o%offz)
      end if
    end associate
  end subroutine

  subroutine save_header(g)
    use iso_c_binding
    
    class(grid), intent(inout) :: g
    
    real(rp), allocatable :: xs(:), ys(:), zs(:)
    
    integer(int32), parameter :: scalar_flag = 1001, vector_flag = 1002
    
    character(2512) :: file_name
    character :: ch
    integer :: io
    
    file_name = ""

    file_name = g%g_out%fname
    
    open(newunit=g%g_out%hu,file = file_name, &
      access='stream',status='replace',form="unformatted",action="write")
      
    write(g%g_out%hu) 1_int32                                         ! endianness
    write(g%g_out%hu) int(storage_size(1_real32) / &
                        CHARACTER_STORAGE_SIZE, &
                      int32)                                      ! size of real used

    xs = g%g_out%x
    ys = g%g_out%y
    zs = g%g_out%z
    
    write(g%g_out%hu) int(g%g_out%nx, int32), int(g%g_out%ny, int32), int(g%g_out%nz, int32)
    write(g%g_out%hu) real(xs,real32)
    write(g%g_out%hu) real(ys,real32)
    write(g%g_out%hu) real(zs,real32)
    
    call g%rewind_header
    
    do
      read(g%hu, iostat=io) ch
      if (io/=0) exit
      write(g%g_out%hu) ch
    end do
    
    close(g%g_out%hu)
    
    
    call g%rewind_header
  end subroutine

  
  subroutine process_next(g,status)
    class(grid), intent(inout) ::g
    integer, intent(inout) :: status
    character(:), allocatable :: title
    
    call g%read_var_from_header(title, status)
   
    if (status /= EOF) then
      call g%read_var_from_eaf(status)
      call g%save_var_to_eaf(status)
    end if
  end subroutine
  
  subroutine open_frame(g, fname)
    class(grid), intent(inout) :: g
    character(*) :: fname
    integer :: io
    character(256) :: msg
    open(newunit=g%fu,file=fname,access="stream",form="unformatted",status="old", &
         action="read",iostat=io, position="rewind", iomsg=msg)
    if (io/=0) then
      stop msg
    end if
  end subroutine
  
  subroutine create_frame(g, fname)
    class(grid), intent(inout) :: g
    character(*) :: fname
    integer :: io
    character(256) :: msg
    open(newunit=g%fu,file=fname,access="stream",form="unformatted",status="replace", &
         action="write",iostat=io, position="rewind", iomsg=msg)
    if (io/=0) then
      stop msg
    end if
  end subroutine
  
  subroutine finalize(g)
    type(grid) :: g
    deallocate(g%vec)
  end subroutine
  
  subroutine rewind_header(g)
    class(grid), intent(inout) :: g
    read(g%hu,pos=g%header_pos)
  end subroutine
  
  subroutine close_header(g)
    class(grid), intent(inout) :: g
    close(g%hu)
  end subroutine
  
  subroutine close_frame(g)
    class(grid), intent(inout) :: g
    close(g%fu)
  end subroutine
  
end module Types





program eaf_subset
  use Types
  use Endianness
  
  implicit none
  
  integer :: nx, ny, nz
 
  integer :: status, arg_len
  
  character(:), allocatable :: in_name, out_name, arg
  
  character(256), allocatable :: fields_arr(:)
  
  integer :: i1=1, i2=1, j1=1, j2=1, k1=1, k2=1, iarg
  
  call GetEndianness

  if (command_argument_count()<2) then
    write(*,*) "Usage: eaf_subset in_name out_name [extent] [field]"
    stop 1
  end if
 
  call get_command_argument(1, length=arg_len)
  allocate(character(arg_len) :: in_name)
  call get_command_argument(1, value=in_name)

  call get_command_argument(2, length=arg_len)
  allocate(character(arg_len) :: out_name)
  call get_command_argument(2, value=out_name)
  
  if (in_name==out_name) then
    write(*,*) "Error: in_name and out_name are identical"
    stop
  end if
  
  nx = -1
  ny = -1
  nz = -1  
  allocate(fields_arr(0))
  
  do iarg = 3,4
    if (command_argument_count() >= iarg)  then

      call get_command_argument(iarg, length=arg_len)
      if (allocated(arg)) deallocate(arg)
      allocate(character(arg_len) :: arg)
      call get_command_argument(iarg, value=arg)
      
      if (arg(1:1)>='0' .and. arg(1:1)<='9') then
        call parse_extent(arg)
      else
        call parse_fields(arg)
      end if
    end if
  end do

  if (i1<1 .or. j1<1 .or. k1<1) then 
    write(*,*) "Error, one of the extents starts below 1."
    write(*,*) "x extent start:", i1
    write(*,*) "y extent start:", j1
    write(*,*) "z extent start:", k1
    stop 1
  end if
  
  if (i1>i2 .or. j1>j2 .or. k1>k2) then 
    write(*,*) "Error, one of the extents is negative."
    write(*,*) "x extent:", i1, i2
    write(*,*) "y extent:", j1, j2
    write(*,*) "z extent:", k1, k2
    stop 1
  end if
    
  
  call process_series(in_name, out_name)

contains

  subroutine parse_extent(extent)
    character(*), intent(in) :: extent
    integer :: io
    
    read(extent,*, iostat=io) i1,i2,j1,j2,k1,k2
    
    if (io/=0) then
      write(*,*) "Eror parsing the extent from '",trim(extent),"'."
      stop 1
    end if
    
    nx = i2 - i1 + 1
    ny = j2 - j1 + 1
    nz = k2 - k1 + 1    
  end subroutine
  
  
  subroutine parse_fields(fields)
    character(*), intent(in) :: fields
    
    fields_arr = split(fields,",")
  end subroutine

  function split(str,sep) result(res)
    character(256), allocatable :: res(:)
    character(*), intent(in) :: str
    character, intent(in) :: sep
    integer :: i, n, pos, istr
    n = 1
    do i = 1, len(str)
      if (str(i:i)==sep) then
        n = n + 1
      end if
    end do
    
    allocate(res(n))
    res = ""
    
    pos = 1
    do istr = 1,n
      i = pos
      do
        if (str(i:i)==sep) then
          res(istr) = str(pos:i-1)
          pos = i + 1
          exit
        else if (i==len(str)) then
          res(istr) = str(pos:i)
          pos = i
          exit
        end if
        i = i + 1
      end do
    end do
  end function
  
  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function    

  subroutine process_series(series_name, out_name)
    character(*), intent(in) :: series_name, out_name
    character(:), allocatable :: fname, base_name
    type(grid), target :: g
    logical :: ex
    integer :: file_num
    
    g = grid(series_name//".efh")
    
    if (i1<1 .or. j1<1 .or. k1<1) then 
      write(*,*) "Error, one of the extents starts below 1."
      write(*,*) "x extent start:", i1
      write(*,*) "y extent start:", j1
      write(*,*) "z extent start:", k1
      stop 1
    end if
    
    if (i1>i2 .or. j1>j2 .or. k1>k2) then 
      write(*,*) "Error, one of the extents is negative."
      write(*,*) "x extent:", i1, i2
      write(*,*) "y extent:", j1, j2
      write(*,*) "z extent:", k1, k2
      stop 1
    end if
    
    if (i2>g%nx) then
      write(*,*) "Error, requested x extent out of range."
      write(*,*) "requested: ", i1, i2
      write(*,*) "input file:", 1, g%nx
      stop 1
    end if
    
    if (j2>g%ny) then
      write(*,*) "Error, requested y extent out of range."
      write(*,*) "requested: ", j1, j2
      write(*,*) "input file:", 1, g%ny
      stop 1
    end if
    
    if (k2>g%nz) then
      write(*,*) "Error, requested z extent out of range."
      write(*,*) "requested: ", k1, k2
      write(*,*) "input file:", 1, g%nz
      stop 1
    end if
    
    allocate(g%g_out)
    
    g%g_out%fname = out_name//".efh"
    g%g_out%nx = nx
    g%g_out%ny = ny
    g%g_out%nz = nz
    g%g_out%offx = i1 - 1
    g%g_out%offy = j1 - 1
    g%g_out%offz = k1 - 1
    
    g%g_out%x = g%x(i1:i2)
    g%g_out%y = g%y(j1:j2)
    g%g_out%z = g%z(k1:k2)
    
    call g%save_header
    
    file_num = 0
    do
      base_name = series_name//"-"//itoa(file_num)
      fname = base_name//".eaf"
      
      inquire(file=fname, exist=ex)
      if (.not.ex) exit
      
      write(*,'(a)',advance="no") achar(27)//"[2K"//achar(27)//"[1G"
      write(*,'(a)',advance="no") fname
      
      call g%rewind_header
      
      call g%open_frame(fname)
      
      base_name = out_name//"-"//itoa(file_num)
      call g%g_out%create_frame(base_name//".eaf")
      
      do
        call g%process_next(status)
        if (status==eof) exit
      end do
      
      call g%close_frame
      call g%g_out%close_frame
      
      file_num = file_num + 1
    end do
    
    call g%close_header
    write(*,*)
  end subroutine
end program
