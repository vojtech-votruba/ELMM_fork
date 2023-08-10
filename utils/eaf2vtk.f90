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
  use Endianness
  
  implicit none

  type grid
    integer nx,ny,nz
    integer hu, fu, vtku
    integer header_pos
    real(rp),allocatable :: x(:),y(:),z(:)
    character(1024) :: fname
    real(rp), pointer, contiguous :: vec(:,:,:,:)
    real(rp), pointer, contiguous :: sc(:,:,:) => null() ! to save space in memory
  contains
    procedure :: process_next
    procedure :: read_var_from_header
    procedure :: read_var_from_eaf
    procedure :: save_header
    procedure :: save_var
    procedure :: open_frame
    procedure :: rewind_header
    procedure :: close_header, close_frame, close_vtk
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
      stop "Running eaf2vtk on a computer with different endianness not yet supported."
    end if
    
    read(g%hu) itmp
    if (itmp /= storage_size(1_rp)/CHARACTER_STORAGE_SIZE) then
      write(*,*) "eaf2vtk compiled for different real precision.Storage sizes:"
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

  subroutine save_header(g, fname)
    class(grid),intent(inout) :: g
    character(*), intent(in) :: fname
    character(70) :: str

    open(newunit=g%vtku,file = fname, &
      access='stream',status='replace',form="unformatted",action="write")

    write(g%vtku) "# vtk DataFile Version 2.0", lf
    write(g%vtku) "CLMM output file", lf
    write(g%vtku) "BINARY", lf
    write(g%vtku) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),*) g%nx, g%ny, g%nz
    write(g%vtku) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') g%nx,"float"
    write(g%vtku) str, lf
    write(g%vtku) BigEnd(g%x), lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') g%ny,"float"
    write(g%vtku) str, lf
    write(g%vtku) BigEnd(g%y), lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') g%nz,"float"
    write(g%vtku) str, lf
    write(g%vtku) BigEnd(g%z), lf
    str="POINT_DATA"
    write(str(12:),*) g%nx*g%ny*g%nz
    write(g%vtku) str, lf
  end subroutine

  subroutine save_var(g, title, status)
    class(grid),intent(inout) :: g
    character(*), intent(in) :: title
    integer, intent(in) :: status
 
    if (status == SCALAR) then
      g%sc = BigEnd(g%sc)
      write(g%vtku) "SCALARS "//title//" float", lf
      write(g%vtku) "LOOKUP_TABLE default", lf
      write(g%vtku) g%sc, lf
    else
      g%vec = BigEnd(g%vec)
      write(g%vtku) "VECTORS "//title//" float", lf
      write(g%vtku) g%vec, lf
    end if
  end subroutine

  
  subroutine process_next(g,status)
    class(grid), intent(inout) ::g
    integer, intent(inout) :: status
    character(:), allocatable :: title
    
    call g%read_var_from_header(title, status)
   
    if (status /= EOF) then
      call g%read_var_from_eaf(status)
      call g%save_var(title, status)
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
  
  subroutine close_vtk(g)
    class(grid), intent(inout) :: g
    close(g%vtku)
  end subroutine
end module Types


module file_names
  implicit none
  
  interface
    function open_dir(dir_name) result(res) bind(C, name="open_dir")
      use iso_c_binding
      type(c_ptr) :: res
      character(kind=c_char, len=1), intent(in) :: dir_name(*)
    end function
    
    subroutine close_dir(dir) bind(C, name="close_dir")
      use iso_c_binding
      type(c_ptr), value :: dir
    end subroutine
    
    subroutine next_file(dir, file_name, name_len) bind(C, name="next_file")
      use iso_c_binding
      type(c_ptr), value :: dir
      character(kind=c_char, len=1), intent(out) :: file_name(256)
      integer(c_int), intent(out) :: name_len
    end subroutine
  end interface
end module


program eaf2vtk
  use iso_c_binding
  use Types
  use file_names
  
  implicit none
  
  character(:), allocatable :: file_name
  
  integer :: arg_len, status
  
  logical :: ex
  integer :: file_num
  
  type(c_ptr) :: dir_ptr = c_null_ptr
  integer(c_int) :: name_len
  
  call GetEndianness()
  
  if (command_argument_count()<1) then
    write(*,*) "Usage:"
    write(*,*) "eaf2vtk file_name.vtk"    
    write(*,*)
    write(*,*) "   for one file"
    write(*,*)
    write(*,*)
    write(*,*) "eaf2vtk series_name"    
    write(*,*)
    write(*,*) "   for the complete series. File names 'frame-series_name-number.eaf' are assumed."
    write(*,*)
    write(*,*)
    write(*,*) "eaf2vtk all"    
    write(*,*)
    write(*,*) "   for all series in the directory."
    stop 1
  end if
  
  call get_command_argument(1, length=arg_len)
  
  allocate(character(arg_len) :: file_name)
  call get_command_argument(1, value=file_name)
  
  if (file_name=='all') then
  
    deallocate(file_name)
    allocate(character(256) :: file_name)
    file_name(:) = ""
    
    dir_ptr = open_dir("."//c_null_char)
    
    do
      call next_file(dir_ptr, file_name, name_len)
      if (name_len > 4) then
        if (file_name(name_len-3:name_len)==".efh") then
          write(*,*) trim(file_name),":"
          call process_series(file_name(:name_len-4))
        end if
      else if (name_len<=0) then
        exit
      end if
    end do
    
    call close_dir(dir_ptr);
    
  else if (file_name(max(arg_len-3,1):arg_len)==".eaf") then
    call process_file(file_name)
  else
    call process_series("frame-"//file_name)    
  end if
  
  deallocate(file_name)
  
contains

  subroutine process_series(series_name)
    character(*), intent(in) :: series_name
    character(:), allocatable :: fname, base_name
    type(grid), target :: g
  
    
    g = grid(series_name//".efh")
    
    file_num = 0
    do
      base_name = series_name//"-"//itoa(file_num)
      fname = base_name//".eaf"
      
      inquire(file=fname, exist=ex)
      if (.not.ex) exit
      
      write(*,'(a)',advance="no") achar(27)//"[2K"//achar(27)//"[1G"
      write(*,'(a)',advance="no") fname
      
      call g%rewind_header
      
      call g%save_header(base_name//".vtk")
      call g%open_frame(fname)
      do
        call g%process_next(status)
        if (status==eof) exit
      end do
      
      call g%close_frame
      call g%close_vtk
      
      file_num = file_num + 1
    end do
    
    call g%close_header
    write(*,*)
  end subroutine
  
  subroutine process_file(file_name)
    character(*), intent(in) :: file_name
    character(:), allocatable :: base_name
    type(grid), target :: g
  
    base_name = file_name(:arg_len-4)
    g = grid(base_name(:scan(base_name, "-", back=.true.)-1)//".efh")
    call g%save_header(base_name//".vtk")
    call g%open_frame(file_name)
    do
      call g%process_next(status)
      if (status==eof) exit
    end do
    
    call g%close_header
    call g%close_frame
    call g%close_vtk
  end subroutine
  
  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function
end program
