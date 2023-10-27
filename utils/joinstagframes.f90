module Kinds
  use iso_fortran_env
  
  integer,parameter :: rp = real64, knd = rp
end module



module Types
  use Kinds
  
  implicit none

  type i3
    integer :: i,j,k
  end type

  type r3
    real(knd) :: x,y,z
  end type

  type irange
    type(i3) min,max
  end type

  type rrange
    type(r3) min,max
  end type

  type TSaveFlags
    logical :: Pr = .false.
    logical :: U = .false.
    logical :: V = .false.
    logical :: W = .false.
    logical :: Viscosity = .false.
    logical :: Temperature = .false.
    logical :: Moisture = .false.
    logical :: Scalar = .false.
    integer :: num_scalars = 0
  end type

  type grid
    integer :: iim = 1, jim = 1, kim = 1
    type(joined_grids), pointer :: global => null()
    integer :: Prnx = 0, Prny = 0, Prnz = 0
    integer :: Unx = 0, Uny = 0, Unz = 0
    integer :: Vnx = 0, Vny = 0, Vnz = 0
    integer :: Wnx = 0, Wny = 0, Wnz = 0
    integer :: offx = 0,offy = 0,offz = 0
    integer :: unit
    integer :: inPrnx = 0, inPrny = 0, inPrnz = 0
    integer :: inUnx = 0, inUny = 0, inUnz = 0
    integer :: inVnx = 0, inVny = 0, inVnz = 0
    integer :: inWnx = 0, inWny = 0, inWnz = 0
    
    type(rrange) :: range
    integer   :: minUi, maxUi, minUj, maxUj, minUk, maxUk
    integer   :: minVi, maxVi, minVj, maxVj, minVk, maxVk
    integer   :: minWi, maxWi, minWj, maxWj, minWk, maxWk
    integer   :: minPri, maxPri, minPrj, maxPrj, minPrk, maxPrk
    integer   :: size
    real(knd),allocatable,dimension(:) :: xU,yV,zW
    real(knd),allocatable,dimension(:) :: xPr,yPr,zPr
      !temperature and scalars use Prsize and min/maxPr_

    character(:), allocatable :: fname, header_fname, mask_fname, series_name, im_path
    type(TSaveFlags) :: save_flags
    
    real(rp) :: time
    
    real(rp), allocatable, dimension(:,:,:) :: U,V, W, Pr, Viscosity, Temperature, Moisture
    real(rp), allocatable, dimension(:,:,:,:) :: Scalar
    
    integer, allocatable, dimension(:,:,:) :: Utype, Vtype, Wtype,Prtype
  contains
    procedure :: read_header_file
    procedure :: read_data_file
    procedure :: read_mask_file
    procedure :: data_fname
    procedure :: deallocate_arrays
  end type
  
  type joined_grids
    character(:), allocatable :: series_name
    integer :: nxims = 1
    integer :: nyims = 1
    integer :: nzims = 1
    
    integer :: min_iim, min_jim, min_kim
    integer :: max_iim, max_jim, max_kim
    
    type(grid), allocatable :: gs(:,:,:)
    
    integer :: Prnx = 0, Prny = 0, Prnz = 0
    integer :: Unx = 0, Uny = 0, Unz = 0
    integer :: Vnx = 0, Vny = 0, Vnz = 0
    integer :: Wnx = 0, Wny = 0, Wnz = 0
    
    real(knd),allocatable,dimension(:) :: xU,yV,zW
    real(knd),allocatable,dimension(:) :: xPr,yPr,zPr
    
    real(rp) :: time
    
    real(rp), allocatable, dimension(:,:,:) :: U,V, W, Pr, Viscosity, Temperature, Moisture
    real(rp), allocatable, dimension(:,:,:,:) :: Scalar
    
    type(TSaveFlags) :: save_flags
  contains
    procedure :: get_shape_from_directories
    procedure :: read_and_join_data
    procedure :: move_data
    procedure :: save_header
    procedure :: save_data
  end type
  
  interface grid
    procedure grid_init
  end interface
  
  interface joined_grids
    procedure joined_grids_init
  end interface
  
  character,parameter :: lf = achar(10)

contains
  

  subroutine read_header_file(g)
    class(grid),intent(inout) :: g
    integer :: io
    integer(int32) :: one, bits

    open(newunit=g%unit,file=g%header_fname,access="stream",form="unformatted",status="old",action="read",iostat=io)

    if (io/=0) then
      g%Prnx = 0
      g%Prny = 0
      g%Prnz = 0
      return
    end if

    read(g%unit) one
    if (one/=1_int32) then
      error stop "Data files probably created with the opposite endianness."
    end if
    
    read(g%unit) bits
    if (bits/=storage_size(1._rp)) then
      write(*,*) "joinstagframes real kind bits:",storage_size(1._rp)
      write(*,*) "data file real kind bits:",bits
      error stop
    end if
    
    read(g%unit) bits
     if (bits/=storage_size(1)) then
      write(*,*) "joinstagframes integer kind bits:",storage_size(1)
      write(*,*) "data file real kind bits:",bits
      error stop
    end if
    
    read(g%unit) bits
    if (bits/=storage_size(.true.)) then
      write(*,*) "joinstagframes logical kind bits:",storage_size(.true.)
      write(*,*) "data file real kind bits:",bits
      error stop
    end if

    
    read(g%unit) g%save_flags%U, g%save_flags%V, g%save_flags%W
    read(g%unit) g%save_flags%Pr, g%save_flags%Viscosity
    read(g%unit) g%save_flags%Temperature, g%save_flags%Moisture, &
                  g%save_flags%Scalar, g%save_flags%num_scalars

    read(g%unit) g%minUi, g%maxUi, g%minUj, g%maxUj, g%minUk, g%maxUk
    read(g%unit) g%minVi, g%maxVi, g%minVj, g%maxVj, g%minVk, g%maxVk
    read(g%unit) g%minWi, g%maxWi, g%minWj, g%maxWj, g%minWk, g%maxWk
    read(g%unit) g%minPri, g%maxPri, g%minPrj, g%maxPrj, g%minPrk, g%maxPrk

    g%Prnx = g%maxPri - g%minPri + 1
    g%Prny = g%maxPrj - g%minPrj + 1
    g%Prnz = g%maxPrk - g%minPrk + 1

    g%Unx = g%maxUi - g%minUi + 1
    g%Uny = g%maxUj - g%minUj + 1
    g%Unz = g%maxUk - g%minUk + 1
    g%Vnx = g%maxVi - g%minVi + 1
    g%Vny = g%maxVj - g%minVj + 1
    g%Vnz = g%maxVk - g%minVk + 1
    g%Wnx = g%maxWi - g%minWi + 1
    g%Wny = g%maxWj - g%minWj + 1
    g%Wnz = g%maxWk - g%minWk + 1
    
    g%inPrnx = g%Prnx - 2
    g%inPrny = g%Prny - 2
    g%inPrnz = g%Prnz - 2
    g%inUnx =  g%inPrnx
    g%inUny =  g%inPrny
    g%inUnz =  g%inPrnz
    g%inVnx =  g%inPrnx
    g%inVny =  g%inPrny
    g%inVnz =  g%inPrnz
    g%inWnx =  g%inPrnx
    g%inWny =  g%inPrny
    g%inWnz =  g%inPrnz
    
    g%size = g%Prnx * g%Prny * g%Prnz

    allocate(g%xU(0:g%Unx-1))
    allocate(g%yV(0:g%Vny-1))
    allocate(g%zW(0:g%Wnz-1))
    allocate(g%xPr(1:g%Prnx))
    allocate(g%yPr(1:g%Prny))
    allocate(g%zPr(1:g%Prnz))

    read(g%unit) g%xU,g%yV,g%zW,g%xPr,g%yPr,g%zPr

    close(g%unit)
    
    if (g%save_flags%U) allocate(g%U(0:g%Unx-1, 1:g%Uny, 1:g%Unz))
    if (g%save_flags%V) allocate(g%V(1:g%Vnx, 0:g%Vny-1, 1:g%Vnz))
    if (g%save_flags%W) allocate(g%W(1:g%Wnx, 1:g%Wny, 0:g%Wnz-1))
    if (g%save_flags%Pr) allocate(g%Pr(g%Prnx, g%Prny, g%Prnz))
    if (g%save_flags%Viscosity) allocate(g%Viscosity(g%Prnx, g%Prny, g%Prnz))
    if (g%save_flags%Temperature) allocate(g%Temperature(g%Prnx, g%Prny, g%Prnz))
    if (g%save_flags%Moisture) allocate(g%Moisture(g%Prnx, g%Prny, g%Prnz))
    if (g%save_flags%Scalar) allocate(g%Scalar(g%Prnx, g%Prny, g%Prnz, g%save_flags%num_scalars))
    
  end subroutine
  
  subroutine read_mask_file(g)
    class(grid),intent(inout) :: g
    integer :: io

    allocate(g%Utype(g%Unx, g%Uny, g%Unz))
    allocate(g%Vtype(g%Vnx, g%Vny, g%Vnz))
    allocate(g%Wtype(g%Wnx, g%Wny, g%Wnz))
    allocate(g%Prtype(g%Prnx, g%Prny, g%Prnz))
    
    open(newunit=g%unit,file=g%mask_fname,access="stream",form="unformatted",status="old",action="read",iostat=io)
    read(g%unit) g%Utype
    read(g%unit) g%Vtype
    read(g%unit) g%Wtype
    read(g%unit) g%Prtype
    close(g%unit)
  end subroutine
  

  
  subroutine read_data_file(g, file_n, stat)
    class(grid),intent(inout) :: g
    integer, intent(in) :: file_n
    integer, intent(out) :: stat
    integer :: io
 
    open(newunit=g%unit,file=g%data_fname(file_n),access="stream",form="unformatted",status="old",action="read",iostat=io)
    if (io/=0) then
      stat = 1
      return
    end if
    read(g%unit) g%time
    if (g%save_flags%U) read(g%unit) g%U
    if (g%save_flags%V) read(g%unit) g%V
    if (g%save_flags%W) read(g%unit) g%W
    if (g%save_flags%Pr) read(g%unit) g%Pr
    if (g%save_flags%Viscosity) read(g%unit) g%Viscosity
    if (g%save_flags%Temperature) read(g%unit) g%Temperature
    if (g%save_flags%Moisture) read(g%unit) g%Moisture
    if (g%save_flags%Scalar) read(g%unit) g%Scalar
    close(g%unit)

    stat = 0
  end subroutine
  
  elemental subroutine deallocate_arrays(g)
    class(grid),intent(inout) :: g
    if (g%save_flags%U) deallocate(g%U)
    if (g%save_flags%V) deallocate(g%V)
    if (g%save_flags%W) deallocate(g%W)
    if (g%save_flags%Pr) deallocate(g%Pr)
    if (g%save_flags%Viscosity) deallocate(g%Viscosity)
    if (g%save_flags%Temperature) deallocate(g%Temperature)
    if (g%save_flags%Moisture) deallocate(g%Moisture)
    if (g%save_flags%Scalar) deallocate(g%Scalar)
    if (allocated(g%Utype)) deallocate(g%Utype)
    if (allocated(g%Vtype)) deallocate(g%Vtype)
    if (allocated(g%Wtype)) deallocate(g%Wtype)
    if (allocated(g%Prtype)) deallocate(g%Prtype)
  end subroutine
  
  pure function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function
  
  pure function data_fname(g, file_n)
    character(:), allocatable :: data_fname
    class(grid), intent(in) :: g
    integer, intent(in) :: file_n
  
    data_fname = g%im_path//"stagframe-"//g%series_name//"-"//itoa(file_n)//".unf"
  end function
  
  
  
  subroutine get_shape_from_directories(gs)
    use iso_c_binding
    
    class(joined_grids), intent(inout) :: gs
    character(kind=c_char, len=:), allocatable :: dirname
    integer :: i
    INTERFACE
      SUBROUTINE file_info(filename,mode,exist,time) BIND(C,name="file_info")
        USE iso_c_binding
        CHARACTER(kind=C_CHAR),INTENT(in) :: filename(*)
        INTEGER(C_INT),INTENT(out) :: mode, exist, time
      END SUBROUTINE
    END INTERFACE
    integer(c_int) :: mode, ex, time

    
    i = 1
    do
      dirname =  "im-"//itoa(i)//"-1-1/" // c_null_char
      call file_info(dirname,mode,ex,time)      
      if (ex==0) exit
      i = i + 1
    end do
    gs%nxims = i-1
    
    i = 1
    do
      dirname =  "im-1-"//itoa(i)//"-1/" // c_null_char
      call file_info(dirname,mode,ex,time)
      if (ex==0) exit
      i = i + 1
    end do
    gs%nyims = i-1
    
    i = 1
    do
      dirname =  "im-1-1-"//itoa(i)//"/" // c_null_char
      call file_info(dirname,mode,ex,time)
      if (ex==0) exit
      i = i + 1
    end do
    gs%nzims = i-1
  end subroutine
  
  
  function grid_init(gs, iim, jim, kim) result(g)
    type(grid) :: g
    class(joined_grids), intent(in), target :: gs
    integer, intent(in) :: iim, jim, kim
    
    g%global => gs
    g%series_name = gs%series_name
    g%iim = iim
    g%jim = jim
    g%kim = kim
    g%im_path = "im-"//itoa(g%iim)//"-"//itoa(g%jim)//"-"//itoa(g%kim)//"/"
    g%header_fname = g%im_path//"stagframe-"//g%series_name//"-"//"header.unf"
    g%mask_fname = g%im_path//"stagframe-"//g%series_name//"-"//"mask.unf"
  end function

  

  function joined_grids_init(series_name) result(gs)
    type(joined_grids) :: gs
    character(*), intent(in) :: series_name
    integer :: iim, jim, kim
    
    gs%series_name = series_name
      
    call gs%get_shape_from_directories

    allocate(gs%gs(gs%nxims, gs%nyims, gs%nzims))
    
    do kim = 1, gs%nzims
      do jim = 1, gs%nyims
        do iim = 1, gs%nxims
          gs%gs(iim, jim, kim) = grid(gs, iim, jim, kim)
          call gs%gs(iim, jim, kim)%read_header_file
        end do
      end do
    end do
    
    
    gs%min_iim = gs%nxims+1
    do iim = 1, gs%nxims
      if (any(gs%gs(iim,:,:)%size > 0)) then
        gs%min_iim = iim
        exit
      end if
    end do
    
    if (gs%min_iim > gs%nxims) then
      write(*,*) "Error, no header files found or all empty."
      stop 2
    end if
    
    gs%min_jim = gs%nyims+1
    do jim = 1, gs%nyims
      if (any(gs%gs(:,jim,:)%size > 0)) then
        gs%min_jim = jim
        exit
      end if
    end do
    gs%min_kim = gs%nzims+1
    do kim = 1, gs%nzims
      if (any(gs%gs(:,:,kim)%size > 0)) then
        gs%min_kim = kim
        exit
      end if
    end do
    
    gs%max_iim = 0
    do iim = gs%nxims, 1, -1
      if (any(gs%gs(iim,:,:)%size > 0)) then
        gs%max_iim = iim
        exit
      end if
    end do
    gs%max_jim = 0
    do jim = gs%nyims, 1, -1
      if (any(gs%gs(:,jim,:)%size > 0)) then
        gs%max_jim = jim
        exit
      end if
    end do
    gs%max_kim = 0
    do kim = gs%nzims, 1, -1
      if (any(gs%gs(:,:, kim)%size > 0)) then
        gs%max_kim = kim
        exit
      end if
    end do

    gs%gs(gs%min_iim, :, :)%offx = 0
    do iim = gs%min_iim+1, gs%max_iim
      gs%gs(iim, :, :)%offx = gs%gs(iim-1, :, :)%offx + gs%gs(iim-1, gs%min_jim, gs%min_kim)%inPrnx
    end do
    gs%Prnx = gs%gs(gs%max_iim,gs%min_jim,gs%min_kim)%offx + gs%gs(gs%max_iim,gs%min_jim,gs%min_kim)%inPrnx
    gs%Unx = gs%gs(gs%max_iim,gs%min_jim,gs%min_kim)%offx + gs%gs(gs%max_iim,gs%min_jim,gs%min_kim)%inUnx
    gs%Vnx = gs%gs(gs%max_iim,gs%min_jim,gs%min_kim)%offx + gs%gs(gs%max_iim,gs%min_jim,gs%min_kim)%inVnx
    gs%Wnx = gs%gs(gs%max_iim,gs%min_jim,gs%min_kim)%offx + gs%gs(gs%max_iim,gs%min_jim,gs%min_kim)%inWnx
    
    gs%gs(:, gs%min_jim, :)%offy = 0
    do jim = gs%min_jim+1, gs%max_jim
      gs%gs(:, jim, :)%offy = gs%gs(:, jim-1, :)%offy + gs%gs(gs%min_iim, jim-1, gs%min_kim)%inPrny
    end do
    gs%Prny = gs%gs(gs%min_iim,gs%max_jim,gs%min_kim)%offy + gs%gs(gs%min_iim,gs%max_jim,gs%min_kim)%inPrny
    gs%Uny = gs%gs(gs%min_iim,gs%max_jim,gs%min_kim)%offy + gs%gs(gs%min_iim,gs%max_jim,gs%min_kim)%inUny
    gs%Vny = gs%gs(gs%min_jim,gs%max_jim,gs%min_kim)%offy + gs%gs(gs%min_iim,gs%max_jim,gs%min_kim)%inVny
    gs%Wny = gs%gs(gs%min_iim,gs%max_jim,gs%min_kim)%offy + gs%gs(gs%min_iim,gs%max_jim,gs%min_kim)%inWny
    
    gs%gs(:, :, gs%min_kim)%offz = 0
    do kim = gs%min_kim+1, gs%max_kim
      gs%gs(:, :, kim)%offz = gs%gs(:, :, kim-1)%offz + gs%gs(gs%min_iim, gs%min_jim, kim-1)%inPrnz
    end do
    gs%Prnz = gs%gs(gs%min_iim,gs%min_jim,gs%max_kim)%offz + gs%gs(gs%min_iim,gs%min_jim,gs%max_kim)%inPrnz
    gs%Unz = gs%gs(gs%min_iim,gs%min_jim,gs%max_kim)%offz + gs%gs(gs%min_iim,gs%min_jim,gs%max_kim)%inUnz
    gs%Vnz = gs%gs(gs%min_jim,gs%min_jim,gs%max_kim)%offz + gs%gs(gs%min_iim,gs%min_jim,gs%max_kim)%inVnz
    gs%Wnz = gs%gs(gs%min_iim,gs%min_jim,gs%max_kim)%offz + gs%gs(gs%min_iim,gs%min_jim,gs%max_kim)%inWnz
    
    write(*,*) "Prn:", gs%Prnx, gs%Prny, gs%Prnz
    write(*,*) "Un:", gs%Unx, gs%Uny, gs%Unz
    write(*,*) "Vn:", gs%Vnx, gs%Vny, gs%Vnz
    write(*,*) "Wn:", gs%Wnx, gs%Wny, gs%Wnz
    
    gs%save_flags = gs%gs(gs%min_iim, gs%min_jim, gs%min_kim)%save_flags
    
    if (gs%save_flags%U) allocate(gs%U(1:gs%Unx, 1:gs%Uny, 1:gs%Unz))
    if (gs%save_flags%V) allocate(gs%V(1:gs%Vnx, 1:gs%Vny, 1:gs%Vnz))
    if (gs%save_flags%W) allocate(gs%W(1:gs%Wnx, 1:gs%Wny, 1:gs%Wnz))
    if (gs%save_flags%Pr) allocate(gs%Pr(1:gs%Prnx, 1:gs%Prny, 1:gs%Prnz))
    if (gs%save_flags%Viscosity) allocate(gs%Viscosity(1:gs%Prnx, 1:gs%Prny, 1:gs%Prnz))
    if (gs%save_flags%Temperature) allocate(gs%Temperature(1:gs%Prnx, 1:gs%Prny, 1:gs%Prnz))
    if (gs%save_flags%Moisture) allocate(gs%Moisture(1:gs%Prnx, 1:gs%Prny, 1:gs%Prnz))
    if (gs%save_flags%Scalar) allocate(gs%Scalar(1:gs%Prnx, 1:gs%Prny, 1:gs%Prnz, gs%save_flags%num_scalars))
    
    allocate(gs%xU(1:gs%Unx))
    allocate(gs%yV(1:gs%Vny))
    allocate(gs%zW(1:gs%Wnz))
    allocate(gs%xPr(1:gs%Prnx))
    allocate(gs%yPr(1:gs%Prny))
    allocate(gs%zPr(1:gs%Prnz))
    
    do iim = gs%min_iim, gs%max_iim
      associate (g=>gs%gs(iim, gs%min_jim, gs%min_kim))
        gs%xU(g%offx + 1 : g%offx + g%inPrnx) = g%xU(2 : g%inPrnx+1)
      end associate
    end do
    do jim = gs%min_jim, gs%max_jim
      associate (g=>gs%gs(gs%min_iim, jim, gs%min_kim))
        gs%yV(g%offy + 1 : g%offy + g%inPrny) = g%yV(2 : g%inPrny+1)
      end associate
    end do
    do kim = gs%min_kim, gs%max_kim
      associate (g=>gs%gs(gs%min_iim, gs%min_jim, kim))
        gs%zW(g%offz + 1 : g%offz + g%inPrnz) = g%zW(2 : g%inPrnz+1)
      end associate
    end do
    do iim = gs%min_iim, gs%max_iim
      associate (g=>gs%gs(iim, gs%min_jim, gs%min_kim))
        gs%xPr(g%offx + 1 : g%offx + g%inPrnx) = g%xPr(2 : g%inPrnx+1)
      end associate
    end do
    do jim = gs%min_jim, gs%max_jim
      associate (g=>gs%gs(gs%min_iim, jim, gs%min_kim))
        gs%yPr(g%offy + 1 : g%offy + g%inPrny) = g%yPr(2 : g%inPrny+1)
      end associate
    end do
    do kim = gs%min_kim, gs%max_kim
      associate (g=>gs%gs(gs%min_iim, gs%min_jim, kim))
        gs%zPr(g%offz + 1 : g%offz + g%inPrnz) = g%zPr(2 : g%inPrnz+1)
      end associate
    end do   
  end function joined_grids_init
  
  subroutine move_data(gs, g)
    class(joined_grids), intent(inout) :: gs
    class(grid), intent(in) :: g   
    integer :: nx, ny, nz, ox, oy, oz
     
    ox = g%offx
    oy = g%offy
    oz = g%offz

    nx = g%inUnx
    ny = g%inUny
    nz = g%inUnz
    if (gs%save_flags%U) gs%U(1+ox:nx+ox, 1+oy:ny+oy, 1+oz:nz+oz) = g%U(2:nx+1, 2:ny+1, 2:nz+1)
    nx = g%inVnx
    ny = g%inVny
    nz = g%inVnz
    if (gs%save_flags%V) gs%V(1+ox:nx+ox, 1+oy:ny+oy, 1+oz:nz+oz) = g%V(2:nx+1, 2:ny+1, 2:nz+1)
    nx = g%inWnx
    ny = g%inWny
    nz = g%inWnz
    if (gs%save_flags%W) gs%W(1+ox:nx+ox, 1+oy:ny+oy, 1+oz:nz+oz) = g%W(2:nx+1, 2:ny+1, 2:nz+1)
    
    nx = g%inPrnx
    ny = g%inPrny
    nz = g%inPrnz
    if (gs%save_flags%Pr) gs%Pr(1+ox:nx+ox, 1+oy:ny+oy, 1+oz:nz+oz) = g%Pr(2:nx+1, 2:ny+1, 2:nz+1)
    if (gs%save_flags%Viscosity) gs%Viscosity(1+ox:nx+ox, 1+oy:ny+oy, 1+oz:nz+oz) = g%Viscosity(2:nx+1, 2:ny+1, 2:nz+1)
    if (gs%save_flags%Temperature) gs%Temperature(1+ox:nx+ox, 1+oy:ny+oy, 1+oz:nz+oz) = g%Temperature(2:nx+1, 2:ny+1, 2:nz+1)
    if (gs%save_flags%Moisture) gs%Moisture(1+ox:nx+ox, 1+oy:ny+oy, 1+oz:nz+oz) = g%Moisture(2:nx+1, 2:ny+1, 2:nz+1)
    if (gs%save_flags%Scalar) gs%Scalar(1+ox:nx+ox, 1+oy:ny+oy, 1+oz:nz+oz,:) = g%Scalar(2:nx+1, 2:ny+1, 2:nz+1,:)
  end subroutine
  
  
  
  subroutine read_and_join_data(gs, first, last)
    class(joined_grids), intent(inout) :: gs
    integer, intent(in) :: first, last
    integer :: stat
    integer :: iframe, iim, jim, kim
    
    call gs%save_header
    
    outer: do iframe = first, last
      write(*,*) "frame:", iframe
      do kim = gs%min_kim, gs%max_kim
        do jim = gs%min_jim, gs%max_jim
          do iim = gs%min_iim, gs%max_iim
            associate(g=>gs%gs(iim, jim, kim))
              call g%read_data_file(iframe, stat)
              if (stat/=0) exit outer
              call gs%move_data(g)
              gs%time = g%time
            end associate
            
          end do
        end do
      end do

      call gs%save_data(iframe)
    end do outer
  
    call gs%gs%deallocate_arrays

  end subroutine
  
  
  subroutine save_header(gs)
    class(joined_grids), intent(in) ::gs
    character(:), allocatable :: file_name
    integer :: unit
    
    file_name = "stagframe-"//gs%series_name//"-header.unf"

    open(newunit=unit, file=file_name, access='stream', form='unformatted', action='write', status='replace')
    write(unit) 1_int32 !endianess can be infered from this
    write(unit) int(storage_size(1._knd),int32)  !save number of bits of the used real kind
    write(unit) int(storage_size(1),int32)  !save number of bits of the used (default) integer kind
    write(unit) int(storage_size(.true.),int32)  !save number of bits of the used (default) logical kind

    write(unit) gs%save_flags%U, gs%save_flags%V, gs%save_flags%W
    write(unit) gs%save_flags%Pr, gs%save_flags%Viscosity
    write(unit) gs%save_flags%Temperature, gs%save_flags%Moisture, &
                  gs%save_flags%Scalar, gs%save_flags%num_scalars

    write(unit) 1, gs%Unx, 1, gs%Uny, 1, gs%Unz
    write(unit) 1, gs%Vnx, 1, gs%Vny, 1, gs%Vnz
    write(unit) 1, gs%Wnx, 1, gs%Wny, 1, gs%Wnz
    write(unit) 1, gs%Prnx, 1, gs%Prny, 1, gs%Prnz
    write(unit) gs%xU,gs%yV,gs%zW,gs%xPr,gs%yPr,gs%zPr

    close(unit)
  end subroutine
  
  subroutine save_data(gs, iframe)
    class(joined_grids), intent(in) ::gs
    integer, intent(in) :: iframe
    character(:), allocatable :: file_name
    integer :: unit
    
    file_name = "stagframe-"//gs%series_name//"-"//itoa(iframe)//".unf"
    
    open(newunit=unit, file=file_name, access="stream", form="unformatted", status="replace", action="write")
    write(unit) gs%time
    if (gs%save_flags%U) write(unit) gs%U
    if (gs%save_flags%V) write(unit) gs%V
    if (gs%save_flags%W) write(unit) gs%W
    if (gs%save_flags%Pr) write(unit) gs%Pr
    if (gs%save_flags%Viscosity) write(unit) gs%Viscosity
    if (gs%save_flags%Temperature) write(unit) gs%Temperature
    if (gs%save_flags%Moisture) write(unit) gs%Moisture
    if (gs%save_flags%Scalar) write(unit) gs%Scalar
    close(unit)
  end subroutine
  
end module Types

program joinstagframes
  use Types

  character(:), allocatable :: domain
  character(80) :: arg
  integer :: arg_len, first = 0, last = huge(1)
  
  type(joined_grids), target :: global
  
  if (command_argument_count()<1) then
    write(*,*) "Usage: joinstagframes label [first, last]"
    stop 1
  end if
 
  call get_command_argument(1, length=arg_len)
  allocate(character(arg_len) :: domain)
  
  call get_command_argument(1, value=domain)
  
  if (command_argument_count()==3) then
    call get_command_argument(2, value=arg)
    read(arg,*) first
    call get_command_argument(3, value=arg)
    read(arg,*) last
  end if
  
  global = joined_grids(domain)

  call global%read_and_join_data(first, last)
end program
