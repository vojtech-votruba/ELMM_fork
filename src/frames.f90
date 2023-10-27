module StaggeredFrames
  use iso_c_binding, only: c_ptr
  use Kinds
  use Pthreads
  use Frames_common

  implicit none

  private
  public i3, r3, irange, rrange, TFrameTimes, TSaveFlags, &
         TStaggeredFrameDomain, SaveStaggeredFrames, FinalizeStaggeredFrames, AddDomain


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

  type, extends(TFrameBase) :: TStaggeredFrameDomain
    !domain for saving full timed data
    !fluxes should be computed by postprocessing, preferably using the same numerical methods
    !as in the primary solver
    private

    type(rrange) :: range
    integer   :: minUi, maxUi, minUj, maxUj, minUk, maxUk
    integer   :: minVi, maxVi, minVj, maxVj, minVk, maxVk
    integer   :: minWi, maxWi, minWj, maxWj, minWk, maxWk
    integer   :: sizeU, sizeV, sizeW
    real(knd),allocatable,dimension(:) :: xU,yV,zW
      !temperature and scalars use Prsize and min/maxPr_

    type(TSaveFlags) :: save_flags

    integer :: buffer_size = 0
    real(knd),allocatable :: buffer(:)

    character(4)  :: suffix = ".unf"
  contains
     procedure :: Fill => TStaggeredFrameDomain_Fill
     procedure :: DoSave => TStaggeredFrameDomain_DoSave
     procedure :: SetRest => TStaggeredFrameDomain_SetRest
     procedure :: SaveHeader => TStaggeredFrameDomain_SaveHeader
     procedure :: SaveMask => TStaggeredFrameDomain_SaveMask
  endtype TStaggeredFrameDomain

  interface TStaggeredFrameDomain
    module procedure TStaggeredFrameDomain_Init
  end interface
  
  interface AddDomain
    module procedure AddStaggeredFrameDomain
  end interface

  interface assignment(=)
    module procedure assign3to1
  end interface

  type(TStaggeredFrameDomain),allocatable :: StaggeredFrameDomains(:)

  integer :: stagframe_unit = 999000

contains

  subroutine assign3to1(A1,A3)
    real(knd),intent(out) :: A1(1:)
    real(knd),intent(in)  :: A3(1:,1:,1:)
    integer :: j,k
    integer :: s1,s2,s3 !sizes
    integer :: off2,off3 !offsets

    !if s1*s2*s3 /= size(A1) let it fail by the Fortran runtime library

    s1 = size(A3,1)
    s2 = size(A3,2)
    s3 = size(A3,3)

    off3 = 0
    do k=1,size(A3,3)
      off2 = 0
      do j=1,size(A3,2)
        A1(off2+off3+1:off2+off3+s1) = A3(:,j,k)
        off2 = off2 + s1
      end do
      off3 = off3 + s1*s2
    end do
  end subroutine assign3to1
  
  
  
  
  
  

  subroutine add_element_sfd(a,e)
    type(TStaggeredFrameDomain),allocatable,intent(inout) :: a(:)
    type(TStaggeredFrameDomain),intent(in) :: e
    type(TStaggeredFrameDomain),allocatable :: tmp(:)

    if (.not.allocated(a)) then
      a = [e]
    else
      call move_alloc(a,tmp)
      allocate(a(size(tmp)+1))
      a(1:size(tmp)) = tmp
      a(size(tmp)+1) = e
    end if
  end subroutine
  
  
  subroutine AddStaggeredFrameDomain(D)
  type(TStaggeredFrameDomain),intent(in) :: D

    call add_element_sfd(StaggeredFrameDomains, D)

  end subroutine
  
  subroutine SaveStaggeredFrames(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    real(knd),intent(in) :: time
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),contiguous,intent(in) :: Pr(-1:,-1:,-1:), Viscosity(-1:,-1:,-1:), &
                                       Temperature(-2:,-2:,-2:), &
                                       Moisture(-2:,-2:,-2:), Scalar(-2:,-2:,-2:,1:)
    integer :: i
    
    if (allocated(StaggeredFrameDomains)) then
      do i=1,size(StaggeredFrameDomains)
        call StaggeredFrameDomains(i)%Save(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
      end do
    end if
  end subroutine
  
  subroutine FinalizeStaggeredFrames
    integer :: i
    if (allocated(StaggeredFrameDomains)) then
      do i=1,size(StaggeredFrameDomains)
        call StaggeredFrameDomains(i)%Finalize
      end do
      deallocate(StaggeredFrameDomains)
    end if
  end subroutine
  
  
  
  
  
  
  


  function TStaggeredFrameDomain_Init(label,range,time_params,save_flags) result(D)
    use Parameters, only: output_dir
    type(TStaggeredFrameDomain) :: D
    character(*) :: label
    type(rrange),intent(in) :: range
    type(TFrameTimes),intent(in) :: time_params
    type(TSaveFlags),intent(in) :: save_flags

    D%base_name = trim(output_dir)//"stagframe-"//label

    D%frame_times = time_params

    D%range = range

    D%save_flags = save_flags

    D%unit = stagframe_unit + 1
    stagframe_unit = D%unit

  end function


  subroutine TStaggeredFrameDomain_SetRest(D, num_scalars)

    use Boundaries, only: GridCoords
    use Parameters, only: xU, yV, zW, xPr, yPr, zPr, Prtype, Utype, Vtype, Wtype

    class(TStaggeredFrameDomain),intent(inout) :: D
    integer,intent(in) :: num_scalars

    call GridCoords(D%minPri, D%minPrj, D%minPrk, D%range%min%x, D%range%min%y, D%range%min%z)

    call GridCoords(D%maxPri, D%maxPrj, D%maxPrk, D%range%max%x, D%range%max%y, D%range%max%z)
    
    
    !FIXME?: this assumes the user wants the domain specified by D%range + side buffers for flux evaluations. Make this behaviour optional?
    D%minPri = D%minPri - 1
    D%minPrj = D%minPrj - 1
    D%minPrk = D%minPrk - 1
    D%maxPri = D%maxPri + 1
    D%maxPrj = D%maxPrj + 1
    D%maxPrk = D%maxPrk + 1

    D%minUi = D%minPri - 1
    D%minUj = D%minPrj
    D%minUk = D%minPrk

    D%minVi = D%minPri
    D%minVj = D%minPrj - 1
    D%minVk = D%minPrk

    D%minWi = D%minPri
    D%minWj = D%minPrj
    D%minWk = D%minPrk - 1

    D%maxUi = D%maxPri
    D%maxUj = D%maxPrj
    D%maxUk = D%maxPrk

    D%maxVi = D%maxPri
    D%maxVj = D%maxPrj
    D%maxVk = D%maxPrk

    D%maxWi = D%maxPri
    D%maxWj = D%maxPrj
    D%maxWk = D%maxPrk

    D%sizeU = max( (D%maxUi - D%minUi + 1) * (D%maxUj - D%minUj + 1) * (D%maxUk - D%minUk + 1) , 0 )
    D%sizeV = max( (D%maxVi - D%minVi + 1) * (D%maxVj - D%minVj + 1) * (D%maxVk - D%minVk + 1) , 0 )
    D%sizeW = max( (D%maxWi - D%minWi + 1) * (D%maxWj - D%minWj + 1) * (D%maxWk - D%minWk + 1) , 0 )
    D%sizePr = max( (D%maxPri - D%minPri + 1) * (D%maxPrj - D%minPrj + 1) * (D%maxPrk - D%minPrk + 1) , 0 )

    allocate(D%xU(D%maxUi-D%minUi+1))  !not needed if Fortran 2003 reallocation enabled, but problems with zero sized arrays in gfortran 4.7.1
    D%xU(:) = xU(D%minUi:D%maxUi)
    allocate(D%yV(D%maxVj-D%minVj+1))
    D%yV(:) = yV(D%minVj:D%maxVj)
    allocate(D%zW(D%maxWk-D%minWk+1))
    D%zW(:) = zW(D%minWk:D%maxWk)
    allocate(D%xPr(D%maxPri-D%minPri+1))
    D%xPr(:) = xPr(D%minPri:D%maxPri)
    allocate(D%yPr(D%maxPrj-D%minPrj+1))
    D%yPr(:) = yPr(D%minPrj:D%maxPrj)
    allocate(D%zPr(D%maxPrk-D%minPrk+1))
    D%zPr(:) = zPr(D%minPrk:D%maxPrk)

    D%buffer_size = 1 !first position is a time-stamp

    if (D%save_flags%U) D%buffer_size = D%buffer_size + D%sizeU
    if (D%save_flags%V) D%buffer_size = D%buffer_size + D%sizeV
    if (D%save_flags%W) D%buffer_size = D%buffer_size + D%sizeW
    if (D%save_flags%Pr) D%buffer_size = D%buffer_size + D%sizePr
    if (D%save_flags%Viscosity) D%buffer_size = D%buffer_size + D%sizePr
    if (D%save_flags%Temperature) D%buffer_size = D%buffer_size + D%sizePr
    if (D%save_flags%Moisture) D%buffer_size = D%buffer_size + D%sizePr
    if (D%save_flags%Scalar) D%buffer_size = D%buffer_size + D%sizePr * num_scalars

    if (D%save_flags%Scalar) then
      D%save_flags%num_scalars = num_scalars
    else !probably redundant
      D%save_flags%num_scalars = 0
    end if

    call D%SaveHeader

    call D%SaveMask(Utype, Vtype, Wtype, Prtype)

    allocate(D%buffer(0:D%buffer_size-1))

    D%ranges_set = .true.

  end subroutine TStaggeredFrameDomain_SetRest



  subroutine TStaggeredFrameDomain_Fill(D, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    !Fill the output buffer for asynchronous output
    class(TStaggeredFrameDomain),intent(inout) :: D
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),contiguous,intent(in) :: Pr(-1:,-1:,-1:), Viscosity(-1:,-1:,-1:), &
                                       Temperature(-2:,-2:,-2:), Moisture(-2:,-2:,-2:), &
                                       Scalar(-2:,-2:,-2:,1:)
    integer :: offset
    integer :: i

    if (.not.D%ranges_set) call D%SetRest(num_scalars = size(Scalar,4))

    !assume allocated by SetRest, otherwise error by the runtime library
    
    !offset = 0
    D%buffer(0) = D%times(D%frame_number)
    offset = 1

    if (D%save_flags%U) then
      D%buffer(offset:offset+D%sizeU-1) = &
          U(D%minUi:D%maxUi, D%minUj:D%maxUj, D%minUk:D%maxUk)

      offset = offset + D%sizeU
    end if

    if (D%save_flags%V) then
      D%buffer(offset:offset+D%sizeV-1) = &
          V(D%minVi:D%maxVi, D%minVj:D%maxVj, D%minVk:D%maxVk)

      offset = offset + D%sizeV
    end if

    if (D%save_flags%W) then
      D%buffer(offset:offset+D%sizeW-1) = &
          W(D%minWi:D%maxWi, D%minWj:D%maxWj, D%minWk:D%maxWk)

      offset = offset + D%sizeW
    end if

    if (D%save_flags%Pr) then
      !FIXME: problem near boundaries because of the boundary buffer
      D%buffer(offset:offset+D%sizePr-1) = &
          Pr(D%minPri:D%maxPri, D%minPrj:D%maxPrj, D%minPrk:D%maxPrk)

      offset = offset + D%sizePr
    end if

    if (D%save_flags%Viscosity) then
      D%buffer(offset:offset+D%sizePr-1) = &
          Viscosity(D%minPri:D%maxPri, D%minPrj:D%maxPrj, D%minPrk:D%maxPrk)

      offset = offset + D%sizePr
    end if

    if (D%save_flags%Temperature) then
      D%buffer(offset:offset+D%sizePr-1) = &
          Temperature(D%minPri:D%maxPri, D%minPrj:D%maxPrj, D%minPrk:D%maxPrk)

      offset = offset + D%sizePr
    end if

    if (D%save_flags%Moisture) then
      D%buffer(offset:offset+D%sizePr-1) = &
          Moisture(D%minPri:D%maxPri, D%minPrj:D%maxPrj, D%minPrk:D%maxPrk)

      offset = offset + D%sizePr
    end if

    if (D%save_flags%Scalar) then
      do i=1,size(Scalar,4)
        D%buffer(offset:offset+D%sizePr-1) = &
            Scalar(D%minPri:D%maxPri, D%minPrj:D%maxPrj, D%minPrk:D%maxPrk,i)

        offset = offset + D%sizePr
      end do
    end if

  end subroutine TStaggeredFrameDomain_Fill



  recursive subroutine SaveBuffer(Dptr) bind(C)
    use iso_c_binding, only: c_f_pointer
    type(c_ptr),value :: Dptr
    type(TStaggeredFrameDomain),pointer :: D
    character(2512) :: file_name

    call c_f_pointer(Dptr, D)

    write(file_name,'(a,i0,a)') trim(D%base_name)//"-",D%frame_number,trim(D%suffix)

!     write(*,*) "Writing frame: ",file_name, " size: ",D%buffer_size

    open(unit=D%unit, file=file_name, access='stream', form='unformatted', action='write', status='replace')

    write(D%unit) D%buffer

    close(D%unit)

  end subroutine



  subroutine TStaggeredFrameDomain_DoSave(D)
    use iso_c_binding, only: c_loc, c_funloc
    use stop_procedures, only: error_stop
    class(TStaggeredFrameDomain),target,asynchronous,intent(inout) :: D
    integer :: err
    
    err = 1

    select type (Dnp => D)
      type is (TStaggeredFrameDomain)
        call pthread_create_opaque(Dnp%threadptr, &
                                   c_funloc(SaveBuffer), &
                                   c_loc(Dnp), &
                                   err)
      class default
        call error_stop("Error: wrog type of D in TStaggeredFrameDomain_DoSave.")
    end select

    if (err==0) then
      D%in_progress = .true.
    else
      write (*,*) "Error in creating frame thread. Will run again synchronously. Code:",err
      
      select type (Dnp => D)
        type is (TStaggeredFrameDomain)
          call SaveBuffer(c_loc(Dnp))
        class default
          call error_stop("Error: wrog type of D in TStaggeredFrameDomain_DoSave.")
      end select
    end if
  end subroutine TStaggeredFrameDomain_DoSave



  subroutine TStaggeredFrameDomain_SaveHeader(D)
    class(TStaggeredFrameDomain),intent(in) :: D
    character(2512) :: file_name

    file_name = trim(D%base_name)//"-header"//trim(D%suffix)

    open(unit=D%unit, file=file_name, access='stream', form='unformatted', action='write', status='replace')

    write(D%unit) 1_int32 !endianess can be infered from this
    write(D%unit) int(storage_size(1._knd),int32)  !save number of bits of the used real kind
    write(D%unit) int(storage_size(1),int32)  !save number of bits of the used (default) integer kind
    write(D%unit) int(storage_size(.true.),int32)  !save number of bits of the used (default) logical kind

    write(D%unit) D%save_flags%U, D%save_flags%V, D%save_flags%W
    write(D%unit) D%save_flags%Pr, D%save_flags%Viscosity
    write(D%unit) D%save_flags%Temperature, D%save_flags%Moisture, &
                  D%save_flags%Scalar, D%save_flags%num_scalars

    write(D%unit) D%minUi, D%maxUi, D%minUj, D%maxUj, D%minUk, D%maxUk
    write(D%unit) D%minVi, D%maxVi, D%minVj, D%maxVj, D%minVk, D%maxVk
    write(D%unit) D%minWi, D%maxWi, D%minWj, D%maxWj, D%minWk, D%maxWk
    write(D%unit) D%minPri, D%maxPri, D%minPrj, D%maxPrj, D%minPrk, D%maxPrk

    write(D%unit) D%xU,D%yV,D%zW,D%xPr,D%yPr,D%zPr

    close(D%unit)

  end subroutine


  subroutine TStaggeredFrameDomain_SaveMask(D, Utype,Vtype, Wtype, Prtype)
    class(TStaggeredFrameDomain),intent(in) :: D
    integer,intent(in) :: Prtype(0:,0:,0:)
    integer,dimension(-2:,-2:,-2:),intent(in) :: Utype, Vtype, Wtype
    character(2512) :: file_name

    file_name = trim(D%base_name)//"-mask"//trim(D%suffix)

    open(unit=D%unit, file=file_name, access='stream', form='unformatted', action='write', status='replace')

    write(D%unit) Utype(D%minUi:D%maxUi, D%minUj:D%maxUj, D%minUk:D%maxUk)
    write(D%unit) Vtype(D%minVi:D%maxVi, D%minVj:D%maxVj, D%minVk:D%maxVk)
    write(D%unit) Wtype(D%minWi:D%maxWi, D%minWj: D%maxWj, D%minWk:D%maxWk)
    write(D%unit) Prtype(D%minPri:D%maxPri, D%minPrj:D%maxPrj, D%minPrk:D%maxPrk)

    close(D%unit)

  end subroutine


end module StaggeredFrames






