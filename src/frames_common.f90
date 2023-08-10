module Frames_common

  use Kinds
  use iso_c_binding, only: c_ptr
  use Pthreads

  type TFrameTimes
    real(knd) :: start, end
    integer :: nframes
  end type
  
  type bounding_box
    real(knd) :: xmin = -huge(1._knd)/2, xmax = huge(1._knd)
    real(knd) :: ymin = -huge(1._knd)/2, ymax = huge(1._knd)
    real(knd) :: zmin = -huge(1._knd)/2, zmax = huge(1._knd)
  end type

  


  type, abstract :: TFrameBase

    type(TFrameTimes) :: frame_times

    integer   :: frame_number = -1 !number of the previous frame

    real(knd), allocatable :: times(:)

    integer   :: minPri, maxPri, minPrj, maxPrj, minPrk, maxPrk

    integer   :: sizePr

    real(knd),allocatable,dimension(:) :: xPr,yPr,zPr
    
    type(bounding_box), allocatable :: bbox
    
    logical :: active = .true.

    logical :: ranges_set = .false.

    character(2512) :: base_name = ""

    integer :: unit

    logical :: in_progress = .false.

    type(c_ptr) :: threadptr
  contains
    procedure(dosave_interface), deferred :: DoSave
    procedure(fill_interface), deferred :: Fill
    procedure :: SaveTimes => TFrameBase_SaveTimes
    procedure :: Save => TFrameBase_Save
    procedure :: Wait => TFrameBase_Wait
    procedure :: Finalize => TFrameBase_Finalize
  end type
  
  abstract interface
    subroutine fill_interface(D, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
      import
      class(TFrameBase),intent(inout) :: D
      real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
      real(knd),contiguous,intent(in) :: Pr(-1:,-1:,-1:), Viscosity(-1:,-1:,-1:), &
                              Temperature(-1:,-1:,-1:), Moisture(-1:,-1:,-1:), &
                              Scalar(-1:,-1:,-1:,1:)
    end subroutine
    subroutine dosave_interface(D)
      import
      class(TFrameBase),target,asynchronous,intent(inout) :: D
    end subroutine
  end interface

  interface Add
    module procedure TimeSeries_Add
  end interface
  
contains

  subroutine TimeSeries_Add(array,idx,val)
    real(knd),allocatable,intent(inout) :: array(:)
    integer,intent(in)   :: idx
    real(knd),intent(in) :: val
    real(knd),allocatable :: tmp(:)

    !assume allocated, otherwise error by the runtime library

    if (idx < lbound(array,1)) then
      write(*,*) "Error. Tried to write under the start of the series."
      stop
    end if
    if (idx > ubound(array,1)) then
      allocate( tmp(lbound(array,1):lbound(array,1)+(idx-lbound(array,1))*2) )
      tmp(lbound(array,1):ubound(array,1)) = array
      deallocate(array)
      call move_alloc(tmp, array)
    end if
    
    array(idx) = val

  end subroutine TimeSeries_Add
  
  
  subroutine TFrameBase_SaveTimes(D)
    class(TFrameBase),intent(in) :: D
    character(2512) :: file_name
    integer :: i

    file_name = trim(D%base_name)//"-times.txt"

    open(unit=D%unit, file=file_name, access='sequential', form='formatted', action='write', status='replace')

    do i = 0, D%frame_number
      write(D%unit,'(i8,2x,es12.5)') i, D%times(i)
    end do
    close(D%unit)

  end subroutine
  
  
  
  subroutine TFrameBase_Save(D, time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    class(TFrameBase),target,asynchronous,intent(inout) :: D
    real(knd),intent(in) :: time
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),contiguous,intent(in) :: Pr(-1:,-1:,-1:), &
                                       Temperature(-1:,-1:,-1:), Viscosity(-1:,-1:,-1:), &
                                       Moisture(-1:,-1:,-1:), Scalar(-1:,-1:,-1:,1:)

    associate(start   => D%frame_times%start,&
              end     => D%frame_times%end,&
              nframes => D%frame_times%nframes)

      if ( (time >= start) .and. (time <= end + (end-start)/(nframes-1)) &
        .and. (time >= start + (D%frame_number+1)*(end-start) / (nframes-1)) ) then

        D%frame_number = D%frame_number + 1

        if (.not.allocated(D%times)) allocate(D%times(D%frame_number:D%frame_number+100))

        call Add(D%times,D%frame_number,time)

        if (D%in_progress) call D%Wait

        call D%Fill(U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)

        call D%DoSave

      end if

    end associate
  end subroutine TFrameBase_Save
  
  
  subroutine TFrameBase_Wait(D)
    class(TFrameBase),asynchronous,intent(inout) :: D
    integer :: err

    call pthread_join_opaque(D%threadptr,err)

    if (err/=0) then
      write (*,*) "Error in joining staggered frame thread. Will try to continue anyway. Code:",err
    end if

    D%in_progress = .false.
  end subroutine



  subroutine TFrameBase_Finalize(D)
    class(TFrameBase),intent(inout) :: D

    if (D%in_progress) call D%Wait

    call D%SaveTimes

  end subroutine


end module
