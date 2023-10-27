module Puffs

  use Parameters
  
  use ArrayUtilities
  
  use VolumeSources, only: ScalarFlVolumesContainer, ScalarFlVolume
  
  implicit none
  
  private
  
  public PuffSource, PuffSources, PreparePuffs, DoPuffs, InitPuffSources, enable_puffs
  
  type PuffSource
    real(knd) :: start_time,end_time,period, duration, last_start = -huge(1._knd)
    logical   :: clear !clear the scalar field? (all scalars for which source below exists)
    type(ScalarFlVolumesContainer),allocatable :: sources(:) !must be allocated to 0 if empty
    logical, private :: restart = .false.
    real(knd), private :: dt_fraction = 0
  end type

  type(PuffSource),allocatable :: PuffSources(:) !index means individual puff source objects
  
  integer :: enable_puffs = 1

  
contains

  subroutine PreparePuffs(Scalar, RK_stage, RK_stages, time, dt)
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(inout) :: Scalar
    integer, intent(in) :: RK_stage, RK_stages
    real(knd), intent(in) :: time, dt
    integer :: i
  
    if (allocated(PuffSources)) then
      do i = 1,size(PuffSources)
        call aux(PuffSources(i))
      end do
    end if
    
  contains
  
    subroutine aux(p)
      type(PuffSource), intent(inout) :: p
      integer :: j
      !we assume period >> dt and duration >> dt in all following
          
      if (RK_stage==1) then
        if (time+dt > p%last_start + p%period .and. time+dt >= p%start_time .and. time <= p%end_time) then
          p%restart = .true.
        else
          p%restart = .false.
        end if
        

        if ((.not.p%restart) .and. (time <= p%last_start + p%duration) .and. (time >= p%start_time)) then
          p%dt_fraction = max(min(p%last_start + p%duration, time + dt, p%end_time) - time, 0._knd) / dt
        else if (p%restart) then
          p%dt_fraction = time + dt - max(p%last_start + p%period, time, p%start_time) / dt
        else
          p%dt_fraction = 0
        end if
        
        if (p%clear.and.p%restart.and.RK_stage==1) then
          do j=1,size(p%sources)
            associate(src => p%sources(j))
              if (src%scalar_number>0 .and. src%scalar_number<=num_of_scalars) then
                  call set(Scalar(:,:,:,src%scalar_number), 0)
              end if
            end associate
          end do
        end if
      end if
      
      if (p%restart.and.RK_stage==RK_stages) then
        if (p%last_start>=p%start_time) then
          p%last_start = p%last_start + p%period
        else
          p%last_start = p%start_time
        end if
      end if
    end subroutine
  end subroutine

  subroutine DoPuffs(Scalar_derivative, sc)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: Scalar_derivative
    integer, intent(in) :: sc
    integer :: i
      
    if (allocated(PuffSources)) then
      do i = 1,size(PuffSources)
        call DoPuff(PuffSources(i), Scalar_derivative, sc)
      end do
    end if
  end subroutine
  
  
  subroutine DoPuff(p, Scalar_derivative, sc)
    type(PuffSource),intent(inout) :: p
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: Scalar_derivative
    integer, intent(in) :: sc
    integer :: i,j
  
    if (p%dt_fraction>0) then
      do i=1,size(p%sources)

        associate(src => p%sources(i))
          if (src%scalar_number==sc) then
          
            do j=1,size(src%volumes)
              associate (xi => src%volumes(j)%xi, &
                         yj => src%volumes(j)%yj, &
                         zk => src%volumes(j)%zk)
                  Scalar_derivative(xi,yj,zk) = Scalar_derivative(xi,yj,zk) + &
                                                  p%dt_fraction * src%volumes(j)%flux
              end associate
            end do
            
          end if
        end associate

      end do
    end if

  end subroutine
  
  
  subroutine InitPuffSources
   use Boundaries, only: GridCoords, InDomain
#if CUSTOM_PUFF_SOURCES
   interface
     subroutine  CustomPuffSources
     end subroutine
   end interface
   
   call CustomPuffSources
#endif
  end subroutine



end module Puffs
