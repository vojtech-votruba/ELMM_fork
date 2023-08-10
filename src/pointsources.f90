module PointSources
  !In this implementation point sources reduce to a single cell volume source
  use Parameters
  use VolumeSources
  
  implicit none

  type :: ScalarPointSource
    integer   :: scalar_number
    real(knd) :: position(3)
    real(knd) :: flux !total flux in unit/time_unit
  contains
    procedure :: point_source
  end type
  
  type(ScalarPointSource), allocatable :: ScalarPointSources(:)
  
  contains

    function point_source(self) result(res)
      use Boundaries
      use ArrayUtilities
      type(ScalarFlVolumesContainer) :: res
      class(ScalarPointSource), intent(in) :: self
      real(knd) :: rpos(3)
      integer :: ipos(3)
      
      res%scalar_number = self%scalar_number
      
      allocate(res%volumes(0))
      
      rpos = self%position
      
      if (InDomain(rpos)) then
      
        call GridCoords(ipos,rpos)

        res%volumes =  [ScalarFlVolume(ipos, self%flux / volumePr(ipos))]
      end if

    end function
    
    subroutine InitPointSources
      integer :: i
      type(ScalarFlVolumesContainer) :: point
      
      if (.not.allocated(ScalarPointSources)) then
      
        allocate(ScalarPointSources(0))
        
      else if (num_of_scalars>0) then
      
        do i=1,size(ScalarPointSources)
          point = ScalarPointSources(i)%point_source()
          if (.not.empty(point)) call Add(ScalarFlVolumes,point)
        end do
        
      end if
      
    contains
      logical function empty(src)
        type(ScalarFlVolumesContainer), intent(in) :: src
        empty = .not.allocated(src%volumes)
        if (.not.empty) empty = size(src%volumes)==0
      end function
    end subroutine
    
end module PointSources
