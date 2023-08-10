module LineSources
  !In this implementation line sources reduce to a set of single cell volume sources
  use Parameters
  use VolumeSources
  
  implicit none

  public :: ScalarLineSource, ScalarLineSources, InitLineSources

  private

  type :: ScalarLineSource
    integer   :: scalar_number
    real(knd) :: start(3), end(3)
    real(knd) :: flux !total flux in unit/time_unit, distribuion is uniform
    integer   :: number_of_points = 10000 !number of virtual point sources
  contains
    procedure :: point_sources
  end type
  
  type(ScalarLineSource), allocatable :: ScalarLineSources(:)

  
  !TODO: line sources with nonuniform intensity could be useful, use a similar (or extended?)
  ! class with a function pointer
  
  contains

    function point_sources(self) result(res)
      use Boundaries
      use ArrayUtilities
      type(ScalarFlVolumesContainer) :: res
      class(ScalarLineSource), intent(in) :: self
      real(knd) :: rpos(3), flux_per_point
      integer :: ipos(3), last_ipos(3)
      integer :: i
      
      res%scalar_number = self%scalar_number
      
      allocate(res%volumes(0))
      
      last_ipos = -huge(1)
      
      flux_per_point = self%flux * dist(self%end, self%start) / self%number_of_points

      do i = 0, self%number_of_points-1
        rpos = self%start + ((i+1/2._knd)*(self%end-self%start))/self%number_of_points
        
        if (InDomain(rpos)) then
        
          call GridCoords(ipos,rpos)

          if (all(ipos==last_ipos)) then

            associate(vs => res%volumes)
              vs(size(vs))%flux = vs(size(vs))%flux + flux_per_point / volumePr(ipos)
            end associate

          else
            res%volumes = [res%volumes , &
                           ScalarFlVolume(ipos, flux_per_point / volumePr(ipos))]
          end if
          
          last_ipos = ipos
          
        end if
      end do

    end function
    
    subroutine InitLineSources
      integer :: i
      type(ScalarFlVolumesContainer) :: points
      
      if (.not.allocated(ScalarLineSources)) then
      
        allocate(ScalarLineSources(0))
        
      else if (num_of_scalars>0) then
      
        do i = 1, size(ScalarLineSources)
          points = ScalarLineSources(i)%point_sources()
          if (.not.empty(points)) call Add(ScalarFlVolumes,points)
        end do
        
      end if
      
    contains
      logical function empty(src)
        type(ScalarFlVolumesContainer), intent(in) :: src
        empty = .not.allocated(src%volumes)
        if (.not.empty) empty = size(src%volumes)==0
      end function
    end subroutine
    
end module LineSources