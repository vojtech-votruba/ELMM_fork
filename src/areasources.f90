module AreaSources
  !In this implementation area sources reduce to a set of single cell volume sources above the surface
  use Parameters
  use VolumeSources
  use GeometricShapes2D
  
  implicit none

  public :: ScalarAreaSource, ScalarAreaSources, InitAreaSources

  private

  type :: ScalarAreaSource
    class(GeometricShape2D), allocatable :: GeometricShape
    integer   :: scalar_number
    real(knd) :: flux !total flux in unit/time_unit, distribution is uniform
  contains
    procedure :: point_sources
  end type
  
  type(ScalarAreaSource), allocatable :: ScalarAreaSources(:)

  interface ScalarAreaSource
    module procedure ScalarAreaSource_Init
  end interface
  
contains

  function ScalarAreaSource_Init(gs, scalar_number, flux) result(res)
    type(ScalarAreaSource) :: res
    class(GeometricShape2D), intent(in) :: gs
    integer, intent(in) :: scalar_number
    real(knd), intent(in) :: flux

    allocate(res%GeometricShape, source=gs)
    res%scalar_number = scalar_number
    res%flux = flux
  end function

  function point_sources(self) result(res)
    use Boundaries
    use ArrayUtilities
#ifdef PAR
    use custom_par
#endif
    type(ScalarFlVolumesContainer) :: res
    class(ScalarAreaSource), intent(in) :: self
    integer :: i, j, k
    real(knd) :: total_volume
    
    res%scalar_number = self%scalar_number
    
    allocate(res%volumes(0))

    total_volume = 0
    
    do j = 1, Prny
      do i = 1, Prnx
        if (self%GeometricShape%Inside(xPr(i), yPr(j))) then
          if (Prtype(i,j,Prnz)>0) cycle
          if (Prtype(i,j,1)<=0 .and. &
              (Btype(Bo)>=BC_MPI_BOUNDS_MIN .and. Btype(Bo)<=BC_MPI_BOUNDS_MAX)) cycle

          do k = 1, Prnz
            if (Prtype(i,j,k)<=0) then
              total_volume = total_volume + dxmin * dymin * dzPr(k)
              !the flux inScalarFlVolume is the local concentration derivative dc/dt
              !but self%flux is the total amount m released per unit time. m = c * total_volume
              res%volumes = [res%volumes , &
                             ScalarFlVolume([i, j, k], self%flux)]
              exit
            end if
          end do

        end if
      end do
    end do

#ifdef PAR
    total_volume = par_co_sum(total_volume)
#endif
    if (total_volume>0) res%volumes%flux = res%volumes%flux / total_volume

  end function
  
  subroutine InitAreaSources
    integer :: i
    type(ScalarFlVolumesContainer) :: points
    
    if (.not.allocated(ScalarAreaSources)) then
    
      allocate(ScalarAreaSources(0))
      
    else if (num_of_scalars>0) then
    
      do i = 1, size(ScalarAreaSources)
        points = ScalarAreaSources(i)%point_sources()
        if (.not.empty(points)) call Add(ScalarFlVolumes, points)
      end do
      
    end if
    
  contains
    logical function empty(src)
      type(ScalarFlVolumesContainer), intent(in) :: src
      empty = .not.allocated(src%volumes)
      if (.not.empty) empty = size(src%volumes)==0
    end function
  end subroutine
    
end module AreaSources
