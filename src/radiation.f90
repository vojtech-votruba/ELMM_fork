module Radiation
  use Parameters

  implicit none
  
  !in future it may point to other modules implementing some actual radiation schemes

  logical :: enable_radiation_profile = .false.
  
  real(knd), allocatable :: rad_flux_convergence_prof(:)

contains

  subroutine apply_radiation_heat_profile(Temperature)
    !a simple profile of radiation heating (flux convergence)
    real(knd), contiguous, intent(inout) :: Temperature(-1:,-1:,-1:)
    integer :: i ,j, k

    !$omp parallel do collapse(3)
    do k = 1, Prnz
      do j = 1, Prny
        do i = 1, Prnx
          Temperature(i,j,k) = Temperature(i,j,k) + rad_flux_convergence_prof(k)
        end do
      end do
    end do
  end subroutine

end module Radiation

