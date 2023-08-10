module Large_scale_processes
  use Parameters

  implicit none
  
  !in future it may point to other modules implementing some actual radiation schemes
  
  logical :: enable_large_scale_processes = .false.

  logical :: enable_temperature_advection_profile = .false.
  logical :: enable_moisture_advection_profile = .false.
  
  real(knd), allocatable :: temperature_advection_prof(:)
  real(knd), allocatable :: moisture_advection_prof(:)
  
  
  logical :: enable_subsidence_profile = .false.
  
  real(knd), allocatable :: subsidence_profile(:)

contains

  subroutine apply_temperature_advection_profile(Temperature)
    real(knd), contiguous, intent(inout) :: Temperature(-1:,-1:,-1:)
    integer :: i ,j, k

    !$omp parallel do collapse(3)
    do k = 1, Prnz
      do j = 1, Prny
        do i = 1, Prnx
          Temperature(i,j,k) = Temperature(i,j,k) + temperature_advection_prof(k)
        end do
      end do
    end do
  end subroutine

  subroutine apply_moisture_advection_profile(Moisture)
    real(knd), contiguous, intent(inout) :: Moisture(-1:,-1:,-1:)
    integer :: i ,j, k

    !$omp parallel do collapse(3)
    do k = 1, Prnz
      do j = 1, Prny
        do i = 1, Prnx
          Moisture(i,j,k) = Moisture(i,j,k) + moisture_advection_prof(k)
        end do
      end do
    end do
  end subroutine
  
  
  
  subroutine apply_subsidence_profile(Array_adv, Array)
    !the subsidence profile contains velocities <w>
    !positive means upwards
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(inout) :: Array_adv, Array
    integer :: i, j, k
    
    !$omp parallel do private(i,j,k)
    do k = 1, Prnz
      do j = 1, Prny
        do i = 1, Prnx
          Array_adv(i,j,k) = Array_adv(i,j,k) - subsidence_profile(k) * (Array(i,j,k+1) - Array(i,j,k-1)) / (zPr(k+1)-zPr(k-1))
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

 
  
  subroutine temperature_large_scale_terms(Temperature)
    real(knd), contiguous, intent(inout) :: Temperature(-1:,-1:,-1:)

    if (enable_temperature_advection_profile) then
      call apply_temperature_advection_profile(Temperature)
    end if

    !subsidence treated separately, because it requires the current solution values
  end subroutine


  subroutine moisture_large_scale_terms(Moisture)
    real(knd), contiguous, intent(inout) :: Moisture(-1:,-1:,-1:)

    if (enable_moisture_advection_profile) then
      call apply_moisture_advection_profile(Moisture)
    end if

    !subsidence treated separately, because it requires the current solution values
  end subroutine



end module Large_scale_processes

