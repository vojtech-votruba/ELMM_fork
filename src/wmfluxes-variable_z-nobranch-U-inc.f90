!These macros generate lines too long for Fortran.
!This file is therefore used only to generate the 
!files specific to each component which is than adjusted by hand 
!to have lines short enough.
!Command to generate specific files:
!    gfortran -E -cpp -Dcomp=1  wmfluxes-variable_z-nobranch-inc.f90 > wmfluxes-variable_z-nobranch-U-inc.f90
















      !$omp do
      do i = 1, size(UxmWMpoints)

            wrk(UxmWMpoints(i)%xi,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) = &
            wrk(UxmWMpoints(i)%xi,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) + &
  nu((UxmWMpoints(i)%xi-1)+1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) * &
  (U((UxmWMpoints(i)%xi-1)+1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk)-U(UxmWMpoints(i)%xi-1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk)) * &
  recdxmin2

            wrk(UxmWMpoints(i)%xi,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) = &
            wrk(UxmWMpoints(i)%xi,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) + &
              UxmWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(UxmWMpoints)

            wrk(UxmWMpoints(i)%xi-1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) = &
            wrk(UxmWMpoints(i)%xi-1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) - &
  nu((UxmWMpoints(i)%xi-1)+1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) * &
  (U((UxmWMpoints(i)%xi-1)+1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk)-U(UxmWMpoints(i)%xi-1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk)) * &
  recdxmin2 / dzPr(UxmWMpoints(i)%zk-1)

            wrk(UxmWMpoints(i)%xi-1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) = &
            wrk(UxmWMpoints(i)%xi-1,UxmWMpoints(i)%yj,UxmWMpoints(i)%zk) - &
              UxmWMpoints(i)%fluxm

      end do
      !$omp end do

      !$omp do
      do i = 1, size(UxpWMpoints)

            wrk(UxpWMpoints(i)%xi+1,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) = &
            wrk(UxpWMpoints(i)%xi+1,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) + &
  nu((UxpWMpoints(i)%xi)+1,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) * &
  (U((UxpWMpoints(i)%xi)+1,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk)-U(UxpWMpoints(i)%xi,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk)) * &
  recdxmin2

            wrk(UxpWMpoints(i)%xi+1,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) = &
            wrk(UxpWMpoints(i)%xi+1,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) + &
              UxpWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(UxpWMpoints)

            wrk(UxpWMpoints(i)%xi,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) = &
            wrk(UxpWMpoints(i)%xi,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) - &
  nu((UxpWMpoints(i)%xi)+1,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) * &
  (U((UxpWMpoints(i)%xi)+1,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk)-U(UxpWMpoints(i)%xi,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk)) * &
  recdxmin2

            wrk(UxpWMpoints(i)%xi,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) = &
            wrk(UxpWMpoints(i)%xi,UxpWMpoints(i)%yj,UxpWMpoints(i)%zk) - &
              UxpWMpoints(i)%fluxm

      end do
      !$omp end do

      !$omp do
      do i = 1, size(UymWMpoints)

            wrk(UymWMpoints(i)%xi,UymWMpoints(i)%yj,UymWMpoints(i)%zk) = &
            wrk(UymWMpoints(i)%xi,UymWMpoints(i)%yj,UymWMpoints(i)%zk) + &
  (nu((UymWMpoints(i)%xi)+1,(UymWMpoints(i)%yj-1)+1,UymWMpoints(i)%zk)+&
  nu((UymWMpoints(i)%xi)+1,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk)+&
  nu(UymWMpoints(i)%xi,(UymWMpoints(i)%yj-1)+1,UymWMpoints(i)%zk)+&
  nu(UymWMpoints(i)%xi,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk))*&
  (U(UymWMpoints(i)%xi,(UymWMpoints(i)%yj-1)+1,UymWMpoints(i)%zk)-U(UymWMpoints(i)%xi,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk))*&
  recdymin2/4

            wrk(UymWMpoints(i)%xi,UymWMpoints(i)%yj,UymWMpoints(i)%zk) = &
            wrk(UymWMpoints(i)%xi,UymWMpoints(i)%yj,UymWMpoints(i)%zk) + &
              UymWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(UymWMpoints)

            wrk(UymWMpoints(i)%xi,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk) = &
            wrk(UymWMpoints(i)%xi,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk) - &
  (nu((UymWMpoints(i)%xi)+1,(UymWMpoints(i)%yj-1)+1,UymWMpoints(i)%zk)+&
  nu((UymWMpoints(i)%xi)+1,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk)+&
  nu(UymWMpoints(i)%xi,(UymWMpoints(i)%yj-1)+1,UymWMpoints(i)%zk)+&
  nu(UymWMpoints(i)%xi,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk))*&
  (U(UymWMpoints(i)%xi,(UymWMpoints(i)%yj-1)+1,UymWMpoints(i)%zk)-U(UymWMpoints(i)%xi,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk))*&
  recdymin2/4

            wrk(UymWMpoints(i)%xi,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk) = &
            wrk(UymWMpoints(i)%xi,UymWMpoints(i)%yj-1,UymWMpoints(i)%zk) - &
              UymWMpoints(i)%fluxm

      end do
      !$omp end do

      !$omp do
      do i = 1, size(UypWMpoints)

            wrk(UypWMpoints(i)%xi,UypWMpoints(i)%yj+1,UypWMpoints(i)%zk) = &
            wrk(UypWMpoints(i)%xi,UypWMpoints(i)%yj+1,UypWMpoints(i)%zk) + &
  (nu((UypWMpoints(i)%xi)+1,(UypWMpoints(i)%yj)+1,UypWMpoints(i)%zk)+&
  nu((UypWMpoints(i)%xi)+1,UypWMpoints(i)%yj,UypWMpoints(i)%zk)+&
  nu(UypWMpoints(i)%xi,(UypWMpoints(i)%yj)+1,UypWMpoints(i)%zk)+&
  nu(UypWMpoints(i)%xi,UypWMpoints(i)%yj,UypWMpoints(i)%zk))*&
  (U(UypWMpoints(i)%xi,(UypWMpoints(i)%yj)+1,UypWMpoints(i)%zk)-U(UypWMpoints(i)%xi,UypWMpoints(i)%yj,UypWMpoints(i)%zk))*&
  recdymin2/4

            wrk(UypWMpoints(i)%xi,UypWMpoints(i)%yj+1,UypWMpoints(i)%zk) = &
            wrk(UypWMpoints(i)%xi,UypWMpoints(i)%yj+1,UypWMpoints(i)%zk) + &
              UypWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(UypWMpoints)

            wrk(UypWMpoints(i)%xi,UypWMpoints(i)%yj,UypWMpoints(i)%zk) = &
            wrk(UypWMpoints(i)%xi,UypWMpoints(i)%yj,UypWMpoints(i)%zk) - &
  (nu((UypWMpoints(i)%xi)+1,(UypWMpoints(i)%yj)+1,UypWMpoints(i)%zk)+&
  nu((UypWMpoints(i)%xi)+1,UypWMpoints(i)%yj,UypWMpoints(i)%zk)+&
  nu(UypWMpoints(i)%xi,(UypWMpoints(i)%yj)+1,UypWMpoints(i)%zk)+&
  nu(UypWMpoints(i)%xi,UypWMpoints(i)%yj,UypWMpoints(i)%zk))*&
  (U(UypWMpoints(i)%xi,(UypWMpoints(i)%yj)+1,UypWMpoints(i)%zk)-U(UypWMpoints(i)%xi,UypWMpoints(i)%yj,UypWMpoints(i)%zk))*&
  recdymin2/4

            wrk(UypWMpoints(i)%xi,UypWMpoints(i)%yj,UypWMpoints(i)%zk) = &
            wrk(UypWMpoints(i)%xi,UypWMpoints(i)%yj,UypWMpoints(i)%zk) - &
              UypWMpoints(i)%fluxm

      end do
      !$omp end do

      !$omp do
      do i = 1, size(UzmWMpoints)

            wrk(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk) = &
            wrk(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk) + &
  (nu((UzmWMpoints(i)%xi)+1,UzmWMpoints(i)%yj,(UzmWMpoints(i)%zk-1)+1)+nu((UzmWMpoints(i)%xi)+1,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1)+nu(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,(UzmWMpoints(i)%zk-1)+1)+nu(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1))*(U(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,(UzmWMpoints(i)%zk-1)+1)-U(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1))/dzW(UzmWMpoints(i)%zk-1)/4 / dzPr(UzmWMpoints(i)%zk)

            wrk(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk) = &
            wrk(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk) + &
              UzmWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(UzmWMpoints)

            wrk(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1) = &
            wrk(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1) - &
  (nu((UzmWMpoints(i)%xi)+1,UzmWMpoints(i)%yj,(UzmWMpoints(i)%zk-1)+1)+nu((UzmWMpoints(i)%xi)+1,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1)+nu(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,(UzmWMpoints(i)%zk-1)+1)+nu(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1))*(U(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,(UzmWMpoints(i)%zk-1)+1)-U(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1))/dzW(UzmWMpoints(i)%zk-1)/4 / dzPr(UzmWMpoints(i)%zk-1)

            wrk(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1) = &
            wrk(UzmWMpoints(i)%xi,UzmWMpoints(i)%yj,UzmWMpoints(i)%zk-1) - &
              UzmWMpoints(i)%fluxm

      end do
      !$omp end do
      !$omp do
      do i = 1, size(UzpWMpoints)

            wrk(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk+1) = &
            wrk(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk+1) + &
  (nu((UzpWMpoints(i)%xi)+1,UzpWMpoints(i)%yj,(UzpWMpoints(i)%zk)+1)+nu((UzpWMpoints(i)%xi)+1,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk)+nu(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,(UzpWMpoints(i)%zk)+1)+nu(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk))*(U(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,(UzpWMpoints(i)%zk)+1)-U(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk))/dzW(UzpWMpoints(i)%zk)/4 / dzPr(UzpWMpoints(i)%zk+1)

            wrk(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk+1) = &
            wrk(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk+1) + &
              UzpWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(UzpWMpoints)

            wrk(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk) = &
            wrk(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk) - &
  (nu((UzpWMpoints(i)%xi)+1,UzpWMpoints(i)%yj,(UzpWMpoints(i)%zk)+1)+nu((UzpWMpoints(i)%xi)+1,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk)+nu(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,(UzpWMpoints(i)%zk)+1)+nu(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk))*(U(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,(UzpWMpoints(i)%zk)+1)-U(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk))/dzW(UzpWMpoints(i)%zk)/4 / dzPr(UzpWMpoints(i)%zk)

            wrk(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk) = &
            wrk(UzpWMpoints(i)%xi,UzpWMpoints(i)%yj,UzpWMpoints(i)%zk) - &
              UzpWMpoints(i)%fluxm

      end do
      !$omp end do
      
