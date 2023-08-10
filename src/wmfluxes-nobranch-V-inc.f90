!Automatically generated from wmfluxes-nobranch-inc.f90
!Do not edit directly if the original template has to be changed.

      !$omp do
      do i = 1, size(VxmWMpoints)

            wrk(VxmWMpoints(i)%xi,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk) = &
            wrk(VxmWMpoints(i)%xi,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk) + &
  (nu((VxmWMpoints(i)%xi-1)+1,(VxmWMpoints(i)%yj)+1,VxmWMpoints(i)%zk)+&
   nu((VxmWMpoints(i)%xi-1)+1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk)+&
   nu(VxmWMpoints(i)%xi-1,(VxmWMpoints(i)%yj)+1,VxmWMpoints(i)%zk)+&
   nu(VxmWMpoints(i)%xi-1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk))&
  *(V((VxmWMpoints(i)%xi-1)+1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk)-V(VxmWMpoints(i)%xi-1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk))&
  *recdxmin2/4

            wrk(VxmWMpoints(i)%xi,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk) = &
            wrk(VxmWMpoints(i)%xi,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk) + &
              VxmWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(VxmWMpoints)

            wrk(VxmWMpoints(i)%xi-1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk) = &
            wrk(VxmWMpoints(i)%xi-1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk) - &
  (nu((VxmWMpoints(i)%xi-1)+1,(VxmWMpoints(i)%yj)+1,VxmWMpoints(i)%zk)+&
   nu((VxmWMpoints(i)%xi-1)+1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk)+&
   nu(VxmWMpoints(i)%xi-1,(VxmWMpoints(i)%yj)+1,VxmWMpoints(i)%zk)+&
   nu(VxmWMpoints(i)%xi-1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk))&
  *(V((VxmWMpoints(i)%xi-1)+1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk)-V(VxmWMpoints(i)%xi-1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk))&
  *recdxmin2/4

            wrk(VxmWMpoints(i)%xi-1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk) = &
            wrk(VxmWMpoints(i)%xi-1,VxmWMpoints(i)%yj,VxmWMpoints(i)%zk) - &
              VxmWMpoints(i)%fluxp

      end do
      !$omp end do

      !$omp do
      do i = 1, size(VxpWMpoints)
      
            wrk(VxpWMpoints(i)%xi+1,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk) = &
            wrk(VxpWMpoints(i)%xi+1,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk) + &
  (nu((VxpWMpoints(i)%xi)+1,(VxpWMpoints(i)%yj)+1,VxpWMpoints(i)%zk)+&
   nu((VxpWMpoints(i)%xi)+1,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk)+&
   nu(VxpWMpoints(i)%xi,(VxpWMpoints(i)%yj)+1,VxpWMpoints(i)%zk)+&
   nu(VxpWMpoints(i)%xi,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk))&
  *(V((VxpWMpoints(i)%xi)+1,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk)-V(VxpWMpoints(i)%xi,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk))&
  *recdxmin2/4

            wrk(VxpWMpoints(i)%xi+1,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk) = &
            wrk(VxpWMpoints(i)%xi+1,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk) + &
              VxpWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(VxpWMpoints)

            wrk(VxpWMpoints(i)%xi,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk) = &
            wrk(VxpWMpoints(i)%xi,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk) - &
  (nu((VxpWMpoints(i)%xi)+1,(VxpWMpoints(i)%yj)+1,VxpWMpoints(i)%zk)+&
   nu((VxpWMpoints(i)%xi)+1,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk)+&
   nu(VxpWMpoints(i)%xi,(VxpWMpoints(i)%yj)+1,VxpWMpoints(i)%zk)+&
   nu(VxpWMpoints(i)%xi,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk))&
  *(V((VxpWMpoints(i)%xi)+1,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk)-V(VxpWMpoints(i)%xi,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk))&
  *recdxmin2/4

            wrk(VxpWMpoints(i)%xi,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk) = &
            wrk(VxpWMpoints(i)%xi,VxpWMpoints(i)%yj,VxpWMpoints(i)%zk) - &
              VxpWMpoints(i)%fluxp

      end do
      !$omp end do

      !$omp do
      do i = 1, size(VymWMpoints)

            wrk(VymWMpoints(i)%xi,VymWMpoints(i)%yj,VymWMpoints(i)%zk) = &
            wrk(VymWMpoints(i)%xi,VymWMpoints(i)%yj,VymWMpoints(i)%zk) + &
  nu(VymWMpoints(i)%xi,(VymWMpoints(i)%yj-1)+1,VymWMpoints(i)%zk) * &
  (V(VymWMpoints(i)%xi,(VymWMpoints(i)%yj-1)+1,VymWMpoints(i)%zk)-V(VymWMpoints(i)%xi,VymWMpoints(i)%yj-1,VymWMpoints(i)%zk))&
  * recdymin2

            wrk(VymWMpoints(i)%xi,VymWMpoints(i)%yj,VymWMpoints(i)%zk) = &
            wrk(VymWMpoints(i)%xi,VymWMpoints(i)%yj,VymWMpoints(i)%zk) + &
              VymWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(VymWMpoints)

            wrk(VymWMpoints(i)%xi,VymWMpoints(i)%yj-1,VymWMpoints(i)%zk) = &
            wrk(VymWMpoints(i)%xi,VymWMpoints(i)%yj-1,VymWMpoints(i)%zk) - &
  nu(VymWMpoints(i)%xi,(VymWMpoints(i)%yj-1)+1,VymWMpoints(i)%zk) * &
  (V(VymWMpoints(i)%xi,(VymWMpoints(i)%yj-1)+1,VymWMpoints(i)%zk)-V(VymWMpoints(i)%xi,VymWMpoints(i)%yj-1,VymWMpoints(i)%zk))&
  * recdymin2

            wrk(VymWMpoints(i)%xi,VymWMpoints(i)%yj-1,VymWMpoints(i)%zk) = &
            wrk(VymWMpoints(i)%xi,VymWMpoints(i)%yj-1,VymWMpoints(i)%zk) - &
              VymWMpoints(i)%fluxp

      end do
      !$omp end do

      !$omp do
      do i = 1, size(VypWMpoints)
            wrk(VypWMpoints(i)%xi,VypWMpoints(i)%yj+1,VypWMpoints(i)%zk) = &
            wrk(VypWMpoints(i)%xi,VypWMpoints(i)%yj+1,VypWMpoints(i)%zk) + &
  nu(VypWMpoints(i)%xi,(VypWMpoints(i)%yj)+1,VypWMpoints(i)%zk) * &
  (V(VypWMpoints(i)%xi,(VypWMpoints(i)%yj)+1,VypWMpoints(i)%zk)-V(VypWMpoints(i)%xi,VypWMpoints(i)%yj,VypWMpoints(i)%zk))&
  * recdymin2

            wrk(VypWMpoints(i)%xi,VypWMpoints(i)%yj+1,VypWMpoints(i)%zk) = &
            wrk(VypWMpoints(i)%xi,VypWMpoints(i)%yj+1,VypWMpoints(i)%zk) + &
              VypWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(VypWMpoints)

            wrk(VypWMpoints(i)%xi,VypWMpoints(i)%yj,VypWMpoints(i)%zk) = &
            wrk(VypWMpoints(i)%xi,VypWMpoints(i)%yj,VypWMpoints(i)%zk) - &
  nu(VypWMpoints(i)%xi,(VypWMpoints(i)%yj)+1,VypWMpoints(i)%zk) * &
  (V(VypWMpoints(i)%xi,(VypWMpoints(i)%yj)+1,VypWMpoints(i)%zk)-V(VypWMpoints(i)%xi,VypWMpoints(i)%yj,VypWMpoints(i)%zk))&
  * recdymin2

            wrk(VypWMpoints(i)%xi,VypWMpoints(i)%yj,VypWMpoints(i)%zk) = &
            wrk(VypWMpoints(i)%xi,VypWMpoints(i)%yj,VypWMpoints(i)%zk) - &
              VypWMpoints(i)%fluxp

      end do
      !$omp end do

      !$omp do
      do i = 1, size(VzmWMpoints)

            wrk(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk) = &
            wrk(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk) + &
  (nu(VzmWMpoints(i)%xi,(VzmWMpoints(i)%yj)+1,(VzmWMpoints(i)%zk-1)+1)+&
   nu(VzmWMpoints(i)%xi,(VzmWMpoints(i)%yj)+1,VzmWMpoints(i)%zk-1)+&
   nu(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,(VzmWMpoints(i)%zk-1)+1)+&
   nu(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk-1))&
  *(V(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,(VzmWMpoints(i)%zk-1)+1)-V(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk-1))&
  *recdzmin2/4

            wrk(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk) = &
            wrk(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk) + &
              VzmWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(VzmWMpoints)

            wrk(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk-1) = &
            wrk(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk-1) - &
  (nu(VzmWMpoints(i)%xi,(VzmWMpoints(i)%yj)+1,(VzmWMpoints(i)%zk-1)+1)+&
   nu(VzmWMpoints(i)%xi,(VzmWMpoints(i)%yj)+1,VzmWMpoints(i)%zk-1)+&
   nu(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,(VzmWMpoints(i)%zk-1)+1)+&
   nu(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk-1))&
  *(V(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,(VzmWMpoints(i)%zk-1)+1)-V(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk-1))&
  *recdzmin2/4

            wrk(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk-1) = &
            wrk(VzmWMpoints(i)%xi,VzmWMpoints(i)%yj,VzmWMpoints(i)%zk-1) - &
              VzmWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(VzpWMpoints)

            wrk(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk+1) = &
            wrk(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk+1) + &
  (nu(VzpWMpoints(i)%xi,(VzpWMpoints(i)%yj)+1,(VzpWMpoints(i)%zk)+1)+&
   nu(VzpWMpoints(i)%xi,(VzpWMpoints(i)%yj)+1,VzpWMpoints(i)%zk)+&
   nu(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,(VzpWMpoints(i)%zk)+1)+&
   nu(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk))&
  *(V(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,(VzpWMpoints(i)%zk)+1)-V(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk))&
  *recdzmin2/4

            wrk(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk+1) = &
            wrk(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk+1) + &
              VzpWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(VzpWMpoints)

            wrk(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk) = &
            wrk(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk) - &
  (nu(VzpWMpoints(i)%xi,(VzpWMpoints(i)%yj)+1,(VzpWMpoints(i)%zk)+1)+&
   nu(VzpWMpoints(i)%xi,(VzpWMpoints(i)%yj)+1,VzpWMpoints(i)%zk)+&
   nu(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,(VzpWMpoints(i)%zk)+1)+&
   nu(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk))&
  *(V(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,(VzpWMpoints(i)%zk)+1)-V(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk))&
  *recdzmin2/4

            wrk(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk) = &
            wrk(VzpWMpoints(i)%xi,VzpWMpoints(i)%yj,VzpWMpoints(i)%zk) - &
              VzpWMpoints(i)%fluxp

      end do
      !$omp end do
     
