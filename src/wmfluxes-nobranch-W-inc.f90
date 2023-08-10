!Automatically generated from wmfluxes-nobranch-inc.f90
!Do not edit directly if the original template has to be changed.

      !$omp do
      do i = 1, size(WxmWMpoints)

            wrk(WxmWMpoints(i)%xi,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk) = &
            wrk(WxmWMpoints(i)%xi,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk) + &
  (nu((WxmWMpoints(i)%xi-1)+1,WxmWMpoints(i)%yj,(WxmWMpoints(i)%zk)+1)+&
   nu((WxmWMpoints(i)%xi-1)+1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk)+&
   nu(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,(WxmWMpoints(i)%zk)+1)+&
   nu(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk))&
  *(W((WxmWMpoints(i)%xi-1)+1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk)-W(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk))&
  *recdxmin2/4

            wrk(WxmWMpoints(i)%xi,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk) = &
            wrk(WxmWMpoints(i)%xi,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk) + &
              WxmWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(WxmWMpoints)

            wrk(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk) = &
            wrk(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk) - &
  (nu((WxmWMpoints(i)%xi-1)+1,WxmWMpoints(i)%yj,(WxmWMpoints(i)%zk)+1)+&
   nu((WxmWMpoints(i)%xi-1)+1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk)+&
   nu(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,(WxmWMpoints(i)%zk)+1)+&
   nu(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk))&
  *(W((WxmWMpoints(i)%xi-1)+1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk)-W(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk))&
  *recdxmin2/4

            wrk(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk) = &
            wrk(WxmWMpoints(i)%xi-1,WxmWMpoints(i)%yj,WxmWMpoints(i)%zk) - &
              WxmWMpoints(i)%fluxp

      end do
      !$omp end do

      !$omp do
      do i = 1, size(WxpWMpoints)

            wrk(WxpWMpoints(i)%xi+1,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk) = &
            wrk(WxpWMpoints(i)%xi+1,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk) + &
  (nu((WxpWMpoints(i)%xi)+1,WxpWMpoints(i)%yj,(WxpWMpoints(i)%zk)+1)+&
   nu((WxpWMpoints(i)%xi)+1,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk)+&
   nu(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,(WxpWMpoints(i)%zk)+1)+&
   nu(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk))&
  *(W((WxpWMpoints(i)%xi)+1,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk)-W(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk))&
  *recdxmin2/4

            wrk(WxpWMpoints(i)%xi+1,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk) = &
            wrk(WxpWMpoints(i)%xi+1,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk) + &
              WxpWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(WxpWMpoints)

            wrk(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk) = &
            wrk(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk) - &
  (nu((WxpWMpoints(i)%xi)+1,WxpWMpoints(i)%yj,(WxpWMpoints(i)%zk)+1)+&
   nu((WxpWMpoints(i)%xi)+1,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk)+&
   nu(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,(WxpWMpoints(i)%zk)+1)+&
   nu(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk))&
  *(W((WxpWMpoints(i)%xi)+1,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk)-W(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk))&
  *recdxmin2/4

            wrk(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk) = &
            wrk(WxpWMpoints(i)%xi,WxpWMpoints(i)%yj,WxpWMpoints(i)%zk) - &
              WxpWMpoints(i)%fluxp

      end do
      !$omp end do

      !$omp do
      do i = 1, size(WymWMpoints)

            wrk(WymWMpoints(i)%xi,WymWMpoints(i)%yj,WymWMpoints(i)%zk) = &
            wrk(WymWMpoints(i)%xi,WymWMpoints(i)%yj,WymWMpoints(i)%zk) + &
  (nu(WymWMpoints(i)%xi,(WymWMpoints(i)%yj-1)+1,(WymWMpoints(i)%zk)+1)+&
   nu(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,(WymWMpoints(i)%zk)+1)+&
   nu(WymWMpoints(i)%xi,(WymWMpoints(i)%yj-1)+1,WymWMpoints(i)%zk)+&
   nu(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,WymWMpoints(i)%zk))&
  *(W(WymWMpoints(i)%xi,(WymWMpoints(i)%yj-1)+1,WymWMpoints(i)%zk)-W(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,WymWMpoints(i)%zk))&
  *recdymin2/4

            wrk(WymWMpoints(i)%xi,WymWMpoints(i)%yj,WymWMpoints(i)%zk) = &
            wrk(WymWMpoints(i)%xi,WymWMpoints(i)%yj,WymWMpoints(i)%zk) + &
              WymWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(WymWMpoints)

            wrk(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,WymWMpoints(i)%zk) = &
            wrk(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,WymWMpoints(i)%zk) - &
  (nu(WymWMpoints(i)%xi,(WymWMpoints(i)%yj-1)+1,(WymWMpoints(i)%zk)+1)+&
   nu(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,(WymWMpoints(i)%zk)+1)+&
   nu(WymWMpoints(i)%xi,(WymWMpoints(i)%yj-1)+1,WymWMpoints(i)%zk)+&
   nu(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,WymWMpoints(i)%zk))&
  *(W(WymWMpoints(i)%xi,(WymWMpoints(i)%yj-1)+1,WymWMpoints(i)%zk)-W(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,WymWMpoints(i)%zk))&
  *recdymin2/4

            wrk(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,WymWMpoints(i)%zk) = &
            wrk(WymWMpoints(i)%xi,WymWMpoints(i)%yj-1,WymWMpoints(i)%zk) - &
              WymWMpoints(i)%fluxp

      end do
      !$omp end do

      !$omp do
      do i = 1, size(WypWMpoints)

            wrk(WypWMpoints(i)%xi,WypWMpoints(i)%yj+1,WypWMpoints(i)%zk) = &
            wrk(WypWMpoints(i)%xi,WypWMpoints(i)%yj+1,WypWMpoints(i)%zk) + &
  (nu(WypWMpoints(i)%xi,(WypWMpoints(i)%yj)+1,(WypWMpoints(i)%zk)+1)+&
   nu(WypWMpoints(i)%xi,WypWMpoints(i)%yj,(WypWMpoints(i)%zk)+1)+&
   nu(WypWMpoints(i)%xi,(WypWMpoints(i)%yj)+1,WypWMpoints(i)%zk)+&
   nu(WypWMpoints(i)%xi,WypWMpoints(i)%yj,WypWMpoints(i)%zk))&
  *(W(WypWMpoints(i)%xi,(WypWMpoints(i)%yj)+1,WypWMpoints(i)%zk)-W(WypWMpoints(i)%xi,WypWMpoints(i)%yj,WypWMpoints(i)%zk))&
  *recdymin2/4

            wrk(WypWMpoints(i)%xi,WypWMpoints(i)%yj+1,WypWMpoints(i)%zk) = &
            wrk(WypWMpoints(i)%xi,WypWMpoints(i)%yj+1,WypWMpoints(i)%zk) + &
              WypWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(WypWMpoints)

            wrk(WypWMpoints(i)%xi,WypWMpoints(i)%yj,WypWMpoints(i)%zk) = &
            wrk(WypWMpoints(i)%xi,WypWMpoints(i)%yj,WypWMpoints(i)%zk) - &
  (nu(WypWMpoints(i)%xi,(WypWMpoints(i)%yj)+1,(WypWMpoints(i)%zk)+1)+&
   nu(WypWMpoints(i)%xi,WypWMpoints(i)%yj,(WypWMpoints(i)%zk)+1)+&
   nu(WypWMpoints(i)%xi,(WypWMpoints(i)%yj)+1,WypWMpoints(i)%zk)+&
   nu(WypWMpoints(i)%xi,WypWMpoints(i)%yj,WypWMpoints(i)%zk))&
  *(W(WypWMpoints(i)%xi,(WypWMpoints(i)%yj)+1,WypWMpoints(i)%zk)-W(WypWMpoints(i)%xi,WypWMpoints(i)%yj,WypWMpoints(i)%zk))&
  *recdymin2/4

            wrk(WypWMpoints(i)%xi,WypWMpoints(i)%yj,WypWMpoints(i)%zk) = &
            wrk(WypWMpoints(i)%xi,WypWMpoints(i)%yj,WypWMpoints(i)%zk) - &
              WypWMpoints(i)%fluxp

      end do
      !$omp end do

      !$omp do
      do i = 1, size(WzmWMpoints)

            wrk(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk) = &
            wrk(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk) + &
  nu(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,(WzmWMpoints(i)%zk-1)+1)*&
  (W(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,(WzmWMpoints(i)%zk-1)+1)-W(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk-1))&
  *recdzmin2

            wrk(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk) = &
            wrk(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk) + &
              WzmWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(WzmWMpoints)

            wrk(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk-1) = &
            wrk(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk-1) - &
  nu(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,(WzmWMpoints(i)%zk-1)+1)*&
  (W(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,(WzmWMpoints(i)%zk-1)+1)-W(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk-1))&
  *recdzmin2

            wrk(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk-1) = &
            wrk(WzmWMpoints(i)%xi,WzmWMpoints(i)%yj,WzmWMpoints(i)%zk-1) - &
              WzmWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(WzpWMpoints)

            wrk(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk+1) = &
            wrk(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk+1) + &
  nu(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,(WzpWMpoints(i)%zk)+1)*&
  (W(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,(WzpWMpoints(i)%zk)+1)-W(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk))&
  *recdzmin2

            wrk(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk+1) = &
            wrk(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk+1) + &
              WzpWMpoints(i)%fluxp

      end do
      !$omp end do
      !$omp do
      do i = 1, size(WzpWMpoints)

            wrk(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk) = &
            wrk(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk) - &
  nu(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,(WzpWMpoints(i)%zk)+1)*&
  (W(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,(WzpWMpoints(i)%zk)+1)-W(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk))&
  *recdzmin2

            wrk(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk) = &
            wrk(WzpWMpoints(i)%xi,WzpWMpoints(i)%yj,WzpWMpoints(i)%zk) - &
              WzpWMpoints(i)%fluxp

      end do
      !$omp end do
      
