!These macros generate lines too long for Fortran.
!This file is therefore used only to generate the 
!files specific to each component which is than adjusted by hand 
!to havelines short enough.
!Command to generate specific files:
!    gfortran -E -cpp -Dcomp=1  wmfluxes-nobranch-inc.f90 > wmfluxes-nobranch-U-inc.f90

#if comp == 1

#define xm UxmWMpoints 
#define xp UxpWMpoints 
#define ym UymWMpoints 
#define yp UypWMpoints 
#define zm UzmWMpoints 
#define zp UzpWMpoints

#define flx(i,j,k) nu((i)+1,j,k) * (U((i)+1,j,k)-U(i,j,k)) * recdxmin2
#define fly(i,j,k) (nu((i)+1,(j)+1,k)+nu((i)+1,j,k)+nu(i,(j)+1,k)+nu(i,j,k))*(U(i,(j)+1,k)-U(i,j,k))*recdymin2/4
#define flz(i,j,k) (nu((i)+1,j,(k)+1)+nu((i)+1,j,k)+nu(i,j,(k)+1)+nu(i,j,k))*(U(i,j,(k)+1)-U(i,j,k))*recdzmin2/4

#elif comp == 2

#define xm VxmWMpoints 
#define xp VxpWMpoints 
#define ym VymWMpoints 
#define yp VypWMpoints 
#define zm VzmWMpoints 
#define zp VzpWMpoints 

#define flx(i,j,k) (nu((i)+1,(j)+1,k)+nu((i)+1,j,k)+nu(i,(j)+1,k)+nu(i,j,k))*(V((i)+1,j,k)-V(i,j,k))*recdxmin2/4
#define fly(i,j,k) nu(i,(j)+1,k) * (V(i,(j)+1,k)-V(i,j,k)) * recdymin2
#define flz(i,j,k) (nu(i,(j)+1,(k)+1)+nu(i,(j)+1,k)+nu(i,j,(k)+1)+nu(i,j,k))*(V(i,j,(k)+1)-V(i,j,k))*recdzmin2/4

#else

#define xm WxmWMpoints 
#define xp WxpWMpoints 
#define ym WymWMpoints 
#define yp WypWMpoints 
#define zm WzmWMpoints 
#define zp WzpWMpoints 

#define flx(i,j,k) (nu((i)+1,j,(k)+1)+nu((i)+1,j,k)+nu(i,j,(k)+1)+nu(i,j,k))*(W((i)+1,j,k)-W(i,j,k))*recdxmin2/4
#define fly(i,j,k) (nu(i,(j)+1,(k)+1)+nu(i,j,(k)+1)+nu(i,(j)+1,k)+nu(i,j,k))*(W(i,(j)+1,k)-W(i,j,k))*recdymin2/4
#define flz(i,j,k) nu(i,j,(k)+1)*(W(i,j,(k)+1)-W(i,j,k))*recdzmin2

#endif

      !$omp do
      do i = 1, size(xm)
#define p xm(i)
            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) + &
  flx(p%xi-1,p%yj,p%zk)

            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) + &
              p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do i = 1, size(xm)
#define p xm(i)
            wrk(p%xi-1,p%yj,p%zk) = &
            wrk(p%xi-1,p%yj,p%zk) - &
  flx(p%xi-1,p%yj,p%zk)

            wrk(p%xi-1,p%yj,p%zk) = &
            wrk(p%xi-1,p%yj,p%zk) - &
              p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do i = 1, size(xp)
#define p xp(i)
            wrk(p%xi+1,p%yj,p%zk) = &
            wrk(p%xi+1,p%yj,p%zk) + &
  flx(p%xi,p%yj,p%zk)

            wrk(p%xi+1,p%yj,p%zk) = &
            wrk(p%xi+1,p%yj,p%zk) + &
              p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do i = 1, size(xp)
#define p xp(i)
            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) - &
  flx(p%xi,p%yj,p%zk)

            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) - &
              p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do i = 1, size(ym)
#define p ym(i)
            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) + &
  fly(p%xi,p%yj-1,p%zk)

            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) + &
              p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do i = 1, size(ym)
#define p ym(i)
            wrk(p%xi,p%yj-1,p%zk) = &
            wrk(p%xi,p%yj-1,p%zk) - &
  fly(p%xi,p%yj-1,p%zk)

            wrk(p%xi,p%yj-1,p%zk) = &
            wrk(p%xi,p%yj-1,p%zk) - &
              p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do i = 1, size(yp)
#define p yp(i)
            wrk(p%xi,p%yj+1,p%zk) = &
            wrk(p%xi,p%yj+1,p%zk) + &
  fly(p%xi,p%yj,p%zk)

            wrk(p%xi,p%yj+1,p%zk) = &
            wrk(p%xi,p%yj+1,p%zk) + &
              p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do i = 1, size(yp)
#define p yp(i)
            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) - &
  fly(p%xi,p%yj,p%zk)

            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) - &
              p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do i = 1, size(zm)
#define p zm(i)
            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) + &
  flz(p%xi,p%yj,p%zk-1)

            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) + &
              p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do i = 1, size(zm)
#define p zm(i)
            wrk(p%xi,p%yj,p%zk-1) = &
            wrk(p%xi,p%yj,p%zk-1) - &
  flz(p%xi,p%yj,p%zk-1)

            wrk(p%xi,p%yj,p%zk-1) = &
            wrk(p%xi,p%yj,p%zk-1) - &
              p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do i = 1, size(zp)
#define p zp(i)
            wrk(p%xi,p%yj,p%zk+1) = &
            wrk(p%xi,p%yj,p%zk+1) + &
  flz(p%xi,p%yj,p%zk)

            wrk(p%xi,p%yj,p%zk+1) = &
            wrk(p%xi,p%yj,p%zk+1) + &
              p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do i = 1, size(zp)
#define p zp(i)
            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) - &
  flz(p%xi,p%yj,p%zk)

            wrk(p%xi,p%yj,p%zk) = &
            wrk(p%xi,p%yj,p%zk) - &
              p%fluxp
#undef p
      end do
      !$omp end do
      
#undef flx
#undef fly
#undef flz
#undef xm
#undef xp
#undef ym
#undef yp
#undef zm
#undef zp
