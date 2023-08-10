#if comp == 1

#define xm UxmWMpoints 
#define xp UxpWMpoints 
#define ym UymWMpoints 
#define yp UypWMpoints 
#define zm UzmWMpoints 
#define zp UzpWMpoints 

#elif comp == 2

#define xm VxmWMpoints 
#define xp VxpWMpoints 
#define ym VymWMpoints 
#define yp VypWMpoints 
#define zm VzmWMpoints 
#define zp VzpWMpoints 

#else

#define xm WxmWMpoints 
#define xp WxpWMpoints 
#define ym WymWMpoints 
#define yp WypWMpoints 
#define zm WzmWMpoints 
#define zp WzpWMpoints 

#endif

      !$omp do
      do i = 1, size(xm)
#define p xm(i)
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
              p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do i = 1, size(xp)
#define p xp(i)
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
              p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do i = 1, size(ym)
#define p ym(i)
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
              p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do i = 1, size(yp)
#define p yp(i)
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
              p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do i = 1, size(zm)
#define p zm(i)
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
              p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do i = 1, size(zp)
#define p zp(i)
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
              p%fluxp
#undef p
      end do
      !$omp end do
      
#undef xm
#undef xp
#undef ym
#undef yp
#undef zm
#undef zp
