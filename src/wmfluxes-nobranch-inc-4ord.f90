#if comp == 1

#define wrk U2

#if dir==1
#define xm UxmWMpoints
#define xp UxpWMpoints
#elif dir==2
#define ym UymWMpoints 
#define yp UypWMpoints
#else
#define zm UzmWMpoints 
#define zp UzpWMpoints
#endif

#elif comp == 2

#define wrk V2

#if dir==1
#define xm VxmWMpoints 
#define xp VxpWMpoints 
#elif dir==2
#define ym VymWMpoints 
#define yp VypWMpoints 
#else
#define zm VzmWMpoints 
#define zp VzpWMpoints 
#endif

#else

#define wrk W2

#if dir==1
#define xm WxmWMpoints 
#define xp WxpWMpoints 
#elif dir==2
#define ym WymWMpoints 
#define yp WypWMpoints 
#else
#define zm WzmWMpoints 
#define zp WzpWMpoints 
#endif

#endif



      !deconvolution of flux
      ! [Fl]i = -1/24*Fl(i) + 13*Fl(i)/12 -1/24*Fl(i)

      
      
#if dir==1



      !$omp do
      do ip = 1, size(xm)
#define p xm(ip)
        i = p%xi
        j = p%yj
        k = p%zk

        wrk(i,j,k) = &
          wrk(i,j,k) + (D1*Fl(i-1,j,k) + D0*Fl(i,j,k) + D1*Fl(i+1,j,k))

        wrk(i,j,k) = &
          wrk(i,j,k) + p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do ip = 1, size(xm)
#define p xm(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i-1,j,k) = &
          wrk(i-1,j,k) - (D1*Fl(i-1,j,k) + D0*Fl(i,j,k) + D1*Fl(i+1,j,k))

        wrk(i-1,j,k) = &
          wrk(i-1,j,k) - p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do ip = 1, size(xp)
#define p xp(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i+1,j,k) = &
          wrk(i+1,j,k) + (D1*Fl(i,j,k) + D0*Fl(i+1,j,k) + D1*Fl(i+2,j,k))

        wrk(i+1,j,k) = &
          wrk(i+1,j,k) + p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do ip = 1, size(xp)
#define p xp(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i,j,k) = &
          wrk(i,j,k) - (D1*Fl(i,j,k) + D0*Fl(i+1,j,k) + D1*Fl(i+2,j,k))

        wrk(i,j,k) = &
          wrk(i,j,k) - p%fluxp
#undef p
      end do
      !$omp end do
      
      
      
#elif dir==2

      

      !$omp do
      do ip = 1, size(ym)
#define p ym(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i,j,k) = &
          wrk(i,j,k) + (D1*Fl(i,j-1,k) + D0*Fl(i,j,k) + D1*Fl(i,j+1,k))

        wrk(i,j,k) = &
          wrk(i,j,k) + p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do ip = 1, size(ym)
#define p ym(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i,j-1,k) = &
          wrk(i,j-1,k) - (D1*Fl(i,j-1,k) + D0*Fl(i,j,k) + D1*Fl(i,j+1,k))

        wrk(i,j-1,k) = &
          wrk(i,j-1,k) - p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do ip = 1, size(yp)
#define p yp(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i,j+1,k) = &
          wrk(i,j+1,k) + (D1*Fl(i,j,k) + D0*Fl(i,j+1,k) + D1*Fl(i,j+2,k))

        wrk(i,j+1,k) = &
          wrk(i,j+1,k) + p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do ip = 1, size(yp)
#define p yp(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i,j,k) = &
          wrk(i,j,k) - (D1*Fl(i,j,k) + D0*Fl(i,j+1,k) + D1*Fl(i,j+2,k))

        wrk(i,j,k) = &
          wrk(i,j,k) - p%fluxp
#undef p
      end do
      !$omp end do
      
      
!dir==3  
#else

      

      !$omp do
      do ip = 1, size(zm)
#define p zm(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i,j,k) = &
          wrk(i,j,k) + (D1*Fl(i,j,k-1) + D0*Fl(i,j,k) + D1*Fl(i,j,k+1))

        wrk(i,j,k) = &
          wrk(i,j,k) + p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do ip = 1, size(zm)
#define p zm(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i,j,k-1) = &
          wrk(i,j,k-1) - (D1*Fl(i,j,k-1) + D0*Fl(i,j,k) + D1*Fl(i,j,k+1))

        wrk(i,j,k-1) = &
          wrk(i,j,k-1) - p%fluxp
#undef p
      end do
      !$omp end do

      !$omp do
      do ip = 1, size(zp)
#define p zp(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i,j,k+1) = &
          wrk(i,j,k+1) + (D1*Fl(i,j,k) + D0*Fl(i,j,k+1) + D1*Fl(i,j,k+2))

        wrk(i,j,k+1) = &
          wrk(i,j,k+1) + p%fluxp
#undef p
      end do
      !$omp end do
      !$omp do
      do ip = 1, size(zp)
#define p zp(ip)
        i = p%xi
        j = p%yj
        k = p%zk
        
        wrk(i,j,k) = &
          wrk(i,j,k) - (D1*Fl(i,j,k) + D0*Fl(i,j,k+1) + D1*Fl(i,j,k+2))

        wrk(i,j,k) = &
          wrk(i,j,k) - p%fluxp
#undef p
      end do
      !$omp end do
      
      
      
#endif
!dir==x,y,z



      
#undef flx
#undef fly
#undef flz
#undef xm
#undef xp
#undef ym
#undef yp
#undef zm
#undef zp

#undef wrk
