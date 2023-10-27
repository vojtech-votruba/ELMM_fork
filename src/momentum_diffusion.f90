module MomentumDiffusion
  use Parameters
  use ArrayUtilities
  use Tiling, only: tilenx, tileny, tilenz
  use Boundaries, only: BoundUVW

  implicit none
  
  
  real(knd), dimension(:,:,:), allocatable :: Apu, ApV, ApW, work

contains




  subroutine MomentumDiffusion_2ord(U2,V2,W2,U,V,W)
    use Parameters, nu => Viscosity
    use Wallmodels
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2,V2,W2
    real(knd) :: recdxmin2, recdymin2, recdzmin2
    integer :: i,j,k,bi,bj,bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 6
       
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2

    !$omp parallel private(i,j,k,bi,bj,bk)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           if (Uflx_mask(i+1,j,k)) &
             U2(i,j,k) = U2(i,j,k) + &
              nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) *recdxmin2
           if (Uflx_mask(i,j,k)) &
             U2(i,j,k) = U2(i,j,k) - &
               nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
           if (Ufly_mask(i,j+1,k)) &
             U2(i,j,k) = U2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
           if (Ufly_mask(i,j,k)) &
             U2(i,j,k) = U2(i,j,k) - &
               0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
           if (Uflz_mask(i,j,k+1)) &
             U2(i,j,k) = U2(i,j,k) + &
               0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
           if (Uflz_mask(i,j,k)) &
             U2(i,j,k) = U2(i,j,k) - &
               0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 1
#define wrk U2
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           if (Vflx_mask(i+1,j,k)) &
             V2(i,j,k) = V2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
           if (Vflx_mask(i,j,k)) &
             V2(i,j,k) = V2(i,j,k) - &
               0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
           if (Vfly_mask(i,j+1,k)) &
             V2(i,j,k) = V2(i,j,k) + &
               nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
           if (Vfly_mask(i,j,k)) &
             V2(i,j,k) = V2(i,j,k) - &
               nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
           if (Vflz_mask(i,j,k+1)) &
             V2(i,j,k) = V2(i,j,k) + &
               0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
           if (Vflz_mask(i,j,k)) &
             V2(i,j,k) = V2(i,j,k) - &
               0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 2
#define wrk V2
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           if (Wflx_mask(i+1,j,k)) &
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
           if (Wflx_mask(i,j,k)) &
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
           if (Wfly_mask(i,j+1,k)) &
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
           if (Wfly_mask(i,j,k)) &
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
           if (Wflz_mask(i,j,k+1)) &
             W2(i,j,k) = W2(i,j,k) + &
               nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
           if (Wflz_mask(i,j,k)) &
             W2(i,j,k) = W2(i,j,k) - &
               nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 3
#define wrk W2
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
    !$omp end parallel

  end subroutine MomentumDiffusion_2ord














  subroutine MomentumDiffusion_nobranch_2ord(U2,V2,W2,U,V,W)
    use Parameters, nu => Viscosity
    use Wallmodels
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2,V2,W2
    real(knd) :: recdxmin2, recdymin2, recdzmin2
    integer :: i,j,k,bi,bj,bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 3
       
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)
       
    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2

    !$omp parallel private(i,j,k,bi,bj,bk)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
             U2(i,j,k) = U2(i,j,k) + &
              nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) * recdxmin2
             U2(i,j,k) = U2(i,j,k) - &
               nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
             U2(i,j,k) = U2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
             U2(i,j,k) = U2(i,j,k) - &
               0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
             U2(i,j,k) = U2(i,j,k) + &
               0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
             U2(i,j,k) = U2(i,j,k) - &
               0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 1
#define wrk U2
#include "wmfluxes-nobranch-U-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
             V2(i,j,k) = V2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
             V2(i,j,k) = V2(i,j,k) - &
               0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
             V2(i,j,k) = V2(i,j,k) + &
               nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
             V2(i,j,k) = V2(i,j,k) - &
               nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
             V2(i,j,k) = V2(i,j,k) + &
               0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
             V2(i,j,k) = V2(i,j,k) - &
               0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 2
#define wrk V2
#include "wmfluxes-nobranch-V-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
             W2(i,j,k) = W2(i,j,k) + &
               nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
             W2(i,j,k) = W2(i,j,k) - &
               nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 3
#define wrk W2
#include "wmfluxes-nobranch-W-inc.f90"
#undef wrk
#undef comp
    !$omp end parallel

  end subroutine MomentumDiffusion_nobranch_2ord









  subroutine MomentumDiffusion_variable_z_nobranch_2ord(U2,V2,W2,U,V,W)
    use Parameters, nu => Viscosity
    use Wallmodels
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2,V2,W2
    real(knd) :: recdxmin2, recdymin2, tmp
    integer :: i,j,k,bi,bj,bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 3
       
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)
       
    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2

    !$omp parallel private(i,j,k,bi,bj,bk,tmp)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
             U2(i,j,k) = U2(i,j,k) + &
              nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) * recdxmin2
             U2(i,j,k) = U2(i,j,k) - &
               nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
             U2(i,j,k) = U2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
             U2(i,j,k) = U2(i,j,k) - &
               0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
             tmp = 0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * &
                   (U(i,j,k+1)-U(i,j,k)) / dzW(k)
             tmp = tmp - 0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * &
                   (U(i,j,k)-U(i,j,k-1)) / dzW(k-1)
             tmp = tmp / dzPr(k)
             U2(i,j,k) = U2(i,j,k) + tmp
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 1
#define wrk U2
#include "wmfluxes-variable_z-nobranch-U-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
             V2(i,j,k) = V2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
             V2(i,j,k) = V2(i,j,k) - &
               0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
             V2(i,j,k) = V2(i,j,k) + &
               nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
             V2(i,j,k) = V2(i,j,k) - &
               nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
             tmp = 0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * &
                   (V(i,j,k+1)-V(i,j,k)) / dzW(k)
             tmp = tmp - 0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * &
                   (V(i,j,k)-V(i,j,k-1)) / dzW(k-1)
             tmp = tmp / dzPr(k)
             V2(i,j,k) = V2(i,j,k) + tmp
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 2
#define wrk V2
#include "wmfluxes-variable_z-nobranch-V-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
             tmp = nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) / dzPr(k+1)
             tmp = tmp - nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) / dzPr(k)
             tmp = tmp / dzW(k)
             W2(i,j,k) = W2(i,j,k) + tmp
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 3
#define wrk W2
#include "wmfluxes-variable_z-nobranch-W-inc.f90"
#undef wrk
#undef comp
    !$omp end parallel

  end subroutine MomentumDiffusion_variable_z_nobranch_2ord









  subroutine MomentumDiffusion_nobranch_4ord(U2,V2,W2,U,V,W)
    use Parameters, nu => Viscosity
    use Wallmodels
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2,V2,W2
    real(knd), dimension(:,:,:), allocatable :: Fl
    real(knd) :: recdxmin2, recdymin2, recdzmin2
    integer :: i,j,k,bi,bj,bk
    integer :: tnx, tny, tnz
    
    real(knd), parameter :: C1 = 9._knd / 8, C3 = 1._knd / (8*3)
    real(knd), parameter :: D0 = 13._knd / 12, D1 = -1._knd / 24

    integer, parameter :: narr = 3
    
    !for the include files
    integer :: ip
       
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)
       
    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2
    
    if (.not.allocated(Fl)) &
      allocate(Fl(-1:max(Unx,Vnx,Wnx)+2, &
                  -1:max(Uny,Vny,Wny)+2, &
                  -1:max(Unz,Vnz,Wnz)+2))

    !$omp parallel private(i,j,k,bi,bj,bk)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 0, Unx+2, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx+2)
           Fl(i,j,k) = nu(i,j,k) * (C1*(U(i,j,k)-U(i-1,j,k)) - C3*(U(i+1,j,k)-U(i-2,j,k))) * recdxmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           U2(i,j,k) = U2(i,j,k) + (C1*(Fl(i+1,j,k)-Fl(i,j,k)) - C3*(Fl(i+2,j,k)-Fl(i-1,j,k)))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    
#define comp 1
#define dir 1
#include "wmfluxes-nobranch-inc-4ord.f90"
#undef dir
#undef comp
    
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 0, Uny+2, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny+2)
         do i = bi, min(bi+tnx-1,Unx)
           Fl(i,j,k) = 0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * &
               (C1*(U(i,j,k)-U(i,j-1,k)) - C3*(U(i,j+1,k)-U(i,j-2,k))) * recdymin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           U2(i,j,k) = U2(i,j,k) + (C1*(Fl(i,j+1,k)-Fl(i,j,k)) - C3*(Fl(i,j+2,k)-Fl(i,j-1,k)))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
               
#define comp 1
#define dir 2
#include "wmfluxes-nobranch-inc-4ord.f90"
#undef dir
#undef comp
                   
    !$omp do schedule(runtime) collapse(3)
    do bk = 0, Unz+2, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz+2)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           Fl(i,j,k) = 0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * &
               (C1*(U(i,j,k)-U(i,j,k-1)) - C3*(U(i,j,k+1)-U(i,j,k-2))) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           U2(i,j,k) = U2(i,j,k) + (C1*(Fl(i,j,k+1)-Fl(i,j,k)) - C3*(Fl(i,j,k+2)-Fl(i,j,k-1)))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do

#define comp 1
#define dir 3
#include "wmfluxes-nobranch-inc-4ord.f90"
#undef dir
#undef comp
        
    
    
    

    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 0, Vnx+2, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx+2)
           Fl(i,j,k) = 0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * &
               (C1*(V(i,j,k)-V(i-1,j,k)) - C3*(V(i+1,j,k)-V(i-2,j,k))) * recdxmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           V2(i,j,k) = V2(i,j,k) + (C1*(Fl(i+1,j,k)-Fl(i,j,k)) - C3*(Fl(i+2,j,k)-Fl(i-1,j,k))) 
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do

#define comp 2
#define dir 1
#include "wmfluxes-nobranch-inc-4ord.f90"
#undef dir
#undef comp
                       
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 0, Vny+2, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny+2)
         do i = bi, min(bi+tnx-1,Vnx)
           Fl(i,j,k) = nu(i,j,k) * (C1*(V(i,j,k)-V(i,j-1,k)) - C3*(V(i,j+1,k)-V(i,j-2,k))) * recdymin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           V2(i,j,k) = V2(i,j,k) + (C1*(Fl(i,j+1,k)-Fl(i,j,k)) - C3*(Fl(i,j+2,k)-Fl(i,j-1,k)))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 2
#define dir 2
#include "wmfluxes-nobranch-inc-4ord.f90"
#undef dir
#undef comp
                                  
    !$omp do schedule(runtime) collapse(3)
    do bk = 0, Vnz+2, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz+2)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           Fl(i,j,k) = 0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * &
               (C1*(V(i,j,k)-V(i,j,k-1)) - C3*(V(i,j,k+1)-V(i,j,k-2))) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           V2(i,j,k) = V2(i,j,k) + (C1*(Fl(i,j,k+1)-Fl(i,j,k)) - C3*(Fl(i,j,k+2)-Fl(i,j,k-1)))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 2
#define dir 3
#include "wmfluxes-nobranch-inc-4ord.f90"
#undef dir
#undef comp
                                  

    
    
    


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 0, Wnx+2, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx+2)
           Fl(i,j,k) = 0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * &
               (C1*(W(i,j,k)-W(i-1,j,k)) - C3*(W(i+1,j,k)-W(i-2,j,k))) * recdxmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           W2(i,j,k) = W2(i,j,k) + (C1*(Fl(i+1,j,k)-Fl(i,j,k)) - C3*(Fl(i+2,j,k)-Fl(i-1,j,k)))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 3
#define dir 1
#include "wmfluxes-nobranch-inc-4ord.f90"
#undef dir
#undef comp
                                  
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 0, Wny+2, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny+2)
         do i = bi, min(bi+tnx-1,Wnx)
           Fl(i,j,k) = 0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * &
               (C1*(W(i,j,k)-W(i,j-1,k)) - C3*(W(i,j+1,k)-W(i,j-2,k))) * recdymin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           W2(i,j,k) = W2(i,j,k) + (C1*(Fl(i,j+1,k)-Fl(i,j,k)) - C3*(Fl(i,j+2,k)-Fl(i,j-1,k)))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 3
#define dir 2
#include "wmfluxes-nobranch-inc-4ord.f90"
#undef dir
#undef comp
                                  
    !$omp do schedule(runtime) collapse(3)
    do bk = 0, Wnz+2, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz+2)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           Fl(i,j,k) = nu(i,j,k) * (C1*(W(i,j,k)-W(i,j,k-1)) - C3*(W(i,j,k+1)-W(i,j,k-2))) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           W2(i,j,k) = W2(i,j,k) + (C1*(Fl(i,j,k+1)-Fl(i,j,k)) - C3*(Fl(i,j,k+2)-Fl(i,j,k-1)))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 3
#define dir 3
#include "wmfluxes-nobranch-inc-4ord.f90"
#undef dir
#undef comp



    !$omp end parallel

  end subroutine MomentumDiffusion_nobranch_4ord











  subroutine MomentumDiffusion_4ord_5point(U2,V2,W2,U,V,W)
    use Parameters, nu => Viscosity
    use Wallmodels
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2,V2,W2
    real(knd), dimension(:,:,:), allocatable :: Fl
    real(knd) :: recdxmin2, recdymin2, recdzmin2
    integer :: i,j,k,bi,bj,bk
    integer :: tnx, tny, tnz
    
    real(knd), parameter :: C1 = 15._knd / 12, C3 = 1._knd / 12

    integer, parameter :: narr = 3
    
    !for the include files
    integer :: ip
    
    ! Discretization leading to the 5-point scheme for the second derivative, 
    ! split into a flux divergence.
    ! F(i) = (15 (U(i)-U(i-1)) - 1 (U(i+1) - U(i-2))) / (12 dx)
    ! dU = (F(i+1) - F(i)) / dx
       
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)
       
    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2
    
    if (.not.allocated(Fl)) &
      allocate(Fl(-1:max(Unx,Vnx,Wnx)+2, &
                  -1:max(Uny,Vny,Wny)+2, &
                  -1:max(Unz,Vnz,Wnz)+2))

    !$omp parallel private(i,j,k,bi,bj,bk)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 0, Unx+2, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx+2)
           Fl(i,j,k) = nu(i,j,k) * (C1*(U(i,j,k)-U(i-1,j,k)) - C3*(U(i+1,j,k)-U(i-2,j,k))) * recdxmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           U2(i,j,k) = U2(i,j,k) + (Fl(i+1,j,k)-Fl(i,j,k))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    
#define comp 1
#define dir 1
#include "wmfluxes-nobranch-inc-4ord-5point.f90"
#undef dir
#undef comp
    
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 0, Uny+2, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny+2)
         do i = bi, min(bi+tnx-1,Unx)
           Fl(i,j,k) = 0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * &
               (C1*(U(i,j,k)-U(i,j-1,k)) - C3*(U(i,j+1,k)-U(i,j-2,k))) * recdymin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           U2(i,j,k) = U2(i,j,k) + (Fl(i,j+1,k)-Fl(i,j,k))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
               
#define comp 1
#define dir 2
#include "wmfluxes-nobranch-inc-4ord-5point.f90"
#undef dir
#undef comp
                   
    !$omp do schedule(runtime) collapse(3)
    do bk = 0, Unz+2, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz+2)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           Fl(i,j,k) = 0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * &
               (C1*(U(i,j,k)-U(i,j,k-1)) - C3*(U(i,j,k+1)-U(i,j,k-2))) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           U2(i,j,k) = U2(i,j,k) + (Fl(i,j,k+1)-Fl(i,j,k))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do

#define comp 1
#define dir 3
#include "wmfluxes-nobranch-inc-4ord-5point.f90"
#undef dir
#undef comp
        
    
    
    

    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 0, Vnx+2, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx+2)
           Fl(i,j,k) = 0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * &
               (C1*(V(i,j,k)-V(i-1,j,k)) - C3*(V(i+1,j,k)-V(i-2,j,k))) * recdxmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           V2(i,j,k) = V2(i,j,k) + (Fl(i+1,j,k)-Fl(i,j,k))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do

#define comp 2
#define dir 1
#include "wmfluxes-nobranch-inc-4ord-5point.f90"
#undef dir
#undef comp
                       
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 0, Vny+2, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny+2)
         do i = bi, min(bi+tnx-1,Vnx)
           Fl(i,j,k) = nu(i,j,k) * (C1*(V(i,j,k)-V(i,j-1,k)) - C3*(V(i,j+1,k)-V(i,j-2,k))) * recdymin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           V2(i,j,k) = V2(i,j,k) + (Fl(i,j+1,k)-Fl(i,j,k))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 2
#define dir 2
#include "wmfluxes-nobranch-inc-4ord-5point.f90"
#undef dir
#undef comp
                                  
    !$omp do schedule(runtime) collapse(3)
    do bk = 0, Vnz+2, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz+2)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           Fl(i,j,k) = 0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * &
               (C1*(V(i,j,k)-V(i,j,k-1)) - C3*(V(i,j,k+1)-V(i,j,k-2))) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           V2(i,j,k) = V2(i,j,k) + (Fl(i,j,k+1)-Fl(i,j,k))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 2
#define dir 3
#include "wmfluxes-nobranch-inc-4ord-5point.f90"
#undef dir
#undef comp
                                  

    
    
    


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 0, Wnx+2, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx+2)
           Fl(i,j,k) = 0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * &
               (C1*(W(i,j,k)-W(i-1,j,k)) - C3*(W(i+1,j,k)-W(i-2,j,k))) * recdxmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           W2(i,j,k) = W2(i,j,k) + (Fl(i+1,j,k)-Fl(i,j,k))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 3
#define dir 1
#include "wmfluxes-nobranch-inc-4ord-5point.f90"
#undef dir
#undef comp
                                  
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 0, Wny+2, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny+2)
         do i = bi, min(bi+tnx-1,Wnx)
           Fl(i,j,k) = 0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * &
               (C1*(W(i,j,k)-W(i,j-1,k)) - C3*(W(i,j+1,k)-W(i,j-2,k))) * recdymin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           W2(i,j,k) = W2(i,j,k) + (Fl(i,j+1,k)-Fl(i,j,k))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 3
#define dir 2
#include "wmfluxes-nobranch-inc-4ord-5point.f90"
#undef dir
#undef comp
                                  
    !$omp do schedule(runtime) collapse(3)
    do bk = 0, Wnz+2, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz+2)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           Fl(i,j,k) = nu(i,j,k) * (C1*(W(i,j,k)-W(i,j,k-1)) - C3*(W(i,j,k+1)-W(i,j,k-2))) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           W2(i,j,k) = W2(i,j,k) + (Fl(i,j,k+1)-Fl(i,j,k))
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
             
#define comp 3
#define dir 3
#include "wmfluxes-nobranch-inc-4ord-5point.f90"
#undef dir
#undef comp



    !$omp end parallel

  end subroutine MomentumDiffusion_4ord_5point



















! Remaining parts of implicit diffusion




  subroutine ImplicitDiffusion_ForwEul(U, V, W, U2, V2, W2, U3, V3, W3, coef)
    use Parameters, nu => Viscosity
    use Wallmodels
    real(knd), intent(in),  dimension(-2:,-2:,-2:), contiguous :: U,V,W
    real(knd), intent(in),  dimension(-2:,-2:,-2:), contiguous :: U2,V2,W2
    real(knd), intent(out), dimension(-2:,-2:,-2:), contiguous :: U3,V3,W3
    real(knd), intent(in) :: coef

    real(knd) :: Ap, recdxmin2, recdymin2, recdzmin2
    integer   :: i,j,k


    Ap = coef

    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2

    !$omp parallel private(i,j,k)

    !$omp do
    do k = 1, Unz
      do j = 1, Uny
        do i = 1, Unx
            U3(i,j,k) = 0
            if (Uflx_mask(i+1,j,k)) &
              U3(i,j,k) = U3(i,j,k) + &
                nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) *recdxmin2
            if (Uflx_mask(i,j,k)) &
              U3(i,j,k) = U3(i,j,k) - &
                nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
            if (Ufly_mask(i,j+1,k)) &
              U3(i,j,k) = U3(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
            if (Ufly_mask(i,j,k)) &
              U3(i,j,k) = U3(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
            if (Uflz_mask(i,j,k+1)) &
              U3(i,j,k) = U3(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
            if (Uflz_mask(i,j,k)) &
              U3(i,j,k) = U3(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
        end do
      end do
    end do
    !$omp end do
#define comp 1
#define wrk U3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
    !$omp do
    do k = 1, Unz
      do j = 1, Uny
        do i = 1, Unx
          U3(i,j,k) = U3(i,j,k) * Ap
          U3(i,j,k) = U3(i,j,k) + U(i,j,k) + U2(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait


    !$omp do
    do k = 1, Vnz
      do j = 1, Vny
        do i = 1, Vnx
            V3(i,j,k) = 0
            if (Vflx_mask(i+1,j,k)) &
              V3(i,j,k) = V3(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
            if (Vflx_mask(i,j,k)) &
              V3(i,j,k) = V3(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
            if (Vfly_mask(i,j+1,k)) &
              V3(i,j,k) = V3(i,j,k) + &
                nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
            if (Vfly_mask(i,j,k)) &
              V3(i,j,k) = V3(i,j,k) - &
                nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
            if (Vflz_mask(i,j,k+1)) &
              V3(i,j,k) = V3(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
            if (Vflz_mask(i,j,k)) &
              V3(i,j,k) = V3(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
        end do
      end do
    end do
    !$omp end do
#define comp 2
#define wrk V3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
    !$omp do
    do k = 1, Vnz
      do j = 1, Vny
        do i = 1, Vnx
          V3(i,j,k) = V3(i,j,k) * Ap
          V3(i,j,k) = V3(i,j,k) + V(i,j,k) + V2(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait


    !$omp do
    do k = 1, Wnz
      do j = 1, Wny
        do i = 1, Wnx
            W3(i,j,k) = 0
            if (Wflx_mask(i+1,j,k)) &
              W3(i,j,k) = W3(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
            if (Wflx_mask(i,j,k)) &
              W3(i,j,k) = W3(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
            if (Wfly_mask(i,j+1,k)) &
              W3(i,j,k) = W3(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
            if (Wfly_mask(i,j,k)) &
              W3(i,j,k) = W3(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
            if (Wflz_mask(i,j,k+1)) &
              W3(i,j,k) = W3(i,j,k) + &
                nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
            if (Wflz_mask(i,j,k)) &
              W3(i,j,k) = W3(i,j,k) - &
                nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
        end do
      end do
    end do
    !$omp end do
#define comp 3
#define wrk W3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
    !$omp do
    do k = 1, Wnz
      do j = 1, Wny
        do i = 1, Wnx
          W3(i,j,k) = W3(i,j,k) * Ap
          W3(i,j,k) = W3(i,j,k) + W(i,j,k) + W2(i,j,k)
        end do
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine ImplicitDiffusion_ForwEul






  subroutine ImplicitDiffusion_Iterations(U, V, W, U2, V2, W2, U3, V3, W3, coef)
    use Parameters, nu => Viscosity
    use Wallmodels
    !$ use omp_lib
#ifdef PAR
    use custom_par, only: par_co_max
#endif
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2, V2, W2, U3, V3, W3
    real(knd), intent(in) :: coef
    real(knd) :: recdxmin2,recdymin2,recdzmin2                                                               !reciprocal values of dx**2
    real(knd) :: Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
    integer :: i,j,k,bi,bj,bk,l
    integer, save :: called = 0
    integer :: tnx, tny, tnz, tnx2, tny2, tnz2
    
    integer, parameter :: narr = 3, narr2 = 5
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    tnx2 = tilenx(narr2)
    tny2 = tileny(narr2)
    tnz2 = tilenz(narr2)


    if (called==0) then
      allocate(Apu(1:Unx,1:Uny,1:Unz))
      allocate(ApV(1:Vnx,1:Vny,1:Vnz))
      allocate(ApW(1:Wnx,1:Wny,1:Wnz))
      allocate(work(0:max(Unx+1,Vnx+1,Wnx+1), &
                  0:max(Uny+1,Vny+1,Wny+1), &
                  0:max(Unz+1,Uny+1,Wnz+1)))
      called = 1
    end if


    Ap = coef / 2
    S = 0
    l = 0

    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2

    !$omp parallel private(i,j,k,bi,bj,bk)

    !The explicit part, which doesn't have to be changed inside the loop
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
            work(i,j,k) = 0
            if (Uflx_mask(i+1,j,k)) &
              work(i,j,k) = work(i,j,k) + &
                nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) *recdxmin2
            if (Uflx_mask(i,j,k)) &
              work(i,j,k) = work(i,j,k) - &
                nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
            if (Ufly_mask(i,j+1,k)) &
              work(i,j,k) = work(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
            if (Ufly_mask(i,j,k)) &
              work(i,j,k) = work(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
            if (Uflz_mask(i,j,k+1)) &
              work(i,j,k) = work(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
            if (Uflz_mask(i,j,k)) &
              work(i,j,k) = work(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 1
#define wrk work
#include "wmfluxes-inc.f90"
#undef comp
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
              U2(i,j,k) = U2(i,j,k) + Ap * work(i,j,k)
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do nowait
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
            work(i,j,k) = 0
            if (Vflx_mask(i+1,j,k)) &
              work(i,j,k) = work(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
            if (Vflx_mask(i,j,k)) &
              work(i,j,k) = work(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
            if (Vfly_mask(i,j+1,k)) &
              work(i,j,k) = work(i,j,k) + &
                nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
            if (Vfly_mask(i,j,k)) &
              work(i,j,k) = work(i,j,k) - &
                nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
            if (Vflz_mask(i,j,k+1)) &
              work(i,j,k) = work(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
            if (Vflz_mask(i,j,k)) &
              work(i,j,k) = work(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 2
#include "wmfluxes-inc.f90"
#undef comp
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           V2(i,j,k) = V2(i,j,k) + Ap * work(i,j,k)
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
            work(i,j,k) = 0
            if (Wflx_mask(i+1,j,k)) &
              work(i,j,k) = work(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
            if (Wflx_mask(i,j,k)) &
              work(i,j,k) = work(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
            if (Wfly_mask(i,j+1,k)) &
              work(i,j,k) = work(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
            if (Wfly_mask(i,j,k)) &
              work(i,j,k) = work(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
            if (Wflz_mask(i,j,k+1)) &
              work(i,j,k) = work(i,j,k) + &
                nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
            if (Wflz_mask(i,j,k)) &
              work(i,j,k) = work(i,j,k) - &
                nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 3
#include "wmfluxes-inc.f90"
#undef comp
    !$omp do schedule(runtime)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           W2(i,j,k) = W2(i,j,k) + Ap * work(i,j,k)
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do

    !Auxiliary coefficients to better efficiency in loops
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
            ApU(i,j,k) = 0
            if (Uflx_mask(i+1,j,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                nu(i+1,j,k) * recdxmin2
            if (Uflx_mask(i,j,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                nu(i,j,k) * recdxmin2
            if (Ufly_mask(i,j+1,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * recdymin2
            if (Ufly_mask(i,j,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * recdymin2
            if (Uflz_mask(i,j,k+1)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * recdzmin2
            if (Uflz_mask(i,j,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp end parallel

!        ApU = 1._knd/(1._knd+Ap*ApU)
    call multiply_and_add_scalar(ApU, Ap, 1._knd)
    call reciprocal(ApU, 1._knd)

    !$omp parallel private(i,j,k,bi,bj,bk)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
            ApV(i,j,k) =0
            if (Vflx_mask(i+1,j,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * recdxmin2
            if (Vflx_mask(i,j,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * recdxmin2
            if (Vfly_mask(i,j+1,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                nu(i,j+1,k) * recdymin2
            if (Vfly_mask(i,j,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                nu(i,j,k) * recdymin2
            if (Vflz_mask(i,j,k+1)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * recdzmin2
            if (Vflz_mask(i,j,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp end parallel

!        ApV = 1._knd/(1._knd+Ap*ApV)
    call multiply_and_add_scalar(ApV, Ap, 1._knd)
    call reciprocal(ApV, 1._knd)

    !$omp parallel private(i,j,k,bi,bj,bk,Suavg,Svavg,Swavg)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
            ApW(i,j,k) = 0
            if (Wflx_mask(i+1,j,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * recdxmin2
            if (Wflx_mask(i,j,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * recdxmin2
            if (Wfly_mask(i,j+1,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * recdymin2
            if (Wfly_mask(i,j,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * recdymin2
            if (Wflz_mask(i,j,k+1)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                nu(i,j,k+1) * recdzmin2
            if (Wflz_mask(i,j,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                nu(i,j,k) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp end parallel


!        ApW = 1._knd/(1._knd+Ap*ApW)
    call multiply_and_add_scalar(ApW, Ap, 1._knd)
    call reciprocal(ApW, 1._knd)


    do l = 1, maxCNiter               !Gauss-Seidel iteration for Crank-Nicolson result
      call BoundUVW(U3, V3, W3)

      S = 0
      Su = 0
      Sv = 0
      Sw = 0
      !$omp parallel private(i,j,k,bi,bj,bk,p) reduction(max:Su,Sv,Sw)
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Unz, tnz2
       do bj = 1, Uny, tny2
        do bi = 1, Unx, tnx2
         do k = bk, min(bk+tnz2-1,Unz)
          do j = bj, min(bj+tny2-1,Uny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Unx), 2
            if (Utype(i,j,k)<=0) then
              work(i,j,k) = 0
              if (Uflx_mask(i+1,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i+1,j,k) * U3(i+1,j,k) * recdxmin2
              if (Uflx_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j,k) * U3(i-1,j,k) * recdxmin2
              if (Ufly_mask(i,j+1,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * U3(i,j+1,k) * recdymin2
              if (Ufly_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * U3(i,j-1,k) * recdymin2
              if (Uflz_mask(i,j,k+1)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * U3(i,j,k+1) * recdzmin2
              if (Uflz_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * U3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 1
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Unz, tnz2
       do bj = 1, Uny, tny2
        do bi = 1, Unx, tnx2
         do k = bk, min(bk+tnz2-1,Unz)
          do j = bj, min(bj+tny2-1,Uny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Unx), 2
            if (Utype(i,j,k)<=0) then
              p = Ap * work(i,j,k) + U2(i,j,k) + U(i,j,k)
              p = p * ApU(i,j,k)
              Su = max(Su,abs(p-U3(i,j,k)))
              U3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do

      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Vnz, tnz2
       do bj = 1, Vny, tny2
        do bi = 1, Vnx, tnx2
         do k = bk, min(bk+tnz2-1,Vnz)
          do j = bj, min(bj+tny2-1,Vny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Vnx), 2
            if (Vtype(i,j,k)<=0) then
              work(i,j,k) = 0
              if (Vflx_mask(i+1,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * V3(i+1,j,k) * recdxmin2
              if (Vflx_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * V3(i-1,j,k) * recdxmin2
              if (Vfly_mask(i,j+1,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j+1,k) * V3(i,j+1,k) * recdymin2
              if (Vfly_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j,k) * V3(i,j-1,k) * recdymin2
              if (Vflz_mask(i,j,k+1)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * V3(i,j,k+1) * recdzmin2
              if (Vflz_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * V3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 2
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Vnz, tnz2
       do bj = 1, Vny, tny2
        do bi = 1, Vnx, tnx2
         do k = bk, min(bk+tnz2-1,Vnz)
          do j = bj, min(bj+tny2-1,Vny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Vnx), 2
            if (Vtype(i,j,k)<=0) then
              p = Ap * work(i,j,k) + V2(i,j,k) + V(i,j,k)
              p = p * ApV(i,j,k)
              Sv = max(Sv,abs(p-V3(i,j,k)))
              V3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
         
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Wnz, tnz2
       do bj = 1, Wny, tny2
        do bi = 1, Wnx, tnx2
         do k = bk, min(bk+tnz2-1,Wnz)
          do j = bj, min(bj+tny2-1,Wny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Wnx), 2
            if (Wtype(i,j,k)<=0) then
              work(i,j,k) = 0
              if (Wflx_mask(i+1,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * W3(i+1,j,k) * recdxmin2
              if (Wflx_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * W3(i-1,j,k) * recdxmin2
              if (Wfly_mask(i,j+1,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * W3(i,j+1,k) * recdymin2
              if (Wfly_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * W3(i,j-1,k) * recdymin2
              if (Wflz_mask(i,j,k+1)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j,k+1) * W3(i,j,k+1) * recdzmin2
              if (Wflz_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j,k) * W3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 3
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Wnz, tnz2
       do bj = 1, Wny, tny2
        do bi = 1, Wnx, tnx2
         do k = bk, min(bk+tnz2-1,Wnz)
          do j = bj, min(bj+tny2-1,Wny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Wnx), 2
            if (Wtype(i,j,k)<=0) then
              p = Ap * work(i,j,k) + W2(i,j,k) + W(i,j,k)
              p = p * ApW(i,j,k)
              Sw = max(Sw,abs(p-W3(i,j,k)))
              W3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do


      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Unz, tnz2
       do bj = 1, Uny, tny2
        do bi = 1, Unx, tnx2
         do k = bk, min(bk+tnz2-1,Unz)
          do j = bj, min(bj+tny2-1,Uny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Unx), 2
            if (Utype(i,j,k)<=0) then
              work(i,j,k) = 0
              if (Uflx_mask(i+1,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i+1,j,k) * U3(i+1,j,k) * recdxmin2
              if (Uflx_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j,k) * U3(i-1,j,k) * recdxmin2
              if (Ufly_mask(i,j+1,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * U3(i,j+1,k) * recdymin2
              if (Ufly_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * U3(i,j-1,k) * recdymin2
              if (Uflz_mask(i,j,k+1)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * U3(i,j,k+1) * recdzmin2
              if (Uflz_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * U3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 1
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Unz, tnz2
       do bj = 1, Uny, tny2
        do bi = 1, Unx, tnx2
         do k = bk, min(bk+tnz2-1,Unz)
          do j = bj, min(bj+tny2-1,Uny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Unx), 2
            if (Utype(i,j,k)<=0) then
              p = Ap * work(i,j,k) + U2(i,j,k) + U(i,j,k)
              p = p * ApU(i,j,k)
              Su = max(Su,abs(p-U3(i,j,k)))
              U3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
         
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Vnz, tnz2
       do bj = 1, Vny, tny2
        do bi = 1, Vnx, tnx2
         do k = bk, min(bk+tnz2-1,Vnz)
          do j = bj, min(bj+tny2-1,Vny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Vnx), 2
            if (Vtype(i,j,k)<=0) then
              work(i,j,k) = 0
              if (Vflx_mask(i+1,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * V3(i+1,j,k) * recdxmin2
              if (Vflx_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * V3(i-1,j,k) * recdxmin2
              if (Vfly_mask(i,j+1,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j+1,k) * V3(i,j+1,k) * recdymin2
              if (Vfly_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j,k) * V3(i,j-1,k) * recdymin2
              if (Vflz_mask(i,j,k+1)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * V3(i,j,k+1) * recdzmin2
              if (Vflz_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * V3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 2
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Vnz, tnz2
       do bj = 1, Vny, tny2
        do bi = 1, Vnx, tnx2
         do k = bk, min(bk+tnz2-1,Vnz)
          do j = bj, min(bj+tny2-1,Vny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Vnx), 2
            if (Vtype(i,j,k)<=0) then
              p = Ap * work(i,j,k) + V2(i,j,k) + V(i,j,k)
              p = p * ApV(i,j,k)
              Sv = max(Sv,abs(p-V3(i,j,k)))
              V3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
         
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Wnz, tnz2
       do bj = 1, Wny, tny2
        do bi = 1, Wnx, tnx2
         do k = bk, min(bk+tnz2-1,Wnz)
          do j = bj, min(bj+tny2-1,Wny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Wnx), 2
            if (Wtype(i,j,k)<=0) then
              work(i,j,k) = 0
              if (Wflx_mask(i+1,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * W3(i+1,j,k) * recdxmin2
              if (Wflx_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * W3(i-1,j,k) * recdxmin2
              if (Wfly_mask(i,j+1,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * W3(i,j+1,k) * recdymin2
              if (Wfly_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * W3(i,j-1,k) * recdymin2
              if (Wflz_mask(i,j,k+1)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j,k+1) * W3(i,j,k+1) * recdzmin2
              if (Wflz_mask(i,j,k)) &
                work(i,j,k) = work(i,j,k) + &
                  nu(i,j,k) * W3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Wnz, tnz2
       do bj = 1, Wny, tny2
        do bi = 1, Wnx, tnx2
         do k = bk, min(bk+tnz2-1,Wnz)
          do j = bj, min(bj+tny2-1,Wny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Wnx), 2
            if (Wtype(i,j,k)<=0) then
              p = Ap * work(i,j,k) + W2(i,j,k) + W(i,j,k)
              p = p * ApW(i,j,k)
              Sw = max(Sw,abs(p-W3(i,j,k)))
              W3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
      !$omp end parallel

      S = max(Su,Sv,Sw)
#ifdef PAR
      S = par_co_max(S)
#endif
      if (master) write(*,*) "CN ", l, S

      if (S<=epsCN) exit
    end do


  end subroutine ImplicitDiffusion_Iterations





end module MomentumDiffusion
