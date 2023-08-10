module Subgrid_MixedTimeScale
  use Parameters
  use Filters
  use Tiling

  implicit none

  private
  public SGS_MixedTimeScale

contains

  subroutine SGS_MixedTimeScale(U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W    
    real(knd), parameter :: C_mts = 0.05
    real(knd), dimension(:,:,:), allocatable, save :: Uf, Vf, Wf
    integer :: i, j, k, bi, bj, bk
    real(knd) :: rec_width, T_s, k_es
    integer :: tnx, tny, tnz

    integer, parameter :: narr = 7

    if (gridtype==GRID_VARIABLE_Z) then
      call error_stop("MTS subgrid only implemented for uniform grid.")
    end if
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    rec_width = 1 / (dxmin*dymin*dzmin)**(1._knd/3._knd)

    !allocation on assignment on the first call
    if (.not.allocated(Uf)) then
      Uf = U
      Vf = V
      Wf = W
    end if

    call FilterTopHat(Uf, U, Utype)
    call FilterTopHat(Vf, V, Vtype)
    call FilterTopHat(Wf, W, Wtype)

    !$omp parallel do private(k_es, T_s, i, j, k, bi, bj, bk) schedule(runtime) collapse(3)
    do bk = 1, Prnz, tnz
      do bj = 1, Prny, tny
        do bi = 1, Prnx, tnx
          do k = bk, min(bk+tnz-1,Prnz)
            do j = bj, min(bj+tny-1,Prny)
              do i = bi ,min(bi+tnx-1,Prnx)
                k_es = (U(i,j,k) - Uf(i,j,k))**2 + &
                       (V(i,j,k) - Vf(i,j,k))**2 + &
                       (W(i,j,k) - Wf(i,j,k))**2

                T_s = TimeScale(i, j, k, k_es)

                Viscosity(i,j,k) = C_mts * k_es * T_s
         
                Viscosity(i,j,k) = Viscosity(i,j,k) + molecular_viscosity
              end do
            end do
          end do
        end do
      end do
    end do

  contains

    pure function TimeScale(i,j,k,k_es) result(res)
      real(knd) :: res
      integer,intent(in) :: i,j,k
      real(knd), intent(in) :: k_es
      real(knd), parameter :: Ct = 10
      real(knd) :: D(3,3), magD
  
      call GradientTensorUG(D, i, j, k)
      
      !symmetric part
      D = (D + transpose(D)) / 2
      !deviatoric part
      D = D - (D(1,1) + D(2,2) + D(3,3)) / 3
      
      magD = sqrt(2 * sum(D**2))
                
      res = sqrt(k_es) * rec_width

      res = res + magD / Ct
                  
      res = max(res, epsilon(1._knd))
                  
      res = 1 / res
    end function

    pure subroutine GradientTensorUG(g,i,j,k)
      real(knd),intent(out) :: g(3,3)
      integer,intent(in) :: i,j,k

      g(1,1) = (U(i,j,k)-U(i-1,j,k)) / dxmin
      g(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k)) / (4._knd*dymin)
      g(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) / (4._knd*dzmin)

      g(2,2) = (V(i,j,k)-V(i,j-1,k)) / dymin
      g(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k)) / (4._knd*dxmin)
      g(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) / (4._knd*dzmin)

      g(3,3) = (W(i,j,k)-W(i,j,k-1)) / dzmin
      g(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1)) / (4._knd*dxmin)
      g(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) / (4._knd*dymin)
    end subroutine GradientTensorUG

  end subroutine

  
end module Subgrid_MixedTimeScale
