module Subgrid_MixedTimeScale
  use Parameters
  use Filters
  use Tiling

  implicit none

  private
  public SGS_MixedTimeScale
  
  real(knd), public :: C_MixedTimeScale = 0.05_knd

contains

  subroutine SGS_MixedTimeScale(U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W    
!     real(knd), parameter :: C_mts = 0.05
    real(knd), dimension(:,:,:), allocatable, save :: Uf, Vf, Wf
    integer :: i, j, k, bi, bj, bk
    real(knd) :: rec_width, T_s, k_es
    integer :: tnx, tny, tnz

    integer, parameter :: narr = 7
    real(knd) :: C_mts
    C_mts = C_MixedTimeScale

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

    call FilterTopHat_no_Utype(Uf, U)
    call FilterTopHat_no_Utype(Vf, V)
    call FilterTopHat_no_Utype(Wf, W)

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

  subroutine FilterTopHat_no_Utype(Uf, U)
    use ArrayUtilities, only: set
    real(knd), dimension(-2:,-2:,-2:), intent(out)  :: Uf
    real(knd), dimension(-2:,-2:,-2:), intent(in)  :: U
    integer :: i, j, k, bi, bj, bk
    integer :: mini, minj, mink
    integer :: maxi, maxj, maxk

    real(knd) :: tmp(-2:max(ubound(U,1), ubound(U,2), ubound(U,3)))

    mini = lbound(U,1) + 3
    maxi = ubound(U,1) - 3
    minj = lbound(U,2) + 3
    maxj = ubound(U,2) - 3
    mink = lbound(U,3) + 3
    maxk = ubound(U,3) - 3

    !Shift of the region boundaries at the walls.
    !We do not want to get the outside points into the stencil.
    !It causes too large difference between the original and filtered value.
    if (Btype(We)==BC_NOSLIP.or.Btype(We)==BC_DIRICHLET) then
      mini = mini + 1
    end if
    if (Btype(Ea)==BC_NOSLIP.or.Btype(Ea)==BC_DIRICHLET) then
      maxi = maxi - 1
    end if
    if (Btype(So)==BC_NOSLIP.or.Btype(So)==BC_DIRICHLET) then
      minj = minj + 1
    end if
    if (Btype(No)==BC_NOSLIP.or.Btype(No)==BC_DIRICHLET) then
      maxj = maxj - 1
    end if
    if (Btype(Bo)==BC_NOSLIP.or.Btype(Bo)==BC_DIRICHLET) then
      mink = mink + 1
    end if
    if (Btype(To)==BC_NOSLIP.or.Btype(To)==BC_DIRICHLET) then
      maxk = maxk - 1
    end if
    
    !$omp parallel private(i, j, k, tmp) shared(U, mini, maxi, minj, maxj, mink, maxk)
    
    !filter by the separable kernel in all three directions
    
    !$omp do collapse(2) schedule(dynamic,4)
    do k = mink - 1, maxk + 1
      do j = minj - 1, maxj + 1
        tmp(:ubound(U,1)) = U(:,j,k)
        do i = mini, maxi

            Uf(i,j,k) = 0.25 * (tmp(i+1) + 2 * tmp(i) + tmp(i-1))

        end do
      end do
    end do

    !$omp do schedule(dynamic)
    do k = mink - 1, maxk + 1
      do i = mini, maxi
        tmp(:ubound(U,2)) = Uf(i,:,k)
        do j = minj, maxj
            Uf(i,j,k) = 0.25 * (tmp(j+1) + 2 * tmp(j) + tmp(j-1))
        end do
      end do
    end do

    !$omp do collapse(2) schedule(guided)
    do j = minj, maxj
      do i = mini, maxi
        tmp(:ubound(U,3)) = Uf(i,j,:)
        do k = mink, maxk

            Uf(i,j,k) = 0.25 * (tmp(k+1) + 2 * tmp(k) + tmp(k-1))

        end do
      end do
    end do
    !$omp end parallel
    
    !copy the portions which are not filtered
    
    !$omp do collapse(3)
    do k = -2, mink-1
      do j = -2, ubound(U,2)
        do i = -2, ubound(U,1)
          Uf(i,j,k) = U(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait

    !$omp do collapse(3)
    do k = maxk+1, ubound(U, 3)
      do j = -2, ubound(U,2)
        do i = -2, ubound(U,1)
          Uf(i,j,k) = U(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait
    
    !$omp do collapse(3)
    do k = 1, maxk
      do j = -2, minj-1
        do i = -2, ubound(U,1)
          Uf(i,j,k) = U(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait
    
    !$omp do collapse(3)
    do k = 1, maxk
      do j = maxj+1, ubound(U,2)
        do i = -2, ubound(U,1)
          Uf(i,j,k) = U(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait   

    !$omp do collapse(3)
    do k = 1, maxk
      do j = 1, maxj
        do i = -2, mini-1
          Uf(i,j,k) = U(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait
    
    !$omp do collapse(3)
    do k = 1, maxk
      do j = 1, maxj
        do i = maxi+1, ubound(U,1)
          Uf(i,j,k) = U(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait
        
  end subroutine 

  
end module Subgrid_MixedTimeScale
