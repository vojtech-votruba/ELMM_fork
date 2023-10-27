module Filters

  use Parameters
  use Tiling

  implicit none

  private

  public filtertype, filter_ratios, Filter, FilterTopHat, FilterTopHatSimple, Filter3ord

  integer :: filtertype = 1

  real(knd) :: filter_ratios(0:1) = [ 1.0_knd, 2._knd]

 
contains

  subroutine FilterTopHat(Uf, U, Utype)
    use ArrayUtilities, only: set
    real(knd), dimension(-2:,-2:,-2:), intent(out)  :: Uf
    real(knd), dimension(-2:,-2:,-2:), intent(in)  :: U
    integer, dimension(-2:,-2:,-2:), intent(in)  :: Utype
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
          if (all(Utype(i-1:i+1,j,k)<=0)) then
            Uf(i,j,k) = 0.25 * (tmp(i+1) + 2 * tmp(i) + tmp(i-1))
          else
            Uf(i,j,k) = tmp(i)
          end if
        end do
      end do
    end do

    !$omp do schedule(dynamic)
    do k = mink - 1, maxk + 1
      do i = mini, maxi
        tmp(:ubound(U,2)) = Uf(i,:,k)
        do j = minj, maxj
          if (all(Utype(i,j-1:j+1,k)<=0)) &
            Uf(i,j,k) = 0.25 * (tmp(j+1) + 2 * tmp(j) + tmp(j-1))
        end do
      end do
    end do

    !$omp do collapse(2) schedule(guided)
    do j = minj, maxj
      do i = mini, maxi
        tmp(:ubound(U,3)) = Uf(i,j,:)
        do k = mink, maxk
          if (all(Utype(i,j,k-1:k+1)<=0)) then
            Uf(i,j,k) = 0.25 * (tmp(k+1) + 2 * tmp(k) + tmp(k-1))
          end if
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


  subroutine FilterTopHatSimple(Uf, U)
    use ArrayUtilities, only: set
    real(knd), dimension(0:,0:,0:), intent(out)  :: Uf
    real(knd), dimension(0:,0:,0:), intent(in)  :: U
    integer :: i, j, k, bi, bj, bk
    integer :: mini, minj, mink
    integer :: maxi, maxj, maxk

    real(knd) :: tmp(0:max(ubound(U,1), ubound(U,2), ubound(U,3)))

    mini = 1
    maxi = ubound(U,1) - 1
    minj = 1
    maxj = ubound(U,2) - 1
    mink = 1
    maxk = ubound(U,3) - 1

    Uf = U

    !$omp parallel private(i, j, k, tmp) shared(U, mini, maxi, minj, maxj, mink, maxk)
       
    !filter by the separable kernel in all three directions
    
    !$omp do collapse(2) schedule(dynamic,4)
    do k = mink, maxk
      do j = minj, maxj
        tmp(:ubound(U,1)) = U(:,j,k)
        do i = mini, maxi
            Uf(i,j,k) = 0.25 * (tmp(i+1) + 2 * tmp(i) + tmp(i-1))
        end do
      end do
    end do

    !$omp do schedule(dynamic)
    do k = mink, maxk
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
  end subroutine FilterTopHatSimple


  
  subroutine Filter3ord(U, Utype, dir)
    !3D filter with 3 vanishing moments width 2 delta x
    !Vasilyev, Lund, Moin, 1998, http://dx.doi.org/10.1006/jcph.1998.6060
    real(knd), dimension(-2:,-2:,-2:), intent(inout)  :: U
    integer, dimension(-2:,-2:,-2:), intent(in)  :: Utype
    integer, intent(in) :: dir
    real(knd), parameter :: w(-3:3) = [ -1._knd/32, 0._knd, 9._knd/32, 1._knd/2, 9._knd/32, 0._knd, -1._knd/32 ]
    real(knd) :: S, q
    integer :: i, j, k, ii, jj, kk
    integer :: mini, minj, mink
    integer :: maxi, maxj, maxk

    real(knd) :: tmp(-2:ubound(U,dir))
    
    mini = lbound(U,1) + 3
    maxi = ubound(U,1) - 3
    minj = lbound(U,2) + 3
    maxj = ubound(U,2) - 3
    mink = lbound(U,3) + 3
    maxk = ubound(U,3) - 3

    !The condition is intentionally ==0( not <=0)
    !  to avoid points nearest to the boundary.
    
    if (dir==1) then
      !$omp parallel do private(i, j, k, ii, S, q, tmp) shared(U, Utype, mini, maxi, minj, maxj, mink, maxk) schedule(runtime)
      do k = mink, maxk
       do j = minj, maxj
         tmp = U(:,j,k)
         do i = 1, maxi
           if (Utype(i,j,k) == 0) then
             S = 0
             q = 0
             do ii = i-3, i+3
               if (Utype(ii,j,k)<=0) then
                 S = S + w(ii-i) * tmp(ii)
                 q = q + w(ii-i)
               end if
             end do
             if (abs(q)>0.74) U(i,j,k) = S / q
           end if
         end do
       end do
      end do
      !$omp end parallel do
    else if (dir==2) then
      !$omp parallel do private(i, j, k, jj, S, q, tmp) shared(U, Utype, mini, maxi, minj, maxj, mink, maxk) schedule(runtime)
      do k = mink, maxk
       do i = mini, maxi
         tmp = U(i,:,k)
         do j = 1, maxj
           if (Utype(i,j,k) == 0) then
             S = 0
             q = 0
             do jj = j-3, j+3
               if (Utype(i,jj,k)<=0) then
                 S = S + w(jj-j) * tmp(jj)
                 q = q + w(jj-j)
               end if
             end do
             if (abs(q)>.74) U(i,j,k) = S / q
           end if
         end do
       end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(i, j, k, kk, S, q, tmp) shared(U, Utype, mini, maxi, minj, maxj, mink, maxk) schedule(runtime)
      do j = minj, maxj
       do i = mini, maxi
         tmp = U(i,j,:)
         do k = 1, maxk
           if (Utype(i,j,k) == 0) then               
             S = 0
             q = 0
             do kk = k-3, k+3
               if (Utype(i,j,kk)<=0) then
                 S = S + w(kk-k) * tmp(kk)
                 q = q + w(kk-k)
               end if
             end do
             if (abs(q)>.74) U(i,j,k) = S / q
           end if
         end do
       end do
      end do
      !$omp end parallel do
    end if
  end subroutine Filter3ord



  subroutine Filter(U, Utype)    !Calls a selected filter  U <- filt(U)
    real(knd), dimension(-2:,-2:,-2:), intent(inout) :: U
    integer, dimension(-2:,-2:,-2:), intent(in)  :: Utype

    if (filtertype==1) then
      if (Prnx==1) then
        call Filter3ord(U, Utype, 2)
        call Filter3ord(U, Utype, 3)
      else if (Prny==1) then
        call Filter3ord(U, Utype, 1)
        call Filter3ord(U, Utype, 3)
      else if (Prnz==1) then
        call Filter3ord(U, Utype, 3)
        call Filter3ord(U, Utype, 3)
      else
        call Filter3ord(U, Utype, 1)
        call Filter3ord(U, Utype, 2)
        call Filter3ord(U, Utype, 3)
      end if
    end if

  end subroutine Filter


end module Filters
