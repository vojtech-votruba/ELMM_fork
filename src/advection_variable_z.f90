module MomentumAdvection_variable_z
  use Parameters
  use Boundaries
  use ArrayUtilities
  use Tiling

  implicit none

contains





  subroutine CDUdiv_variable_z(U2, U, V, W)
    real(knd), contiguous, intent(out) :: U2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd) :: Ax, Ay
    integer :: i, j, k, bi, bj, bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 4
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    Ax = 0.25_knd / dxmin
    Ay = 0.25_knd / dymin

    !$omp parallel do private(i, j, k, bi, bj, bk) schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
      do bj = 1, Uny, tny
        do bi = 1, Unx, tnx
          do k = bk, min(bk+tnz-1, Unz)
            do j = bj, min(bj+tny-1, Uny)
              do i = bi, min(bi+tnx-1, Unx)
                U2(i,j,k) = - ((Ax*(U(i+1,j,k) + U(i,j,k)) * (U(i+1,j,k) + U(i,j,k)) &
                              - Ax*(U(i,j,k) + U(i-1,j,k)) * (U(i,j,k) + U(i-1,j,k))) &
                             + (Ay*(U(i,j+1,k) + U(i,j,k)) * (V(i+1,j,k) + V(i,j,k)) &
                              - Ay*(U(i,j,k) + U(i,j-1,k)) * (V(i+1,j-1,k) + V(i,j-1,k))) &
                             + ( (U(i,j,k+1) + U(i,j,k)) * (W(i+1,j,k) + W(i,j,k)) &
                              - (U(i,j,k) + U(i,j,k-1)) * (W(i+1,j,k-1) + W(i,j,k-1))) / dzPr(k)/4)
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine






  subroutine CDVdiv_variable_z(V2, U, V, W)
    real(knd), contiguous, intent(out) :: V2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd) :: Ax, Ay
    integer :: i, j, k, bi, bj, bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 4
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    Ax = 0.25_knd / dxmin
    Ay = 0.25_knd / dymin

    !$omp parallel do private(i, j, k, bi, bj, bk) schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
      do bj = 1, Vny, tny
        do bi = 1, Vnx, tnx
          do k = bk, min(bk+tnz-1, Vnz)
            do j = bj, min(bj+tny-1, Vny)
              do i = bi, min(bi+tnx-1, Vnx)
                V2(i,j,k) = - ((Ay*(V(i,j+1,k) + V(i,j,k)) * (V(i,j+1,k) + V(i,j,k)) &
                               -Ay*(V(i,j,k) + V(i,j-1,k)) * (V(i,j,k) + V(i,j-1,k))) &
                              +(Ax*(V(i+1,j,k) + V(i,j,k)) * (U(i,j+1,k) + U(i,j,k)) &
                               -Ax*(V(i,j,k) + V(i-1,j,k)) * (U(i-1,j+1,k) + U(i-1,j,k))) &
                              +( (V(i,j,k+1) + V(i,j,k)) * (W(i,j+1,k) + W(i,j,k)) &
                               -(V(i,j,k) + V(i,j,k-1)) * (W(i,j+1,k-1) + W(i,j,k-1))) / dzPr(k)/4)
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine





  subroutine CDWdiv_variable_z(W2, U, V, W)
    real(knd), contiguous, intent(out) :: W2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd) :: Ax, Ay
    integer :: i, j, k, bi, bj, bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 4
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    Ax = 0.25_knd / dxmin
    Ay = 0.25_knd / dymin

    !$omp parallel do private(i, j, k, bi, bj, bk) schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
      do bj = 1, Wny, tny
        do bi = 1, Wnx, tnx
          do k = bk, min(bk+tnz-1, Wnz)
            do j = bj, min(bj+tny-1, Wny)
              do i = bi, min(bi+tnx-1, Wnx)
                W2(i,j,k) = - (((W(i,j,k+1) + W(i,j,k)) * (W(i,j,k+1) + W(i,j,k)) &
                              - (W(i,j,k) + W(i,j,k-1)) * (W(i,j,k) + W(i,j,k-1))) / dzW(k)/4 &
                             + (Ay*(W(i,j+1,k) + W(i,j,k)) * (V(i,j,k+1) + V(i,j,k)) &
                              - Ay*(W(i,j,k) + W(i,j-1,k)) * (V(i,j-1,k+1) + V(i,j-1,k))) &
                             + (Ax*(W(i+1,j,k) + W(i,j,k)) * (U(i,j,k+1) + U(i,j,k)) &
                              - Ax*(W(i,j,k) + W(i-1,j,k)) * (U(i-1,j,k+1) + U(i-1,j,k))))
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine











  subroutine CDUadv_variable_z(U2, U, V, W)
    real(knd), contiguous, intent(out) :: U2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd) :: Ax, Ay, Vadv, Wadv
    integer :: i, j, k, bi, bj, bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 4
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    Ax = 0.5_knd / dxmin
    Ay = 0.125_knd / dymin

    !$omp parallel do private(i, j, k, bi, bj, bk, Vadv, Wadv) schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
      do bj = 1, Uny, tny
        do bi = 1, Unx, tnx
          do k = bk, min(bk+tnz-1, Unz)
            do j = bj, min(bj+tny-1, Uny)
              do i = bi, min(bi+tnx-1, Unx)
                Vadv = ( V(i,j,k) + V(i+1,j,k) + V(i,j-1,k) + V(i+1,j-1,k) )
                Wadv = ( W(i,j,k) + W(i+1,j,k) + W(i,j,k-1) + W(i+1,j,k-1) ) / 4
                U2(i,j,k) = U2(i,j,k) &
                           - (Ax*(U(i+1,j,k)-U(i-1,j,k)) * U(i,j,k) &
                           +  Ay*(U(i,j+1,k)-U(i,j-1,k)) * Vadv&
                           +  (U(i,j,k+1)-U(i,j,k-1)) / (zPr(k+1)-zPr(k-1)) * Wadv )
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine






  subroutine CDVadv_variable_z(V2, U, V, W)
    real(knd), contiguous, intent(out) :: V2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd) :: Ax, Ay, Uadv, Wadv
    integer :: i, j, k, bi, bj, bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 4
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    Ax = 0.125_knd / dxmin
    Ay = 0.5_knd / dymin

    !$omp parallel do private(i, j, k, bi, bj, bk, Uadv, Wadv) schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
      do bj = 1, Vny, tny
        do bi = 1, Vnx, tnx
          do k = bk, min(bk+tnz-1, Vnz)
            do j = bj, min(bj+tny-1, Vny)
              do i = bi, min(bi+tnx-1, Vnx)
                Uadv = ( U(i,j,k) + U(i,j+1,k) + U(i-1,j,k) + U(i-1,j+1,k) )
                Wadv = ( W(i,j,k) + W(i,j+1,k) + W(i,j,k-1) + W(i,j+1,k-1) ) / 4
                V2(i,j,k) = V2(i,j,k) &
                           - (Ax*(V(i+1,j,k)-V(i-1,j,k)) * Uadv&
                           +  Ay*(V(i,j+1,k)-V(i,j-1,k)) * V(i,j,k) &
                           +  (V(i,j,k+1)-V(i,j,k-1)) / (zPr(k+1)-zPr(k-1)) * Wadv )
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine





  subroutine CDWadv_variable_z(W2, U, V, W)
    real(knd), contiguous, intent(out) :: W2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd) :: Ax, Ay, Uadv, Vadv
    integer :: i, j, k, bi, bj, bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 4
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    Ax = 0.125_knd / dxmin
    Ay = 0.125_knd / dymin

    !$omp parallel do private(i, j, k, bi, bj, bk, Uadv, Vadv) schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
      do bj = 1, Wny, tny
        do bi = 1, Wnx, tnx
          do k = bk, min(bk+tnz-1, Wnz)
            do j = bj, min(bj+tny-1, Wny)
              do i = bi, min(bi+tnx-1, Wnx)
                Uadv = ( U(i,j,k) + U(i,j,k+1) + U(i-1,j,k) + U(i-1,j,k+1) )
                Vadv = ( V(i,j,k) + V(i,j,k+1) + V(i,j-1,k) + V(i,j-1,k+1) )
                W2(i,j,k) = W2(i,j,k) &
                           - (Ax*(W(i+1,j,k)-W(i-1,j,k)) * Uadv&
                           +  Ay*(W(i,j+1,k)-W(i,j-1,k)) * Vadv&
                           +  (W(i,j,k+1)-W(i,j,k-1)) / (zW(k+1)-zW(k-1)) * W(i,j,k) )
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine















  subroutine CDU_variable_z(U2, U, V, W)
    real(knd), contiguous, intent(out) :: U2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)

    call set(U2, 0)
    call CDUdiv_variable_z(U2, U, V, W)
    call CDUadv_variable_z(U2, U, V, W)
    call multiply(U2, 0.5_knd)
  end subroutine






  subroutine CDV_variable_z(V2, U, V, W)
    real(knd), contiguous, intent(out) :: V2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)

    call set(V2, 0)
    call CDVdiv_variable_z(V2, U, V, W)
    call CDVadv_variable_z(V2, U, V, W)
    call multiply(V2, 0.5_knd)
  end subroutine





  subroutine CDW_variable_z(W2, U, V, W)
    real(knd), contiguous, intent(out) :: W2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)

    call set(W2, 0)
    call CDWdiv_variable_z(W2, U, V, W)
    call CDWadv_variable_z(W2, U, V, W)
    call multiply(W2, 0.5_knd)
  end subroutine



end module MomentumAdvection_variable_z
