module ScalarDiffusion

  use Parameters
  use ArrayUtilities
  use ScalarBoundaries
  use Wallmodels
  use ImmersedBoundary, only: Scalar_ImmersedBoundaries
  use Tiling, only: tilenx, tileny, tilenz
  
  
  implicit none
  
  private
  
  public ScalarDiffusion_explicit, ScalarDiffusion_nobranch, ScalarDiffusion_implicit, &
         AddScalarDiffVector, ComputeTDiff, &
         ScalarDiffusion_Deallocate, &
         boundary_interface, &
         constPr_sgs

  real(knd), parameter :: constPr_sgs = 0.6 !constant (default) value of Pr_sgs, which may be further refined elsewhere

  real(knd), allocatable :: Scal3(:,:,:)
  real(knd), allocatable :: Ap(:,:,:)

  abstract interface
    subroutine boundary_interface(array)
      use Parameters
      real(knd), intent(inout) :: array(-1:,-1:,-1:)
    end subroutine
  end interface


contains



  subroutine ScalarDiffusion_implicit(Scal2, Scal, sctype, boundary_procedure, coef)
#ifdef PAR
    use custom_par, only: par_co_max
#endif
    real(knd), contiguous, intent(in)    :: Scal(-1:,-1:,-1:)
    real(knd), contiguous, intent(inout) :: Scal2(-1:,-1:,-1:)
    real(knd), intent(in) :: coef
    integer, intent(in) :: sctype
    procedure(boundary_interface) :: boundary_procedure
    integer   :: nx, ny, nz, i, j, k, bi, bj, bk, l
    real(knd) :: p, S
    real(knd) :: A, Ax, Ay, Az
    integer :: tnx, tny, tnz, tnx2, tny2, tnz2
    
    integer, parameter :: narr = 3, narr2 = 5
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    tnx2 = tilenx(narr2)
    tny2 = tileny(narr2)
    tnz2 = tilenz(narr2)

    
    

    if (.not.allocated(Scal3)) then
      allocate(Scal3(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
      allocate(Ap(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    end if

    nx = Prnx
    ny = Prny
    nz = Prnz


    if (molecular_viscosity > 0) then


      A = coef
      Ax = 1 / (dxmin**2)
      Ay = 1 / (dymin**2)
      Az = 1 / (dzmin**2)

      !$omp parallel private (i,j,k,bi,bj,bk)

      !initital value using forward Euler
      if (gridtype==GRID_UNIFORM) then
        !$omp do schedule(runtime) collapse(3)
        do bk = 1, Prnz, tnz
         do bj = 1, Prny, tny
          do bi = 1, Prnx, tnx
           do k = bk, min(bk+tnz-1,Prnz)
            do j = bj, min(bj+tny-1,Prny)
             do i = bi, min(bi+tnx-1,Prnx)
               if (Prtype(i,j,k)<=0) then
                 Scal3(i,j,k) = ((TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal(i+1,j,k)-Scal(i,j,k))-&
                   (TDiff(i,j,k)+TDiff(i-1,j,k)) * (Scal(i,j,k)-Scal(i-1,j,k)))*Ax

                 Scal3(i,j,k) = Scal3(i,j,k)  + &
                   ((TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal(i,j+1,k)-Scal(i,j,k))-&
                   (TDiff(i,j,k)+TDiff(i,j-1,k)) * (Scal(i,j,k)-Scal(i,j-1,k)))*Ay

                 Scal3(i,j,k) = Scal3(i,j,k)  + &
                   ((TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal(i,j,k+1)-Scal(i,j,k))-&
                   (TDiff(i,j,k)+TDiff(i,j,k-1)) * (Scal(i,j,k)-Scal(i,j,k-1)))*Az
               else
                 Scal3(i,j,k) = 0
               end if
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp end do
      else
        !$omp do schedule(runtime) collapse(3)
        do bk = 1, Prnz, tnz
         do bj = 1, Prny, tny
          do bi = 1, Prnx, tnx
           do k = bk, min(bk+tnz-1,Prnz)
            do j = bj, min(bj+tny-1,Prny)
             do i = bi, min(bi+tnx-1,Prnx)
               Scal3(i,j,k) = (((TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal(i+1,j,k)-Scal(i,j,k)) / dxU(i)-&
                 (TDiff(i,j,k)+TDiff(i-1,j,k)) * (Scal(i,j,k)-Scal(i-1,j,k)) / dxU(i-1)) / (dxPr(i)) + &
                ((TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal(i,j+1,k)-Scal(i,j,k)) / dyV(j)-&
                 (TDiff(i,j,k)+TDiff(i,j-1,k)) * (Scal(i,j,k)-Scal(i,j-1,k)) / dyV(j-1)) / (dyPr(j)) + &
                ((TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal(i,j,k+1)-Scal(i,j,k)) / dzW(k)-&
                 (TDiff(i,j,k)+TDiff(i,j,k-1)) * (Scal(i,j,k)-Scal(i,j,k-1)) / dzW(k-1)) / (dzPr(k)))
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp end do
      end if
      !$omp end parallel

      Ax = 1 / (4 * dxmin**2)
      Ay = 1 / (4 * dymin**2)
      Az = 1 / (4 * dzmin**2)

      !Scal2 = Scal + Scal3 * A
      call assign(Scal2, Scal)
      call add_multiplied(Scal2, Scal3, A)

      call boundary_procedure(Scal2)

      call Scalar_ImmersedBoundaries(Scal2)
      call Scalar_ImmersedBoundaries(Scal3)

      !$omp parallel private(i,j,k,bi,bj,bk)
      if (gridtype==GRID_UNIFORM) then
       !$omp do schedule(runtime) collapse(3)
       do bk = 1, Prnz, tnz
        do bj = 1, Prny, tny
         do bi = 1, Prnx, tnx
          do k = bk, min(bk+tnz-1,Prnz)
           do j = bj, min(bj+tny-1,Prny)
            do i = bi, min(bi+tnx-1,Prnx)
              Ap(i,j,k) = 1 / (1._knd/A+(((TDiff(i+1,j,k)+TDiff(i,j,k)) + &
                                   (TDiff(i,j,k)+TDiff(i-1,j,k)))*Ax + &
                                   ((TDiff(i,j+1,k)+TDiff(i,j,k)) + &
                                   (TDiff(i,j,k)+TDiff(i,j-1,k)))*Ay + &
                                   ((TDiff(i,j,k+1)+TDiff(i,j,k)) + &
                                   (TDiff(i,j,k)+TDiff(i,j,k-1)))*Az))
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do

     else

       !$omp do schedule(runtime) collapse(3)
       do bk = 1, Prnz, tnz
        do bj = 1, Prny, tny
         do bi = 1, Prnx, tnx
          do k = bk, min(bk+tnz-1,Prnz)
           do j = bj, min(bj+tny-1,Prny)
            do i = bi, min(bi+tnx-1,Prnx)
              Ap(i,j,k) = 1 / (1._knd/A+(((TDiff(i+1,j,k)+TDiff(i,j,k)) / dxU(i) + &
                                   (TDiff(i,j,k)+TDiff(i-1,j,k)) / dxU(i-1)) / (4._knd*dxPr(i)) + &
                                   ((TDiff(i,j+1,k)+TDiff(i,j,k)) / dyV(j) + &
                                   (TDiff(i,j,k)+TDiff(i,j-1,k)) / dyV(j-1)) / (4._knd*dyPr(j)) + &
                                   ((TDiff(i,j,k+1)+TDiff(i,j,k)) / dzW(k) + &
                                   (TDiff(i,j,k)+TDiff(i,j,k-1)) / dzW(k-1)) / (4._knd*dzPr(k))))
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do

     end if
     !$omp end parallel

     do l = 1, maxCNiter
       S = 0
       call boundary_procedure(Scal2)

       if (gridtype==GRID_UNIFORM) then
        !$omp parallel private(i,j,k,p) reduction(max:S)
        !$omp do schedule(runtime) collapse(3)
        do bk = 1, Prnz, tnz2
         do bj = 1, Prny, tny2
          do bi = 1, Prnx, tnx2
           do k = bk, min(bk+tnz2-1,Prnz)
            do j = bj, min(bj+tny2-1,Prny)
             do i = bi+mod(bi+j+k-1, 2), min(bi+tnx2-1,Prnx), 2
               if (Prtype(i,j,k)<=0) then
                 p = (Scal(i,j,k)/A) + (Scal3(i,j,k)/4._knd + &
                  ((TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal2(i+1,j,k))-&
                   (TDiff(i,j,k)+TDiff(i-1,j,k)) * (-Scal2(i-1,j,k)))*Ax + &
                  ((TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal2(i,j+1,k))-&
                   (TDiff(i,j,k)+TDiff(i,j-1,k)) * (-Scal2(i,j-1,k)))*Ay + &
                  ((TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal2(i,j,k+1))-&
                   (TDiff(i,j,k)+TDiff(i,j,k-1)) * (-Scal2(i,j,k-1)))*Az&
                  )
                  p = p * Ap(i,j,k)
                  S = max(S, abs(p-Scal2(i,j,k)))
                  Scal2(i,j,k) = p
               end if
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp end do
        !$omp do schedule(runtime) collapse(3)
        do bk = 1, Prnz, tnz2
         do bj = 1, Prny, tny2
          do bi = 1, Prnx, tnx2
           do k = bk, min(bk+tnz2-1,Prnz)
            do j = bj, min(bj+tny2-1,Prny)
             do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Prnx), 2
               if (Prtype(i,j,k)<=0) then
                 p = (Scal(i,j,k)/A) + (Scal3(i,j,k)/4._knd + &
                  ((TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal2(i+1,j,k))-&
                   (TDiff(i,j,k)+TDiff(i-1,j,k)) * (-Scal2(i-1,j,k)))*Ax + &
                  ((TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal2(i,j+1,k))-&
                   (TDiff(i,j,k)+TDiff(i,j-1,k)) * (-Scal2(i,j-1,k)))*Ay + &
                  ((TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal2(i,j,k+1))-&
                   (TDiff(i,j,k)+TDiff(i,j,k-1)) * (-Scal2(i,j,k-1)))*Az&
                  )
                  p = p * Ap(i,j,k)
                  S = max(S, abs(p-Scal2(i,j,k)))
                  Scal2(i,j,k) = p
               end if
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp end do
        !$omp endparallel
       else
        !$omp parallel private(i,j,k,p) reduction(max:S)
        !$omp do schedule(runtime) collapse(3)
        do bk = 1, Prnz, tnz2
         do bj = 1, Prny, tny2
          do bi = 1, Prnx, tnx2
           do k = bk, min(bk+tnz2-1,Prnz)
            do j = bj, min(bj+tny2-1,Prny)
             do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Prnx), 2
               p = (Scal(i,j,k)/A) + (Scal3(i,j,k) + &
                ((TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal2(i+1,j,k)) / dxU(i)-&
                 (TDiff(i,j,k)+TDiff(i-1,j,k)) * (-Scal2(i-1,j,k)) / dxU(i-1)) / (dxPr(i)) + &
                ((TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal2(i,j+1,k)) / dyV(j)-&
                 (TDiff(i,j,k)+TDiff(i,j-1,k)) * (-Scal2(i,j-1,k)) / dyV(j-1)) / (dyPr(j)) + &
                ((TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal2(i,j,k+1)) / dzW(k)-&
                 (TDiff(i,j,k)+TDiff(i,j,k-1)) * (-Scal2(i,j,k-1)) / dzW(k-1)) / (dzPr(k)) &
                )/4._knd
                p = p * Ap(i,j,k)
                S = max(S, abs(p-Scal2(i,j,k)))
                Scal2(i,j,k) = p
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp end do
        !$omp do schedule(runtime) collapse(3)
        do bk = 1, Prnz, tnz2
         do bj = 1, Prny, tny2
          do bi = 1, Prnx, tnx2
           do k = bk, min(bk+tnz2-1, Prnz)
            do j = bj, min(bj+tny2-1, Prny)
             do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Prnx), 2
               p = (Scal(i,j,k)/A) + (Scal3(i,j,k) + &
                ((TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal2(i+1,j,k)) / dxU(i)-&
                 (TDiff(i,j,k)+TDiff(i-1,j,k)) * (-Scal2(i-1,j,k)) / dxU(i-1)) / (dxPr(i)) + &
                ((TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal2(i,j+1,k)) / dyV(j)-&
                 (TDiff(i,j,k)+TDiff(i,j-1,k)) * (-Scal2(i,j-1,k)) / dyV(j-1)) / (dyPr(j)) + &
                ((TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal2(i,j,k+1)) / dzW(k)-&
                 (TDiff(i,j,k)+TDiff(i,j,k-1)) * (-Scal2(i,j,k-1)) / dzW(k-1)) / (dzPr(k)) &
                )/4._knd
                p = p * Ap(i,j,k)
                S = max(S, abs(p-Scal2(i,j,k)))
                Scal2(i,j,k) = p
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp end do
        !$omp endparallel
       end if

#ifdef PAR
       S = par_co_max(S)
#endif
       if (master) write(*,*) "CN scalar", l, S

       if (S<=epsCN) exit
     end do


    else


     call assign(Scal2, Scal)

     call boundary_procedure(Scal2)

     call Scalar_ImmersedBoundaries(Scal2)


    end if
  endsubroutine ScalarDiffusion_implicit








  subroutine ScalarDiffusion_explicit(Scal2, Scal)
    real(knd), contiguous, intent(inout) :: Scal2(-1:,-1:,-1:)
    real(knd), contiguous, intent(in)    :: Scal(-1:,-1:,-1:)
    integer :: i, j, k, bi, bj, bk
    real(knd) :: Ax, Ay, Az
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 6
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    Ax = 1 / (2 * dxmin**2)
    Ay = 1 / (2 * dymin**2)
    Az = 1 / (2 * dzmin**2)

   !$omp parallel do private(i, j, k, bi, bj, bk) schedule(runtime) collapse(3)
    do bk = 1, Prnz, tnz
     do bj = 1, Prny, tny
      do bi = 1, Prnx, tnx
       do k = bk, min(bk+tnz-1, Prnz)
        do j = bj, min(bj+tny-1, Prny)
         do i = bi, min(bi+tnx-1, Prnx)
           if (Scflx_mask(i,j,k)) &
             Scal2(i,j,k) = Scal2(i,j,k) + &
                     (TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal(i+1,j,k)-Scal(i,j,k)) * Ax
           if (Scflx_mask(i-1,j,k)) &
             Scal2(i,j,k) = Scal2(i,j,k) - &
                     (TDiff(i,j,k)+TDiff(i-1,j,k)) * (Scal(i,j,k)-Scal(i-1,j,k)) * Ax

           if (Scfly_mask(i,j,k)) &
             Scal2(i,j,k) = Scal2(i,j,k) + &
                     (TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal(i,j+1,k)-Scal(i,j,k)) * Ay
           if (Scfly_mask(i,j-1,k)) &
             Scal2(i,j,k) = Scal2(i,j,k) - &
                     (TDiff(i,j,k)+TDiff(i,j-1,k)) * (Scal(i,j,k)-Scal(i,j-1,k)) * Ay

           if (Scflz_mask(i,j,k)) &
             Scal2(i,j,k) = Scal2(i,j,k) + &
                     (TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal(i,j,k+1)-Scal(i,j,k)) * Az
           if (Scflz_mask(i,j,k-1)) &
             Scal2(i,j,k) = Scal2(i,j,k) - &
                     (TDiff(i,j,k)+TDiff(i,j,k-1)) * (Scal(i,j,k)-Scal(i,j,k-1)) * Az
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end parallel do
  end subroutine ScalarDiffusion_explicit




  subroutine ScalarDiffusion_nobranch(Scal2, Scal)
    real(knd), contiguous, intent(inout) :: Scal2(-1:,-1:,-1:)
    real(knd), contiguous, intent(in)    :: Scal(-1:,-1:,-1:)
    integer :: i, j, k, bi, bj, bk
    integer :: xi, yj, zk
    real(knd) :: Ax, Ay, Az
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 3
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    Ax = 1 / (2 * dxmin**2)
    Ay = 1 / (2 * dymin**2)
    Az = 1 / (2 * dzmin**2)

    !$omp parallel private(i, j, k, bi, bj, bk, xi, yj, zk)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Prnz, tnz
     do bj = 1, Prny, tny
      do bi = 1, Prnx, tnx
       do k = bk, min(bk+tnz-1, Prnz)
        do j = bj, min(bj+tny-1, Prny)
         do i = bi, min(bi+tnx-1, Prnx)
             Scal2(i,j,k) = Scal2(i,j,k) + &
                     (TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal(i+1,j,k)-Scal(i,j,k)) * Ax
             Scal2(i,j,k) = Scal2(i,j,k) - &
                     (TDiff(i,j,k)+TDiff(i-1,j,k)) * (Scal(i,j,k)-Scal(i-1,j,k)) * Ax

             Scal2(i,j,k) = Scal2(i,j,k) + &
                     (TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal(i,j+1,k)-Scal(i,j,k)) * Ay
             Scal2(i,j,k) = Scal2(i,j,k) - &
                     (TDiff(i,j,k)+TDiff(i,j-1,k)) * (Scal(i,j,k)-Scal(i,j-1,k)) * Ay

             Scal2(i,j,k) = Scal2(i,j,k) + &
                     (TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal(i,j,k+1)-Scal(i,j,k)) * Az
             Scal2(i,j,k) = Scal2(i,j,k) - &
                     (TDiff(i,j,k)+TDiff(i,j,k-1)) * (Scal(i,j,k)-Scal(i,j,k-1)) * Az
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do

    !$omp do
    do i = 1, size(Scflx_points)
      xi = Scflx_points(i)%xi
      yj = Scflx_points(i)%yj
      zk = Scflx_points(i)%zk
      Scal2(xi,yj,zk) = Scal2(xi,yj,zk) - &
                        (TDiff(xi+1,yj,zk)+TDiff(xi,yj,zk)) * (Scal(xi+1,yj,zk)-Scal(xi,yj,zk)) * Ax
    end do
    !$omp end do
    !$omp do
    do i = 1, size(Scflx_points)
      xi = Scflx_points(i)%xi
      yj = Scflx_points(i)%yj
      zk = Scflx_points(i)%zk
      Scal2(xi+1,yj,zk) = Scal2(xi+1,yj,zk) + &
                          (TDiff(xi+1,yj,zk)+TDiff(xi,yj,zk)) * (Scal(xi+1,yj,zk)-Scal(xi,yj,zk)) * Ax
    end do
    !$omp end do nowait
    !$omp do
    do i = 1, size(Scfly_points)
      xi = Scfly_points(i)%xi
      yj = Scfly_points(i)%yj
      zk = Scfly_points(i)%zk
      Scal2(xi,yj,zk) = Scal2(xi,yj,zk) - &
                        (TDiff(xi,yj+1,zk)+TDiff(xi,yj,zk)) * (Scal(xi,yj+1,zk)-Scal(xi,yj,zk)) * Ay
    end do
    !$omp end do
    !$omp do
    do i = 1, size(Scfly_points)
      xi = Scfly_points(i)%xi
      yj = Scfly_points(i)%yj
      zk = Scfly_points(i)%zk
      Scal2(xi,yj+1,zk) = Scal2(xi,yj+1,zk) + &
                          (TDiff(xi,yj+1,zk)+TDiff(xi,yj,zk)) * (Scal(xi,yj+1,zk)-Scal(xi,yj,zk)) * Ay
    end do
    !$omp end do nowait
    !$omp do
    do i = 1, size(Scflz_points)
      xi = Scflz_points(i)%xi
      yj = Scflz_points(i)%yj
      zk = Scflz_points(i)%zk
      Scal2(xi,yj,zk) = Scal2(xi,yj,zk) - &
                        (TDiff(xi,yj,zk+1)+TDiff(xi,yj,zk)) * (Scal(xi,yj,zk+1)-Scal(xi,yj,zk)) * Az
    end do
    !$omp end do
    !$omp do
    do i = 1, size(Scflz_points)
      xi = Scflz_points(i)%xi
      yj = Scflz_points(i)%yj
      zk = Scflz_points(i)%zk
      Scal2(xi,yj,zk+1) = Scal2(xi,yj,zk+1) + &
                          (TDiff(xi,yj,zk+1)+TDiff(xi,yj,zk)) * (Scal(xi,yj,zk+1)-Scal(xi,yj,zk)) * Az
    end do
    !$omp end do
    !$omp end parallel
  end subroutine ScalarDiffusion_nobranch








  subroutine AddScalarDiffVector(ScU, ScV, ScW, Scal, weight, probes_flux, px, py, pz)
    real(knd), contiguous, intent(inout) :: ScU(:,:,:)
    real(knd), contiguous, intent(inout) :: ScV(:,:,:)
    real(knd), contiguous, intent(inout) :: ScW(:,:,:)
    real(knd), contiguous, intent(in)    :: Scal(-1:,-1:,-1:)
    real(knd), intent(in) :: weight
    real(knd), intent(out) :: probes_flux(:,:) !component, position
    integer, intent(in) :: px(:), py(:), pz(:)
    integer :: i, j, k, probe

    !$omp parallel private (i,j,k)
    !$omp do
    do k = 1, Unz
     do j = 1, Uny
      do i = 1, Unx
        ScU(i,j,k) = ScU(i,j,k) + &
                   weight * (TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal(i+1,j,k)-Scal(i,j,k)) / dxmin
      end do
     end do
    end do
    !$omp end do nowait
    !$omp do
    do probe = 1, size(px)
        i = px(probe)
        j = py(probe)
        k = pz(probe)
        probes_flux(1, probe) = probes_flux(1, probe) + (TDiff(i+1,j,k)+TDiff(i,j,k)) * (Scal(i+1,j,k)-Scal(i,j,k)) / dxmin
    end do
    !$omp end do nowait
    !$omp do
    do k = 1, Vnz
     do j = 1, Vny
      do i = 1, Vnx
        ScV(i,j,k) = ScV(i,j,k) + &
                   weight * (TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal(i,j+1,k)-Scal(i,j,k)) / dymin
      end do
     end do
    end do
    !$omp end do nowait
    !$omp do
    do probe = 1, size(px)
        i = px(probe)
        j = py(probe)
        k = pz(probe)
        probes_flux(2, probe) = probes_flux(2, probe) + (TDiff(i,j+1,k)+TDiff(i,j,k)) * (Scal(i,j+1,k)-Scal(i,j,k)) / dymin
    end do
    !$omp end do nowait
    !$omp do
    do k = 1, Wnz
     do j = 1, Wny
      do i = 1, Wnx
        ScW(i,j,k) = ScW(i,j,k) + &
                   weight * (TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal(i,j,k+1)-Scal(i,j,k)) / dzmin
      end do
     end do
    end do
    !$omp end do nowait
    !$omp do
    do probe = 1, size(px)
        i = px(probe)
        j = py(probe)
        k = pz(probe)
        probes_flux(3, probe) = probes_flux(3, probe) + (TDiff(i,j,k+1)+TDiff(i,j,k)) * (Scal(i,j,k+1)-Scal(i,j,k)) / dzmin
    end do
    !$omp end do nowait
    !$omp end parallel
  endsubroutine AddScalarDiffVector
  

  subroutine ComputeTDiff(U, V, W)
    real(knd), contiguous, intent(in) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    integer:: i,j,k
    real(knd), parameter :: Pr_sgs = constPr_sgs !if variable, then implement as an internal function due to problems with inlining

    if (molecular_viscosity > 0) then
      !$omp parallel do private(i,j,k)
      do k = 1, Prnz
        do j = 1, Prny
          do i = 1, Prnx
            TDiff(i,j,k) = (Viscosity(i,j,k) - molecular_viscosity) / Pr_sgs + molecular_diffusivity
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(i,j,k)
      do k = 1, Prnz
        do j = 1, Prny
          do i = 1, Prnx
            TDiff(i,j,k) = Viscosity(i,j,k) / Pr_sgs
          end do
        end do
      end do
      !$omp end parallel do
    end if
  end subroutine ComputeTDiff


  subroutine ScalarDiffusion_Deallocate
    if (allocated(Ap)) deallocate(Ap)
    if (allocated(Scal3)) deallocate(Scal3)
  end subroutine

end module ScalarDiffusion
