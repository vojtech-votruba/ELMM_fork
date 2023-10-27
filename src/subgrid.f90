module Subgrid

  use Parameters
  use Subgrid_MixedTimeScale, only: SGS_MixedTimeScale, C_MixedTimeScale

  implicit none

  private
  public :: sgstype, SubgridModel, TKEDissipation, C_MixedTimeScale

!   real(knd),parameter :: CSmag = 0.122_knd

  integer :: sgstype

  integer, parameter, public :: SmagorinskyModel = 1, SigmaModel = 2, VremanModel = 3, &
                                StabSubgridModel = 4, MixedTimeScaleModel = 5, WALE = 6

  real(knd), public :: C_Smagorinsky = 0.122_knd, &
                       C_Sigma = 1.04_knd, &
                       C_Vreman = 0.041_knd, &
                       C_StabSubgrid = 1.04_knd, &
                       C_WALE = 0.58_knd
                     
  
  contains

    subroutine SGS_Smag(U,V,W,filter_ratio)  !Standard Smagorinsky model with implicit filtering
     use ArrayUtilities, only: add
     real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
     real(knd), intent(in) :: filter_ratio
     integer :: i,j,k

      !$omp parallel do private(i,j,k)
      do k = 1, Prnz
       do j = 1, Prny
        do i = 1, Prnx
          Viscosity(i,j,k) = NuSmag(i, j, k, U, V, W, filter_ratio)
        end do
       end do
      end do
      !$omp end parallel do
      if (molecular_viscosity > 0) then
        call add(Viscosity, molecular_viscosity)
      end if

    endsubroutine SGS_Smag




    real(knd) function NuSmag(i, j, k, U, V, W, filter_ratio)    !subgrid viscosity for Smagorinsky with implicit filtering
      integer, intent(in) :: i, j, k
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      real(knd), intent(in) :: filter_ratio
      real(knd) :: S(1:3,1:3)
      real(knd) :: width,Sbar
      real(knd) :: CSmag
      CSmag = C_Smagorinsky

      width = filter_ratio * (dxPr(i)*dyPr(j)*dzPr(k))**(1._knd/3._knd)
      call StrainIJ(i, j, k, U, V, W, S)
      Sbar = Strainu(S)
      NuSmag = Sbar*(width*CSmag)**2
    endfunction NuSmag





    subroutine SGS_StabSmag(U,V,W,Temperature,filter_ratio)  !Smagorinsky with a stability correction Brown et al. (1994)
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Temperature
      real(knd), intent(in) :: filter_ratio
      real(knd) :: Ri, l, l0
      real(knd) :: width, Sbar
      real(knd), parameter :: CS = 0.17_knd
      real(knd) :: S(1:3,1:3)
      integer :: i,j,k

      do k = 1, Prnz
       do j = 1, Prny
        do i = 1 , Prnx
          width = filter_ratio * (dxPr(i)*dyPr(j)*dzPr(k))**(1._knd/3._knd)

          call StrainIJ(i, j, k, U, V, W, S)
          Sbar = Strainu(S)

          Ri = Rig(i,j,k,U,V,temperature)


          l0 = CS*width
          l = WallDamp(l0,z0B,zPr(k))

          Viscosity(i,j,k)  = Sbar*Fm(Ri)*l
          TDiff(i,j,k) = Sbar*Fh(Ri)*l
        end do
       end do
      end do
      if (molecular_viscosity > 0) then
        Viscosity  = Viscosity  + molecular_viscosity
        TDiff = TDiff + molecular_diffusivity
      end if
    endsubroutine SGS_StabSmag


    pure real(knd) function Fm(Ri)       !Adjustment function for stability
      real(knd), intent(in) :: Ri           !Pointwise Richarson number
      real(knd), parameter :: Ric = 0.25_knd

      if (Ri>=Ric) then
       Fm  =  0
      elseif (Ri>0) then
       Fm  =  (1._knd-Ri/Ric)**4
      else
       Fm  =  sqrt(1._knd - 16._knd*Ri)
      end if
    end function Fm

    pure real(knd) function Fh(Ri)       !Adjustment function for stability
      real(knd), intent(in) :: Ri           !Pointwise Richarson number
      real(knd), parameter :: Ric = 0.25_knd

      if (Ri>=Ric) then
       Fh  =  0
      elseif (Ri>0) then
       Fh  =  (1./0.7_knd)  *  (1._knd - Ri/Ric)**4  *  (1._knd - 1.2_knd*Ri)
      else
       Fh  =  sqrt(1._knd - 40._knd*Ri) / 0.7_knd
      end if
    end function Fh

    pure real(knd) function WallDamp(l0, z0, z) !Wall damping for Smagorinsky model
      real(knd), intent(in) :: l0, z0, z

      WallDamp= 1._knd / (&
                    1._knd/l0**2  +  1._knd/(0.4_knd*(z+z0))**2 &
                    )
    end function WallDamp


    pure real(knd) function Rig(i, j, k, U, V, temperature)
      integer, intent(in) :: i, j, k
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: temperature
      real(knd) :: num, denom

      num = (grav_acc/temperature_ref) * &
            (temperature(i,j,k+1)-temperature(i,j,k-1)) / (zPr(k+1)-zPr(k-1))
      
      denom = ((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) / (2._knd*(zPr(k+1)-zPr(k-1))))**2
      
      denom = denom + &
              ((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) / (2._knd*(zPr(k+1)-zPr(k-1))))**2
      
      if (abs(denom)>tiny(1._knd)*100) then
       Rig = num/denom
      else
       Rig = 0
      end if
    endfunction Rig






    pure real(knd) function Strainu(S)      !Magnitude of the strain rate tensor.
      real(knd), intent(in) :: S(1:3,1:3)
      integer ::ii,jj

      Strainu = 0
      do jj = 1, 3
       do ii = 1, 3
        Strainu = Strainu + S(ii,jj)*S(ii,jj)
       end do
      end do
      Strainu = SQRT(2._knd*Strainu)
    endfunction Strainu



     subroutine StrainIJ(i, j, k, U, V, W, S)     !Computes components of the strain rate tensor.
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      real(knd), intent(out) :: S(1:3,1:3)
      integer, intent(in) :: i, j, k
      real(knd) :: D(1:3,1:3)
      integer :: ii, jj

      D = 0

      D(1,1) = (U(i,j,k)-U(i-1,j,k)) / dxmin
      D(2,2) = (V(i,j,k)-V(i,j-1,k)) / dymin
      D(3,3) = (W(i,j,k)-W(i,j,k-1)) / dzPr(k)
      D(1,2) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k)) / (4 * dymin)
      D(1,3) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) / (2 * (zPr(k+1)-zPr(k-1)))
      D(2,1) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k)) / (4 * dxmin)
      D(2,3) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) / (2 * (zPr(k+1)-zPr(k-1)))
      D(3,1) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1)) / (4 * dxmin)
      D(3,2) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) / (4 * dymin)

      do jj = 1, 3
       do ii = 1, 3
        S(ii,jj) = (D(ii,jj) + D(jj,ii)) / 2
       end do
      end do
    endsubroutine StrainIJ
    
    
    
    
     function TKEDissipation(i, j, k, U, V, W) result(res)
      !includes resolved and subgrid
      ! epsilon =  2*nu*Sij*Sij where Sij = (di_uj + dj_ui)/2
      real(knd) :: res
      integer, intent(in) :: i,j,k
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
      real(knd) :: S(1:3,1:3)
      integer :: ii, jj

      call StrainIJ(i, j, k, U, V, W, S)
      res = 0
      do jj = 1, 3
        do ii = 1, 3
          res = res + S(ii,jj)*S(ii,jj)
        end do
      end do
      res = 2 * Viscosity(i,j,k) * res
    end function




    subroutine SGS_Vreman(U,V,W,filter_ratio)   !Vreman subgrid model (Physics of Fluids, 2004)
      use Tiling, only: tilenx, tileny, tilenz
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      real(knd), intent(in) :: filter_ratio
      real(knd) :: aa, bb
      real(knd), dimension(1:3,1:3) :: a, b
      real(knd) :: dx2, dy2, dz2, dz_2k
!       real(knd),parameter ::c = 0.041
      real(knd) :: c2
      integer :: i, j, k, bi, bj, bk, ii, jj
      integer :: tnx, tny, tnz
      
      real(knd), parameter :: C1 = 9._knd / 8, C3 = 1._knd / (8*3)
      real(knd), parameter :: D1 = 2._knd / 3, D3 = 1._knd / 12

      integer,parameter :: narr = 4
      real(knd) :: c
      c = C_Vreman
    
      tnx = tilenx(narr)
      tny = tileny(narr)
      tnz = tilenz(narr)

      c2 = c * filter_ratio**2

      dx2 = dxmin**2
      dy2 = dymin**2
      dz2 = dzmin**2

      !NOTE: Consider consolidating the common part, but it would be branching in a tight loop.
      ! But there might be enough computing to offset that.
      if (gridtype==GRID_VARIABLE_Z) then
        !$omp parallel do private(aa,bb,a,b,i,j,k,bi,bj,bk,ii,jj) schedule(runtime) collapse(3)
        do bk = 1, Prnz, tnz
         do bj = 1, Prny, tny
          do bi = 1, Prnx, tnx
           do k = bk, min(bk+tnz-1, Prnz)
            do j = bj, min(bj+tny-1, Prny)
             do i = bi, min(bi+tnx-1, Prnx)
            
                dz_2k = 2 * (zPr(k+1) - zPr(k-1))
                dz2 = dzPr(k)**2
                
                
                a(1,1) = (U(i,j,k)-U(i-1,j,k)) / dxmin
                a(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k)) / (4._knd*dymin)
                a(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) / dz_2k

                a(2,2) = (V(i,j,k)-V(i,j-1,k)) / dymin
                a(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k)) / (4._knd*dxmin)
                a(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) / dz_2k

                a(3,3) = (W(i,j,k)-W(i,j,k-1)) / dzPr(k)
                a(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1)) / (4._knd*dxmin)
                a(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) / (4._knd*dymin)

                do jj = 1, 3
                  do ii = 1, 3
                    b(ii,jj) = dx2 * a(1,ii) * a(1,jj) + &
                               dy2 * a(2,ii) * a(2,jj) + &
                               dz2 * a(3,ii) * a(3,jj)
                  end do
                end do

                bb =      b(1,1)*b(2,2) - b(1,2)**2
                bb = bb + b(1,1)*b(3,3) - b(1,3)**2
                bb = bb + b(2,2)*b(3,3) - b(2,3)**2

                aa = 0

                do jj = 1, 3
                  do ii = 1, 3
                    aa = aa + (a(ii,jj)**2)
                  end do
                end do



                Viscosity(i,j,k) = c2 * sqrt(bb/aa)

                Viscosity(i,j,k) = max(0._knd, Viscosity(i,j,k)) + molecular_viscosity
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp end parallel do
      else
        if (discretization_order==4) then
          !$omp parallel do private(aa,bb,a,b,i,j,k,bi,bj,bk,ii,jj) schedule(runtime) collapse(3)
          do bk = 1, Prnz, tnz
           do bj = 1, Prny, tny
            do bi = 1, Prnx, tnx
             do k = bk, min(bk+tnz-1, Prnz)
              do j = bj, min(bj+tny-1, Prny)
               do i = bi, min(bi+tnx-1, Prnx)
                  a(1,1) = (C1*(U(i,j,k)-U(i-1,j,k)) - C3*(U(i+1,j,k)-U(i-2,j,k))) / dxmin
                  a(2,1) = ( &
                             D1*(U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k)) - &
                             D3*(U(i,j+2,k)+U(i-1,j+2,k)-U(i,j-2,k)-U(i-1,j-2,k)) &
                           ) / (2*dymin)
                  a(3,1) = ( &
                             D1*(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) - &
                             D1*(U(i,j,k+2)+U(i-1,j,k+2)-U(i,j,k-2)-U(i-1,j,k-2)) &
                           ) / (2*dzmin)

                  a(2,2) = (C1*(V(i,j,k)-V(i,j-1,k)) - C3*(V(i,j+1,k)-V(i,j-2,k))) / dymin
                  a(1,2) = ( &
                             D1*(V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k)) - &
                             D3*(V(i+2,j,k)+V(i+2,j-1,k)-V(i-2,j,k)-V(i-2,j-1,k)) &
                           ) / (2*dxmin)
                  a(3,2) = ( &
                             D1*(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) - &
                             D3*(V(i,j,k+2)+V(i,j-1,k+2)-V(i,j,k-2)-V(i,j-1,k-2)) &
                           ) / (2*dzmin)

                  a(3,3) = (C1*(W(i,j,k)-W(i,j,k-1)) - C3*(W(i,j,k)-W(i,j,k-1))) / dzmin
                  a(1,3) = ( &
                             D1*(W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1)) - &
                             D3*(W(i+2,j,k)+W(i+2,j,k-1)-W(i-2,j,k)-W(i-2,j,k-1)) &
                           ) / (2*dxmin)
                  a(2,3) = ( &
                             D1*(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) - &
                             D3*(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) &
                           ) / (2*dymin)

                  do jj = 1, 3
                    do ii = 1, 3
                      b(ii,jj) = dx2 * a(1,ii) * a(1,jj) + &
                                 dy2 * a(2,ii) * a(2,jj) + &
                                 dz2 * a(3,ii) * a(3,jj)
                    end do
                  end do

                  bb =      b(1,1)*b(2,2) - b(1,2)**2
                  bb = bb + b(1,1)*b(3,3) - b(1,3)**2
                  bb = bb + b(2,2)*b(3,3) - b(2,3)**2

                  aa = 0

                  do jj = 1, 3
                    do ii = 1, 3
                      aa = aa + (a(ii,jj)**2)
                    end do
                  end do



                  Viscosity(i,j,k) = c2 * sqrt(bb/aa)

                  Viscosity(i,j,k) = max(0._knd, Viscosity(i,j,k)) + molecular_viscosity
               end do
              end do
             end do
            end do
           end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(aa,bb,a,b,i,j,k,bi,bj,bk,ii,jj) schedule(runtime) collapse(3)
          do bk = 1, Prnz, tnz
           do bj = 1, Prny, tny
            do bi = 1, Prnx, tnx
             do k = bk, min(bk+tnz-1, Prnz)
              do j = bj, min(bj+tny-1, Prny)
               do i = bi, min(bi+tnx-1, Prnx)
                  a(1,1) = (U(i,j,k)-U(i-1,j,k)) / dxmin
                  a(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k)) / (4._knd*dymin)
                  a(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) / (4._knd*dzmin)

                  a(2,2) = (V(i,j,k)-V(i,j-1,k)) / dymin
                  a(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k)) / (4._knd*dxmin)
                  a(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) / (4._knd*dzmin)

                  a(3,3) = (W(i,j,k)-W(i,j,k-1)) / dzmin
                  a(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1)) / (4._knd*dxmin)
                  a(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) / (4._knd*dymin)

                  do jj = 1, 3
                    do ii = 1, 3
                      b(ii,jj) = dx2 * a(1,ii) * a(1,jj) + &
                                 dy2 * a(2,ii) * a(2,jj) + &
                                 dz2 * a(3,ii) * a(3,jj)
                    end do
                  end do

                  bb =      b(1,1)*b(2,2) - b(1,2)**2
                  bb = bb + b(1,1)*b(3,3) - b(1,3)**2
                  bb = bb + b(2,2)*b(3,3) - b(2,3)**2

                  aa = 0

                  do jj = 1, 3
                    do ii = 1, 3
                      aa = aa + (a(ii,jj)**2)
                    end do
                  end do



                  Viscosity(i,j,k) = c2 * sqrt(bb/aa)

                  Viscosity(i,j,k) = max(0._knd, Viscosity(i,j,k)) + molecular_viscosity
               end do
              end do
             end do
            end do
           end do
          end do
          !$omp end parallel do
        end if
      end if
    endsubroutine SGS_Vreman



    subroutine SGS_Sigma(U, V, W, filter_ratio)
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      real(knd), intent(in) :: filter_ratio
      if (gridtype==GRID_VARIABLE_Z) then
        call SGS_Sigma_variable_z(U, V, W, filter_ratio)
      else
        call SGS_Sigma_uniform(U, V, W, filter_ratio)
      end if
    end subroutine SGS_Sigma

    subroutine SGS_Sigma_uniform(U, V, W, filter_ratio)
      !from Nicoud, Toda, Cabrit, Bose, Lee, http://dx.doi.org/10.1063/1.3623274
      use Tiling, only: tilenx, tileny, tilenz
      use ieee_exceptions
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      real(knd), intent(in) :: filter_ratio
!       real(knd), parameter :: Csig = 1.04_knd
      integer   :: i, j, k, bi, bj, bk
      real(knd) :: width, C, D, g(3,3), s1, s2, s3
      integer, parameter :: sigma_knd = knd
      integer :: tnx, tny, tnz

      integer,parameter :: narr = 4
      
      logical :: saved_fpe_mode(size(ieee_all))
      
      real(knd) :: Csig
      Csig = C_Sigma
      
      call ieee_get_halting_mode(ieee_all, saved_fpe_mode)
      call ieee_set_halting_mode(ieee_all, .false.)
    
      tnx = tilenx(narr)
      tny = tileny(narr)
      tnz = tilenz(narr)

      width = filter_ratio * (dxmin*dymin*dzmin)**(1._knd/3._knd)
      C = (Csig*width)**2

      !TODO: profile and perhaps consolidate into a single loop.
      if (discretization_order==4) then
        !$omp parallel do private(g,s1,s2,s3,D,i,j,k,bi,bj,bk) schedule(runtime) collapse(3)
        do bk = 1, Prnz, tnz
         do bj = 1, Prny, tny
          do bi = 1, Prnx, tnx
           do k = bk, min(bk+tnz-1, Prnz)
            do j = bj, min(bj+tny-1, Prny)
             do i = bi, min(bi+tnx-1, Prnx)

              call GradientTensorUG4(g, i, j, k)

              call Sigmas(s1,s2,s3,g)

              D = (s3 * (s1 - s2) * (s2 - s3)) / s1**2

              D = max(0._knd, D)

              Viscosity(i,j,k) = C * D + molecular_viscosity

             end do
            end do
           end do
          end do
         end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do private(g,s1,s2,s3,D,i,j,k,bi,bj,bk) schedule(runtime) collapse(3)
        do bk = 1, Prnz, tnz
         do bj = 1, Prny, tny
          do bi = 1, Prnx, tnx
           do k = bk, min(bk+tnz-1, Prnz)
            do j = bj, min(bj+tny-1, Prny)
             do i = bi, min(bi+tnx-1, Prnx)

              call GradientTensorUG(g, i, j, k)

              call Sigmas(s1,s2,s3,g)

              D = (s3 * (s1 - s2) * (s2 - s3)) / s1**2

              D = max(0._knd, D)

              Viscosity(i,j,k) = C * D + molecular_viscosity

             end do
            end do
           end do
          end do
         end do
        end do
        !$omp end parallel do
      end if
      
      call ieee_set_halting_mode(ieee_all, saved_fpe_mode)
    contains

      pure subroutine GradientTensorUG(g ,i, j, k)
        real(knd), intent(out) :: g(3,3)
        integer, intent(in) :: i, j, k

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


      pure subroutine GradientTensorUG4(g ,i, j, k)
        real(knd), intent(out) :: g(3,3)
        integer, intent(in) :: i, j, k
        real(knd), parameter :: C1 = 9._knd / 8, C3 = 1._knd / (8*3)
        real(knd), parameter :: D1 = 2._knd / 3, D3 = 1._knd / 12

        g(1,1) = (C1*(U(i,j,k)-U(i-1,j,k)) - C3*(U(i+1,j,k)-U(i-2,j,k))) / dxmin
        g(2,1) = ( &
                    D1*(U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k)) - &
                    D3*(U(i,j+2,k)+U(i-1,j+2,k)-U(i,j-2,k)-U(i-1,j-2,k)) &
                 ) / (2*dymin)
        g(3,1) = ( &
                    D1*(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) - &
                    D1*(U(i,j,k+2)+U(i-1,j,k+2)-U(i,j,k-2)-U(i-1,j,k-2)) &
                 ) / (2*dzmin)

        g(2,2) = (C1*(V(i,j,k)-V(i,j-1,k)) - C3*(V(i,j+1,k)-V(i,j-2,k))) / dymin
        g(1,2) = ( &
                    D1*(V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k)) - &
                    D3*(V(i+2,j,k)+V(i+2,j-1,k)-V(i-2,j,k)-V(i-2,j-1,k)) &
                 ) / (2*dxmin)
        g(3,2) = ( &
                    D1*(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) - &
                    D3*(V(i,j,k+2)+V(i,j-1,k+2)-V(i,j,k-2)-V(i,j-1,k-2)) &
                 ) / (2*dzmin)

        g(3,3) = (C1*(W(i,j,k)-W(i,j,k-1)) - C3*(W(i,j,k)-W(i,j,k-1))) / dzmin
        g(1,3) = ( &
                    D1*(W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1)) - &
                    D3*(W(i+2,j,k)+W(i+2,j,k-1)-W(i-2,j,k)-W(i-2,j,k-1)) &
                 ) / (2*dxmin)
        g(2,3) = ( &
                    D1*(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) - &
                    D3*(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) &
                 ) / (2*dymin)
      end subroutine GradientTensorUG4


      pure subroutine Sigmas(s1,s2,s3,grads)
        !from Hasan, Basser, Parker, Alexander, http://dx.doi.org/10.1006/jmre.2001.2400
        !via Nicoud, Toda, Cabrit, Bose, Lee, http://dx.doi.org/10.1063/1.3623274
        use ieee_arithmetic

        real(knd), intent(out) :: s1, s2, s3
        real(knd), intent(in)  :: grads(3,3)

        real(sigma_knd) :: trG2, i1, i2, i3, a1, a2, a3, c, G(3,3)

        !G = matmul(transpose(grads),grads)
        G(1,1) = grads(1,1)**2           + grads(2,1)**2           + grads(3,1)**2
        G(1,2) = grads(1,1) * grads(1,2) + grads(2,1) * grads(2,2) + grads(3,1) * grads(3,2)
        G(1,3) = grads(1,1) * grads(1,3) + grads(2,1) * grads(2,3) + grads(3,1) * grads(3,3)
        G(2,1) = grads(1,1) * grads(1,2) + grads(2,1) * grads(2,2) + grads(3,1) * grads(3,2)
        G(2,2) = grads(1,2)**2           + grads(2,2)**2           + grads(3,2)**2
        G(2,3) = grads(1,2) * grads(1,3) + grads(2,2) * grads(2,3) + grads(3,2) * grads(3,3)
        G(3,1) = grads(1,1) * grads(1,3) + grads(2,1) * grads(2,3) + grads(3,1) * grads(3,3)
        G(3,2) = grads(1,2) * grads(1,3) + grads(2,2) * grads(2,3) + grads(3,2) * grads(3,3)
        G(3,3) = grads(1,3)**2           + grads(2,3)**2           + grads(3,3)**2

        trG2 = dot_product(G(:,1),G(:,1)) +&
               dot_product(G(:,2),G(:,2)) +&
               dot_product(G(:,3),G(:,3))

        i1 = G(1,1) + G(2,2) + G(3,3)

        i2 = (i1**2 - trG2)
        i2 = i2 / 2

        i3 = det3x3(G)

        a1 = max((i1**2)/9 - i2/3,0._sigma_knd)

        a2 = (i1**3)/27 - i1*i2/6 +i3/2

        !This requires no FPE trapping is in progress!
        c =  a2 / sqrt(a1**3)

        !If c is NaN let it be 1.
        c = max(-1._sigma_knd, min(1._sigma_knd,c))
        a3 = acos(c) / 3

        c = 2*sqrt(a1)

        s1 = real( sqrt( i1/3 + c*cos(a3) ) , knd )

        s2 = real( sqrt( max( i1/3 - c*cos(pi/3 + a3) , 0._sigma_knd ) ) , knd)

        s3 = real( sqrt( max( i1/3 - c*cos(pi/3 - a3) , 0._sigma_knd ) ) , knd)

      end subroutine Sigmas


      pure function det3x3(A) result (res)
        real(sigma_knd), intent(in) :: A(3,3)
        real(sigma_knd) :: res

        res =   A(1,1) * A(2,2) * A(3,3)  &
              - A(1,1) * A(2,3) * A(3,2)  &
              - A(1,2) * A(2,1) * A(3,3)  &
              + A(1,2) * A(2,3) * A(3,1)  &
              + A(1,3) * A(2,1) * A(3,2)  &
              - A(1,3) * A(2,2) * A(3,1)
      end function det3x3

    end subroutine SGS_Sigma_uniform
    
    
    
    
    subroutine SGS_Sigma_variable_z(U, V, W, filter_ratio)
      !from Nicoud, Toda, Cabrit, Bose, Lee, http://dx.doi.org/10.1063/1.3623274
      use Tiling, only: tilenx, tileny, tilenz
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      real(knd), intent(in) :: filter_ratio
!       real(knd), parameter :: Csig = 1.04_knd
      integer   :: i, j, k, bi, bj, bk
      real(knd) :: width, dz23, C, D, g(3,3), s1, s2, s3
      integer, parameter :: sigma_knd = knd
      integer :: tnx, tny, tnz

      integer,parameter :: narr = 4
      real(knd) :: Csig
      Csig = C_Sigma
    
      tnx = tilenx(narr)
      tny = tileny(narr)
      tnz = tilenz(narr)

      width = filter_ratio * (dxmin*dymin)**(1._knd/3._knd)
      C = (Csig*width)**2

      !$omp parallel do private(g,s1,s2,s3,D,i,j,k,bi,bj,bk) schedule(runtime) collapse(3)
      do bk = 1, Prnz, tnz
       do bj = 1, Prny, tny
        do bi = 1, Prnx, tnx
         do k = bk, min(bk+tnz-1, Prnz)
          dz23 = dzPr(k)**(2._knd/3._knd)
          do j = bj, min(bj+tny-1, Prny)
           do i = bi, min(bi+tnx-1, Prnx)

            call GradientTensorUG(g, i, j, k)

            call Sigmas(s1,s2,s3,g)

            D = (s3 * (s1 - s2) * (s2 - s3)) / s1**2

            D = max(0._knd, D)

            Viscosity(i,j,k) = (C*dz23) * D + molecular_viscosity

           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end parallel do

    contains

      pure subroutine GradientTensorUG(g ,i, j, k)
        real(knd), intent(out) :: g(3,3)
        integer, intent(in) :: i, j, k

        g(1,1) = (U(i,j,k)-U(i-1,j,k)) / dxmin
        g(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k)) / (4._knd*dymin)
        g(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) / (2*(zPr(k+1)-zPr(k-1)))

        g(2,2) = (V(i,j,k)-V(i,j-1,k)) / dymin
        g(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k)) / (4._knd*dxmin)
        g(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) / (2*(zPr(k+1)-zPr(k-1)))

        g(3,3) = (W(i,j,k)-W(i,j,k-1)) / dzPr(k)
        g(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1)) / (4._knd*dxmin)
        g(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) / (4._knd*dymin)
      end subroutine GradientTensorUG


      pure subroutine Sigmas(s1,s2,s3,grads)
        !from Hasan, Basser, Parker, Alexander, http://dx.doi.org/10.1006/jmre.2001.2400
        !via Nicoud, Toda, Cabrit, Bose, Lee, http://dx.doi.org/10.1063/1.3623274
        use ieee_arithmetic

        real(knd), intent(out) :: s1, s2, s3
        real(knd), intent(in)  :: grads(3,3)

        real(sigma_knd) :: trG2, i1, i2, i3, a1, a2, a3, c, G(3,3)

        !G = matmul(transpose(grads),grads)
        G(1,1) = grads(1,1)**2           + grads(2,1)**2           + grads(3,1)**2
        G(1,2) = grads(1,1) * grads(1,2) + grads(2,1) * grads(2,2) + grads(3,1) * grads(3,2)
        G(1,3) = grads(1,1) * grads(1,3) + grads(2,1) * grads(2,3) + grads(3,1) * grads(3,3)
        G(2,1) = grads(1,1) * grads(1,2) + grads(2,1) * grads(2,2) + grads(3,1) * grads(3,2)
        G(2,2) = grads(1,2)**2           + grads(2,2)**2           + grads(3,2)**2
        G(2,3) = grads(1,2) * grads(1,3) + grads(2,2) * grads(2,3) + grads(3,2) * grads(3,3)
        G(3,1) = grads(1,1) * grads(1,3) + grads(2,1) * grads(2,3) + grads(3,1) * grads(3,3)
        G(3,2) = grads(1,2) * grads(1,3) + grads(2,2) * grads(2,3) + grads(3,2) * grads(3,3)
        G(3,3) = grads(1,3)**2           + grads(2,3)**2           + grads(3,3)**2

        trG2 = dot_product(G(:,1),G(:,1)) +&
               dot_product(G(:,2),G(:,2)) +&
               dot_product(G(:,3),G(:,3))

        i1 = G(1,1) + G(2,2) + G(3,3)

        i2 = (i1**2 - trG2)
        i2 = i2 / 2

        i3 = det3x3(G)

        a1 = max((i1**2)/9 - i2/3,0._sigma_knd)

        a2 = (i1**3)/27 - i1*i2/6 +i3/2

        !This requires no FPE trapping is in progress!
        c =  a2 / sqrt(a1**3)

        !If c is NaN let it be 1.
        c = max(-1._sigma_knd, min(1._sigma_knd,c))
        a3 = acos(c) / 3

        c = 2*sqrt(a1)

        s1 = real( sqrt( i1/3 + c*cos(a3) ) , knd )

        s2 = real( sqrt( max( i1/3 - c*cos(pi/3 + a3) , 0._sigma_knd ) ) , knd)

        s3 = real( sqrt( max( i1/3 - c*cos(pi/3 - a3) , 0._sigma_knd ) ) , knd)

      end subroutine Sigmas


      pure function det3x3(A) result (res)
        real(sigma_knd), intent(in) :: A(3,3)
        real(sigma_knd) :: res

        res =   A(1,1) * A(2,2) * A(3,3)  &
              - A(1,1) * A(2,3) * A(3,2)  &
              - A(1,2) * A(2,1) * A(3,3)  &
              + A(1,2) * A(2,3) * A(3,1)  &
              + A(1,3) * A(2,1) * A(3,2)  &
              - A(1,3) * A(2,2) * A(3,1)
      end function det3x3

    end subroutine SGS_Sigma_variable_z
    
    
    
    
    subroutine SGS_Sigma_stability(U,V,W,filter_ratio)
      !from Nicoud, Toda, Cabrit, Bose, Lee, http://dx.doi.org/10.1063/1.3623274
#ifdef PAR
      use custom_par
#endif
      use Tiling, only: tilenx, tileny, tilenz
      use VerticalProfiles, only: enable_profiles, current_profiles
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      real(knd), intent(in) :: filter_ratio
!       real(knd), parameter :: Csig = 1.04_knd
      integer   :: i, j, k, bi, bj, bk
      real(knd) :: width, C, Pr_sgs, D, g(3,3), s1, s2, s3
      integer, parameter :: sigma_knd = knd

      real(knd) :: tempfl_p(0:Prnz), momfl_p(0:Prnz)
      
      logical :: enable_stability_correction
      integer :: tnx, tny, tnz

      integer,parameter :: narr = 4
      real(knd) :: Csig
      Csig = C_StabSubgrid
    
      tnx = tilenx(narr)
      tny = tileny(narr)
      tnz = tilenz(narr)
      
      enable_stability_correction = enable_buoyancy .and. enable_profiles

      if (enable_stability_correction) then
        tempfl_p = current_profiles%tempfl + current_profiles%tempflsgs
        momfl_p = sqrt((current_profiles%uw + current_profiles%uwsgs)**2 + &
                       (current_profiles%vw + current_profiles%vwsgs)**2)
#ifdef PAR
        tempfl_p = par_co_sum(tempfl_p, comm=comm_plane_xy)
        momfl_p = par_co_sum(momfl_p, comm=comm_plane_xy)
#endif
        tempfl_p = tempfl_p / (gPrnx*gPrny)
        momfl_p = momfl_p / (gPrnx*gPrny)
      end if      

      width = filter_ratio * (dxmin*dymin*dzmin)**(1._knd/3._knd)

      !$omp parallel do private(g,s1,s2,s3,D,i,j,k,bi,bj,bk,C,Pr_sgs) schedule(runtime) collapse(3)
      do bk = 1, Prnz, tnz
       do bj = 1, Prny, tny
        do bi = 1, Prnx, tnx
         do k = bk, min(bk+tnz-1, Prnz)
         
          if (enable_stability_correction) then
            call stability_correction_L_MO(k, width, C, Pr_sgs)
            C = (C * Csig * width)**2
          else
            C = (Csig*width)**2
          end if
            
          do j = bj, min(bj+tny-1, Prny)
           do i = bi, min(bi+tnx-1, Prnx)

            call GradientTensorUG(g, i, j, k)

            call Sigmas(s1,s2,s3,g)

            D = (s3 * (s1 - s2) * (s2 - s3)) / s1**2

            D = max(0._knd, D)

            Viscosity(i,j,k) = C * D

            TDiff(i,j,k) = Viscosity(i,j,k) / Pr_sgs + molecular_diffusivity
            Viscosity(i,j,k) = Viscosity(i,j,k) + molecular_viscosity

           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end parallel do

    contains

      pure subroutine GradientTensorUG(g, i, j, k)
        real(knd), intent(out) :: g(3,3)
        integer, intent(in) :: i, j, k

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


      pure subroutine Sigmas(s1, s2, s3, grads)
        !from Hasan, Basser, Parker, Alexander, http://dx.doi.org/10.1006/jmre.2001.2400
        !via Nicoud, Toda, Cabrit, Bose, Lee, http://dx.doi.org/10.1063/1.3623274
        use ieee_arithmetic

        real(knd), intent(out) :: s1, s2, s3
        real(knd), intent(in)  :: grads(3,3)

        real(sigma_knd) :: trG2, i1, i2, i3, a1, a2, a3, c, G(3,3)

        G = real( matmul(transpose(grads),grads) ,sigma_knd)

        trG2 = dot_product(G(:,1),G(:,1)) +&
               dot_product(G(:,2),G(:,2)) +&
               dot_product(G(:,3),G(:,3))

        i1 = G(1,1) + G(2,2) + G(3,3)

        i2 = (i1**2 - trG2)
        i2 = i2 / 2

        i3 = det3x3(G)

        a1 = max((i1**2)/9 - i2/3,0._sigma_knd)

        a2 = (i1**3)/27 - i1*i2/6 +i3/2

        !This requires no FPE trapping is in progress!
        c =  a2 / sqrt(a1**3)

        !If c is NaN let it be 1.
        c = max(-1._sigma_knd, min(1._sigma_knd,c))
        a3 = acos(c) / 3

        c = 2*sqrt(a1)

        s1 = real( sqrt( i1/3 + c*cos(a3) ) , knd )

        s2 = real( sqrt( max( i1/3 - c*cos(pi/3 + a3) , 0._sigma_knd ) ) , knd)

        s3 = real( sqrt( max( i1/3 - c*cos(pi/3 - a3) , 0._sigma_knd ) ) , knd)

      end subroutine Sigmas


      pure function det3x3(A) result (res)
        real(sigma_knd), intent(in) :: A(3,3)
        real(sigma_knd) :: res

        res =   A(1,1) * A(2,2) * A(3,3)  &
              - A(1,1) * A(2,3) * A(3,2)  &
              - A(1,2) * A(2,1) * A(3,3)  &
              + A(1,2) * A(2,3) * A(3,1)  &
              + A(1,3) * A(2,1) * A(3,2)  &
              - A(1,3) * A(2,2) * A(3,1)
      end function det3x3
      
      subroutine stability_correction_L_MO(k, width, C, Pr_sgs)
        real(knd), intent(out) :: C, Pr_sgs
        integer, intent(in) :: k
        real(knd), intent(in) :: width
        real(knd) :: tempfl, momfl, L, wL
        
        tempfl = (tempfl_p(k) + tempfl_p(k-1))/2
        momfl = (momfl_p(k) + momfl_p(k-1)) / 2
        if (tempfl < 0 .and. momfl < 0) then
          L = - momfl**(3._knd/2._knd) * temperature_ref / &
               ( 0.4_knd* grav_acc * tempfl)
          wL = width / L
          C = 1 / (1 + wL)
          !own fit of graph in Zeid et al. 2010, JFM 665
          Pr_sgs = 0.5 + 0.85 * (1+tanh(log(1.2_knd*wL)))/2
        else
          C = 1
          Pr_sgs = 0.5
        end if
      end subroutine

    end subroutine SGS_Sigma_stability

    subroutine SGS_WALE(U,V,W, filter_ratio)
      !WALE subgrid model - Flow, Turbul. Combust, Nicoud F., Ducros F. (1999)
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      real(knd), intent(in) :: filter_ratio
      integer :: i, j, k
      integer :: ii, jj, kk, ll
      real (knd) :: width
      real(knd) :: S(1:3,1:3)
      real(knd) :: Omega(1:3,1:3)
      real (knd) :: S_sqr, Omega_sqr, IVs
      real (knd) :: SdSd, OP1, OP2


      width = filter_ratio * (dxmin*dymin*dzmin)**(1._knd/3._knd)

      do k = 1, Prnz
        do j = 1, Prny
          do i = 1, Prnx
            call StrainIJ(i, j, k, U, V, W, S)
            call OmegaIJ(Omega, i, j, k)
            S_sqr = 0
            Omega_sqr = 0
            IVs = 0

            do ii = 1, 3
              do jj = 1, 3
                S_sqr = S_sqr + S(ii, jj)*S(ii, jj)
                Omega_sqr = Omega_sqr + Omega(ii, jj) * Omega(ii, jj)
              end do
            end do

            do ii = 1, 3
              do jj = 1, 3
                do kk = 1, 3
                  do ll = 1, 3
                    IVs = IVs + S(ii, kk) * S(kk, jj) * Omega(jj, ll) * Omega(ll, ii)
                  end do
                end do
              end do
            end do

            SdSd = (1._knd/6._knd)*(S_sqr*S_sqr + Omega_sqr*Omega_sqr) + (2._knd/3._knd)*(S_sqr * Omega_sqr) + 2._knd*IVs
            SdSd = max(SdSd, 0._knd)

            OP1 = (SdSd)**(3._knd/2._knd)
            OP2 = (S_sqr)**(5._knd/2._knd) + (SdSd)**(5._knd/4._knd)
            if (OP1 == 0._knd .AND. OP2 == 0._knd) then
              Viscosity(i,j,k) = 0._knd
            else
              Viscosity(i,j,k) = (C_WALE * width)**(2._knd) * OP1/OP2
            endif
            Viscosity(i,j,k) = Viscosity(i,j,k) + molecular_viscosity
          end do
        end do
      end do

      contains
        pure subroutine OmegaIJ(Omega, i, j, k)
          real(knd), intent(out) :: Omega(1:3,1:3)
          integer, intent(in) :: i, j, k
          real(knd) :: D(1:3,1:3)
          integer :: ii, jj

          D = 0

          D(1,1) = (U(i,j,k)-U(i-1,j,k)) / dxmin
          D(2,2) = (V(i,j,k)-V(i,j-1,k)) / dymin
          D(3,3) = (W(i,j,k)-W(i,j,k-1)) / dzPr(k)
          D(1,2) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k)) / (4 * dymin)
          D(1,3) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) / (2 * (zPr(k+1)-zPr(k-1)))
          D(2,1) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k)) / (4 * dxmin)
          D(2,3) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) / (2 * (zPr(k+1)-zPr(k-1)))
          D(3,1) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1)) / (4 * dxmin)
          D(3,2) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)) / (4 * dymin)

          do jj = 1, 3
            do ii = 1, 3
              Omega(ii,jj) = (D(ii,jj) - D(jj,ii)) / 2
            end do
          end do

        end subroutine OmegaIJ

    endsubroutine SGS_WALE
    
    subroutine SubgridModel(U, V, W)
      !dispatch of subgrid models
      use Filters, only: filtertype, filter_ratios
      real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
      
      if (sgstype==SmagorinskyModel) then
                        call SGS_Smag(U, V, W, filter_ratios(filtertype))
      else if (sgstype==SigmaModel) then
                        call SGS_Sigma(U, V, W, filter_ratios(filtertype))
      else if (sgstype==VremanModel) then
                        call SGS_Vreman(U, V, W, filter_ratios(filtertype))
      else if (sgstype==StabSubgridModel) then
                        call SGS_Sigma_stability(U, V, W, filter_ratios(filtertype))
      else if (sgstype==MixedTimeScaleModel) then
                        call SGS_MixedTimeScale(U, V, W)
      else if (sgstype==WALE) then
                        call SGS_WALE(U,V,W, filter_ratios(filtertype))
      else

        Viscosity = molecular_viscosity

      end if

    end subroutine

end module Subgrid
