module Sponge

  use Parameters
#ifdef PAR
    use custom_par
#endif

  implicit none
  
  private
  
  public :: enable_top_sponge, enable_top_sponge_scalar, &
            enable_in_sponge_x, enable_out_sponge_x, enable_out_sponge_y, &
            SpongeTop, SpongeOut, SpongeTopScalar, &
            top_sponge_bottom, sponge_to_profiles, &
            U_sponge_avg, V_sponge_avg, W_sponge_avg

  logical :: enable_top_sponge = .false.
  logical :: enable_top_sponge_scalar = .false.
  logical :: enable_in_sponge_x = .false.
  logical :: enable_out_sponge_x = .false.
  logical :: enable_out_sponge_y = .false.

  logical :: sponge_to_profiles = .false.

  real(knd), dimension(:), allocatable :: U_sponge_avg, V_sponge_avg, W_sponge_avg

  real(knd) :: top_sponge_bottom = huge(1._knd)

  real(knd), dimension(:), allocatable :: DF, avg
  
contains

  elemental function DampF(x) result(res)
    real(knd) :: res
    real(knd), intent(in) :: x

    if (x<=0) then
      res = 1
    else if (x>=1) then
      res = 0
    else
      res = (1 - 0.1_knd*x**2) * &
              ( 1 - (1 - exp(10._knd*x**2)) / (1 - exp(10._knd)) )
    end if
  end function
  
  elemental function ScalarDampF(x) result(res)
    real(knd) :: res
    real(knd), intent(in)::x

    if (x <= 0) then
      res = 1
    else if (x >= 1) then
      res = 0
    else
      res = (1 - 0.1_knd*x**2) * &
              ( 1 - (1 - exp(10._knd*x**2)) / (1 - exp(10._knd)) )
    end if
  end function
  

  subroutine SpongeTop(U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    integer   :: i, j, k, lo_U, lo_V, lo_W
    real(knd) :: ze, zs, zb, p
    
    if (top_sponge_bottom<zW(Prnz)) then
      lo_U = Unz - min(Unz, int((zW(Prnz) - top_sponge_bottom) / dzmin) ) + 1
      lo_V = Vnz - min(Vnz, int((zW(Prnz) - top_sponge_bottom) / dzmin) ) + 1
      lo_W = Wnz - min(Wnz, int((zW(Prnz) - top_sponge_bottom) / dzmin) ) + 1
      zs = top_sponge_bottom
      ze = gzmax
   
      if (.not.allocated(DF)) allocate(DF(  min(lo_U, lo_V, lo_W)  :  max(Unz, Vnz, Wnz)))
      if (.not.allocated(avg)) allocate(avg( min(lo_U, lo_V, lo_W)  :  max(Unz, Vnz, Wnz)))

      !$omp parallel private(i, j, k, p, zb)
      
      if (sponge_to_profiles) then
        avg(lo_U:Unz) = U_sponge_avg(lo_U:Unz)
      else
        !$omp do
        do k = lo_U, Unz
          avg(k) = 0
        end do
        
        !$omp do
        do k = lo_U, Unz
          p = 0

          do j = 1, Uny
            do i = 1, Unx
              p = p + U(i,j,k)
            end do
          end do
          avg(k) = p
        end do
        
#ifdef PAR
        avg = par_co_sum(avg, comm_plane_xy)
#endif

        !$omp do
        do k = lo_U, Unz
          avg(k) = avg(k) / (gUnx*gUny)
        end do
      end if

      !$omp do
      do k = lo_U, Unz
        zb = (zPr(k)-zs) / (ze-zs)
        DF(k) = DampF(zb)
      end do

      !$omp do
      do k = lo_U, Unz
        do j = -1, Uny + 1
          do i = -1, Unx + 1
            U(i,j,k) = avg(k) + DF(k) * (U(i,j,k) - avg(k))
          end do
        end do
      end do


      if (sponge_to_profiles) then
        avg(lo_V:Vnz) = V_sponge_avg(lo_V:Vnz)
      else
        !$omp do
        do k = lo_V, Vnz
          avg(k) = 0
        end do

        !$omp do
        do k = lo_V, Vnz
          p = 0

          do j = 1, Vny
            do i = 1, Vnx
              p = p + V(i,j,k)
            end do
          end do
          avg(k) = p
        end do

#ifdef PAR
        avg = par_co_sum(avg, comm_plane_xy)
#endif

        !$omp do
        do k = lo_V, Vnz
          avg(k) = avg(k) / (gVnx*gVny)
        end do
      end if

      !$omp do
      do k = lo_V, Vnz
        zb = (zPr(k)-zs) / (ze-zs)
        DF(k) = DampF(zb)
      end do

      !$omp do
      do k = lo_V, Vnz
        do j = -1, Vny + 1
          do i = -1, Vnx + 1
            V(i,j,k) = avg(k) + DF(k) * (V(i,j,k) - avg(k))
          end do
        end do
      end do


      if (sponge_to_profiles) then
        avg(lo_W:Wnz) = W_sponge_avg(lo_W:Wnz)
      else
        !$omp do
        do k = lo_W, Wnz
          avg(k) = 0
        end do

        !$omp do
        do k = lo_W, Wnz
          p = 0

          do j = 1, Wny
            do i = 1, Wnx
              p = p + W(i,j,k)
            end do
          end do
          avg(k) = p
        end do

#ifdef PAR
        avg = par_co_sum(avg, comm_plane_xy)
#endif

        !$omp do
        do k = lo_W, Wnz
          avg(k) = avg(k) / (gWnx*gWny)
        end do
      end if

      !$omp do
      do k = lo_W, Wnz
        zb = (zW(k)-zs) / (ze-zs)
        DF(k) = DampF(zb)
      end do

      !$omp do
      do k = lo_W, Wnz
        do j = -1, Wny + 1
          do i = -1, Wnx + 1
            W(i,j,k) = avg(k) + DF(k) * (W(i,j,k) - avg(k))
          end do
        end do
      end do
      
      !$omp end parallel
    
    end if

  end subroutine SpongeTop
  
  
  
  subroutine SpongeTopScalar(Phi)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: Phi
    integer :: i, j, k, bufn
    real(knd) :: ze, zs, zb, p

    if (top_sponge_bottom<zW(Prnz)) then

      bufn = min(Prnz, int((zW(Prnz) - top_sponge_bottom) / dzmin) )
      zs = top_sponge_bottom
      ze = gzmax

      if (.not.allocated(DF)) allocate(DF(  min(Unz, Vnz, Wnz) - bufn  :  max(Unz, Vnz, Wnz)))
      if (.not.allocated(avg)) allocate(avg( min(Unz, Vnz, Wnz) - bufn  :  max(Unz, Vnz, Wnz)))

      !$omp parallel private(i,j,k,p,zb)
      
      !$omp do
      do k = Prnz-bufn, Prnz
        avg(k) = 0
      end do
      !$omp end do

      !$omp do
      do k = Prnz-bufn, Prnz
        p = 0
        do j = 1, Prny
          do i = 1, Prnx
            p = p + Phi(i,j,k)
          end do
        end do
        avg(k) = p
      end do
      !$omp end do

#ifdef PAR
        avg = par_co_sum(avg, comm_plane_xy)
#endif

      !$omp do
      do k = Prnz-bufn, Prnz
        avg(k) = avg(k) / (gPrnx*gPrny)
      end do
      !$omp end do

      !$omp do
      do k = Prnz-bufn, Prnz
        zb = (zPr(k)-zs) / (ze-zs)
        DF(k) = ScalarDampF(zb)
      end do
      !$omp end do

      !$omp do
      do k = Prnz-bufn, Prnz
        do j = -1, Prny+1
          do i = -1, Prnx+1
            Phi(i,j,k) = avg(k) + DF(k) * (Phi(i,j,k)-avg(k))
          end do
        end do
      end do
      !$omp end do

      !$omp end parallel
    end if

  endsubroutine SpongeTopScalar
  
  
  
  subroutine SpongeOut(U, V, W, temperature)
    real(knd), contiguous, intent(inout), dimension(-2:,-2:,-2:) :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: temperature
    
    !NOTE: currently having both of them enabled will likely lead to strange results
    if (enable_in_sponge_x) call SpongeIn_X(U, V, W, temperature)
    if (enable_out_sponge_x) call SpongeOut_X(U, V, W, temperature)
    if (enable_out_sponge_y) call SpongeOut_Y(U, V, W, temperature)
  end subroutine

  subroutine SpongeIn_X(U, V, W, temperature)
    real(knd), contiguous, intent(inout), dimension(-2:,-2:,-2:) :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: temperature
    integer   :: i, j, k, bufn
    integer   :: hi, lo, loU, hiU
    real(knd) :: p, xe, xs, xb, DF

    !size end extent of the buffer region
    bufn = min(max(5,Prnx/50), Prnx/4)
    xs = xU(0)
    xe = xU(bufn + 2)

    !extent of the probe region where the local average is taken from
    loU = bufn+1
    hiU = max(Unx/3, 50)

    lo = bufn+1
    hi = max(Prnx/3, 50)

    !$omp parallel private(i, j, k, p, xb, DF)

    !$omp do collapse(2)
    do k = 1, Unz
      do j = 1, Uny
        p = 0
        do i = loU, hiU
          p = p + U(i,j,k)
        end do
        p = p / (hiU - loU + 1)
        do i = Unx - bufn, Unx + 1
          xb = (xe-xU(i)) / (xe-xs)
          DF = DampF(xb)
          U(i,j,k) = p + DF * (U(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    !$omp do collapse(2)
    do k = 1, Vnz
      do j = 1, Vny
        p = 0
        do i = lo, hi
          p = p + V(i,j,k)
        end do
        p = p / (hi - lo + 1)
        do i = Vnx-bufn, Vnx + 1
          xb = (xe-xPr(i)) / (xe-xs)
          DF = DampF(xb)
          V(i,j,k) = p + DF * (V(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    !$omp do collapse(2)
    do k = 1, Wnz
      do j = 1, Wny
        p = 0
        do i = lo, hi
          p = p + W(i,j,k)
        end do
        p = p / (hi - lo + 1)
        do i = Wnx - bufn, Wnx + 1
          xb = (xe-xPr(i)) / (xe-xs)
          DF = DampF(xb)
          W(i,j,k) = p + DF * (W(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    if (enable_buoyancy) then
      !$omp do collapse(2)
      do k = 1, Prnz
        do j = 1, Prny
          p = 0
          do i = lo, hi
            p = p + temperature(i,j,k)
          end do
          p = p / (hi - lo + 1)
          do i = Prnx - bufn, Prnx + 1
            xb = (xe-xPr(i)) / (xe-xs)
            DF = DampF(xb)
            temperature(i,j,k) = p + DF * (temperature(i,j,k) - p)
          end do
        end do
      end do
      !$omp end do
    end if
    !$omp end parallel

  end subroutine SpongeIn_X


  subroutine SpongeOut_X(U, V, W, temperature)
    real(knd), contiguous, intent(inout), dimension(-2:,-2:,-2:) :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: temperature
    integer   :: i, j, k, bufn
    integer   :: hi, lo, loU, hiU
    real(knd) :: p, xe, xs, xb, DF

    !size end extent of the buffer region
    bufn = min(max(5,Prnx/50), Prnx/4)
    xs = xU(Prnx - bufn-2)
    xe = xU(Prnx)

    !extent of the probe region where the local average is taken from
    loU = max(2*Unx/3, Unx-50)
    hiU = Unx-bufn-1

    lo = max(2*Prnx/3, Prnx-50)
    hi = Prnx-bufn-1

    !$omp parallel private(i, j, k, p, xb, DF)

    !$omp do collapse(2)
    do k = 1, Unz
      do j = 1, Uny
        p = 0
        do i = loU, hiU
          p = p + U(i,j,k)
        end do
        p = p / (hiU - loU + 1)
        do i = Unx - bufn, Unx + 1
          xb = (xU(i)-xs) / (xe-xs)
          DF = DampF(xb)
          U(i,j,k) = p + DF * (U(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    !$omp do collapse(2)
    do k = 1, Vnz
      do j = 1, Vny
        p = 0
        do i = lo, hi
          p = p + V(i,j,k)
        end do
        p = p / (hi - lo + 1)
        do i = Vnx-bufn, Vnx + 1
          xb = (xPr(i)-xs) / (xe-xs)
          DF = DampF(xb)
          V(i,j,k) = p + DF * (V(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    !$omp do collapse(2)
    do k = 1, Wnz
      do j = 1, Wny
        p = 0
        do i = lo, hi
          p = p + W(i,j,k)
        end do
        p = p / (hi - lo + 1)
        do i = Wnx - bufn, Wnx + 1
          xb = (xPr(i)-xs) / (xe-xs)
          DF = DampF(xb)
          W(i,j,k) = p + DF * (W(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    if (enable_buoyancy) then
      !$omp do collapse(2)
      do k = 1, Prnz
        do j = 1, Prny
          p = 0
          do i = lo, hi
            p = p + temperature(i,j,k)
          end do
          p = p / (hi - lo + 1)
          do i = Prnx - bufn, Prnx + 1
            xb = (xPr(i)-xs) / (xe-xs)
            DF = DampF(xb)
            temperature(i,j,k) = p + DF * (temperature(i,j,k) - p)
          end do
        end do
      end do
      !$omp end do
    end if
    !$omp end parallel

  end subroutine SpongeOut_X


  subroutine SpongeOut_Y(U, V, W, temperature)
    use custom_par
    real(knd), contiguous, intent(inout), dimension(-2:,-2:,-2:) :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: temperature
    integer   :: i, j, k, bufn
    integer   :: hi, lo, loV, hiV
    real(knd) :: p, ye, ys, yb, DF

    !TODO: Properly parallelize for a sponge zone across several images
    if (jim==nyims) then
      !size end extent of the buffer region
      bufn = min(max(5,Prnx/50), gPrny/4, Prny/4)
      ys = yV(Prny - bufn-2)
      ye = yV(Prny)

      !extent of the probe region where the local average is taken from
      loV = max(2*gPrny/3-offset_to_global_y, gPrny-50-offset_to_global_y, 1)
      hiV = Vny-bufn-1

      lo = max(2*gPrny/3-offset_to_global_y, gPrny-50-offset_to_global_y, 1)
      hi = Prny-bufn-1

      !$omp parallel private(i, j, k, p, yb, DF)

      !$omp do collapse(2)
      do k = 1, Unz
        do i = 1, Unx
          p = 0
          do j = lo, hi
            p = p + U(i,j,k)
          end do
          p = p / (hi - lo + 1)
          do j = Uny - bufn, Uny + 1
            yb = (yPr(j)-ys) / (ye-ys)       
            DF = DampF(yb)
            U(i,j,k) = p + DF * (U(i,j,k) - p)
          end do
        end do
      end do
      !$omp end do

      !$omp do collapse(2)
      do k = 1, Vnz
        do i = 1, Vnx
          p = 0
          do j = loV, hiV
            p = p + V(i,j,k)
          end do
          p = p / (hiV - loV + 1)
          do j = Vny-bufn, Vny + 1
            yb = (yV(j)-ys) / (ye-ys)
            DF = DampF(yb)
            V(i,j,k) = p + DF * (V(i,j,k) - p)
          end do           
        end do
      end do
      !$omp end do

      !$omp do collapse(2)
      do k = 1, Wnz
        do i = 1, Wnx
          p = 0
          do j = lo, hi
            p = p + W(i,j,k)
          end do
          p = p / (hi - lo + 1)
          do j = Wny - bufn, Wny + 1
            yb = (yPr(j)-ys) / (ye-ys)
            DF = DampF(yb)
            W(i,j,k) = p + DF * (W(i,j,k) - p)
          end do
        end do
      end do
      !$omp end do

      if (enable_buoyancy) then
        !$omp do collapse(2)
        do k = 1, Prnz
          do i = 1, Prny
            p = 0
            do j = lo, hi
              p = p + temperature(i,j,k)
            end do
            p = p / (hi - lo + 1)
            do j = Prny - bufn, Prny + 1
              yb = (yPr(j)-ys) / (ye-ys)
              DF = DampF(yb)
              temperature(i,j,k) = p + DF * (temperature(i,j,k) - p)
            end do
          end do
        end do
        !$omp end do
      end if
      !$omp end parallel
    end if

  end subroutine SpongeOut_Y



end module
