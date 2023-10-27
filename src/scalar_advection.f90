module ScalarAdvection

  use Parameters
  use ArrayUtilities
  use Wallmodels
  
  implicit none
  
  private
  
  public AdvScalar, AddScalarAdvVector, ScalarAdvection_Deallocate
  
  real(knd), allocatable :: Slope(:,:,:)

  
contains



  subroutine AdvScalar(Scal2, Scal, U, V, W, dt, temperature_flux_profile)
    real(knd), contiguous, intent(out) :: Scal2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: Scal(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd),             intent(in)  :: dt
    real(knd), contiguous, intent(out), optional :: temperature_flux_profile(0:)
    real(knd), allocatable, save :: temperature_flux_profileLoc(:)

    if (.not.allocated(temperature_flux_profileLoc)) then
      allocate(temperature_flux_profileLoc(0:Prnz))
    end if

    
    if (discretization_order==4) then
      call KappaScalar_4ord(Scal2, Scal, &
                            U, V, W, &
                            temperature_flux_profileLoc)
    else
#ifdef MODIFIED_KAPPA
      call KappaScalar_mod_delta(Scal2, Scal, &
                                  U, V, W, &
                                  dt, &
                                  temperature_flux_profileLoc)
#else
      call KappaScalar(Scal2, Scal, &
                        U, V, W, &
                        temperature_flux_profileLoc)
#endif
    end if

    if (present(temperature_flux_profile)) then
      if (size(temperature_flux_profile)==size(temperature_flux_profileLoc)) &
        temperature_flux_profile(0:Prnz) = temperature_flux_profileLoc
    end if

  endsubroutine AdvScalar

  subroutine CDSScalar(Scal2, Scal, U, V, W)
    real(knd), contiguous, intent(out) :: Scal2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: Scal(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    integer :: nx, ny, nz, i, j, k
    real(knd) :: Ax, Ay, Az

    nx = Prnx
    ny = Prny
    nz = Prnz

    Scal2 = 0

    Ax = dxmin / 2
    Ay = dymin / 2
    Az = dzmin / 2


    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Scal2(i,j,k) = Scal2(i,j,k) - ( (Az * (Scal(i,j,k+1)+Scal(i,j,k)) * (W(i,j,k)) &
                                      -    Az * (Scal(i,j,k)+Scal(i,j,k-1)) * (W(i,j,k-1))) &
                                      +   (Ay * (Scal(i,j+1,k)+Scal(i,j,k)) * (V(i,j,k)) &
                                      -    Ay * (Scal(i,j,k)+Scal(i,j-1,k)) * (V(i,j-1,k))) &
                                      +   (Ax * (Scal(i+1,j,k)+Scal(i,j,k)) * (U(i,j,k)) &
                                      -    Ax * (Scal(i,j,k)+Scal(i-1,j,k)) * (U(i-1,j,k))))
        end do
      end do
    end do
  end subroutine CDSScalar



  subroutine KappaScalar(Scal2, Scal, &
                         U, V, W, &
                         temperature_flux_profile)
    !Kappa scheme with flux limiter
    !Hunsdorfer et al. 1995, JCP
    real(knd), contiguous, intent(out) :: Scal2(-2:,-2:,-2:) 
    real(knd), contiguous, intent(in)  :: Scal(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd), contiguous, intent(out) :: temperature_flux_profile(0:)
    integer   :: i, j, k, l
    real(knd) :: Ax, Ay, Az              !Auxiliary variables to store muliplication constants for efficiency
    real(knd) :: vel, sl, sr, flux
    real(knd), parameter :: eps = 1e-8

    if (.not.allocated(Slope)) then
      allocate(Slope(-1:Prnx+2, -1:Prny+2, -1:Prnz+2))
    end if

    Ax = 1 / dxmin
    Ay = 1 / dymin
    Az = 1 / dzmin


    call set(Scal2, 0._knd)
    call set(Slope, 0._knd)

    !$omp parallel private(i,j,k,l,vel,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 1, Prny
      do i = 0, Prnx
       if (U(i,j,k)>0) then
        sr = Scal(i+1,j,k) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i-1,j,k)
       else
        sr = Scal(i,j,k) - Scal(i+1,j,k)
        sl = Scal(i+1,j,k) - Scal(i+2,j,k)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)))
      end do
     end do
    end do
    !$omp end do

    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 1, Prny
      do i = 0, Prnx
        if (Scflx_mask(i,j,k)) then
          if (U(i,j,k)>0) then
           flux = U(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i-1,j,k)) * Slope(i,j,k)/2._knd)
          else
           flux = U(i,j,k) * (Scal(i+1,j,k) + (Scal(i+1,j,k)-Scal(i+2,j,k)) * Slope(i,j,k)/2._knd)
          end if

          Scal2(i,j,k) = Scal2(i,j,k) - Ax * flux
          Scal2(i+1,j,k) = Scal2(i+1,j,k) + Ax * flux
        end if
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel


    call set(Slope, 0._knd)

    !$omp parallel private(i,j,k,l,vel,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 0, Prny
      do i = 1, Prnx
       if (V(i,j,k)>0) then
        sr = Scal(i,j+1,k) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i,j-1,k)
       else
        sr = Scal(i,j,k) - Scal(i,j+1,k)
        sl = Scal(i,j+1,k) - Scal(i,j+2,k)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)))
      end do
     end do
    end do
    !$omp end do


    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 0, Prny
      do i = 1, Prnx
        if (Scfly_mask(i,j,k)) then
          if (V(i,j,k)>0) then
           flux = V(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j-1,k)) * Slope(i,j,k)/2._knd)
          else
           flux = V(i,j,k) * (Scal(i,j+1,k) + (Scal(i,j+1,k)-Scal(i,j+2,k)) * Slope(i,j,k)/2._knd)
          end if

          Scal2(i,j,k) = Scal2(i,j,k) - Ay * flux
          Scal2(i,j+1,k) = Scal2(i,j+1,k) + Ay * flux
        end if
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel


    call set(Slope, 0._knd)

    !$omp parallel private(i,j,k,l,vel,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 0, Prnz
     do j = 1, Prny
      do i = 1, Prnx
       if (W(i,j,k)>0) then
        sr = Scal(i,j,k+1) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i,j,k-1)
       else
        sr = Scal(i,j,k) - Scal(i,j,k+1)
        sl = Scal(i,j,k+1) - Scal(i,j,k+2)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)))
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel

    call set(temperature_flux_profile, 0._knd)

    !$omp parallel private(i,j,k,l,vel,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    do l = 0, 1  !odd-even separation to avoid a race condition
      !$omp do reduction(+:temperature_flux_profile) schedule(runtime)
      do k = 0+l, Prnz, 2
       do j = 1, Prny
        do i = 1, Prnx
          if (Scflz_mask(i,j,k)) then
            vel = W(i,j,k)
            if (vel>0) then
             flux = vel * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j,k-1)) * Slope(i,j,k)/2._knd)
            else
             flux = vel * (Scal(i,j,k+1) + (Scal(i,j,k+1)-Scal(i,j,k+2)) * Slope(i,j,k)/2._knd)
            end if

            temperature_flux_profile(k) = temperature_flux_profile(k) + flux

            Scal2(i,j,k) = Scal2(i,j,k) - Az * flux
            Scal2(i,j,k+1) = Scal2(i,j,k+1) + Az * flux
          end if
        end do
       end do
      end do
      !$omp end do
    end do
    !$omp end parallel

  contains

    real(knd) pure function  FluxLimiter(r)
      real(knd), intent(in) :: r

      FluxLimiter = max(0._knd, &
                        min(2 * r, &
                            min(2._knd, &
                                (1 + 2*r) / 3) ) )
    end function

  endsubroutine KappaScalar




  subroutine KappaScalar_4ord(Scal2, Scal, &
                         U, V, W, &
                         temperature_flux_profile)
    !version with advection velocities following the 4th-order discrete divergence-free condition
    ! c.f. Hokpunna, Manhart, 2010, eq. 11, https://dx.doi.org/10.1016/j.jcp.2010.05.042
    ! [U]i = -1/24 U(i+1) + 13/12 U(i) - 1/24 U(i-1)
    ! ([U]i - [U]i-1) / dx = 9/8*(U(i)-U(i-1)) / dx - 1/8*(U(i+1)-U(i-2))/(3*dx)
                         
    !Kappa scheme with flux limiter
    !Hunsdorfer et al. 1995, JCP
    real(knd), contiguous, intent(out) :: Scal2(-2:,-2:,-2:) 
    real(knd), contiguous, intent(in)  :: Scal(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd), contiguous, intent(out) :: temperature_flux_profile(0:)
    integer   :: i, j, k, l
    real(knd) :: Ax, Ay, Az              !Auxiliary variables to store muliplication constants for efficiency
    real(knd) :: sl, sr, flux
    real(knd) :: Uadv, Vadv, Wadv
    real(knd), parameter :: D0 = 13._knd / 12, D1 = -1._knd / 24
    real(knd), parameter :: eps = 1e-8


    if (.not.allocated(Slope)) then
      allocate(Slope(-1:Prnx+2, -1:Prny+2, -1:Prnz+2))
    end if

    Ax = 1 / dxmin
    Ay = 1 / dymin
    Az = 1 / dzmin


    call set(Scal2, 0._knd)
    call set(Slope, 0._knd)

    !$omp parallel private(i,j,k,l,Uadv,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 1, Prny
      do i = 0, Prnx
       Uadv = D1 * U(i-1,j,k) + D0 * U(i,j,k) + D1 * U(i+1,j,k)
       if (Uadv > 0) then
        sr = Scal(i+1,j,k) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i-1,j,k)
       else
        sr = Scal(i,j,k) - Scal(i+1,j,k)
        sl = Scal(i+1,j,k) - Scal(i+2,j,k)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)))
      end do
     end do
    end do
    !$omp end do

    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 1, Prny
      do i = 0, Prnx
        if (Scflx_mask(i,j,k)) then
          Uadv = D1 * U(i-1,j,k) + D0 * U(i,j,k) + D1 * U(i+1,j,k)
          if (Uadv > 0) then
           flux = Uadv * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i-1,j,k)) * Slope(i,j,k)/2._knd)
          else
           flux = Uadv * (Scal(i+1,j,k) + (Scal(i+1,j,k)-Scal(i+2,j,k)) * Slope(i,j,k)/2._knd)
          end if

          Scal2(i,j,k) = Scal2(i,j,k) - Ax * flux
          Scal2(i+1,j,k) = Scal2(i+1,j,k) + Ax * flux
        end if
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel


    call set(Slope, 0._knd)

    !$omp parallel private(i,j,k,l,Vadv,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 0, Prny
      do i = 1, Prnx
       Vadv = D1 * V(i,j-1,k) + D0 * V(i,j,k) + D1 * V(i,j+1,k)
       if (Vadv > 0) then
        sr = Scal(i,j+1,k) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i,j-1,k)
       else
        sr = Scal(i,j,k) - Scal(i,j+1,k)
        sl = Scal(i,j+1,k) - Scal(i,j+2,k)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)))
      end do
     end do
    end do
    !$omp end do


    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 0, Prny
      do i = 1, Prnx
        if (Scfly_mask(i,j,k)) then
          Vadv = D1 * V(i,j-1,k) + D0 * V(i,j,k) + D1 * V(i,j+1,k)
          if (Vadv > 0) then
           flux = Vadv * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j-1,k)) * Slope(i,j,k)/2._knd)
          else
           flux = Vadv * (Scal(i,j+1,k) + (Scal(i,j+1,k)-Scal(i,j+2,k)) * Slope(i,j,k)/2._knd)
          end if

          Scal2(i,j,k) = Scal2(i,j,k) - Ay * flux
          Scal2(i,j+1,k) = Scal2(i,j+1,k) + Ay * flux
        end if
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel


    call set(Slope, 0._knd)

    !$omp parallel private(i,j,k,l,Wadv,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 0, Prnz
     do j = 1, Prny
      do i = 1, Prnx
       Wadv = D1 * W(i,j,k-1) + D0 * W(i,j,k) + D1 * W(i,j,k+1)
       if (Wadv > 0) then
        sr = Scal(i,j,k+1) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i,j,k-1)
       else
        sr = Scal(i,j,k) - Scal(i,j,k+1)
        sl = Scal(i,j,k+1) - Scal(i,j,k+2)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)))
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel

    call set(temperature_flux_profile, 0._knd)

    !$omp parallel private(i,j,k,l,Wadv,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    do l = 0, 1  !odd-even separation to avoid a race condition
      !$omp do reduction(+:temperature_flux_profile) schedule(runtime)
      do k = 0+l, Prnz, 2
       do j = 1, Prny
        do i = 1, Prnx
          if (Scflz_mask(i,j,k)) then
            Wadv = D1 * W(i,j,k-1) + D0 * W(i,j,k) + D1 * W(i,j,k+1)
            if (Wadv > 0) then
             flux = Wadv * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j,k-1)) * Slope(i,j,k)/2._knd)
            else
             flux = Wadv * (Scal(i,j,k+1) + (Scal(i,j,k+1)-Scal(i,j,k+2)) * Slope(i,j,k)/2._knd)
            end if

            temperature_flux_profile(k) = temperature_flux_profile(k) + flux

            Scal2(i,j,k) = Scal2(i,j,k) - Az * flux
            Scal2(i,j,k+1) = Scal2(i,j,k+1) + Az * flux
          end if
        end do
       end do
      end do
      !$omp end do
    end do
    !$omp end parallel

  contains

    real(knd) pure function  FluxLimiter(r)
      real(knd), intent(in) :: r

      FluxLimiter = max(0._knd, &
                        min(2 * r, &
                            min(2._knd, &
                                (1 + 2*r) / 3) ) )
    end function

  endsubroutine KappaScalar_4ord





  subroutine KappaScalar_mod_delta(Scal2, Scal, &
                                   U, V, W, &
                                   dt, &
                                   temperature_flux_profile) 
    !Kappa scheme with flux limiter
    !Hunsdorfer et al. 1995, JCP
    !delta modified for smaller Courant numbers according to Hunsdorfer et al.
    real(knd), contiguous, intent(out) :: Scal2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: Scal(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd),             intent(in)  :: dt
    real(knd), contiguous, intent(out) :: temperature_flux_profile(0:)
    integer   :: i, j, k, l
    real(knd) :: Ax, Ay, Az              !Auxiliary variables to store muliplication constants for efficiency
    real(knd) :: vel, sl, sr, flux
    real(knd), parameter :: eps = 1e-8

    if (.not.allocated(Slope)) then
      allocate(Slope(-1:Prnx+2, -1:Prny+2, -1:Prnz+2))
    end if

    Ax = 1 / dxmin
    Ay = 1 / dymin
    Az = 1 / dzmin


    call set(Scal2, 0._knd)
    call set(Slope, 0._knd)

    !$omp parallel private(i,j,k,l,vel,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 1, Prny
      do i = 0, Prnx
       if (U(i,j,k)>0) then
        sr = Scal(i+1,j,k) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i-1,j,k)
       else
        sr = Scal(i,j,k) - Scal(i+1,j,k)
        sl = Scal(i+1,j,k) - Scal(i+2,j,k)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)), &
                                  abs(U(i,j,k)) * dt / dxmin)
      end do
     end do
    end do
    !$omp end do

    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 1, Prny
      do i = 0, Prnx
        if (Scflx_mask(i,j,k)) then
          if (U(i,j,k)>0) then
           flux = U(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i-1,j,k)) * Slope(i,j,k)/2._knd)
          else
           flux = U(i,j,k) * (Scal(i+1,j,k) + (Scal(i+1,j,k)-Scal(i+2,j,k)) * Slope(i,j,k)/2._knd)
          end if

          Scal2(i,j,k) = Scal2(i,j,k) - Ax * flux
          Scal2(i+1,j,k) = Scal2(i+1,j,k) + Ax * flux
        end if
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel


    call set(Slope, 0._knd)

    !$omp parallel private(i,j,k,l,vel,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 0, Prny
      do i = 1, Prnx
       if (V(i,j,k)>0) then
        sr = Scal(i,j+1,k) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i,j-1,k)
       else
        sr = Scal(i,j,k) - Scal(i,j+1,k)
        sl = Scal(i,j+1,k) - Scal(i,j+2,k)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)), &
                                  abs(V(i,j,k)) * dt / dymin)
      end do
     end do
    end do
    !$omp end do


    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 0, Prny
      do i = 1, Prnx
        if (Scfly_mask(i,j,k)) then
          if (V(i,j,k)>0) then
           flux = V(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j-1,k)) * Slope(i,j,k)/2._knd)
          else
           flux = V(i,j,k) * (Scal(i,j+1,k) + (Scal(i,j+1,k)-Scal(i,j+2,k)) * Slope(i,j,k)/2._knd)
          end if

          Scal2(i,j,k) = Scal2(i,j,k) - Ay * flux
          Scal2(i,j+1,k) = Scal2(i,j+1,k) + Ay * flux
        end if
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel


    call set(Slope, 0._knd)

    !$omp parallel private(i,j,k,l,vel,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 0, Prnz
     do j = 1, Prny
      do i = 1, Prnx
       if (W(i,j,k)>0) then
        sr = Scal(i,j,k+1) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i,j,k-1)
       else
        sr = Scal(i,j,k) - Scal(i,j,k+1)
        sl = Scal(i,j,k+1) - Scal(i,j,k+2)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)), &
                                  abs(W(i,j,k)) * dt / dzmin)
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel

    call set(temperature_flux_profile, 0._knd)

    !$omp parallel private(i,j,k,l,vel,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    do l = 0, 1  !odd-even separation to avoid a race condition
      !$omp do reduction(+:temperature_flux_profile) schedule(runtime)
      do k = 0+l, Prnz, 2
       do j = 1, Prny
        do i = 1, Prnx
          if (Scflz_mask(i,j,k)) then
            vel = W(i,j,k)
            if (vel>0) then
             flux = vel * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j,k-1)) * Slope(i,j,k)/2._knd)
            else
             flux = vel * (Scal(i,j,k+1) + (Scal(i,j,k+1)-Scal(i,j,k+2)) * Slope(i,j,k)/2._knd)
            end if

            temperature_flux_profile(k) = temperature_flux_profile(k) + flux

            Scal2(i,j,k) = Scal2(i,j,k) - Az * flux
            Scal2(i,j,k+1) = Scal2(i,j,k+1) + Az * flux
          end if
        end do
       end do
      end do
      !$omp end do
    end do
    !$omp end parallel

  contains

    real(knd) pure function  FluxLimiter(r, C)
      real(knd), intent(in) :: r !slope ratio
      real(knd), intent(in) :: C !Courant number
      real(knd), parameter :: eps = 1e-6_knd

      FluxLimiter = max(0._knd, &
                        min(2 * r, &
                            min(2 * max(1._knd, (1 - C) / (C + eps)), &
                                (1 + 2 * r) / 3) ) )
    end function

  endsubroutine KappaScalar_mod_delta
  
  
  
  
  
  
  
  
  
  
  
  
  subroutine CDS4Scalar(Scal2, Scal, &
                         U, V, W, &
                         temperature_flux_profile)
    !advection velocities following the 4th-order discrete divergence-free condition
    ! c.f. Hokpunna, Manhart, 2010, eq. 11, https://dx.doi.org/10.1016/j.jcp.2010.05.042
    ! [U]i = -1/24 U(i+1) + 13/12 U(i) - 1/24 U(i-1)
    ! ([U]i - [U]i-1) / dx = 9/8*(U(i)-U(i-1)) / dx - 1/8*(U(i+1)-U(i-2))/(3*dx)
                         
    ! Central scheme, e.g. in Wicker, Skamarock (2002), eq. 4b or Wesseling (2001) p. 150
    ! F = Uadv * (7/12 (C(j) + C(j+1)) - 1/12 (C(j-1) + C(j+2)))
    
    real(knd), contiguous, intent(out) :: Scal2(-2:,-2:,-2:) 
    real(knd), contiguous, intent(in)  :: Scal(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd), contiguous, intent(out) :: temperature_flux_profile(0:)
    integer   :: i, j, k, l
    real(knd) :: Ax, Ay, Az              !Auxiliary variables to store muliplication constants for efficiency
    real(knd) :: sl, sr, flux
    real(knd) :: Uadv, Vadv, Wadv
    real(knd), parameter :: D0 = 13._knd / 12, D1 = -1._knd / 24
    real(knd), parameter :: C0 = 7._knd / 12, C1 = 1._knd / 12
    real(knd), parameter :: eps = 1e-8

    Ax = 1 / dxmin
    Ay = 1 / dymin
    Az = 1 / dzmin


    call set(Scal2, 0._knd)

    !$omp parallel private(i,j,k,l,Uadv,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 1, Prny
      do i = 0, Prnx
        if (Scflx_mask(i,j,k)) then
          Uadv = D1 * U(i-1,j,k) + D0 * U(i,j,k) + D1 * U(i+1,j,k)

          flux = Uadv * (C0 * (Scal(i,j,k) + Scal(i+1,j,k)) - C1 * (Scal(i-1,j,k) + Scal(i+2,j,k)))

          Scal2(i,j,k) = Scal2(i,j,k) - Ax * flux
          Scal2(i+1,j,k) = Scal2(i+1,j,k) + Ax * flux
        end if
      end do
     end do
    end do
    !$omp end do

    !$omp do schedule(runtime)
    do k = 1, Prnz
     do j = 0, Prny
      do i = 1, Prnx
        if (Scfly_mask(i,j,k)) then
          Vadv = D1 * V(i,j-1,k) + D0 * V(i,j,k) + D1 * V(i,j+1,k)

          flux = Vadv * (C0 * (Scal(i,j,k) + Scal(i,j+1,k)) - C1 * (Scal(i,j-1,k) + Scal(i,j+2,k)))

          Scal2(i,j,k) = Scal2(i,j,k) - Ay * flux
          Scal2(i,j+1,k) = Scal2(i,j+1,k) + Ay * flux
        end if
      end do
     end do
    end do
    !$omp end do nowait
    !$omp end parallel


    call set(temperature_flux_profile, 0._knd)

    !$omp parallel private(i,j,k,l,Wadv,sl,sr,flux) shared(Slope,Scal,Scal2,temperature_flux_profile)
    do l = 0, 1  !odd-even separation to avoid a race condition
      !$omp do reduction(+:temperature_flux_profile) schedule(runtime)
      do k = 0+l, Prnz, 2
       do j = 1, Prny
        do i = 1, Prnx
          if (Scflz_mask(i,j,k)) then
            Wadv = D1 * W(i,j,k-1) + D0 * W(i,j,k) + D1 * W(i,j,k+1)

            flux = Wadv * (C0 * (Scal(i,j,k) + Scal(i,j,k+1)) - C1 * (Scal(i,j,k-1) + Scal(i,j,k+2)))

            temperature_flux_profile(k) = temperature_flux_profile(k) + flux

            Scal2(i,j,k) = Scal2(i,j,k) - Az * flux
            Scal2(i,j,k+1) = Scal2(i,j,k+1) + Az * flux
          end if
        end do
       end do
      end do
      !$omp end do
    end do
    !$omp end parallel

  endsubroutine CDS4Scalar

  
  






  subroutine AddScalarAdvVector(ScU, ScV, ScW, Scal, U, V, W, weight, probes_flux, px, py, pz) !Kappa scheme with flux limiter
    real(knd), contiguous, intent(inout) :: ScU(:,:,:) !Hunsdorfer et al. 1995, JCP
    real(knd), contiguous, intent(inout) :: ScV(:,:,:)
    real(knd), contiguous, intent(inout) :: ScW(:,:,:)
    real(knd), contiguous, intent(in)    :: Scal(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)    :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd), intent(in) :: weight
    real(knd), intent(out) :: probes_flux(:,:) !component, position
    integer, intent(in) :: px(:), py(:), pz(:)
    integer :: i, j, k, probe
    real(knd) :: sl, sr
    real(knd), parameter ::eps = 1e-8

    if (.not.allocated(Slope)) then
      allocate(Slope(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    end if


    call set(Slope,0._knd)

    !$omp parallel private(i,j,k,sl,sr) shared(Slope,Scal,ScU)
    !$omp do
    do k = 1, Unz
     do j = 1, Uny
      do i = 1, Unx
       if (U(i,j,k)>0) then
        sr = Scal(i+1,j,k) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i-1,j,k)
       else
        sr = Scal(i,j,k) - Scal(i+1,j,k)
        sl = Scal(i+1,j,k) - Scal(i+2,j,k)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)))
      end do
     end do
    end do
    !$omp end do

    !$omp do
    do k = 1, Unz
     do j = 1, Uny
      do i = 1, Unx
       if (U(i,j,k)>0) then
        ScU(i,j,k) = ScU(i,j,k) + &
                   weight * U(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i-1,j,k)) * Slope(i,j,k)/2._knd)
       else
        ScU(i,j,k) = ScU(i,j,k) + &
                   weight * U(i,j,k) * (Scal(i+1,j,k) + (Scal(i+1,j,k)-Scal(i+2,j,k)) * Slope(i,j,k)/2._knd)
       end if
      end do
     end do
    end do
    !$omp end do nowait

    !$omp do
    do probe = 1, size(px)
        i = px(probe)
        j = py(probe)
        k = pz(probe)
        if (U(i,j,k)>0) then
          probes_flux(1, probe) = U(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i-1,j,k)) * Slope(i,j,k)/2._knd)
        else
          probes_flux(1, probe) = U(i,j,k) * (Scal(i+1,j,k) + (Scal(i+1,j,k)-Scal(i+2,j,k)) * Slope(i,j,k)/2._knd)
       end if
    end do
    !$omp end do nowait
    !$omp end parallel

    call set(Slope,0._knd)

    !$omp parallel private(i,j,k,sl,sr) shared(Slope,Scal,ScV)
    !$omp do
    do k = 1, Vnz
     do j = 1, Vny
      do i = 1, Vnx
       if (V(i,j,k)>0) then
        sr = Scal(i,j+1,k) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i,j-1,k)
       else
        sr = Scal(i,j,k) - Scal(i,j+1,k)
        sl = Scal(i,j+1,k) - Scal(i,j+2,k)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)))
      end do
     end do
    end do
    !$omp end do


    !$omp do
    do k = 1, Vnz
     do j = 1, Vny
      do i = 1, Vnx
       if (V(i,j,k)>0) then
        ScV(i,j,k) = ScV(i,j,k) + &
                   weight * V(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j-1,k)) * Slope(i,j,k)/2._knd)
       else
        ScV(i,j,k) = ScV(i,j,k) + &
                   weight * V(i,j,k) * (Scal(i,j+1,k) + (Scal(i,j+1,k)-Scal(i,j+2,k)) * Slope(i,j,k)/2._knd)
       end if
      end do
     end do
    end do
    !$omp end do nowait

    !$omp do
    do probe = 1, size(px)
        i = px(probe)
        j = py(probe)
        k = pz(probe)
        if (V(i,j,k)>0) then
          probes_flux(2, probe) = V(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j-1,k)) * Slope(i,j,k)/2._knd)
        else
          probes_flux(2, probe) = V(i,j,k) * (Scal(i,j+1,k) + (Scal(i,j+1,k)-Scal(i,j+2,k)) * Slope(i,j,k)/2._knd)
       end if
    end do
    !$omp end do nowait
    !$omp end parallel


    call set(Slope ,0._knd)

    !$omp parallel private(i,j,k,sl,sr) shared(Slope,Scal,ScW)
    !$omp do
    do k = 1, Wnz
     do j = 1, Wny
      do i = 1, Wnx
       if (W(i,j,k)>0) then
        sr = Scal(i,j,k+1) - Scal(i,j,k)
        sl = Scal(i,j,k) - Scal(i,j,k-1)
       else
        sr = Scal(i,j,k) - Scal(i,j,k+1)
        sl = Scal(i,j,k+1) - Scal(i,j,k+2)
       end if
       Slope(i,j,k) = FluxLimiter((sr + eps*sign(1._knd, sl)) / (sl + eps*sign(1._knd, sl)))
      end do
     end do
    end do
    !$omp end do

    !$omp do
    do k = 1, Wnz
     do j = 1, Wny
      do i = 1, Wnx
       if (W(i,j,k)>0) then
        ScW(i,j,k) = ScW(i,j,k) + &
                   weight * W(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j,k-1)) * Slope(i,j,k)/2._knd)
       else
        ScW(i,j,k) = ScW(i,j,k) + &
                   weight * W(i,j,k) * (Scal(i,j,k+1) + (Scal(i,j,k+1)-Scal(i,j,k+2)) * Slope(i,j,k)/2._knd)
       end if
      end do
     end do
    end do
    !$omp end do

    !$omp do
    do probe = 1, size(px)
        i = px(probe)
        j = py(probe)
        k = pz(probe)
        if (W(i,j,k)>0) then
          probes_flux(3, probe) = W(i,j,k) * (Scal(i,j,k) + (Scal(i,j,k)-Scal(i,j,k-1)) * Slope(i,j,k)/2._knd)
        else
          probes_flux(3, probe) = W(i,j,k) * (Scal(i,j,k+1) + (Scal(i,j,k+1)-Scal(i,j,k+2)) * Slope(i,j,k)/2._knd)
       end if
    end do
    !$omp end do nowait
    !$omp end parallel

  contains

    real(knd) pure function  FluxLimiter(r)
      real(knd), intent(in) :: r
      FluxLimiter = max(0._knd, &
                        min(2._knd*r, &
                            min(2._knd, &
                                (1 + 2 * r) / 3)))
    end function

  endsubroutine AddScalarAdvVector
  





  subroutine ScalarAdvection_Deallocate
    if (allocated(Slope)) deallocate(Slope)
  end subroutine ScalarAdvection_Deallocate



end module ScalarAdvection
