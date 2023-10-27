module Scalars
 use Parameters
 use ArrayUtilities
 use Wallmodels
 use Boundaries
 use ScalarBoundaries
 use ScalarAdvection
 use ScalarDiffusion
 use Sponge, only: enable_top_sponge_scalar, SpongeTopScalar

implicit none
  private
  public ScalarRK3, constPr_sgs, Rig, partdiam, partrho, percdistrib, &
     ComputeTDiff, &
     TScalarProfile, &
     TemperatureProfileObj, MoistureProfileObj, SubsidenceProfileObj, &
     Interpolate1DProfile, InitScalarProfile, &
     InitScalar, &
     SubsidenceGradient, InitSubsidenceProfile, &
     AddScalarAdvVector, AddScalarDiffVector, &
     Scalars_Deallocate


  real(knd), dimension(:), allocatable :: partdiam, partrho, percdistrib !diameter of particles <=0 for gas

  real(knd) :: SubsidenceGradient = 0

  type TScalarProfileSection
    real(knd) :: top, height
    real(knd) :: jump
    real(knd) :: gradient
  end type


  type TScalarProfile
    real(knd), allocatable :: points(:,:)
    type(TScalarProfileSection), allocatable :: sections(:)
    logical   :: randomize
    real(knd) :: randomizeTop
    real(knd) :: randomizeAmplitude
  end type

  type(TScalarProfile) :: TemperatureProfileObj
  type(TScalarProfile) :: MoistureProfileObj

  type(TScalarProfile) :: SubsidenceProfileObj

  !module variables to enable deallocation before program end
  real(knd), dimension(:,:,:),  allocatable :: Temperature_adv, Temperature_2
  real(knd), dimension(:,:,:),  allocatable :: Moisture_adv,    Moisture_2
  real(knd), dimension(:,:,:,:), allocatable, save :: Scalar_adv,      Scalar_2

  abstract interface
    subroutine extra_interface
    end subroutine
  end interface

contains


  subroutine Scalars_Deallocate
    if (allocated(Temperature_adv)) deallocate(Temperature_adv)
    if (allocated(Moisture_adv)) deallocate(Moisture_adv)
    if (allocated(Scalar_adv)) deallocate(Scalar_adv)
    if (allocated(Temperature_2)) deallocate(Temperature_2)
    if (allocated(Moisture_2)) deallocate(Moisture_2)
    if (allocated(Scalar_2)) deallocate(Scalar_2)
    call ScalarAdvection_Deallocate
    call ScalarDiffusion_Deallocate
  end subroutine

  subroutine ScalarRK3(U, V, W, Pr, &
                       Temperature, Moisture, Scalar, &
                       RK_stage, dt, &
                       temperature_flux_profile, moisture_flux_profile)
    use RK3
    use VolumeSources, only: ScalarVolumeSources
    use Puffs, only: DoPuffs, PreparePuffs
    use WaterThermodynamics, only: compute_liquid_water_content
    real(knd), contiguous, intent(in)    :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)    :: Pr(-1:,-1:,-1:)
    real(knd), contiguous, intent(inout) :: Temperature(-2:,-2:,-2:), Moisture(-2:,-2:,-2:)
    real(knd), contiguous, intent(inout) :: Scalar(-2:,-2:,-2:,1:)
    real(knd), intent(in)                :: dt
    real(knd), contiguous, intent(out)   :: temperature_flux_profile(:)
    real(knd), contiguous, intent(out)   :: moisture_flux_profile(:)
    integer,  intent(in)                 :: RK_stage
    integer :: sc
    integer, save :: called = 0


    if (called==0) then

      called = 1

      allocate(Scalar_adv, mold=Scalar)
      allocate(Scalar_2, mold=Scalar)

      allocate(Temperature_adv, mold=Temperature)
      allocate(Temperature_2, mold=Temperature)

      allocate(Moisture_adv, mold=Moisture)
      allocate(Moisture_2, mold=Moisture)
    end if

    call DivergenceWM(U, V, W)

    if (enable_buoyancy) then

      call stage(Temperature, Temperature_2, Temperature_adv, &
                 ScalarTypeTemperature, TempBtype, &
                 BoundTemperature, TemperatureExtra, &
                 temperature_flux_profile)
      block
        !TODO:move into a procedure reduce the number of involved points. 
        !FIXME Ultimately, fix the spurious fluxes so that this is not needed.
        integer :: i, j, k, ii, jj, kk, n
        real(knd) :: t
        do k = -1, Prnz+2
          do j = -1, Prny+2
            do i = -1, Prnx+2
              if (Prtype(i,j,k)>0) then
                t = 0
                n = 0
                if (i>-1 .and. Prtype(i-1,j,k)<=0) then
                  t = t + Temperature(i-1,j,k) 
                  n = n + 1
                end if
                if (i<Prnx+2 .and. Prtype(i+1,j,k)<=0) then
                  t = t + Temperature(i+1,j,k) 
                  n = n + 1
                end if
                if (j>-1 .and. Prtype(i,j-1,k)<=0) then
                  t = t + Temperature(i,j-1,k) 
                  n = n + 1
                end if
                if (j<Prny+2 .and. Prtype(i,j+1,k)<=0) then
                  t = t + Temperature(i,j+1,k) 
                  n = n + 1
                end if
                if (k>-1 .and. Prtype(i,j,k-1)<=0) then
                  t = t + Temperature(i,j,k-1) 
                  n = n + 1
                end if
                if (k<Prnz+2 .and. Prtype(i,j,k+1)<=0) then
                  t = t + Temperature(i,j,k+1) 
                  n = n + 1
                end if
                if (n>0) then
                  Temperature(i,j,k) = t / n
                else
                  t = 0
                  n = 0
                  if (i>-1) then
                    t = t + Temperature(i-1,j,k) 
                    n = n + 1
                  end if
                  if (i<Prnx+2) then
                    t = t + Temperature(i+1,j,k) 
                    n = n + 1
                  end if
                  if (j>-1) then
                    t = t + Temperature(i,j-1,k) 
                    n = n + 1
                  end if
                  if (j<Prny+2) then
                    t = t + Temperature(i,j+1,k) 
                    n = n + 1
                  end if
                  if (k>-1) then
                    t = t + Temperature(i,j,k-1) 
                    n = n + 1
                  end if
                  if (k<Prnz+2) then
                    t = t + Temperature(i,j,k+1) 
                    n = n + 1
                  end if
                end if
              end if
            end do
          end do
        end do
      end block
    end if

    if (enable_moisture) then

      call stage(Moisture, Moisture_2, Moisture_adv, &
                 ScalarTypeMoisture, MoistBtype, &
                 BoundMoisture, MoistureExtra, &
                 moisture_flux_profile)
      
      where (Prtype(-2:Prnx+3,-2:Prny+3,-2:Prnz+3)>0) Moisture = moisture_ref
      
    end if
    
    if (enable_liquid) then
      call compute_liquid_water_content(Temperature, Moisture, Pr)
    end if

    call PreparePuffs(Scalar, RK_stage, RK_stages, time_stepping%time, time_stepping%dt)

    do sc = 1, num_of_scalars
      call stage(Scalar(:,:,:,sc), Scalar_2(:,:,:,sc), Scalar_adv(:,:,:,sc), &
                 ScalarTypePassive, ScalBtype, &
                 BoundScalar, ScalarExtra)
    end do

    if (computedeposition>0) call Deposition(Scalar_2,2._knd*RK_alpha(RK_stage)*dt)

    if (computegravsettling>0) call GravSettling(Scalar_2,2._knd*RK_alpha(RK_stage)*dt)


    contains

      subroutine stage(Array, Array2, Array_adv, &
                       scalar_type, btype, &
                       boundary_procedure, extra_procedure, &
                       flux_profile)
        use Large_scale_processes, only: enable_subsidence_profile, apply_subsidence_profile
        real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: Array, Array2, Array_adv
        integer, intent(in) :: scalar_type, btype(6)
        procedure(boundary_interface) :: boundary_procedure
        procedure(extra_interface) :: extra_procedure
        real(knd), intent(out), contiguous, optional :: flux_profile(0:)
        integer :: i,j,k


        call boundary_procedure(Array)

        if (RK_stage>1) then

          call assign(Array2, Array_adv)
          call multiply(Array2, RK_rho(RK_stage)*dt)

        else

          call set(Array2, 0._knd)

        end if

        call AdvScalar(Array_adv, Array, U, V, W, dt, flux_profile)
        
        !$omp parallel do private(i,j,k)
        do k = 0, Prnz+1
          do j = 0, Prny+1
            do i = 0, Prnx+1
              if (Prtype(i,j,k)>0) Array_adv(i,j,k) = 0
            end do
          end do
        end do
        !$omp end parallel do

        if (explicit_scalar_diffusion) then
          if (discretization_order==4) then
            call ScalarDiffusion_4ord_5point(Array_adv, Array)
          else
            call ScalarDiffusion_nobranch(Array_adv, Array)
          end if
        end if

        if (enable_subsidence_profile) call apply_subsidence_profile(Array_adv, Array)

        call extra_procedure

        call add_multiplied(Array2, Array_adv, RK_beta(RK_stage)*dt)

        call add(Array, Array2)

        if (.not.explicit_scalar_diffusion) then
          call boundary_procedure(Array)

          call ScalarDiffusion_implicit(Array2, Array, &
                          scalar_type, boundary_procedure, 2._knd*RK_alpha(RK_stage)*dt)
          call assign(Array, Array2)
        end if

        if (enable_top_sponge_scalar) call SpongeTopScalar(Array)

        call boundary_procedure(Array)

      end subroutine
      
      subroutine TemperatureExtra
        use DiabaticProcesses
        use Large_scale_processes
        use VolumeSources, only: TemperatureVolumeSources
        integer :: i, first, last

        call TemperatureVolumeSources(Temperature_adv)
        
        if (enable_diabatic_processes) call diabatic_heat_terms(Temperature_adv)

        if (enable_large_scale_processes) call temperature_large_scale_terms(Temperature_adv)

        if (TempBtype(To)==BC_AUTOMATIC_FLUX) then
          first = min(Prnz*5/6, Prnz-5)
          last = Prnz-5
          sideTemp(To) = sum(temperature_flux_profile(first:last)) / (last-first+1)

          do i = 1, size(WMPoints)
            if (WMPoints(i)%zk==Prnz) WMPoints(i)%temperature_flux = - sideTemp(To)
          end do
        end if

        do i = 1, size(WMPoints)
          associate (p => WMPoints(i))
            Temperature_adv(p%xi,p%yj,p%zk) = Temperature_adv(p%xi,p%yj,p%zk) + &
                                          p%temperature_flux * p%area_factor

            Temperature_adv(p%xi,p%yj,p%zk) = Temperature_adv(p%xi,p%yj,p%zk) + p%div * Temperature(p%xi,p%yj,p%zk)

          end associate
        end do
      end subroutine

      subroutine MoistureExtra
        use Large_scale_processes
        use VolumeSources, only: MoistureVolumeSources
        integer :: i, first, last

        call MoistureVolumeSources(Moisture_adv)
        
        if (enable_large_scale_processes) call moisture_large_scale_terms(Moisture_adv)

        if (MoistBtype(To)==BC_AUTOMATIC_FLUX) then
          first = min(Prnz*5/6, Prnz-5)
          last = Prnz-5
          sideMoist(To) = sum(moisture_flux_profile(first:last)) / (last-first+1)
        end if

        do i = 1, size(WMPoints)
          associate (p => WMPoints(i))
            Moisture_adv(p%xi,p%yj,p%zk) = Moisture_adv(p%xi,p%yj,p%zk) + &
                                          p%moisture_flux * p%area_factor

            Moisture_adv(p%xi,p%yj,p%zk) = Moisture_adv(p%xi,p%yj,p%zk) + p%div * Moisture(p%xi,p%yj,p%zk)

          end associate
        end do
      end subroutine

      subroutine ScalarExtra
        use VolumeSources, only: ScalarVolumeSources
        integer :: i

        call ScalarVolumeSources(Scalar_adv(:,:,:,sc), sc, RK_stage==1)

        call DoPuffs(Scalar_adv(:,:,:,sc), sc)

        do i = 1, size(WMPoints)
          associate (p => WMPoints(i))

            Scalar_adv(p%xi,p%yj,p%zk,sc) = Scalar_adv(p%xi,p%yj,p%zk,sc) + p%div * Scalar(p%xi,p%yj,p%zk,sc)

          end associate
        end do
      end subroutine

  end subroutine ScalarRK3







  pure real(knd) function AirDensity(press, temp)
    real(knd), intent(in) :: press, temp

    AirDensity = press / (287.05_knd * temp)
  endfunction AirDensity


  pure real(knd) function AirDynVisc(temp)
    real(knd), intent(in) :: temp

    AirDynVisc = 1.85e-5_knd
  endfunction AirDynVisc



  pure real(knd) function CorrFactor(dp, press, temp)
    real(knd), intent(in) :: dp, press, temp
    real(dbl) :: l

    l = MeanFreePath(press, temp)
    CorrFactor = real(1 + (2*l/dp) * (1.257_knd + 0.4_knd * exp(-0.55_knd*dp/l)), knd)
  endfunction CorrFactor


  pure real(dbl) function MeanFreePath(press, temp)
    real(knd), intent(in) :: press, temp

    MeanFreePath = 2.24e-5_dbl * temp / press
  endfunction MeanFreePath


  pure real(dbl) function BrownDiffusivity(dp, press, temp)
    real(knd), intent(in) :: dp, press, temp
    real(dbl) :: C

    C = Corrfactor(dp, press, temp)
    BrownDiffusivity = BoltzC * temp * C / (3 * pi * AirDynVisc(temp) * dp)
  endfunction BrownDiffusivity


  pure real(knd) function BrownEff(dp, press, temp)
    real(knd), intent(in) :: dp, press, temp
    real(knd) :: Sc

    Sc = real( (AirDynVisc(temp) / AirDensity(press, temp)) / &
               BrownDiffusivity(dp, press, temp), knd)
    BrownEff = Sc**(-0.54_knd)
  endfunction BrownEff


  pure real(knd) function ImpactEff(dp, press, temp, ustar, visc)
    real(knd), intent(in) :: dp, press, temp, ustar, visc
    real(knd) :: St

    St = SedimVelocity2(dp, press, temp) * ustar**2 / visc
    ImpactEff = St**2 / (400 + St**2)
  endfunction ImpactEff


  pure real(knd) function SedimVelocity2(dp, press, temp)
    real(knd), intent(in) :: dp, temp, press
    real(knd) :: C

    C = CorrFactor(dp, press, temp)
    SedimVelocity2 = AirDensity(press, temp)*dp**2 * 9.81_knd * C / (18 * AirDynVisc(temp))
  endfunction SedimVelocity2

  pure real(dbl) function SedimVelocity(dp, rhop, press, temp)
    real(knd), intent(in) :: dp, rhop, temp, press
    real(dbl) :: C, us, rho, mu

    rho = AirDensity(press, temp)
    mu = AirDynVisc(temp)
    C = CorrFactor(dp, press, temp)
    us = 1 + (0.42_dbl * C**2 * rho * rhop / (108 * mu**2)) * dp**3 * (1-rho/rhop) * 9.81_dbl
    us = sqrt(us)
    us = (12._dbl * mu / (0.42_dbl  *C * rho * dp)) * (us - 1._dbl)
    SedimVelocity = us
  endfunction SedimVelocity


  pure real(knd) function DyerH(zL)
    real(knd), intent(in) :: zL

    if (zL >= 0) then
     DyerH = 1._knd + 5._knd * zl
    else
     DyerH = 1._knd / sqrt(1._knd - 16._knd * zl)
    end if
  endfunction DyerH



  pure real(knd) function AerResist(z, z0, zL, ustar, visc)
    real(knd), intent(in) :: z, z0, zL, ustar, visc
    real(knd), parameter :: yplcrit = 11.225

    if (z>z0.and.z0>0) then
      AerResist = (log(z/z0) - DyerH(zl)) / (0.4_knd * ustar)
    else
      if (z * ustar / visc < yplcrit) then
        AerResist = z * ustar / visc
      else
       AerResist = (log(abs(ustar*z/visc)) / 0.4_knd + 5.2_knd) / ustar
      end if
    end if
  endfunction AerResist

  pure real(knd) function DepositionVelocity3(dp, rhop, press, temp, z, z0, zL, ustar) !EMRAS recommended values
    real(knd), intent(in) :: dp, rhop, press, temp, z, z0, zL, ustar

    if (dp < 1e-6) then
      DepositionVelocity3 = 0.5e-4
    else if (dp < 2e-6) then
      DepositionVelocity3 = 1.5e-4
    else if (dp < 10e-6) then
      DepositionVelocity3 = 10e-4
    else
      DepositionVelocity3 = 80e-4
    end if
  endfunction DepositionVelocity3



  pure real(knd) function DepositionVelocity2(dp, press, temp, z, z0, zL, ustar)
    real(knd), intent(in) :: dp, press, temp, z, z0, zL, ustar
    real(knd) :: SurfResist, St, visc

    visc = real(AirDynVisc(temp) / AirDensity(press,temp), knd)
    St = SedimVelocity2(dp, press, temp) * ustar**2 / (visc)
    SurfResist = 1 / (3 * ustar * (BrownEff(dp, press, temp) + ImpactEff(dp, press, temp, ustar, visc)))
    DepositionVelocity2 = SedimVelocity2(dp, press, temp)+1 / &
                          (AerResist(z, z0, zL, ustar, visc) + SurfResist)
  endfunction DepositionVelocity2


  pure real(knd) function DepositionVelocity(dp, rhop, press, temp, z, z0, zL, ustar) !Kharchenko
    real(knd), intent(in) :: dp, rhop, press, temp, z, z0, zL, ustar
    real(dbl) :: visc, us, Intz, Intexp, BD, tp
    real(dbl), parameter :: zexp = 0.01

    us = SedimVelocity(dp, rhop, press, temp)
    visc = AirDynVisc(temp) / AirDensity(press, temp)
    tp = (us/9.81_knd)  *ustar**2 / visc
    BD = BrownDiffusivity(dp, press, temp)

    if (zl>=0) then
      Intz = -log(z/zexp+6*(zl-zexp*zl/z)) / Karman
    else
      Intz = ((sqrt(1-9*zl)-1) * (sqrt(1-9*zexp*zl/z)+1))
      Intz = Intz/((sqrt(1-9*zl)+1) * (sqrt(1-9*zexp*zl/z)-1))
      Intz = -Intz/Karman
    end if

    Intexp = -367.8_knd
    Intexp = Intexp + 16.4 * log(visc/BD)
    Intexp = Intexp - 0.73 * log(100*dp) * log(1e4*BD) - (log(100*dp))**2 / 2
    Intexp = Intexp + 0.13 * log(0.03/z0)
    Intexp = Intexp + 0.25 * log(0.2/ustar) * (1 - 0.2*log(0.03/z0))
    Intexp = Intexp - 0.03 * log(tp)*log(0.03/z0)
    Intexp = Intexp - 32.7 * log(100*dp)
    Intexp = -exp(Intexp)

    DepositionVelocity = real(us/(1._knd-exp((us/ustar) * (Intexp+Intz))), knd)
  endfunction DepositionVelocity


  pure real(knd) function DepositionFlux(WMP, conc, partdiam, rhop)
    type(WMPoint), intent(in) :: WMP
    real(knd), intent(in) :: conc, partdiam, rhop
    real(knd) :: press, temp, depvel

    press = 101300
    temp = temperature_ref
    if ( (.2_knd*WMP%distz)**2 > (WMP%distx)**2+(WMP%disty)**2 .and. WMP%distz>0 ) then
      depvel = DepositionVelocity(partdiam, rhop, press, temp, WMP%distz, WMP%z0, 0._knd, WMP%ustar)
    else
      depvel = 0
    end if
    DepositionFlux = abs(depvel)*abs(conc)

  endfunction DepositionFlux

  subroutine Deposition(Scal, coef)
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(inout) :: Scal
    real(knd), intent(in) :: coef
    integer   :: i, j
    real(knd) :: deptmp
    
    do j = 1, size(WMPoints)
      associate(p => WMPoints(j))
        if (allocated(p%depscalar)) then
          if (partdistrib>0) then
            do i = 1, partdistrib
              deptmp = abs(DepositionFlux(p, Scal(p%xi,p%yj,p%zk,1)*percdistrib(i), partdiam(i), partrho(i))) &
                  * coef * dxPr(p%xi) * dyPr(p%yj)
              p%depscalar(1) = p%depscalar(1) + deptmp / (dxPr(p%xi)*dyPr(p%yj)*dzPr(p%zk))
              Scal(p%xi,p%yj,p%zk,1) = Scal(p%xi,p%yj,p%zk,1) - deptmp / (dxPr(p%xi)*dyPr(p%yj)*dzPr(p%zk))
            end do
          else
            do i = 1, num_of_scalars
              deptmp = abs(DepositionFlux(p, Scal(p%xi,p%yj,p%zk,i), partdiam(i), partrho(i))) &
                  * coef * dxPr(p%xi) * dyPr(p%yj)
              p%depscalar(i) = p%depscalar(i) + deptmp / (dxPr(p%xi)*dyPr(p%yj)*dzPr(p%zk))
              Scal(p%xi,p%yj,p%zk,i) = Scal(p%xi,p%yj,p%zk,i) - deptmp / (dxPr(p%xi)*dyPr(p%yj)*dzPr(p%zk))
            end do
          end if
        end if
      end associate
    end do
    
  endsubroutine Deposition


  subroutine Gravsettling(Scal, coef)
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(inout) :: Scal
    integer :: i, j, k, l
    real(knd), dimension(Prnx,Prny,Prnz) :: flux
    real(knd) :: coef, press, temp, us

    press = 101300
    temp = temperature_ref
    if (partdistrib==0) then
      do l = 1, num_of_scalars
        do k = 1, Prnz
          do j = 1, Prny
            do i = 1, Prnx
              us = real( SedimVelocity(partdiam(l), partrho(l), press, temp), knd)
              flux(i,j,k) = us * Scal(i,j,k+1,l) * coef * dxPr(i) * dyPr(j)
            end do
          end do
        end do
        do k = 1, Prnz-1
          do j = 1, Prny
            do i = 1, Prnx
              Scal(i,j,k+1,l) = Scal(i,j,k+1,l) - flux(i,j,k) / (dxPr(i)*dyPr(j)*dzPr(k+1))
              Scal(i,j,k  ,l) = Scal(i,j,k  ,l) + flux(i,j,k) / (dxPr(i)*dyPr(j)*dzPr(k))
            end do
          end do
        end do
      end do
    end if
  endsubroutine Gravsettling


  pure real(knd) function Rig(i, j, k, U, V, temperature)
    integer, intent(in) :: i, j, k
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: temperature
    real(knd) :: num, denom

    num = (grav_acc/temperature_ref) * (temperature(i,j,k+1)-temperature(i,j,k-1)) / &
          (zPr(k+1)-zPr(k-1))

    denom = ((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) / &
             (2 * (zPr(k+1)-zPr(k-1))))**2
    denom = denom + ((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) / &
                     (2 * (zPr(k+1)-zPr(k-1))))**2

    if (abs(denom)>1E-5_knd*abs(num) .and. abs(denom)>epsilon(1._knd)) then
      Rig = num / denom
    else
      Rig = 100000 * sign(1.0_knd, num) * sign(1.0_knd, denom)
    end if
  endfunction Rig



  subroutine Interpolate1DProfile(prof, prof_obj, lo, w_grid)
    use Interpolation
    real(knd), contiguous, intent(inout) :: prof(lo:)
    type(TScalarProfile), intent(inout) :: prof_obj
    integer, intent(in) :: lo
    logical, optional :: w_grid
    logical :: w_grid_loc
    real(knd), allocatable :: coefs(:,:)
    integer :: k, i_int, up
    real(knd) :: z
    
    w_grid_loc = .false.
    if (present(w_grid)) w_grid_loc = w_grid
    
    up = ubound(prof, 1)
    
    if (allocated(prof_obj%points)) then
    
      allocate(coefs(0:3,size(prof_obj%points,1)))
      call linear_interpolation(prof_obj%points(:,1), prof_obj%points(:,2), coefs)
      
      i_int = 0
      do k = lo, up
        if (w_grid_loc) then
          z = zW(k)
        else
          z = zPr(k)
        end if
        
        prof(k) = linear_interpolation_eval(z, prof_obj%points(:,1), coefs, i_int)
      end do
    else
      prof = 0
    end if
  end subroutine

  subroutine InitScalarProfile(ScalarIn, ScalarProfile, default_value)
    use Interpolation
    real(knd), contiguous, intent(inout) :: ScalarIn(-2:,-2:)
    type(TScalarProfile), intent(inout) :: ScalarProfile
    real(knd), intent(in) :: default_value
    integer   :: SectionToUse(lbound(ScalarIn,2):ubound(ScalarIn,2))
    integer   :: section, nSections, s
    integer   :: i, j, k
    real(knd) :: temp
    real(knd), allocatable :: coefs(:,:)
    
    if (allocated(ScalarProfile%points)) then
    
      allocate(coefs(0:3,size(ScalarProfile%points,1)))
      call linear_interpolation(ScalarProfile%points(:,1), ScalarProfile%points(:,2), coefs)
      
      i = 0
      do k = -1, Prnz+2
        ScalarIn(:,k) = linear_interpolation_eval(zPr(k), ScalarProfile%points(:,1), coefs, i)
      end do
        
    else

      nSections = size(ScalarProfile%Sections)

      if (nSections > 0) then

         if (size(ScalarProfile%Sections)>0) then
           if (ScalarProfile%Sections(1)%jump<=0) ScalarProfile%Sections(1)%jump = default_value
         end if

        ScalarProfile%Sections(1)%height = ScalarProfile%Sections(1)%top

        do i = 2, nSections
          if (ScalarProfile%Sections(i)%top < ScalarProfile%Sections(i-1)%top) &
            ScalarProfile%Sections(i)%top = ScalarProfile%Sections(i-1)%top

          ScalarProfile%Sections(i)%height = ScalarProfile%Sections(i)%top - ScalarProfile%Sections(i-1)%top
        end do

        if (ScalarProfile%Sections(nSections)%top < zW(Prnz+2) ) &
            ScalarProfile%Sections(nSections)%top = zW(Prnz+2)

        section = 1

        do k = -1, Prnz+2
          do while (zPr(k) > ScalarProfile%Sections(section)%top) !should be safe, because last top is adjusted above
            section = section + 1
          end do
          SectionToUse(k) = section
        end do

      else

        SectionToUse = 0

      end if

      do k = -1, Prnz+2

        s = SectionToUse(k)

        if (s==0) then

          temp = default_value

        else

          temp = 0

          do i = 1, s-1
            temp = temp + ScalarProfile%Sections(i)%jump
            temp = temp + ScalarProfile%Sections(i)%height * ScalarProfile%Sections(i)%gradient
          end do

          temp = temp + ScalarProfile%Sections(s)%jump

          if (s>1) then
            temp = temp + (zPr(k) - ScalarProfile%Sections(s-1)%top) * ScalarProfile%Sections(s)%gradient
          else
            temp = temp + zPr(k) * ScalarProfile%Sections(s)%gradient
          end if

        end if

        do j = -1, Prny+2
            ScalarIn(j,k) = temp
        end do

      end do
    
    end if
  end subroutine InitScalarProfile


  subroutine InitScalar(ScalarIn, ScalarProfile, Sc)
    use rng_par_zig
    !$ use omp_lib
    real(knd), contiguous, intent(in)  :: ScalarIn(-2:,-2:)
    type(TScalarProfile), intent(in)  :: ScalarProfile
    real(knd), contiguous, intent(out) :: Sc(-2:,-2:,-2:)
    real(knd) :: p
    integer   :: i, j, k, tid

    if (ScalarProfile%randomize) then
      tid = 0
      !$omp parallel private(i,j,k,p,tid)
      !$ tid = omp_get_thread_num()
      
      !$omp do
      do k = 0, Prnz+1
        do j = 0, Prny+1
           do i = 0, Prnx+1

             if (zPr(k) <= ScalarProfile%randomizeTop) then
               call rng_uni(p, tid)
               p = p - 0.5
             else
               p = 0
             end if

             Sc(i,j,k)=ScalarIn(j,k) + ScalarProfile%randomizeAmplitude * 2 * p

           end do
        end do
      end do
      !$omp end parallel

    else
      !$omp parallel workshare
      forall(i=0:Prnx+1) Sc(i,:,:) = ScalarIn(:,:)
      !$omp end parallel workshare
    end if

  end subroutine InitScalar


  subroutine InitSubsidenceProfile
    use Large_scale_processes
    !the subsidence profile contains velocities <w>
    !positive means upwards
    integer :: k

    if (allocated(SubsidenceProfileObj%points)) then
    
      enable_subsidence_profile = .true. 
    
      allocate(subsidence_profile(1:Prnz))
      call Interpolate1DProfile(subsidence_profile, SubsidenceProfileObj, 1)  
      
    else if (SubsidenceGradient/=0) then
    
      enable_subsidence_profile = .true. 
    
      allocate(subsidence_profile(0:Prnz))
      subsidence_profile = [ (- zPr(k) * SubsidenceGradient, k=0, Prnz) ]
      
    else
    
      enable_subsidence_profile = .false. 
          
    end if
  end subroutine InitSubsidenceProfile


end module Scalars
