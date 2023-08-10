module TimeSteps

  use Parameters
  use Dynamics
  use Boundaries, only: BoundUVW, Bound_Q
  use Pressure, only: PressureCorrection, pressure_solution, POISSON_SOLVER_NONE
  use Outputs, only: current_profiles
  use Scalars, only: ScalarRK3
  use Turbinlet, only: default_turbulence_generator, GetInletFromFile
  use Sponge, only: enable_top_sponge, enable_out_sponge_x, enable_out_sponge_y, &
                    SpongeTop, SpongeOut

  implicit none


  private
  public TMarchRK3


contains


  subroutine TMarchRK3(U, V, W, Pr, Temperature, Moisture, Scalar, delta)
    use RK3
#ifdef PAR
    use custom_par
    use domains_bc_par
#endif
    real(knd), allocatable, intent(inout) :: U(:,:,:), V(:,:,:) ,W(:,:,:), Pr(:,:,:)
    real(knd), allocatable, intent(inout) :: Temperature(:,:,:), Moisture(:,:,:), Scalar(:,:,:,:)
    real(knd), intent(out) :: delta

    integer :: RK_stage
    integer, save :: called = 0
    integer(int64), save :: trate
    integer(int64), save :: time1, time2

#ifdef  CUSTOM_TIMESTEP_PROCEDURE
    interface
      subroutine CustomTimeStepProcedure
      end subroutine
    end interface
#endif

    time_stepping%effective_time = time_stepping%time

#ifdef PAR
    !uses previous dt (it should be fixed anyway, but...)
    call par_exchange_domain_bounds(U, V, W, Temperature, Moisture, Scalar, &
                                    time_stepping%time, time_stepping%dt)
#endif

    if (called==0) then
      called = 1

      !just to allocate it and make it defined in all points
      U2 = U
      V2 = V
      W2 = W

      Ustar = U
      Vstar = V
      Wstar = W


      call BoundUVW(U, V, W)

      if (enable_buoyancy) call BoundTemperature(temperature)

      call IBMomentum(U, V, W)


      if (enable_ibm_mass_sources) allocate(Q(0:Prnx+1, 0:Prny+1, 0:Prnz+1))

      if (debugparam>1) call system_clock(count_rate=trate)
    end if


    if (any(Btype==BC_TURBULENT_INLET)) then
      call default_turbulence_generator%time_step(Uin, Vin, Win, time_stepping%dt)
    else if (Btype(We)==BC_INLET_FROM_FILE) then
      call GetInletFromFile(time_stepping%time)
    end if

    call TimeStepLength(U, V, W, time_stepping)

    if (master) then
      if (time_stepping%variable_time_steps) then
        write (*,'(a,f12.6,a,es12.4)') "  dt: ", time_stepping%dt
      else
        if (time_stepping%enable_CFL_check) then
#ifdef PAR
          if (enable_multiple_domains) then
            write (*,'(a,i0,a,f6.3)') "  domain: ", domain_index,"   CFL: ", time_stepping%CFL
          else
            write (*,'(a,f6.3)') "  CFL: ", time_stepping%CFL
          end if
#else
          write (*,'(a,f6.3)') "  CFL: ", time_stepping%CFL
#endif          
        end if
      end if
    end if
    
#ifdef  CUSTOM_TIMESTEP_PROCEDURE
    call CustomTimeStepProcedure
#endif


    do RK_stage = 1, RK_stages

      if (debugparam>1.and.master) call system_clock(count=time1)

      time_stepping%effective_time = time_stepping%time + 2 * sum(RK_alpha(1:RK_stage-1)) * time_stepping%dt      

#ifdef PAR
      call par_update_pr_gradient(time_stepping%effective_time)
#endif

      call SubgridStresses(U, V, W, Pr, Temperature, Moisture)


      call Convection(U, V, W, &
                      U2, V2, W2, &
                      Ustar, Vstar, Wstar, &
                      Temperature, Moisture, Scalar, &
                      Pr, &
                      RK_beta, RK_rho, RK_stage, time_stepping%dt)

      call ScalarRK3(U, V, W, Pr, &
                     Temperature, Moisture, Scalar, &
                     RK_stage, time_stepping%dt, &
                     current_profiles%tempfl, current_profiles%moistfl)

      call OtherTerms(U, V, W, &
                      U2, V2, W2, &
                      Pr, &
                      2*RK_alpha(RK_stage)*time_stepping%dt)

      if (enable_top_sponge)  then

          call SpongeTop(U2, V2, W2)
      end if

      if (enable_out_sponge_x .or. enable_out_sponge_y) then

          call SpongeOut(U2, V2, W2, temperature)
      end if

      if (enable_fixed_flow_rate) call CorrectFlowRate(U2, V2, W2)

      call BoundUVW(U2, V2, W2)


      time_stepping%effective_time = time_stepping%time + 2 * sum(RK_alpha(1:RK_stage)) * time_stepping%dt      

#ifdef PAR
      !Does nudging of the solution to the parent boundary solution in the relaxation zones.
      call par_domain_bound_relaxation(U2, V2, W2, Temperature, Moisture, Scalar, time_stepping%effective_time)
#endif


#ifdef PAR
      if (RK_stage==RK_stages) call par_domain_two_way_nesting_feedback(U2, V2, W2, Temperature, Moisture, Scalar, &
                                            time_stepping%time, time_stepping%dt)
#endif

      call IBMomentum(U2, V2, W2)

      call BoundUVW(U2, V2, W2)

      if (enable_ibm_mass_sources) then
          call IBMassSources(Q, U2, V2, W2)
      end if



      if (pressure_solution%poisson_solver > POISSON_SOLVER_NONE) then
        call PressureCorrection(U2, V2, W2, Pr, Q, 2*RK_alpha(RK_stage)*time_stepping%dt)
      end if



      if (RK_stage==1) delta = 0

      if ( debuglevel>0 .or. steady==1 ) then

        if (Unx*Uny*Unz > 0) &
          delta = delta + sum(abs(U(1:Unx,1:Uny,1:Unz) - U2(1:Unx,1:Uny,1:Unz))) / (Unx*Uny*Unz)

        if (Vnx*Vny*Vnz > 0) &
          delta = delta + sum(abs(V(1:Vnx,1:Vny,1:Vnz) - V2(1:Vnx,1:Vny,1:Vnz))) / (Vnx*Vny*Vnz)

        if (Wnx*Wny*Wnz > 0) &
          delta = delta + sum(abs(W(1:Wnx,1:Wny,1:Wnz) - W2(1:Wnx,1:Wny,1:Wnz))) / (Wnx*Wny*Wnz)

        if (RK_stage==RK_stages) then
#ifdef PAR
          delta = par_co_sum(delta)
#endif
          if (master) write(*,*) "delta",delta
        end if
      end if

      call exchange_alloc(U, U2)
      call exchange_alloc(V, V2)
      call exchange_alloc(W, W2)



      if (enable_out_sponge_x .or. enable_out_sponge_y) then

        call SpongeOut(U, V, W, temperature)

      end if


      call NullInterior(U, V, W)


      call BoundUVW(U, V, W)

      call IBMomentum(U, V, W)

      call BoundUVW(U, V, W)


      if (debugparam>1.and.master) then
        call system_clock(count=time2)
        write (*,*) "ET of part 1", (time2-time1) / real(trate, dbl)
        time1 = time2
      end if

    end do
    

  end subroutine TMarchRK3



















  subroutine OtherTerms(U, V, W, U2, V2, W2, Pr, coef)
    real(knd), contiguous,  intent(inout) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd), contiguous,  intent(inout) :: Pr(-1:,-1:,-1:)
    real(knd), allocatable, intent(inout) :: U2(:,:,:), V2(:,:,:), W2(:,:,:)
    real(knd), intent(in) :: coef

    real(knd) :: S

    integer :: i, j, k
    integer, save :: called=0

    if (called==0) then
      allocate(U3(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)))
      allocate(V3(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)))
      allocate(W3(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)))
      called=1
    end if


    call PressureGrad(Pr, U2, V2, W2, coef)

    if (.not.explicit_diffusion) then

        !semi-implicit diffusion

        Re_gt_0: if (molecular_viscosity > 0) then

          !Diffusion using Crank Nicolson
          !first approximation using forward Euler
          !iteration SOR or Gauss-Seidel

          call ImplicitDiffusion_ForwEul(U, V, W, U2, V2, W2, U3, V3, W3, coef)


          call BoundUVW(U3, V3, W3)

          call IBMomentum(U3, V3, W3)

          !Performs the diffusion terms

          if (gridtype==GRID_UNIFORM) then
            call ImplicitDiffusion_Iterations(U, V, W, U2, V2, W2, U3, V3, W3, coef)
          else
            call error_stop("Non-uniform grid support was dropped.")
          end if

          call exchange_alloc(U2, U3)
          call exchange_alloc(V2, V3)
          call exchange_alloc(W2, W3)


        else  Re_gt_0  !Re<=0

          U2 = U + U2
          V2 = V + V2
          W2 = W + W2

        end if   Re_gt_0


        call BoundUVW(U2, V2, W2)

        call IBMomentum(U2, V2, W2)

        if (debuglevel>=2) then  !Compute and output the mean friction in the domain.
          S = 0
          do k = 1, Unz
           do j = 1, Uny
            do i = 1, Unx
              S = S-((Viscosity(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) / dxPr(i+1) - &
              Viscosity(i,j,k) * (U(i,j,k)-U(i-1,j,k)) / dxPr(i)) / dxU(i) + &
                (0.25_knd * (Viscosity(i+1,j+1,k)+Viscosity(i+1,j,k)+Viscosity(i,j+1,k)+Viscosity(i,j,k))* &
                      (U(i,j+1,k)-U(i,j,k)) / dyV(j) - &
                0.25_knd * (Viscosity(i+1,j,k)+Viscosity(i+1,j-1,k)+Viscosity(i,j,k)+Viscosity(i,j-1,k))* &
                      (U(i,j,k)-U(i,j-1,k)) / dyV(j-1)) / dyPr(j) + &
                 (0.25_knd * (Viscosity(i+1,j,k+1)+Viscosity(i+1,j,k)+Viscosity(i,j,k+1)+Viscosity(i,j,k))* &
                      (U(i,j,k+1)-U(i,j,k)) / dzW(k) - &
                0.25_knd * (Viscosity(i+1,j,k)+Viscosity(i+1,j,k-1)+Viscosity(i,j,k)+Viscosity(i,j,k-1))* &
                      (U(i,j,k)-U(i,j,k-1)) / dzW(k-1)) / dzPr(k))
            end do
           end do
          end do

          S = S / (Unx*Uny*Unz)
          write(*,*) "Mean friction:", S
        end if

    else !explicit diffusion

        call add(U2, U)
        call add(V2, V)
        call add(W2, W)

    end if

  end subroutine OtherTerms








  
  subroutine IBMomentum(U, V, W)
    use ImmersedBoundary, only: Up => UIBPoints, &
                                Vp => VIBPoints, &
                                Wp => WIBPoints, &
                                Interpolate => TIBPoint_Interpolate

    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    integer :: i

    if (size(Up) + size(Vp) + size(Wp)>0) then
      if (.not. allocated(Uwm)) allocate(Uwm(size(Up)))
      if (.not. allocated(Vwm)) allocate(Vwm(size(Vp)))
      if (.not. allocated(Wwm)) allocate(Wwm(size(Wp)))

      !$omp parallel
      !$omp do
      do i = 1, ubound(Up, 1)
          Uwm(i) = Interpolate(Up(i), U, -2)
      end do
      !$omp end do
      !$omp do
      do i = 1, ubound(Up, 1)
          U(Up(i)%xi, Up(i)%yj, Up(i)%zk) = Uwm(i)
      end do
      !$omp end do nowait

      !$omp do
      do i = 1, ubound(Vp, 1)
          Vwm(i) = Interpolate(Vp(i), V, -2)
      end do
      !$omp end do
      !$omp do
      do i = 1, ubound(Vp, 1)
          V(Vp(i)%xi, Vp(i)%yj, Vp(i)%zk) = Vwm(i)
      end do
      !$omp end do nowait

      !$omp do
      do i = 1, ubound(Wp, 1)
          Wwm(i) = Interpolate(Wp(i), W, -2)
      end do
      !$omp end do
      !$omp do
      do i = 1, ubound(Wp, 1)
          W(Wp(i)%xi, Wp(i)%yj, Wp(i)%zk) = Wwm(i)
      end do
      !$omp end do
      !$omp end parallel

    end if
  end subroutine IBMomentum






  subroutine IBMassSources(Q, U, V, W)
    use ImmersedBoundary, only: Up => UIBPoints, &
                                Vp => VIBPoints, &
                                Wp => WIBPoints
    use vtkarray

    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
    real(knd), intent(out) :: Q(0:,0:,0:)
    integer :: i, xi, yj, zk
    !deconvolution of the flux velocity
    real(knd), parameter :: C0 = 26._knd / 24, C1 = -1._knd / 24
    real(knd) :: tmp

    call set(Q, 0)
    
    if (discretization_order==4) then
    
      !$omp parallel
      !$omp do private(xi, yj, zk, tmp)
      do i = 1, ubound(Up, 1)
        xi = Up(i)%xi
        yj = Up(i)%yj
        zk = Up(i)%zk
        
        tmp = (C1*U(xi-1,yj,zk) + C0*U(xi,yj,zk) + C1*U(xi+1,yj,zk)) / dxmin
        
        !$omp atomic
        Q(xi,yj,zk)   = Q(xi,yj,zk)   + tmp
        !$omp atomic
        Q(xi+1,yj,zk) = Q(xi+1,yj,zk) - tmp
      end do
      !$omp end do

      !$omp do private(xi, yj, zk, tmp)
      do i = 1, ubound(Vp, 1)
        xi = Vp(i)%xi
        yj = Vp(i)%yj
        zk = Vp(i)%zk

        tmp = (C1*V(xi,yj-1,zk) + C0*V(xi,yj,zk) + C1*V(xi,yj+1,zk)) / dymin
        
        !$omp atomic
        Q(xi,yj,zk)   = Q(xi,yj,zk)   + tmp
        !$omp atomic
        Q(xi,yj+1,zk) = Q(xi,yj+1,zk) - tmp
      end do
      !$omp end do


      !$omp do private(xi, yj, zk, tmp)
      do i = 1, ubound(Wp, 1)
        xi = Wp(i)%xi
        yj = Wp(i)%yj
        zk = Wp(i)%zk

        tmp = (C1*W(xi,yj,zk-1) + C0*W(xi,yj,zk) + C1*W(xi,yj,zk+1)) / dzmin
        
        !$omp atomic
        Q(xi,yj,zk)   = Q(xi,yj,zk)   + tmp
        !$omp atomic
        Q(xi,yj,zk+1) = Q(xi,yj,zk+1) - tmp
      end do
      !$omp end do
      !$omp end parallel
      
    else

      !$omp parallel
      !$omp do private(xi, yj, zk)
      do i = 1, ubound(Up, 1)
        xi = Up(i)%xi
        yj = Up(i)%yj
        zk = Up(i)%zk
        !$omp atomic
        Q(xi,yj,zk)   = Q(xi,yj,zk)   + U(xi,yj,zk) / dxPr(xi)
        !$omp atomic
        Q(xi+1,yj,zk) = Q(xi+1,yj,zk) - U(xi,yj,zk) / dxPr(xi+1)
      end do
      !$omp end do

      !$omp do private(xi, yj, zk)
      do i = 1, ubound(Vp, 1)
        xi = Vp(i)%xi
        yj = Vp(i)%yj
        zk = Vp(i)%zk

        !$omp atomic
        Q(xi,yj,zk)   = Q(xi,yj,zk)   + V(xi,yj,zk) / dyPr(yj)
        !$omp atomic
        Q(xi,yj+1,zk) = Q(xi,yj+1,zk) - V(xi,yj,zk) / dyPr(yj+1)
      end do
      !$omp end do


      !$omp do private(xi, yj, zk)
      do i = 1, ubound(Wp, 1)
        xi = Wp(i)%xi
        yj = Wp(i)%yj
        zk = Wp(i)%zk

        !$omp atomic
        Q(xi,yj,zk)   = Q(xi,yj,zk)   + W(xi,yj,zk) / dzPr(zk)
        !$omp atomic
        Q(xi,yj,zk+1) = Q(xi,yj,zk+1) - W(xi,yj,zk) / dzPr(zk+1)
      end do
      !$omp end do
      !$omp end parallel
      
    end if

    call Bound_Q(Q)


  end subroutine IBMassSources




  subroutine NullInterior(U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    integer :: i

    !$omp parallel private(i)
    !$omp do
    do i = 1, nUnull
      U(Unull(1,i),Unull(2,i),Unull(3,i)) = 0
    end do
    !$omp end do nowait

    !$omp do
    do i = 1, nVnull
      V(Vnull(1,i),Vnull(2,i),Vnull(3,i)) = 0
    end do
    !$omp end do nowait

    !$omp do
    do i = 1, nWnull
      W(Wnull(1,i),Wnull(2,i),Wnull(3,i)) = 0
    end do
    !$omp end do
    !$omp end parallel

  end subroutine NullInterior


end module TimeSteps










