module Pressure

  use Parameters
  use Boundaries
  use PoissonResidues
  use PoissonSolvers

  implicit none

  private
  public InitPressureCorrection, PressureCorrection, InitHydrostaticPressure, &
          pressure_solution, pressure_solution_control
  
  real(knd) :: bound_area(6)
  real(knd) :: bound_cell_area(6)
  real(knd) :: correction_area
  
  integer, parameter, public :: POISSON_SOLVER_NONE = 0, &
                                POISSON_SOLVER_SOR = 1, &
                                POISSON_SOLVER_POISFFT = 2, &
                                POISSON_SOLVER_MULTIGRID = 3
                                
  ! incremental - the Poisson equation solves for the inceremental pressure update, momentum equations use last stage pressure gradient
  ! pressure - the Poisson equation solves for the pressure field itself, no pressure gradient in momentum equations
  integer, parameter, public :: PROJECTION_METHOD_INCREMENTAL = 0, &
                                PROJECTION_METHOD_PRESSURE = 1

  type pressure_solution_control
    logical :: check_mass_flux = .false.
    logical :: report_mass_flux = .false.
    logical :: correct_mass_flux(6) = .false.
    integer :: poisson_solver = POISSON_SOLVER_POISFFT
    integer :: projection_method = PROJECTION_METHOD_INCREMENTAL
    logical :: check_divergence = .false.
    logical :: check_poisson_residue = .false.
    real(knd) :: top_pressure = 0    !mean pressure at the top boundary - calculated
    real(knd) :: bottom_pressure = 0
  end type
  
  type(pressure_solution_control) :: pressure_solution


contains

  subroutine InitPressureCorrection
    use custom_par
    
    integer :: bound_n_free(6)
    integer :: i, j, k

    bound_cell_area(We) = dymin * dzmin
    bound_cell_area(Ea) = dymin * dzmin
    bound_cell_area(So) = dxmin * dzmin
    bound_cell_area(No) = dxmin * dzmin
    bound_cell_area(Bo) = dymin * dxmin
    bound_cell_area(To) = dymin * dxmin
    
    
    if (gridtype == GRID_UNIFORM) then
      !no reason to save a few CPU cycles, but we save some numerical accuracy when summing up integers
        
      bound_n_free = 0

      if (iim==1     .and. Btype(We)/=BC_NOSLIP .and. Btype(We)/=BC_PERIODIC) &
        bound_n_free(We) = count(Utype(0,1:Uny,1:Unz)<=0)

      if (iim==nxims .and. Btype(Ea)/=BC_NOSLIP .and. Btype(Ea)/=BC_PERIODIC) &
        bound_n_free(Ea) = count(Utype(Prnx,1:Uny,1:Unz)<=0)

      if (jim==1     .and. Btype(So)/=BC_NOSLIP .and. Btype(So)/=BC_PERIODIC) &
        bound_n_free(So) = count(Vtype(1:Vnx,0,1:Vnz)<=0)

      if (jim==nyims .and. Btype(No)/=BC_NOSLIP .and. Btype(No)/=BC_PERIODIC) &
        bound_n_free(No) = count(Vtype(1:Vnx,Prny,1:Vnz)<=0)

      if (kim==1     .and. Btype(Bo)/=BC_NOSLIP .and. Btype(Bo)/=BC_PERIODIC) &
        bound_n_free(Bo) = count(Wtype(1:Wnx,1:Wny,0)<=0)

      if (kim==nzims .and. Btype(To)/=BC_NOSLIP .and. Btype(To)/=BC_PERIODIC) &
        bound_n_free(To) = count(Wtype(1:Wnx,1:Wny,Prnz)<=0)

#ifdef PAR
      bound_n_free = par_co_sum(bound_n_free)
#endif

      bound_area = bound_cell_area * bound_n_free

    else
    
      bound_area = 0
      
      if (iim==1     .and. Btype(We)/=BC_NOSLIP .and. Btype(We)/=BC_PERIODIC) then
        do k = 1, Prnz
          do j = 1, Prny
            if (Utype(0,j,k)<=0) bound_area(We) = bound_area(We) + dyPr(j) * dzPr(k)
          end do
        end do
      end if
    
      if (iim==nxims .and. Btype(Ea)/=BC_NOSLIP .and. Btype(Ea)/=BC_PERIODIC) then
        do k = 1, Prnz
          do j = 1, Prny
            if (Utype(Prnx,j,k)<=0) bound_area(Ea) = bound_area(Ea) + dyPr(j) * dzPr(k)
          end do
        end do
      end if
    
      if (jim==1     .and. Btype(So)/=BC_NOSLIP .and. Btype(So)/=BC_PERIODIC) then
        do k = 1, Prnz
          do i = 1, Prnx
            if (Vtype(i,0,k)<=0) bound_area(So) = bound_area(So) + dxPr(i) * dzPr(k)
          end do
        end do
      end if
    
      if (jim==nyims .and. Btype(No)/=BC_NOSLIP .and. Btype(No)/=BC_PERIODIC) then
        do k = 1, Prnz
          do i = 1, Prnx
            if (Vtype(i,Prny,k)<=0) bound_area(No) = bound_area(No) + dxPr(i) * dzPr(k)
          end do
        end do
      end if
    
      if (kim==1     .and. Btype(Bo)/=BC_NOSLIP .and. Btype(Bo)/=BC_PERIODIC) then
        do j = 1, Prny
          do i = 1, Prnx
            if (Wtype(i,j,0)<=0) bound_area(Bo) = bound_area(Bo) + dxPr(i) * dyPr(j)
          end do
        end do
      end if
    
      if (kim==nzims .and. Btype(To)/=BC_NOSLIP .and. Btype(To)/=BC_PERIODIC) then
        do j = 1, Prny
          do i = 1, Prnx
            if (Wtype(i,j,Prnz)<=0) bound_area(To) = bound_area(To) + dxPr(i) * dyPr(j)
          end do
        end do
      end if
    
#ifdef PAR
      bound_area = par_co_sum(bound_area)
#endif

    end if
    
    correction_area = sum(bound_area, mask=pressure_solution%correct_mass_flux)

  end subroutine


  subroutine PressureCorrection(U,V,W,Pr,Q,coef)                    !Pressure correction
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout)     :: U, V, W !Phi is computed in Poisson eq. with div of U in RHS
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(inout)        :: Pr      !Depend ing on active projection method Phi becomes new pressure
    real(knd), dimension(:,:,:), allocatable, intent(in) :: Q       !or is added to last pressure
    real(knd), intent(in) :: coef
                                                           !U,V,W velocity field for correction
    real(knd), save, allocatable :: Phi(:,:,:), RHS(:,:,:) !Pr pressure
                                                           !coef cofficient from Runge Kutta, Q mass sources from immersed boundary
    real(TIM) :: dt2,dt3                                      !RHS right hand side of eq. with divergence of U
                                                           !Phi computed pseudopressure, saved as first guess for next time
    real(knd) :: mass_flux
    real(knd) :: res
    integer, save :: called = 0
    integer(int64), save :: trate
    integer(int64), save :: time1, time2, time3, time4

    if (called==0) then
      allocate(Phi(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
      allocate(RHS(0:Prnx+1,0:Prny+1,0:Prnz+1))
      Phi = 0

      if (debugparam>1) call system_clock(count_rate=trate)
    end if
    called = called + 1


    if (debugparam>1 .and. called>1) call system_clock(count=time1)

    dt2 = coef
    dt3 = coef / 2


    call PrePoisson(U,V,W,Q,RHS,dt2,mass_flux)


    if (debugparam>1 .and. called>1) call system_clock(count=time2)

    if (pressure_solution%report_mass_flux .and. pressure_solution%check_mass_flux) then
      if (master) write(*,*) "total mass flux:", mass_flux
    end if

    if (pressure_solution%poisson_solver==POISSON_SOLVER_SOR) then

        call PoissSOR(Phi,RHS)

    else if (pressure_solution%poisson_solver==POISSON_SOLVER_POISFFT) then

        if (gridtype==GRID_VARIABLE_Z) then
          call Poiss_PoisFFT_variable_z(Phi,RHS)
        else
          call Poiss_PoisFFT(Phi,RHS)
        end if

    else if (pressure_solution%poisson_solver==POISSON_SOLVER_MULTIGRID) then

        if (Prny==1) then
          call PoissMG2d(Phi,RHS)
        else
          call PoissMG(Phi,RHS)
        end if

    end if


    call Bound_Phi(Phi)
    
    
    if (pressure_solution%check_poisson_residue) then
      call Poiss_Residue(Phi, RHS, res)
      
      if (master) write(*,*) "Poisson equation residue:", res
    end if


    if (debugparam>1 .and. called>1) then
     call system_clock(count=time3)
     if (master) write(*,*) "ET of part 2", real(time3-time2)/real(trate)
    endif


    call PostPoisson(U,V,W,Pr,Q,Phi,dt2,dt3)


    if (debugparam>1 .and. called>1) then
     call system_clock(count=time4)
     if (master) write(*,*) "ET of part 3", real((time4-time1)-(time3-time2))/real(trate)
    endif
    
  end subroutine PressureCorrection





  subroutine PrePoisson(U,V,W,Q,RHS,dt2,mass_flux)
    use custom_par
    real(knd), contiguous, intent(inout) :: U(-2:,-2:,-2:)
    real(knd), contiguous, intent(inout) :: V(-2:,-2:,-2:)
    real(knd), contiguous, intent(inout) :: W(-2:,-2:,-2:)
    real(knd), contiguous, intent(out)   :: RHS(0:,0:,0:)
    real(knd), allocatable, intent(in) :: Q(:,:,:)
    real(knd), intent(out)   :: mass_flux
    real(knd), intent(in)    :: dt2
    real(knd) :: dt2_rec
    integer   :: i, j, k
    integer   ::  side
    real(knd) :: flux(6)
    real(knd) :: df_side(6), df
    real(knd), parameter :: C1 = 9._knd / 8, C3 = 1._knd / (8*3)

    dt2_rec = 1._knd / dt2

    call BoundUVW(U, V, W)

    flux = 0
    if (pressure_solution%check_mass_flux) then
    
      if (gridtype==GRID_VARIABLE_Z) then
        if (iim==1 .and. Btype(We)/=BC_NOSLIP .and. Btype(We)/=BC_PERIODIC) then
          do k = 1, Unz
            do j = 1, Uny
              if (Utype(0,j,k)<=0) flux(We) = flux(We) + U(0,j,k) * dymin * dzPr(k)
            end do
          end do
        end if
        if (iim==nxims .and. Btype(Ea)/=BC_NOSLIP .and. Btype(Ea)/=BC_PERIODIC) then
          do k = 1, Unz
            do j = 1, Uny
              if (Utype(Prnx,j,k)<=0) flux(Ea) = flux(Ea) - U(Prnx,j,k) * dymin * dzPr(k)
            end do
          end do
        end if
        if (jim==1 .and. Btype(So)/=BC_NOSLIP .and. Btype(So)/=BC_PERIODIC) then
          do k = 1, Vnz
            do i = 1, Vnx
              if (Vtype(i,0,k)<=0) flux(So) = flux(So) + V(i,0,k) * dxmin * dzPr(k)
            end do
          end do
        end if
        if (jim==nyims .and. Btype(No)/=BC_NOSLIP .and. Btype(No)/=BC_PERIODIC) then
          do k = 1, Vnz
            do i = 1, Vnx
              if (Vtype(i,Prny,k)<=0) flux(No) = flux(No) - V(i,Prny,k) * dxmin * dzPr(k)
            end do
          end do
        end if
        if (kim==1 .and. Btype(Bo)/=BC_NOSLIP .and. Btype(Bo)/=BC_PERIODIC) then
          do j = 1, Wny
            do i = 1, Wnx
              if (Wtype(i,j,0)<=0) flux(Bo) = flux(Bo) + W(i,j,0) * dxmin * dymin
            end do
          end do
        end if
        if (kim==nzims .and. Btype(To)/=BC_NOSLIP .and. Btype(To)/=BC_PERIODIC) then
          do j = 1, Wny
            do i = 1, Wnx
              if (Wtype(i,j,Prnz)<=0) flux(To) = flux(To) - W(i,j,Prnz) * dxmin * dymin
            end do
          end do
        end if
        
      else
      
        if (iim==1 .and. Btype(We)/=BC_NOSLIP .and. Btype(We)/=BC_PERIODIC) then
          do k = 1, Unz
            do j = 1, Uny
              if (Utype(0,j,k)<=0) flux(We) = flux(We) + U(0,j,k)
            end do
          end do
        end if
        if (iim==nxims .and. Btype(Ea)/=BC_NOSLIP .and. Btype(Ea)/=BC_PERIODIC) then
          do k = 1, Unz
            do j = 1, Uny
              if (Utype(Prnx,j,k)<=0) flux(Ea) = flux(Ea) - U(Prnx,j,k)
            end do
          end do
        end if
        if (jim==1 .and. Btype(So)/=BC_NOSLIP .and. Btype(So)/=BC_PERIODIC) then
          do k = 1, Vnz
            do i = 1, Vnx
              if (Vtype(i,0,k)<=0) flux(So) = flux(So) + V(i,0,k)
            end do
          end do
        end if
        if (jim==nyims .and. Btype(No)/=BC_NOSLIP .and. Btype(No)/=BC_PERIODIC) then
          do k = 1, Vnz
            do i = 1, Vnx
              if (Vtype(i,Prny,k)<=0) flux(No) = flux(No) - V(i,Prny,k)
            end do
          end do
        end if
        if (kim==1 .and. Btype(Bo)/=BC_NOSLIP .and. Btype(Bo)/=BC_PERIODIC) then
          do j = 1, Wny
            do i = 1, Wnx
              if (Wtype(i,j,0)<=0) flux(Bo) = flux(Bo) + W(i,j,0)
            end do
          end do
        end if
        if (kim==nzims .and. Btype(To)/=BC_NOSLIP .and. Btype(To)/=BC_PERIODIC) then
          do j = 1, Wny
            do i = 1, Wnx
              if (Wtype(i,j,Prnz)<=0) flux(To) = flux(To) - W(i,j,Prnz)
            end do
          end do
        end if
        
        flux = flux * bound_cell_area
        
      end if


#ifdef PAR
      !TODO: only the relevant planes
      flux = par_co_sum(flux)
#endif

      df = sum(flux)
      
      mass_flux = df
      
      do side = We, To
        if (pressure_solution%correct_mass_flux(side)) then
          !part of the flux that goes to that side
          df_side(side) = df * bound_area(side) / correction_area
          
          !the velocity to be added is flux / the area
          df_side(side) = df_side(side) / bound_area(side)          
        end if
      end do
    
      if (pressure_solution%correct_mass_flux(We)) then
        do k = 1, Unz
          do j = 1, Uny
            if (Utype(0,j,k)<=0) U(-2:0,j,k) = U(-2:0,j,k) - df_side(We)
          end do
        end do
      end if
      if (pressure_solution%correct_mass_flux(Ea)) then
        do k = 1, Unz
          do j = 1, Uny
            if (Utype(Prnx,j,k)<=0) U(Prnx:Unx+3,j,k) = U(Prnx:Unx+3,j,k) + df_side(Ea)
          end do
        end do
      end if
      if (pressure_solution%correct_mass_flux(So)) then
          do k = 1, Vnz
            do i = 1, Vnx
              if (Vtype(i,0,k)<=0) V(i,-2:0,k) = V(i,-2:0,k) - df_side(So)
            end do
          end do
      end if
      if (pressure_solution%correct_mass_flux(No)) then
          do k = 1, Vnz
            do i = 1, Vnx
              if (Vtype(i,Prny,k)<=0) V(i,Prny:Vny+3,k) = V(i,Prny:Vny+3,k) + df_side(No)
            end do
          end do
      end if
      if (pressure_solution%correct_mass_flux(Bo)) then
          do j = 1, Wny
            do i = 1, Wnx
              if (Wtype(i,j,0)<=0) W(i,j,-2:0) = W(i,j,-2:0) - df_side(Bo)
            end do
          end do
      end if
      if (pressure_solution%correct_mass_flux(To)) then
          do j = 1, Wny
            do i = 1, Wnx
              if (Wtype(i,j,Prnz)<=0) W(i,j,Prnz:Wnz+3) = W(i,j,Prnz:Wnz+3) + df_side(To)
            end do
          end do
      end if

    end if !check mass flux    

    if (discretization_order == 4) then
    
      !$omp parallel private(i,j,k)
      !$omp do
      do k = 1, Prnz            !divergence of U -> RHS
        do j = 1, Prny
          do i = 1, Prnx
             RHS(i,j,k) = ( C1*(U(i,j,k)-U(i-1,j,k)) - C3*(U(i+1,j,k)-U(i-2,j,k)) ) / dxmin &
                        + ( C1*(V(i,j,k)-V(i,j-1,k)) - C3*(V(i,j+1,k)-V(i,j-2,k)) ) / dymin &
                        + ( C1*(W(i,j,k)-W(i,j,k-1)) - C3*(W(i,j,k+1)-W(i,j,k-2)) ) / dzmin
          end do
        end do
      end do
      !$omp end do
      
      if (allocated(Q)) then
        !$omp do
        do k = 1, Prnz
          do j = 1, Prny
            do i = 1, Prnx
               RHS(i,j,k) = RHS(i,j,k) - Q(i,j,k)
            end do
          end do
        end do
        !$omp end do
      end if
      
      !$omp do
      do k = 1, Prnz
        do j = 1, Prny
          do i = 1, Prnx
             RHS(i,j,k) = RHS(i,j,k) * dt2_rec
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
      
    else
    
      if (gridtype==GRID_VARIABLE_Z) then
      
        if (allocated(Q)) then
          !$omp parallel do private(i,j,k)
          do k = 1, Prnz            !divergence of U -> RHS
            do j = 1, Prny
              do i = 1, Prnx
                RHS(i,j,k) = (U(i,j,k) - U(i-1,j,k)) / dxmin &
                            + (V(i,j,k) - V(i,j-1,k)) / dymin &
                            + (W(i,j,k) - W(i,j,k-1)) / dzPr(k) &
                            - Q(i,j,k)
                RHS(i,j,k) = RHS(i,j,k) * dt2_rec
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(i,j,k)
          do k = 1, Prnz            !divergence of U -> RHS
            do j = 1, Prny
              do i = 1, Prnx
                RHS(i,j,k) = (U(i,j,k) - U(i-1,j,k)) / dxmin &
                            + (V(i,j,k) - V(i,j-1,k)) / dymin &
                            + (W(i,j,k) - W(i,j,k-1)) / dzPr(k)
                RHS(i,j,k) = RHS(i,j,k) * dt2_rec
              end do
            end do
          end do
          !$omp end parallel do
        end if

      else
      
        if (allocated(Q)) then
          !$omp parallel do private(i,j,k)
          do k = 1, Prnz            !divergence of U -> RHS
            do j = 1, Prny
              do i = 1, Prnx
                RHS(i,j,k) = (U(i,j,k) - U(i-1,j,k)) / dxmin &
                            + (V(i,j,k) - V(i,j-1,k)) / dymin &
                            + (W(i,j,k) - W(i,j,k-1)) / dzmin &
                            - Q(i,j,k)
                RHS(i,j,k) = RHS(i,j,k) * dt2_rec
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(i,j,k)
          do k = 1, Prnz            !divergence of U -> RHS
            do j = 1, Prny
              do i = 1, Prnx
                RHS(i,j,k) = (U(i,j,k) - U(i-1,j,k)) / dxmin &
                            + (V(i,j,k) - V(i,j-1,k)) / dymin &
                            + (W(i,j,k) - W(i,j,k-1)) / dzmin
                RHS(i,j,k) = RHS(i,j,k) * dt2_rec
              end do
            end do
          end do
          !$omp end parallel do
        end if
        
      end if
      
    end if

  end subroutine PrePoisson




  subroutine PostPoisson(U,V,W,Pr,Q,Phi,dt2,dt3)
#ifdef PAR
    use custom_par, only: kim, nzims, &
                          par_co_max, par_broadcast_from_last_z, par_co_sum_plane_xy
    use exchange_par, only: par_exchange_UVW
#endif
    real(knd), contiguous, intent(inout) :: U(-2:,-2:,-2:)
    real(knd), contiguous, intent(inout) :: V(-2:,-2:,-2:)
    real(knd), contiguous, intent(inout) :: W(-2:,-2:,-2:)
    real(knd), contiguous, intent(inout) :: Pr(-1:,-1:,-1:)
    real(knd), allocatable, intent(in) :: Q(:,:,:)
    real(knd), intent(inout) :: Phi(-1:,-1:,-1:)
    real(knd), intent(in)    :: dt2,dt3
    real(knd) :: Phi_ref,Au,Av,Aw,dxmin2,dymin2,dzmin2,S,p
    integer   :: i,j,k
    real(knd), parameter :: C1 = 9._knd / 8, C3 = 1._knd / (8*3)


    Au = dt2/dxmin
    Av = dt2/dymin
    Aw = dt2/dzmin

    dxmin2 = dxmin**2
    dymin2 = dymin**2
    dzmin2 = dzmin**2

    !$omp parallel private (i,j,k)
    if (discretization_order == 4) then
      !$omp do
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            U(i,j,k) = U(i,j,k) - Au * ( C1 * (Phi(i+1,j,k) - Phi(i  ,j,k)) - &
                                         C3 * (Phi(i+2,j,k) - Phi(i-1,j,k)) )
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            V(i,j,k) = V(i,j,k) - Av * ( C1 * (Phi(i,j+1,k) - Phi(i,j  ,k)) - &
                                         C3 * (Phi(i,j+2,k) - Phi(i,j-1,k)) )
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Wnz
        do j = 1, Wny
          do i = 1, Wnx
            W(i,j,k) = W(i,j,k) - Aw * ( C1 * (Phi(i,j,k+1) - Phi(i,j,k  )) - &
                                         C3 * (Phi(i,j,k+2) - Phi(i,j,k-1)) )
          end do
        end do
      end do
      !$omp end do nowait

    else
      
      !$omp do
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            U(i,j,k) = U(i,j,k) - Au * (Phi(i+1,j,k) - Phi(i,j,k))
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            V(i,j,k) = V(i,j,k) - Av * (Phi(i,j+1,k) - Phi(i,j,k))
          end do
        end do
      end do
      !$omp end do nowait
      if (gridtype==GRID_VARIABLE_Z) then
        !$omp do
        do k = 1, Wnz
          do j = 1, Wny
            do i = 1, Wnx
              W(i,j,k) = W(i,j,k) - dt2 * (Phi(i,j,k+1) - Phi(i,j,k)) / dzW(k)
            end do
          end do
        end do
        !$omp end do nowait
      else
        !$omp do
        do k = 1, Wnz
          do j = 1, Wny
            do i = 1, Wnx
              W(i,j,k) = W(i,j,k) - Aw * (Phi(i,j,k+1) - Phi(i,j,k))
            end do
          end do
        end do
        !$omp end do nowait
      end if
    end if

    if (pressure_solution%projection_method==PROJECTION_METHOD_PRESSURE) then
        !$omp do
        do k = 1, Prnz
          do j = 1, Prny
            do i = 1, Prnx
              Pr(i,j,k) = Phi(i,j,k)
            end do
          end do
        end do
        !$omp end do
    else
      if (explicit_diffusion) then
        !$omp do
        do k = 1, Prnz
          do j = 1, Prny
            do i = 1, Prnx
              Pr(i,j,k) = Pr(i,j,k) + Phi(i,j,k)
            end do
          end do
        end do
        !$omp end do
      else
        !$omp do
        do k = 1, Prnz
          do j = 1, Prny
            do i = 1, Prnx
              Pr(i,j,k) = Pr(i,j,k) + Phi(i,j,k) - &
                          dt3 * Viscosity(i,j,k) * (((Phi(i+1,j,k)-Phi(i,j,k)) - &
                                                    (Phi(i,j,k)-Phi(i-1,j,k)))/dxmin2 + &
                                                    ((Phi(i,j+1,k)-Phi(i,j,k)) - &
                                                    (Phi(i,j,k)-Phi(i,j-1,k)))/dymin2 + &
                                                    ((Phi(i,j,k+1)-Phi(i,j,k)) - &
                                                    (Phi(i,j,k)-Phi(i,j,k-1)))/dzmin2)
            end do
          end do
        end do
        !$omp end do
      end if
    end if

    Phi_ref = 0
#ifdef PAR
    !images in top plane compute the reference pressure
    if (kim==nzims) then
      !$omp do reduction(+:Phi_ref)
      do j = 1, Prny
        do i = 1, Prnx
          Phi_ref = Phi_ref + Pr(i,j,Prnz)
        end do
      end do
      !$omp end do

      !$omp single
      Phi_ref = par_co_sum_plane_xy(Phi_ref)
      Phi_ref = Phi_ref / (gPrnx * gPrny)
      !$omp end single
    end if

    !all kim==nzims broadcast to images with smaller kim
    !$omp single
    call par_broadcast_from_last_z(Phi_ref)
    !$omp end single

#else

    !$omp do reduction(+:Phi_ref)
    do j = 1, Prny
      do i = 1, Prnx
        Phi_ref = Phi_ref + Pr(i,j,Prnz)
      end do
    end do
    !$omp end do

    !$omp single
    Phi_ref = Phi_ref / (Prnx*Prny)
    !$omp end single

#endif

    !$omp do reduction(+:Phi_ref)
    do k = 1, Prnz
      do j = 1, Prny
        do i = 1, Prnx
          Pr(i,j,k) = Pr(i,j,k) - (Phi_ref)
        end do
      end do
    end do
    !$omp end do   

    !$omp end parallel

#ifdef PAR    
    call par_exchange_UVW(U, V, W)
#endif

    call Bound_Pr(Pr)

    if (pressure_solution%check_divergence) then
      call BoundUVW(U,V,W)
      
      S = 0
      
      if (discretization_order == 4) then
      
        if (allocated(Q)) then
          !$omp parallel do private(i,j,k) reduction(max:S)
          do k = 1, Prnz            !divergence of U -> RHS
            do j = 1, Prny
              do i = 1, Prnx
                 p = ( C1*(U(i,j,k)-U(i-1,j,k)) - C3*(U(i+1,j,k)-U(i-2,j,k)) ) / dxmin &
                   + ( C1*(V(i,j,k)-V(i,j-1,k)) - C3*(V(i,j+1,k)-V(i,j-2,k)) ) / dymin &
                   + ( C1*(W(i,j,k)-W(i,j,k-1)) - C3*(W(i,j,k+1)-W(i,j,k-2)) ) / dzmin
                 S = max(S, abs(p - Q(i,j,k)))              
              end do
            end do
          end do
          !$omp end parallel do
        else
          !$omp parallel do private(i,j,k) reduction(max:S)
          do k = 1, Prnz            !divergence of U -> RHS
            do j = 1, Prny
              do i = 1, Prnx
                 p = ( C1*(U(i,j,k)-U(i-1,j,k)) - C3*(U(i+1,j,k)-U(i-2,j,k)) ) / dxmin &
                   + ( C1*(V(i,j,k)-V(i,j-1,k)) - C3*(V(i,j+1,k)-V(i,j-2,k)) ) / dymin &
                   + ( C1*(W(i,j,k)-W(i,j,k-1)) - C3*(W(i,j,k+1)-W(i,j,k-2)) ) / dzmin
                 S = max(S, abs(p))
              end do
            end do
          end do
          !$omp end parallel do
        end if
        
      else
        if (gridtype==GRID_VARIABLE_Z) then
          if (allocated(Q)) then
            !$omp parallel do private(i,j,k,p) reduction(max:S)
            do k = 1, Prnz            !divergence of U -> RHS
              do j = 1, Prny
                do i = 1, Prnx
                  p =   (U(i,j,k) - U(i-1,j,k)) / (dxmin) &
                      + (V(i,j,k) - V(i,j-1,k)) / (dymin) &
                      + (W(i,j,k) - W(i,j,k-1)) / (dzPr(k)) &
                      - Q(i,j,k)
                  S = max(S,abs(p))
                end do
              end do
            end do
            !$omp end parallel do
          else
            !$omp parallel do private(i,j,k,p) reduction(max:S)
            do k = 1, Prnz            !divergence of U -> RHS
              do j = 1, Prny
                do i = 1, Prnx
                  p =   (U(i,j,k) - U(i-1,j,k)) / (dxmin) &
                      + (V(i,j,k) - V(i,j-1,k)) / (dymin) &
                      + (W(i,j,k) - W(i,j,k-1)) / (dzPr(k))

                  S = max(S,abs(p))
                end do
              end do
            end do
            !$omp end parallel do
          end if
        else
          if (allocated(Q)) then
            !$omp parallel do private(i,j,k,p) reduction(max:S)
            do k = 1, Prnz            !divergence of U -> RHS
              do j = 1, Prny
                do i = 1, Prnx
                  p =   (U(i,j,k) - U(i-1,j,k)) / (dxmin) &
                      + (V(i,j,k) - V(i,j-1,k)) / (dymin) &
                      + (W(i,j,k) - W(i,j,k-1)) / (dzmin) &
                      - Q(i,j,k)
                  S = max(S,abs(p))
                end do
              end do
            end do
            !$omp end parallel do
          else
            !$omp parallel do private(i,j,k,p) reduction(max:S)
            do k = 1, Prnz            !divergence of U -> RHS
              do j = 1, Prny
                do i = 1, Prnx
                  p =   (U(i,j,k) - U(i-1,j,k)) / (dxmin) &
                      + (V(i,j,k) - V(i,j-1,k)) / (dymin) &
                      + (W(i,j,k) - W(i,j,k-1)) / (dzmin)

                  S = max(S,abs(p))
                end do
              end do
            end do
            !$omp end parallel do
          end if
        end if
        
      end if
      
#ifdef PAR
      S = par_co_max(S)
#endif
      
      if (master) write(*,*) "max divergence:", S
    end if

  end subroutine PostPoisson


  subroutine InitHydrostaticPressure(Pr, Temperature, Moisture)
    use PhysicalProperties
    real(knd), contiguous, intent(out) :: Pr(-1:,-1:,-1:)
    real(knd), contiguous, intent(in) :: Temperature(-1:,-1:,-1:), Moisture(-1:,-1:,-1:)
    real(knd) :: t_virt, t_virt_prev, p
    integer :: i, j, k

    interface theta_v
      procedure :: theta_v_ijk
      procedure :: theta_v_Tq
    end interface

#ifdef PAR
    call InitHydrostaticPressure_par(Pr, Temperature, Moisture)
    return
#endif
    
    ! in Pascals
    ! the pressure at the ground is used in configurations
    p = pressure_solution%bottom_pressure

    do k = 1, Prnz
      t_virt = sum(TempIn(1:Prny,k))/Prny

      if (enable_moisture) t_virt = theta_v(t_virt, sum(MoistIn(1:Prny,k))/Prny)

      p = p - rho_air_ref * grav_acc*dzPr(k) * &
              (1 - &
                      ( temperature_ref - t_virt ) &
                      / temperature_ref &
              )
    end do
    !this will be used as a fixed reference in the simulation
    pressure_solution%top_pressure = p
    if (master) write(*,*) "top_pressure:", pressure_solution%top_pressure

    ! reference_pressure_z(k) does not contain any thermal stratification.
    ! Its effect is explicitly simulated and contained in the Pr(i,j,k) field.
    allocate(reference_pressure_z(1:Prnz))
    
    reference_pressure_z(Prnz) = pressure_solution%top_pressure + rho_air_ref * grav_acc*(zW(Prnz)-zPr(Prnz))

    
    
    do k = Prnz-1, 1, -1
      reference_pressure_z(k) = reference_pressure_z(k+1) + rho_air_ref * grav_acc*dzW(k)
    end do



    ! Pr(i,j,k) is the physical pressure in Pascals divided by rho_air_ref
    do j = 1, Vny+1
      do i = 1, Unx+1
        t_virt_prev = theta_v(i,j,Prnz)
        Pr(i,j,Prnz) = - grav_acc*(zW(Prnz+1)-zPr(Prnz)) * &
                ( t_virt_prev - temperature_ref ) &
                / temperature_ref
        do k = Prnz-1, 1, -1
          t_virt = theta_v(i,j,k)
          Pr(i,j,k) = Pr(i,j,k+1) - &
                 grav_acc*dzW(k) * &
                ( (t_virt+t_virt_prev)/2._knd - temperature_ref ) &
                / temperature_ref
          t_virt_prev = t_virt
        end do
      end do
    end do


  contains

      pure function theta_v_ijk(i,j,k) result(res)
        real(knd) :: res
        integer, intent(in) :: i,j,k

        if (enable_moisture) then
          res = Temperature(i,j,k) * (1._knd + 0.61_knd * Moisture(i,j,k))
        else
          res = Temperature(i,j,k)
        end if
      end function

      pure function theta_v_tq(T, q) result(res)
        real(knd) :: res
        real(knd), intent(in) :: T, q

        if (enable_moisture) then
          res = T * (1._knd + 0.61_knd * q)
        else
          res = T
        end if
      end function

  end subroutine InitHydrostaticPressure

#ifdef PAR
  subroutine InitHydrostaticPressure_par(Pr, Temperature, Moisture)
    use custom_par
    use PhysicalProperties
    real(knd), contiguous, intent(out) :: Pr(-1:,-1:,-1:)
    real(knd), contiguous, intent(in) :: Temperature(-1:,-1:,-1:), Moisture(-1:,-1:,-1:)
    real(knd) :: t_virt, t_virt_prev, p
    integer :: i, j, k
    integer :: nz

    interface theta_v
      procedure :: theta_v_ijk
      procedure :: theta_v_Tq
    end interface

    
    ! in Pascals
    ! the pressure at the ground is used in configurations
    if (kim==1) then
      p = pressure_solution%bottom_pressure
    else
      call par_recv(p, b_im)
    end if

    do k = 1, Prnz
      t_virt = par_co_sum_plane_xy(sum(TempIn(1:Prny,k)))/gPrny

      if (enable_moisture) t_virt = theta_v(t_virt, par_co_sum_plane_xy(sum(MoistIn(1:Prny,k)))/gPrny)

      p = p - rho_air_ref * grav_acc*dzPr(k) * &
              (1 - &
                      ( temperature_ref - t_virt ) &
                      / temperature_ref &
              )
    end do
    !this will be used as a fixed reference in the simulation

    if (kim == nzims) then
      pressure_solution%top_pressure = p
    else
      call par_send(p, t_im)
    end if
    
    call par_broadcast_from_last_z(pressure_solution%top_pressure)
    
    
    if (master) write(*,*) "top_pressure:", pressure_solution%top_pressure
  
    
    
    ! reference_pressure_z(k) does not contain any thermal stratification.
    ! Its effect is explicitly simulated and contained in the Pr(i,j,k) field.
    allocate(reference_pressure_z(1:Wnz+1))
    

    if (kim == nzims) then    
      reference_pressure_z(Prnz) = pressure_solution%top_pressure + rho_air_ref * grav_acc*(zW(Prnz)-zPr(Prnz))
      nz = Prnz - 1
    else
      call par_recv(reference_pressure_z(Prnz+1), t_im)
      nz = Prnz
    end if
    
    do k = nz, 1, -1
      reference_pressure_z(k) = reference_pressure_z(k+1) + rho_air_ref * grav_acc*dzW(k)
    end do
    
    if (kim>1) call par_send(reference_pressure_z(1), b_im)



    ! Pr(i,j,k) is the physical pressure in Pascals divided by rho_air_ref
    do j = 1, Vny+1
      do i = 1, Unx+1
      
        if (kim == nzims) then
          t_virt_prev = theta_v(i,j,Prnz)
          Pr(i,j,Prnz) = - grav_acc*(zW(Prnz+1)-zPr(Prnz)) * &
                  ( t_virt_prev - temperature_ref ) &
                  / temperature_ref
          nz = Prnz -1
        else
          call par_recv(Pr(i,j,Prnz+1), t_im)
          call par_recv(t_virt_prev, t_im)
          nz = Prnz
        end if
          
        do k = nz, 1, -1
          t_virt = theta_v(i,j,k)
          Pr(i,j,k) = Pr(i,j,k+1) - &
                 grav_acc*dzW(k) * &
                ( (t_virt+t_virt_prev)/2._knd - temperature_ref ) &
                / temperature_ref
          t_virt_prev = t_virt
        end do
        
        if (kim > 1) then
          call par_send(Pr(i,j,1), b_im)
          call par_send(t_virt_prev, b_im)
        end if
        
      end do
    end do


  contains

      pure function theta_v_ijk(i,j,k) result(res)
        real(knd) :: res
        integer, intent(in) :: i,j,k

        if (enable_moisture) then
          res = Temperature(i,j,k) * (1._knd + 0.61_knd * Moisture(i,j,k))
        else
          res = Temperature(i,j,k)
        end if
      end function

      pure function theta_v_tq(T, q) result(res)
        real(knd) :: res
        real(knd), intent(in) :: T, q

        if (enable_moisture) then
          res = T * (1._knd + 0.61_knd * q)
        else
          res = T
        end if
      end function

  end subroutine InitHydrostaticPressure_par
#endif  



end module Pressure
