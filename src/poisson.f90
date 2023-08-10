module PoissonResidues
  use Parameters
  use Boundaries
  use custom_par
    
  implicit none
  
  private
  
  public Poiss_Residue

contains

  subroutine Poiss_Residue(Phi, RHS, res)
    real(knd), intent(in) :: Phi(-1:,-1:,-1:)
    real(knd), intent(in) :: RHS(0:,0:,0:)
    real(knd), intent(out) :: res
    
    call Residue_variable_z(Phi, RHS, 1/dxmin**2, 1/dxmin**2, 1/dymin**2, 1/dymin**2, zPr(0:), zW(-1:), res)
  end subroutine

 
  subroutine Residue_variable_z(Phi,RHS, &
                    Aw,Ae,As,An,z,z_u,   &
                    res)
    real(knd), intent(in) :: Phi(-1:,-1:,-1:)
    real(knd), intent(in) :: RHS(0:,0:,0:)
    real(knd), intent(in) :: Aw,Ae
    real(knd), intent(in) :: As,An
    real(knd), intent(in) :: z(0:), z_u(-1:)
    real(knd), intent(out) :: res
    integer :: i, j, k
    real(knd) :: p, Ap, Ab, At
    
    res = 0
    
    do k = 1, Prnz
      do j = 1, Prny
          do i = 1, Prnx
            Ab = 1 / ((z(k) - z(k-1)) * (z_u(k)-z_u(k-1)))
            At = 1 / ((z(k+1) - z(k)) * (z_u(k)-z_u(k-1)))
            Ap = Aw + Ae + As + An + Ab + At
            
            p = 0
            p = p + Phi(i-1,j,k) * Aw
            p = p + Phi(i+1,j,k) * Ae
            p = p + Phi(i,j-1,k) * As
            p = p + Phi(i,j+1,k) * An
            p = p + Phi(i,j,k-1) * Ab
            p = p + Phi(i,j,k+1) * At
            p = p - RHS(i,j,k)
            p = abs(-p + Ap*Phi(i,j,k))
            res = max(res, abs(p))
          end do
        end do
    end do
#ifdef PAR
    res = par_co_max(res)
#endif
  end subroutine Residue_variable_z
   
end module PoissonResidues


module PoissonSolvers

  use Parameters
  use Multigrid,   only: PoissMG
  use Multigrid2D, only: PoissMG2d

  implicit none


contains

  subroutine Poiss_PoisFFT(Phi,RHS)
#ifdef DPREC
    use PoisFFT, PoisFFT_Solver => PoisFFT_Solver3D_DP
#else
    use PoisFFT, PoisFFT_Solver => PoisFFT_Solver3D_SP
#endif
#ifdef PAR
    use custom_par
#endif

    type(PoisFFT_Solver),save :: Solver
    real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Phi
    real(knd),dimension(0:,0:,0:),intent(in) :: RHS
    logical, save :: called = .false.

    integer(int64), save :: trate
    integer(int64)       :: t1, t2

    integer :: discretization
    
    if (discretization_order == 4) then
      discretization = PoisFFT_FiniteDifference4
    else
      discretization = PoisFFT_FiniteDifference2
    end if

    if (.not.called) then
#ifdef PAR
      Solver =  PoisFFT_Solver([Prnx,Prny,Prnz], &
                               [gxmax-gxmin,gymax-gymin,gzmax-gzmin], &
                               PoissonBtype, &
                               discretization, &
                               gPrns,offsets_to_global,poisfft_comm)
#else
      Solver =  PoisFFT_Solver([Prnx,Prny,Prnz], &
                               [gxmax-gxmin,gymax-gymin,gzmax-gzmin], &
                               PoissonBtype, &
                               discretization)
#endif
      called = .true.

      call system_clock(count_rate=trate)

    end if


    call system_clock(count=t1)


    call Execute(Solver,Phi,RHS)


    call system_clock(count=t2)
    if (master) then
      poisson_solver_time = poisson_solver_time + real(t2-t1,dbl)/real(trate,dbl)
      if (master.and.debugparam > 1) write(*,*) "solver cpu time", real(t2-t1)/real(trate)
    end if

  end subroutine Poiss_PoisFFT





  subroutine Poiss_PoisFFT_variable_z(Phi,RHS)
#ifdef DPREC
    use PoisFFT, PoisFFT_Solver => PoisFFT_Solver3D_nonuniform_z_DP
#else
    use PoisFFT, PoisFFT_Solver => PoisFFT_Solver3D_nonuniform_z_SP
#endif
#ifdef PAR
    use custom_par
#endif

    type(PoisFFT_Solver),save :: Solver
    real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Phi
    real(knd),dimension(0:,0:,0:),intent(in) :: RHS
    logical, save :: called = .false.

    integer(int64), save :: trate
    integer(int64)       :: t1, t2

    integer :: discretization

    if (discretization_order == 4) then
      discretization = PoisFFT_FiniteDifference4
    else
      discretization = PoisFFT_FiniteDifference2
    end if

    if (.not.called) then
#ifdef PAR
      Solver =  PoisFFT_Solver([Prnx,Prny,Prnz], &
                               [gxmax-gxmin,gymax-gymin], gzPr(1:gPrnz), gzW(0:gPrnz),  &
                               PoissonBtype, &
                               discretization, &
                               gPrns,offsets_to_global,poisfft_comm)
#else
      Solver =  PoisFFT_Solver([Prnx,Prny,Prnz], &
                               [gxmax-gxmin,gymax-gymin], zPr(1:Prnz), zW(0:Prnz), &
                               PoissonBtype, &
                               discretization)
#endif
      called = .true.

      call system_clock(count_rate=trate)

    end if


    call system_clock(count=t1)


    call Execute(Solver,Phi,RHS)


    call system_clock(count=t2)
    if (master) then
      poisson_solver_time = poisson_solver_time + real(t2-t1,dbl)/real(trate,dbl)
      if (master.and.debugparam > 1) write(*,*) "solver cpu time", real(t2-t1)/real(trate)
    end if

  end subroutine Poiss_PoisFFT_variable_z





  subroutine PoissSOR(Phi,RHS) 
#ifdef PAR
    use Boundaries
    use custom_par
#endif
    !Solves Poisson equation using Successive over-relaxation

    real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Phi
    real(knd),dimension(0:,0:,0:),intent(in) :: RHS
    integer,save :: called=0
    integer :: i,j,k,l
    real(knd) :: S,P,Ap
    real(knd),dimension(:),allocatable,save :: Aw,Ae,As,An,Ab,At


    write (*,*) "Computing Poisson equation"
    S=0

    if (called==0) then                             !coefficients based on grid spacing computed only once
      allocate(Aw(1:Prnx),Ae(1:Prnx))
      allocate(As(1:Prny),An(1:Prny))
      allocate(Ab(1:Prnz),At(1:Prnz))
      forall(i=1:Prnx)
        Ae(i)=1._knd/(dxU(i)*dxPr(i))
        Aw(i)=1._knd/(dxU(i-1)*dxPr(i))
      end forall
      forall(j=1:Prny)
        An(j)=1._knd/(dyV(j)*dyPr(j))
        As(j)=1._knd/(dyV(j-1)*dyPr(j))
      end forall
      forall(k=1:Prnz)
        At(k)=1._knd/(dzW(k)*dzPr(k))
        Ab(k)=1._knd/(dzW(k-1)*dzPr(k))
      end forall
      called=1
    end if

    l=0
    S=huge(1.0_knd)
     
    do while (l<=maxPoissoniter.and.S>epsPoisson)
      l=l+1
      S=0
#ifdef PAR
      call Bound_Phi(Phi)
#endif
      !$OMP PARALLEL PRIVATE(i,j,k,p) REDUCTION(max:S)
      !$OMP DO
      do k=1,Prnz
        do j=1,Prny
          do i=1+mod(j+k,2),Prnx,2
            p=0
            Ap=0
            if (i > 1.or.(Btype(We)>=BC_MPI_BOUNDS_MIN.and.Btype(We)<=BC_MPI_BOUNDS_MAX)) then
                      p = p + Phi(i-1,j,k)*Aw(i)
                      Ap=Ap+Aw(i)
            else if (Btype(We)==BC_PERIODIC) then
                      p = p + Phi(Prnx,j,k)*Aw(i)
                      Ap=Ap+Aw(i)
            end if
            if (i < Prnx.or.(Btype(Ea)>=BC_MPI_BOUNDS_MIN.and.Btype(Ea)<=BC_MPI_BOUNDS_MAX)) then
                      p = p + Phi(i+1,j,k)*Ae(i)
                      Ap=Ap+Ae(i)
            else if (Btype(Ea)==BC_PERIODIC) then
                      p = p + Phi(1,j,k)*Ae(i)
                      Ap=Ap+Ae(i)
            end if
            if (j > 1.or.(Btype(So)>=BC_MPI_BOUNDS_MIN.and.Btype(So)<=BC_MPI_BOUNDS_MAX)) then
                      p = p + Phi(i,j-1,k)*As(j)
                      Ap=Ap+As(j)
            else if (Btype(So)==BC_PERIODIC) then
                      p = p + Phi(i,Prny,k)*As(j)
                      Ap=Ap+As(j)
            end if
            if (j < Prny.or.(Btype(No)>=BC_MPI_BOUNDS_MIN.and.Btype(No)<=BC_MPI_BOUNDS_MAX)) then
                      p = p + Phi(i,j+1,k)*An(j)
                      Ap=Ap+An(j)
            else if (Btype(No)==BC_PERIODIC) then
                      p = p + Phi(i,1,k)*An(j)
                      Ap=Ap+An(j)
            end if
            if (k > 1.or.(Btype(Bo)>=BC_MPI_BOUNDS_MIN.and.Btype(Bo)<=BC_MPI_BOUNDS_MAX)) then
                      p = p + Phi(i,j,k-1)*Ab(k)
                      Ap=Ap+Ab(k)
            else if (Btype(Bo)==BC_PERIODIC) then
                      p = p + Phi(i,j,Prnz)*Ab(k)
                      Ap=Ap+Ab(k)
            end if
            if (k < Prnz.or.(Btype(To)>=BC_MPI_BOUNDS_MIN.and.Btype(To)<=BC_MPI_BOUNDS_MAX)) then
                      p = p + Phi(i,j,k+1)*At(k)
                      Ap=Ap+At(k)
            else if (Btype(To)==BC_PERIODIC) then
                      p = p + Phi(i,j,1)*At(k)
                      Ap=Ap+At(k)
            end if
            p=p-RHS(i,j,k)
            p=p/Ap
            S=max(abs(p-Phi(i,j,k)),S)
            Phi(i,j,k)=p
          end do
        end do
      end do
      !$OMP ENDDO
      !$OMP DO
      do k=1,Prnz
        do j=1,Prny
          do i=1+mod(j+k+1,2),Prnx,2
            p=0
            Ap=0
            if (i > 1) then
                      p = p + Phi(i-1,j,k)*Aw(i)
                      Ap=Ap+Aw(i)
            else if (Btype(We)==BC_PERIODIC) then
                      p = p + Phi(Prnx,j,k)*Aw(i)
                      Ap=Ap+Aw(i)
            end if
            if (i < Prnx) then
                      p = p + Phi(i+1,j,k)*Ae(i)
                      Ap=Ap+Ae(i)
            else if (Btype(We)==BC_PERIODIC) then
                      p = p + Phi(1,j,k)*Ae(i)
                      Ap=Ap+Ae(i)
            end if
            if (j > 1) then
                      p = p + Phi(i,j-1,k)*As(j)
                      Ap=Ap+As(j)
            else if (Btype(No)==BC_PERIODIC) then
                      p = p + Phi(i,Prny,k)*As(j)
                      Ap=Ap+As(j)
            end if
            if (j < Prny) then
                      p = p + Phi(i,j+1,k)*An(j)
                      Ap=Ap+An(j)
            else if (Btype(No)==BC_PERIODIC) then
                      p = p + Phi(i,1,k)*An(j)
                      Ap=Ap+An(j)
            end if
            if (k > 1) then
                      p = p + Phi(i,j,k-1)*Ab(k)
                      Ap=Ap+Ab(k)
            else if (Btype(To)==BC_PERIODIC) then
                      p = p + Phi(i,j,Prnz)*Ab(k)
                      Ap=Ap+Ab(k)
            end if
            if (k < Prnz) then
                      p = p + Phi(i,j,k+1)*At(k)
                      Ap=Ap+At(k)
            else if (Btype(To)==BC_PERIODIC) then
                      p = p + Phi(i,j,1)*At(k)
                      Ap=Ap+At(k)
            end if
            p=p-RHS(i,j,k)
            p=p/Ap
            S=max(abs(p-Phi(i,j,k)),S)
            Phi(i,j,k)=p
          end do
        end do
      end do
      !$OMP ENDDO
      !$OMP ENDPARALLEL   
      p=abs(maxval(Phi(1:Prnx,1:Prny,1:Prnz)))
      if (p>0) S=S/p
      
#ifdef PAR
      S = par_co_max(S)
#endif

      if (MOD(l,10)==0)  write (*,*) "   Poisson iter: ",l,S
    end do

  end subroutine PoissSOR


end module PoissonSolvers
