module Multigrid2d

  use Parameters
  use Stop_procedures

  implicit none

  private

  public POISSMG2d,SetMGParams2d

  type TMGArr
    real(knd),allocatable,dimension(:,:) :: Arr
  end type

  type TCoefs
    real(knd) :: dx,dy,dz,Aw,Ae,As,An,Ab,At
    integer :: nx,nz
  end type


  type(TMGArr),allocatable,dimension(:) :: PhiMG,RHSMG,ResMG
  type(TCoefs),allocatable,dimension(:) :: CoefMg

  integer :: bnx,bnz !Base cube dimensions for Multigrid, Prnx = bnx*2**level
  integer :: LMG !depth of multigrid
  integer :: minmglevel !innermost MG level
  integer :: mgncgc !type of cycling 1..V cycle, 2..W cycle
  integer :: mgnpre
  integer :: mgnpost
  integer :: mgmaxinnerGSiter
  real(knd) :: mgepsinnerGS


contains

  subroutine SetMGParams2d(llmg, lminmglevel, lbnx,  lbnz,&
                           lmgncgc, lmgnpre, lmgnpost, lmgmaxinnerGSiter, lmgepsinnerGS)

    integer,intent(in) ::   llmg, lminmglevel, lbnx, lbnz, lmgncgc, lmgnpre, lmgnpost, lmgmaxinnerGSiter
    real(knd),intent(in) :: lmgepsinnerGS

    lmg = llmg
    minmglevel = lminmglevel
    bnx = lbnx
    bnz = lbnz
    mgncgc = lmgncgc
    mgnpre = lmgnpre
    mgnpost = lmgnpost
    mgmaxinnerGSiter = lmgmaxinnerGSiter
    mgepsinnerGS = lmgepsinnerGS
  end subroutine SetMGParams2d


  pure subroutine Prolongate(AFine,ACoarse,level)
    integer,intent(in) :: level
    real(knd),dimension(-1:,-1:),intent(in) :: ACoarse
    real(knd),dimension(-1:,-1:),intent(inout) :: AFine
    integer :: i,k,nx,nz

    nx = bnx*2**level !level means from which grid we interpolate
    nz = bnz*2**level

    do k = 0,nz
      do i = 0,nx
                      AFine(2*i,2*k) = AFine(2*i,2*k)+ACoarse(i,k)
                      AFine(2*i+1,2*k) = AFine(2*i+1,2*k)+0.5*ACoarse(i,k)
                      AFine(2*i-1,2*k) = AFine(2*i-1,2*k)+0.5*ACoarse(i,k)
                      AFine(2*i,2*k+1) = AFine(2*i,2*k+1)+0.5*ACoarse(i,k)
                      AFine(2*i,2*k-1) = AFine(2*i,2*k-1)+0.5*ACoarse(i,k)
                      AFine(2*i+1,2*k+1) = AFine(2*i+1,2*k+1)+0.25*ACoarse(i,k)
                      AFine(2*i-1,2*k-1) = AFine(2*i-1,2*k-1)+0.25*ACoarse(i,k)
                      AFine(2*i-1,2*k+1) = AFine(2*i-1,2*k+1)+0.25*ACoarse(i,k)
                      AFine(2*i+1,2*k-1) = AFine(2*i+1,2*k-1)+0.25*ACoarse(i,k)
     end do
    end do

    call Bound_Phi_MG(Afine,2*nx,2*nz)
  end subroutine Prolongate



  subroutine Restrict(ACoarse,AFine,level)
    integer,intent(in) :: level
    real(knd),dimension(-1:,-1:),intent(out) :: ACoarse
    real(knd),dimension(-1:,-1:),intent(inout) :: AFine
    real(knd) :: p,S,w
    integer :: i,k,ii,kk,nx,nz
    real(knd) :: weight(-1:1,-1:1)

    nx = bnx*2**level !level means on which grid we restrict
    nz = bnz*2**level

    call Bound_Phi_MG(Afine,2*nx,2*nz)

    weight = 1
    weight(-1,:) = weight(-1,:) / 2._knd
    weight( 1,:) = weight( 1,:) / 2._knd
    weight(:,-1) = weight(:,-1) / 2._knd
    weight(:, 1) = weight(:, 1) / 2._knd

    !$omp parallel do private(i,k,ii,kk,w,S,p)
    do k = 0,nz
      do i = 0,nx

        p = 0
        S = 0

        do kk = max(2*k-1, 0), min(2*k+1, 2*nz)
          do ii = max(2*i-1, 0), min(2*i+1, 2*nx)
            w = weight(ii-2*i, kk-2*k)
            S = S + w * AFine(ii,kk)
            p = p + w
          end do
        end do

        ACoarse(i,k) = S/p

      end do
    end do
    !$omp end parallel do
    if (Btype(Ea)==BC_PERIODIC) ACoarse(0,:) = ACoarse(nx,:)
    if (Btype(To)==BC_PERIODIC) ACoarse(:,0) = ACoarse(:,nz)
  end subroutine Restrict







  pure subroutine BOUND_Phi_MG(Phi,nx,nz)
    real(knd),intent(inout) :: Phi(-1:,-1:)
    integer,intent(in) :: nx,nz
    integer :: i,k


     if (Btype(We)==BC_PERIODIC) then
      do k = 0,nz
        Phi(-1,k) = Phi(nx-1,k)
      end do
    else
      do k = 0,nz
        Phi(-1,k) = Phi(0,k)
      end do
     end if

     if (Btype(Ea)==BC_PERIODIC) then
      do k = 0,nz
        Phi(nx+1,k) = Phi(1,k)
      end do
     else
      do k = 0,nz
        Phi(nx+1,k) = Phi(nx,k)
      end do
     end if

     if (Btype(Bo)==BC_PERIODIC) then
       do i=-1,nx+1                      !Periodic BC
        Phi(i,-1) = Phi(i,nz-1)
       end do
     else
       do i=-1,nx+1                      !Other BCs
        Phi(i,-1) = Phi(i,0)
       end do
     end if

     if (Btype(To)==BC_PERIODIC) then
       do i=-1,nx+1                      !Periodic BC
        Phi(i,nz+1) = Phi(i,1)
       end do
     else
       do i=-1,nx+1                      !Other BCs
        Phi(i,nz+1) = Phi(i,nz)
       end do
     end if
  end subroutine BOUND_Phi_MG



  pure function ind(level,nulx,nulz,i,k)
    integer :: ind
    integer,intent(in) :: level,i,k,nulx,nulz
    ind = (1-nulx)+i+(k-nulz)*(CoefMG(level)%nx+(1-nulx))
  end function ind




  subroutine MG_GE(level)
    use Lapack
    integer,intent(in) :: level
    integer :: i,k,l,info
    integer,allocatable,dimension(:),save :: ige,ipivot,work2
    real(knd),allocatable,dimension(:),save :: xge,bge,work,R,C,ferr,berr
    real(knd),allocatable,dimension(:,:),save :: age,af
    real(knd) :: Ap,rcond
    integer,save :: nx,nz,nulx,nulz,nxyz,called = 0
    character(1),save :: equed

    if (called==0) then
      nx = CoefMG(level)%nx
      nz = CoefMG(level)%nz

      if (Btype(Ea)==BC_PERIODIC) then
        nulx = 1
      else
        nulx = 0
      end if
      if (Btype(To)==BC_PERIODIC) then
        nulz = 1
      else
        nulz = 0
      end if
      nxyz = (nx+(1-nulx))*(nz+(1-nulz))
      allocate(xge(nxyz))
      allocate(bge(nxyz))
      allocate(ige(nxyz))
      allocate(ipivot(nxyz))
      allocate(age(nxyz,nxyz))
      allocate(af(nxyz,nxyz))
      allocate(work(4*nxyz))
      allocate(work2(nxyz))
      allocate(R(nxyz))
      allocate(C(nxyz))
      allocate(ferr(nxyz))
      allocate(berr(nxyz))

      called = 1


      age = 0
      bge = 0
      xge = 0
      ige = 0
      do k = nulz,nz
        do i = nulx,nx
          l = ind(level,nulx,nulz,i,k)
          Ap = 0
          if (i>nulx) then
                    age(l,ind(level,nulx,nulz,i-1,k)) = CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          else if (Btype(We)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,nx,k)) = CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          end if
          if (i<nx) then
                    age(l,ind(level,nulx,nulz,i+1,k)) = CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          else if (Btype(Ea)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,nulx,k)) = CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          end if
          if (k>nulz) then
                    age(l,ind(level,nulx,nulz,i,k-1)) = CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          else if (Btype(Bo)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,i,nz)) = CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          end if
          if (k<nz) then
                    age(l,ind(level,nulx,nulz,i,k+1)) = CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          else if (Btype(To)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,i,nulz)) = CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          end if

          age(l,l)=-Ap
          bge(l) = RHSMG(level)%Arr(i,k)
        end do
      end do

      l = ind(level,nulx,nulz,nx,nz)
      age(1,:) = 0
      age(1,1) = age(2,2)
      bge(1) = 0

      call gesvx("E","N",nxyz,1,age,nxyz,af,nxyz,ipivot,EQUED,R,C,bge,nxyz,xge,nxyz,rcond, &
                ferr,berr,work,work2,info)

      write (*,*) "RCOND",rcond
      if (info/=0) then
      call error_stop
      write (*,*) "info",info
      end if

    else

      do k = nulz,nz
        do i = nulx,nx
          l = ind(level,nulx,nulz,i,k)
          bge(l) = RHSMG(level)%Arr(i,k)
        end do
      end do
      bge(1) = 0

      call gesvx("F","N",nxyz,1,age,nxyz,af,nxyz,ipivot,EQUED,R,C,bge,nxyz,xge,nxyz,rcond, &
                ferr,berr,work,work2,info)

      if (info/=0) then
        write (*,*) "info",info
        call error_stop
      end if
    end if


    do k = nulz,nz
      do i = nulx,nx
        l = ind(level,nulx,nulz,i,k)
        PhiMG(level)%Arr(i,k) = xge(l)
      end do
    end do

    if (Btype(Ea)==BC_PERIODIC) then
      PhiMG(level)%Arr(0,:) = PhiMG(level)%Arr(nx,:)
    end if
    if (Btype(To)==BC_PERIODIC) then
      PhiMG(level)%Arr(:,0) = PhiMG(level)%Arr(:,nz)
    end if

  end subroutine MG_GE


  subroutine MG_LU(level)
    use Lapack
    integer,intent(in) :: level
    integer :: i,k,l,info
    integer,allocatable,dimension(:),save :: ige,ipivot,work2
    real(knd),allocatable,dimension(:),save :: xge,bge,work,R,C,ferr,berr
    real(knd),allocatable,dimension(:,:),save :: age,af
    real(knd) :: Ap
    integer,save :: nx,nz,nulx,nulz,nxyz,called = 0

    if (called==0) then
      nx = CoefMG(level)%nx
      nz = CoefMG(level)%nz

      if (Btype(Ea)==BC_PERIODIC) then
        nulx = 1
      else
        nulx = 0
      end if
      if (Btype(To)==BC_PERIODIC) then
        nulz = 1
      else
        nulz = 0
      end if
      nxyz = (nx+(1-nulx))*(nz+(1-nulz))
      allocate(xge(nxyz))
      allocate(bge(nxyz))
      allocate(ige(nxyz))
      allocate(ipivot(nxyz))
      allocate(age(nxyz,nxyz))
      allocate(af(nxyz,nxyz))
      allocate(work(4*nxyz))
      allocate(work2(nxyz))
      allocate(R(nxyz))
      allocate(C(nxyz))
      allocate(ferr(nxyz))
      allocate(berr(nxyz))

      called = 1


      age = 0
      xge = 0
      ige = 0
      do k = nulz,nz
        do i = nulx,nx
          l = ind(level,nulx,nulz,i,k)
          Ap = 0
          if (i>nulx) then
                    age(l,ind(level,nulx,nulz,i-1,k)) = CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          else if (Btype(We)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,nx,k)) = CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          end if
          if (i<nx) then
                    age(l,ind(level,nulx,nulz,i+1,k)) = CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          else if (Btype(Ea)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,nulx,k)) = CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          end if
          if (k>nulz) then
                    age(l,ind(level,nulx,nulz,i,k-1)) = CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          else if (Btype(Bo)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,i,nz)) = CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          end if
          if (k<nz) then
                    age(l,ind(level,nulx,nulz,i,k+1)) = CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          else if (Btype(To)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,i,nulz)) = CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          end if

          age(l,l)=-Ap
        end do
      end do

      l = ind(level,nulx,nulz,nx,nz)
      age(1,:) = 0
      age(1,1) = age(2,2)

      call  getrf(nxyz,nxyz,age,nxyz,ipivot,info)

      if (info/=0) then
        call error_stop
        write (*,*) "info",info
      end if
    end if


    do k = nulz,nz
      do i = nulx,nx
        l = ind(level,nulx,nulz,i,k)
        bge(l) = RHSMG(level)%Arr(i,k)
      end do
    end do
    bge(1) = 0

    call getrs("N",nxyz,1,age,nxyz,ipivot,bge,nxyz,info)

    if (info/=0) then
      write (*,*) "info",info
      call error_stop
    end if


    do k = nulz,nz
      do i = nulx,nx
        l = ind(level,nulx,nulz,i,k)
        PhiMG(level)%Arr(i,k) = bge(l)
      end do
    end do

    if (Btype(Ea)==BC_PERIODIC) then
      PhiMG(level)%Arr(0,:) = PhiMG(level)%Arr(nx,:)
    end if
    if (Btype(To)==BC_PERIODIC) then
      PhiMG(level)%Arr(:,0) = PhiMG(level)%Arr(:,nz)
    end if

  end subroutine MG_LU





  subroutine MG_INV(level)
    use Lapack
    integer,intent(in) :: level
    integer :: i,k,l,info,ldwork
    integer,allocatable,dimension(:) :: ipivot
    real(knd),allocatable,dimension(:),save :: xge,bge,work
    real(knd),allocatable,dimension(:,:),save :: age
    real(knd) :: Ap
    integer,save :: nx,nz,nulx,nulz,nxyz,called = 0

    if (called==0) then
      nx = CoefMG(level)%nx
      nz = CoefMG(level)%nz

      if (Btype(Ea)==BC_PERIODIC) then
        nulx = 1
      else
        nulx = 0
      end if
      if (Btype(To)==BC_PERIODIC) then
        nulz = 1
      else
        nulz = 0
      end if
      nxyz = (nx+(1-nulx))*(nz+(1-nulz))
      allocate(xge(nxyz))
      allocate(bge(nxyz))
      allocate(ipivot(nxyz))
      allocate(age(nxyz,nxyz))


      called = 1


      age = 0
      xge = 0
      do k = nulz,nz
        do i = nulx,nx
          l = ind(level,nulx,nulz,i,k)
          Ap = 0
          if (i>nulx) then
                    age(l,ind(level,nulx,nulz,i-1,k)) = CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          else if (Btype(We)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,nx,k)) = CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          end if
          if (i<nx) then
                    age(l,ind(level,nulx,nulz,i+1,k)) = CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          else if (Btype(Ea)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,nulx,k)) = CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          end if
          if (k>nulz) then
                    age(l,ind(level,nulx,nulz,i,k-1)) = CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          else if (Btype(Bo)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,i,nz)) = CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          end if
          if (k<nz) then
                    age(l,ind(level,nulx,nulz,i,k+1)) = CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          else if (Btype(To)==BC_PERIODIC) then
                    age(l,ind(level,nulx,nulz,i,nulz)) = CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          end if

          age(l,l)=-Ap
        end do
      end do

      l = ind(level,nulx,nulz,nx,nz)
      age(1,:) = 0
      age(1,1) = age(2,2)

      call  getrf(nxyz,nxyz,age,nxyz,ipivot,info)

      if (info/=0) then
        write (*,*) "info",info
        call error_stop
      end if

      allocate(work(1))

      call getri(nxyz,age,nxyz,ipivot,work,-1,info)

      if (info/=0) then
        write (*,*) "info",info
        call error_stop
      end if

      ldwork = nint(work(1))

      deallocate(work)
      allocate(work(ldwork))

      call getri(nxyz,age,nxyz,ipivot,work,ldwork,info)

      if (info/=0) then
        write (*,*) "info",info
        call error_stop
      end if

    end if


    do k = nulz,nz
      do i = nulx,nx
        l = ind(level,nulx,nulz,i,k)
        bge(l) = RHSMG(level)%Arr(i,k)
      end do
    end do
    bge(1) = 0


    xge = 0
    do k = 1,nxyz !rows of X{j}
      do i = 1,nxyz !columns of A(i,j) if (age(i,j)>100*tiny(1._knd))
        xge(i) = xge(i)+age(i,k)*bge(k)
      end do
    end do

    do k = nulz,nz
      do i = nulx,nx
        l = ind(level,nulx,nulz,i,k)
        PhiMG(level)%Arr(i,k) = xge(l)
      end do
    end do

    if (Btype(Ea)==BC_PERIODIC) then
      PhiMG(level)%Arr(0,:) = PhiMG(level)%Arr(nx,:)
    end if
    if (Btype(To)==BC_PERIODIC) then
      PhiMG(level)%Arr(:,0) = PhiMG(level)%Arr(:,nz)
    end if

  end subroutine MG_INV









  subroutine MG_GS(level,niter)
    integer,intent(in) :: level,niter
    integer :: i,k,l
    real(knd) :: p,Ap

    do l = 1,niter
      !$omp parallel private(i,k,p,Ap)
      !$omp do
      do k = 0,CoefMG(level)%nz
        do i = 0+mod(k,2),CoefMG(level)%nx,2
          p = 0
          Ap = 0
          if (i>0) then
                    p = p+PhiMG(level)%Arr(i-1,k)*CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          else if (Btype(We)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(CoefMG(level)%nx-1,k)*CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          end if
          if (i<CoefMG(level)%nx) then
                    p = p+PhiMG(level)%Arr(i+1,k)*CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          else if (Btype(Ea)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(1,k)*CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          end if
          if (k>0) then
                    p = p+PhiMG(level)%Arr(i,k-1)*CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          else if (Btype(Bo)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(i,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          end if
          if (k<CoefMG(level)%nz) then
                    p = p+PhiMG(level)%Arr(i,k+1)*CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          else if (Btype(To)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(i,1)*CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          end if
          p = p-RHSMG(level)%Arr(i,k)

          p = p/Ap
          PhiMG(level)%Arr(i,k) = p
        end do
      end do
      !$omp end do
      !$omp do
      do k = 0,CoefMG(level)%nz
        do i = 0+mod(k+1,2),CoefMG(level)%nx,2
          p = 0
          Ap = 0
          if (i>0) then
                    p = p+PhiMG(level)%Arr(i-1,k)*CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          else if (Btype(We)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(CoefMG(level)%nx-1,k)*CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          end if
          if (i<CoefMG(level)%nx) then
                    p = p+PhiMG(level)%Arr(i+1,k)*CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          else if (Btype(Ea)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(1,k)*CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          end if
          if (k>0) then
                    p = p+PhiMG(level)%Arr(i,k-1)*CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          else if (Btype(Bo)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(i,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          end if
          if (k<CoefMG(level)%nz) then
                    p = p+PhiMG(level)%Arr(i,k+1)*CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          else if (Btype(To)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(i,1)*CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          end if
          p = p-RHSMG(level)%Arr(i,k)
          p = p/Ap

          PhiMG(level)%Arr(i,k) = p
        end do
      end do
      !$omp end do
      !$omp end parallel

    end do
  end subroutine MG_GS

  subroutine MG_res(level,R)
    integer,intent(in) :: level
    real(knd),intent(out) :: R
    integer :: i,k
    real(knd),save :: p,Ap

    call BOUND_Phi_MG(PhiMG(level)%Arr,CoefMG(level)%nx,CoefMG(level)%nz)
    do k = 0,CoefMG(level)%nz
      do i = 0,CoefMG(level)%nx
        p = 0
        Ap = 0
        if (i>0) then
                  p = p+PhiMG(level)%Arr(i-1,k)*CoefMG(level)%Aw
                  Ap = Ap+CoefMG(level)%Aw
        else if (Btype(We)==BC_PERIODIC) then
                  p = p+PhiMG(level)%Arr(CoefMG(level)%nx-1,k)*CoefMG(level)%Aw
                  Ap = Ap+CoefMG(level)%Aw
        end if
        if (i<CoefMG(level)%nx) then
                  p = p+PhiMG(level)%Arr(i+1,k)*CoefMG(level)%Ae
                  Ap = Ap+CoefMG(level)%Ae
        else if (Btype(Ea)==BC_PERIODIC) then
                  p = p+PhiMG(level)%Arr(1,k)*CoefMG(level)%Ae
                  Ap = Ap+CoefMG(level)%Ae
        end if
        if (k>0) then
                  p = p+PhiMG(level)%Arr(i,k-1)*CoefMG(level)%Ab
                  Ap = Ap+CoefMG(level)%Ab
        else if (Btype(Bo)==BC_PERIODIC) then
                  p = p+PhiMG(level)%Arr(i,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                  Ap = Ap+CoefMG(level)%Ab
        end if
        if (k<CoefMG(level)%nz) then
                  p = p+PhiMG(level)%Arr(i,k+1)*CoefMG(level)%At
                  Ap = Ap+CoefMG(level)%At
        else if (Btype(To)==BC_PERIODIC) then
                  p = p+PhiMG(level)%Arr(i,1)*CoefMG(level)%At
                  Ap = Ap+CoefMG(level)%At
        end if
        p = p-RHSMG(level)%Arr(i,k)
        p=-p +Ap*PhiMG(level)%Arr(i,k)
        ResMG(level)%Arr(i,k) = p
      end do
    end do

    R = MAXVAL(ResMG(level)%Arr(0:CoefMG(level)%nx,0:CoefMG(level)%nz))

  end subroutine MG_res


  subroutine MG_clear(level)
    integer,intent(in) :: level
    integer :: l

    do l = level,minmglevel,-1
      PhiMG(l)%Arr = 0
      RHSMG(l)%Arr = 0
      ResMG(l)%Arr = 0
    end do
  end subroutine



  recursive subroutine MG_CGC(level,eps,ncgc,npre,npost,R)
    integer,intent(in) :: level,ncgc,npre,npost
    real(knd),intent(in) :: eps
    real(knd),intent(out) :: R
    real(knd) :: R1
    integer :: k

    if (level == minmglevel) then !1
      call MG_GE(level)
      call MG_res(level,R)

      if (R>mgepsinnerGS**2) then
        R1 = R
        do k = 1,mgmaxinnerGSiter
          call MG_GS(level, 10)
          call MG_res(level,R)

          if (R < mgepsinnerGS**2) exit
          if (R1<R*1.05_knd) exit
          R1 = R
        end do
      end if
    else
      do k = 1, ncgc !number of recurrent calls
        call MG_GS(level,npre)
        call MG_res(level,R)

        call MG_clear(level-1)

        call Restrict(RHSMG(level-1)%Arr,ResMG(level)%Arr,level-1)

        call MG_CGC(level-1,eps,ncgc,npre,npost,R)

        call Prolongate(PhiMG(level)%Arr,PhiMG(level-1)%Arr,level-1)

        call MG_GS(level,npost)
        call MG_res(level,R)
      end do
    end if

  end subroutine MG_CGC


  subroutine POISSMG2d(Phi,RHS)        !Solves Poisson equation using Successive over-relaxation
    real(knd),dimension(0:,0:,0:),intent(inout) :: Phi
    real(knd),dimension(1:,1:,1:),intent(in) :: RHS
    integer :: i,k,l,nx,nz,sx,sz
    real(knd) :: mgeps,R,Phiref
    integer,save :: called = 0

    mgeps = epsPoisson
    Phi = 0

    if (Btype(Ea)==BC_PERIODIC) then
      sx = 1
    else
      sx = 0
    end if

    if (Btype(To)==BC_PERIODIC) then
      sz = 1
    else
      sz = 0
    end if


    if (called==0) then
      allocate(PhiMG(minmglevel:LMG),RHSMG(minmglevel:LMG),ResMG(minmglevel:LMG),CoefMG(minmglevel:LMG))
      do l = LMG,minmglevel,-1
        CoefMG(l)%nx = bnx*2**l
        CoefMG(l)%nz = bnz*2**l
        CoefMG(l)%dx = dxmin*2**(LMG-l)
        CoefMG(l)%dz = dzmin*2**(LMG-l)
        allocate(PhiMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%nz+1))
        allocate(RHSMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%nz+1))
        allocate(ResMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%nz+1))
        CoefMG(l)%Ae = 1._knd/(CoefMG(l)%dx*CoefMG(l)%dx)
        CoefMG(l)%Aw = 1._knd/(CoefMG(l)%dx*CoefMG(l)%dx)
        CoefMG(l)%At = 1._knd/(CoefMG(l)%dz*CoefMG(l)%dz)
        CoefMG(l)%Ab = 1._knd/(CoefMG(l)%dz*CoefMG(l)%dz)
      end do
      called = 1
    end if

    nx = bnx*2**LMG
    nz = bnz*2**LMG

    if (nx-sx/=Prnx-1.or.nz-sz/=Prnz-1) then
      write (*,*) "Incorrect dimensions, multigrid, vs. grid defined in grid.conf:"
      write (*,*) 0+sx,":",nx,"--",1,":",Prnx
      write (*,*) 0+sz,":",nz,"--",1,":",Prnz
      call error_stop
    end if

    PhiMG(LMG)%Arr = 0
    PhiMG(LMG)%Arr(0+sx:nx,0+sz:nz) = Phi(1:Prnx,1,1:Prnz)
    RHSMG(LMG)%Arr(0+sx:nx,0+sz:nz) = RHS(1:Prnx,1,1:Prnz)

    if (Btype(Ea)==BC_PERIODIC) then
      RHSMG(LMG)%Arr(0,:) = RHSMG(LMG)%Arr(Prnx,:)
      PhiMG(LMG)%Arr(0,:) = PhiMG(LMG)%Arr(Prnx,:)
    end if
    if (Btype(To)==BC_PERIODIC) then
      RHSMG(LMG)%Arr(:,0) = RHSMG(LMG)%Arr(:,Prnz)
      PhiMG(LMG)%Arr(:,0) = PhiMG(LMG)%Arr(:,Prnz)
    end if

    do l = 1,maxPoissoniter
      write (*,*) "MG iteration:",l
      Phiref = sum(PhiMG(LMG)%Arr(0:CoefMG(LMG)%nx,0:CoefMG(LMG)%nz))/&
              ((CoefMG(LMG)%nx+1)*(CoefMG(LMG)%nz+1))
      do k = 0,CoefMG(LMG)%nz
        do i = 0,CoefMG(LMG)%nx
          PhiMG(LMG)%Arr(i,k) = PhiMG(LMG)%Arr(i,k)-Phiref
        end do
      end do
      call Bound_Phi_MG(PhiMG(LMG)%Arr,CoefMG(LMG)%nx,CoefMG(LMG)%nz)
      call MG_CGC(LMG,mgeps,mgncgc,mgnpre,mgnpost,R)

      if (R<mgeps)   exit
      write (*,*) "MG residuum",R
    end do
    write (*,*) "MG residuum",R

    forall(l = 1:Prny) Phi(1:Prnx,l,1:Prnz) = PhiMG(LMG)%Arr(0+sx:nx,0+sz:nz)

  end subroutine POISSMG2d





end module Multigrid2d
