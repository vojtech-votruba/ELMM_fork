module Multigrid

  use Parameters
  use Stop_procedures

  implicit none

  private

  public PoissMG, SetMGParams

  type TMGArr
    real(knd),allocatable,dimension(:,:,:) :: Arr
  end type

  type TCoefs
    real(knd) :: dx,dy,dz,Aw,Ae,As,An,Ab,At
    integer :: nx,ny,nz
  end type


  type(TMGArr),allocatable,dimension(:) :: PhiMG,RHSMG,ResMG
  type(TCoefs),allocatable,dimension(:) :: CoefMg

  integer :: bnx,bny,bnz !Base cube dimensions for Multigrid, Prnx = bnx*2**level
  integer :: LMG !depth of multigrid
  integer :: minmglevel !innermost MG level
  integer :: mgncgc !type of cycling 1..V cycle, 2..W cycle
  integer :: mgnpre
  integer :: mgnpost
  integer :: mgmaxinnerGSiter
  real(knd) :: mgepsinnerGS

  integer,dimension(0:8) :: nxa,nya,nza
  real(knd),dimension(0:8) :: Aw,Ae,As,An,Ab,At

contains


  subroutine SetMGParams(llmg, lminmglevel, lbnx, lbny, lbnz,&
                           lmgncgc, lmgnpre, lmgnpost, lmgmaxinnerGSiter, lmgepsinnerGS)

   integer,intent(in) ::   llmg, lminmglevel, lbnx, lbny, lbnz, lmgncgc, lmgnpre, lmgnpost, lmgmaxinnerGSiter
   real(knd),intent(in) :: lmgepsinnerGS

   lmg = llmg
   minmglevel = lminmglevel
   bnx = lbnx
   bny = lbny
   bnz = lbnz
   mgncgc = lmgncgc
   mgnpre = lmgnpre
   mgnpost = lmgnpost
   mgmaxinnerGSiter = lmgmaxinnerGSiter
   mgepsinnerGS = lmgepsinnerGS
  end subroutine SetMGParams

  subroutine Prolongate(AFine,ACoarse,level)
    integer,intent(in) :: level
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: ACoarse
    real(knd),dimension(-1:,-1:,-1:),intent(inout) :: AFine
    integer :: i,j,k,ii,jj,kk,nx,ny,nz
    real(knd),dimension(-1:1,-1:1,-1:1),save :: Cf
    integer,save :: called = 0


    nx = bnx*2**level !level means from which grid we interpolate
    ny = bny*2**level
    nz = bnz*2**level


    if (called==0) then
      do kk=-1,1
        do jj=-1,1
          do ii=-1,1
            if (ii==0.and.jj==0.and.kk==0) then
                Cf(ii,jj,kk) = 1._knd
            else if (abs(ii)+abs(jj)+abs(kk)==1) then
                Cf(ii,jj,kk) = 1/2._knd
            else if (abs(ii)+abs(jj)+abs(kk)==2) then
                Cf(ii,jj,kk) = 1/4._knd
            else
                Cf(ii,jj,kk) = 1/8._knd
            end if
          end do
        end do
      end do
      called = 0
    end if

    !$omp parallel private (i,j)
    do kk=-1,1
      do jj=-1,1
        do ii=-1,1
          !$omp do
          do k = 0,nz
            do j = 0,ny
              do i = 0,nx
                AFine(2*i+ii,2*j+jj,2*k+kk) = AFine(2*i+ii,2*j+jj,2*k+kk) + Cf(ii,jj,kk) * ACoarse(i,j,k)
              end do
            end do
          end do
          !$omp end do
        end do
      end do
    end do
    !$omp end parallel

    call Bound_Phi_MG(Afine,2*nx,2*ny,2*nz)

  end subroutine Prolongate


  subroutine Restrict(ACoarse,AFine,level)
    integer,intent(in) :: level
    real(knd),dimension(-1:,-1:,-1:),intent(out) :: ACoarse
    real(knd),dimension(-1:,-1:,-1:),intent(inout) :: AFine
    real(knd) :: p,S,w
    integer :: i,j,k,ii,jj,kk,nx,ny,nz
    real(knd) :: weight(-1:1,-1:1,-1:1)


    nx = bnx*2**level !level means on which grid we restrict
    ny = bny*2**level
    nz = bnz*2**level

    weight = 1
    weight(-1,:,:) = weight(-1,:,:) / 2._knd
    weight( 1,:,:) = weight( 1,:,:) / 2._knd
    weight(:,-1,:) = weight(:,-1,:) / 2._knd
    weight(:, 1,:) = weight(:, 1,:) / 2._knd
    weight(:,:,-1) = weight(:,:,-1) / 2._knd
    weight(:,:, 1) = weight(:,:, 1) / 2._knd

    call Bound_Phi_MG(Afine,2*nx,2*ny,2*nz)

    !$omp parallel do private(i,j,k,ii,jj,kk,w,S,p)
    do k = 0,nz
      do j = 0,ny
        do i = 0,nx

            p = 0
            S = 0

            do kk = max(2*k-1, 0), min(2*k+1, 2*nz)
            do jj = max(2*j-1, 0), min(2*j+1, 2*ny)
              do ii = max(2*i-1, 0), min(2*i+1, 2*nx)
              w = weight(ii-2*i, jj-2*j, kk-2*k)
              S = S + w * AFine(ii,jj,kk)
              p = p + w
              end do
            end do
            end do

            ACoarse(i,j,k) = S/p

        end do
      end do
    end do
    !$omp end parallel do
    if (Btype(Ea)==BC_PERIODIC) ACoarse(0,:,:) = ACoarse(nx,:,:)
    if (Btype(No)==BC_PERIODIC) ACoarse(:,0,:) = ACoarse(:,ny,:)
    if (Btype(To)==BC_PERIODIC) ACoarse(:,:,0) = ACoarse(:,:,nz)

 end subroutine Restrict







  pure subroutine BOUND_Phi_MG(Phi,nx,ny,nz)
    real(knd),intent(inout) :: Phi(-1:,-1:,-1:)
    integer,intent(in) :: nx,ny,nz
    integer :: i,j,k


     if (Btype(We)==BC_PERIODIC) then
      do k = 0,nz
       do j = 0,ny                      !Periodic BC
        Phi(-1,j,k) = Phi(nx-1,j,k)
       end do
      end do
    else
      do k = 0,nz
       do j = 0,ny                      !Other BCs
        Phi(-1,j,k) = Phi(0,j,k)
       end do
      end do
     end if

     if (Btype(Ea)==BC_PERIODIC) then
      do k = 0,nz
       do j = 0,ny                      !Periodic BC
        Phi(nx+1,j,k) = Phi(1,j,k)
       end do
      end do
     else
      do k = 0,nz
       do j = 0,ny                      !Other BCs
        Phi(nx+1,j,k) = Phi(nx,j,k)
       end do
      end do
     end if

     if (Btype(So)==BC_PERIODIC) then
     do k = 0,nz
       do i=-1,nx+1                      !Periodic BC
        Phi(i,-1,k) = Phi(i,ny-1,k)
       end do
      end do
     else
      do k = 0,nz
       do i=-1,nx+1                      !Other BCs
        Phi(i,-1,k) = Phi(i,0,k)
       end do
      end do
     end if

     if (Btype(No)==BC_PERIODIC) then
     do k = 0,nz
       do i=-1,nx+1                      !Periodic BC
        Phi(i,ny+1,k) = Phi(i,1,k)
       end do
      end do
     else
      do k = 0,nz
       do i=-1,nx+1                      !Other BCs
        Phi(i,ny+1,k) = Phi(i,ny,k)
       end do
      end do
     end if

     if (Btype(Bo)==BC_PERIODIC) then
      do j=-1,ny+1
       do i=-1,nx+1                      !Periodic BC
        Phi(i,j,-1) = Phi(i,j,nz-1)
       end do
      end do
     else
      do j=-1,ny+1
       do i=-1,nx+1                      !Other BCs
        Phi(i,j,-1) = Phi(i,j,0)
       end do
      end do
     end if

     if (Btype(To)==BC_PERIODIC) then
      do j=-1,ny+1
       do i=-1,nx+1                      !Periodic BC
        Phi(i,j,nz+1) = Phi(i,j,1)
       end do
      end do
     else
      do j=-1,ny+1
       do i=-1,nx+1                      !Other BCs
        Phi(i,j,nz+1) = Phi(i,j,nz)
       end do
      end do
     end if
  end subroutine BOUND_Phi_MG



  pure function ind(level,nulx,nuly,nulz,i,j,k)
    integer :: ind
    integer,intent(in) :: level,i,j,k,nulx,nuly,nulz
    ind = (1-nulx)+i+(j-nuly)*(CoefMG(level)%nx+(1-nulx))+(k-nulz)*(CoefMG(level)%nx+1-nulx)*(CoefMG(level)%ny+1-nuly)
  end function ind




  subroutine MG_GE(level)
    use Lapack
    integer,intent(in) :: level
    integer :: i,j,k,l,info
    integer,allocatable,dimension(:),save :: ige,ipivot,work2
    real(knd),allocatable,dimension(:),save :: xge,bge,work,R,C,ferr,berr
    real(knd),allocatable,dimension(:,:),save :: age,af
    real(knd) :: Ap,rcond
    integer,save :: nx,ny,nz,nulx,nuly,nulz,nxyz,called = 0
    character(1),save :: equed

    if (called==0) then
      nx = CoefMG(level)%nx
      ny = CoefMG(level)%ny
      nz = CoefMG(level)%nz

      if (Btype(Ea)==BC_PERIODIC) then
        nulx = 1
      else
        nulx = 0
      end if
      if (Btype(No)==BC_PERIODIC) then
        nuly = 1
      else
        nuly = 0
      end if
      if (Btype(To)==BC_PERIODIC) then
        nulz = 1
      else
        nulz = 0
      end if
      nxyz = (nx+(1-nulx))*(ny+(1-nuly))*(nz+(1-nulz))
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
        do j = nuly,ny
          do i = nulx,nx
            l = ind(level,nulx,nuly,nulz,i,j,k)
            Ap = 0
            if (i>nulx) then
                      age(l,ind(level,nulx,nuly,nulz,i-1,j,k)) = CoefMG(level)%Aw
                      Ap = Ap+CoefMG(level)%Aw
            else if (Btype(We)==BC_PERIODIC) then
                      age(l,ind(level,nulx,nuly,nulz,nx,j,k)) = CoefMG(level)%Aw
                      Ap = Ap+CoefMG(level)%Aw
            end if
            if (i<nx) then
                      age(l,ind(level,nulx,nuly,nulz,i+1,j,k)) = CoefMG(level)%Ae
                      Ap = Ap+CoefMG(level)%Ae
            else if (Btype(Ea)==BC_PERIODIC) then
                      age(l,ind(level,nulx,nuly,nulz,nulx,j,k)) = CoefMG(level)%Ae
                      Ap = Ap+CoefMG(level)%Ae
            end if
            if (j>nuly) then
                      age(l,ind(level,nulx,nuly,nulz,i,j-1,k)) = CoefMG(level)%As
                      Ap = Ap+CoefMG(level)%As
            else if (Btype(So)==BC_PERIODIC) then
                      age(l,ind(level,nulx,nuly,nulz,i,ny,k)) = CoefMG(level)%As
                      Ap = Ap+CoefMG(level)%As
            end if
            if (j<ny) then
                      age(l,ind(level,nulx,nuly,nulz,i,j+1,k)) = CoefMG(level)%An
                      Ap = Ap+CoefMG(level)%An
            else if (Btype(No)==BC_PERIODIC) then
                      age(l,ind(level,nulx,nuly,nulz,i,nuly,k)) = CoefMG(level)%An
                      Ap = Ap+CoefMG(level)%An
            end if
            if (k>nulz) then
                      age(l,ind(level,nulx,nuly,nulz,i,j,k-1)) = CoefMG(level)%Ab
                      Ap = Ap+CoefMG(level)%Ab
            else if (Btype(Bo)==BC_PERIODIC) then
                      age(l,ind(level,nulx,nuly,nulz,i,j,nz)) = CoefMG(level)%Ab
                      Ap = Ap+CoefMG(level)%Ab
            end if
            if (k<nz) then
                      age(l,ind(level,nulx,nuly,nulz,i,j,k+1)) = CoefMG(level)%At
                      Ap = Ap+CoefMG(level)%At
            else if (Btype(To)==BC_PERIODIC) then
                      age(l,ind(level,nulx,nuly,nulz,i,j,nulz)) = CoefMG(level)%At
                      Ap = Ap+CoefMG(level)%At
            end if

            age(l,l)=-Ap
            bge(l) = RHSMG(level)%Arr(i,j,k)
          end do
        end do
      end do

      l = ind(level,nulx,nuly,nulz,nx,ny,nz)
      age(1,:) = 0
      age(1,1) = age(2,2)
      bge(1) = 0

      call gesvx("E","N",nxyz,1,age,nxyz,af,nxyz,ipivot,EQUED,R,C,bge,nxyz,xge,nxyz, &
                rcond,ferr,berr,work,work2,info)

      write (*,*) "RCOND",rcond
      if (info/=0) then
        write (*,*) "info",info
        call error_stop
      end if
    else


    do k = nulz,nz
      do j = nuly,ny
          do i = nulx,nx
            l = ind(level,nulx,nuly,nulz,i,j,k)
            bge(l) = RHSMG(level)%Arr(i,j,k)
          end do
      end do
    end do
    bge(1) = 0

    call gesvx("F","N",nxyz,1,age,nxyz,af,nxyz,ipivot,EQUED,R,C,bge,nxyz,xge,nxyz, &
                rcond,ferr,berr,work,work2,info)

    if (info/=0) then
      write (*,*) "info",info
      call error_stop
    end if
  end if


    do k = nulz,nz
      do j = nuly,ny
        do i = nulx,nx
          l = ind(level,nulx,nuly,nulz,i,j,k)
          PhiMG(level)%Arr(i,j,k) = xge(l)
        end do
      end do
    end do

    if (Btype(Ea)==BC_PERIODIC) then
      PhiMG(level)%Arr(0,:,:) = PhiMG(level)%Arr(nx,:,:)
    end if
    if (Btype(No)==BC_PERIODIC) then
      PhiMG(level)%Arr(:,0,:) = PhiMG(level)%Arr(:,ny,:)
    end if
    if (Btype(To)==BC_PERIODIC) then
      PhiMG(level)%Arr(:,:,0) = PhiMG(level)%Arr(:,:,nz)
    end if

  end subroutine MG_GE


  subroutine MG_LU(level)
    use Lapack
    integer,intent(in) :: level
    integer :: i,j,k,l,info
    integer,allocatable,dimension(:),save :: ige,ipivot,work2
    real(knd),allocatable,dimension(:),save :: xge,bge,work,R,C,ferr,berr
    real(knd),allocatable,dimension(:,:),save :: age,af
    real(knd) :: Ap
    integer,save :: nx,ny,nz,nulx,nuly,nulz,nxyz,called = 0

    if (called==0) then
     nx = CoefMG(level)%nx
     ny = CoefMG(level)%ny
     nz = CoefMG(level)%nz

     if (Btype(Ea)==BC_PERIODIC) then
      nulx = 1
     else
      nulx = 0
     end if
     if (Btype(No)==BC_PERIODIC) then
      nuly = 1
     else
      nuly = 0
     end if
     if (Btype(To)==BC_PERIODIC) then
      nulz = 1
     else
      nulz = 0
     end if
     nxyz = (nx+(1-nulx))*(ny+(1-nuly))*(nz+(1-nulz))
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
           do j = nuly,ny
              do i = nulx,nx
                l = ind(level,nulx,nuly,nulz,i,j,k)
                Ap = 0
                if (i>nulx) then
                          age(l,ind(level,nulx,nuly,nulz,i-1,j,k)) = CoefMG(level)%Aw
                          Ap = Ap+CoefMG(level)%Aw
                else if (Btype(We)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,nx,j,k)) = CoefMG(level)%Aw
                          Ap = Ap+CoefMG(level)%Aw
                end if
                if (i<nx) then
                          age(l,ind(level,nulx,nuly,nulz,i+1,j,k)) = CoefMG(level)%Ae
                          Ap = Ap+CoefMG(level)%Ae
                else if (Btype(Ea)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,nulx,j,k)) = CoefMG(level)%Ae
                          Ap = Ap+CoefMG(level)%Ae
                end if
                if (j>nuly) then
                          age(l,ind(level,nulx,nuly,nulz,i,j-1,k)) = CoefMG(level)%As
                          Ap = Ap+CoefMG(level)%As
                else if (Btype(So)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,i,ny,k)) = CoefMG(level)%As
                          Ap = Ap+CoefMG(level)%As
                end if
                if (j<ny) then
                          age(l,ind(level,nulx,nuly,nulz,i,j+1,k)) = CoefMG(level)%An
                          Ap = Ap+CoefMG(level)%An
                else if (Btype(No)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,i,nuly,k)) = CoefMG(level)%An
                          Ap = Ap+CoefMG(level)%An
                end if
                if (k>nulz) then
                          age(l,ind(level,nulx,nuly,nulz,i,j,k-1)) = CoefMG(level)%Ab
                          Ap = Ap+CoefMG(level)%Ab
                else if (Btype(Bo)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,i,j,nz)) = CoefMG(level)%Ab
                          Ap = Ap+CoefMG(level)%Ab
                end if
                if (k<nz) then
                          age(l,ind(level,nulx,nuly,nulz,i,j,k+1)) = CoefMG(level)%At
                          Ap = Ap+CoefMG(level)%At
                else if (Btype(To)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,i,j,nulz)) = CoefMG(level)%At
                          Ap = Ap+CoefMG(level)%At
                end if

                age(l,l)=-Ap
               end do
           end do
       end do

       l = ind(level,nulx,nuly,nulz,nx,ny,nz)
       age(1,:) = 0
       age(1,1) = age(2,2)

       call  getrf(nxyz,nxyz,age,nxyz,ipivot,info)

       if (info/=0) then
          write (*,*) "info",info
          call error_stop
       end if
    end if


    do k = nulz,nz
      do j = nuly,ny
          do i = nulx,nx
            l = ind(level,nulx,nuly,nulz,i,j,k)
            bge(l) = RHSMG(level)%Arr(i,j,k)
          end do
      end do
    end do
    bge(1) = 0

    call getrs("N",nxyz,1,age,nxyz,ipivot,bge,nxyz,info)

    if (info/=0) then
      write (*,*) "info",info
      call error_stop
    end if


    do k = nulz,nz
      do j = nuly,ny
        do i = nulx,nx
          l = ind(level,nulx,nuly,nulz,i,j,k)
          PhiMG(level)%Arr(i,j,k) = bge(l)
        end do
      end do
    end do

    if (Btype(Ea)==BC_PERIODIC) then
      PhiMG(level)%Arr(0,:,:) = PhiMG(level)%Arr(nx,:,:)
    end if
    if (Btype(No)==BC_PERIODIC) then
      PhiMG(level)%Arr(:,0,:) = PhiMG(level)%Arr(:,ny,:)
    end if
    if (Btype(To)==BC_PERIODIC) then
      PhiMG(level)%Arr(:,:,0) = PhiMG(level)%Arr(:,:,nz)
    end if

  end subroutine MG_LU





  subroutine MG_INV(level)
    use Lapack
    integer,intent(in) :: level
    integer :: i,j,k,l,info,ldwork
    integer,allocatable,dimension(:) :: ipivot
    real(knd),allocatable,dimension(:),save :: xge,bge,work
    real(knd),allocatable,dimension(:,:),save :: age
    real(knd) :: Ap
    integer,save :: nx,ny,nz,nulx,nuly,nulz,nxyz,called = 0

    if (called==0) then
     nx = CoefMG(level)%nx
     ny = CoefMG(level)%ny
     nz = CoefMG(level)%nz

     if (Btype(Ea)==BC_PERIODIC) then
      nulx = 1
     else
      nulx = 0
     end if
     if (Btype(No)==BC_PERIODIC) then
      nuly = 1
     else
      nuly = 0
     end if
     if (Btype(To)==BC_PERIODIC) then
      nulz = 1
     else
      nulz = 0
     end if
     nxyz = (nx+(1-nulx))*(ny+(1-nuly))*(nz+(1-nulz))
     allocate(xge(nxyz))
     allocate(bge(nxyz))
     allocate(ipivot(nxyz))
     allocate(age(nxyz,nxyz))


     called = 1


     age = 0
   !   bge = 0
     xge = 0
        do k = nulz,nz
           do j = nuly,ny
              do i = nulx,nx
                l = ind(level,nulx,nuly,nulz,i,j,k)
                Ap = 0
                if (i>nulx) then
                          age(l,ind(level,nulx,nuly,nulz,i-1,j,k)) = CoefMG(level)%Aw
                          Ap = Ap+CoefMG(level)%Aw
                else if (Btype(We)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,nx,j,k)) = CoefMG(level)%Aw
                          Ap = Ap+CoefMG(level)%Aw
                end if
                if (i<nx) then
                          age(l,ind(level,nulx,nuly,nulz,i+1,j,k)) = CoefMG(level)%Ae
                          Ap = Ap+CoefMG(level)%Ae
                else if (Btype(Ea)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,nulx,j,k)) = CoefMG(level)%Ae
                          Ap = Ap+CoefMG(level)%Ae
                end if
                if (j>nuly) then
                          age(l,ind(level,nulx,nuly,nulz,i,j-1,k)) = CoefMG(level)%As
                          Ap = Ap+CoefMG(level)%As
                else if (Btype(So)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,i,ny,k)) = CoefMG(level)%As
                          Ap = Ap+CoefMG(level)%As
                end if
                if (j<ny) then
                          age(l,ind(level,nulx,nuly,nulz,i,j+1,k)) = CoefMG(level)%An
                          Ap = Ap+CoefMG(level)%An
                else if (Btype(No)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,i,nuly,k)) = CoefMG(level)%An
                          Ap = Ap+CoefMG(level)%An
                end if
                if (k>nulz) then
                          age(l,ind(level,nulx,nuly,nulz,i,j,k-1)) = CoefMG(level)%Ab
                          Ap = Ap+CoefMG(level)%Ab
                else if (Btype(Bo)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,i,j,nz)) = CoefMG(level)%Ab
                          Ap = Ap+CoefMG(level)%Ab
                end if
                if (k<nz) then
                          age(l,ind(level,nulx,nuly,nulz,i,j,k+1)) = CoefMG(level)%At
                          Ap = Ap+CoefMG(level)%At
                else if (Btype(To)==BC_PERIODIC) then
                          age(l,ind(level,nulx,nuly,nulz,i,j,nulz)) = CoefMG(level)%At
                          Ap = Ap+CoefMG(level)%At
                end if

                age(l,l)=-Ap
   !              bge(l) = RHSMG(level)%Arr(i,j,k)
               end do
           end do
       end do

       l = ind(level,nulx,nuly,nulz,nx,ny,nz)
       age(1,:) = 0
       age(1,1) = age(2,2)
   !     bge(1) = 0

       call  getrf(nxyz,nxyz,age,nxyz,ipivot,info)

       if (info/=0) then
        call error_stop
        write (*,*) "info",info
       end if

       allocate(work(1))

       call getri(nxyz,age,nxyz,ipivot,work,-1,info)

       if (info/=0) then
        call error_stop
        write (*,*) "info",info
       end if

       ldwork = nint(work(1))

       deallocate(work)
       allocate(work(ldwork))

       call getri(nxyz,age,nxyz,ipivot,work,ldwork,info)

       if (info/=0) then
        call error_stop
        write (*,*) "info",info
       end if

    end if


       do k = nulz,nz
          do j = nuly,ny
             do i = nulx,nx
               l = ind(level,nulx,nuly,nulz,i,j,k)
               bge(l) = RHSMG(level)%Arr(i,j,k)
             end do
          end do
       end do
       bge(1) = 0


   !     do j = 1,nxyz !rows of X{j}
   !      xge(j) = 0
   !      do i = 1,nxyz !columns of A(i,j)
   !       xge(j) = xge(i,j)+age(j,i)*b(i)
   !      end do
   !     end do

       xge = 0
       do j = 1,nxyz !rows of X{j}
        do i = 1,nxyz !columns of A(i,j) if (age(i,j)>100*tiny(1._knd))
         xge(i) = xge(i)+age(i,j)*bge(j)
        end do
       end do



        do k = nulz,nz
           do j = nuly,ny
              do i = nulx,nx
                l = ind(level,nulx,nuly,nulz,i,j,k)
                PhiMG(level)%Arr(i,j,k) = xge(l)
              end do
           end do
        end do

    if (Btype(Ea)==BC_PERIODIC) then
     PhiMG(level)%Arr(0,:,:) = PhiMG(level)%Arr(nx,:,:)
    end if
    if (Btype(No)==BC_PERIODIC) then
     PhiMG(level)%Arr(:,0,:) = PhiMG(level)%Arr(:,ny,:)
    end if
    if (Btype(To)==BC_PERIODIC) then
     PhiMG(level)%Arr(:,:,0) = PhiMG(level)%Arr(:,:,nz)
    end if

  end subroutine MG_INV








  subroutine MG_GS(level,niter)
    integer,intent(in) :: level,niter
    integer :: i,j,k,l
    real(knd) :: p,Ap


    do l = 1,niter
      !$omp parallel private(i,j,k,p,Ap)
      !$omp do
      do k = 0,CoefMG(level)%nz
          do j = 0,CoefMG(level)%ny
              do i = 0+mod(j+k,2),CoefMG(level)%nx,2
                p = 0
                Ap = 0
                if (i>0) then
                          p = p+PhiMG(level)%Arr(i-1,j,k)*CoefMG(level)%Aw
                          Ap = Ap+CoefMG(level)%Aw
                else if (Btype(We)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(CoefMG(level)%nx-1,j,k)*CoefMG(level)%Aw
                          Ap = Ap+CoefMG(level)%Aw
                end if
                if (i<CoefMG(level)%nx) then
                          p = p+PhiMG(level)%Arr(i+1,j,k)*CoefMG(level)%Ae
                          Ap = Ap+CoefMG(level)%Ae
                else if (Btype(Ea)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(1,j,k)*CoefMG(level)%Ae
                          Ap = Ap+CoefMG(level)%Ae
                end if
                if (j>0) then
                          p = p+PhiMG(level)%Arr(i,j-1,k)*CoefMG(level)%As
                          Ap = Ap+CoefMG(level)%As
                else if (Btype(So)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(i,CoefMG(level)%ny-1,k)*CoefMG(level)%As
                          Ap = Ap+CoefMG(level)%As
                end if
                if (j<CoefMG(level)%ny) then
                          p = p+PhiMG(level)%Arr(i,j+1,k)*CoefMG(level)%An
                          Ap = Ap+CoefMG(level)%An
                else if (Btype(No)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(i,1,k)*CoefMG(level)%An
                          Ap = Ap+CoefMG(level)%An
                end if
                if (k>0) then
                          p = p+PhiMG(level)%Arr(i,j,k-1)*CoefMG(level)%Ab
                          Ap = Ap+CoefMG(level)%Ab
                else if (Btype(Bo)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(i,j,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                          Ap = Ap+CoefMG(level)%Ab
                end if
                if (k<CoefMG(level)%nz) then
                          p = p+PhiMG(level)%Arr(i,j,k+1)*CoefMG(level)%At
                          Ap = Ap+CoefMG(level)%At
                else if (Btype(To)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(i,j,1)*CoefMG(level)%At
                          Ap = Ap+CoefMG(level)%At
                end if
                p = p-RHSMG(level)%Arr(i,j,k)

                p = p/Ap
                PhiMG(level)%Arr(i,j,k) = p
              end do
          end do
      end do
      !$omp end do
      !$omp do
      do k = 0,CoefMG(level)%nz
          do j = 0,CoefMG(level)%ny
              do i = 0+mod(j+k+1,2),CoefMG(level)%nx,2
                p = 0
                Ap = 0
                if (i>0) then
                          p = p+PhiMG(level)%Arr(i-1,j,k)*CoefMG(level)%Aw
                          Ap = Ap+CoefMG(level)%Aw
                else if (Btype(We)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(CoefMG(level)%nx-1,j,k)*CoefMG(level)%Aw
                          Ap = Ap+CoefMG(level)%Aw
                end if
                if (i<CoefMG(level)%nx) then
                          p = p+PhiMG(level)%Arr(i+1,j,k)*CoefMG(level)%Ae
                          Ap = Ap+CoefMG(level)%Ae
                else if (Btype(Ea)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(1,j,k)*CoefMG(level)%Ae
                          Ap = Ap+CoefMG(level)%Ae
                end if
                if (j>0) then
                          p = p+PhiMG(level)%Arr(i,j-1,k)*CoefMG(level)%As
                          Ap = Ap+CoefMG(level)%As
                else if (Btype(So)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(i,CoefMG(level)%ny-1,k)*CoefMG(level)%As
                          Ap = Ap+CoefMG(level)%As
                end if
                if (j<CoefMG(level)%ny) then
                          p = p+PhiMG(level)%Arr(i,j+1,k)*CoefMG(level)%An
                          Ap = Ap+CoefMG(level)%An
                else if (Btype(No)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(i,1,k)*CoefMG(level)%An
                          Ap = Ap+CoefMG(level)%An
                end if
                if (k>0) then
                          p = p+PhiMG(level)%Arr(i,j,k-1)*CoefMG(level)%Ab
                          Ap = Ap+CoefMG(level)%Ab
                else if (Btype(Bo)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(i,j,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                          Ap = Ap+CoefMG(level)%Ab
                end if
                if (k<CoefMG(level)%nz) then
                          p = p+PhiMG(level)%Arr(i,j,k+1)*CoefMG(level)%At
                          Ap = Ap+CoefMG(level)%At
                else if (Btype(To)==BC_PERIODIC) then
                          p = p+PhiMG(level)%Arr(i,j,1)*CoefMG(level)%At
                          Ap = Ap+CoefMG(level)%At
                end if
                p = p-RHSMG(level)%Arr(i,j,k)
                p = p/Ap

                PhiMG(level)%Arr(i,j,k) = p
              end do
          end do
      end do
      !$omp end do
      !$omp end parallel

    end do

  end subroutine MG_GS

  subroutine MG_res(level,R)
    integer,intent(in) :: level
    real(knd),intent(out) :: R
    integer :: i,j,k
    real(knd),save :: p,Ap


    !$omp parallel do private(p,Ap)
    do k = 0,CoefMG(level)%nz
      do j = 0,CoefMG(level)%ny
        do i = 0,CoefMG(level)%nx
          p = 0
          Ap = 0
          if (i>0) then
                    p = p+PhiMG(level)%Arr(i-1,j,k)*CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          else if (Btype(We)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(CoefMG(level)%nx-1,j,k)*CoefMG(level)%Aw
                    Ap = Ap+CoefMG(level)%Aw
          end if
          if (i<CoefMG(level)%nx) then
                    p = p+PhiMG(level)%Arr(i+1,j,k)*CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          else if (Btype(Ea)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(1,j,k)*CoefMG(level)%Ae
                    Ap = Ap+CoefMG(level)%Ae
          end if
          if (j>0) then
                    p = p+PhiMG(level)%Arr(i,j-1,k)*CoefMG(level)%As
                    Ap = Ap+CoefMG(level)%As
          else if (Btype(So)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(i,CoefMG(level)%ny-1,k)*CoefMG(level)%As
                    Ap = Ap+CoefMG(level)%As
          end if
          if (j<CoefMG(level)%ny) then
                    p = p+PhiMG(level)%Arr(i,j+1,k)*CoefMG(level)%An
                    Ap = Ap+CoefMG(level)%An
          else if (Btype(No)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(i,1,k)*CoefMG(level)%An
                    Ap = Ap+CoefMG(level)%An
          end if
          if (k>0) then
                    p = p+PhiMG(level)%Arr(i,j,k-1)*CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          else if (Btype(Bo)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(i,j,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                    Ap = Ap+CoefMG(level)%Ab
          end if
          if (k<CoefMG(level)%nz) then
                    p = p+PhiMG(level)%Arr(i,j,k+1)*CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          else if (Btype(To)==BC_PERIODIC) then
                    p = p+PhiMG(level)%Arr(i,j,1)*CoefMG(level)%At
                    Ap = Ap+CoefMG(level)%At
          end if
          p = p-RHSMG(level)%Arr(i,j,k)
          p=-p +Ap*PhiMG(level)%Arr(i,j,k)
          ResMG(level)%Arr(i,j,k) = p
        end do
      end do
    end do
    !$omp end parallel do

    R = MAXVAL(ResMG(level)%Arr(0:CoefMG(level)%nx,0:CoefMG(level)%ny,0:CoefMG(level)%nz))

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


    if (level == minmglevel) then
      R = huge(mgepsinnerGS)/2.


      call MG_GE(level)
      call MG_res(level,R)


      if (R>mgepsinnerGS**2) then
        R1 = R
        do k = 1,mgmaxinnerGSiter

          call MG_GS(level, 10)
          call MG_Res(level, R)

          if (R < mgepsinnerGS) exit
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


  subroutine PoissMG(Phi,RHS)        !Solves Poisson equation using Successive over-relaxation
    real(knd),dimension(0:,0:,0:),intent(inout) :: Phi
    real(knd),dimension(1:,1:,1:),intent(in) :: RHS
    integer :: l,nx,ny,nz,sx,sy,sz
    real(knd) :: mgeps,R
    integer,save :: called = 0

    mgeps = epsPoisson
    Phi = 0

    if (Btype(Ea)==BC_PERIODIC) then
      sx = 1
    else
      sx = 0
    end if
    if (Btype(No)==BC_PERIODIC) then
      sy = 1
    else
      sy = 0
    end if
    if (Btype(To)==BC_PERIODIC) then
      sz = 1
    else
      sz = 0
    end if


    if (called==0) then
      allocate(PhiMG(0:8),RHSMG(0:8),ResMG(0:8),CoefMG(0:8))

      do l = 8,0,-1
        if (l>LMG.or.l<minmglevel) then
          CoefMG(l)%nx = 0
          CoefMG(l)%ny = 0
          CoefMG(l)%nz = 0
          allocate(PhiMG(l)%Arr(1,1,1))
          allocate(RHSMG(l)%Arr(1,1,1))
          allocate(ResMG(l)%Arr(1,1,1))
        else
          CoefMG(l)%nx = bnx*2**l
          CoefMG(l)%ny = bny*2**l
          CoefMG(l)%nz = bnz*2**l
          CoefMG(l)%dx = dxmin*2**(LMG-l)
          CoefMG(l)%dy = dymin*2**(LMG-l)
          CoefMG(l)%dz = dzmin*2**(LMG-l)
          allocate(PhiMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%ny+1,-1:CoefMG(l)%nz+1))
          allocate(RHSMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%ny+1,-1:CoefMG(l)%nz+1))
          allocate(ResMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%ny+1,-1:CoefMG(l)%nz+1))
          CoefMG(l)%Ae = 1._knd/(CoefMG(l)%dx*CoefMG(l)%dx)
          CoefMG(l)%Aw = 1._knd/(CoefMG(l)%dx*CoefMG(l)%dx)
          CoefMG(l)%An = 1._knd/(CoefMG(l)%dy*CoefMG(l)%dy)
          CoefMG(l)%As = 1._knd/(CoefMG(l)%dy*CoefMG(l)%dy)
          CoefMG(l)%At = 1._knd/(CoefMG(l)%dz*CoefMG(l)%dz)
          CoefMG(l)%Ab = 1._knd/(CoefMG(l)%dz*CoefMG(l)%dz)
        end if
      end do

      called = 1
    end if

    nx = bnx*2**LMG
    ny = bny*2**LMG
    nz = bnz*2**LMG

    if (nx-sx/=Prnx-1.or.ny-sy/=Prny-1.or.nz-sz/=Prnz-1) then
       write (*,*) "Incorrect dimensions, multigrid, vs. grid defined in grid.conf:"
       write (*,*) 0+sx, ":", nx, "--", 1, ":", Prnx
       write (*,*) 0+sy, ":", ny, "--", 1, ":", Prny
       write (*,*) 0+sz, ":", nz, "--", 1, ":", Prnz
       call error_stop
    end if


    PhiMG(LMG)%Arr = 0
    PhiMG(LMG)%Arr(0+sx:nx,0+sy:ny,0+sz:nz) = Phi(1:Prnx,1:Prny,1:Prnz)
    RHSMG(LMG)%Arr(0+sx:nx,0+sy:ny,0+sz:nz) = RHS(1:Prnx,1:Prny,1:Prnz)

    if (Btype(Ea)==BC_PERIODIC) then
      RHSMG(LMG)%Arr(0,:,:) = RHSMG(LMG)%Arr(Prnx,:,:)
      PhiMG(LMG)%Arr(0,:,:) = PhiMG(LMG)%Arr(Prnx,:,:)
    end if
    if (Btype(No)==BC_PERIODIC) then
      RHSMG(LMG)%Arr(:,0,:) = RHSMG(LMG)%Arr(:,Prny,:)
      PhiMG(LMG)%Arr(:,0,:) = PhiMG(LMG)%Arr(:,Prny,:)
    end if
    if (Btype(To)==BC_PERIODIC) then
      RHSMG(LMG)%Arr(:,:,0) = RHSMG(LMG)%Arr(:,:,Prnz)
      PhiMG(LMG)%Arr(:,:,0) = PhiMG(LMG)%Arr(:,:,Prnz)
    end if





    nxa=[ CoefMG(0)%nx, CoefMG(1)%nx, CoefMG(2)%nx, CoefMG(3)%nx, CoefMG(4)%nx,&
          CoefMG(5)%nx, CoefMG(6)%nx, CoefMG(7)%nx, CoefMG(8)%nx ]
    nya=[ CoefMG(0)%ny, CoefMG(1)%ny, CoefMG(2)%ny, CoefMG(3)%ny, CoefMG(4)%ny,&
          CoefMG(5)%ny, CoefMG(6)%ny, CoefMG(7)%ny, CoefMG(8)%ny ]
    nza=[ CoefMG(0)%nz, CoefMG(1)%nz, CoefMG(2)%nz, CoefMG(3)%nz, CoefMG(4)%nz,&
          CoefMG(5)%nz, CoefMG(6)%nz, CoefMG(7)%nz, CoefMG(8)%nz ]

    Aw=[ CoefMG(0)%Aw, CoefMG(1)%Aw, CoefMG(2)%Aw, CoefMG(3)%Aw, CoefMG(4)%Aw,&
         CoefMG(5)%Aw, CoefMG(6)%Aw, CoefMG(7)%Aw, CoefMG(8)%Aw ]
    Ae=[ CoefMG(0)%Ae, CoefMG(1)%Ae, CoefMG(2)%Ae, CoefMG(3)%Ae, CoefMG(4)%Ae,&
         CoefMG(5)%Ae, CoefMG(6)%Ae, CoefMG(7)%Ae, CoefMG(8)%Ae ]
    As=[ CoefMG(0)%As, CoefMG(1)%As, CoefMG(2)%As, CoefMG(3)%As, CoefMG(4)%As,&
         CoefMG(5)%As, CoefMG(6)%As, CoefMG(7)%As, CoefMG(8)%As ]
    An=[ CoefMG(0)%An, CoefMG(1)%An, CoefMG(2)%An, CoefMG(3)%An, CoefMG(4)%An,&
         CoefMG(5)%An, CoefMG(6)%An, CoefMG(7)%An, CoefMG(8)%An ]
    Ab=[ CoefMG(0)%Ab, CoefMG(1)%Ab, CoefMG(2)%Ab, CoefMG(3)%Ab, CoefMG(4)%Ab,&
         CoefMG(5)%Ab, CoefMG(6)%Ab, CoefMG(7)%Ab, CoefMG(8)%Ab ]
    At=[ CoefMG(0)%At, CoefMG(1)%At, CoefMG(2)%At, CoefMG(3)%At, CoefMG(4)%At,&
         CoefMG(5)%At, CoefMG(6)%At, CoefMG(7)%At, CoefMG(8)%At ]

    do l = 1,maxPoissoniter
      write (*,*) "MG iteration:",l

      call MG_CGC(LMG,mgeps,mgncgc,mgnpre,mgnpost,R)

      write (*,*) "MG residuum",R
      if (R<mgeps)   exit
    end do

    Phi(1:Prnx,1:Prny,1:Prnz) = PhiMG(LMG)%Arr(0+sx:nx,0+sy:ny,0+sz:nz)

  end subroutine PoissMG


end module Multigrid
