module exchange_par
  use custom_par, only: domain_comm, &
                        w_rank, e_rank, s_rank, n_rank, b_rank, t_rank, neigh_ranks ,&
                        nxims, nyims, nzims, iim, jim, kim, myrank

  use custom_mpi, only: MPI_knd, MPI_STATUS_SIZE, MPI_REQUEST_NULL, &
                        MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE

  use Parameters, only: We, Ea, So, No, Bo, To, BC_MPI_PERIODIC

  use Kinds
  
  use exchange_mpi_derived_types
  
  
  implicit none
  
  private
  
  public par_exchange_boundaries, par_exchange_boundaries_xz, par_exchange_boundaries_yz, &
         par_exchange_Pr, &
         par_exchange_Q, par_exchange_Sc_x, par_exchange_Sc_y, par_exchange_Sc_z, &
         par_exchange_U_x, par_exchange_U_y, par_exchange_U_z, &
         par_exchange_UVW, &
         par_init_exchange

  interface
    subroutine MPI_Recv(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)
      import
      real(knd) :: BUF(*)
      integer   :: COUNT, DATATYPE, SOURCE, TAG, COMM
      integer   :: STATUS(MPI_STATUS_SIZE), IERROR
    end subroutine
    subroutine MPI_Send(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      import
      real(knd) :: BUF(*)
      integer   :: COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    end subroutine
    subroutine MPI_IRecv(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR)
      import
      real(knd) :: BUF(*)
      integer   :: COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR
    end subroutine
    subroutine MPI_ISend(BUF, COUNT, DATATYPE, DEST, TAG, COMM, REQUEST, IERROR)
      import
      real(knd) :: BUF(*)
      integer   :: COUNT, DATATYPE, DEST, TAG, COMM, REQUEST, IERROR
    end subroutine
    subroutine MPI_Waitall(COUNT, ARRAY_OF_REQUESTS, ARRAY_OF_STATUSES, IERROR)
      import
      integer :: COUNT, ARRAY_OF_REQUESTS(*)
      integer :: ARRAY_OF_STATUSES(MPI_STATUS_SIZE,*), IERROR
    end subroutine
  end interface
  
      
contains

  subroutine par_init_exchange
    use Parameters, only: Prnx, Prny, Prnz, &
                          Unx, Uny, Unz, &
                          Vnx, Vny, Vnz, &
                          Wnx, Wny, Wnz
                          
    call init_mpi_derived_types(send_mpi_types(:,1), recv_mpi_types(:,1), Unx, Uny, Unz, -2, 3)
    call init_mpi_derived_types(send_mpi_types(:,2), recv_mpi_types(:,2), Vnx, Vny, Vnz, -2, 3)
    call init_mpi_derived_types(send_mpi_types(:,3), recv_mpi_types(:,3), Wnx, Wny, Wnz, -2, 3)  
    call init_mpi_derived_types(send_mpi_types(:,4), recv_mpi_types(:,4), Prnx, Prny, Prnz, -1, 2)
    call init_mpi_derived_types(send_mpi_types(:,6), recv_mpi_types(:,6), Prnx, Prny, Prnz, -1, 2)
    call init_mpi_derived_types_Q(send_mpi_types(:,7), recv_mpi_types(:,7), Prnx, Prny, Prnz)
  end subroutine

  
  subroutine par_exchange_boundaries(Phi, Btype, component, dir)
    real(knd), intent(inout), contiguous :: Phi(:,:,:)
    integer, intent(in) :: Btype(6)
    integer, intent(in) :: component
    integer, intent(in), optional :: dir
    integer :: ierr
    logical :: xdir, ydir, zdir
    integer :: requests(12)
    
    requests = MPI_REQUEST_NULL
    
    ierr = 0
    
    if (present(dir)) then
      xdir = .false.
      ydir = .false.
      zdir = .false.

      select case(dir)
        case(1)
          xdir = .true.
        case(2)
          ydir = .true.
        case(3)
          zdir = .true.
      end select
    else
      xdir = .true.
      ydir = .true.
      zdir = .true.
    end if

    !internal boundaries
    if (xdir) then
      call send_w
      call send_e

      call recv_e
      call recv_w
    end if

    if (ydir) then
      call send_s
      call send_n
      
      call recv_n
      call recv_s
    end if

    if (zdir) then
      call send_b
      call send_t

      call recv_t
      call recv_b
    end if

    !global domain boundaries
    if (xdir) then
      if (Btype(We)==BC_MPI_PERIODIC.or.Btype(Ea)==BC_MPI_PERIODIC) then
        if (iim==1) then
          call send(Phi, We)
        else if (iim==nxims) then
          call recv(Phi, Ea)
        end if     
        if (iim==nxims) then
          call send(Phi, Ea)
        else if (iim==1) then
          call recv(Phi, We)
        end if
      end if
    end if

    if (ydir) then  
      if (Btype(So)==BC_MPI_PERIODIC.or.Btype(No)==BC_MPI_PERIODIC) then
        if (jim==1) then
          call send(Phi, So)
        else if (jim==nyims) then
          call recv(Phi, No)
        end if
        if (jim==nyims) then
          call send(Phi, No)
        else if (jim==1) then
          call recv(Phi, So)
        end if
      end if
    end if

    if (zdir) then
      if (Btype(Bo)==BC_MPI_PERIODIC.or.Btype(To)==BC_MPI_PERIODIC) then
        if (kim==1) then
          call send(Phi, Bo)
        else if (kim==nzims) then
          call recv(Phi, To)
        end if
        if (kim==nzims) then
          call send(Phi, To)
        else if (kim==1) then
          call recv(Phi, Bo)
        end if
      end if
    end if
    
    call MPI_Waitall(size(requests), requests, MPI_STATUSES_IGNORE, ierr)
    
  contains
  
  
    subroutine send(a, side)
      real(knd), contiguous, intent(in) :: a(:,:,:)
      integer, intent(in) :: side

      call MPI_ISend(a, 1, send_mpi_types(side, component), neigh_ranks(side), &
                     1000+component, domain_comm, requests(side), ierr)
      if (ierr/=0) stop "error sending MPI message."    
    end subroutine
    
    subroutine recv(a, side)
      real(knd), intent(out) :: a(:,:,:)
      integer, intent(in) :: side

      call MPI_IRecv(a, 1, recv_mpi_types(side, component), neigh_ranks(side), &
                     1000+component, domain_comm, requests(side+6), ierr)
      if (ierr/=0) stop "error receiving MPI message."
    end subroutine
    

    subroutine send_w
      if (iim>1) then
        call send(Phi, We)
      end if
    end subroutine
    subroutine recv_w
      if (iim>1) then
        call recv(Phi, We)
      end if
    end subroutine
    subroutine send_e
      if (iim<nxims) then
        call send(Phi, Ea)
      end if
    end subroutine       
    subroutine recv_e
      if (iim<nxims) then
        call recv(Phi, Ea)
      end if
     end subroutine
    subroutine send_s
      if (jim>1) then
        call send(Phi, So)
      end if
    end subroutine
    subroutine recv_s
      if (jim>1) then
        call recv(Phi, So)
      end if
    end subroutine
    subroutine send_n
      if (jim<nyims) then
        call send(Phi, No)
      end if
    end subroutine
    subroutine recv_n
      if (jim<nyims) then
        call recv(Phi, No)
      end if
    end subroutine
    subroutine send_b
      if (kim>1) then
        call send(Phi, Bo)
      end if
    end subroutine
    subroutine recv_b
      if (kim>1) then
        call recv(Phi, Bo)
      end if
    end subroutine
    subroutine send_t
      if (kim<nzims) then
        call send(Phi, To)
      end if
    end subroutine
    subroutine recv_t
      if (kim<nzims) then
        call recv(Phi, To)
      end if
    end subroutine
    
  end subroutine par_exchange_boundaries
  
  
  subroutine par_exchange_UVW(U, V, W)
    use Parameters, only: Btype
    real(knd), intent(inout), contiguous :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    call par_exchange_boundaries(U, Btype, 1)
    call par_exchange_boundaries(V, Btype, 2)
    call par_exchange_boundaries(W, Btype, 3)
  end subroutine

  subroutine par_exchange_U_x(U, component)
    use Parameters, only: Btype
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
    integer, intent(in) :: component
    call par_exchange_boundaries(U, Btype, component, dir=1)
  end subroutine
  
  subroutine par_exchange_U_y(U, component)
    use Parameters, only: Btype
    integer, intent(in) :: component
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
    call par_exchange_boundaries(U, Btype, component, dir=2)
  end subroutine
  
  subroutine par_exchange_U_z(U, component)
    use Parameters, only: Btype
    integer, intent(in) :: component
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
    call par_exchange_boundaries(U, Btype, component, dir=3)
  end subroutine
  
  subroutine par_exchange_Sc_x(U, SBtype)
    use Parameters, only: Prnx, Prny, Prnz
    real(knd), intent(inout) :: U(-1:,-1:,-1:)
    integer, intent(in) :: SBtype(6)
    call par_exchange_boundaries(U, SBtype, 4, dir=1)
  end subroutine
  
  subroutine par_exchange_Sc_y(U, SBType)
    use Parameters, only: Prnx, Prny, Prnz
    real(knd), intent(inout) :: U(-1:,-1:,-1:)
    integer, intent(in) :: SBtype(6)
    call par_exchange_boundaries(U, SBtype, 4, dir=2)
  end subroutine
  
  subroutine par_exchange_Sc_z(U, SBType)
    use Parameters, only: Prnx, Prny, Prnz
    real(knd), intent(inout) :: U(-1:,-1:,-1:)
    integer, intent(in) :: SBtype(6)
    call par_exchange_boundaries(U, SBtype, 4, dir=3)
  end subroutine
   
 
  
  
  
  
  subroutine par_exchange_Pr(Phi)
    use Parameters, only: We, Ea, So, No, Bo, To, BC_MPI_PERIODIC, Prnx, Prny, Prnz, Btype
    real(knd), intent(inout), contiguous :: Phi(0:,0:,0:)
    integer :: ierr
    integer :: nx, ny, nz
    integer :: requests(12)
    
    requests = MPI_REQUEST_NULL

    !internal boundaries
    call send_w
    call send_e
    call send_s
    call send_n
    call send_b
    call send_t

    call recv_w
    call recv_e
    call recv_s
    call recv_n
    call recv_b
    call recv_t
   

    !global domain boundaries
    if (Btype(We)==BC_MPI_PERIODIC.or.Btype(Ea)==BC_MPI_PERIODIC) then
      if (iim==1) then
        call send(Phi, We)
      else if (iim==nxims) then
        call recv(Phi, Ea)
      end if     
    end if

    if (Btype(So)==BC_MPI_PERIODIC.or.Btype(No)==BC_MPI_PERIODIC) then
      if (jim==1) then
        call send(Phi, So)
      else if (jim==nyims) then
        call recv(Phi, No)
      end if
    end if
          
    if (Btype(Bo)==BC_MPI_PERIODIC.or.Btype(To)==BC_MPI_PERIODIC) then
      if (kim==1) then
        call send(Phi, Bo)
      else if (kim==nzims) then
        call recv(Phi, To)
      end if
    end if

    call MPI_Waitall(size(requests), requests, MPI_STATUSES_IGNORE, ierr)    
    
  contains
  
  
    subroutine send(a, side)
      real(knd), contiguous, intent(in) :: a(:,:,:)
      integer, intent(in) :: side

      call MPI_ISend(a, 1, send_mpi_types(side, 6), neigh_ranks(side), &
                     1000+8, domain_comm, requests(side), ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    
    subroutine recv(a, side)
      real(knd), intent(out) :: a(:,:,:)
      integer, intent(in) :: side

      call MPI_IRecv(a, 1, recv_mpi_types(side, 6), neigh_ranks(side), &
                     1000+8, domain_comm, requests(side+6), ierr)
      if (ierr/=0) stop "error receiving MPI message."
    end subroutine
    

    subroutine send_w
      if (iim>1) then
        call send(Phi, We)
      end if
    end subroutine
     
    subroutine recv_w
      if (iim>1) then
        call recv(Phi, We)
      end if
    end subroutine
     
    subroutine send_e
      if (iim<nxims) then
        call send(Phi, Ea)
      end if
    end subroutine
    
    subroutine recv_e
      if (iim<nxims) then
        call recv(Phi, Ea)
      end if
    end subroutine
    
    subroutine send_s
      if (jim>1) then
        call send(Phi, So)
      end if
    end subroutine

    subroutine recv_s
      if (jim>1) then
        call recv(Phi, So)
      end if
    end subroutine

    subroutine send_n
      if (jim<nyims) then
        call send(Phi, No)
      end if
    end subroutine
    
    subroutine recv_n
      if (jim<nyims) then
        call recv(Phi, No)
      end if
    end subroutine
    
    subroutine send_b
      if (kim>1) then
        call send(Phi, Bo)
      end if
    end subroutine
    
    subroutine recv_b
      if (kim>1) then
        call recv(Phi, Bo)
      end if
    end subroutine
    
    subroutine send_t
      if (kim<nzims) then
        call send(Phi, To)
      end if
    end subroutine
    
    subroutine recv_t
      if (kim<nzims) then
        call recv(Phi, To)
      end if
    end subroutine
    
  end subroutine par_exchange_Pr
  
 

 
  
  
  subroutine par_exchange_Q(Phi)
    use Parameters, only: We, Ea, So, No, Bo, To, BC_MPI_PERIODIC, &
                          nx=>Prnx, ny=>Prny, nz=>Prnz, Btype
    real(knd), intent(inout), contiguous :: Phi(0:,0:,0:)
    integer :: ierr
    integer :: requests(12)
    real(knd) :: tmp_w(ny,nz), tmp_e(ny,nz), &
                 tmp_s(nx+2,nz), tmp_n(nx+2,nz), &
                 tmp_b(nx+2,ny+2), tmp_t(nx+2,ny+2)
    logical :: update(6)
    
    requests = MPI_REQUEST_NULL

    ierr = 0
    
    update = .false.
    
    !internal boundaries
    call send_e
    call send_w
    call send_n
    call send_s
    call send_t
    call send_b
    
    call recv_w
    call recv_e
    call recv_s
    call recv_n
    call recv_b
    call recv_t

    !global domain boundaries
    if (Btype(We)==BC_MPI_PERIODIC.or.Btype(Ea)==BC_MPI_PERIODIC) then
      if (iim==1) then
        call recv(tmp_w, We)
      else if (iim==nxims) then
        call send(Phi, Ea)
      end if     
      if (iim==nxims) then
        call recv(tmp_e, Ea)
      else if (iim==1) then
        call send(Phi, We)
      end if
    end if

    if (Btype(So)==BC_MPI_PERIODIC.or.Btype(No)==BC_MPI_PERIODIC) then
      if (jim==1) then
        call recv(tmp_s, So)
      else if (jim==nyims) then
        call send(Phi, No)
      end if
      if (jim==nyims) then
        call recv(tmp_n, No)
      else if (jim==1) then
        call send(Phi, So)
      end if
    end if

    if (Btype(Bo)==BC_MPI_PERIODIC.or.Btype(To)==BC_MPI_PERIODIC) then
      if (kim==1) then
        call recv(tmp_b, Bo)
      else if (kim==nzims) then
        call send(Phi, To)
      end if
      if (kim==nzims) then
        call recv(tmp_t, To)
      else if (kim==1) then
        call send(Phi, Bo)
      end if
    end if
    
    call MPI_Waitall(6, requests(7:12), MPI_STATUSES_IGNORE, ierr)
    
    if (update(We)) Phi(1,1:ny,1:nz) = Phi(1,1:ny,1:nz) + tmp_w
    if (update(Ea)) Phi(nx,1:ny,1:nz) = Phi(nx,1:ny,1:nz) + tmp_e
    if (update(So)) Phi(0:nx+1,1,1:nz) = Phi(0:nx+1,1,1:nz) + tmp_s
    if (update(No)) Phi(0:nx+1,ny,1:nz) = Phi(0:nx+1,ny,1:nz) + tmp_n
    if (update(Bo)) Phi(0:nx+1,0:ny+1,1) = Phi(0:nx+1,0:ny+1,1) + tmp_b
    if (update(To)) Phi(0:nx+1,0:ny+1,nz) = Phi(0:nx+1,0:ny+1,nz) + tmp_t
    
    call MPI_Waitall(6, requests(1:6), MPI_STATUSES_IGNORE, ierr)

  contains
  
  
    subroutine send(a, side)
      real(knd), contiguous, intent(in) :: a(:,:,:)
      integer, intent(in) :: side

      call MPI_ISend(a, 1, send_mpi_types(side, 7), neigh_ranks(side), &
                     1000+7 , domain_comm, requests(side), ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    
    subroutine recv(a, side)
      real(knd), contiguous, intent(out) :: a(:,:)
      integer, intent(in) :: side

      call MPI_IRecv(a, size(a), MPI_KND, neigh_ranks(side), &
                     1000+7 , domain_comm, requests(side+6), ierr)
      if (ierr/=0) stop "error receiving MPI message."
      
      update(side) = .true.
    end subroutine
    

    subroutine recv_w
      if (iim>1) then
        call recv(tmp_w, We)
      end if
    end subroutine
    subroutine send_w
      if (iim>1) then
        call send(Phi, We)
      end if
    end subroutine
    subroutine recv_e
      if (iim<nxims) then
        call recv(tmp_e, Ea)
      end if
    end subroutine       
    subroutine send_e
      if (iim<nxims) then
        call send(Phi, Ea)
      end if
     end subroutine
    subroutine recv_s
      if (jim>1) then
        call recv(tmp_s, So)
      end if
    end subroutine
    subroutine send_s
      if (jim>1) then
        call send(Phi, So)
      end if
    end subroutine
    subroutine recv_n
      if (jim<nyims) then
        call recv(tmp_n, No)
      end if
    end subroutine
    subroutine send_n
      if (jim<nyims) then
        call send(Phi, No)
      end if
    end subroutine
    subroutine recv_b
      if (kim>1) then
        call recv(tmp_b, Bo)
      end if
    end subroutine
    subroutine send_b
      if (kim>1) then
        call send(Phi, Bo)
      end if
    end subroutine
    subroutine recv_t
      if (kim<nzims) then
        call recv(tmp_t, To)
      end if
    end subroutine
    subroutine send_t
      if (kim<nzims) then
        call send(Phi, To)
      end if
    end subroutine
    
  end subroutine par_exchange_Q
  
  
  
  
  
  
  

  
  
  subroutine par_exchange_boundaries_xz(Phi, nx, nz, Btype, lbx, lbz, widthx, widthz)
    use Parameters, only: We, Ea, Bo, To, BC_MPI_PERIODIC
    real(knd), intent(inout), contiguous :: Phi(lbx:, lbz:)
    integer, intent(in) :: nx, nz
    integer, intent(in) :: Btype(6)
    integer, intent(in) :: lbx, lbz, widthx, widthz
    logical :: oddx, oddz, evenx, evenz
    integer :: ierr, tag
    logical :: xdir, zdir
    
    ierr = 0; tag = 1500
    
    oddx = mod(iim,2) == 1
    evenx = .not. oddx
    
    oddz = mod(kim,2) == 1
    evenz = .not. oddz

    !internal boundaries

    tag = tag + 1
    if (oddx) then
      call send_w
    else
      call recv_e
    end if
    if (evenx) then
      call send_w
    else
      call recv_e
    end if

    tag = tag + 1
    if (oddx) then
      call send_e
    else
      call recv_w
    end if
    if (evenx) then
      call send_w
    else
      call recv_w
    end if


    tag = tag + 1
    if (oddz) then
      call send_b
    else
      call recv_t
    end if
    if (evenz) then
      call send_b
    else
      call recv_t
    end if

    tag = tag + 1
    if (oddz) then
      call send_t
    else
      call recv_b
    end if
    if (evenz) then
      call send_t
    else
      call recv_b
    end if

    !global domain boundaries

    tag = tag + 1
    if (Btype(We)==BC_MPI_PERIODIC.or.Btype(Ea)==BC_MPI_PERIODIC) then
      if (iim==1) then
        call send(Phi(1:0+widthx,1:nz), w_rank)
      else if (iim==nxims) then
        call recv(Phi(nx+1:nx+widthx,1:nz), e_rank)
      end if
      if (iim==nxims) then
        call send(Phi(nx+1-widthx:nx,1:nz), e_rank)
      else if (iim==1) then
        call recv(Phi(1-widthx:0,1:nz), w_rank)
      end if
    end if



    tag = tag + 1
    if (Btype(Bo)==BC_MPI_PERIODIC.or.Btype(To)==BC_MPI_PERIODIC) then
      if (kim==1) then
        call send(Phi(1-widthx:nx+widthx,1:0+widthz), b_rank)
      else if (kim==nzims) then
        call recv(Phi(1-widthx:nx+widthx,nz+1:nz+widthz), t_rank)
      end if
      if (kim==nzims) then
        call send(Phi(1-widthx:nx+widthx,nz+1-widthz:nz), t_rank)
      else if (kim==1) then
        call recv(Phi(1-widthx:nx+widthx,1-widthz:0), b_rank)
      end if
    end if

    
  contains
  
  
    subroutine send(a,to)
      real(knd), intent(in) :: a(:,:)
      integer, intent(in) :: to
      !Must be domain_comm, for `to` and `from` derived from `x_rank` to be valid!
      call MPI_Send(a, size(a) , MPI_KND, to, tag, domain_comm, ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    
    subroutine recv(a,from)
      real(knd), intent(out) :: a(:,:)
      integer, intent(in) :: from

      call MPI_Recv(a, size(a) , MPI_KND, from, tag, domain_comm, MPI_STATUS_IGNORE, ierr)
      if (ierr/=0) stop "error receiving MPI message."
    end subroutine
    

    subroutine send_w
      if (iim>1) then
        call send(Phi(1:0+widthx,1:nz), w_rank)
      end if
    end subroutine
    subroutine recv_w
      if (iim>1) then
        call recv(Phi(1-widthx:0,1:nz), w_rank)
      end if
    end subroutine
    subroutine send_e
      if (iim<nxims) then
        call send(Phi(nx+1-widthx:nx,1:nz), e_rank)
      end if
    end subroutine
    subroutine recv_e
      if (iim<nxims) then
        call recv(Phi(nx+1:nx+widthx,1:nz), e_rank)
      end if
    end subroutine
    subroutine send_b
      if (kim>1) then
        call send(Phi(1-widthx:nx+widthx,1:0+widthz), b_rank)
      end if
    end subroutine
    subroutine recv_b
      if (kim>1) then
        call recv(Phi(1-widthx:nx+widthx,1-widthz:0), b_rank)
      end if
    end subroutine
    subroutine send_t
      if (kim<nzims) then
        call send(Phi(1-widthx:nx+widthx,nz+1-widthz:nz), t_rank)
      end if
    end subroutine
    subroutine recv_t
      if (kim<nzims) then
        call recv(Phi(1-widthx:nx+widthx,nz+1:nz+widthz), t_rank)
      end if
    end subroutine
    
  end subroutine par_exchange_boundaries_xz
  
  
  subroutine par_exchange_boundaries_yz(Phi, ny, nz, Btype, lby, lbz, widthy, widthz)
    use Parameters, only: So, No, Bo, To, BC_MPI_PERIODIC
    real(knd), intent(inout), contiguous :: Phi(lby:, lbz:)
    integer, intent(in) :: ny, nz
    integer, intent(in) :: Btype(6)
    integer, intent(in) :: lby, lbz, widthy, widthz
    logical :: oddy, oddz, eveny, evenz
    integer :: ierr, tag
    logical :: ydir, zdir
    
    ierr = 0; tag = 1600
    
    oddy = mod(jim,2) == 1
    eveny = .not. oddy
    
    oddz = mod(kim,2) == 1
    evenz = .not. oddz

    !internal boundaries

    tag = tag + 1
    if (oddy) then
      call send_s
    else
      call recv_n
    end if
    if (eveny) then
      call send_s
    else
      call recv_n
    end if

    tag = tag + 1
    if (oddy) then
      call send_n
    else
      call recv_s
    end if
    if (eveny) then
      call send_n
    else
      call recv_s
    end if


    tag = tag + 1
    if (oddz) then
      call send_b
    else
      call recv_t
    end if
    if (evenz) then
      call send_b
    else
      call recv_t
    end if

    tag = tag + 1
    if (oddz) then
      call send_t
    else
      call recv_b
    end if
    if (evenz) then
      call send_t
    else
      call recv_b
    end if

    !global domain boundaries

    tag = tag + 1
    if (Btype(So)==BC_MPI_PERIODIC.or.Btype(No)==BC_MPI_PERIODIC) then
      if (jim==1) then
        call send(Phi(1:0+widthy,1:nz), s_rank)
      else if (jim==nyims) then
        call recv(Phi(ny+1:ny+widthy,1:nz), n_rank)
      end if
      if (jim==nyims) then
        call send(Phi(ny+1-widthy:ny,1:nz), n_rank)
      else if (jim==1) then
        call recv(Phi(1-widthy:0,1:nz), s_rank)
      end if
    end if



    tag = tag + 1
    if (Btype(Bo)==BC_MPI_PERIODIC.or.Btype(To)==BC_MPI_PERIODIC) then
      if (kim==1) then
        call send(Phi(1-widthy:ny+widthy,1:0+widthz), b_rank)
      else if (kim==nzims) then
        call recv(Phi(1-widthy:ny+widthy,nz+1:nz+widthz), t_rank)
      end if
      if (kim==nzims) then
        call send(Phi(1-widthy:ny+widthy,nz+1-widthz:nz), t_rank)
      else if (kim==1) then
        call recv(Phi(1-widthy:ny+widthy,1-widthz:0), b_rank)
      end if
    end if

    
  contains
  
  
    subroutine send(a,to)
      real(knd), intent(in) :: a(:,:)
      integer, intent(in) :: to
      !Must be domain_comm, for `to` and `from` derived from `x_rank` to be valid!
      call MPI_Send(a, size(a) , MPI_KND, to, tag, domain_comm, ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    
    subroutine recv(a,from)
      real(knd), intent(out) :: a(:,:)
      integer, intent(in) :: from

      call MPI_Recv(a, size(a) , MPI_KND, from, tag, domain_comm, MPI_STATUS_IGNORE, ierr)
      if (ierr/=0) stop "error receiving MPI message."
    end subroutine
    

    subroutine send_s
      if (jim>1) then
        call send(Phi(1:0+widthy,1:nz), s_rank)
      end if
    end subroutine
    subroutine recv_s
      if (jim>1) then
        call recv(Phi(1-widthy:0,1:nz), s_rank)
      end if
    end subroutine
    subroutine send_n
      if (jim<nyims) then
        call send(Phi(ny+1-widthy:ny,1:nz), n_rank)
      end if
    end subroutine
    subroutine recv_n
      if (jim<nyims) then
        call recv(Phi(ny+1:ny+widthy,1:nz), n_rank)
      end if
    end subroutine
    subroutine send_b
      if (kim>1) then
        call send(Phi(1-widthy:ny+widthy,1:0+widthz), b_rank)
      end if
    end subroutine
    subroutine recv_b
      if (kim>1) then
        call recv(Phi(1-widthy:ny+widthy,1-widthz:0), b_rank)
      end if
    end subroutine
    subroutine send_t
      if (kim<nzims) then
        call send(Phi(1-widthy:ny+widthy,nz+1-widthz:nz), t_rank)
      end if
    end subroutine
    subroutine recv_t
      if (kim<nzims) then
        call recv(Phi(1-widthy:ny+widthy,nz+1:nz+widthz), t_rank)
      end if
    end subroutine
    
  end subroutine par_exchange_boundaries_yz
  
  
end module exchange_par
