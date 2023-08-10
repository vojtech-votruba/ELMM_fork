module Boundaries
  use Parameters
#ifdef PAR
  use exchange_par
#endif

implicit none


  private
  public volumePr, &
         GridCoords, GridCoords_interp, GridCoords_interp_U, GridCoords_interp_V, GridCoords_interp_W, &
         InDomain,&
         BoundUVW, BoundU, Bound_Phi, Bound_Pr, Bound_Q,&
         ShearInlet, ParabolicInlet, ConstantInlet

  interface GridCoords
    module procedure GridCoords_scalar
    module procedure GridCoords_vector
  end interface

  interface InDomain
    module procedure InDomain_3r
    module procedure InDomain_vec
  end interface
 
  interface InGlobalDomain
    module procedure InGlobalDomain_3r
    module procedure InGlobalDomain_vec
  end interface
 
  interface volumePr
    module procedure volumePr_arr
    module procedure volumePr_int
  end interface

 contains

  pure function volumePr_arr(indexes) result(res)
    real(knd) :: res
    integer, intent(in) :: indexes(3)
    res = dxPr(indexes(1))*dyPr(indexes(2))*dzPr(indexes(3))
  end function

  elemental function volumePr_int(i,j,k) result(res)
    real(knd) :: res
    integer, intent(in) :: i,j,k
    res = dxPr(i)*dyPr(j)*dzPr(k)
  end function

  !TODO: one common base subroutine just called with different arguments
  !consider binary search for general case (branch mispredictions vs. more iterations)
  elemental subroutine GridCoords_scalar(xi, yj, zk, x, y, z)
    integer, intent(out):: xi, yj, zk
    real(knd), intent(in):: x, y, z
    integer :: i

    
    if (gridtype==GRID_UNIFORM) then

        xi = min( max(nint( (x - xU(0))/dxmin + 0.5_knd ),1) , Prnx)
        yj = min( max(nint( (y - yV(0))/dymin + 0.5_knd ),1) , Prny)
        zk = min( max(nint( (z - zW(0))/dzmin + 0.5_knd ),1) , Prnz)

    else

      xi = Prnx
      do i = 1, Prnx
       if (xU(i)>=x) then
                     xi = i
                     exit
                    end if
      end do

      yj = Prny
      do i = 1, Prny
       if (yV(i)>=y) then
                     yj = i
                     exit
                    end if
      end do
      zk = Prnz
      do i = 1, Prnz
       if (zW(i)>=z) then
                     zk = i
                     exit
                    end if
      end do

    end if
  end subroutine GridCoords_scalar

  subroutine GridCoords_vector(ri, r)
    integer, intent(out):: ri(3)
    real(knd), intent(in):: r(3)
    
    call GridCoords_scalar(ri(1), ri(2), ri(3), r(1), r(2), r(3))
  end subroutine



  !finds the index from which interpolate by bilinear interpolation
  !the point (x,y,z) should lie between [xi,yj,zk] and [xi+1,yj+1,zk+1]
  elemental subroutine GridCoords_interp(xi, yj, zk, x, y, z)
    integer, intent(out):: xi, yj, zk
    real(knd), intent(in):: x, y, z
    integer :: i
    
    xi = min( max(floor( (x - xU(0))/dxmin + 0.5_knd ),0) , Prnx)
    yj = min( max(floor( (y - yV(0))/dymin + 0.5_knd ),0) , Prny)
    zk = min( max(floor( (z - zW(0))/dzmin + 0.5_knd ),0) , Prnz)
  end subroutine

  !finds the index from which interpolate by bilinear interpolation
  !the point (x,y,z) should lie between [xi,yj,zk] and [xi+1,yj+1,zk+1]
  elemental subroutine GridCoords_interp_U(xi, yj, zk, x, y, z)
    integer, intent(out) :: xi, yj, zk
    real(knd), intent(in) :: x, y, z
    integer :: i

     xi = min( max(floor( (x - xU(0))/dxmin ),0) , Prnx)
     yj = min( max(floor( (y - yV(0))/dymin + 0.5_knd ),0) , Prny)
     zk = min( max(floor( (z - zW(0))/dzmin + 0.5_knd ),0) , Prnz)
  end subroutine


  !finds the index from which interpolate by bilinear interpolation
  !the point (x,y,z) should lie between [xi,yj,zk] and [xi+1,yj+1,zk+1]
  elemental subroutine GridCoords_interp_V(xi, yj, zk, x, y, z)
    integer, intent(out) :: xi, yj, zk
    real(knd), intent(in) :: x, y, z
    integer :: i

    xi = min( max(floor( (x - xU(0))/dxmin + 0.5_knd ),0) , Prnx)
    yj = min( max(floor( (y - yV(0))/dymin ),0), Prny)
    zk = min( max(floor( (z - zW(0))/dzmin + 0.5_knd ),0) , Prnz)
  end subroutine


  !finds the index from which interpolate by bilinear interpolation
  !the point (x,y,z) should lie between [xi,yj,zk] and [xi+1,yj+1,zk+1]
  elemental subroutine GridCoords_interp_W(xi, yj, zk, x, y, z)
    integer, intent(out) :: xi, yj, zk
    real(knd), intent(in) :: x, y, z
    integer :: i

    xi = min( max(floor( (x - xU(0))/dxmin + 0.5_knd ),0) , Prnx)
    yj = min( max(floor( (y - yV(0))/dymin + 0.5_knd ),0) , Prny)
    zk = min( max(floor( (z - zW(0))/dzmin ),0) , Prnz)
  end subroutine




  logical function InDomain_3r(x, y, z) result(res)
    real(knd), intent(in) :: x, y, z
    res = (x>=xU(0) .and. x<=xU(Prnx)) .and. &
          (y>=yV(0) .and. y<=yV(Prny)) .and. &
          (z>=zW(0) .and. z<=zW(Prnz))
  end function

  logical function InDomain_vec(r)
    real(knd), intent(in) :: r(3)
    InDomain_vec = InDomain_3r(r(1),r(2),r(3))
  end function


  logical function InGlobalDomain_3r(x, y, z) result(res)
    real(knd), intent(in) :: x, y, z
    res = (x>=gxmin .and. x<=gxmax) .and. &
          (y>=gymin .and. y<=gymax) .and. &
          (z>=gzmin .and. z<=gzmax)
  end function

  logical function InGlobalDomain_vec(r)
    real(knd), intent(in) :: r(3)
    InGlobalDomain_vec = InGlobalDomain_3r(r(1),r(2),r(3))
  end function


  recursive subroutine BoundU(component, U, Uin, regime)
    integer, intent(in)                   :: component
    real(knd), contiguous, intent(inout)  :: U(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)     :: Uin(-2:,-2:)
    integer, optional, intent(in)         :: regime
    integer :: reg, i, j, k, nx, ny, nz
    !intermediate regime is when doing BC for one of the terms for dU/dt
    !currently used when explicit filtering is on only
    integer, parameter :: intermediate = 2

    if (present(regime)) then
      reg = regime
    else
      reg = 0
    end if

    if (component==1) then
      nx = Unx
      ny = Uny
      nz = Unz
    else if (component==2) then
      nx = Vnx
      ny = Vny
      nz = Vnz
    else
      nx = Wnx
      ny = Wny
      nz = Wnz
    end if

    !$omp parallel private(i,j,k)

    !!! corners and edges for periodic conditions
    if (Btype(Ea)==BC_PERIODIC.and.Btype(No)==BC_PERIODIC.and.Btype(To)==BC_PERIODIC) then
    !$omp sections
    !$omp section
     do k = nz+1, nz+3
      do j = ny+1, ny+3
       do i = nx+1, nx+3
        U(i,j,k) = U(i-nx,j-ny,k-nz)
       end do
      end do
     end do
    !$omp section
     do k = -2, 0
      do j = ny+1, ny+3
       do i = nx+1, nx+3
        U(i,j,k) = U(i-nx,j-ny,k+nz)
       end do
      end do
     end do
    !$omp section
     do k = nz+1, nz+3
      do j = -2, 0
       do i = nx+1, nx+3
        U(i,j,k) = U(i-nx,j+ny,k-nz)
       end do
      end do
     end do
    !$omp section
     do k = nz+1, nz+3
      do j = ny+1, ny+3
       do i = -2, 0
        U(i,j,k) = U(i+nx,j-ny,k-nz)
       end do
      end do
     end do
    !$omp section
     do k = -2, 0
      do j = -2, 0
       do i = nx+1, nx+3
        U(i,j,k) = U(i-nx,j+ny,k+nz)
       end do
      end do
     end do
    !$omp section
     do k = nz+1, nz+3
      do j = -2, 0
       do i = -2, 0
        U(i,j,k) = U(i+nx,j+ny,k-nz)
       end do
      end do
     end do
    !$omp section
     do k = -2, 0
      do j = ny+1, ny+3
       do i = -2, 0
        U(i,j,k) = U(i+nx,j-ny,k+nz)
       end do
      end do
     end do
    !$omp section
     do k = -2, 0
      do j = -2, 0
       do i = -2, 0
        U(i,j,k) = U(i+nx,j+ny,k+nz)
       end do
      end do
     end do
     !$omp end sections
    end if

    if (Btype(Ea)==BC_PERIODIC.and.Btype(No)==BC_PERIODIC) then
     !$omp sections
     !$omp section
     do k = 1, nz
      do j = ny+1, ny+3
       do i = nx+1, nx+3
        U(i,j,k) = U(i-nx,j-ny,k)
       end do
      end do
     end do
     !$omp section
     do k = 1, nz
      do j = -2, 0
       do i = nx+1, nx+3
        U(i,j,k) = U(i-nx,j+ny,k)
       end do
      end do
     end do
     !$omp section
     do k = 1, nz
      do j = ny+1, ny+3
       do i = -2, 0
        U(i,j,k) = U(i+nx,j-ny,k)
       end do
      end do
     end do
     !$omp section
     do k = 1, nz
      do j = -2, 0
       do i = -2, 0
        U(i,j,k) = U(i+nx,j+ny,k)
       end do
      end do
     end do
     !$omp end sections
    end if


    if (Btype(Ea)==BC_PERIODIC.and.Btype(To)==BC_PERIODIC) then
     !$omp sections
     !$omp section
     do k = nz+1, nz+3
      do j = 1, ny
       do i = nx+1, nx+3
        U(i,j,k) = U(i-nx,j,k-nz)
       end do
      end do
     end do
     !$omp section
     do k = -2, 0
      do j = 1, ny
       do i = nx+1, nx+3
        U(i,j,k) = U(i-nx,j,k+nz)
       end do
      end do
     end do
     !$omp section
     do k = nz+1, nz+3
      do j = 1, ny
       do i = -2, 0
        U(i,j,k) = U(i+nx,j,k-nz)
       end do
      end do
     end do
     !$omp section
     do k = -2, 0
      do j = 1, ny
       do i = -2, 0
        U(i,j,k) = U(i+nx,j,k+nz)
       end do
      end do
     end do
     !$omp end sections
    end if


    if (Btype(No)==BC_PERIODIC.and.Btype(To)==BC_PERIODIC) then
     !$omp sections
     !$omp section
     do k = nz+1, nz+3
      do j = ny+1, ny+3
       do i = 1, nx
        U(i,j,k) = U(i,j-ny,k-nz)
       end do
      end do
     end do
     !$omp section
     do k = -2, 0
      do j = ny+1, ny+3
       do i = 1, nx
        U(i,j,k) = U(i,j-ny,k+nz)
       end do
      end do
     end do
     do k = nz+1, nz+3
      do j = -2, 0
       do i = 1, nx
        U(i,j,k) = U(i,j+ny,k-nz)
       end do
      end do
     end do
     !$omp section
     do k = -2, 0
      do j = -2, 0
       do i = 1, nx
        U(i,j,k) = U(i,j+ny,k+nz)
       end do
      end do
     end do
     !$omp end sections
    end if 

    !$omp end parallel
    
    
#ifdef PAR
    call par_exchange_U_x(U, component)
#endif

    !$omp parallel sections
    !$omp section
    if (Btype(We)==BC_DIRICHLET.and.reg/=intermediate) then
      if (component==1) then
         do k = 1, nz
          do j = 1, ny                       !Dirichlet inlet
           U(0,j,k) = Uin(j,k)
           U(-1,j,k) = Uin(j,k) + (Uin(j,k)-U(1,j,k))
           U(-2,j,k) = Uin(j,k) + (Uin(j,k)-U(2,j,k))
          end do
         end do
      else
         do k = 1, nz
          do j = 1, ny                       !Dirichlet inlet
           U(0,j,k) = -U(1,j,k)
           U(-1,j,k) = -U(2,j,k)
           U(-2,j,k) = -U(3,j,k)
          end do
         end do
      end if
    else if (Btype(We)==BC_NOSLIP .or. (component==1.and.Btype(We)==BC_FREESLIP) .or. &
             (Btype(We)==BC_DIRICHLET.and.reg==intermediate)) then
      if (component==1) then
         do k = 1, nz
          do j = 1, ny                       !Solid wall
           U(0,j,k) = 0
           U(-1,j,k) = -U(1,j,k)
           U(-2,j,k) = -U(2,j,k)
          end do
         end do
      else
         do k = 1, nz
          do j = 1, ny                       !Solid wall
           U(0,j,k) = -U(1,j,k)
           U(-1,j,k) = -U(2,j,k)
           U(-2,j,k) = -U(3,j,k)
          end do
         end do
      end if
    else if (Btype(We)==BC_NEUMANN.or.(component/=1.and.Btype(We)==BC_FREESLIP)) then
      if (component==1) then
         do k = 1, nz
          do j = 1, ny                       !Neumann inlet
           U(0,j,k) = U(1,j,k)
           U(-1,j,k) = U(1,j,k)
           U(-2,j,k) = U(1,j,k)
          end do
         end do
      else
         do k = 1, nz
          do j = 1, ny                       !Neumann inlet
           U(0,j,k) = U(1,j,k)
           U(-1,j,k) = U(1,j,k)
           U(-2,j,k) = U(1,j,k)
          end do
         end do
      end if
    else if (Btype(We)==BC_PERIODIC) then  !Periodic BC
      do k = 1, nz
       do j = 1, ny
        U(0,j,k) = U(nx,j,k)
        U(-1,j,k) = U(nx-1,j,k)
        U(-2,j,k) = U(nx-2,j,k)
       end do
      end do
    else if (Btype(We)==BC_TURBULENT_INLET.or.Btype(We)==BC_INLET_FROM_FILE) then
      if (reg/=intermediate) then
        if (component==1) then
          do k = -1, nz+2
           do j = -1, ny+2
            U(0,j,k) = Uin(j,k)
            U(-1,j,k) = Uin(j,k) + (Uin(j,k)-U(1,j,k))
            U(-2,j,k) = Uin(j,k) + (Uin(j,k)-U(2,j,k))
           end do
          end do
        else
          do k = -1, nz+2
           do j = -1, ny+2
            U(0,j,k) = Uin(j,k) + (Uin(j,k)-U(1,j,k))
            U(-1,j,k) = Uin(j,k) + (Uin(j,k)-U(2,j,k))
            U(-2,j,k) = Uin(j,k) + (Uin(j,k)-U(3,j,k))
           end do
          end do
        end if
      else
        if (component==1) then
           do k = -1, nz+2
            do j = -1, ny+2
             U(0,j,k) = 0
             U(-1,j,k) = -U(1,j,k)
             U(-2,j,k) = -U(2,j,k)
            end do
           end do
        else
           do k = -1, nz+2
            do j = -1, ny+2
             U(0,j,k) = -U(1,j,k)
             U(-1,j,k) = -U(2,j,k)
             U(-2,j,k) = -U(3,j,k)
            end do
           end do
        end if
      end if
    end if


    !$omp section
    if (Btype(Ea)==BC_DIRICHLET.and.reg/=intermediate) then
      if (component==1) then
        do k = 1, nz
         do j = 1, ny                       !Dirichlet inlet
          U(nx+1,j,k) = Uin(j,k)
          U(nx+2,j,k) = Uin(j,k) + (Uin(j,k)-U(nx,j,k))
          U(nx+3,j,k) = Uin(j,k) + (Uin(j,k)-U(nx-1,j,k))
         end do
        end do
      else
        do k = 1, nz
         do j = 1, ny                       !Dirichlet inlet
          U(nx+1,j,k) = -U(nx,j,k)
          U(nx+2,j,k) = -U(nx-1,j,k)
          U(nx+3,j,k) = -U(nx-2,j,k)
         end do
        end do
      end if
    else if (Btype(Ea)==BC_NOSLIP .or. (component==1.and.Btype(Ea)==BC_FREESLIP) .or. &
              (Btype(Ea)==BC_DIRICHLET.and.reg==intermediate)) then
      if (component==1) then
        do k = 1, nz
         do j = 1, ny                       !Solid wall
          U(nx+1,j,k) = 0
          U(nx+2,j,k) = -U(nx,j,k)
          U(nx+3,j,k) = -U(nx-1,j,k)
         end do
        end do
      else
        do k = 1, nz
         do j = 1, ny                       !Solid wall
          U(nx+1,j,k) = -U(nx,j,k)
          U(nx+2,j,k) = -U(nx-1,j,k)
          U(nx+3,j,k) = -U(nx-2,j,k)
         end do
        end do
      end if
    else if (Btype(Ea)==BC_NEUMANN.or.(component/=1.and.Btype(Ea)==BC_FREESLIP)) then   !Neumann outlet
      do k = 1, nz
       do j = 1, ny
        U(nx+1,j,k) = U(nx,j,k)
        U(nx+2,j,k) = U(nx,j,k)
        U(nx+3,j,k) = U(nx,j,k)
       end do
      end do
    else if (Btype(Ea)==BC_PERIODIC) then  !Periodic BC
      do k = 1, nz
       do j = 1, ny
        U(nx+1,j,k) = U(1,j,k)
        U(nx+2,j,k) = U(2,j,k)
        U(nx+3,j,k) = U(3,j,k)
       end do
      end do
    else if (Btype(Ea)==BC_TURBULENT_INLET) then
      if (reg/=intermediate) then
        if (component==1) then
          do k = 1, nz
           do j = 1, ny
            U(nx+1,j,k) = Uin(j,k)
            U(nx+2,j,k) = Uin(j,k) + (Uin(j,k)-U(nx,j,k))
            U(nx+3,j,k) = Uin(j,k) + (Uin(j,k)-U(nx-1,j,k))
           end do
          end do
        else
          do k = 1, nz
           do j = 1, ny
            U(nx+1,j,k) = Uin(j,k) + (Uin(j,k)-U(nx,j,k))
            U(nx+2,j,k) = Uin(j,k) + (Uin(j,k)-U(nx-1,j,k))
            U(nx+3,j,k) = Uin(j,k) + (Uin(j,k)-U(nx-2,j,k))
           end do
          end do
        end if
      else
        if (component==1) then
           do k = 1, nz
            do j = 1, ny
             U(nx+1,j,k) = 0
             U(nx+2,j,k) = -U(nx,j,k)
             U(nx+3,j,k) = -U(nx-1,j,k)
            end do
           end do
        else
           do k = 1, nz
            do j = 1, ny
             U(nx+1,j,k) = -U(nx,j,k)
             U(nx+2,j,k) = -U(nx-1,j,k)
             U(nx+3,j,k) = -U(nx-2,j,k)
            end do
           end do
        end if
      end if
    end if
    !$omp end parallel sections
    
#ifdef PAR    
    call par_exchange_U_y(U, component)
#endif

    !$omp parallel sections
    !$omp section
    if (Btype(So)==BC_DIRICHLET.and.reg/=intermediate) then
      if (component==2) then
        do k = 1, nz
         do i = -2, nx+3                       !Dirichlet inlet
          U(i,0,k) = sideU(component,So)
          U(i,-1,k) = sideU(component,So) + (sideU(component,So)-U(i,1,k))
          U(i,-2,k) = sideU(component,So) + (sideU(component,So)-U(i,2,k))
         end do
        end do
      else
        do k = 1, nz
         do i = -2, nx+3                       !Dirichlet inlet
          U(i,0,k) = sideU(component,So) + (sideU(component,So)-U(i,1,k))
          U(i,-1,k) = sideU(component,So) + (sideU(component,So)-U(i,2,k))
          U(i,-2,k) = sideU(component,So) + (sideU(component,So)-U(i,3,k))
         end do
        end do
      end if
    else if (Btype(So)==BC_NOSLIP .or. (component==2.and.Btype(So)==BC_FREESLIP) .or. &
               (Btype(So)==BC_DIRICHLET.and.reg==intermediate)) then
      if (component==2) then
        do k = 1, nz
         do i = -2, nx+3                       !Solid wall
          U(i,0,k) = 0
          U(i,-1,k) = -U(i,1,k)
          U(i,-2,k) = -U(i,2,k)
         end do
        end do
      else
        do k = 1, nz
         do i = -2, nx+3                       !Solid wall
          U(i,0,k) = -U(i,1,k)
          U(i,-1,k) = -U(i,2,k)
          U(i,-2,k) = -U(i,3,k)
         end do
        end do
      end if
    else if (Btype(So)==BC_NEUMANN.or.(component/=2.and.Btype(So)==BC_FREESLIP)) then
      do k = 1, nz
       do i = -2, nx+3                       !Neumann inlet
        U(i,0,k) = U(i,1,k)
        U(i,-1,k) = U(i,1,k)
        U(i,-2,k) = U(i,1,k)
       end do
      end do
    else if (Btype(So)==BC_PERIODIC) then  !Periodic BC
      do k = 1, nz
       do i = -2, nx+3
        U(i,0,k) = U(i,ny,k)
        U(i,-1,k) = U(i,ny-1,k)
        U(i,-2,k) = U(i,ny-2,k)
       end do
      end do
    else if (Btype(So)==BC_TURBULENT_INLET.or.Btype(So)==BC_INLET_FROM_FILE) then
      if (reg/=intermediate) then
        if (component==2) then        
          do k = -1, nz+2
           do i = -1, nx+2
            U(i,0,k) = Uin(i,k)
            U(i,-1,k) = Uin(i,k) + (Uin(i,k)-U(i,1,k))
            U(i,-2,k) = Uin(i,k) + (Uin(i,k)-U(i,2,k))
           end do
          end do
        else
          do k = -1, nz+2
           do i = -1, nx+2
            U(i,0,k) = Uin(i,k) + (Uin(i,k)-U(i,1,k))
            U(i,-1,k) = Uin(i,k) + (Uin(i,k)-U(i,2,k))
            U(i,-2,k) = Uin(i,k) + (Uin(i,k)-U(i,3,k))
           end do
          end do
        end if
      else
        if (component==2) then
           do k = -1, nz+2
            do i = -1, nx+2
             U(i,0,k) = 0
             U(i,-1,k) = -U(i,1,k)
             U(i,-2,k) = -U(i,2,k)
            end do
           end do
        else
           do k = -1, nz+2
            do i = -1, nx+2
             U(i,0,k) = -U(i,1,k)
             U(i,-1,k) = -U(i,2,k)
             U(i,-2,k) = -U(i,3,k)
            end do
           end do
        end if
      end if
    end if


    !$omp section
    if (Btype(No)==BC_DIRICHLET.and.reg/=intermediate) then
      if (component==2) then
        do k = 1, nz
         do i = -2, nx+3                       !Dirichlet inlet
          U(i,ny+1,k) = sideU(component,No)
          U(i,ny+2,k) = sideU(component,No) + (sideU(component,No)-U(i,ny,k))
          U(i,ny+3,k) = sideU(component,No) + (sideU(component,No)-U(i,ny-1,k))
         end do
        end do
      else
        do k = 1, nz
         do i = -2, nx+3                       !Dirichlet inlet
          U(i,ny+1,k) = sideU(component,No) + (sideU(component,No)-U(i,ny,k))
          U(i,ny+2,k) = sideU(component,No) + (sideU(component,No)-U(i,ny-1,k))
          U(i,ny+3,k) = sideU(component,No) + (sideU(component,No)-U(i,ny-2,k))
         end do
        end do
      end if
    else if (Btype(No)==BC_NOSLIP .or. (component==2.and.Btype(No)==BC_FREESLIP) .or. &
             (Btype(No)==BC_DIRICHLET.and.reg==intermediate)) then
      if (component==2) then
        do k = 1, nz
         do i = -2, nx+3                       !Solid wall
          U(i,ny+1,k) = 0
          U(i,ny+2,k) = -U(i,ny,k)
          U(i,ny+3,k) = -U(i,ny-1,k)
         end do
        end do
      else
        do k = 1, nz
         do i = -2, nx+3                       !Solid wall
          U(i,ny+1,k) = -U(i,ny,k)
          U(i,ny+2,k) = -U(i,ny-1,k)
          U(i,ny+3,k) = -U(i,ny-2,k)
         end do
        end do
      end if
    else if (Btype(No)==BC_NEUMANN .or. (component/=2.and.Btype(No)==BC_FREESLIP)) then
      do k = 1, nz
       do i = -2, nx+3                       !Neumann inlet
        U(i,ny+1,k) = U(i,ny,k)
        U(i,ny+2,k) = U(i,ny,k)
        U(i,ny+3,k) = U(i,ny,k)
       end do
      end do
    else if (Btype(No)==BC_PERIODIC) then  !Periodic BC
      do k = 1, nz
       do i = -2, nx+3
        U(i,ny+1,k) = U(i,1,k)
        U(i,ny+2,k) = U(i,2,k)
        U(i,ny+3,k) = U(i,3,k)
       end do
      end do
    else if (Btype(No)==BC_TURBULENT_INLET.or.Btype(No)==BC_INLET_FROM_FILE) then
      if (reg/=intermediate) then
        if (component==2) then
          do k = -1, nz+2
           do i = -1, nx+2
            U(i,ny+1,k) = Uin(i,k)
            U(i,ny+2,k) = Uin(i,k) + (Uin(i,k)-U(i,ny,k))
            U(i,ny+3,k) = Uin(i,k) + (Uin(i,k)-U(i,ny-1,k))
           end do
          end do
        else
          do k = -1, nz+2
           do i = -1, nx+2
            U(i,ny+1,k) = Uin(i,k) + (Uin(i,k)-U(i,ny,k))
            U(i,ny+2,k) = Uin(i,k) + (Uin(i,k)-U(i,ny-1,k))
            U(i,ny+3,k) = Uin(i,k) + (Uin(i,k)-U(i,ny-2,k))
           end do
          end do
        end if
      else
        if (component==2) then
           do k = -1, nz+2
            do i = -1, nx+2
             U(i,ny+1,k) = 0
             U(i,ny+2,k) = -U(i,ny,k)
             U(i,ny+3,k) = -U(i,ny-1,k)
            end do
           end do
        else
           do k = -1, nz+2
            do i = -1, nx+2
             U(i,ny+1,k) = -U(i,ny,k)
             U(i,ny+2,k) = -U(i,ny-1,k)
             U(i,ny+3,k) = -U(i,ny-2,k)
            end do
           end do
        end if
      end if
    end if
    !$omp end parallel sections

#ifdef PAR
    call par_exchange_U_z(U, component)
#endif

    !$omp parallel sections
    !$omp section
    if (Btype(Bo)==BC_DIRICHLET .and. reg/=intermediate) then
      if (component==3) then
        do j = -2, ny+3
         do i = -2, nx+3                       !Dirichlet inlet
          U(i,j,0) = sideU(component,Bo)
          U(i,j,-1) = sideU(component,Bo) + (sideU(component,Bo)-U(i,j,1))
          U(i,j,-2) = sideU(component,Bo) + (sideU(component,Bo)-U(i,j,2))
         end do
        end do
      else
        do j = -2, ny+3
         do i = -2, nx+3                       !Dirichlet inlet
          U(i,j,0) = sideU(component,Bo) + (sideU(component,Bo)-U(i,j,1))
          U(i,j,-1) = sideU(component,Bo) + (sideU(component,Bo)-U(i,j,2))
          U(i,j,-2) = sideU(component,Bo) + (sideU(component,Bo)-U(i,j,3))
         end do
        end do
      end if
    else if (Btype(Bo)==BC_NOSLIP .or. (component==3.and.Btype(Bo)==BC_FREESLIP) .or. &
               (Btype(Bo)==BC_DIRICHLET.and.reg==intermediate)) then
      if (component==3) then
        do j = -2, ny+3
         do i = -2, nx+3                       !Solid wall
          U(i,j,0) = 0
          U(i,j,-1) = -U(i,j,1)
          U(i,j,-2) = -U(i,j,2)
         end do
        end do
      else
        do j = -2, ny+3
         do i = -2, nx+3                       !Solid wall
          U(i,j,0) = -U(i,j,1)
          U(i,j,-1) = -U(i,j,2)
          U(i,j,-2) = -U(i,j,3)
         end do
        end do
      end if
    else if (Btype(Bo)==BC_NEUMANN.or.(component/=3.and.Btype(Bo)==BC_FREESLIP)) then
      do j = -2, ny+3
       do i = -2, nx+3                       !Neumann inlet
        U(i,j,0) = U(i,j,1)
        U(i,j,-1) = U(i,j,1)
        U(i,j,-2) = U(i,j,1)
       end do
      end do
    else if (Btype(Bo)==BC_PERIODIC) then  !Periodic BC
     do j = -2, ny+3
      do i = -2, nx+3
       U(i,j,0) = U(i,j,nz)
       U(i,j,-1) = U(i,j,nz-1)
       U(i,j,-2) = U(i,j,nz-2)
      end do
     end do
    end if

    !$omp section
    if (Btype(To)==BC_DIRICHLET.and.reg/=intermediate) then
      if (component==3) then
        do j = -2, ny+3
         do i = -2, nx+3                       !Dirichlet inlet
          U(i,j,nz+1) = sideU(component,To)
          U(i,j,nz+2) = sideU(component,To) + (sideU(component,To)-U(i,j,nz))
          U(i,j,nz+3) = sideU(component,To) + (sideU(component,To)-U(i,j,nz-1))
         end do
        end do
      else
        do j = -2, ny+3
         do i = -2, nx+3                       !Dirichlet inlet
          U(i,j,nz+1) = sideU(component,To) + (sideU(component,To)-U(i,j,nz))
          U(i,j,nz+2) = sideU(component,To) + (sideU(component,To)-U(i,j,nz-1))
          U(i,j,nz+3) = sideU(component,To) + (sideU(component,To)-U(i,j,nz-2))
         end do
        end do
      end if
    else if (Btype(To)==BC_NOSLIP .or. &
             (component==3.and.(Btype(To)==BC_FREESLIP .or. &
                                Btype(To)==BC_AUTOMATIC_FLUX)) .or. &
             (Btype(To)==BC_DIRICHLET.and.reg==intermediate) ) then
      if (component==3) then
        do j = -2, ny+3
         do i = -2, nx+3                       !Solid wall
          U(i,j,nz+1) = 0
          U(i,j,nz+2) = -U(i,j,nz)
          U(i,j,nz+3) = -U(i,j,nz-1)
         end do
        end do
      else
        do j = -2, ny+3
         do i = -2, nx+3                       !Solid wall
          U(i,j,nz+1) = -U(i,j,nz)
          U(i,j,nz+2) = -U(i,j,nz-1)
          U(i,j,nz+3) = -U(i,j,nz-2)
         end do
        end do
      end if
    else if (Btype(To)==BC_NEUMANN .or. &
             (component/=3.and.(Btype(To)==BC_FREESLIP .or. &
                                Btype(To)==BC_AUTOMATIC_FLUX))) then
      do j = -2, ny+3
       do i = -2, nx+3                       !Neumann inlet
        U(i,j,nz+1) = U(i,j,nz)
        U(i,j,nz+2) = U(i,j,nz)
        U(i,j,nz+3) = U(i,j,nz)
       end do
      end do
    else if (Btype(To)==BC_PERIODIC) then  !Periodic BC
     do j = -2, ny+3
      do i = -2, nx+3
       U(i,j,nz+1) = U(i,j,1)
       U(i,j,nz+2) = U(i,j,2)
       U(i,j,nz+3) = U(i,j,3)
      end do
     end do
    end if
    !$omp end parallel sections
  end subroutine BoundU


  subroutine BoundUVW(U, V, W, regime)
#ifdef PAR
    use domains_bc_par
#endif
    real(knd), contiguous, intent(inout)     :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    integer, optional, intent(in) :: regime
#ifdef PAR
    integer :: reg

    reg = 1
    if (present(regime)) reg=regime

    if (reg/=2) call par_update_domain_bounds_UVW(U, V, W, time_stepping%effective_time)

#endif

    call BoundU(1, U, Uin, regime)
    call BoundU(2, V, Vin, regime)
    call BoundU(3, W, Win, regime)
  end subroutine



  subroutine Bound_Phi(Phi)
    real(knd), intent(inout) :: Phi(-1:,-1:,-1:)
    integer :: i, j, k, nx, ny, nz

    nx = Prnx
    ny = Prny
    nz = Prnz

#ifdef PAR
    call par_exchange_boundaries(Phi, Btype, 4)
#endif
    
    if (Btype(We)==BC_PERIODIC) then
     do k = 1, nz
      do j = 1, ny                      !Periodic BC
       Phi(-1:0,j,k) = Phi(nx-1:nx,j,k)
      end do
     end do
    else if ((Btype(We)<BC_MPI_BOUNDS_MIN) .or. (Btype(We)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do j = 1, ny                      !Other BCs
       Phi(-1:0,j,k) = Phi(2:1:-1,j,k)
      end do
     end do
    end if

    if (Btype(Ea)==BC_PERIODIC) then
     do k = 1, nz
      do j = 1, ny                      !Periodic BC
       Phi(nx+1:nx+2,j,k) = Phi(1:2,j,k)
      end do
     end do
    else if ((Btype(Ea)<BC_MPI_BOUNDS_MIN) .or. (Btype(Ea)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do j = 1, ny                      !Other BCs
       Phi(nx+1:nx+2,j,k) = Phi(nx:nx-1:-1,j,k)
      end do
     end do
    end if

    if (Btype(So)==BC_PERIODIC) then
     do k = 1, nz
      do i = 1, nx                      !Periodic BC
       Phi(i,-1:0,k) = Phi(i,ny-1:ny,k)
      end do
     end do
    else if ((Btype(So)<BC_MPI_BOUNDS_MIN) .or. (Btype(So)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do i = 1, nx                      !Other BCs
       Phi(i,-1:0,k) = Phi(i,2:1:-1,k)
      end do
     end do
    end if

    if (Btype(No)==BC_PERIODIC) then
    do k = 1, nz
      do i = 1, nx                      !Periodic BC
       Phi(i,ny+1:ny+2,k) = Phi(i,1:2,k)
      end do
     end do
    else if ((Btype(No)<BC_MPI_BOUNDS_MIN) .or. (Btype(No)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do i = 1, nx                      !Other BCs
       Phi(i,ny+1:ny+2,k) = Phi(i,ny:ny-1:-1,k)
      end do
     end do
    end if

    if (Btype(Bo)==BC_PERIODIC) then
     do j = 1, ny
      do i = 1, nx                      !Periodic BC
       Phi(i,j,-1:0) = Phi(i,j,nz-1:nz)
      end do
     end do
    else if ((Btype(Bo)<BC_MPI_BOUNDS_MIN) .or. (Btype(Bo)>BC_MPI_BOUNDS_MAX)) then
     do j = 1, ny
      do i = 1, nx                      !Other BCs
       Phi(i,j,-1:0) = Phi(i,j,2:1:-1)
      end do
     end do
    end if

    if (Btype(To)==BC_PERIODIC) then
     do j = 1, ny
      do i = 1, nx                      !Periodic BC
       Phi(i,j,nz+1:nz+2) = Phi(i,j,1:2)
      end do
     end do
    else if ((Btype(To)<BC_MPI_BOUNDS_MIN) .or. (Btype(To)>BC_MPI_BOUNDS_MAX)) then
     do j = 1, ny
      do i = 1, nx                      !Other BCs
       Phi(i,j,nz+1:nz+2) = Phi(i,j,nz:nz-1:-1)
      end do
     end do
    end if
  end subroutine Bound_Phi

  
  subroutine Bound_Pr(Pr)
    real(knd), intent(inout) :: Pr(-1:,-1:,-1:)
    integer :: i, j, k, nx, ny, nz

    nx = Prnx
    ny = Prny
    nz = Prnz
    
#ifdef PAR
    call par_exchange_Pr(Pr)
#endif

    !extrapolate so that the 4th order gradient reduces to 2nd order
    if (Btype(We)==BC_PERIODIC) then
     do k = 1, nz
      do j = 1, ny
       Pr(0,j,k) = Pr(nx,j,k)
      end do
     end do
    else if ((Btype(We)<BC_MPI_BOUNDS_MIN) .or. (Btype(We)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do j = 1, ny
       Pr(0,j,k) = Pr(3,j,k) - 3 * (Pr(2,j,k) - Pr(1,j,k))
      end do
     end do
    end if
    
    if (Btype(Ea)==BC_PERIODIC) then
     do k = 1, nz
      do j = 1, ny
       Pr(nx+1:nx+2,j,k) = Pr(1:2,j,k)
      end do
     end do
    else if ((Btype(Ea)<BC_MPI_BOUNDS_MIN) .or. (Btype(Ea)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do j = 1, ny
       Pr(nx+1,j,k) = Pr(nx-2,j,k) - 3 * (Pr(nx-1,j,k) - Pr(nx,j,k))
      end do
     end do
    end if
    
    if (Btype(So)==BC_PERIODIC) then
     do k = 1, nz
      do i = 1, nx
       Pr(i,0,k) = Pr(i,ny,k)
      end do
     end do
    else if ((Btype(So)<BC_MPI_BOUNDS_MIN) .or. (Btype(So)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do i = 1, nx
       Pr(i,0,k) = Pr(i,j,k) - 3 * (Pr(i,2,k) - Pr(i,1,k))
      end do
     end do
    end if
    
    if (Btype(No)==BC_PERIODIC) then
     do k = 1, nz
      do i = 1, nx
       Pr(i,ny+1:ny+2,k) = Pr(i,1:2,k)
      end do
     end do
    else if ((Btype(No)<BC_MPI_BOUNDS_MIN) .or. (Btype(No)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do i = 1, nx
       Pr(i,ny+1,k) = Pr(i,ny-2,k) - 3 * (Pr(i,ny-1,k) - Pr(i,ny,k))
      end do
     end do
    end if
    
    if (Btype(Bo)==BC_PERIODIC) then
     do j = 1, ny
      do i = 1, nx
       Pr(i,j,0) = Pr(i,j,nz)
      end do
     end do
    else if ((Btype(Bo)<BC_MPI_BOUNDS_MIN) .or. (Btype(Bo)>BC_MPI_BOUNDS_MAX)) then
     do j = 1, ny
      do i = 1, nx
       Pr(i,j,0) = Pr(i,j,3) - 3 * (Pr(i,j,2) - Pr(i,j,1))
      end do
     end do
    end if
    
    if (Btype(To)==BC_PERIODIC) then
     do j = 1, ny
      do i = 1, nx
       Pr(i,j,nz+1:nz+2) = Pr(i,j,1:2)
      end do
     end do
    else if ((Btype(To)<BC_MPI_BOUNDS_MIN) .or. (Btype(To)>BC_MPI_BOUNDS_MAX)) then
     do j = 1, ny
      do i = 1, nx
       Pr(i,j,nz+1) = Pr(i,j,nz-2) - 3 * (Pr(i,j,nz-1) - Pr(i,j,nz))
      end do
     end do
    end if

  end subroutine Bound_Pr



  subroutine Bound_Q(Phi)
    real(knd), intent(inout) :: Phi(0:,0:,0:)
    integer :: i, j, k, nx, ny, nz

    nx = Prnx
    ny = Prny
    nz = Prnz

#ifdef PAR
    call par_exchange_Q(Phi)
#endif

    !The above filled the buffers, but we have to add
    ! the content of the buffers to the inside points.

    if (Btype(We)==BC_PERIODIC) then
     do k = 1, nz
      do j = 1, ny                      !Periodic BC
       Phi(nx,j,k) = Phi(0,j,k) + Phi(nx,j,k)
      end do
     end do
    else if ((Btype(We)<BC_MPI_BOUNDS_MIN) .or. (Btype(We)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do j = 1, ny                      !Other BCs
       Phi(1,j,k) = Phi(1,j,k) + Phi(0,j,k)
      end do
     end do
    end if

    if (Btype(Ea)==BC_PERIODIC) then
     do k = 1, nz
      do j = 1, ny                      !Periodic BC
       Phi(1,j,k) = Phi(1,j,k) + Phi(nx+1,j,k)
      end do
     end do
    else if ((Btype(Ea)<BC_MPI_BOUNDS_MIN) .or. (Btype(Ea)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do j = 1, ny                      !Other BCs
       Phi(nx,j,k) = Phi(nx,j,k) + Phi(nx+1,j,k)
      end do
     end do
    end if

    if (Btype(So)==BC_PERIODIC) then
     do k = 1, nz
      do i = 1, nx                      !Periodic BC
       Phi(i,ny,k) = Phi(i,ny,k) + Phi(i,0,k)
      end do
     end do
    else if ((Btype(So)<BC_MPI_BOUNDS_MIN) .or. (Btype(So)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do i = 1, nx                      !Other BCs
       Phi(i,1,k) = Phi(i,1,k) + Phi(i,0,k)
      end do
     end do
    end if

    if (Btype(No)==BC_PERIODIC) then
     do k = 1, nz
      do i = 1, nx                      !Periodic BC
       Phi(i,1,k) = Phi(i,1,k) + Phi(i,ny+1,k)
      end do
     end do
    else if ((Btype(No)<BC_MPI_BOUNDS_MIN) .or. (Btype(No)>BC_MPI_BOUNDS_MAX)) then
     do k = 1, nz
      do i = 1, nx                      !Other BCs
       Phi(i,ny,k) = Phi(i,ny,k) + Phi(i,ny+1,k)
      end do
     end do
    end if

    if (Btype(Bo)==BC_PERIODIC) then
     do j = 1, ny
      do i = 1, nx                      !Periodic BC
       Phi(i,j,nz) = Phi(i,j,nz) + Phi(i,j,0)
      end do
     end do
    else if ((Btype(Bo)<BC_MPI_BOUNDS_MIN) .or. (Btype(Bo)>BC_MPI_BOUNDS_MAX)) then
     do j = 1, ny
      do i = 1, nx                      !Other BCs
       Phi(i,j,1) = Phi(i,j,1) + Phi(i,j,0)
      end do
     end do
    end if

    if (Btype(To)==BC_PERIODIC) then
     do j = 1, ny
      do i = 1, nx                      !Periodic BC
       Phi(i,j,1) = Phi(i,j,1) + Phi(i,j,nz+1)
      end do
     end do
    else if ((Btype(To)<BC_MPI_BOUNDS_MIN) .or. (Btype(To)>BC_MPI_BOUNDS_MAX)) then
     do j = 1, ny
      do i = 1, nx                      !Other BCs
       Phi(i,j,nz) = Phi(i,j,nz) + Phi(i,j,nz+1)
      end do
     end do
    end if
  end subroutine Bound_Q















  subroutine ConstantInlet

    Uin = Uinlet * cos(windangle/180._knd*pi)
    Vin = Uinlet * sin(windangle/180._knd*pi)
    Win = 0

  end subroutine

  subroutine ShearInlet(G)
    real(knd) :: G
    integer :: j, k

    do k = lbound(Uin,2), ubound(Uin,2)
     do j = lbound(Uin,1), ubound(Uin,1)
       Uin(j,k) = G * (zPr(k) - ( (zW(Wnz+1)+zW(0)) / 2._knd ) )
     end do
    end do

    Vin = 0
    Win = 0

  end subroutine

  subroutine ParabolicInlet
    integer :: j, k
    real(knd) :: lz

    lz = zW(Prnz) - zW(0)

    do k = lbound(Uin,2), ubound(Uin,2)
     do j = lbound(Uin,1), ubound(Uin,1)
       Uin(j,k) = 1.5_knd * Uinlet * &
                    (1 - &
                      ( &
                        (lz/2._knd - (zPr(k)-zW(0))) / (lz/2._knd) &
                      )**2  &
                    )
     end do
    end do

    Vin = 0
    Win = 0

  end subroutine


end module Boundaries
