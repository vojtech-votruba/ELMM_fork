module ScalarBoundaries
  use Parameters
#ifdef PAR
  use exchange_par
#endif

  implicit none

  private

  public BoundTemperature, BoundMoisture, BoundScalar, &
         BoundViscosity
  
contains



  subroutine CommonBase(arr,Sc_btype,side,in,BsideArr,BsideFlArr)
    real(knd),intent(inout) :: arr(-2:,-2:,-2:)
    integer  ,intent(in)    :: Sc_btype(6)
    real(knd),intent(in)    :: side(6)
    real(knd),intent(in),optional :: in(-2:,-2:)
    real(knd),intent(in),allocatable,optional :: BsideArr(:,:), BsideFlArr(:,:)
    integer :: i,j,k,nx,ny,nz

    nx = Prnx
    ny = Prny
    nz = Prnz

#ifdef PAR
    call par_exchange_Sc_x(arr,Sc_btype)
#endif
   
    if (Sc_btype(We)==BC_DIRICHLET) then
      !TODO: use these types as scalar boundary condiitons!
      if (present(in) .and. &
          (Btype(We)==BC_TURBULENT_INLET.or.Btype(We)==BC_INLET_FROM_FILE)) then
        do k = 1,nz
          do j = 1,ny
            arr(-2:0,j,k)  = in(j,k) - (arr(3:1:-1,j,k)-in(j,k))
          end do
        end do
      else
        do k = 1,nz
          do j = 1,ny
            arr(-2:0,j,k)  = side(We) - (arr(3:1:-1,j,k)-side(We))
          end do
        end do
      end if 
    else if (Sc_btype(We)==BC_PERIODIC) then
      do k = 1,nz
        do j = 1,ny
         arr(-2:0,j,k)  = arr(nx-2:nx,j,k)
        end do
      end do
    else if (Sc_btype(We)==BC_NEUMANN) then
      do k = 1,nz
        do j = 1,ny
          arr(0,j,k)  = arr(1,j,k) - side(We)*dxmin
          arr(-1,j,k) = arr(1,j,k) - side(We)*2*dxmin
          arr(-2,j,k) = arr(1,j,k) - side(We)*3*dxmin
        end do
      end do
    else if (Sc_btype(We)==BC_CONSTFLUX) then
      do k = 1,nz
        do j = 1,ny
          arr(0,j,k)  = arr(1,j,k) + &
                side(We)*dxU(0)/((TDiff(1,j,k)+TDiff(0,j,k))/(2._knd))
          arr(-1,j,k) = arr(0,j,k) - (arr(1,j,k)-arr(0,j,k))
          arr(-2,j,k) = arr(0,j,k) - 2*(arr(1,j,k)-arr(0,j,k))
        end do
      end do
    else if (((Sc_btype(We)<BC_MPI_BOUNDS_MIN) .or. (Sc_btype(We)>BC_MPI_BOUNDS_MAX)) .and. &
             ((Sc_btype(We)<BC_DOMAIN_BOUNDS_MIN) .or. (Sc_btype(We)>BC_DOMAIN_BOUNDS_MAX))) then
      do k = 1,nz
        do j = 1,ny
          arr(-2:0,j,k)  = arr(3:1:-1,j,k)
        end do
      end do
    end if

    if (Sc_btype(Ea)==BC_DIRICHLET) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1:nx+3,j,k) = side(Ea) - (arr(nx:nx-2:-1,j,k)-side(Ea))
        end do
      end do
    else if (Sc_btype(Ea)==BC_PERIODIC) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1:nx+3,j,k) = arr(1:3,j,k)
        end do
      end do
    else if (Sc_btype(Ea)==BC_NEUMANN) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1,j,k) = arr(nx,j,k) + side(Ea)*dxmin
          arr(nx+2,j,k) = arr(nx,j,k) + side(Ea)*2*dxmin
          arr(nx+3,j,k) = arr(nx,j,k) + side(Ea)*3*dxmin
        end do
      end do
    else if (Sc_btype(Ea)==BC_CONSTFLUX) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1,j,k) = arr(nx,j,k) - &
              side(Ea)*dxU(nx)/((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2._knd))
          arr(nx+2,j,k) = arr(nx+1,j,k) - (arr(nx,j,k)-arr(nx+1,j,k))
          arr(nx+2,j,k) = arr(nx+1,j,k) - 2*(arr(nx,j,k)-arr(nx+1,j,k))
        end do
      end do
    else if (((Sc_btype(Ea)<BC_MPI_BOUNDS_MIN) .or. (Sc_btype(Ea)>BC_MPI_BOUNDS_MAX)) .and. &
             ((Sc_btype(Ea)<BC_DOMAIN_BOUNDS_MIN) .or. (Sc_btype(Ea)>BC_DOMAIN_BOUNDS_MAX))) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1:nx+3,j,k) = arr(nx:nx-2:-1,j,k)
        end do
      end do
    end if

#ifdef PAR
    call par_exchange_Sc_y(arr,Sc_btype)
#endif
    
    if (Sc_btype(So)==BC_DIRICHLET) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,-2:0,k)  = side(So) - (arr(i,3:1:-1,k)-side(So))
        end do
      end do
    else if (Sc_btype(So)==BC_PERIODIC) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,-2:0,k)  = arr(i,ny-2:ny,k)
        end do
      end do
    else if (Sc_btype(So)==BC_NEUMANN) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,0,k)  = arr(i,1,k) - side(So)*dymin
          arr(i,-1,k) = arr(i,1,k) - side(So)*2*dymin
          arr(i,-2,k) = arr(i,1,k) - side(So)*3*dymin
        end do
      end do
    else if (Sc_btype(So)==BC_CONSTFLUX) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,0,k)  = arr(i,1,k) + &
                side(So)*dyV(0)/((TDiff(i,1,k)+TDiff(i,0,k))/(2._knd))
          arr(i,-1,k) = arr(i,0,k) - (arr(i,1,k)-arr(i,0,k))
          arr(i,-2,k) = arr(i,0,k) - 2*(arr(i,1,k)-arr(i,0,k))
        end do
      end do
    else if ((Sc_btype(So)<BC_MPI_BOUNDS_MIN) .or. (Sc_btype(So)>BC_MPI_BOUNDS_MAX) .and. &
             ((Sc_btype(So)<BC_DOMAIN_BOUNDS_MIN) .or. (Sc_btype(So)>BC_DOMAIN_BOUNDS_MAX))) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,-2:0,k)  = arr(i,3:1:-1,k)
        end do
      end do
    end if

    if (Sc_btype(No)==BC_DIRICHLET) then
      do k = 1,nz
        do i=-1,nx+2
           arr(i,ny+1:ny+3,k) = side(No) - (arr(i,ny:ny-2:-1,k)-side(No))
         end do
       end do
    else if (Sc_btype(No)==BC_PERIODIC) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,ny+1:ny+3,k) = arr(i,1:3,k)
        end do
      end do
    else if (Sc_btype(No)==BC_NEUMANN) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,ny+1,k) = arr(i,ny,k) + side(No)*dymin
          arr(i,ny+2,k) = arr(i,ny,k) + side(No)*2*dymin
          arr(i,ny+3,k) = arr(i,ny,k) + side(No)*3*dymin
        end do
      end do
    else if (Sc_btype(No)==BC_CONSTFLUX) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,ny+1,k) = arr(i,ny+1,k) - &
              side(No)*dyV(ny)/((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2._knd))
          arr(i,ny+2,k) = arr(i,ny+1,k) - (arr(i,ny,k)-arr(i,ny+1,k))
          arr(i,ny+3,k) = arr(i,ny+1,k) - 2*(arr(i,ny,k)-arr(i,ny+1,k))
        end do
      end do
    else if (((Sc_btype(No)<BC_MPI_BOUNDS_MIN) .or. (Sc_btype(No)>BC_MPI_BOUNDS_MAX)) .and. &
             ((Sc_btype(No)<BC_DOMAIN_BOUNDS_MIN) .or. (Sc_btype(No)>BC_DOMAIN_BOUNDS_MAX))) then    
      do k = 1,nz
        do i=-1,ny+2
          arr(i,ny+1:ny+3,k) = arr(i,ny:ny-2:-1,k)
        end do
      end do
    end if

#ifdef PAR
    call par_exchange_Sc_z(arr,Sc_btype)
#endif
    
    if (Sc_btype(Bo)==BC_DIRICHLET) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,-2:0) = side(Bo) - (arr(i,j,3:1:-1)-side(Bo))
        end do
      end do
    else if (Sc_btype(Bo)==BC_PERIODIC) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,-2:0)  = arr(i,j,nz:nz-2:-1)
        end do
      end do
    else if (Sc_btype(Bo)==BC_NEUMANN) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,0)  = arr(i,j,1) - side(Bo)*dzW(0)
          arr(i,j,-1) = arr(i,j,1) - side(Bo)*(zPr(1)-zPr(-1))
          arr(i,j,-2) = arr(i,j,1) - side(Bo)*(zPr(1)-zPr(-2))
        end do
      end do
    else if (Sc_btype(Bo)==BC_CONSTFLUX) then
      do j=-1,ny+2
        do i=-1,nx+2
            arr(i,j,0) = arr(i,j,1) + &
                side(Bo)*dzW(0)/((TDiff(i,j,1)+TDiff(i,j,0))/(2._knd))
          arr(i,j,-1) = arr(i,j,0) - (arr(i,j,1)-arr(i,j,0))
          arr(i,j,-2) = arr(i,j,0) - 2*(arr(i,j,1)-arr(i,j,0))
        end do
      end do
    else if (((Sc_btype(Bo)<BC_MPI_BOUNDS_MIN) .or. (Sc_btype(Bo)>BC_MPI_BOUNDS_MAX)) .and. &
             ((Sc_btype(Bo)<BC_DOMAIN_BOUNDS_MIN) .or. (Sc_btype(Bo)>BC_DOMAIN_BOUNDS_MAX))) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,-2:0)  = arr(i,j,3:1:-1)
        end do
      end do
    end if

    if (Sc_btype(To)==BC_DIRICHLET) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1:nz+3) = side(To) - (arr(i,j,nz:nz-2:-1)-side(To))
        end do
      end do
    else if (Sc_btype(To)==BC_PERIODIC) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1:3) = arr(i,j,1:3)
        end do
      end do
    else if (Sc_btype(To)==BC_NEUMANN) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1) = arr(i,j,nz) + side(To)*dzW(nz+1)
          arr(i,j,nz+2) = arr(i,j,nz) + side(To)*(zPr(nz+2)-zPr(nz))
          arr(i,j,nz+3) = arr(i,j,nz) + side(To)*(zPr(nz+3)-zPr(nz))
        end do
      end do
    else if (Sc_btype(To)==BC_CONSTFLUX.or.Sc_btype(To)==BC_AUTOMATIC_FLUX) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1) = arr(i,j,nz) - &
              side(To)*dzW(nz)/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._knd))
          arr(i,j,nz+2) = arr(i,j,nz+1) - (arr(i,j,nz)-arr(i,j,nz+1))
          arr(i,j,nz+3) = arr(i,j,nz+1) - 2*(arr(i,j,nz)-arr(i,j,nz+1))
        end do
      end do
    else if (((Sc_btype(To)<BC_MPI_BOUNDS_MIN) .or. (Sc_btype(To)>BC_MPI_BOUNDS_MAX)) .and. &
             ((Sc_btype(To)<BC_DOMAIN_BOUNDS_MIN) .or. (Sc_btype(To)>BC_DOMAIN_BOUNDS_MAX))) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1:nz+3) = arr(i,j,nz:nz-2:-1)
        end do
      end do
     end if
  end subroutine CommonBase




  subroutine BoundScalar(SCAL)
    real(knd),intent(inout) :: Scal(-2:,-2:,-2:)

    call CommonBase(SCAL,ScalBtype,sideScal)

  end subroutine BoundScalar






  subroutine BoundTemperature(Temperature)
#ifdef PAR
    use domains_bc_par
#endif
    real(knd),intent(inout) :: Temperature(-2:,-2:,-2:)
    
#ifdef PAR
    call par_update_domain_bounds_temperature(Temperature, time_stepping%effective_time)
#endif

    call CommonBase(Temperature,TempBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr)

  end subroutine BoundTemperature


  subroutine BoundMoisture(Moisture)
#ifdef PAR
    use domains_bc_par
#endif
    real(knd),intent(inout) :: Moisture(-2:,-2:,-2:)

#ifdef PAR
    call par_update_domain_bounds_moisture(Moisture, time_stepping%effective_time)
#endif

    call CommonBase(Moisture,MoistBtype,sideMoist, MoistIn,BsideMArr,BsideMFLArr)

  end subroutine BoundMoisture


  subroutine BoundViscosity(Nu)
    real(knd),contiguous,intent(inout) :: Nu(-1:,-1:,-1:)
    integer :: i,j,k,nx,ny,nz

    nx = Prnx
    ny = Prny
    nz = Prnz
    
#ifdef PAR
    call par_exchange_visc_x(Nu, Btype)
#endif

    if (Btype(Ea)==BC_PERIODIC) then
      do k = 1,nz
       do j = 1,ny
         Nu(0,j,k) = Nu(nx,j,k)
         Nu(nx+1,j,k) = Nu(1,j,k)
       end do
      end do
    else if ((Btype(Ea)<BC_MPI_BOUNDS_MIN) .or. (Btype(Ea)>BC_MPI_BOUNDS_MAX)) then
      do k = 1,nz
       do j = 1,ny
         Nu(0,j,k) = Nu(1,j,k)
         Nu(nx+1,j,k) = Nu(nx,j,k)
       end do
      end do
    end if

#ifdef PAR
    call par_exchange_visc_y(Nu, Btype)
#endif

    if (Btype(No)==BC_PERIODIC) then
      do k = 1,nz
       do i = 0,nx+1
         Nu(i,0,k) = Nu(i,ny,k)
         Nu(i,ny+1,k) = Nu(i,1,k)
       end do
      end do
    else if ((Btype(No)<BC_MPI_BOUNDS_MIN) .or. (Btype(No)>BC_MPI_BOUNDS_MAX)) then
      do k = 1,nz
       do i = 0,nx+1
         Nu(i,0,k) = Nu(i,1,k)
         Nu(i,ny+1,k) = Nu(i,ny,k)
       end do
      end do
    end if
    
#ifdef PAR
    call par_exchange_visc_z(Nu, Btype)
#endif

    if (Btype(To)==BC_PERIODIC) then
      do j = 0,ny+1
       do i = 0,nx+1
         Nu(i,j,0) = Nu(i,j,nz)
         Nu(i,j,nz+1) = Nu(i,j,1)
       end do
      end do
    else if ((Btype(To)<BC_MPI_BOUNDS_MIN) .or. (Btype(To)>BC_MPI_BOUNDS_MAX)) then
      do j = 0,ny+1
       do i = 0,nx+1
         Nu(i,j,0) = Nu(i,j,1)
         Nu(i,j,nz+1) = Nu(i,j,nz)
       end do
      end do
    end if
    
  end subroutine BoundViscosity




end module ScalarBoundaries
