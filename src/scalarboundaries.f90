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



  subroutine CommonBase(arr,btype,side,in,BsideArr,BsideFlArr)
    real(knd),intent(inout) :: arr(-1:,-1:,-1:)
    integer  ,intent(in)    :: btype(6)
    real(knd),intent(in)    :: side(6)
    real(knd),intent(in),optional :: in(-1:,-1:)
    real(knd),intent(in),allocatable,optional :: BsideArr(:,:), BsideFlArr(:,:)
    integer :: i,j,k,nx,ny,nz

    nx = Prnx
    ny = Prny
    nz = Prnz

#ifdef PAR
    call par_exchange_Sc_x(arr,btype)
#endif
    
    if (btype(We)==BC_DIRICHLET) then
      do k = 1,nz
        do j = 1,ny
          if (present(in)) then
            arr(0,j,k)  = in(j,k) - (arr(1,j,k)-in(j,k))
            arr(-1,j,k) = in(j,k) - (arr(2,j,k)-in(j,k))
          else
            arr(0,j,k)  = side(We) - (arr(1,j,k)-side(We))
            arr(-1,j,k) = side(We) - (arr(2,j,k)-side(We))
          end if
        end do
      end do
    else if (btype(We)==BC_PERIODIC) then
      do k = 1,nz
        do j = 1,ny
         arr(0,j,k)  = arr(nx,j,k)
         arr(-1,j,k) = arr(nx-1,j,k)
        end do
      end do
    else if (btype(We)==BC_NEUMANN) then
      do k = 1,nz
        do j = 1,ny
          arr(0,j,k)  = arr(1,j,k) - side(We)*dxU(0)
          arr(-1,j,k) = arr(1,j,k) - side(We)*(dxU(0)+dxU(-1))
        end do
      end do
    else if (btype(We)==BC_CONSTFLUX) then
      do k = 1,nz
        do j = 1,ny
          arr(0,j,k)  = arr(1,j,k) + &
                side(We)*dxU(0)/((TDiff(1,j,k)+TDiff(0,j,k))/(2._knd))
          arr(-1,j,k) = arr(0,j,k) - (arr(1,j,k)-arr(0,j,k))
        end do
      end do
    else if (((btype(We)<BC_MPI_BOUNDS_MIN) .or. (btype(We)>BC_MPI_BOUNDS_MAX)) .and. &
             ((btype(We)<BC_DOMAIN_BOUNDS_MIN) .or. (btype(We)>BC_DOMAIN_BOUNDS_MAX))) then
      do k = 1,nz
        do j = 1,ny
          arr(0,j,k)  = arr(1,j,k)
          arr(-1,j,k) = arr(2,j,k)
        end do
      end do
    end if

    if (btype(Ea)==BC_DIRICHLET) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1,j,k) = side(Ea) - (arr(nx,j,k)-side(Ea))
          arr(nx+2,j,k) = side(Ea) - (arr(nx-1,j,k)-side(Ea))
        end do
      end do
    else if (btype(Ea)==BC_PERIODIC) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1,j,k) = arr(1,j,k)
          arr(nx+2,j,k) = arr(2,j,k)
        end do
      end do
    else if (btype(Ea)==BC_NEUMANN) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1,j,k) = arr(nx,j,k) + side(Ea)*dxU(nx+1)
          arr(nx+2,j,k) = arr(nx,j,k) + side(Ea)*(dxU(nx+1)+dxU(nx+2))
        end do
      end do
    else if (btype(Ea)==BC_CONSTFLUX) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1,j,k) = arr(nx,j,k) - &
              side(Ea)*dxU(nx)/((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2._knd))
          arr(nx+2,j,k) = arr(nx+1,j,k) - (arr(nx,j,k)-arr(nx+1,j,k))
        end do
      end do
    else if (((btype(Ea)<BC_MPI_BOUNDS_MIN) .or. (btype(Ea)>BC_MPI_BOUNDS_MAX)) .and. &
             ((btype(Ea)<BC_DOMAIN_BOUNDS_MIN) .or. (btype(Ea)>BC_DOMAIN_BOUNDS_MAX))) then
      do k = 1,nz
        do j = 1,ny
          arr(nx+1,j,k) = arr(nx,j,k)
          arr(nx+2,j,k) = arr(nx-1,j,k)
        end do
      end do
    end if

#ifdef PAR
    call par_exchange_Sc_y(arr,btype)
#endif
    
    if (btype(So)==BC_DIRICHLET) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,0,k)  = side(So) - (arr(i,1,k)-side(So))
          arr(i,-1,k) = side(So) - (arr(i,2,k)-side(So))
        end do
      end do
    else if (btype(So)==BC_PERIODIC) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,0,k)  = arr(i,ny,k)
          arr(i,-1,k) = arr(i,ny-1,k)
        end do
      end do
    else if (btype(So)==BC_NEUMANN) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,0,k)  = arr(i,1,k) - side(So)*dyV(0)
          arr(i,-1,k) = arr(i,1,k) - side(So)*(dyV(0)+dyV(-1))
        end do
      end do
    else if (btype(So)==BC_CONSTFLUX) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,0,k)  = arr(i,1,k) + &
                side(So)*dyV(0)/((TDiff(i,1,k)+TDiff(i,0,k))/(2._knd))
          arr(i,-1,k) = arr(i,0,k) - (arr(i,1,k)-arr(i,0,k))
        end do
      end do
    else if ((btype(So)<BC_MPI_BOUNDS_MIN) .or. (btype(So)>BC_MPI_BOUNDS_MAX) .and. &
             ((btype(So)<BC_DOMAIN_BOUNDS_MIN) .or. (btype(So)>BC_DOMAIN_BOUNDS_MAX))) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,0,k)  = arr(i,1,k)
          arr(i,-1,k) = arr(i,2,k)
        end do
      end do
    end if

    if (btype(No)==BC_DIRICHLET) then
      do k = 1,nz
        do i=-1,nx+2
           arr(i,ny+1,k) = side(No) - (arr(i,ny,k)-side(No))
           arr(i,ny+2,k) = side(No) - (arr(i,ny-1,k)-side(No))
         end do
       end do
    else if (btype(No)==BC_PERIODIC) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,ny+1,k) = arr(i,1,k)
          arr(i,ny+2,k) = arr(i,2,k)
        end do
      end do
    else if (btype(No)==BC_NEUMANN) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,ny+1,k) = arr(i,ny,k) + side(No)*dyV(ny+1)
          arr(i,ny+2,k) = arr(i,ny,k) + side(No)*(dyV(ny+1)+dyV(ny+2))
        end do
      end do
    else if (btype(No)==BC_CONSTFLUX) then
      do k = 1,nz
        do i=-1,nx+2
          arr(i,ny+1,k) = arr(i,ny+1,k) - &
              side(No)*dyV(ny)/((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2._knd))
          arr(i,ny+2,k) = arr(i,ny+1,k) - (arr(i,ny,k)-arr(i,ny+1,k))
        end do
      end do
    else if (((btype(No)<BC_MPI_BOUNDS_MIN) .or. (btype(No)>BC_MPI_BOUNDS_MAX)) .and. &
             ((btype(No)<BC_DOMAIN_BOUNDS_MIN) .or. (btype(No)>BC_DOMAIN_BOUNDS_MAX))) then    
      do k = 1,nz
        do i=-1,ny+2
          arr(i,ny+1,k) = arr(i,ny,k)
          arr(i,ny+2,k) = arr(i,ny-1,k)
        end do
      end do
    end if

#ifdef PAR
    call par_exchange_Sc_z(arr,btype)
#endif
    
    if (btype(Bo)==BC_DIRICHLET) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1) = side(Bo) - (arr(i,j,1)-side(Bo))
          arr(i,j,nz+2) = side(Bo) - (arr(i,j,2)-side(Bo))
        end do
      end do
    else if (btype(Bo)==BC_PERIODIC) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,0)  = arr(i,j,nz)
          arr(i,j,-1) = arr(i,j,nz-1)
        end do
      end do
    else if (btype(Bo)==BC_NEUMANN) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,0)  = arr(i,j,1) - side(Bo)*dzW(0)
          arr(i,j,-1) = arr(i,j,1) - side(Bo)*(dzW(0)+dzW(-1))
        end do
      end do
    else if (btype(Bo)==BC_CONSTFLUX) then
      do j=-1,ny+2
        do i=-1,nx+2
            arr(i,j,0) = arr(i,j,1) + &
                side(Bo)*dzW(0)/((TDiff(i,j,1)+TDiff(i,j,0))/(2._knd))
          arr(i,j,-1) = arr(i,j,0)-(arr(i,j,1)-arr(i,j,0))
        end do
      end do
    else if (((btype(Bo)<BC_MPI_BOUNDS_MIN) .or. (btype(Bo)>BC_MPI_BOUNDS_MAX)) .and. &
             ((btype(Bo)<BC_DOMAIN_BOUNDS_MIN) .or. (btype(Bo)>BC_DOMAIN_BOUNDS_MAX))) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,0)  = arr(i,j,1)
          arr(i,j,-1) = arr(i,j,2)
        end do
      end do
    end if

    if (btype(To)==BC_DIRICHLET) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1) = side(To) - (arr(i,j,nz)-side(To))
          arr(i,j,nz+2) = side(To) - (arr(i,j,nz-1)-side(To))
        end do
      end do
    else if (btype(To)==BC_PERIODIC) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1) = arr(i,j,1)
          arr(i,j,nz+2) = arr(i,j,2)
        end do
      end do
    else if (btype(To)==BC_NEUMANN) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1) = arr(i,j,nz) + side(To)*dzW(nz+1)
          arr(i,j,nz+2) = arr(i,j,nz) + side(To)*(dzW(nz+1)+dzW(nz+2))
        end do
      end do
    else if (btype(To)==BC_CONSTFLUX.or.btype(To)==BC_AUTOMATIC_FLUX) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1) = arr(i,j,nz) - &
              side(To)*dzW(nz)/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._knd))
          arr(i,j,nz+2) = arr(i,j,nz+1) - (arr(i,j,nz)-arr(i,j,nz+1))
        end do
      end do
    else if (((btype(To)<BC_MPI_BOUNDS_MIN) .or. (btype(To)>BC_MPI_BOUNDS_MAX)) .and. &
             ((btype(To)<BC_DOMAIN_BOUNDS_MIN) .or. (btype(To)>BC_DOMAIN_BOUNDS_MAX))) then
      do j=-1,ny+2
        do i=-1,nx+2
          arr(i,j,nz+1) = arr(i,j,nz)
          arr(i,j,nz+2) = arr(i,j,nz-1)
        end do
      end do
     end if
  end subroutine CommonBase




  subroutine BoundScalar(SCAL)
    real(knd),intent(inout) :: SCAL(-1:,-1:,-1:)

    call CommonBase(SCAL,ScalBtype,sideScal)

  end subroutine BoundScalar






  subroutine BoundTemperature(Temperature)
#ifdef PAR
    use domains_bc_par
#endif
    real(knd),intent(inout) :: Temperature(-1:,-1:,-1:)
    
#ifdef PAR
    call par_update_domain_bounds_temperature(Temperature, time_stepping%effective_time)
#endif

    call CommonBase(Temperature,TempBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr)

  end subroutine BoundTemperature


  subroutine BoundMoisture(Moisture)
#ifdef PAR
    use domains_bc_par
#endif
    real(knd),intent(inout) :: Moisture(-1:,-1:,-1:)

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
    call par_exchange_Sc_x(Nu, Btype)
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
    call par_exchange_Sc_y(Nu, Btype)
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
    call par_exchange_Sc_z(Nu, Btype)
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
