module ImmersedBoundaryWM
  use Parameters
  use SolidBodies
  use GeometricShapes, only: Terrain, Ray
  use WallModels
  
  implicit none

  private
  public GetSolidBodiesWM, GetSolidBodiesWM_UVW, GetWMFluxes
 
contains
  
  subroutine GetSolidBodiesWM
    type(WMPoint)            :: p
    type(SolidBody), pointer :: SB
    integer                  :: neighbours(3,6)
    real(knd)     :: dist,nearx,neary,nearz
    integer       :: i,j,k,m,n,o,r
    integer       :: nb    
    real(knd) :: min_wall_distance
    
    !This must not be stricter than the threshold in FindInsideCells.
    min_wall_distance = (dxmin*dymin*dzmin)**(1._knd/3)/20

    if (wallmodeltype==0) return

    allocate(p%depscalar(num_of_scalars))

    !six triplets [1,0,0], [-1,0,0], [0,1,0],...
    neighbours = 0
    neighbours(1,1) =  1
    neighbours(1,2) = -1
    neighbours(2,3) =  1
    neighbours(2,4) = -1
    neighbours(3,5) =  1
    neighbours(3,6) = -1

    !Pobably not worth parallelizing for just one pass. It can speed it up,
    !but we need a specific order and we call the sort() in the end
    ! which would have more work to do.
    !Also, gfortran 4.8 has a bug which causes crashes when this loop is parallelized.

    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       if (Prtype(i,j,k)<=0) then
        dist = huge(dist)
        nb = 0
         
        do r=1,6
           m=neighbours(1,r)
           n=neighbours(2,r)
           o=neighbours(3,r)

           if ((Prtype(i+m,j+n,k+o)>0).and.Prtype(i+m,j+n,k+o)/=nb.and.(sum(abs([m,n,o]))==1)) then
             call SetCurrentSB(SB,Prtype(i+m,j+n,k+o))
             call SB%Closest(nearx,neary,nearz,xPr(i),yPr(j),zPr(k))
             if (norm2([nearx-xPr(i),neary-yPr(j),nearz-zPr(k)])<dist) then
              dist = norm2([nearx-xPr(i),neary-yPr(j),nearz-zPr(k)])
              nb = Prtype(i+m,j+n,k+o)
             end if
           end if
        end do

        if (nb>0) then
                    
          call SetCurrentSB(SB,nb)
          
          p%xi = i
          p%yj = j
          p%zk = k
          p%distx = nearx-xPr(i)
          p%disty = neary-yPr(j)
          p%distz = nearz-zPr(k)
          
          dist = norm2([p%distx,p%disty,p%distz])
          
          if (dist==0) then
            write(*,*) "Prtype(i,j,k)", Prtype(i,j,k)
            write(*,*) "i,j,k", i,j,k
            write(*,*) "xPr(i), yPr(j), zPr(k)", xPr(i), yPr(j), zPr(k)
            call SB%Closest(nearx,neary,nearz,xPr(i),yPr(j),zPr(k))
            write(*,*) "nearx, neary, nearz", nearx, neary, nearz
            call error_stop("zero distance WM point")
          end if
          
          if (dist<min_wall_distance) then
            write(*,*) "ijk", p%xi, p%yj, p%zk
            write(*,*) "dist dx, dy, dz", p%distx, p%disty, p%distz
            write(*,*) "grid dx, dy, dz", dxmin, dymin, dzmin
            call error_stop("Error, WM point is too close to the wall!")
          end if
                
          p%ustar = 1
          
          !HACK
          if (hypot(p%distx,p%disty)<abs(p%distz)*5) then
            !Santamouris - Environmental Design of Urban Buildings
            p%albedo = 0.3 !red brick 
            p%emissivity = 0.9 !red brick 
          else
            p%albedo = SB%albedo
            p%emissivity = SB%emissivity
          end if

          select type (geomshape => SB%GeometricShape)
            type is (Terrain)
              if (geomshape%PrPoints(i,j)%rough) then
                p%z0 = geomshape%PrPoints(i,j)%z0
              else
                p%z0 = 0
              end if
            class default
              if (SB%rough) then
                p%z0 = SB%z0
              else
                p%z0 = 0
              end if
          end select
          
          p%z0H = p%z0

          call AddWMPoint(p)
        end if
       end if
      end do
     end do
    end do
    
  end subroutine GetSolidBodiesWM
  
  
  
  
  subroutine GetSolidBodiesWM_UVW
    integer                  :: neighbours(3,MINUSX:PLUSZ)
    real(knd), target        :: r_neighbours(3,MINUSX:PLUSZ)
    real(knd) :: min_wall_distance
    
    !This must not be stricter than the threshold in FindInsideCells.
    min_wall_distance = (dxmin*dymin*dzmin)**(1._knd/3)/20

    if (wallmodeltype==0) return

    !six triplets [1,0,0], [-1,0,0], [0,1,0],...
    neighbours = 0
    neighbours(1, MINUSX) = -1
    neighbours(1, PLUSX)  =  1
    neighbours(2, MINUSY) = -1
    neighbours(2, PLUSY)  =  1
    neighbours(3, MINUSZ) = -1
    neighbours(3, PLUSZ)  =  1
    
    r_neighbours = real(neighbours, knd)
    
    call helper(1, Unx, Uny, Unz, xU(-2:), yPr, zPr, Utype)

    call helper(2, Vnx, Vny, Vnz, xPr, yV(-2:), zPr, Vtype)

    call helper(3, Wnx, Wny, Wnz, xPr, yPr, zW(-2:), Wtype)

  contains
  
    subroutine helper(component, nx, ny, nz, x, y, z, Xtype)
      integer, intent(in) :: component, nx, ny, nz
      real(knd), intent(in) :: x(-2:), y(-2:), z(-2:)
      integer, intent(in) :: Xtype(-2:,-2:,-2:)
      type(WMPointUVW) :: p
      type(SolidBody), pointer :: SB
      real(knd), pointer, contiguous :: dirvec(:)
      real(knd)     :: nearx, neary, nearz, t, dist, distvec(3)
      integer       :: i, j, k, m, n, o, dir
      
      !See above in GetSolidBodiesWM why no OpenMP.

      do k = 1, nz
       do j = 1, ny
        do i = 1, nx
        
          if (Xtype(i,j,k) < 0) then

            do dir = MINUSX, PLUSZ

              dirvec => r_neighbours(:,dir)

              m = neighbours(1,dir)
              n = neighbours(2,dir)
              o = neighbours(3,dir)
              
              if ((Xtype(i+m,j+n,k+o)>0)) then
              
                call SetCurrentSB(SB, Xtype(i+m,j+n,k+o))
                
                call SB%Closest(nearx ,neary, nearz, x(i), y(j), z(k))
                
                distvec = [nearx - x(i), neary - y(j), nearz - z(k)]
                
                dist = norm2(distvec)

                if (.not. right_direction(distvec, dirvec)) then
                  t = SB%ClosestOnLineOut( x(i+m), y(j+n), z(k+o), &
                                                   x(i),   y(j),   z(k) )
                  nearx = x(i+m) + t * ( x(i) - x(i+m) )
                  neary = y(j+n) + t * ( y(j) - y(j+n) )
                  nearz = z(k+o) + t * ( z(k) - z(k+o) )
                  distvec = [nearx - x(i), neary - y(j), nearz - z(k)]
                  if (norm2(distvec)<dist) then
                    !Potentially problematic
                    !write(*,'(a)',advance="no") '.'
                    if (norm2(distvec) < min_wall_distance) then
                      !Really problematic
                      write(*,*) "Warning, inconsistent processing of geometry detected in GetSolidBodiesWM_UVW,"
                      write(*,*) "point",i,j,k,"dir",dir
                      write(*,*) "distance on line out",norm2(distvec)
                      write(*,*) "distance",dist
                      write(*,*) "minimal allowed distance", min_wall_distance
                      call SB%Closest(nearx ,neary, nearz, x(i), y(j), z(k))
                  
                      distvec = [nearx - x(i), neary - y(j), nearz - z(k)]
                    end if
                  end if
                end if

                p%xi = i
                p%yj = j
                p%zk = k
                
                p%distx = distvec(1)
                p%disty = distvec(2)
                p%distz = distvec(3)
                
                p%ustar = 1
                
                select type (geomshape => SB%GeometricShape)
                  type is (Terrain)
                  
                    select case (component)
                      case (1)
                        if (geomshape%UPoints(i,j)%rough) then
                          p%z0 = geomshape%UPoints(i,j)%z0
                        else
                          p%z0 = 0
                        end if
                      case (2)
                        if (geomshape%VPoints(i,j)%rough) then
                          p%z0 = geomshape%VPoints(i,j)%z0
                        else
                          p%z0 = 0
                        end if
                      case (3)
                        if (geomshape%PrPoints(i,j)%rough) then
                          p%z0 = geomshape%PrPoints(i,j)%z0
                        else
                          p%z0 = 0
                        end if
                    end select
                    
                  class default
                  
                    if (SB%rough) then
                      p%z0 = SB%z0
                    else
                      p%z0 = 0
                    end if
                    
                end select
                
                if (norm2(distvec) < min_wall_distance) then
                  write(*,*) "component", component
                  write(*,*) "ijk", p%xi, p%yj, p%zk
                  write(*,*) "x, y, z", x(i),y(j),z(k)
                  write(*,*) "wall x, y, z", nearx, neary, nearz
                  write(*,*) "dist dx, dy, dz", distvec
                  write(*,*) "grid dx, dy, dz", dxmin, dymin, dzmin
                  call error_stop("Error, WM point UVW is too close to the wall!")
                end if
                
                p%z0H = p%z0
                
                call AddWMPointUVW(p, component, dir)
              end if

            end do
             
          end if

        end do
       end do
      end do

    end subroutine
  
    logical function right_direction(a, b)
      !checks if the angle between two vectors is small enough
      real(knd), intent(in) :: a(3), b(3)
      right_direction = ( dot_product(a, b) / (norm2(a) * norm2(b)) ) > 0.4_knd
    end function

  end subroutine GetSolidBodiesWM_UVW
  
  
  subroutine GetWMFluxes
    use SolarRadiation
    use PhysicalProperties
    use PlantBody_type
    
    real(knd) :: inc_radiation_flux, radiation_balance, angle_to_sun, &
                 total_heat_flux, sensible_heat_flux, latent_heat_flux
    real(knd) :: out_norm(3), xr(3), distv(3), svf
    
    type(Ray) :: r
    
    integer :: i
    
    enable_radiation = .true.
    
    if (enable_buoyancy .and. enable_radiation) then
      do i=1,size(WMPoints)
      
        associate(p => WMPoints(i))

          distv = [p%distx, p%disty, p%distz]

          out_norm = - distv / norm2(distv)

          angle_to_sun = max(0._knd, dot_product(out_norm, vector_to_sun))

          xr = [xPr(p%xi), yPr(p%yj), zPr(p%zk)] + distv * 0.9_knd

          if (angle_to_sun>0._knd) then
            r = Ray(xr, vector_to_sun)
            if ((SolidBodiesList%any(DoIntersect)) .or. &
                (PlantBodiesList%any(DoIntersectPlants))) then
              angle_to_sun = 0
            end if
          end if

          !temporary
          svf = sky_view_factor(xr)

          inc_radiation_flux = angle_to_sun * solar_direct_flux() + &
                               solar_diffuse_flux()*svf + &
                               in_lw_radiation()*svf

          radiation_balance = inc_radiation_flux * (1-p%albedo) - &
                              out_lw_radiation(p%emissivity, temperature_ref) * svf

          total_heat_flux = radiation_balance - &
                            (0.1 * radiation_balance) !crude guess of the storage flux

          if (enable_moisture) then
            latent_heat_flux = total_heat_flux * p%evaporative_fraction
            sensible_heat_flux = total_heat_flux - latent_heat_flux
            p%moisture_flux = latent_heat_flux / (rho_air_ref * Lv_water_ref)
          else
            sensible_heat_flux = total_heat_flux
          end if

          p%temperature_flux = sensible_heat_flux / (rho_air_ref * Cp_air_ref)

        end associate  
        
      end do
      
      call SaveFluxes
    end if
    
    
    contains
    
      real(knd) function sky_view_factor(xyz)
        real(knd),intent(in) :: xyz(3)
        integer :: i,nfree

        
        nfree = 0
        
        do i=1,svf_nrays
          r = Ray(xyz, svf_vecs(:,i))
        
          if (.not.((SolidBodiesList%any(DoIntersect)) .or. &
                    (PlantBodiesList%any(DoIntersectPlants)))) then
            nfree = nfree + 1
          end if
        end do
        
        sky_view_factor = real(nfree,knd) / svf_nrays
        
      end function
  
      logical function DoIntersect(item) result(res)
        type(SolidBody), intent(in) :: item

        res = item%IntersectsRay(r)
      end function
      
      logical function DoIntersectPlants(item) result(res)
        type(PlantBody), intent(in) :: item

        res = item%IntersectsRay(r)
      end function
     
      subroutine SaveFluxes
        use VTKArray
        real(knd),allocatable :: temperature_flux(:,:,:)
        allocate(temperature_flux(1:Prnx,1:Prny,1:Prnz))
        
        temperature_flux = 0
        
        do i=1,size(WMPoints)
          associate(p => WMPoints(i))
            temperature_flux(p%xi,p%yj,p%zk) = p%temperature_flux
          end associate
        end do

!         call VtkArraySimple("tempfl.vtk",temperature_flux)
        
      end subroutine
      
  end subroutine GetWMFluxes

end module ImmersedBoundaryWM


















module IBPoint_types
  use Kinds
  
  implicit none

  type TInterpolationPoint
    integer   :: xi                !coordinates of the interpolation points
    integer   :: yj
    integer   :: zk
    real(knd) :: coef              !interpolation coefficients for the interpolation points
  endtype TInterpolationPoint

  type :: TVelIBPoint
    integer   :: component      !1..U, 2..V, 3..W
    integer   :: xi             !coordinates of the grid point
    integer   :: yj
    integer   :: zk
    real(knd) :: distx          !vector to the nearest boundary point
    real(knd) :: disty
    real(knd) :: distz
    integer   :: dirx           !integer form of the above vector (~sign(distx))
    integer   :: diry
    integer   :: dirz
    integer   :: interp         !kind of interpolation 0.. none (boundarypoint), 1..linear, 2..bilinear, 3..trilinear
    integer   :: interpdir      !direction of interpolation in the linear case, for bilinear it is normal direction to the interpolation plane
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
 !  contains
 !    procedure Create         => TVelIBPoint_Create
  end type TVelIBPoint



  type :: TScalFlIBPoint
    integer                :: xi                       !coordinates of the grid point
    integer                :: yj
    integer                :: zk
    real(knd)              :: dist                     !distance to the boundary
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
    integer                :: interp                   !kind of interpolation 1.. none (1 point outside), 2..linear, 4..bilinear  other values not allowed
    real(knd)              :: temperature_flux = 0      !desired temperature flux
    real(knd)              :: moisture_flux = 0      !desired temperature flux
    integer :: n_WMPs = 0 !number of associated wall model points feeding the  !  contains
 !    procedure Create         => TScalFlIBPoint_Create
  end type TScalFlIBPoint
end module IBPoint_types


module VelIBPoint_list
  use IBPoint_types
#define TYPEPARAM type(TVelIBPoint)
#include "list-inc-def.f90"
contains
#include "list-inc-proc.f90"
#undef TYPEPARAM
end module VelIBPoint_list


module ScalFlIBPoint_list
  use IBPoint_types
#define TYPEPARAM type(TScalFlIBPoint)
#include "list-inc-def.f90"
contains
#include "list-inc-proc.f90"
#undef TYPEPARAM
end module ScalFlIBPoint_list






module ImmersedBoundary

  use Parameters
  use IBPoint_types
  use VelIBPoint_list, only: VelIBPointList => list
  use ScalFlIBPoint_list, only: ScalFlIBPointList => list
  use SolidBodies
  use ImmersedBoundaryWM

  implicit none

  private

  public TIBPoint, TIBPoint_Interpolate, TIBPoint_Viscosity, &
         UIBPoints, VIBPoints, WIBPoints, ScalFlIBPoints, &
         GetSolidBodiesBC, InitIBPFluxes, &
         Scalar_ImmersedBoundaries
         !InitSolidBodies imported from SolidBodies




  interface Create
    module procedure TVelIBPoint_Create
    module procedure TScalFlIBPoint_Create
  end interface Create



  type TIBPoint
    integer   :: xi
    integer   :: yj
    integer   :: zk
    real(knd) :: dist
    real(knd) :: temperature_flux = 0
    real(knd) :: moisture_flux = 0
    integer   :: interp
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
    integer :: n_WMPs = 0 !number of associated wall model points feeding the fluxes
 !   contains
 !     procedure Interpolate      => TIBPoint_Interpolate
 !     procedure InterpolateTDiff => TIBPoint_Interpolate_TDiff
 !     procedure ScalFlSource     => TIBPoint_ScalFlSource
 !     procedure MomentumSource   => TIBPoint_MomentumSource
 !     procedure Viscosity        => TIBPoint_Viscosity
  end type TIBPoint

  type(VelIBPointList)  :: UIBPointsList, VIBPointsList, WIBPointsList
  type(ScalFlIBPointList)  :: ScalFlIBPointsList

  type(TIBPoint),dimension(:),allocatable,save :: UIBPoints, VIBPoints, WIBPoints
  type(TIBPoint),dimension(:),allocatable,save :: ScalFlIBPoints

  interface assignment (=)
    module procedure VelIBPtoIBP
  end interface

  interface assignment (=)
    module procedure ScalFlIBPtoIBP
  end interface
  
  interface IB_interpolation_coefs
    module procedure IB_interpolation_coefs_1st_order
  end interface

  !number of least square interpolation points
  !must be higher than the interpolation polynomial order
  integer, parameter :: n_ls = 2
   
contains

  subroutine  Scalar_ImmersedBoundaries(Scal)
    real(knd), contiguous :: Scal(-1:,-1:,-1:)
    integer :: i

    !$omp parallel do
    do i = 1, size(ScalFlIBPoints)
      Scal(ScalFlIBPoints(i)%xi, ScalFlIBPoints(i)%yj, ScalFlIBPoints(i)%zk) = &
        TIBPoint_Interpolate(ScalFlIBPoints(i), Scal, -1)
    end do
    !$omp end parallel do
  end subroutine

  subroutine VelIBPtoIBP(IBP,VelIBP)
    type(TIBPoint),intent(out)    :: IBP
    type(TVelIBPoint),intent(in)  :: VelIBP

    IBP%xi = VelIBP%xi
    IBP%yj = VelIBP%yj
    IBP%zk = VelIBP%zk
    IBP%interp = VelIBP%interp

    allocate(IBP%IntPoints(size(VelIBP%IntPoints)))

    IBP%IntPoints = VelIBP%IntPoints
  end subroutine VelIBPtoIBP

  subroutine ScalFlIBPtoIBP(IBP,ScalFlIBP)
    type(TIBPoint),intent(out)    :: IBP
    type(TScalFlIBPoint),intent(in)  :: ScalFlIBP

    IBP%xi = ScalFlIBP%xi
    IBP%yj = ScalFlIBP%yj
    IBP%zk = ScalFlIBP%zk
    IBP%dist = ScalFlIBP%dist
    IBP%temperature_flux = ScalFlIBP%temperature_flux
    IBP%moisture_flux = ScalFlIBP%moisture_flux
    IBP%interp = ScalFlIBP%interp

    allocate(IBP%IntPoints(size(ScalFlIBP%IntPoints)))

    IBP%IntPoints = ScalFlIBP%IntPoints
    
    IBP%n_WMPs = ScalFlIBP%n_WMPs
  end subroutine ScalFlIBPtoIBP


  
  pure recursive function TIBPoint_Interpolate(IBP,U,lb) result(Uint)
    real(knd) :: Uint
    type(TIBPoint),intent(in) :: IBP
    integer,intent(in) :: lb
    real(knd),dimension(lb:,lb:,lb:),intent(in) :: U
    integer :: i

    Uint = 0

    do i=1,IBP%interp
      Uint = Uint + IBP%IntPoints(i)%coef * U(IBP%IntPoints(i)%xi,&
                                              IBP%IntPoints(i)%yj,&
                                              IBP%IntPoints(i)%zk)
    end do
  end function TIBPoint_Interpolate


  pure recursive function TIBPoint_InterpolateTDiff(IBP,U) result(Uint)
    real(knd) :: Uint
    type(TIBPoint),intent(in) :: IBP
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: U
    integer :: i,n

    n = 0
    Uint = 0

    do i=1,IBP%interp
      if (abs(IBP%IntPoints(i)%coef-0._knd)>epsilon(1._knd)) then
        Uint = Uint + U(IBP%IntPoints(i)%xi,&
                        IBP%IntPoints(i)%yj,&
                        IBP%IntPoints(i)%zk)
        n = n + 1
      end if
    end do

    Uint = Uint + U(IBP%xi,IBP%yj,IBP%zk)
    Uint = Uint / (n+1)
  end function TIBPoint_InterpolateTDiff

  
  pure recursive function TIBPoint_Viscosity(IBP,Viscosity)  result(src)  !Virtual scalar source for the Immersed Boundary Method with prescribed scalar flux on the boundary
    real(knd) :: src

    type(TIBPoint),intent(in) :: IBP
    real(knd),intent(in)      :: Viscosity(-1:,-1:,-1:)

    src = TIBPoint_Interpolate(IBP,Viscosity,-1)

  end function TIBPoint_Viscosity




  recursive subroutine TVelIBPoint_Create(IBP,xi,yj,zk,xU,yU,zU,Utype,component)
    type(TVelIBPoint),intent(out)               :: IBP
    integer,intent(in)                          :: xi,yj,zk
    real(knd),dimension(-2:),intent(in)         :: xU,yU,zU
    integer,dimension(-2:,-2:,-2:),intent(in)   :: Utype
    integer,intent(in)                          :: component

    type(SolidBody),pointer :: SB
    integer :: dirx,diry,dirz
    real(knd) :: x,y,z
    real(knd) :: x2,y2,z2
    real(knd) :: xnear,ynear,znear
    real(knd) :: tx,ty,tz
    logical :: freexm,freeym,freezm
    logical :: freexp,freeyp,freezp
    logical :: free1xm,free1ym,free1zm
    logical :: free1xp,free1yp,free1zp
    integer :: i

    !real coordinates of the IB forcing point
    x = xU(xi)
    y = yU(yj)
    z = zU(zk)
    call SetCurrentSB(SB,Utype(xi,yj,zk))
    
    IBP%component = component
    
    !integer grid coordinates
    IBP%xi = xi
    IBP%yj = yj
    IBP%zk = zk
    
!HACK: last resort when encountering instabilities near boundaries
! call null_point; return

    !HACK need to make sure we are not interpolating from the boundary buffer at least in
    !some problematic cases which have to be researched and tested.
    block
      integer :: nx, ny, nz

      select case (component)
        case (1)
          nx = Unx
          ny = Uny
          nz = Unz
        case (2)
          nx = Vnx
          ny = Vny
          nz = Vnz
        case default
          nx = Wnx
          ny = Wny
          nz = Wnz
      end select
      
      if (Btype(We)>=BC_DOMAIN_BOUNDS_MIN .and. &
          Btype(We)<=BC_DOMAIN_BOUNDS_MAX .and. &
          xi<=1) then
        call null_point; return
      else if (Btype(Ea)>=BC_DOMAIN_BOUNDS_MIN .and. &
          Btype(Ea)<=BC_DOMAIN_BOUNDS_MAX .and. &
          xi>=nx) then
        call null_point; return
      else if (Btype(So)>=BC_DOMAIN_BOUNDS_MIN .and. &
          Btype(So)<=BC_DOMAIN_BOUNDS_MAX .and. &
          yj<=1) then
        call null_point; return
      else if (Btype(No)>=BC_DOMAIN_BOUNDS_MIN .and. &
          Btype(No)<=BC_DOMAIN_BOUNDS_MAX .and. &
          yj<=ny) then
        call null_point; return
      else if (Btype(Bo)>=BC_DOMAIN_BOUNDS_MIN .and. &
          Btype(Bo)<=BC_DOMAIN_BOUNDS_MAX .and. &
          zk<=1) then
        call null_point; return
      else if (Btype(To)>=BC_DOMAIN_BOUNDS_MIN .and. &
          Btype(To)<=BC_DOMAIN_BOUNDS_MAX .and. &
          zk>=nz) then
        call null_point; return
      end if

    end block


    if (.not. SB%Inside(x,y,z,0._knd)) then
      call null_point
      return
    end if

    call SB%ClosestOut(xnear,ynear,znear,x,y,z)

    if (norm2([x-xnear,y-ynear,z-znear])<=0.1_knd*(dxmin*dymin*dzmin)**(1._knd/3)) then
      call null_point
      return
    end if



    freexm = all([ ( Utype(xi-i,yj  ,zk  )<=0, i = 1, n_ls ) ])
    freeym = all([ ( Utype(xi  ,yj-i,zk  )<=0, i = 1, n_ls ) ])
    freezm = all([ ( Utype(xi  ,yj  ,zk-i)<=0, i = 1, n_ls ) ])
    freexp = all([ ( Utype(xi+i,yj  ,zk  )<=0, i = 1, n_ls ) ])
    freeyp = all([ ( Utype(xi  ,yj+i,zk  )<=0, i = 1, n_ls ) ])
    freezp = all([ ( Utype(xi  ,yj  ,zk+i)<=0, i = 1, n_ls ) ])
    
    free1xm = Utype(xi-1,yj  ,zk  )<=0
    free1ym = Utype(xi  ,yj-1,zk  )<=0
    free1zm = Utype(xi  ,yj  ,zk-1)<=0
    free1xp = Utype(xi+1,yj  ,zk  )<=0
    free1yp = Utype(xi  ,yj+1,zk  )<=0
    free1zp = Utype(xi  ,yj  ,zk+1)<=0
    
    if ((free1xm.and..not.freexm) .or. &
        (free1ym.and..not.freeym) .or. &
        (free1zm.and..not.freezm) .or. &
        (free1xp.and..not.freexp) .or. &
        (free1yp.and..not.freeyp) .or. &
        (free1zp.and..not.freezp)) then
      call null_point;  return
    end if
    
    if (count(Utype(xi-1:xi+1,yj-1:yj+1,zk-1:zk+1)>0)>21) then
      call null_point; return
    end if
    
    if ((free1xm.and.free1xp) .or. (free1ym.and.free1yp) .or. (free1zm.and.free1zp)) then
      call null_point; return
    end if  

    if ( (.not.freexm) .and. (.not.freexp) ) then
      dirx = 0
    else if (freexm) then
      dirx = -1
    else
      dirx = 1
    end if
      
    if ( (.not.freeym) .and. (.not.freeyp) ) then
      diry = 0
    else if (freeym) then
      diry = -1
    else
      diry = 1
    end if
      
    if ( (.not.freezm) .and. (.not.freezp) ) then
      dirz = 0
    else if (freezm) then
      dirz = -1
    else
      dirz = 1
    end if
      
    tx = 0
    ty = 0
    tz = 0
      
    if (dirx/=0) then
      x2 = xU(xi+dirx)
      tx = SB%ClosestOnLineOut(x,y,z,x2,y,z)
      
      if (tx>1 .or. tx<0.05) then
        dirx = 0
        tx = 0
      end if
      
      IBP%distx = tx * (xU(xi+dirx) - xU(xi))
    end if
      
    if (diry/=0) then
      y2 = yU(yj+diry)
      ty = SB%ClosestOnLineOut(x,y,z,x,y2,z)
      
      if (ty>1 .or. ty<0.05) then
        diry = 0
        ty = 0
      end if
      
      IBP%disty = ty * (yU(yj+diry) - yU(yj))
    end if
      
    if (dirz/=0) then
      z2 = zU(zk+dirz)
      tz = SB%ClosestOnLineOut(x,y,z,x,y,z2)
      
      if (tz>1 .or. tz<0.05) then
        dirz = 0
        tz = 0
      end if
      
      IBP%distz = tz * (zU(zk+dirz) - zU(zk))
    end if
    
    if (dirx==0.and.diry==0.and.dirz==0) then
      call null_point
      return
    end if
    
    IBP%dirx = dirx
    IBP%diry = diry
    IBP%dirz = dirz

    IBP%IntPoints = InterpolationPoints(IBP,xU,yU,zU)
    
     
    IBP%interp = size(IBP%IntPoints)
    
  contains
  
    subroutine null_point
      IBP%interp = 0
      allocate(IBP%IntPoints(0))
    end subroutine
    
  end subroutine TVelIBPoint_Create
  
  
  recursive function InterpolationPoints(IBP,xU,yU,zU) result(res)
    type(TInterpolationPoint),allocatable :: res(:)
    type(TVelIBpoint), intent(in)         :: IBP
    real(knd),dimension(-2:), intent(in)  :: xU,yU,zU
    
    type(TInterpolationPoint),target      :: tmp(3*n_ls)
    type(TInterpolationPoint), pointer    :: t(:)
    real(knd) :: xr, yr, zr, x(0:n_ls), y(0:n_ls), z(0:n_ls)
    real(knd) :: b(3), d(3), c
    integer :: i,xi,yj,zk,dirx,diry,dirz
    
!     interface IB_interpolation_coefs
!       module procedure IB_interpolation_coefs_1st_order
!     end interface
    

    xi = IBP%xi
    yj = IBP%yj
    zk = IBP%zk
    dirx = IBP%dirx
    diry = IBP%diry
    dirz = IBP%dirz
    
    !distance from the wall
    d = 1
    
    if (dirx/=0) then
     
      do i = 1, n_ls
        tmp(i)%xi = xi + i*dirx
        tmp(i)%yj = yj
        tmp(i)%zk = zk
      end do
      
      xr = xU(xi) + IBP%distx
      x  = [ ( xU(xi+i*dirx), i = 0, n_ls ) ]
      tmp(1:n_ls)%coef = IB_interpolation_coefs(xr,x)
      
      d(1)  = abs(x(0)-xr)
        
    end if

    if (diry/=0) then

      do i = 1, n_ls
        tmp(n_ls + i)%xi = xi
        tmp(n_ls + i)%yj = yj + i*diry
        tmp(n_ls + i)%zk = zk
      end do
      
      yr = yU(yj) + IBP%disty
      y  = [ ( yU(yj+i*diry) , i = 0, n_ls ) ]
      tmp(n_ls+1:2*n_ls)%coef = IB_interpolation_coefs(yr,y)
      
      d(2) = abs(y(0)-yr)
      
    end if
    
    if (dirz/=0) then
     
      do i = 1, n_ls
        tmp(2*n_ls + i)%xi = xi
        tmp(2*n_ls + i)%yj = yj
        tmp(2*n_ls + i)%zk = zk + i*dirz
      end do

      zr = zU(zk) + IBP%distz
      z  = [ ( zU(zk+i*dirz), i = 0, n_ls )  ]
      tmp(2*n_ls+1:3*n_ls)%coef = IB_interpolation_coefs(zr,z)
      
      d(3) = abs(z(0)-zr)
      
    end if

    !beta 1,2,3 in eq. 17-19 in Peller et al., doi:1.1002/fld.1227
    b = 0
    
    if (dirx/=0) then
      b(1) = product(d)/(d(1))**2
    end if
    
    if (diry/=0) then
      b(2) = product(d)/(d(2))**2
    end if
    
    if (dirz/=0) then
      b(3) = product(d)/(d(3))**2
    end if
    
    !sum of betas
    c = sum(b)

    allocate(res(0))
    
    !NOTE: when well supported by compilers use associate (cf. http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56386)
    if (dirx/=0) then
      t=>tmp(1:n_ls)
      t%coef = t%coef * b(1)/c
      res = [res, t]
    end if

    if (diry/=0) then
      t=>tmp(n_ls+1:2*n_ls)
      t%coef = t%coef * b(2)/c
      res = [res, t]
    end if

    if (dirz/=0) then
      t=>tmp(2*n_ls+1:3*n_ls)
      t%coef = t%coef * b(3)/c
      res = [res, t]
    end if

  end function InterpolationPoints




  pure function IB_interpolation_coefs_2nd_order(xr,x) result(Coefs)
    real(knd),intent(in)  :: xr,x(0:)
    real(knd) :: Coefs(size(x)-1)
    real(knd) :: A1, A2, A4
    integer   :: n

    n = size(x)-1
    if ( n<3 .or. size(x)<=n ) then
      Coefs = 0*x(1:)
    else
      A1 = sum( x(1:) - xr )**2
      A2 = sum( (x(1:)**2 - xr**2) * (x(1:) - xr) )
      A4 = sum( x(1:)**2 - xr**2 )**2
      
      Coefs = ( ( A2 * (x(1:)**2 - xr**2) - A4 * (x(1:) - xr) ) * (x(0)    - xr)&
            +   ( A2 * (x(1:) - xr) - A1 * (x(1:)**2 - xr**2) ) * (x(0)**2 - xr**2) )&
            / (A2**2 - A1*A4)
    end if

  end function


  function IB_interpolation_coefs_1st_order(xr,x) result(Coefs)
    real(knd),intent(in)  :: xr,x(0:)
    real(knd) :: Coefs(size(x)-1)
    real(knd) :: A
    integer   :: n

    n = size(x)-1
    if ( n<2 .or. size(x)<=n ) then
      Coefs = 0*x(1:)
    else
      A = sum( (x(1:) - xr)**2 )

      Coefs = (x(1:) - xr) * (x(0) - xr) / A
    end if

  end function


  subroutine TScalFlIBPoint_Create(IBP,xi,yj,zk)
    type(TScalFlIBPoint),intent(out) :: IBP
    integer,intent(in) :: xi,yj,zk               !grid coordinates of the forcing point
    type(SolidBody),pointer :: SB
    integer :: dirx,diry,dirz,dirx2,diry2,dirz2,nfreedirs,ndirs,i
    real(knd) :: x,y,z,xnear,ynear,znear,distx,disty,distz,t,tx,ty,tz
    logical freep00,free0p0,free00p,freem00,free0m0,free00m

    x = xPr(xi)                                   !physical coordinates of the forcing point
    y = yPr(yj)
    z = zPr(zk)
    call SetCurrentSB(SB,Prtype(xi,yj,zk))
    call SB%ClosestOut(xnear,ynear,znear,x,y,z)

    IBP%xi = xi
    IBP%yj = yj
    IBP%zk = zk

    distx = xnear-x                               !distance to the nearest point on the boundary
    disty = ynear-y
    distz = znear-z
    dirx = nint(sign(1.0_knd,distx))              !direction to the boundary point
    diry = nint(sign(1.0_knd,disty))
    dirz = nint(sign(1.0_knd,distz))


    if (abs(distx)<(xPr(xi+1)-xPr(xi-1))/100._knd) then  !If very close to the boundary, set the boundary here
      distx = 0
      dirx = 0
    end if

    if (abs(disty)<(yPr(yj+1)-yPr(yj-1))/100._knd) then
      disty = 0
      diry = 0
    end if

    if (abs(distz)<(zW(zk+1)-zW(zk-1))/100._knd) then
      distz = 0
      dirz = 0
    end if


    freep00 = (Prtype(xi+1,yj,zk)<=0)   !logicals denoting if the cell in plus x direction is free of SB
    freem00 = (Prtype(xi-1,yj,zk)<=0)
    free0p0 = (Prtype(xi,yj+1,zk)<=0)
    free0m0 = (Prtype(xi,yj-1,zk)<=0)
    free00p = (Prtype(xi,yj,zk+1)<=0)
    free00m = (Prtype(xi,yj,zk-1)<=0)



    if (.not.(freep00.or.freem00)) then
      dirx = 0
      distx = 0
    end if
    if (.not.(free0p0.or.free0m0)) then
      diry = 0
      disty = 0
    end if
    if (.not.(free00p.or.free00m)) then
      dirz = 0
      distz = 0
    end if

    nfreedirs = 0
    ndirs = 0

    if (freep00) nfreedirs = nfreedirs+1
    if (free0p0) nfreedirs = nfreedirs+1
    if (free00p) nfreedirs = nfreedirs+1
    if (freem00) nfreedirs = nfreedirs+1
    if (free0m0) nfreedirs = nfreedirs+1
    if (free00m) nfreedirs = nfreedirs+1

    if (dirx/=0) ndirs = ndirs+1
    if (diry/=0) ndirs = ndirs+1
    if (dirz/=0) ndirs = ndirs+1

    if (nfreedirs>ndirs.or.(freep00.and.freem00).or.(free0p0.and.free0m0).or.(free00p.and.free00m)) then
      IBP%interp = 0    !If more free spaces than directions, or free space in both oposite directions,
      distx = 0         !Assume is an unresolvably small feature and treat as on the boundary.
      disty = 0
      distz = 0
      dirx = 0
      diry = 0
      dirz = 0
      dirx = 0
      diry = 0
      dirz = 0
    end if


    if (dirx==0.and.diry==0.and.dirz==0) then   !If dir>0 treat as on the boundary.

      dirx2 = 0
      diry2 = 0
      dirz2 = 0
      if (freep00) dirx2 = dirx2+1     !Find a free neighbouring cell and interpolate from there.
      if (freem00) dirx2 = dirx2-1
      if (free0p0) diry2 = diry2+1
      if (free0m0) diry2 = diry2-1
      if (free00p) dirz2 = dirz2+1
      if (free00m) dirz2 = dirz2-1

      if (dirx2==0.and.diry2==0.and.dirz2==0) then  !probably a very thin body, just average free points around

        IBP%interp = nfreedirs

        allocate(IBP%IntPoints(nfreedirs))

        i = 1
        IBP%dist = 0
        if (freep00) then
          IBP%IntPoints(i) = TInterpolationPoint(xi + 1, yj, zk, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + xPr(xi+1)-xPr(xi)
        end if
        if (freem00) then
          IBP%IntPoints(i) = TInterpolationPoint(xi - 1, yj, zk, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + xPr(xi)-xPr(xi-1)
        end if
        if (free0p0) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj + 1, zk, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + yPr(yj+1)-yPr(yj)
        end if
        if (free0m0) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj - 1, zk, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + yPr(yj)-yPr(yj-1)
        end if
        if (free00p) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj, zk + 1, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + zPr(zk+1)-zPr(zk)
        end if
        if (free00m) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj, zk - 1, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + zPr(zk)-zPr(zk-1)
        end if
        IBP%dist = IBP%dist / nfreedirs

      else  !some edge

        allocate(IBP%IntPoints(1))

        IBP%IntPoints(1)%xi = IBP%xi+dirx2
        IBP%IntPoints(1)%yj = IBP%yj+diry2
        IBP%IntPoints(1)%zk = IBP%zk+dirz2
        IBP%IntPoints%coef = 1._knd
        IBP%interp = 1
        IBP%dist = sqrt((x-xPr(IBP%xi+dirx2))**2+(y-yPr(IBP%yj+diry2))**2+(z-zPr(IBP%zk+dirz2))**2)

      end if

    elseif (nfreedirs==1.or.ndirs==1) then     !Only one free point outside, interpolate from there.

      allocate(IBP%IntPoints(1))

      IBP%IntPoints(1)%xi = IBP%xi+dirx
      IBP%IntPoints(1)%yj = IBP%yj+diry
      IBP%IntPoints(1)%zk = IBP%zk+dirz
      IBP%IntPoints%coef = 1._knd
      IBP%interp = 1
      IBP%dist = sqrt((x-xPr(IBP%xi+dirx))**2+(y-yPr(IBP%yj+diry))**2+(z-zPr(IBP%zk+dirz))**2)

    elseif (nfreedirs==2.or.ndirs==2) then    !Two free directions, use bilinear interpolation in the plane contaning these two neigbours.

      allocate(IBP%IntPoints(2))

      if (dirx==0) then     !plane yz

        if (abs(disty/distz)<abs(yPr(yj+diry)-y)/abs(zPr(zk+dirz)-z)) then  !Which gridline does the line from this point
          t = (zPr(zk+dirz)-z)/distz                                           !to the boundary intersect? Then use the points
          IBP%IntPoints(1)%xi = IBP%xi
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(disty*t)/abs(yPr(IBP%yj+diry)-yPr(IBP%yj))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk+dirz
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((disty*t)**2+(distz*t)**2)
        else
          t = (yPr(yj+diry)-y)/disty
          IBP%IntPoints(1)%xi = IBP%xi
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(distz*t)/abs(zPr(IBP%zk+dirz)-zPr(IBP%zk))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj+diry
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((disty*t)**2+(distz*t)**2)
        end if

      elseif (diry==0) then     !plane xz

        if (abs(distx/distz)<abs(xPr(xi+dirx)-x)/abs(zPr(zk+dirz)-z)) then
          t = (zPr(zk+dirz)-z)/distz
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(distx*t)/abs(xPr(IBP%xi+dirx)-xPr(IBP%xi))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk+dirz
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(distz*t)**2)
        else
          t = (xPr(xi+dirx)-x)/distx
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(distz*t)/abs(zPr(IBP%zk+dirz)-zPr(IBP%zk))
          IBP%IntPoints(2)%xi = IBP%xi+dirx
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(distz*t)**2)
        end if

      else                 !plane xy

        if (abs(distx/disty)<abs(xPr(xi+dirx)-x)/abs(yPr(yj+diry)-y)) then
          t = (yPr(yj+diry)-y)/disty
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk
          IBP%IntPoints(1)%coef = abs(distx*t)/abs(xPr(IBP%xi+dirx)-xPr(IBP%xi))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj+diry
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(disty*t)**2)
        else
          t = (xPr(xi+dirx)-x)/distx
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk
          IBP%IntPoints(1)%coef = abs(disty*t)/abs(yPr(IBP%yj+diry)-yPr(IBP%yj))
          IBP%IntPoints(2)%xi = IBP%xi+dirx
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(disty*t)**2)
        end if
      end if

    else                         !more than two free directions

      tx = (xPr(xi+dirx)-x)/distx   !a vector pointing from a free point to this forcing point.
      ty = (yPr(yj+diry)-y)/disty
      tz = (zPr(zk+dirz)-z)/distz

      IBP%interp = 4

      allocate(IBP%IntPoints(4))

      if (tx<=ty.and.tx<=tz) then
        !coordinates of interpolation point are therefore
        !xPr(xi+dirx) = x+tx*distx
        !y+tx*disty
        !z+tx*distz
        IBP%dist = sqrt((distx*tx)**2+(disty*tx)**2+(distz*tx)**2)
        IBP%IntPoints(1)%xi = IBP%xi+dirx
        IBP%IntPoints(1)%yj = IBP%yj+diry
        IBP%IntPoints(1)%zk = IBP%zk+dirz
        IBP%IntPoints(1)%coef = abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                                  abs(distz*tx)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(2)%xi = IBP%xi+dirx
        IBP%IntPoints(2)%yj = IBP%yj
        IBP%IntPoints(2)%zk = IBP%zk+dirz
        IBP%IntPoints(2)%coef = (1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                                  abs(distz*tx)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(3)%xi = IBP%xi+dirx
        IBP%IntPoints(3)%yj = IBP%yj+diry
        IBP%IntPoints(3)%zk = IBP%zk
        IBP%IntPoints(3)%coef = abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                                  (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
        IBP%IntPoints(4)%xi = IBP%xi+dirx
        IBP%IntPoints(4)%yj = IBP%yj
        IBP%IntPoints(4)%zk = IBP%zk
        IBP%IntPoints(4)%coef = (1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                                  (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
      elseif (ty<=tx.and.ty<=tz) then
        !the same with ty
        IBP%dist = sqrt((distx*ty)**2+(disty*ty)**2+(distz*ty)**2)
        IBP%IntPoints(1)%xi = IBP%xi+dirx
        IBP%IntPoints(1)%yj = IBP%yj+diry
        IBP%IntPoints(1)%zk = IBP%zk+dirz
        IBP%IntPoints(1)%coef = abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                                  abs(distz*ty)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(2)%xi = IBP%xi
        IBP%IntPoints(2)%yj = IBP%yj+diry
        IBP%IntPoints(2)%zk = IBP%zk+dirz
        IBP%IntPoints(2)%coef = (1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                                  abs(distz*ty)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(3)%xi = IBP%xi+dirx
        IBP%IntPoints(3)%yj = IBP%yj+diry
        IBP%IntPoints(3)%zk = IBP%zk
        IBP%IntPoints(3)%coef = abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                                  (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
        IBP%IntPoints(4)%xi = IBP%xi
        IBP%IntPoints(4)%yj = IBP%yj+diry
        IBP%IntPoints(4)%zk = IBP%zk
        IBP%IntPoints(4)%coef = (1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                                  (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
      else
        !the same with tz
        IBP%dist = sqrt((distx*tz)**2+(disty*tz)**2+(distz*tz)**2)
        IBP%IntPoints(1)%xi = IBP%xi+dirx
        IBP%IntPoints(1)%yj = IBP%yj+diry
        IBP%IntPoints(1)%zk = IBP%zk+dirz
        IBP%IntPoints(1)%coef = abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                                  abs(disty*tz)/abs(yPr(yj+diry)-y)
        IBP%IntPoints(2)%xi = IBP%xi
        IBP%IntPoints(2)%yj = IBP%yj+diry
        IBP%IntPoints(2)%zk = IBP%zk+dirz
        IBP%IntPoints(2)%coef = (1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                                  abs(disty*tz)/abs(yPr(yj+diry)-y)
        IBP%IntPoints(3)%xi = IBP%xi+dirx
        IBP%IntPoints(3)%yj = IBP%yj
        IBP%IntPoints(3)%zk = IBP%zk+dirz
        IBP%IntPoints(3)%coef = abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                                  (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
        IBP%IntPoints(4)%xi = IBP%xi
        IBP%IntPoints(4)%yj = IBP%yj
        IBP%IntPoints(4)%zk = IBP%zk+dirz
        IBP%IntPoints(4)%coef = (1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                                  (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
      end if

    end if

    IBP%temperature_flux = SB%temperature_flux

  end subroutine TScalFlIBPoint_Create


  subroutine MoveIBPointsToArray
    !It would be posible call a generic procedure with the list and array
    ! as an argument, but we want the final arrays not polymorphic.
    integer :: iU, iV, iW, iS

    !$omp parallel sections

    !$omp section
    allocate(UIBPoints(UIBPointsList%Len()))
    iU = 0
    call UIBPointsList%for_each(CopyUIBPoint)

    !$omp section
    allocate(VIBPoints(VIBPointsList%Len()))
    iV = 0
    call VIBPointsList%for_each(CopyVIBPoint)

    !$omp section
    allocate(WIBPoints(WIBPointsList%Len()))
    iW = 0
    call WIBPointsList%for_each(CopyWIBPoint)

    !$omp section
    allocate(ScalFlIBPoints(ScalFlIBPointsList%Len()))
    iS = 0
    call ScalFlIBPointsList%for_each(CopyScalFlIBPoint)

    !$omp end parallel sections

  contains

    subroutine CopyUIBPoint(CurrentIBPoint)
      type(TVelIBPoint) :: CurrentIBPoint

      iU = iU + 1
      UIBPoints(iU) = CurrentIBPoint
    end subroutine

    subroutine CopyVIBPoint(CurrentIBPoint)
      type(TVelIBPoint) :: CurrentIBPoint

      iV = iV + 1
      VIBPoints(iV) = CurrentIBPoint
    end subroutine

    subroutine CopyWIBPoint(CurrentIBPoint)
      type(TVelIBPoint) :: CurrentIBPoint

      iW = iW + 1
      WIBPoints(iW) = CurrentIBPoint
    end subroutine

    subroutine CopyScalFlIBPoint(CurrentIBPoint)
      type(TScalFlIBPoint) :: CurrentIBPoint

      iS = iS + 1
      ScalFlIBPoints(iS) = CurrentIBPoint
    end subroutine

  end subroutine MoveIBPointsToArray


  subroutine AuxNeighbours(Xtype,nx,ny,nz)
    !helper procedure to FindNeighbouringCells
    integer,intent(in)    :: nx,ny,nz
    integer,intent(inout) :: Xtype(-2:,-2:,-2:)
    integer :: i,j,k

    do k = -1, nz+2
      do j = -1, ny+2
        do i = -1, nx+2
          if (Xtype(i,j,k)==0 .and. (any(Xtype(i-1:i+1,j-1:j+1,k-1:k+1)>0))) then
                 Xtype(i,j,k) = -1
          end if
        end do
      end do
    end do

  end subroutine

  subroutine FindNeighbouringCells
    !sets type of the cells closest to the solid body to -1

    !$omp parallel sections
    !$omp section
    call AuxNeighbours(Prtype,Prnx,Prny,Prnz)
    !$omp section
    call AuxNeighbours(Utype,Unx,Uny,Unz)
    !$omp section
    call AuxNeighbours(Vtype,Vnx,Vny,Vnz)
    !$omp section
    call AuxNeighbours(Wtype,Wnx,Wny,Wnz)
    !$omp end parallel sections
  end subroutine FindNeighbouringCells

  
  integer function FindScalarIBCellIndex(xi,yj,zk) result(res)
    !binary search of the immersed boundary cell index with the prescribed coordinates
    integer,intent(in) :: xi,yj,zk
    integer :: idx,d,u,i,idx_d,idx_u,idx_i
    
    idx = f_ijk(xi,yj,zk)
    
    d = 1
    u = size(ScalFlIBPoints)
    
    idx_d = f(d)
    if (idx_d==idx) then
      res = d
      return
    end if
    
    idx_u = f(u)
    if (idx_u==idx) then
      res = u
      return
    end if
    
    res = -1
    
    do while (u-d>1)
      i = (d+u)/2
      idx_i = f(i)
      if (idx_i == idx) then
        res = i
        return
      else if (idx_i<idx) then
        d = i
      else
        u = i
      end if
    end do
    
    contains
      
      pure integer function f_ijk(i,j,k)
        integer,intent(in) :: i,j,k
        f_ijk = (i-1) + Prnx * (j-1) + Prnx * Prny * (k-1)
      end function
      pure integer function f(n)
        integer,intent(in) :: n
        f = (ScalFlIBPoints(n)%xi-1) + &
            Prnx * (ScalFlIBPoints(n)%yj-1) + &
            Prnx * Prny * (ScalFlIBPoints(n)%zk-1)
      end function
  end function FindScalarIBCellIndex
  
!   
!   subroutine BindWMstoIBPs
!     use ArrayUtilities, only: add_element
!     use WallModels, only: WMPoints
!     integer :: neighbours(3,6)
!     integer :: i, j
!     integer :: xn, yn, zn, idx
!     
!     neighbours = 0
!     neighbours(1,1) =  1
!     neighbours(1,2) = -1
!     neighbours(2,3) =  1
!     neighbours(2,4) = -1
!     neighbours(3,5) =  1
!     neighbours(3,6) = -1
!     
!     do i=1,size(WMPoints)
!       associate (p => WMPoints(i)) !gfortran 4.8 bug prevents more items here
!         allocate(p%bound_IBPs(0))
! 
!         do j=1,6
!           xn = p%xi+neighbours(1,j)
!           yn = p%yj+neighbours(2,j)
!           zn = p%zk+neighbours(3,j)
!           if (Prtype(xn,yn,zn)>0.and. &
!               xn>0 .and. xn<=Prnx .and. &
!               yn>0 .and. yn<=Prny .and. &
!               zn>0 .and. zn<=Prnz) then
! 
!             idx = FindScalarIBCellIndex(xn, &
!                                         yn, &
!                                         zn)
!             if (idx>0) then
!               ScalFlIBPoints(idx)%n_WMPs = ScalFlIBPoints(idx)%n_WMPs + 1
!               call add_element(p%bound_IBPs, idx)
!             else
!               call error_stop("This point does not have an IBP!")
!             end if
!           end if
!         end do
!         
!       end associate
!     end do
!   end subroutine BindWMstoIBPs
  
  
  subroutine InitIBPFluxes
  
    !find the corresponding immersed boundary points to wall model points
!     call BindWMstoIBPs
    !compute initial temperature and moisture fluxes
    call GetWMFluxes
  end subroutine
  
!   subroutine SetIBPFluxes
!     use WallModels, only: WMPoints
!     integer :: i,j
!          
!     ScalFlIBPoints%temperature_flux = 0
! 
!     do i=1,size(WMPoints)
!       associate (p => WMPoints(i))
!         if (enable_buoyancy) then
!           do j=1,size(p%bound_IBPs)
!             associate (IBP => ScalFlIBPoints(p%bound_IBPs(j)))
!               IBP%temperature_flux = IBP%temperature_flux + p%temperature_flux/IBP%n_WMPs
!             end associate
!           end do
!         end if
!         
!         if (enable_moisture) then
!           do j=1,size(p%bound_IBPs)
!             associate (IBP => ScalFlIBPoints(p%bound_IBPs(j)))
!               IBP%moisture_flux = IBP%moisture_flux + p%moisture_flux/IBP%n_WMPs
!             end associate
!           end do
!         end if
!         
!         if (p%zk==1) then
!           if (enable_buoyancy) BsideTFlArr(p%xi,p%yj) = p%temperature_flux
!           if (enable_moisture) BsideMFlArr(p%xi,p%yj) = p%moisture_flux
!         end if
!       end associate
!     end do
!   
!   end subroutine SetIBPFluxes
!   
  


  subroutine InitImBoundaries
    type(TVelIBPoint) IBP
    type(TScalFlIBPoint) SIBP
    integer :: i,j,k

!strange runtime errors in gfortran 4.8, should not be a race condition
#if !( (defined __GFORTRAN__) && (__GNUC__==4) && (__GNUC_MINOR__<=8) )
    !$omp parallel private(i,j,k,IBP,SIBP)
#endif
    !$omp do collapse(3) schedule(dynamic, 1000) 
    do k = 1,Unz
     do j = 1,Uny
      do i = 1,Unx
       if (Utype(i,j,k)>0) then
        if (Utype(i+1,j,k)<=0 .or. Utype(i-1,j,k)<=0 .or. Utype(i,j+1,k)<=0 &
          .or. Utype(i,j-1,k)<=0 .or. Utype(i,j,k+1)<=0 .or. Utype(i,j,k-1)<=0)  then
            call  Create(IBP,i,j,k,xU(-2:),yPr(-2:),zPr(-2:),Utype,1)
            !$omp critical
            call  UIBPointsList%add(IBP)
            !$omp end critical
        end if
       end if
      end do
     end do
    end do
    !$omp end do nowait

    !$omp do collapse(3) schedule(dynamic, 1000) 
    do k = 1,Vnz
     do j = 1,Vny
      do i = 1,Vnx
       if (Vtype(i,j,k)>0) then
        if (Vtype(i+1,j,k)<=0 .or. Vtype(i-1,j,k)<=0 .or. Vtype(i,j+1,k)<=0 &
          .or. Vtype(i,j-1,k)<=0 .or. Vtype(i,j,k+1)<=0 .or. Vtype(i,j,k-1)<=0)  then
            call  Create(IBP,i,j,k,xPr(-2:),yV(-2:),zPr(-2:),Vtype,2)
            !$omp critical
            call  VIBPointsList%add(IBP)
            !$omp end critical
        end if
       end if
      end do
     end do
    end do
    !$omp end do nowait

    !$omp do collapse(3) schedule(dynamic, 1000) 
    do k = 1,Wnz
     do j = 1,Wny
      do i = 1,Wnx
       if (Wtype(i,j,k)>0) then
        if (Wtype(i+1,j,k)<=0 .or. Wtype(i-1,j,k)<=0 .or. Wtype(i,j+1,k)<=0 &
          .or. Wtype(i,j-1,k)<=0 .or. Wtype(i,j,k+1)<=0 .or. Wtype(i,j,k-1)<=0)  then
            call  Create(IBP,i,j,k,xPr(-2:),yPr(-2:),zW(-2:),Wtype,3)
            !$omp critical
            call  WIBPointsList%add(IBP)
            !$omp end critical
        end if
       end if
      end do
     end do
    end do
    !$omp end do nowait

    !$omp do collapse(3) schedule(dynamic, 1000) 
    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       if (Prtype(i,j,k)>0) then
        if (Prtype(i+1,j,k)<=0 .or. Prtype(i-1,j,k)<=0 .or. Prtype(i,j+1,k)<=0 &
          .or. Prtype(i,j-1,k)<=0 .or. Prtype(i,j,k+1)<=0 .or. Prtype(i,j,k-1)<=0)  then
            call  Create(SIBP,i,j,k)
            !$omp critical
            call  ScalFlIBPointsList%add(SIBP)
            !$omp end critical
        end if
       end if
      end do
     end do
    end do
    !$omp end do nowait
#if !( (defined __GFORTRAN__) && (__GNUC__==4) && (__GNUC_MINOR__<=8) )
    !$omp end parallel
#endif

    call MoveIBPointsToArray

    !$omp parallel sections
    !$omp section
    call UIBPointsList%finalize
    !$omp section
    call VIBPointsList%finalize
    !$omp section
    call WIBPointsList%finalize
    !$omp section
    call ScalFlIBPointsList%finalize
    !$omp end parallel sections

  end subroutine InitImBoundaries


  subroutine GetSolidBodiesBC
    use custom_par, only: par_sync_out
    
    call par_sync_out("      ...finding cells inside solid bodies.")
    call FindInsideCells

    call par_sync_out("      ...finding cells neighbouring solid bodies.")
    call FindNeighbouringCells

    call par_sync_out("      ...preparing immersed boundary interpolations.")
    call InitImBoundaries

    call par_sync_out("      ...preparing wall model scalar points.")
    call GetSolidBodiesWM

    call par_sync_out("      ...preparing wall model velocity points.")
    call GetSolidBodiesWM_UVW

  end subroutine GetSolidBodiesBC


end module ImmersedBoundary
