module VolumeSourceBody_type
  use Kinds
  use Body_class, only: Body
  
  implicit none

  type, extends(Body) :: VolumeSourceBody
    procedure(temperature_interface),nopass,pointer :: get_temperature_flux => null()
    procedure(temperature_interface),nopass,pointer :: get_moisture_flux    => null()
    procedure(scalar_flux_interface),nopass,pointer :: get_scalar_flux      => null()
    logical :: variable = .false.
  end type VolumeSourceBody

  abstract interface
    function temperature_interface(x,y,z) result(res)
      use Kinds
      real(knd) :: res
      real(knd),intent(in) :: x,y,z
    end function
    function scalar_flux_interface(x,y,z,num_of_scalar) result(res)
      use Kinds
      real(knd) :: res
      real(knd),intent(in) :: x,y,z
      integer,intent(in) :: num_of_scalar
    end function
  end interface
    

  
#define TYPEPARAM type(VolumeSourceBody)
#include "list-inc-def.f90"

  type(List) :: SourceBodiesList

contains
#include "list-inc-proc.f90"
#undef TYPEPARAM

end module



module PlantBody_type
  use Kinds
  use Body_class, only: Body
  
  implicit none

  type, extends(Body) :: PlantBody
    integer :: plant_type !problem specific, used by custom routines
    real(knd) :: albedo = 0.3_knd
    real(knd) :: emmissivity = 0.7_knd
    real(knd) :: evaporative_fraction = 0.6_knd
    procedure(resistance_interface),nopass,pointer  :: get_resistance       => null()
  end type PlantBody
  
  abstract interface
    function resistance_interface(x,y,z) result(res)
      use Kinds
      real(knd) :: res
      real(knd),intent(in) :: x,y,z
    end function
  end interface
  
#define TYPEPARAM type(PlantBody)
#include "list-inc-def.f90"

  type(List) :: PlantBodiesList
  

contains
#include "list-inc-proc.f90"
#undef TYPEPARAM

  
end module




module VolumeSources
  use Parameters
  use Lists
  use Body_class
  use VolumeSourceBody_type, only: VolumeSourceBody, SourceBodiesList, scalar_flux_interface
  use PlantBody_type, only: PlantBody, PlantBodiesList

  implicit none
! 
!   private
! 
!   public UResistanceVolumes, VResistanceVolumes, WResistanceVolumes, &
!          TemperatureFlVolumes, MoistureFlVolumes, ScalarFlVolumes, VolumeSourceBody

  type :: GridVolume
    integer   :: xi, yj, zk
  end type GridVolume

  type,extends(GridVolume) :: FluxVolume
    real(knd) :: flux
  end type FluxVolume
  
  interface FluxVolume
    module procedure FluxVolume_v3
  end interface

  type,extends(FluxVolume) :: TResistanceVolume
    ! f_drag_i = - Cd * a * V * u_i , V = sqrt(sum(u_i**2))
    !flux = Cd * a
  end type TResistanceVolume

  type,extends(FluxVolume) :: TemperatureFlVolume
    ! <T'w'> = temperature_flux, Qh = rho*Cp*temperature_flux
    !flux = temperature_flux
  end type TemperatureFlVolume

  type,extends(FluxVolume) :: MoistureFlVolume
    ! <q'w'> = moisture_flux, Qe = rho*Lv*moisture_flux
    !flux = moisture_flux
  end type MoistureFlVolume

  type,extends(FluxVolume) :: ScalarFlVolume
    ! <c'w'> = scalar_flux
    !flux = flux
    procedure(scalar_flux_interface),nopass,pointer :: get_flux => null()
  end type ScalarFlVolume

  interface ScalarFlVolume
    module procedure ScalarFlVolume_3i
    module procedure ScalarFlVolume_v3
  end interface

  type ScalarFlVolumesContainer
    integer :: scalar_number
    type(ScalarFlVolume), allocatable :: Volumes(:)
  end type

  interface Add
    module procedure AddScalarFlVolumesContainer
  end interface


  type(List) :: UResistanceVolumesList, &
                 VResistanceVolumesList, &
                 WResistanceVolumesList, &
                 TemperatureFlVolumesList, &
                 MoistureFlVolumesList
                 
  type ScalarFlVolumesListContainer
    type(List) :: list
  end type
                 
  type(ScalarFlVolumesListContainer),allocatable :: ScalarFlVolumesLists(:)

  !final volume momentum sinks - resistances
  type(TResistanceVolume),allocatable :: UResistanceVolumes(:), &
                                         VResistanceVolumes(:), &
                                         WResistanceVolumes(:)

  !final scalar quantities sources/sinks
  type(TemperatureFlVolume),allocatable :: TemperatureFlVolumes(:)
  type(MoistureFlVolume)   ,allocatable :: MoistureFlVolumes(:)
  !ScalarFlVolumes do not change in time, get_flux is not associated
  type(ScalarFlVolumesContainer),allocatable :: ScalarFlVolumes(:)
  !VariableScalarFlVolumes van change in time, get_flux is associated and called each time step
  type(ScalarFlVolumesContainer),allocatable :: VariableScalarFlVolumes(:)
  
  contains

    subroutine InitVolumeSourceBodies
#ifdef CUSTOMPB
      interface
        subroutine CustomVolumeSourceBodies
        end subroutine
      end interface

      call CustomVolumeSourceBodies
#endif
    end subroutine InitVolumeSourceBodies
    
    
    function FluxVolume_v3(pos,flux) result(res)
      type(FluxVolume) :: res
      integer,intent(in) :: pos(3)
      real(knd),intent(in) :: flux
      res%xi = pos(1)
      res%yj = pos(2)
      res%zk = pos(3)
      res%flux = flux
    end function

    
    function ScalarFlVolume_v3(pos,flux) result(res)
      type(ScalarFlVolume) :: res
      integer,intent(in) :: pos(3)
      real(knd),intent(in) :: flux
      res%FluxVolume = FluxVolume(pos,flux)
    end function

    
    function ScalarFlVolume_3i(x,y,z,flux) result(res)
      type(ScalarFlVolume) :: res
      integer,intent(in) :: x,y,z
      real(knd),intent(in) :: flux
      res%FluxVolume = FluxVolume([x,y,z],flux)
    end function

    
    subroutine InsideCellsToLists
      !find cells inside the canopy and store them in a list
      call PlantBodiesList%for_each(GetUCells)

      call PlantBodiesList%for_each(GetVCells)

      call PlantBodiesList%for_each(GetWCells)

      allocate(ScalarFlVolumesLists(num_of_scalars))

      call SourceBodiesList%for_each(GetOtherCells_VolumeSources)

      call PlantBodiesList%for_each(GetOtherCells_Plants)


      contains

        subroutine GetUCells(PB)
          type(PlantBody) :: PB
          type(TResistanceVolume) :: elem
          integer :: i,j,k

           if (associated(PB%get_resistance)) then
             do k = 0,Unz+1
              do j = 0,Uny+1
               do i = 0,Unx+1
                  if (Inside(PB,xU(i),yPr(j),zPr(k))) then
                    elem%xi = i
                    elem%yj = j
                    elem%zk = k
                    elem%flux = PB%get_resistance(xU(i),yPr(j),zPr(k))
                    call UResistanceVolumesList%add(elem)
                  end if
               enddo
              enddo
             enddo
           end if
        end subroutine

        subroutine GetVCells(PB)
          type(PlantBody) :: PB
          type(TResistanceVolume) :: elem
          integer :: i,j,k

          if (associated(PB%get_resistance)) then
            do k = 0,Vnz+1
             do j = 0,Vny+1
              do i = 0,Vnx+1
                 if (Inside(PB,xPr(i),yV(j),zPr(k))) then
                   elem%xi = i
                   elem%yj = j
                   elem%zk = k
                   elem%flux = PB%get_resistance(xPr(i),yV(j),zPr(k))
                   call VResistanceVolumesList%add(elem)
                 end if
              enddo
             enddo
            enddo
          end if
        end subroutine

        subroutine GetWCells(PB)
          type(PlantBody) :: PB
          type(TResistanceVolume) :: elem
          integer :: i,j,k

          if (associated(PB%get_resistance)) then
            do k = 0,Wnz+1
             do j = 0,Wny+1
              do i = 0,Wnx+1
                 if (Inside(PB,xPr(i),yPr(j),zW(k))) then
                   elem%xi = i
                   elem%yj = j
                   elem%zk = k
                   elem%flux = PB%get_resistance(xPr(i),yPr(j),zW(k))
                   call WResistanceVolumesList%add(elem)
                 end if
              enddo
             enddo
            enddo
          end if
        end subroutine

        subroutine GetOtherCells_VolumeSources(VSB)
          type(VolumeSourceBody) :: VSB
          type(TemperatureFlVolume) :: Telem
          type(MoistureFlVolume) :: Melem
          type(ScalarFlVolume) :: Selem
          integer :: i,j,k,sc
          
          do k = 0,Prnz+1
           do j = 0,Prny+1
            do i = 0,Prnx+1
               if (VSB%Inside(xPr(i),yPr(j),zPr(k))) then
                 Telem%xi = i
                 Telem%yj = j
                 Telem%zk = k
                 Melem%xi = i
                 Melem%yj = j
                 Melem%zk = k

                 Telem%flux = 0
                 Melem%flux = 0

                 
                 if (enable_buoyancy .and. associated(VSB%get_temperature_flux)) then
                 
                     Telem%flux = Telem%flux + &
                                  VSB%get_temperature_flux(xPr(i),yPr(j),zPr(k))
                 end if
                 
                 if (enable_moisture .and. associated(VSB%get_moisture_flux)) then

                     Melem%flux = Melem%flux + &
                                  VSB%get_moisture_flux(xPr(i),yPr(j),zPr(k))
                   
                 end if

                 if (Telem%flux/=0) call TemperatureFlVolumesList%add(Telem)

                 if (Melem%flux/=0) call MoistureFlVolumesList%add(Melem)
                  
                 
                 if (associated(VSB%get_scalar_flux)) then
                   do sc = 1,num_of_scalars
                     Selem%xi = i
                     Selem%yj = j
                     Selem%zk = k
                     
                     if (VSB%variable) then
                       Selem%get_flux => VSB%get_scalar_flux
                     else
                       Selem%get_flux => null()
                     end if
                     
                     Selem%flux = VSB%get_scalar_flux(xPr(i),yPr(j),zPr(k),sc)
                     
                     if (Selem%flux/=0) call ScalarFlVolumesLists(sc)%list%add(Selem)
                   end do
                 end if
               end if
            enddo
           enddo
          enddo
        end subroutine

        subroutine GetOtherCells_Plants(PB)
          use SolarRadiation, only: enable_radiation
          type(PlantBody) :: PB
          type(TemperatureFlVolume) :: Telem
          type(MoistureFlVolume) :: Melem
          integer :: i, j, k
          
          do k = 0,Prnz+1
           do j = 0,Prny+1
            do i = 0,Prnx+1
               if (PB%Inside(xPr(i),yPr(j),zPr(k))) then
               
                 Telem%xi = i
                 Telem%yj = j
                 Telem%zk = k
                 Melem%xi = i
                 Melem%yj = j
                 Melem%zk = k
                   
                 if (enable_buoyancy .and. enable_radiation) then
                   if (on_border(PB,i,j,k)) then
                     call GetRadiationFluxes(PB,Telem,Melem,xPr(i),yPr(j),zPr(k))
                   else
                     Telem%flux = 0
                     Melem%flux = 0
                   end if
                 else
                   Telem%flux = 0
                   Melem%flux = 0
                 end if

                 if (Telem%flux/=0) call TemperatureFlVolumesList%add(Telem)

                 if (Melem%flux/=0) call MoistureFlVolumesList%add(Melem)

               end if
            enddo
           enddo
          enddo
        end subroutine

        
        
        logical function on_border(PB,xi,yj,zk)
          type(PlantBody), intent(in) :: PB
          integer, intent(in) :: xi,yj,zk
          on_border = insPB(PB,xi,yj,zk) .and. &
                      .not. all([insPB(PB,xi-1,yj,zk), &
                                 insPB(PB,xi+1,yj,zk), &
                                 insPB(PB,xi,yj-1,zk), &
                                 insPB(PB,xi,yj+1,zk), &
                                 insPB(PB,xi,yj,zk-1), &
                                 insPB(PB,xi,yj,zk+1)])
        end function
        
        logical function insPB(PB,xi,yj,zk)
          type(PlantBody), intent(in) :: PB
          integer, intent(in) :: xi,yj,zk
          insPB = PB%Inside(xPr(xi),yPr(yj),zPr(zk))
        end function
    end subroutine InsideCellsToLists


    
    
    
    subroutine GetRadiationFluxes(PB,Telem,Melem,x,y,z)
      use SolarRadiation
      use PhysicalProperties
      use GeometricShapes, only: Ray
      use SolidBodies, only: SolidBodiesList, SolidBody
      class(PlantBody),intent(inout) :: PB
      type(TemperatureFlVolume),intent(inout) :: Telem
      type(MoistureFlVolume),intent(inout)    :: Melem
      real(knd),intent(in) :: x, y, z
      
      type(Ray) :: r
      real(knd) :: out_norm(3)
      real(knd) :: inc_radiation_flux, radiation_balance, angle_to_sun, &
                   total_heat_flux, sensible_heat_flux, latent_heat_flux
      real(knd) :: xr(3), svf

      out_norm = PB%OutwardNormal(x,y,z)

      if (norm2(out_norm)<epsilon(1._knd)) then
        write(*,*) "Error, zero outward normal at file",__FILE__,"line",__LINE__
        call error_stop
      end if

      angle_to_sun = max(0._knd, dot_product(out_norm, vector_to_sun))

      xr = [x,y,z] * out_norm * min(dxmin,dymin,dzmin) / 10
      
      if (angle_to_sun>0._knd) then
        r = Ray(xr, &
                vector_to_sun)

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

      radiation_balance = inc_radiation_flux * (1-PB%albedo) - &
                          out_lw_radiation(PB%emmissivity, temperature_ref) * svf

      total_heat_flux = radiation_balance - &
                        (0.1 * radiation_balance) !crude guess of the storage flux

      !FIXME: surface flux to volume flux APPROXIMATE!                  
      total_heat_flux = total_heat_flux/dot_product([dxmin,dymin,dzmin],abs(out_norm))
                        
      if (enable_moisture) then
        latent_heat_flux = total_heat_flux * PB%evaporative_fraction
        sensible_heat_flux = total_heat_flux - latent_heat_flux
        Melem%flux  = latent_heat_flux / (rho_air_ref * Lv_water_ref)
      else
        sensible_heat_flux = total_heat_flux
        Melem%flux  = 0
      end if
      
      Telem%flux = sensible_heat_flux / (rho_air_ref * Cp_air_ref)

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
      

    end subroutine GetRadiationFluxes

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    subroutine MovePointsToArray
    !It would be posible call a generic procedure with the list and array
    ! as an argument, but we want the final arrays not polymorphic.
      integer :: i, inorm, ivar, j, component

      allocate(UResistanceVolumes(UResistanceVolumesList%Len()))
      allocate(VResistanceVolumes(VResistanceVolumesList%Len()))
      allocate(WResistanceVolumes(WResistanceVolumesList%Len()))
      allocate(TemperatureFlVolumes(TemperatureFlVolumesList%Len()))
      allocate(MoistureFlVolumes(MoistureFlVolumesList%Len()))

      allocate(ScalarFlVolumes(num_of_scalars))
      
      allocate(VariableScalarFlVolumes(num_of_scalars))
      
      do j=1,num_of_scalars
        allocate(ScalarFlVolumes(j)%volumes(0))
        allocate(VariableScalarFlVolumes(j)%volumes(0))
      end do

      i = 0
      component = 1
      call UResistanceVolumesList%for_each(CopyPoint)
      i = 0
      component = 2
      call VResistanceVolumesList%for_each(CopyPoint)
      i = 0
      component = 3
      call WResistanceVolumesList%for_each(CopyPoint)
      i = 0
      call TemperatureFlVolumesList%for_each(CopyPoint)
      i = 0
      call MoistureFlVolumesList%for_each(CopyPoint)

      do j=1,num_of_scalars
        inorm = 0
        ivar = 0
        call ScalarFlVolumesLists(j)%list%for_each(CopyPoint)
      end do
      
!       if (enable_radiation) call SaveFluxes

      contains

        subroutine CopyPoint(elem)
          class(*) :: elem

          select type (elem)
            type is (TResistanceVolume)          
              i = i + 1
              if (component==1) then
                UResistanceVolumes(i) = elem
              elseif (component==2) then
                VResistanceVolumes(i) = elem
              else
                WResistanceVolumes(i) = elem
              endif
            type is (TemperatureFlVolume)
              i = i + 1
              TemperatureFlVolumes(i) = elem
            type is (MoistureFlVolume)
              i = i + 1
              MoistureFlVolumes(i) = elem
            type is (ScalarFlVolume)
              if (associated(elem%get_flux)) then
                call add_element(VariableScalarFlVolumes(j)%volumes, elem)
              else
                call add_element(ScalarFlVolumes(j)%volumes, elem)
              end if
            class default
              call error_stop("Type error in volume source list.")
          end select
        end subroutine
        
        subroutine add_element(a,e)
          type(ScalarFlVolume),allocatable,intent(inout) :: a(:)
          type(ScalarFlVolume),intent(in) :: e
          type(ScalarFlVolume),allocatable :: tmp(:)

          if (.not.allocated(a)) then
            a = [e]
          else
            call move_alloc(a,tmp)
            allocate(a(size(tmp)+1))
            a(1:size(tmp)) = tmp
            a(size(tmp)+1) = e
          end if
        end subroutine

        subroutine SaveFluxes
          use VTKArray
          real(knd),allocatable :: temperature_flux(:,:,:)
          allocate(temperature_flux(0:Prnx+1,0:Prny+1,0:Prnz+1))
          
          temperature_flux = 0
          
          do i=1,size(TemperatureFlVolumes)
            associate(p => TemperatureFlVolumes(i))
              temperature_flux(p%xi,p%yj,p%zk) = p%flux
            end associate
          end do
          
          call VtkArrayBin(trim(output_dir)//"tempflplants.vtk",temperature_flux)
          
        end subroutine
        
    end subroutine MovePointsToArray



    subroutine InitVolumeSources
 
      call InsideCellsToLists
 
      call MovePointsToArray

    end subroutine InitVolumeSources

    

    subroutine ResistanceForce(U2,V2,W2,U,V,W)
      real(knd),dimension(-2:,-2:,-2:),intent(inout) :: U2,V2,W2
      real(knd),dimension(-2:,-2:,-2:),intent(in)    :: U,V,W

      call apply(U2,U,UResistanceVolumes,totU_u)

      call apply(V2,V,VResistanceVolumes,totU_v)

      call apply(W2,W,WResistanceVolumes,totU_w)

      contains

         pure real(knd) function totU_u(i,j,k)
           integer,intent(in) :: i,j,k

           totU_u = hypot( hypot( U(i,j,k), &
                                  (V(i,j,k)+V(i,j-1,k)+V(i-1,j,k)+V(i-1,j-1,k))/4._knd ) , &
                                  (W(i,j,k)+W(i,j,k-1)+W(i-1,j,k)+W(i-1,j,k-1))/4._knd )
         end function

         pure real(knd) function totU_v(i,j,k)
           integer,intent(in) :: i,j,k

           totU_v = hypot( hypot( (U(i,j,k)+U(i-1,j,k)+U(i,j-1,k)+U(i-1,j-1,k))/4._knd, &
                                  V(i,j,k) ) , &
                                  (W(i,j,k)+W(i,j,k-1)+W(i,j-1,k)+W(i,j-1,k-1))/4._knd )
         end function

         pure real(knd) function totU_w(i,j,k)
           integer,intent(in) :: i,j,k

           totU_w = hypot( hypot( (U(i,j,k)+U(i-1,j,k)+U(i,j,k-1)+U(i-1,j,k-1))/4._knd, &
                                  (V(i,j,k)+V(i,j-1,k)+V(i,j,k-1)+V(i,j-1,k-1))/4._knd ) , &
                                  W(i,j,k) )
         end function

         subroutine apply(X2,X,src,fun)
           real(knd),dimension(-2:,-2:,-2:),intent(inout)  :: X2
           real(knd),dimension(-2:,-2:,-2:),intent(in)     :: X
           type(TResistanceVolume),allocatable,intent(in)  :: src(:)
           procedure(totU_u) :: fun
           integer :: i

           if (allocated(src)) then
             do i=1,size(src)
               associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
                 X2(xi,yj,zk) = X2(xi,yj,zk) - src(i)%flux * fun(xi,yj,zk) * X(xi,yj,zk)              
               end associate
             end do
           end if
         end subroutine
         
    end subroutine ResistanceForce

    subroutine TemperatureVolumeSources(Temperature)
      real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Temperature
      integer :: i
!       call FluxKernel(Temperature,TemperatureFlVolumes)
      associate (X=>Temperature, src=> TemperatureFlVolumes)
       do i=1,size(src)
         associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
           X(xi,yj,zk) = X(xi,yj,zk) + src(i)%flux
         end associate
       end do
      end associate      
    end subroutine TemperatureVolumeSources

    subroutine MoistureVolumeSources(Moisture)
      real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Moisture
      integer :: i
!       call FluxKernel(Moisture,MoistureFlVolumes)
      associate (X=>Moisture, src=> MoistureFlVolumes)
       do i=1,size(src)
         associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
           X(xi,yj,zk) = X(xi,yj,zk) + src(i)%flux
         end associate
       end do
      end associate
    end subroutine MoistureVolumeSources

    subroutine ScalarVolumeSources(Scalar, sc, first_stage)
      real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Scalar
      integer, intent(in) :: sc
      logical, intent(in) :: first_stage

      if (allocated(ScalarFlVolumes)) then
        call FluxKernel(Scalar(:,:,:), ScalarFlVolumes(sc)%volumes)
      end if

      if (allocated(VariableScalarFlVolumes)) then
        if (first_stage) then
          call VariableFluxKernel(Scalar(:,:,:), VariableScalarFlVolumes(sc)%volumes, sc)
        else
          call FluxKernel(Scalar(:,:,:), VariableScalarFlVolumes(sc)%volumes)
        end if
      end if

    end subroutine ScalarVolumeSources

    subroutine FluxKernel(X, src)
      real(knd),intent(inout) :: X(-1:,-1:,-1:)
      type(ScalarFlVolume),intent(in) :: src(:)
      integer :: i
      !Assume src is allocated. It must hold if we called MovePointsToArray properly.
      do i=1,size(src)
        associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
          X(xi,yj,zk) = X(xi,yj,zk) + src(i)%flux
        end associate
      end do

    end subroutine FluxKernel
    
    
    subroutine VariableFluxKernel(X, src, sc)
      real(knd),intent(inout) :: X(-1:,-1:,-1:)
      type(ScalarFlVolume),intent(inout) :: src(:)
      integer, intent(in) :: sc
      integer :: i
      !Assume src is allocated. It must hold if we called MovePointsToArray properly.
      !Assume get_flux is associated. It must hold if we called MovePointsToArray properly.
      do i=1,size(src)
        associate (xp => xPr(src(i)%xi), &
                   yp => yPr(src(i)%yj), &
                   zp => zPr(src(i)%zk))
          src(i)%flux = src(i)%get_flux(xp, yp, zp, sc)
        end associate
      end do
      
      do i=1,size(src)
        associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
          X(xi,yj,zk) = X(xi,yj,zk) + src(i)%flux
        end associate
      end do

    end subroutine VariableFluxKernel
    
    
    
    
    subroutine AddScalarFlVolumesContainer(l,r)
      type(ScalarFlVolumesContainer),intent(inout) :: l(:)
      type(ScalarFlVolumesContainer),intent(in)    :: r
      type(ScalarFlVolume),allocatable :: tmp(:)

      associate (sn => r%scalar_number)
        if (size(r%volumes)>0) then
          !NOTE: the shorter version problematic in gfortran4.8 and ifort 14
          !l(sn)%volumes = [l(sn)%volumes, r%volumes]
          allocate(tmp( size(l(sn)%volumes) + size(r%volumes) ))
          
          tmp(1:size(l(sn)%volumes)) = l(sn)%volumes
          tmp(size(l(sn)%volumes)+1:) = r%volumes
          
          call move_alloc(tmp,l(sn)%volumes)
        end if
      end associate
    end subroutine

end module VolumeSources



