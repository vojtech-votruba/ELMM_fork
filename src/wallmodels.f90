module Sort
  use iso_c_binding
#define pure 
  implicit none

  interface
    subroutine qsort(array,elem_count,elem_size,compare) bind(C,name="qsort")
      import
      type(c_ptr),value       :: array
      integer(c_size_t),value :: elem_count
      integer(c_size_t),value :: elem_size
      type(c_funptr),value    :: compare !int(*compare)(const void *, const void *)
    end subroutine qsort !standard C library qsort
  end interface
end module Sort

module WMPoint_types
  use Kinds

  type point
    integer   :: xi
    integer   :: yj
    integer   :: zk
  end type

  type WMpoint   !points in which we apply wall model

    integer   :: xi
    integer   :: yj
    integer   :: zk

    real(knd) :: distx
    real(knd) :: disty
    real(knd) :: distz

    real(knd) :: area_factor ![m^-1] area of the solid wall divided by the volume of the cell(xi,yj,zk)

    real(knd) :: z0 = 0
    real(knd) :: z0H = 0
    real(knd) :: ustar = 1
    logical :: prescribed_ustar = .false.

    real(knd) :: temperature = 0
    real(knd) :: temperature_flux = 0
    logical :: prescribed_temperature = .false.

    real(knd) :: moisture = 0
    real(knd) :: moisture_flux = 0
    logical :: prescribed_moisture = .false.

    real(knd) :: div = 0 !divergence when zeroing out trans-boundary velocities

    real(knd) :: wallu = 0
    real(knd) :: wallv = 0
    real(knd) :: wallw = 0

    real(knd) :: albedo = 0.1 !for shortwave radiation - asphalt
    real(knd) :: emissivity = 0.9 !for longwave radiation - asphalt
    
    real(knd) :: evaporative_fraction = 0

    real(knd),allocatable:: depscalar(:)

  end type WMpoint


  type WMpointUVW   !points in which we apply wall model
  
    real(knd) :: fluxp, fluxm

    integer   :: xi
    integer   :: yj
    integer   :: zk

    real(knd) :: distx
    real(knd) :: disty
    real(knd) :: distz

    real(knd) :: z0 = 0
    real(knd) :: z0H = 0
    real(knd) :: ustar = 1
    logical :: prescribed_ustar = .false.

    real(knd) :: wallu = 0
    real(knd) :: wallv = 0
    real(knd) :: wallw = 0

    real(knd) :: temperature = 0
    real(knd) :: temperature_flux = 0
    logical :: prescribed_temperature = .false.
        
    real(knd) :: moisture = 0
    real(knd) :: moisture_flux = 0
    logical :: prescribed_moisture = .false.    
  end type WMpointUVW

end module

module WMPointLists
  use WMPoint_types
#define TYPEPARAM type(WMPoint)
#include "list-inc-def.f90"
contains
#include "list-inc-proc.f90"
#undef TYPEPARAM
end module

module WMPointUVWLists
  use WMPoint_types
#define TYPEPARAM type(WMPointUVW)
#include "list-inc-def.f90"
contains
#include "list-inc-proc.f90"
#undef TYPEPARAM
end module


module Wallmodels
  use Parameters
  use WMPoint_types
  use WMPointLists, only: WMPointList => list
  use WMPointUVWLists, only: WMPointUVWList => list

  implicit none

  private
  public WMPoint, WMpointUVW, &
         AddWMPoint, AddWMPointUVW, &
         MoveWMPointsToArray, GetOutsideBoundariesWM, InitWMMasks, &
         ComputeViscsWM, ComputeUVWFluxesWM, DivergenceWM, &
         GroundDeposition, &
         GroundUstar, GroundUstarUVW, &
         GroundTFlux, GroundTFluxUVW, &
         GroundMFlux, GroundMFluxUVW, &
         wallmodeltype

  integer, parameter, public :: MINUSX = 1, PLUSX = 2, MINUSY = 3, PLUSY = 4, MINUSZ = 5, PLUSZ = 6

  type(WMPointUVWList) :: UxmWMPsL, UxpWMPsL, UymWMPsL, UypWMPsL, UzmWMPsL, UzpWMPsL, &
                          VxmWMPsL, VxpWMPsL, VymWMPsL, VypWMPsL, VzmWMPsL, VzpWMPsL, &
                          WxmWMPsL, WxpWMPsL, WymWMPsL, WypWMPsL, WzmWMPsL, WzpWMPsL

  type(WMPointList) :: WMPointsList

  type(WMPointUVW), dimension(:), allocatable, public, target :: &
      UxmWMpoints, UxpWMpoints, UymWMpoints, UypWMpoints, UzmWMpoints, UzpWMPoints, &
      VxmWMpoints, VxpWMpoints, VymWMpoints, VypWMpoints, VzmWMpoints, VzpWMpoints, &
      WxmWMpoints, WxpWMpoints, WymWMpoints, WypWMpoints, WzmWMpoints, WzpWMpoints

  type(WMPoint), dimension(:), allocatable, public :: WMPoints

  type(point), dimension(:), allocatable, public :: Scflx_points, Scfly_points, Scflz_points

  logical(slog), dimension(:,:,:), allocatable, public :: Scflx_mask, Scfly_mask, Scflz_mask

  logical(slog), dimension(:,:,:), allocatable, public :: Uflx_mask, Ufly_mask, Uflz_mask
  logical(slog), dimension(:,:,:), allocatable, public :: Vflx_mask, Vfly_mask, Vflz_mask
  logical(slog), dimension(:,:,:), allocatable, public :: Wflx_mask, Wfly_mask, Wflz_mask

  integer :: wallmodeltype
  
  type point_container
    type(WMPointUVW), dimension(:), pointer :: points
  end type
  
  type(point_container), public :: WMPointsUVW(6,3)

contains



  impure elemental subroutine ComputeAreaFactor(p)
    type(WMpoint),intent(inout):: p
    real(knd) :: out_norm(3), d, l, o, S, V
    
    !FIXME: Not accurate!!!
    p%area_factor = 0
    if (all([-p%distx, -p%disty, -p%distz]==0)) return !HACK!!!
    
    V = dxmin*dymin*dzmin
    
    out_norm = [-p%distx, -p%disty, -p%distz] / &
                     norm2([p%distx,p%disty,p%distz])
    out_norm = abs(out_norm)

    !FIXME
    if (abs(out_norm(1))>abs(out_norm(2)).and.abs(out_norm(1))>abs(out_norm(3))) then
      o = out_norm(1)
      d = sqrt(dymin*dzmin)
    else if (abs(out_norm(2))>abs(out_norm(1)).and.abs(out_norm(2))>abs(out_norm(3))) then
      o = out_norm(2)
      d = sqrt(dxmin*dzmin)
    else
      o = out_norm(3)
      d = sqrt(dxmin*dymin)
    end if
    
    l = d / o
    S = d*l
    
    p%area_factor = S / V
  end subroutine  



  subroutine AddWMPoint(WMP)
    type(WMPoint), intent(in) :: WMP

    call WMPointsList%add(WMP)

  end subroutine


  subroutine AddWMPointUVW(WMP, component, direction)
    type(WMPointUVW), intent(in) :: WMP
    integer, intent(in) :: component, direction

    select case (component)
      case (1)
        call add(WMP, direction, UxmWMPsL, UxpWMPsL, UymWMPsL, UypWMPsL, UzmWMPsL, UzpWMPsL)
      case (2)
        call add(WMP, direction, VxmWMPsL, VxpWMPsL, VymWMPsL, VypWMPsL, VzmWMPsL, VzpWMPsL)
      case (3)
        call add(WMP, direction, WxmWMPsL, WxpWMPsL, WymWMPsL, WypWMPsL, WzmWMPsL, WzpWMPsL)
    end select

  contains
    subroutine add(p, dir, xm, xp, ym, yp, zm, zp)
      type(WMPointUVW), intent(in) :: p
      integer, intent(in) :: dir
      type(WMPointUVWList), intent(inout) :: xm, xp, ym, yp, zm, zp
      select case (dir)
        case (MINUSX)
          call xm%add(p)
        case (PLUSX)
          call xp%add(p)
        case (MINUSY)
          call ym%add(p)
        case (PLUSY)
          call yp%add(p)
        case (MINUSZ)
          call zm%add(p)
        case (PLUSZ)
          call zp%add(p)
      end select
    end subroutine
  end subroutine




  subroutine MoveWMPointsToArray

    call MoveWMPointsToArrayPr
    call MoveWMPointsToArrayUVW
  end subroutine




  subroutine MoveWMPointsToArrayPr
    type(WMPoint),pointer :: p
    integer :: i

    allocate(WMPoints(WMPointsList%len()))

    if (size(WMPoints)>0) then

      call WMPointsList%iter_restart

      do i = 1, size(WMPoints)
        p => WMPointsList%iter_next()
        if (associated(p)) then
          WMPoints(i) = p
        else
          call error_stop("Assert error, pointer not associated. File "// &
            __FILE__ &
            //" line ",__LINE__)
        end if
      end do

      call RemoveDuplicateWMPoints(WMPoints)

      call ComputeAreaFactor(WMPoints)

    end if

    call WMPointsList%finalize

  end subroutine MoveWMPointsToArrayPr



  subroutine MoveWMPointsToArrayUVW
    type(WMPointUVW),pointer :: p
    integer :: i

    call helper(UxmWMPsL, UxmWMpoints)
    call helper(UxpWMPsL, UxpWMpoints)
    call helper(UymWMPsL, UymWMpoints)
    call helper(UypWMPsL, UypWMpoints)
    call helper(UzmWMPsL, UzmWMpoints)
    call helper(UzpWMPsL, UzpWMpoints)

    call helper(VxmWMPsL, VxmWMpoints)
    call helper(VxpWMPsL, VxpWMpoints)
    call helper(VymWMPsL, VymWMpoints)
    call helper(VypWMPsL, VypWMpoints)
    call helper(VzmWMPsL, VzmWMpoints)
    call helper(VzpWMPsL, VzpWMpoints)

    call helper(WxmWMPsL, WxmWMpoints)
    call helper(WxpWMPsL, WxpWMpoints)
    call helper(WymWMPsL, WymWMpoints)
    call helper(WypWMPsL, WypWMpoints)
    call helper(WzmWMPsL, WzmWMpoints)
    call helper(WzpWMPsL, WzpWMpoints)
    
    
    WMPointsUVW(1,1)%points => UxmWMpoints
    WMPointsUVW(2,1)%points => UxpWMpoints
    WMPointsUVW(3,1)%points => UymWMpoints
    WMPointsUVW(4,1)%points => UypWMpoints
    WMPointsUVW(5,1)%points => UzmWMpoints
    WMPointsUVW(6,1)%points => UzpWMpoints

    WMPointsUVW(1,2)%points => VxmWMpoints
    WMPointsUVW(2,2)%points => VxpWMpoints
    WMPointsUVW(3,2)%points => VymWMpoints
    WMPointsUVW(4,2)%points => VypWMpoints
    WMPointsUVW(5,2)%points => VzmWMpoints
    WMPointsUVW(6,2)%points => VzpWMpoints

    WMPointsUVW(1,3)%points => WxmWMpoints
    WMPointsUVW(2,3)%points => WxpWMpoints
    WMPointsUVW(3,3)%points => WymWMpoints
    WMPointsUVW(4,3)%points => WypWMpoints
    WMPointsUVW(5,3)%points => WzmWMpoints
    WMPointsUVW(6,3)%points => WzpWMpoints
    
  contains
    subroutine helper(l, arr)
      type(WMPointUVWList), intent(inout) :: l
      type(WMPointUVW), allocatable, intent(out) :: arr(:)

      allocate(arr(l%len()))

      if (size(arr)>0) then

        call l%iter_restart

        do i = 1, size(arr)
          p => l%iter_next()
          if (associated(p)) then
            arr(i) = p
          else
            write(*,*) "Assert error, pointer not associated. File ", &
              __FILE__ &
              ," line ",__LINE__
            call error_stop
          end if
        end do

      end if

      call l%finalize
    end subroutine
  end subroutine MoveWMPointsToArrayUVW



  subroutine RemoveDuplicateWMPoints(WMPoints)
    use iso_c_binding
    use Sort
    type(WMPoint),allocatable,dimension(:),target,intent(inout)  :: WMPoints
    type(WMPoint),allocatable,dimension(:) :: TMP
    integer :: i,n

    !Choose the one closer to a wall. If of the same distance, choose the later one.
    !For wider compatibility we do not use MOLD= or SOURCE= in allocate.

    if (size(WMPoints)>1) then

      allocate(TMP(size(WMPoints)))

      call qsort(c_loc(WMPoints(1)), &
                 size(WMPoints,kind = c_size_t), &
                 storage_size(WMPoints,c_size_t)/storage_size(c_char_'a',c_size_t), &
  !                c_sizeof(WMPoints), &  F08+TS29113 seem to require C interoperable variable as argument.
                 c_funloc(CompareWMPoints))

      TMP(1) = WMPoints(1)
      n = 1
      do i = 2,size(WMPoints)
          if (WMPoints(i-1)%xi/=WMPoints(i)%xi .or. &
              WMPoints(i-1)%yj/=WMPoints(i)%yj .or. &
              WMPoints(i-1)%zk/=WMPoints(i)%zk) then !if the point is not duplicate of previous one

                n = n+1
                TMP(n) = WMPoints(i)

          end if
      end do

      WMPoints = TMP(1:n)

   end if

  end subroutine RemoveDuplicateWMPoints






  function CompareWMPoints(Aptr,Bptr) bind(C,name="CompareWMPoints") result(res)
    use iso_c_binding
    integer(c_int)         :: res
    type(c_ptr),value :: Aptr,Bptr
    type(WMPoint),pointer  :: A,B

    call c_f_pointer(Aptr,A)
    call c_f_pointer(Bptr,B)

    if ((A%xi+(A%yj-1)*Prnx+(A%zk-1)*Prnx*Prny) < (B%xi+(B%yj-1)*Prnx+(B%zk-1)*Prnx*Prny)) then
      res = -1_c_int
    else if ((A%xi+(A%yj-1)*Prnx+(A%zk-1)*Prnx*Prny) > (B%xi+(B%yj-1)*Prnx+(B%zk-1)*Prnx*Prny)) then
      res =  1_c_int
    else if (A%distx**2+A%disty**2+A%distz**2 < B%distx**2+B%disty**2+B%distz**2) then
      res = -1_c_int
    else if (A%distx**2+A%disty**2+A%distz**2 > B%distx**2+B%disty**2+B%distz**2) then
      res =  1_c_int
    else
      res =  0_c_int
    end if

  end function CompareWMPoints









 
  subroutine GetOutsideBoundariesWM(nscalars)
    integer, intent(in) :: nscalars

    if (wallmodeltype/=0) then

      call GetOutsideBoundariesWM_UVW

      call GetOutsideBoundariesWM_Pr(nscalars)

    end if

  end subroutine






  subroutine GetOutsideBoundariesWM_Pr(nscalars)
    integer, intent(in) :: nscalars
    integer       :: i,j,k
    type(WMPoint) :: WMP

    allocate(WMP%depscalar(nscalars))
    WMP%depscalar = 0
    WMP%evaporative_fraction = 0.1

    !type==0 below, therefore we need 
    ! type of free bounderies not to be -1
    if (Btype(We)==BC_NOSLIP) then
      do k = 1,Prnz
       do j = 1,Prny
         if (Prtype(1,j,k)==0) then
           WMP%xi = 1
           WMP%yj = j
           WMP%zk = k
           WMP%distx = xU(0) - xPr(1)
           WMP%disty = 0
           WMP%distz = 0
           WMP%ustar = 1

           WMP%z0 = z0W
           WMP%z0H = WMP%z0
           call AddWMPoint(WMP)
         end if
       end do
      end do
    end if

    if (Btype(Ea)==BC_NOSLIP) then
      do k = 1,Prnz
       do j = 1,Prny
         if (Prtype(Prnx,j,k)==0) then
           WMP%xi = Prnx
           WMP%yj = j
           WMP%zk = k
           WMP%distx = xU(Unx+1) - xPr(Prnx)
           WMP%disty = 0
           WMP%distz = 0
           WMP%ustar = 1

           WMP%z0 = z0E
           WMP%z0H = WMP%z0
           call AddWMPoint(WMP)
         end if
       end do
      end do
    end if

    if (Btype(So)==BC_NOSLIP.or.(Btype(So)==BC_DIRICHLET.and.sideU(2,So)==0)) then
      do k = 1,Prnz
       do i = 1,Prnx
         if (Prtype(i,1,k)==0) then
           WMP%xi = i
           WMP%yj = 1
           WMP%zk = k
           WMP%distx = 0
           WMP%disty = yV(0) - yPr(1)
           WMP%distz = 0
           WMP%ustar = 1

           if (Btype(So)==BC_DIRICHLET) then
             WMP%wallu = sideU(1,So)
             WMP%wallv = 0
             WMP%wallw = sideU(3,So)
           end if

           WMP%z0 = z0S
           WMP%z0H = WMP%z0
           call AddWMPoint(WMP)
         end if
       end do
      end do
    end if

    if (Btype(No)==BC_NOSLIP.or.(Btype(No)==BC_DIRICHLET.and.sideU(2,No)==0)) then
      do k = 1,Prnz
       do i = 1,Prnx
         if (Prtype(i,Prny,k)==0) then
           WMP%xi = i
           WMP%yj = Prny
           WMP%zk = k
           WMP%distx = 0
           WMP%disty = yV(Vny+1) - yPr(Prny)
           WMP%distz = 0
           WMP%ustar = 1

           if (Btype(No)==BC_DIRICHLET) then
             WMP%wallu = sideU(1,No)
             WMP%wallv = 0
             WMP%wallw = sideU(3,No)
           end if

           WMP%z0 = z0N
           WMP%z0H = WMP%z0
           call AddWMPoint(WMP)
         end if
       end do
      end do
    end if

    if (Btype(Bo)==BC_NOSLIP.or.(Btype(Bo)==BC_DIRICHLET.and.sideU(3,Bo)==0)) then
      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,1)==0) then

           WMP%xi = i
           WMP%yj = j
           WMP%zk = 1
           WMP%distx = 0
           WMP%disty = 0
           WMP%distz = zW(0) - zPr(1)
           WMP%ustar = 1

           if (Btype(Bo)==BC_DIRICHLET) then
             WMP%wallu = sideU(1,Bo)
             WMP%wallv = sideU(2,Bo)
             WMP%wallw = 0
           end if

           WMP%z0 = z0B
           WMP%z0H = WMP%z0

           if (TempBtype(Bo)==BC_CONSTFLUX) then
             WMP%temperature_flux = sideTemp(Bo)
           else
             WMP%temperature_flux = 0
           end if

           if (TempBtype(Bo)==BC_DIRICHLET) then
             WMP%temperature = sideTemp(Bo)
             WMP%prescribed_temperature = .true.
           end if

           if (MoistBtype(Bo)==BC_CONSTFLUX) then
             WMP%moisture_flux = sideMoist(Bo)
          else
             WMP%moisture_flux = 0
           end if

           if (MoistBtype(Bo)==BC_DIRICHLET) then
             WMP%moisture = sideMoist(Bo)
             WMP%prescribed_moisture = .true.
           end if

           call AddWMPoint(WMP)

         end if
       end do
      end do

      if (TempBtype(Bo)==BC_CONSTFLUX.or.TempBtype(Bo)==BC_DIRICHLET) then
        TempBtype(Bo) = 0
      end if

      if (MoistBtype(Bo)==BC_CONSTFLUX) then
        MoistBtype(Bo) = 0
      end if

    end if

    if (Btype(To)==BC_NOSLIP.or.(Btype(To)==BC_DIRICHLET.and.sideU(3,To)==0)) then

      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,Prnz)==0) then
           WMP%xi = i
           WMP%yj = j
           WMP%zk = Prnz
           WMP%distx = 0
           WMP%disty = 0
           WMP%distz = zW(Wnz+1) - zPr(Prnz)
           WMP%ustar = 1

           if (Btype(To)==BC_DIRICHLET) then
             WMP%wallu = sideU(1,To)
             WMP%wallv = sideU(2,To)
             WMP%wallw = 0
           end if
            
           WMP%z0 = z0T
           WMP%z0H = WMP%z0

           if (TempBtype(To)==BC_CONSTFLUX) then
             WMP%temperature_flux = - sideTemp(To)
           else
             WMP%temperature_flux = 0
           end if

           if (TempBtype(To)==BC_DIRICHLET) then
             WMP%temperature = sideTemp(To)
             WMP%prescribed_temperature = .true.
           end if

           if (MoistBtype(To)==BC_CONSTFLUX) then
             WMP%moisture_flux = - sideMoist(To)
           end if

           if (MoistBtype(To)==BC_DIRICHLET) then
             WMP%moisture = sideMoist(To)
             WMP%prescribed_moisture = .true.
           end if

           call AddWMPoint(WMP)
         end if
       end do
      end do

      if (TempBtype(To)==BC_CONSTFLUX.or.TempBtype(To)==BC_DIRICHLET) then
        TempBtype(To) = 0
      end if

      if (MoistBtype(To)==BC_CONSTFLUX) then
        MoistBtype(To) = 0
      end if

    end if

  end subroutine GetOutsideBoundariesWM_Pr


  subroutine GetOutsideBoundariesWM_UVW

    call helper(1, Unx, Uny, Unz, xU(-2:), yPr, zPr, Utype)

    call helper(2, Vnx, Vny, Vnz, xPr, yV(-2:), zPr, Vtype)

    call helper(3, Wnx, Wny, Wnz, xPr, yPr, zW(-2:), Wtype)

  contains

    subroutine helper(component, nx, ny, nz, x, y, z, Xtype)
      integer, intent(in) :: component, nx, ny, nz
      real(knd), intent(in) :: x(-2:), y(-2:), z(-2:)
      integer, intent(in) :: Xtype(-2:,-2:,-2:)
      integer       :: i,j,k
      type(WMPointUVW) :: p

      !type==0 below, therefore we need 
      ! type of free bounderies not to be -1
      if (Btype(We)==BC_NOSLIP) then
        do k = 1,nz
         do j = 1,ny
           if (Xtype(1,j,k)==0.and.Xtype(0,j,k)<=0) then
             p%xi = 1
             p%yj = j
             p%zk = k
             p%distx = xU(0) - x(1)
             p%disty = 0
             p%distz = 0
             p%ustar = 1

             p%z0 = z0W
             p%z0H = p%z0
             call AddWMPointUVW(p, component, We)
           end if
         end do
        end do
      end if

      if (Btype(Ea)==BC_NOSLIP) then
        do k = 1,nz
         do j = 1,ny
           if (Xtype(nx,j,k)==0.and.Xtype(nx+1,j,k)<=0) then
             p%xi = nx
             p%yj = j
             p%zk = k
             p%distx = xU(Unx+1) - x(nx)
             p%disty = 0
             p%distz = 0
             p%ustar = 1

             p%z0 = z0E
             p%z0H = p%z0
             call AddWMPointUVW(p, component, Ea)
           end if
         end do
        end do
      end if

      if (Btype(So)==BC_NOSLIP.or.(Btype(So)==BC_DIRICHLET.and.sideU(2,So)==0)) then
        do k = 1,nz
         do i = 1,nx
           if (Xtype(i,1,k)==0.and.Xtype(i,0,k)<=0) then
             p%xi = i
             p%yj = 1
             p%zk = k
             p%distx = 0
             p%disty = yV(0) - y(1)
             p%distz = 0
             p%ustar = 1

             if (Btype(So)==BC_DIRICHLET) then
               p%wallu = sideU(1,So)
               p%wallv = 0
               p%wallw = sideU(3,So)
             end if

             p%z0 = z0S
             p%z0H = p%z0
             call AddWMPointUVW(p, component, So)
           end if
         end do
        end do
      end if

      if (Btype(No)==BC_NOSLIP.or.(Btype(No)==BC_DIRICHLET.and.sideU(2,No)==0)) then
        do k = 1,nz
         do i = 1,nx
           if (Xtype(i,ny,k)==0.and.Xtype(i,ny+1,k)<=0) then
             p%xi = i
             p%yj = ny
             p%zk = k
             p%distx = 0
             p%disty = yV(Vny+1) - y(ny)
             p%distz = 0
             p%ustar = 1

             if (Btype(No)==BC_DIRICHLET) then
               p%wallu = sideU(1,No)
               p%wallv = 0
               p%wallw = sideU(3,No)
             end if

             p%z0 = z0N
             p%z0H = p%z0
             call AddWMPointUVW(p, component, No)
           end if
         end do
        end do
      end if

      if (Btype(Bo)==BC_NOSLIP.or.(Btype(Bo)==BC_DIRICHLET.and.sideU(3,Bo)==0)) then
        do j = 1,ny
         do i = 1,nx
           if (Xtype(i,j,1)==0.and.Xtype(i,j,0)<=0) then

             p%xi = i
             p%yj = j
             p%zk = 1
             p%distx = 0
             p%disty = 0
             p%distz = zW(0) - z(1)
             p%ustar = 1

             if (Btype(Bo)==BC_DIRICHLET) then
               p%wallu = sideU(1,Bo)
               p%wallv = sideU(2,Bo)
               p%wallw = 0
             end if

             p%z0 = z0B
             p%z0H = p%z0

             if (TempBtype(Bo)==BC_CONSTFLUX) then
               p%temperature_flux = sideTemp(Bo)
             else
               p%temperature_flux = 0
             end if

             if (TempBtype(Bo)==BC_DIRICHLET) then
               p%temperature = sideTemp(Bo)
               p%prescribed_temperature = .true.
             end if

             if (MoistBtype(Bo)==BC_CONSTFLUX) then
               p%moisture_flux = sideMoist(Bo)
             else
               p%moisture_flux = 0
             end if

             if (MoistBtype(Bo)==BC_DIRICHLET) then
               p%moisture = sideMoist(Bo)
               p%prescribed_moisture = .true.
             end if

             call AddWMPointUVW(p, component, Bo)

           end if
         end do
        end do
      end if

      if (Btype(To)==BC_NOSLIP.or.(Btype(To)==BC_DIRICHLET.and.sideU(3,To)==0)) then

        do j = 1,ny
         do i = 1,nx
           if (Xtype(i,j,nz)==0.and.Xtype(i,j,nz+1)<=0) then
             p%xi = i
             p%yj = j
             p%zk = nz
             p%distx = 0
             p%disty = 0
             p%distz = zW(Wnz+1) - z(nz)
             p%ustar = 1

             if (Btype(To)==BC_DIRICHLET) then
               p%wallu = sideU(1,To)
               p%wallv = sideU(2,To)
               p%wallw = 0
             end if

             p%z0 = z0T
             p%z0H = p%z0
             
             if (TempBtype(To)==BC_CONSTFLUX) then
               p%temperature_flux = - sideTemp(To)
             else
               p%temperature_flux = 0
             end if

             if (TempBtype(To)==BC_DIRICHLET) then
               p%temperature = sideTemp(To)
               p%prescribed_temperature = .true.
             end if
             
             if (MoistBtype(To)==BC_CONSTFLUX) then
               p%moisture_flux = - sideMoist(To)
             else
               p%moisture_flux = 0
             end if

             if (MoistBtype(To)==BC_DIRICHLET) then
               p%moisture = sideMoist(To)
               p%prescribed_moisture = .true.
             end if
             
             call AddWMPointUVW(p, component, To)
           end if
         end do
        end do
      end if

    end subroutine helper

  end subroutine GetOutsideBoundariesWM_UVW




  
  
  subroutine InitWMMasks
    integer :: i, j, k
    integer :: nx, ny, nz

    allocate(Scflx_mask(0:Prnx,Prny,Prnz))
    allocate(Scfly_mask(Prnx,0:Prny,Prnz))
    allocate(Scflz_mask(Prnx,Prny,0:Prnz))

    !$omp parallel
    !$omp workshare
    Scflx_mask = .true.
    Scfly_mask = .true.
    Scflz_mask = .true.
    !$omp end workshare
    !$omp end parallel

    if (wallmodeltype/=0) then

      !$omp parallel
      !$omp do collapse(3) schedule(dynamic, 100)
      do k = 1, Prnz
        do j = 1, Prny
          do i = 0, Prnx
            if (Prtype(i,j,k)>0.or.Prtype(i+1,j,k)>0) Scflx_mask(i,j,k) = .false.
          end do
        end do
      end do
      !$omp do collapse(3) schedule(dynamic, 100)
      do k = 1, Prnz
        do j = 0, Prny
          do i = 1, Prnx
            if (Prtype(i,j,k)>0.or.Prtype(i,j+1,k)>0) Scfly_mask(i,j,k) = .false.
          end do
        end do
      end do
      !$omp do collapse(3) schedule(dynamic, 100)
      do k = 0, Prnz
        do j = 1, Prny
          do i = 1, Prnx
            if (Prtype(i,j,k)>0.or.Prtype(i,j,k+1)>0) Scflz_mask(i,j,k) = .false.
          end do
        end do
      end do


      !$omp single
      nx = 0
      ny = 0
      nz = 0
      !$omp end single

      !$omp do reduction(+:nx) collapse(3) schedule(dynamic, 100)
      do k = 1, Prnz
        do j = 1, Prny
          do i = 0, Prnx
            if (Prtype(i,j,k)>0.neqv.Prtype(i+1,j,k)>0) nx = nx + 1
         end do
        end do
      end do
      !$omp do reduction(+:ny) collapse(3) schedule(dynamic, 100)
      do k = 1, Prnz
        do j = 0, Prny
          do i = 1, Prnx
            if (Prtype(i,j,k)>0.neqv.Prtype(i,j+1,k)>0) ny = ny + 1
          end do
        end do
      end do
      !$omp do reduction(+:nz) collapse(3) schedule(dynamic, 100)
      do k = 0, Prnz
        do j = 1, Prny
          do i = 1, Prnx
            if (Prtype(i,j,k)>0.neqv.Prtype(i,j,k+1)>0) nz = nz + 1
          end do
        end do
      end do
      !$omp end parallel

      allocate(Scflx_points(nx))
      allocate(Scfly_points(ny))
      allocate(Scflz_points(nz))
      nx = 0
      ny = 0
      nz = 0

      
      do k = 1, Prnz
        do j = 1, Prny
          do i = 0, Prnx
            if (Prtype(i,j,k)>0.neqv.Prtype(i+1,j,k)>0) then
              nx = nx + 1
              Scflx_points(nx) = point(i,j,k)
            end if
          end do
        end do
      end do
      do k = 1, Prnz
        do j = 0, Prny
          do i = 1, Prnx
            if (Prtype(i,j,k)>0.neqv.Prtype(i,j+1,k)>0) then
              ny = ny + 1
              Scfly_points(ny) = point(i,j,k)
            end if
          end do
        end do
      end do
      do k = 0, Prnz
        do j = 1, Prny
          do i = 1, Prnx
            if (Prtype(i,j,k)>0.neqv.Prtype(i,j,k+1)>0) then
              nz = nz + 1
              Scflz_points(nz) = point(i,j,k)
            end if
          end do
        end do
      end do
    else
    
      allocate(Scflx_points(0))
      allocate(Scfly_points(0))
      allocate(Scflz_points(0))

    end if

    call set_masks(Uflx_mask, Ufly_mask, Uflz_mask, &
                   UxmWMpoints, UxpWMpoints, UymWMpoints, UypWMpoints, UzmWMpoints, UzpWMpoints, Unx, Uny, Unz, Utype)

    call set_masks(Vflx_mask, Vfly_mask, Vflz_mask, &
                   VxmWMpoints, VxpWMpoints, VymWMpoints, VypWMpoints, VzmWMpoints, VzpWMpoints, Vnx, Vny, Vnz, Vtype)

    call set_masks(Wflx_mask, Wfly_mask, Wflz_mask, &
                   WxmWMpoints, WxpWMpoints, WymWMpoints, WypWMpoints, WzmWMpoints, WzpWMpoints, Wnx, Wny, Wnz, Wtype)

  contains
  
    subroutine set_masks(flx, fly, flz, xm, xp, ym, yp, zm, zp, nx, ny, nz, Xtype)
      logical(slog), dimension(:,:,:), allocatable, intent(out) :: flx, fly, flz
      type(WMpointUVW),  dimension(:), intent(in) :: xm, xp, ym, yp, zm, zp
      integer, intent(in) :: nx, ny, nz
      integer, intent(in) :: Xtype(-2:,-2:,-2:)
      integer :: i, j, k

      allocate(flx(nx+1,ny,nz))
      allocate(fly(nx,ny+1,nz))
      allocate(flz(nx,ny,nz+1))
  
      flx = .true.
      fly = .true.
      flz = .true.

      if (wallmodeltype/=0) then

        !$omp parallel

        !$omp do private(i)
        do i = 1, size(xm)
          flx(xm(i)%xi,xm(i)%yj,xm(i)%zk) = .false.
        end do
        !$omp end do
        !$omp do private(i)
        do i = 1, size(xp)
          flx(xp(i)%xi+1,xp(i)%yj,xp(i)%zk) = .false.
        end do
        !$omp end do nowait
        !$omp do private(i)
        do i = 1, size(ym)
          fly(ym(i)%xi,ym(i)%yj,ym(i)%zk) = .false.
        end do
        !$omp end do
        !$omp do private(i)
        do i = 1, size(yp)
          fly(yp(i)%xi,yp(i)%yj+1,yp(i)%zk) = .false.
        end do
        !$omp end do nowait
        !$omp do private(i)
        do i = 1, size(zm)
          flz(zm(i)%xi,zm(i)%yj,zm(i)%zk) = .false.
        end do
        !$omp end do
        !$omp do private(i)
        do i = 1, size(zp)
          flz(zp(i)%xi,zp(i)%yj,zp(i)%zk+1) = .false.
        end do
        !$omp end do

        !$omp do private(i,j,k) collapse(3) schedule(dynamic, 100)
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if (Xtype(i,j,k)>0.and.Xtype(i-1,j,k)>0) flx(i,j,k) = .false.
              if (Xtype(i,j,k)>0.and.Xtype(i,j-1,k)>0) fly(i,j,k) = .false.
              if (Xtype(i,j,k)>0.and.Xtype(i,j,k-1)>0) flz(i,j,k) = .false.
            end do
          end do
        end do
        !$omp end do

        !$omp end parallel

      end if

    end subroutine

  end subroutine InitWMMasks



  subroutine DivergenceWM(U, V, W)
    real(knd), contiguous, intent(in)    :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer :: i, xi, yj, zk
    real(knd) :: div

    !$omp parallel do private(i, xi, yj, zk, div)
    do i=1,size(WMPoints)
      xi = WMPoints(i)%xi
      yj = WMPoints(i)%yj
      zk = WMPoints(i)%zk

      div = 0
      if (Prtype(xi+1,yj,zk)<=0) div = div + U(xi,yj,zk) / dxmin
      if (Prtype(xi-1,yj,zk)<=0) div = div - U(xi-1,yj,zk) / dxmin
      if (Prtype(xi,yj+1,zk)<=0) div = div + V(xi,yj,zk) / dymin
      if (Prtype(xi,yj-1,zk)<=0) div = div - V(xi,yj-1,zk) / dymin
      if (Prtype(xi,yj,zk+1)<=0) div = div + W(xi,yj,zk) / dzPr(zk)
      if (Prtype(xi,yj,zk-1)<=0) div = div - W(xi,yj,zk-1) / dzPr(zk)
      WMPoints(i)%div = div
    end do
    !$omp end parallel do

  end subroutine





















  pure subroutine WMFlatUstar(ustar,vel,dist)
    real(knd),intent(inout) :: ustar
    real(knd),intent(in) :: vel,dist
    real(knd),parameter :: eps = 1e-4_knd
    real(knd),parameter :: zpl_lam = 5._knd, zpl_turb = 30._knd
    real(knd) :: ustar2, ustar_lam, z_pl
    integer :: i

    if (molecular_viscosity <= 0) then
      ustar = 0
      return
    end if

    ustar_lam = sqrt(molecular_viscosity * vel / dist)

    if ((dist * ustar_lam / molecular_viscosity) < zpl_lam) then

      ustar = ustar_lam

    else   !turbulent region
      i = 0

      do
        i = i+1
        ustar2 = ustar

        z_pl = dist * ustar2 / molecular_viscosity

        if (z_pl <= zpl_lam) then
          !viscous sublayer
          ustar = sqrt(molecular_viscosity * vel / dist)
        else if (z_pl < zpl_turb) then
          !buffer layer
          ustar = vel / (5*log(z_pl)-3.05)
        else
          !logarithmic layer
          ustar = vel / (log(z_pl)/0.4_knd + 5.5_knd)
        end if

        if  (abs(ustar-ustar2)/abs(ustar)<eps) exit

        if (i>=50) then
                    ustar = 0
                    exit
        end if

      end do

    end if

  end subroutine WMFlatUstar


  subroutine WMFlatStress(ustar,distvect,uvect,walluvect,tan_vel,distance,tan_vect)
    real(knd), intent(inout) :: ustar
    real(knd), intent(in)    :: distvect(3), uvect(3), walluvect(3)
    real(knd), optional, intent(out)   :: tan_vel, distance, tan_vect(3)
    real(knd) :: vel, dist, tan_vect_loc(3)

    call vel_and_dist(tan_vect_loc, dist, uvect, walluvect, distvect)

    vel = norm2(tan_vect_loc)

    if (vel/=0) then

      call WMFlatUstar(ustar,vel,dist)

    end if

    if (ustar<0) ustar = 0

    if (present(tan_vel)) tan_vel = vel
    if (present(distance)) distance = dist
    if (present(tan_vect)) tan_vect = tan_vect_loc

  end subroutine WMFlatStress


  subroutine WMFlatVisc(visc,ustar,distvect,uvect,walluvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: distvect(3), uvect(3), walluvect(3)
    real(knd) :: vel, dist

    if (molecular_viscosity<=0) then
      ustar = 0
      visc = 0
      return
    end if

    call WMFlatStress(ustar,distvect,uvect,walluvect,vel,dist)

    if (vel>0) then

     if (dist*ustar/molecular_viscosity>1) then

       visc = ustar**2 * dist/vel

     else

       visc = molecular_viscosity

     end if

    else if (molecular_viscosity > 0) then

      visc = molecular_viscosity

    else

      visc = 0

    end if

  end subroutine WMFlatVisc




  pure subroutine WMRoughUstar(ustar,vel,dist,z0)
    real(knd),intent(inout) :: ustar
    real(knd),intent(in) :: vel,dist,z0
    real(knd),parameter  :: eps = 1e-4_knd
    real(knd),parameter  :: yplcrit = 11.225_knd

    !if in the laminar region for flat boundary layer, treat as flat and laminar
    if (molecular_viscosity>0) then
      ustar = sqrt(molecular_viscosity * vel / dist)
      if (dist * ustar / molecular_viscosity < yplcrit*1.1_knd) return
    end if
        

    !under z0 the whole concept of rougness parameter breaks down
    if (dist<=z0*2_knd) then
      if (molecular_viscosity > 0) then
        call WMFlatUstar(ustar,vel,dist)
      else
        ustar = 0
      end if
    else 
      ustar = vel * 0.41_knd / log(dist/z0)
    end if

  end subroutine WMRoughUstar


  pure subroutine WMRoughStress(ustar,z0,distvect,uvect,walluvect,tan_vel,distance,tan_vect)
    real(knd), intent(inout) :: ustar
    real(knd), intent(in)    :: z0
    real(knd), intent(in)    :: distvect(3), uvect(3), walluvect(3)
    real(knd), optional, intent(out)   :: tan_vel, distance, tan_vect(3)
    real(knd) :: vel, dist, tan_vect_loc(3)

    call vel_and_dist(tan_vect_loc, dist, uvect, walluvect, distvect)

    vel = norm2(tan_vect_loc)

    if (vel/=0) then

      call WMRoughUstar(ustar,vel,dist,z0)

    end if

    if (ustar<0) ustar = 0

    if (present(tan_vel)) tan_vel = vel
    if (present(distance)) distance = dist
    if (present(tan_vect)) tan_vect = tan_vect_loc

  end subroutine WMRoughStress


  pure subroutine WMRoughVisc(visc,ustar,z0,distvect,uvect,walluvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: z0
    real(knd),intent(in)    :: distvect(3), uvect(3), walluvect(3)
    real(knd) :: vel, dist

    call WMRoughStress(ustar,z0,distvect,uvect,walluvect,vel,dist)

    if (vel>0 .and. ustar**2 * dist/vel > molecular_viscosity) then
      visc = ustar**2 * dist/vel
    else if (molecular_viscosity > 0) then
      visc = molecular_viscosity
    else
      visc = 0
    end if

  end subroutine WMRoughVisc


  pure subroutine vel_and_dist(tan_vect, dist, uvect, walluvect, distvect)
    real(knd), intent(out) :: tan_vect(3),  dist
    real(knd), intent(in)  :: uvect(3), walluvect(3), distvect(3)

    dist = norm2(distvect)

    tan_vect = uvect - walluvect

    tan_vect = tan_vect - dot_product(tan_vect,distvect) * distvect / dist**2  !tangential part
  end subroutine





  pure subroutine WMFlatPrGradUstar(ustar,vel,prgrad,dist)
    real(knd),intent(out) :: ustar
    real(knd),intent(in) :: vel,prgrad,dist
    real(knd),parameter :: eps = 1e-3_knd
    real(knd),parameter :: yplcrit = 11.225_knd
    real(knd),parameter :: ustar_div = 1000
    real(knd) :: ustar1,ustar2,ustar_lam
    integer :: i
    integer,parameter :: maxiter = 30

    !u/u_* = z * u_* / nu  +  dp/dx * z**2 / (2 * u_*)
    ustar_lam = sqrt(abs((molecular_viscosity/dist) * (vel - dist**2 * prgrad/2) ))

    if ((dist*ustar_lam / molecular_viscosity)<yplcrit) then

      ustar = ustar_lam

    else   !turbulent region


      i = 0

      ustar1 = ustar_lam

      do
        i = i+1

        ustar2 = newguess(ustar1)
        if (ustar2 < 0) then
          if (ustar1>max(100*vel,10*Uinlet).or.i>=30) then
            ustar = ustar_lam
            exit
          end if
          ustar1 = ustar1 * 10
        else if (ustar2<tiny(1._knd)) then
          ustar = ustar_lam
          exit
        else if (abs(ustar1-ustar2)/abs(ustar1)<eps) then
          ustar = ustar2
          exit
        else if (i<20) then
          ustar1 = ustar2
        else
          if (abs(ustar1-ustar2)/abs(ustar1)<0.1) then
            ustar = ustar2
          else
            ustar = ustar_lam
          end if
          exit
        end if
      end do

    end if !laminar/turbulent

    contains

      pure function newguess(ustar)
         !linearize the function by letting the ustar in log constant
         !  and solve the quadratic equation for the larger root
         real(knd) :: newguess
         real(knd),intent(in) :: ustar
         real(knd) :: a,b,c,D
         !u/u_* = dp/dx * z / (k (u_*)**2)  + (1/k) * ln(z * u_* / nu) + B

         a = log(ustar*dist / molecular_viscosity) / 0.41_knd + 5.2_knd
         b = - vel
         c = prgrad * dist / 0.41
         !function to find root of is f(ustar) = a*ustar**2 + b*ustar + c
         D = b**2 - 4*a*c

         if (D<0) then !solution does not exist
           newguess = -1
         else
           newguess = (-b+sqrt(D))/(2*a)
         end if
      end function

  end subroutine WMFlatPrGradUstar




  pure subroutine WMFlatPrGradVisc(visc,ustar,distvect,uvect,walluvect,prgradvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: distvect(3),uvect(3),walluvect(3),prgradvect(3)
    real(knd) :: vect(3),vel,dist,prgrad

    dist = sqrt(sum(distvect**2))

    vect = uvect - walluvect

    vect = vect - dot_product(vect,distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    if (vel>=tiny(1._knd)) then
      !in the same direction as tangential velocity vector
      prgrad = dot_product( prgradvect , vect ) / vel

    else
      !tangential to the wall
      prgrad = sqrt(sum((prgradvect - dot_product(prgradvect,distvect) * distvect / dist**2)**2))

    end if

    call WMFlatPrGradUstar(ustar,vel,prgrad,dist)

    if (ustar<0) ustar = 0

    if (vel>0) then
     if (molecular_viscosity <= 0) then
       visc = 0
     else if (dist*ustar / molecular_viscosity > 1) then
       visc = ustar**2 * dist/vel
     else
       visc = molecular_viscosity
     end if

    else if (molecular_viscosity > 0) then

      visc = molecular_viscosity

    else

      visc = 0

    end if

  end subroutine WMFlatPrGradVisc



















  pure real(knd) function PsiM_MO(zeta)
    real(knd),intent(in):: zeta
    real(knd) :: x

    if (zeta<0) then
      x = (1-15._knd*zeta)**(1/4._knd)
      PsiM_MO = log(((1+x**2)/2._knd)*((1+x)/2._knd)**2)-2._knd*atan(x)+pi/2
    else
      PsiM_MO = - 4.8_knd * zeta
    end if
  end function PsiM_MO


  pure real(knd) function PsiH_MO(zeta)
    real(knd),intent(in):: zeta
    real(knd) :: x

    if (zeta<0) then
      x = (1-15._knd*zeta)**(1/4._knd)
      PsiH_MO = 2._knd*log((1+x**2)/2._knd)
    else
      PsiH_MO = - 7.8_knd * zeta
    end if
  end function PsiH_MO


  pure real(knd) function PsiM_MO_mod(zeta) result(res)
    real(knd),intent(in):: zeta
    real(knd) :: x
    !For L defined without k!
    if (zeta<0) then
      x = (1-6._knd*zeta)**(1/4._knd)
      res = log(((1+x**2)/2._knd)*((1+x)/2._knd)**2)-2._knd*atan(x)+pi/2
    else
      !Zilitinkevich, Esau - 2006
      res = - 3 * zeta**(5._knd/6._knd) 
    end if
  end function PsiM_MO_mod


  pure real(knd) function PsiH_MO_mod(zeta) result(res)
    real(knd),intent(in):: zeta
    real(knd) :: x
    !For L defined without k!
    if (zeta<0) then
      x = (1-6._knd*zeta)**(1/4._knd)
      res = 2._knd*log((1+x**2)/2._knd)
    else
      !Zilitinkevich, Esau - 2006
      res = - 2.5_knd * zeta**(4._knd/5._knd)
    end if
  end function PsiH_MO_mod


  pure real(knd) function Obukhov_zL(ustar,temperature_flux,tempref,g,z)
    real(knd),intent(in):: ustar,temperature_flux,tempref,g,z

    Obukhov_zL = z*(0.4_knd*(g/tempref)*temperature_flux)/(-ustar**3)
  end function Obukhov_zL



  pure subroutine WM_MO_Flux_ustar(vel,dist,ustar,z0,temperature_flux,temperature_ref,grav_acc)
    real(knd),intent(inout) :: ustar
    real(knd),parameter  :: eps = 1e-3
    real(knd),parameter  :: yplcrit = 11.225_knd
    real(knd),intent(in) :: vel,dist,z0,temperature_flux
    real(knd),intent(in) :: temperature_ref,grav_acc
    real(knd) :: ustar2,zL,zL2,Psi
    integer :: i

    if (dist<=z0) then

     if (molecular_viscosity > 0) then

      if ((dist*ustar/molecular_viscosity) < yplcrit) then
        ustar = sqrt(molecular_viscosity * vel / dist)
      else
        ustar = vel / ( log(abs(ustar*dist/molecular_viscosity)) / 0.4_knd + 5.2_knd )
      end if
     else
       ustar = 0
     end if

    else

     i = 0
     zL = 0
     Psi = 0

     do
      i = i+1
      ustar2 = ustar

      ustar = ustar + (max(vel*0.4_knd/(log(max((dist/z0)-Psi,1E-5))), 0._knd) - ustar) / 2

      if (ustar<1E-4) then
       zL = -10000
      else
       zL2 = Obukhov_zL(ustar,temperature_flux,temperature_ref,grav_acc,dist)
       zL = zL+(zL2-zL)/2
      end if

      Psi = PsiM_MO(zL)

      if  (abs(ustar-ustar2)/max(abs(ustar),1.e-3_knd)<eps) exit

      if (i>=50) then
                  ustar = 0
                  exit
      end if

     end do

    end if

  end subroutine WM_MO_Flux_ustar



  pure subroutine WM_MO_Dirichlet_ustar_tfl(ustar,temperature_flux,dist,z0,z0H,vel,tempdif)
    real(knd),intent(inout) :: ustar,temperature_flux
    real(knd),intent(in) :: vel,dist,z0,z0H,tempdif
    real(knd),parameter :: eps = 1e-3_knd
    real(knd),parameter :: yplcrit = 11.225_knd
    real(knd),parameter :: k_U = 0.4_knd
    real(knd),parameter :: k_T = 0.47_knd
    real(knd) :: zL, zL0, psi_m, psi_h, ustar_lam, dist_plus
    integer :: i

    call WMRoughUstar(ustar,vel,dist,z0)

    !if in the laminar region for flat boundary layer, treat as flat and laminar
    if (molecular_viscosity>0) then
      ustar_lam = sqrt(molecular_viscosity * vel / dist)
      dist_plus = dist * ustar_lam / molecular_viscosity
    else
      dist_plus = huge(dist_plus)
    end if
    
    if (dist<=z0 .or. dist_plus < 1.1 *  yplcrit) then

      temperature_flux = - molecular_diffusivity * tempdif / dist

    else
       psi_m = 0
       psi_h = 0
       zL0 = -10000._knd
    
       i = 0
       do
         i = i+1
         
         !L does not contain the Karman constant!
         zL =  - dist * grav_acc * temperature_flux / (temperature_ref * ustar**3)
         
         psi_m = PsiM_MO_mod(zL)
         psi_h = PsiH_MO_mod(zL)
         
         ustar = vel * k_U / (log(dist/z0) - psi_m)
         temperature_flux = - tempdif * ustar * k_T / (log(dist/z0H) - psi_h)
         
         if  (i>1 .and. abs(zL-zL0)/max(abs(zL),1.e-4_knd)<eps) exit
         if (i>=50.or.zL>10000) exit
                 
         zl0 = zL
       end do

       if (i>=50.or.zL>10000) then
         ustar = sqrt(molecular_viscosity * vel / dist)
         temperature_flux = - molecular_diffusivity * tempdif / dist
       end if
    end if
  end subroutine WM_MO_Dirichlet_ustar_tfl




  pure subroutine WM_MO_Dirichlet_ustar_tfl_mfl(ustar, temperature_flux, moisture_flux, &
                                                dist, z0, z0H, &
                                                vel, &
                                                temp, surf_temp, &
                                                moist, surf_moist)
    real(knd),intent(inout) :: ustar, temperature_flux, moisture_flux
    real(knd),intent(in) :: vel, dist, z0, z0H
    real(knd),intent(in) :: temp, surf_temp, moist, surf_moist
    real(knd),parameter :: eps = 1e-3_knd
    real(knd),parameter :: yplcrit = 11.225_knd
    real(knd),parameter :: k_U = 0.4_knd
    real(knd),parameter :: k_T = 0.47_knd, k_Q = k_T
    real(knd),parameter :: e_coef = 0.61_knd
    real(knd) :: zL, zL0, psi_m, psi_h, psi_e
    real(knd) :: tempdif, moistdif, virt_flux
    integer :: i

    tempdif = temp - surf_temp
    moistdif = moist - surf_moist
    
    call WMRoughUstar(ustar,vel,dist,z0)

    if (dist<=z0) then

      temperature_flux = - molecular_diffusivity * tempdif / dist
      moisture_flux = - molecular_diffusivity * moistdif / dist

    else
       psi_m = 0
       psi_h = 0
       psi_e = 0
       zL0 = -10000._knd
    
       i = 0
       do
         i = i+1
         
         !assuming no saturation
         virt_flux = temperature_flux * (1 + e_coef * moist) + &
                     e_coef * temp * moisture_flux
         
         !L does not contain the Karman constant!
         zL =  - dist * grav_acc * virt_flux / (temperature_ref * ustar**3)
         
         psi_m = PsiM_MO_mod(zL)
         psi_h = PsiH_MO_mod(zL)
         psi_e = psi_h
         
         ustar = vel * k_U / (log(dist/z0) - psi_m)
         temperature_flux = - tempdif * ustar * k_T / (log(dist/z0H) - psi_h)
         moisture_flux = - moistdif * ustar * k_Q / (log(dist/z0H) - psi_e)
         
         if  (i>1 .and. abs(zL-zL0)/max(abs(zL),1.e-4_knd)<eps) exit
         if (i>=50.or.zL>10000) exit
                 
         zl0 = zL
       end do

       if (i>=50.or.zL>10000) then
         ustar = sqrt(molecular_viscosity * vel / dist)
         temperature_flux = - molecular_diffusivity * tempdif / dist
       end if
    end if
  end subroutine WM_MO_Dirichlet_ustar_tfl_mfl







  pure subroutine WM_MO_FLUX(visc,ustar,temperature_flux,z0,distvect,uvect)
    real(knd),intent(out) :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: z0,temperature_flux
    real(knd),intent(in)    :: distvect(3),uvect(3)
    real(knd) :: vect(3),vel,dist

    dist = norm2(distvect)

    vect = uvect - dot_product(uvect,distvect) * distvect / dist**2  !tangential part

    vel = norm2(vect)

    if (vel/=0) then
      call WM_MO_FLUX_ustar(vel,dist,ustar,z0,temperature_flux,temperature_ref,grav_acc)
      
      if (ustar<0) ustar = 0
      
      if (ustar*ustar*dist/vel > molecular_viscosity) then
        visc = ustar*ustar*dist/vel
      else if (molecular_viscosity > 0) then
        visc = molecular_viscosity
      else
        visc = 0
      end if
    else
      ustar = 0
      if (molecular_viscosity > 0) then
        visc = molecular_viscosity
      else
        visc = 0
      end if
    end if
    

  end subroutine WM_MO_FLUX



  pure subroutine WM_MO_FluxStress(ustar,temperature_flux,z0,distvect,uvect,tan_vect)
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: z0,temperature_flux
    real(knd),intent(in)    :: distvect(3),uvect(3)
    real(knd),intent(out)   :: tan_vect(3)
    real(knd) :: vel,dist

    call vel_and_dist(tan_vect, dist, uvect, [0._knd,0._knd,0._knd], distvect)

    vel = norm2(tan_vect)

    if (vel/=0) then
      call WM_MO_FLUX_ustar(vel,dist,ustar,z0,temperature_flux,temperature_ref,grav_acc)
      if (ustar<0) ustar = 0
    else
      ustar = 0
    end if

  end subroutine WM_MO_FluxStress



  pure subroutine WM_MO_Dirichlet(visc,ustar,temperature_flux,z0,z0H,tempdif,distvect,uvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar,temperature_flux
    real(knd),intent(in)    :: z0,z0H
    real(knd),intent(in)    :: tempdif ! temperature difference surface - nearest point
    real(knd),intent(in)    :: distvect(3),uvect(3)
    real(knd) :: vect(3),vel,dist

    dist = sqrt(sum(distvect**2))

    vect = uvect - dot_product(uvect,distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    call WM_MO_Dirichlet_ustar_tfl(ustar,temperature_flux,dist,z0,z0H,vel,tempdif)

    if (ustar<0) ustar = 0

    if ((vel>0.or.temperature_flux/=0)) then
      visc = max(ustar*ustar*dist/vel, molecular_viscosity)
    else if (molecular_viscosity > 0) then
      visc = molecular_viscosity
    else
      visc = 0
    end if

  end subroutine WM_MO_Dirichlet


  pure subroutine WM_MO_Dirichlet_moist(visc, &
                                        ustar, temperature_flux, moisture_flux, &
                                        z0, z0H, &
                                        temp, surf_temp, moist, surf_moist, &
                                        distvect, uvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar, temperature_flux, moisture_flux
    real(knd),intent(in)    :: z0, z0H
    real(knd),intent(in)    :: temp, surf_temp, moist, surf_moist
    real(knd),intent(in)    :: distvect(3),uvect(3)
    real(knd) :: vect(3),vel,dist

    dist = sqrt(sum(distvect**2))

    vect = uvect - dot_product(uvect, distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    call WM_MO_Dirichlet_ustar_tfl_mfl(ustar, temperature_flux, moisture_flux, &
                                       dist, z0, z0H, &
                                       vel, temp, surf_temp, moist, surf_moist)

    if (ustar<0) ustar = 0

    if ((vel>0.or.temperature_flux/=0)) then
      visc = max(ustar*ustar*dist/vel, molecular_viscosity)
    else if (molecular_viscosity > 0) then
      visc = molecular_viscosity
    else
      visc = 0
    end if

  end subroutine WM_MO_Dirichlet_moist


  pure subroutine WM_MO_DirichletStress(ustar,temperature_flux,z0,z0H,tempdif,distvect,uvect,tan_vect)
    real(knd),intent(inout) :: ustar,temperature_flux
    real(knd),intent(in)    :: z0,z0H
    real(knd),intent(in)    :: tempdif ! temperature difference surface - nearest point
    real(knd),intent(in)    :: distvect(3),uvect(3)
    real(knd),intent(out)   :: tan_vect(3)
    real(knd) :: vel,dist

    call vel_and_dist(tan_vect, dist, uvect, [0._knd,0._knd,0._knd], distvect)

    vel = sqrt(sum(tan_vect**2))

    call WM_MO_Dirichlet_ustar_tfl(ustar,temperature_flux,dist,z0,z0H,vel,tempdif)

    if (ustar<0) ustar = 0

  end subroutine WM_MO_DirichletStress


  pure subroutine WM_MO_DirichletStress_moist(ustar, temperature_flux, moisture_flux, &
                                              z0, z0H, &
                                              temp, surf_temp, moist, surf_moist, &
                                              distvect, uvect, tan_vect)
    real(knd),intent(inout) :: ustar, temperature_flux, moisture_flux
    real(knd),intent(in)    :: z0, z0H
    real(knd),intent(in)    :: temp, surf_temp
    real(knd),intent(in)    :: moist, surf_moist
    real(knd),intent(in)    :: distvect(3), uvect(3)
    real(knd),intent(out)   :: tan_vect(3)
    real(knd) :: vel,dist

    call vel_and_dist(tan_vect, dist, uvect, [0._knd,0._knd,0._knd], distvect)

    vel = sqrt(sum(tan_vect**2))

    call WM_MO_Dirichlet_ustar_tfl_mfl(ustar, temperature_flux, moisture_flux, &
                                      dist, z0, z0H, &
                                      vel, temp, surf_temp, moist, surf_moist)

    if (ustar<0) ustar = 0

  end subroutine WM_MO_DirichletStress_moist



  pure subroutine bound_temperature_flux(Nu)
    real(knd),intent(inout):: Nu(-1:,-1:)
    integer :: i,j,nx,ny

    nx = Prnx
    ny = Prny

    if (Btype(Ea)==BC_PERIODIC) then
      do j = 1,ny
        Nu(0,j) = Nu(nx,j)
        Nu(nx+1,j) = Nu(1,j)
      end do
    else
      do j = 1,ny
        Nu(0,j) = Nu(1,j)
        Nu(nx+1,j) = Nu(nx,j)
      end do
    end if

    if (Btype(No)==BC_PERIODIC) then
      do i = 1,nx
        Nu(i,0) = Nu(i,ny)
        Nu(i,ny+1) = Nu(i,1)
      end do
    else
      do i = 1,nx
        Nu(i,0) = Nu(i,1)
        Nu(i,ny+1) = Nu(i,ny)
      end do
    end if
  end subroutine bound_temperature_flux




  pure recursive subroutine WallPrGradient(prgrad,i,j,k,Pr,Prtype)
    real(knd),intent(out) :: prgrad(3)
    integer,intent(in)    :: i,j,k
    real(knd),intent(in)  :: Pr(-1:,-1:,-1:)
    integer,intent(in)    :: Prtype(0:,0:,0:)
    integer :: n

    prgrad = 0

    n=0
    if (Prtype(i+1,j,k)>0 .and. i<Prnx) then
      prgrad(1) = prgrad(1) + (Pr(i+1,j,k) - Pr(i,j,k))/(dxU(i))
      n = n + 1
    end if
    if (Prtype(i-1,j,k)>0 .and. i>1) then
      prgrad(1) = prgrad(1) + (Pr(i,j,k) - Pr(i-1,j,k))/(dxU(i-1))
      n = n + 1
    end if
    if (n>0) prgrad(1) = prgrad(1)/n

    n=0
    if (Prtype(i,j+1,k)>0 .and. j<Prny) then
      prgrad(2) = prgrad(2) + (Pr(i,j+1,k) - Pr(i,j,k))/(dyV(j))
      n = n + 1
    end if
    if (Prtype(i,j-1,k)>0 .and. j>1) then
      prgrad(2) = prgrad(2) + (Pr(i,j,k) - Pr(i,j-1,k))/(dyV(j-1))
      n = n + 1
    end if
    if (n>0) prgrad(2) = prgrad(2)/n

    n=0
    if (Prtype(i,j,k+1)>0 .and. k<Prnz) then
      prgrad(3) = prgrad(3) + (Pr(i,j,k+1) - Pr(i,j,k))/(dzW(k))
      n = n + 1
    end if
    if (Prtype(i,j,k-1)>0 .and. k>1) then
      prgrad(3) = prgrad(3) + (Pr(i,j,k) - Pr(i,j,k-1))/(dzW(k-1))
      n = n + 1
    end if
    if (n>0) prgrad(3) = prgrad(3)/n

    prgrad = prgrad + [pr_gradient_x, pr_gradient_y, 0._knd]

  end subroutine WallPrGradient



  subroutine ComputeViscsWM(U,V,W,Pr,Temperature,Moisture)
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    real(knd),dimension(-1:,-1:,-1:),   intent(in) :: Pr
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: Temperature
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: Moisture
    integer :: i, xi, yj, zk
    real(knd) :: tdif
    real(knd) :: dist(3), vel(3), wallvel(3), prgrad(3)
    real(knd) :: visc

    !$omp parallel do private(i,xi,yj,zk,tdif, vel, wallvel, dist, prgrad, visc) schedule(guided,5)
    do i = 1,size(WMPoints)

      xi = WMPoints(i)%xi
      yj = WMPoints(i)%yj
      zk = WMPoints(i)%zk

      dist = [WMPoints(i)%distx, WMPoints(i)%disty, WMPoints(i)%distz]

      vel(1) = (U(xi,yj,zk)+U(xi-1,yj,zk))/2._knd
      vel(2) = (V(xi,yj,zk)+V(xi,yj-1,zk))/2._knd
      vel(3) = (W(xi,yj,zk)+W(xi,yj,zk-1))/2._knd

      wallvel = [WMPoints(i)%wallu, WMPoints(i)%wallv, WMPoints(i)%wallw]


      if (WMPoints(i)%z0>0) then

#ifdef CUSTOM_SURFACE_TEMPERATURE_FLUX
         WMPoints(i)%temperature_flux = SurfaceTemperatureFlux(xPr(xi), yPr(yj), zPr(zk), time_stepping%time)
#endif

         if (WMPoints(i)%prescribed_ustar) then
         
           continue
           
         else if (enable_buoyancy .and. WMPoints(i)%prescribed_temperature) then

#ifdef CUSTOM_SURFACE_TEMPERATURE
           WMPoints(i)%temperature = SurfaceTemperature(xPr(xi), yPr(yj), zPr(zk), time_stepping%time)
#endif
           tdif = Temperature(xi,yj,zk) - WMPoints(i)%temperature

           if (enable_moisture .and. WMPoints(i)%prescribed_temperature) then
#ifdef CUSTOM_SURFACE_MOISTURE
             WMPoints(i)%moisture = SurfaceMoisture(xPr(xi), yPr(yj), zPr(zk), time_stepping%time)
#endif
             
             call WM_MO_Dirichlet_moist(visc, WMPoints(i)%ustar, &
                                  WMPoints(i)%temperature_flux, &
                                  WMPoints(i)%moisture_flux, &
                                  WMPoints(i)%z0, WMPoints(i)%z0H, &
                                  Temperature(xi,yj,zk), WMPoints(i)%temperature, &
                                  Moisture(xi,yj,zk), WMPoints(i)%moisture, &
                                  dist, vel)
           else
             call WM_MO_Dirichlet(visc, WMPoints(i)%ustar, &
                                  WMPoints(i)%temperature_flux, &
                                  WMPoints(i)%z0, WMPoints(i)%z0H, tdif, dist, vel)
           end if


         else if (enable_buoyancy .and. WMPoints(i)%temperature_flux>0) then

           call WM_MO_FLUX(visc, WMPoints(i)%ustar, WMPoints(i)%temperature_flux, &
                           WMPoints(i)%z0, dist, vel)

         else

           call WMRoughVisc(visc, &
                            WMPoints(i)%ustar, WMPoints(i)%z0, &
                            dist, vel, wallvel)
         end if

       else

         if (molecular_viscosity <= 0) then
           call error_stop("The wall model requires positive viscosity or roughness length.")
         end if

         if (wallmodeltype == 2) then

           call WallPrGradient(prgrad,xi,yj,zk,Pr,Prtype)

           call WMFlatPrGradVisc(visc, &
                           WMPoints(i)%ustar, &
                           dist, vel, wallvel, prgrad)

         else

           call WMFlatVisc(visc, &
                           WMPoints(i)%ustar, &
                           dist, vel, wallvel)

         end if
       end if

       Viscosity(xi,yj,zk) = visc

    end do
    !$omp end parallel do

  end subroutine ComputeViscsWM






  subroutine ComputeUVWFluxesWM(U, V, W, Pr, Temperature, Moisture)
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    real(knd),dimension(-1:,-1:,-1:),   intent(in) :: Pr
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: Temperature
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: Moisture


    call fluxes(UxmWMpoints, 1, MINUSX, xU, yPr, zPr)
    call fluxes(UxpWMpoints, 1, PLUSX,  xU, yPr, zPr)
    call fluxes(UymWMpoints, 1, MINUSY, xU, yPr, zPr)
    call fluxes(UypWMpoints, 1, PLUSY,  xU, yPr, zPr)
    call fluxes(UzmWMpoints, 1, MINUSZ, xU, yPr, zPr)
    call fluxes(UzpWMpoints, 1, PLUSZ,  xU, yPr, zPr)

    call fluxes(VxmWMpoints, 2, MINUSX, xPr, yV, zPr)
    call fluxes(VxpWMpoints, 2, PLUSX,  xPr, yV, zPr)
    call fluxes(VymWMpoints, 2, MINUSY, xPr, yV, zPr)
    call fluxes(VypWMpoints, 2, PLUSY,  xPr, yV, zPr)
    call fluxes(VzmWMpoints, 2, MINUSZ, xPr, yV, zPr)
    call fluxes(VzpWMpoints, 2, PLUSZ,  xPr, yV, zPr)

    call fluxes(WxmWMpoints, 3, MINUSX, xPr, yPr, zW)
    call fluxes(WxpWMpoints, 3, PLUSX,  xPr, yPr, zW)
    call fluxes(WymWMpoints, 3, MINUSY, xPr, yPr, zW)
    call fluxes(WypWMpoints, 3, PLUSY,  xPr, yPr, zW)
    call fluxes(WzmWMpoints, 3, MINUSZ, xPr, yPr, zW)
    call fluxes(WzpWMpoints, 3, PLUSZ,  xPr, yPr, zW)


  contains

    subroutine fluxes(points, component, direction, x, y, z)
      type(WMPointUVW), intent(inout), target :: points(:)
      integer, intent(in) :: component, direction
      real(knd), intent(in) :: x(:), y(:), z(:)
      integer :: point, xi, yj, zk
      real(knd) :: dist(3), vel(3), wallvel(3), tan_vect(3), mag, drec(6), temp, moist
      type(WMPointUVW), pointer :: p
      real(knd), parameter :: eps = sqrt(epsilon(1._knd))
      
      drec = [1/dxmin, 1/dxmin, 1/dymin, 1/dymin, 1/dzmin, 1/dzmin]

      !$omp parallel do private(point,xi,yj,zk,dist,vel,wallvel,tan_vect,p,mag) schedule(guided,5)
      do point = 1,size(points)

!         associate (p => points(i))  NOTE: ASSOCIATE supported only in OpenMP 4
          p => points(point)
          xi = p%xi
          yj = p%yj
          zk = p%zk

          dist = [p%distx, p%disty, p%distz]

          vel = local_velocity(U,V,W,component,xi,yj,zk)

          wallvel = [p%wallu, p%wallv, p%wallw]

          if (p%z0>0) then

            if (p%prescribed_ustar) then
            
              continue
              
            else if (enable_buoyancy .and. p%prescribed_temperature) then
             
#ifdef CUSTOM_SURFACE_TEMPERATURE
              p%temperature = SurfaceTemperature(x(xi)+dist(1), &
                                                 y(yj)+dist(2), &
                                                 z(zk)+dist(3), &
                                                 time_stepping%time)
#endif

              temp = local_value(Temperature, component, xi, yj, zk)

              if (enable_moisture) then
               
#ifdef CUSTOM_SURFACE_MOISTURE
                p%moisture = SurfaceMoisture(x(xi)+dist(1), &
                                             y(yj)+dist(2), &
                                             z(zk)+dist(3), &
                                             time_stepping%time)
#endif

                moist = local_value(Moisture, component, xi, yj, zk)

                call WM_MO_DirichletStress_moist(p%ustar, p%temperature_flux, p%moisture_flux, &
                                           p%z0, p%z0H, &
                                           temp, p%temperature, moist, p%moisture, &
                                           dist, vel, tan_vect)
              else
               
                call WM_MO_DirichletStress(p%ustar, p%temperature_flux, &
                                           p%z0, p%z0H, temp-p%temperature, &
                                           dist, vel, tan_vect)
              end if

            else if (enable_buoyancy .and. p%temperature_flux>0) then
!TODO: The temperature flux is always 0 now in momentum points, so no stability effect here!

               call WM_MO_FluxStress(p%ustar, p%temperature_flux, &
                                     p%z0, dist, vel, tan_vect)

            else

               call WMRoughStress(p%ustar, p%z0, &
                                  dist, vel, wallvel, tan_vect = tan_vect)
            end if

          else

            if (molecular_viscosity <= 0) then
               call error_stop("The wall model requires positive viscosity or roughness length.")
            end if


            call WMFlatStress(p%ustar, &
                              dist, vel, wallvel, tan_vect = tan_vect)
          end if

          mag = norm2(tan_vect)

          if (mag>eps) then
          
            ! This is a bit ugly, but saves CPU cycles.
            ! The actual fluxes are scaled by the ratio of the cell boundary area
            !and the cell volume, which is the grid size perpendicular to the boundary.
            ! It also enables the wmfluxes include files to be simpler.
           
            if (gridtype==GRID_VARIABLE_Z .and. direction==MINUSZ) then
            
              p%fluxp = (tan_vect(component) / mag) * p%ustar**2
              p%fluxm = p%fluxp
              if (component==1.or.component==2) then
                p%fluxp = p%fluxp / dzPr(zk)
                p%fluxm = p%fluxm / dzPr(zk-1)
              else
                p%fluxp = p%fluxp / dzW(zk)
                p%fluxm = p%fluxm / dzW(zk-1)
              end if
              
            else if (gridtype==GRID_VARIABLE_Z .and. direction==PLUSZ) then
            
              p%fluxp = (tan_vect(component) / mag) * p%ustar**2
              p%fluxm = p%fluxp
              if (component==1.or.component==2) then
                p%fluxp = p%fluxp / dzPr(zk+1)
                p%fluxm = p%fluxm / dzPr(zk)
              else
                p%fluxp = p%fluxp / dzW(zk+1)
                p%fluxm = p%fluxm / dzW(zk)
              end if
              
            else
            
              p%fluxp = (tan_vect(component) / mag) * p%ustar**2 * drec(direction)
              p%fluxm = p%fluxp
             
            end if

            if (mod(direction,2)==1) then
              p%fluxp = - p%fluxp
              p%fluxm = - p%fluxm
            end if

          else

             p%fluxp = 0
             p%fluxm = p%fluxp

          end if

!         end associate

      end do
      !$omp end parallel do

    end subroutine

  end subroutine ComputeUVWFluxesWM


  pure function local_value(C,component,xi,yj,zk) result(val)
    real(knd) :: val
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: C
    integer, intent(in) :: component, xi, yj, zk

    select case (component)
      case (1)
        val =  (C(xi+1,yj,  zk  ) + C(xi,yj,zk)) / 2
      case (2)
        val =  (C(xi,  yj+1,zk  ) + C(xi,yj,zk)) / 2
      case default
        val =  (C(xi,  yj,  zk+1) + C(xi,yj,zk)) / 2
    end select
  end function


  pure function local_velocity(U,V,W,component,xi,yj,zk) result(vel)
    real(knd) :: vel(3)
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    integer, intent(in) :: component, xi, yj, zk

    select case (component)
      case (1)
        vel(1) =   U(xi,  yj,  zk)
        vel(2) = ( V(xi+1,yj,  zk) + &
                   V(xi,  yj,  zk) + &
                   V(xi+1,yj-1,zk) + &
                   V(xi,  yj-1,zk) &
                 ) / 4._knd
        vel(3) = ( W(xi+1,yj, zk  ) + &
                   W(xi,  yj, zk  ) + &
                   W(xi+1,yj, zk-1) + &
                   W(xi,  yj, zk-1) &
                 ) / 4._knd
      case (2)
        vel(1) = ( U(xi  ,yj+1,zk) + &
                   U(xi,  yj,  zk) + &
                   U(xi-1,yj+1,zk) + &
                   U(xi-1,yj,  zk) &
                 ) / 4._knd
        vel(2) =   V(xi, yj, zk)
        vel(3) = ( W(xi, yj+1,zk  ) + &
                   W(xi, yj,  zk  ) + &
                   W(xi, yj+1,zk-1) + &
                   W(xi, yj,  zk-1) &
                 ) / 4._knd
      case default
        vel(1) = ( U(xi  ,yj, zk+1) + &
                   U(xi,  yj, zk  ) + &
                   U(xi-1,yj, zk+1) + &
                   U(xi-1,yj, zk  ) &
                 ) / 4._knd
        vel(2) = ( V(xi, yj,  zk+1) + &
                   V(xi, yj,  zk  ) + &
                   V(xi, yj-1,zk+1) + &
                   V(xi, yj-1,zk  ) &
                 ) / 4._knd
        vel(3) =   W(xi, yj, zk)
    end select
  end function


  real(knd) function GroundUstar()
    integer :: i, n

    if (offset_to_global_z==0) then
      GroundUstar = 0
      n = 0
      !$omp parallel do private(i) reduction(+:GroundUstar,n) schedule(guided,5)
      do i = 1, size(WMPoints)
        if (WMPoints(i)%zk==1) then
          GroundUstar = GroundUstar + WMPoints(i)%ustar
          n = n + 1
        end if
      end do
      !$omp end parallel do

      GroundUstar = GroundUstar / max(n, 1)
    else
      GroundUstar = 0
    end if
  end function GroundUstar


  real(knd) function GroundUstarUVW() result(res)
    integer :: i, j, k, n

    if (any(WMPoints%zk == 1)) then
      res = 0
      n = 0
      !$omp parallel reduction(+:res,n)
      do j = 1,2
        do i = 1,6
            !$omp do private(k) schedule(guided,5)
            do k = 1, size(WMPointsUVW(i,j)%points)
              if (WMPointsUVW(i,j)%points(k)%zk==1) then
                res = res + WMPointsUVW(i,j)%points(k)%ustar
                n = n + 1
              end if
            end do
            !$omp end do nowait
        end do
      end do
      !$omp end parallel

      res = res / max(n, 1)
    else
      res = 0
    end if
  end function GroundUstarUVW


  pure real(knd) function GroundTFlux()
    if (any(WMPoints%zk == 1)) then
      GroundTFlux = sum(WMPoints%temperature_flux, mask = (WMPoints%zk == 1)) / count(WMPoints%zk == 1)
    else
      GroundTFlux = 0
    end if
  end function GroundTFlux


  real(knd) function GroundTFluxUVW() result(res)
    integer :: i, j, n

    if (any(WMPoints%zk == 1)) then
      res = 0
      n = 0
      do j = 1,2
        do i = 1,6
          associate(p => WMPointsUVW(i,j))
            n = n +  count(p%points%zk == 1)
            res = res + sum(p%points%temperature_flux, mask = (p%points%zk == 1))
          end associate
        end do
      end do
      if (n>0) then
        res = res / n
      else
        res = 0
      end if
    else
      res = 0
    end if
  end function GroundTFluxUVW


  pure real(knd) function GroundMFlux()
    if (any(WMPoints%zk == 1)) then
      GroundMFlux = sum(WMPoints%moisture_flux, mask = (WMPoints%zk == 1)) / count(WMPoints%zk == 1)
    else
      GroundMFlux = 0
    end if
  end function GroundMFlux
  
  
  real(knd) function GroundMFluxUVW() result(res)
    integer :: i, j, n

    if (any(WMPoints%zk == 1)) then
      res = 0
      n = 0
      do j = 1,2
        do i = 1,6
          associate(p => WMPointsUVW(i,j))
            n = n +  count(p%points%zk == 1)
            res = res + sum(p%points%moisture_flux, mask = (p%points%zk == 1))
          end associate
        end do
      end do
      if (n>0) then
        res = res / n
      else
        res = 0
      end if
    else
      res = 0
    end if
  end function GroundMFluxUVW





  pure real(knd) function TotalUstar()
    if (size(WMPoints) > 0) then
      TotalUstar = sum(WMPoints%ustar) / size(WMPoints%zk)
    else
      TotalUstar = 0
    end if
  end function TotalUstar





  pure function GroundDeposition() result(depos)
    real(knd), dimension(:) :: depos(1:Prnx,1:Prny,num_of_scalars)

    integer :: i, j

    depos = 0

    do j = 1, size(WMPoints)

      if (allocated(WMPoints(j)%depscalar)) then

        do i = 1, num_of_scalars
          depos(WMPoints(j)%xi,WMPoints(j)%yj,i) = depos(WMPoints(j)%xi,WMPoints(j)%yj,i) + &
                                                   WMPoints(j)%depscalar(i)
        end do

      end if

    end do
  end function GroundDeposition



#if CUSTOM_SURFACE_TEMPERATURE
  real(knd) function SurfaceTemperature(x,y,z,t)
   real(knd),intent(in):: x,y,z
   real(TIM),intent(in):: t
   interface
     function CustomSurfaceTemperature(x,y,z,t) result(res)
       use Kinds
       real(knd) :: res
       real(knd), intent(in) :: x, y, z
       real(tim), intent(in) :: t
     end function
   end interface
   SurfaceTemperature = CustomSurfaceTemperature(x,y,z,t)
  end function
#endif


#if CUSTOM_SURFACE_TEMPERATURE_FLUX
  real(knd) function SurfaceTemperatureFlux(x,y,z,t)
   real(knd),intent(in):: x,y,z
   real(TIM),intent(in):: t
   interface
     function CustomSurfaceTemperatureFlux(x,y,z,t) result(res)
       use Kinds
       real(knd) :: res
       real(knd), intent(in) :: x, y, z
       real(tim), intent(in) :: t
     end function
   end interface
   SurfaceTemperatureFlux = CustomSurfaceTemperatureFlux(x,y,z,t)
  end function
#endif


 end module Wallmodels
