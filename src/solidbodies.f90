module SolidBodies

  use Parameters
  use Body_class
  use GeometricShapes

  implicit none

  private

  public SolidBody, InitSolidBodies, SetCurrentSB, FindInsideCells, &
         obstacles_file, roughness_file, displacement_file
#ifdef CUSTOMSB
  public  AddSolidBody, SolidBodiesList
#endif

  type, extends(Body) :: SolidBody
    logical   :: rough = .false.                             !T rough surface, F flat surface
    real(knd) :: z0 = 0                                      !roughness parameter
    real(knd) :: z0H = 0                                      !roughness parameter
    real(knd) :: temperature_flux = 0
    real(knd) :: moisture_flux = 0
    !other scalar fluxes assumed zero
    
    !building brick - Santamouris - Environmental Design of Urban Buildings
    real(knd) :: emissivity = 0.45 
    real(knd) :: albedo = 0.3
  end type SolidBody
  
#define TYPEPARAM type(SolidBody)
#include "list-inc-def.f90"

  type(List) :: SolidBodiesList

  type(SolidBody), allocatable, target :: SolidBodiesArray(:)

  character(1024) :: obstacles_file = ''
  
  character(1024) :: roughness_file = ''
  
  character(1024) :: displacement_file = ''
  
  !global bounding box containing all obstacles
  !Can be set by a custom code (CustomSolidBodies) to speed up the search for inside cells.
  real(knd), public :: obstacles_bbox(6) = [-huge(1._knd)/2, huge(1._knd), &
                                            -huge(1._knd)/2, huge(1._knd), &
                                            -huge(1._knd)/2, huge(1._knd)]
  
  interface AddSolidBody
    module procedure AddSolidBody_scalar
    module procedure AddSolidBody_array
  end interface
  
  interface SolidBody
    module procedure SolidBody_Init
  end interface

contains

#include "list-inc-proc.f90"
#undef TYPEPARAM

  subroutine SetCurrentSB(SB, n)
    type(SolidBody), pointer, intent(out) :: SB
    integer,intent(in) :: n

    SB => SolidBodiesArray(n)
  end subroutine SetCurrentSB








  subroutine FindInsideCells
#ifdef PAR
    use custom_par
#endif
    !find if the gridpoints lie inside a solid body and write its number
    !do not nullify the .type arrays, they could have been made nonzero by other unit

    real(knd) :: delta
    integer :: nbody

    delta = (dxmin*dymin*dzmin)**(1._knd/3)/20

    do nbody = 1, size(SolidBodiesArray)
      call SetPrtype(SolidBodiesArray(nbody))
      call SetUtype(SolidBodiesArray(nbody))
      call SetVtype(SolidBodiesArray(nbody))
      call SetWtype(SolidBodiesArray(nbody))
    end do
    
    !set information about numbers of occupied points
    n_free_im   = count(Prtype(1:Prnx,1:Prny,1:Prnz)==0)
    n_free_im_U = count(Utype(1:Unx,1:Uny,1:Unz)==0)
    n_free_im_V = count(Vtype(1:Vnx,1:Vny,1:Vnz)==0)
    n_free_im_W = count(Wtype(1:Wnx,1:Wny,1:Wnz)==0)
    
    n_full_im   = Prnx*Prny*Prnz - n_free_im
    n_full_im_U = Unx*Uny*Unz - n_free_im_U
    n_full_im_V = Vnx*Vny*Vnz - n_free_im_V
    n_full_im_W = Wnx*Wny*Wnz - n_free_im_W
    
#ifdef PAR
    n_free_domain   = par_co_sum(n_free_im)
    n_free_domain_U = par_co_sum(n_free_im_U)
    n_free_domain_V = par_co_sum(n_free_im_V)
    n_free_domain_W = par_co_sum(n_free_im_W)
#else    
    n_free_domain   = n_free_im
    n_free_domain_U = n_free_im_U
    n_free_domain_V = n_free_im_V
    n_free_domain_W = n_free_im_W
#endif

    n_full_domain   = gPrnx*gPrny*gPrnz - n_free_domain
    n_full_domain_U = gUnx*gUny*gUnz - n_free_domain_U
    n_full_domain_V = gVnx*gVny*gVnz - n_free_domain_V
    n_full_domain_W = gWnx*gWny*gWnz - n_free_domain_W

  contains

    subroutine SetPrtype(CurrentSB)
      type(SolidBody) :: CurrentSB
      integer :: i,j,k

      if (CurrentSB%numofbody==0) call error_stop("Error, numofbody==0, did you use AddSolidBody()?")
  
      !$omp parallel do private(i,j,k) schedule(dynamic)
      do k = -2, Prnz+3
        if (zPr(k) > obstacles_bbox(To) .or. zPr(k) < obstacles_bbox(Bo)) cycle
        do j = -2, Prny+3
          if (yPr(j) > obstacles_bbox(No) .or. yPr(j) < obstacles_bbox(So)) cycle
          do i = -2, Prnx+3
            if (xPr(i) > obstacles_bbox(Ea) .or. xPr(i) < obstacles_bbox(We)) cycle
            if (CurrentSB%Inside(xPr(i), yPr(j), zPr(k), delta)) &
                      Prtype(i,j,k) = CurrentSB%numofbody
          enddo
        enddo
      enddo
      !$omp end parallel do
    end subroutine

    subroutine SetUtype(CurrentSB)
      type(SolidBody) :: CurrentSB
      integer :: i,j,k

      if (CurrentSB%numofbody==0) call error_stop("Error, numofbody==0, did you use AddSolidBody()?")
  
      !$omp parallel do private(i,j,k) schedule(dynamic)
      do k = -2, Unz+3
        if (zPr(k) > obstacles_bbox(To) .or. zPr(k) < obstacles_bbox(Bo)) cycle 
        do j = -2, Uny+3
          if (yPr(j) > obstacles_bbox(No) .or. yPr(j) < obstacles_bbox(So)) cycle
          do i = -2, Unx+3
            if (xU(i) > obstacles_bbox(Ea) .or. xU(i) < obstacles_bbox(We)) cycle
            if (CurrentSB%Inside(xU(i), yPr(j) ,zPr(k), delta)) &
                        Utype(i,j,k) = CurrentSB%numofbody
          enddo
        enddo
      enddo
      !$omp end parallel do
    end subroutine

    subroutine SetVtype(CurrentSB)
      type(SolidBody) :: CurrentSB
      integer :: i,j,k

      if (CurrentSB%numofbody==0) call error_stop("Error, numofbody==0, did you use AddSolidBody()?")
  
      !$omp parallel do private(i,j,k) schedule(dynamic)
      do k = -2, Vnz+3
        if (zPr(k) > obstacles_bbox(To) .or. zPr(k) < obstacles_bbox(Bo)) cycle
        do j = -2, Vny+3
          if (yV(j) > obstacles_bbox(No) .or. yV(j) < obstacles_bbox(So)) cycle
          do i = -2, Vnx+3
            if (xPr(i) > obstacles_bbox(Ea) .or. xPr(i) < obstacles_bbox(We)) cycle
            if (CurrentSB%Inside(xPr(i), yV(j), zPr(k), delta)) &
                        Vtype(i,j,k) = CurrentSB%numofbody
          enddo
        enddo
      enddo
      !$omp end parallel do
    end subroutine

    subroutine SetWtype(CurrentSB)
      type(SolidBody) :: CurrentSB
      integer :: i,j,k

      !$omp parallel do private(i,j,k) schedule(dynamic)
      do k = -2, Wnz+3
        if (zW(k) > obstacles_bbox(To) .or. zW(k) < obstacles_bbox(Bo)) cycle
        do j = -2, Wny+3
          if (yPr(j) > obstacles_bbox(No) .or. yPr(j) < obstacles_bbox(So)) cycle
          do i = -2, Wnx+3
            if (xPr(i) > obstacles_bbox(Ea) .or. xPr(i) < obstacles_bbox(We)) cycle
            if (CurrentSB%Inside(xPr(i), yPr(j), zW(k), delta)) &
                        Wtype(i,j,k) = CurrentSB%numofbody
          enddo
        enddo
      enddo
      !$omp end parallel do
    end subroutine

  end subroutine FindInsideCells




  subroutine InitSolidBodies

#ifdef CUSTOMSB
    interface
      subroutine CustomSolidBodies
      end subroutine
    end interface
    !An external subroutine, it should use this module and use AddSolidBody to supply
    ! pointers to the new solid bodies.
    call CustomSolidBodies
#endif

    if (len_trim(obstacles_file)>0) then
      call ReadSolidBodiesFromFile(trim(obstacles_file))
    end if

    call MoveSolidBodiesToArray

  end subroutine InitSolidBodies


  subroutine ReadSolidBodiesFromFile(filename)
    use Strings
    character(*),intent(in) :: filename
    character(5) :: suffix
    
    suffix = filename(index(filename,'.',back=.true.):)

    if (suffix=='.obst'.or.suffix=='.geom') then
      call ReadUnion(filename)
    else if (suffix=='.off') then
      call ReadOff(filename)
    else if (suffix=='.xyz'.or.suffix=='.elv') then
      call ReadTerrain(filename)
    else if (suffix=='.ltop') then
      call ReadTopPoints(filename,.false.)
    else if (suffix=='.rtop') then
      call ReadTopPoints(filename,.true.)
    else
      write(*,*) "Unknown file format "//suffix
      call error_stop
    end if

  end subroutine ReadSolidBodiesFromFile
  
  

  subroutine ReadUnion(filename)
    character(*),intent(in) :: filename
    
    !assume z0 the same as for the lower boundary
    call AddSolidBody(SolidBody(Union(filename), z0 = z0B))
  end subroutine ReadUnion
  
  
  
  subroutine ReadOff(filename)
    character(*),intent(in) :: filename
    logical :: ex

    inquire(file=filename,exist=ex)

    if (ex) then
      !assume z0 the same as for the lower boundary
      call AddSolidBody(SolidBody(Polyhedron(filename), z0 = z0B))
    else
      call error_stop("Error, file "//filename//" does not exist.")
    end if

  end subroutine ReadOff
  
  
  subroutine ReadTopPoints(filename, right)
    character(*),intent(in) :: filename
    logical, intent(in) :: right
    logical :: ex, body_opened
    real(knd), allocatable :: points(:,:)
    integer :: u, stat

    inquire(file=filename,exist=ex)

    if (ex) then
      
      open(newunit=u, file=filename)
      
      body_opened = .false.
      do
        call next_line(stat)
        if (stat>0) exit
      end do
      
      close(u)
    else
      call error_stop("Error, file "//filename//" does not exist.")
    end if

  contains
    subroutine next_line(stat)
      integer, intent(out) :: stat
      character(1024) :: line
      integer :: io
      
      read(u,'(a)', iostat=io) line
      
      if (io/=0) then
        if (body_opened) call close_body
        stat = 1
        return
      end if
      
      if (len_trim(line)==0) then  
        if (body_opened) call close_body
        stat = 3
        return
      end if
        
      call next_point(line, stat)
    end subroutine
    
    subroutine next_point(line, stat)
      character(*), intent(in) :: line
      integer, intent(out) :: stat
      real(knd) :: point(3)
      real(knd), allocatable :: tmp(:,:)
      integer :: io
      
      read(line,*,iostat=io) point
      if (io/=0) then
        stat = 2
        return
      end if
     
      if (.not.body_opened) then
        body_opened = .true.
        allocate(points(3,1))
      else
        tmp = points
        deallocate(points)
        allocate(points(3,size(tmp,2)+1))
        points(:,1:size(tmp,2)) = tmp
      end if
      points(:,size(points,2)) = point
      
      stat = 0
    end subroutine
    
    subroutine close_body
      !assume z0 the same as for the lower boundary
      if (right) then
        call AddSolidBody(SolidBody(ConvexPolyhedron_FromTopPoints(points), z0 = z0B))
      else
        call AddSolidBody(SolidBody(ConvexPolyhedron_FromTopPoints(points(:,size(points,2):1:-1)), z0 = z0B))
      end if

      deallocate(points)
      body_opened = .false.
    end subroutine

  end subroutine ReadTopPoints
  
  
  
  subroutine ReadTerrain(filename)
    use ElevationModels, only: uniform_map
    character(*), intent(in) :: filename

    !assume z0 the same as for the lower boundary if not specified otherwise
    
    if (len_trim(roughness_file)>0) then
      if (len_trim(displacement_file)>0) then
        call AddSolidBody(SolidBody(Terrain(uniform_map(filename), &
                                            uniform_map(trim(roughness_file), default = z0B), &
                                            uniform_map(trim(displacement_file), default = 0._knd))))
      else
        call AddSolidBody(SolidBody(Terrain(uniform_map(trim(filename)), &
                                            uniform_map(trim(roughness_file), default = z0B))))
      end if
    else
      call AddSolidBody(SolidBody(Terrain(uniform_map(filename), z0 = z0B)))
    end if

  end subroutine ReadTerrain


  subroutine MoveSolidBodiesToArray
    integer :: i, nbodies
    type(SolidBody), pointer :: tmp

    nbodies = SolidBodiesList%Len()
    allocate(SolidBodiesArray(nbodies))

    call SolidBodiesList%iter_restart
    do i = 1, nbodies
      tmp => SolidBodiesList%iter_next()
      SolidBodiesArray(i) = tmp
    end do

    call SolidBodiesList%finalize
  end subroutine
  
  
  subroutine AddSolidBody_scalar(SB)
    type(SolidBody),intent(in) :: SB
    !NOTE: expecting numofbody value unspecified and not important in the calling code
    type(SolidBody) :: tmp

    tmp = SB
    tmp%numofbody = SolidBodiesList%Len() + 1

    call SolidBodiesList%add(tmp)

  end subroutine

  subroutine AddSolidBody_array(SB)
    type(SolidBody),intent(in) :: SB(:)
    !NOTE: expecting numofbody value unspecified and not important in the calling code
    integer :: i
    type(SolidBody) :: tmp
    
    do i = 1,size(SB)
      tmp = SB(i)
      tmp%numofbody = SolidBodiesList%Len() + 1

      call SolidBodiesList%add(tmp)
    end do

  end subroutine
  
  function SolidBody_Init(gs, z0, z0H, temperature_flux, moisture_flux) result(res)
    type(SolidBody) :: res
    class(GeometricShape), intent(in) :: gs
    real(knd), optional, intent(in) :: z0, z0H, temperature_flux, moisture_flux
    
    allocate(res%GeometricShape, source=gs)
    if (present(z0)) res%z0 = z0
    if (present(z0H)) then
      res%z0H = z0H
    else
      res%z0H = res%z0
    end if
    res%rough = res%z0 > 0
    if (present(temperature_flux)) res%temperature_flux = temperature_flux
    if (present(moisture_flux)) res%moisture_flux = moisture_flux
  end function


end module SolidBodies
