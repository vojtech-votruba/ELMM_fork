module r3_type
  use Kinds

  type r3
    real(knd) :: x,y,z
  end type
  
  interface assignment(=)
    module procedure r_to_r3
    module procedure v3_to_r3
    module procedure r3_to_v3
  end interface
  
contains
  
  elemental subroutine r_to_r3(lhs,rhs)
    type(r3),intent(out) :: lhs
    real(knd),intent(in) :: rhs
    lhs = r3(rhs,rhs,rhs)
  end subroutine
  
  pure subroutine v3_to_r3(lhs,rhs)
    type(r3),intent(out) :: lhs
    real(knd),intent(in) :: rhs(3)
    lhs = r3(rhs(1),rhs(2),rhs(3))
  end subroutine

  pure subroutine r3_to_v3(lhs,rhs)
    real(knd),intent(out) :: lhs(3)
    type(r3),intent(in)   :: rhs
    lhs = [rhs%x,rhs%y,rhs%z]
  end subroutine
  
  pure function v3(rhs) result(lhs)
    real(knd) :: lhs(3)
    type(r3),intent(in)   :: rhs
    lhs = [rhs%x,rhs%y,rhs%z]
  end function
  
end module r3_type


module GeometricShapes
  use iso_c_binding, only: c_ptr
  use Kinds
  use r3_type
  use Parameters
  use CGAL_Polyhedra
  
  implicit none

  private

  public GeometricShape, Line, Ray, Plane, &
         ConvexPolyhedron, ConvexPolyhedron_FromTopPoints, &
         Polyhedron, &
         Sphere, Ellipsoid, CylJacket, Cylinder, TerrainPoint, Terrain, &
         Translation, Scaling, LinearTransform, Union, ArrayFromObst

  type Bbox
    real(knd) :: xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0
  end type

         
  type,abstract :: GeometricShape
    private
    type(Bbox) :: bbox = Bbox()
  contains
    procedure,private :: in_bbox
    procedure :: Inside => GeometricShape_Inside
    procedure(Inside_interface),private,deferred :: InsideEps
    procedure(Closest_interface), deferred :: Closest
    procedure :: ClosestOut => GeometricShape_ClosestOut
    procedure :: OutwardNormal => GeometricShape_OutwardNormal
    procedure :: IntersectsRay => GeometricShape_IntersectsRay
  end type

  abstract interface
    logical function Inside_interface(self,x,y,z,eps)
      import
      class(GeometricShape),intent(in) :: self
      real(knd),intent(in) :: x,y,z
      real(knd),intent(in) ::eps
    end function
    subroutine Closest_interface(self,xnear,ynear,znear,x,y,z)
      import
      class(GeometricShape),intent(in) :: self
      real(knd),intent(out) :: xnear,ynear,znear
      real(knd),intent(in) :: x,y,z
    end subroutine
  end interface

  type,extends(GeometricShape) :: Line
    real(knd) :: xc,yc,zc
    real(knd) :: a,b,c
  contains
    procedure :: InsideEps => Line_Inside
    procedure :: Closest => Line_Closest
  end type

  type,extends(GeometricShape) :: Ray
    real(knd) :: xc,yc,zc
    real(knd) :: a,b,c
  contains
    procedure :: InsideEps => Ray_Inside
    procedure :: Closest => Ray_Closest
  end type
  
  type,extends(GeometricShape) :: Plane
    real(knd) :: a,b,c,d      !ax+by+cz+d/=0 for inner half-space
    logical gl             !T > in ineq. above F < in ineq. above
   contains
    procedure :: InsideEps => Plane_Inside
    procedure :: Closest => Plane_Closest
  end type


  type,extends(GeometricShape) :: ConvexPolyhedron
    integer :: nplanes = 0
    type(Plane),dimension(:),allocatable :: Planes !intersection of half-spaces
  contains
    procedure,private :: InsideEps => ConvexPolyhedron_Inside
    procedure :: Closest => ConvexPolyhedron_Closest
    procedure :: ClosestOut => ConvexPolyhedron_ClosestOut
  end type


  type,extends(GeometricShape) :: Polyhedron
    type(c_ptr) :: cgalptr
    integer     :: n_ref_points = 2
    type(r3),allocatable    :: ref(:)
  contains
    procedure,private :: ReadOff => Polyhedron_ReadOff
    procedure,private :: InitBbox => Polyhedron_InitBbox
    procedure,private :: InitRefPoints => Polyhedron_InitRefPoints
    procedure,private :: InsideEps => Polyhedron_Inside
    procedure :: Closest => Polyhedron_Closest
    procedure :: IntersectsRay => Polyhedron_IntersectsRay
  end type


  type,extends(GeometricShape) :: Sphere
    real(knd) :: xc,yc,zc,r
  contains
    procedure,private :: InsideEps => Sphere_Inside
    procedure :: Closest => Sphere_Closest
    procedure :: IntersectsRay => Sphere_IntersectsRay
  end type


  type,extends(GeometricShape) :: Ellipsoid
    real(knd) :: xc,yc,zc,a,b,c
  contains
    procedure,private :: InsideEps => Ellipsoid_Inside
    procedure :: Closest => Ellipsoid_Closest
    procedure :: OutwardNormal => Ellipsoid_OutwardNormal
    procedure :: IntersectsRay => Ellipsoid_IntersectsRay
  end type


  type,extends(GeometricShape) :: CylJacket
    real(knd) :: xc,yc,zc
    real(knd) :: a,b,c
    real(knd) :: r
  contains
    procedure :: InsideEps => CylJacket_Inside
    procedure :: Closest => CylJacket_Closest
  end type


  type,extends(GeometricShape) :: Cylinder
    type(CylJacket) Jacket
    type(Plane),allocatable :: Plane1 ,Plane2
  contains
    procedure,private :: InsideEps => Cylinder_Inside
    procedure :: Closest => Cylinder_Closest
    procedure :: ClosestOut => Cylinder_ClosestOut
  end type


  type TerrainPoint
    real(knd) :: elev = 0
    logical :: rough = .false.
    real(knd) :: z0
  end type


  type,extends(GeometricShape) :: Terrain
    type(TerrainPoint),dimension(:,:),allocatable :: UPoints,VPoints,PrPoints !allocate with a buffer of width 1 (i.e. 0:Xnx+1)
  contains
    procedure,private :: InsideEps => Terrain_Inside
    procedure,private,nopass,non_overridable :: GridCoords => Terrain_GridCoords
    procedure :: Closest => Terrain_Closest
    procedure :: LocalPlane => Terrain_LocalPlane
  end type
  
  type,extends(GeometricShape) :: Translation
    class(GeometricShape),allocatable :: original
    type(r3) :: shift
  contains
    procedure,private :: Translation_in_bbox 
    procedure,private :: InsideEps => Translation_Inside
    procedure :: Closest => Translation_Closest
    procedure :: ClosestOut => Translation_ClosestOut
    procedure :: IntersectsRay => Translation_IntersectsRay
  end type
   
  type,extends(GeometricShape) :: Scaling
    class(GeometricShape),allocatable :: original
    type(r3) :: factor
  contains
    procedure,private :: in_bbox => Scaling_in_bbox 
    procedure,private :: InsideEps => Scaling_Inside
    procedure :: Closest => Scaling_Closest
    procedure :: ClosestOut => Scaling_ClosestOut
    procedure :: IntersectsRay => Scaling_IntersectsRay
  end type
   
  type,extends(GeometricShape) :: LinearTransform
    class(GeometricShape),allocatable :: original
    real(knd) :: matrix(3,3), inv_matrix(3,3)
  contains
    procedure,private :: in_bbox => LinearTransform_in_bbox 
    procedure,private :: InsideEps => LinearTransform_Inside
    procedure :: Closest => LinearTransform_Closest
    procedure :: ClosestOut => LinearTransform_ClosestOut
    procedure :: IntersectsRay => LinearTransform_IntersectsRay
  end type
   
  type,extends(GeometricShape) :: Union
    class(GeometricShape),allocatable,private :: items(:)
    integer,private :: size
  contains
    procedure,private :: in_bbox => Union_in_bbox 
    procedure,private :: InsideEps => Union_Inside
    procedure :: Closest => Union_Closest
    procedure :: ClosestOut => Union_ClosestOut
    procedure :: IntersectsRay => Union_IntersectsRay
  end type
   
  interface Closest
    module procedure Line_Closest
  end interface
 
  !initializers
  interface Line
    module procedure Line_Init
  end interface

  interface Ray
    module procedure Ray_Init
    module procedure Ray_Init_v3
    module procedure Ray_Init_r3
  end interface

  interface Plane
    module procedure Plane_Init_3r_32
    module procedure Plane_Init_3r_64
    module procedure Plane_Init_v3_32
    module procedure Plane_Init_v3_64
    module procedure Plane_Init_3xv3_32
    module procedure Plane_Init_3xv3_64
  end interface

  interface ConvexPolyhedron
    module procedure ConvexPolyhedron_Init
  end interface

  interface Polyhedron
    module procedure Polyhedron_Init
  end interface
  
  interface Ellipsoid
    module procedure Ellipsoid_Init
  end interface
  
  interface CylJacket
    module procedure CylJacket_Init
  end interface
  
  interface Cylinder
    module procedure Cylinder_Init
  end interface
  
  interface Terrain
    module procedure Terrain_Init_function
    module procedure Terrain_Init_uniform_DEM
  end interface
  
  interface Translation
    module procedure Translation_Init_3r
    module procedure Translation_Init_3r3
    module procedure Translation_Init_v3
  end interface

  interface Scaling
    module procedure Scaling_Init_r
    module procedure Scaling_Init_r3
    module procedure Scaling_Init_v3
  end interface
  
  interface LinearTransform
    module procedure LinearTransform_Init_scale_r
    module procedure LinearTransform_Init_scale_r3
    module procedure LinearTransform_Init_scale_v3
    module procedure LinearTransform_Init_rot
  end interface
  
  interface Union
    !TODO change to a polymorphic function that changes result according to data read
    module procedure Union_Init
    module procedure Union_Init_File
  end interface
  
  interface GeometricShape
    module procedure GeometricShape_FromFile
  end interface
  
  type top_points_container
    real(knd), allocatable :: points(:,:)
  end type
  
  real(knd),parameter :: unit_matrix_3(3,3) = reshape(source=[0,0,1,0,1,0,0,0,1], &
                                                      shape=[3,3])

contains

  !defaults
  
  logical function GeometricShape_Inside(self,x,y,z,eps) result(ins)
    class(GeometricShape),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),optional,intent(in) ::eps
    
    if (present(eps)) then
      ins = self%InsideEps(x,y,z,eps)
    else
      ins = self%InsideEps(x,y,z,10*epsilon(1._knd))
    end if
    
  end function
  
  subroutine GeometricShape_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(GeometricShape),intent(in) :: self
    real(knd),intent(out)  :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%Closest(xnear,ynear,znear,x,y,z)
  end subroutine
    
  function GeometricShape_OutwardNormal(self,x,y,z) result(res)
    real(knd) :: res(3)
    class(GeometricShape),intent(in) :: self
    real(knd),intent(in) :: x,y,z

    res = huge(res)
    call error_stop("OutwardNormal not implemented for this shape")
  end function
    
  logical function GeometricShape_IntersectsRay(self,r) result(intersects)
    class(GeometricShape),intent(in) :: self
    class(Ray),intent(in) :: r
    
    !shapes will be transparent for (solar) rays if they do not override this method.
    intersects = .false.
 
  end function
  

  
  
  !helpers
  
  real(knd) function PointDist(x1,y1,z1,x2,y2,z2)
    real(knd),intent(in) :: x1,y1,z1,x2,y2,z2

    PointDist = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
  end function
  
   real(knd) function LineDist(x,y,z,xl,yl,zl,a,b,c)
    real(knd),intent(in) :: x,y,z,xl,yl,zl,a,b,c
    real(knd) :: t

    if (((a/=0).or.(b/=0)).or.(c/=0)) then
     t = (a*(x-xl)+b*(y-yl)+c*(z-zl))/(a**2+b**2+c**2)
    else
     t = 0
    endif

    LineDist = sqrt((xl+a*t-x)**2+(yl+b*t-y)**2+(zl+c*t-z)**2)
  end function LineDist
 
  
  

  
  logical function in_bbox(self,x,y,z,eps)
    class(GeometricShape),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    
    associate(b=>self%bbox)
      in_bbox = (x+eps>=b%xmin).and. &
                (x-eps<=b%xmax).and. &
                (y+eps>=b%ymin).and. &
                (y-eps<=b%ymax).and. &
                (z+eps>=b%zmin).and. &
                (z-eps<=b%zmax)
    end associate
  end function
  
  
  

  
  



  function Line_Init(xc,yc,zc,a,b,c) result(res)
    type(Line) :: res
    real(knd),intent(in) :: xc, yc, zc, a, b, c
    
    res%xc = xc
    res%yc = yc
    res%zc = zc
    res%a  = a
    res%b  = b
    res%c  = c

  end function
  
  logical function Line_Inside(self,x,y,z,eps) result(ins)
    class(Line),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps

    ins = .false.
  end function
  

  subroutine Line_Closest(self,xnear,ynear,znear,x,y,z)
    class(Line),intent(in) :: self
    real(knd),intent(out)  :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: t

    if (self%a/=0 .or. self%b/=0 .or. self%c/=0) then
      t = ( self%a*(x-self%xc) + self%b*(y-self%yc) + self%c*(z-self%zc) ) / (self%a**2 + self%b**2 + self%c**2)
    else
      t = 0
    endif

    xnear = self%xc + self%a * t
    ynear = self%yc + self%b * t
    znear = self%zc + self%c * t
  end subroutine


  
  
  
  
  
  
  
  function Ray_Init(xc,yc,zc,a,b,c) result(res)
    type(Ray) :: res
    real(knd),intent(in) :: xc, yc, zc, a, b, c
    
    res%xc = xc
    res%yc = yc
    res%zc = zc
    res%a  = a
    res%b  = b
    res%c  = c

  end function
  
  function Ray_Init_v3(c,v) result(res)
    type(Ray) :: res
    real(knd),intent(in) :: c(3),v(3)
    
    res%xc = c(1)
    res%yc = c(2)
    res%zc = c(3)
    res%a  = v(1)
    res%b  = v(2)
    res%c  = v(3)

  end function
  
  function Ray_Init_r3(c,v) result(res)
    type(Ray) :: res
    type(r3),intent(in) :: c,v
    
    res%xc = c%x
    res%yc = c%y
    res%zc = c%z
    res%a  = v%x
    res%b  = v%y
    res%c  = v%z

  end function
  
  
  logical function Ray_Inside(self,x,y,z,eps) result(ins)
    class(Ray),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps

    ins = .false.
  end function
  
  
  subroutine Ray_Closest(self,xnear,ynear,znear,x,y,z)
    class(Ray),intent(in) :: self
    real(knd),intent(out)  :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: t

    if (self%a/=0 .or. self%b/=0 .or. self%c/=0) then
      t = ( self%a*(x-self%xc) + self%b*(y-self%yc) + self%c*(z-self%zc) ) / (self%a**2 + self%b**2 + self%c**2)
    else
      t = 0
    endif

    xnear = self%xc + self%a * max(t,0._knd)
    ynear = self%yc + self%b * max(t,0._knd)
    znear = self%zc + self%c * max(t,0._knd)
  end subroutine


  
  
  
  
  
  
  
  
  
  
  function Plane_Init_3r(xc,yc,zc,a,b,c) result(res)
    !Point and an outward normal
    type(Plane) :: res
    real(knd),intent(in) :: xc, yc, zc, a, b, c
    
    res%a  = a
    res%b  = b
    res%c  = c
    res%gl = .false.
    res%d = - (a*xc + b*yc + c*zc)

  end function
  
  function Plane_Init_3r_32(xc,yc,zc,a,b,c) result(res)
    !Point and an outward normal
    type(Plane) :: res
    real(real32),intent(in) :: xc, yc, zc, a, b, c
    
    res = Plane_Init_3r(real(xc,knd), real(yc,knd), real(zc,knd), &
                        real(a,knd),  real(b,knd),  real(c,knd))
  end function
  
  function Plane_Init_3r_64(xc,yc,zc,a,b,c) result(res)
    !Point and an outward normal
    type(Plane) :: res
    real(real64),intent(in) :: xc, yc, zc, a, b, c
    
    res = Plane_Init_3r(real(xc,knd), real(yc,knd), real(zc,knd), &
                        real(a,knd),  real(b,knd),  real(c,knd))
  end function
  
  function Plane_Init_v3_32(point, vec) result(res)
    !Point and an outward normal
    type(Plane) :: res
    real(real32),intent(in) :: point(3), vec(3)
    
    res = Plane_Init_3r(real(point(1),knd), real(point(2),knd), real(point(3),knd), &
                        real(vec(1),knd),  real(vec(2),knd),  real(vec(3),knd))
  end function
  
  function Plane_Init_v3_64(point, vec) result(res)
    !Point and an outward normal
    type(Plane) :: res
    real(real64),intent(in) :: point(3), vec(3)
    
    res = Plane_Init_3r(real(point(1),knd), real(point(2),knd), real(point(3),knd), &
                        real(vec(1),knd),  real(vec(2),knd),  real(vec(3),knd))
  end function
  
  function Plane_Init_3xv3_knd(point1, point2, point3) result(res)
    !three points in righthand order
    use ArrayUtilities, only: cross_product
    type(Plane) :: res
    real(knd),intent(in) :: point1(3), point2(3), point3(3)
    real(knd) :: vec(3)
    
    vec = cross_product(point3-point2, point1-point2)
    
    if (norm2(vec)<epsilon(vec)) then
      write(*,*) "Error, in line ",__LINE__," in file ",__FILE__
      write(*,*) "Points",point1, point2, point3,"lie on a single line." 
      call error_stop
    end if
    
    res = Plane(point2, vec)
  end function
  
  function Plane_Init_3xv3_32(point1, point2, point3) result(res)
    !three points in righthand order
    use ArrayUtilities, only: cross_product
    type(Plane) :: res
    real(real32),intent(in) :: point1(3), point2(3), point3(3)
    
    res = Plane_Init_3xv3_knd(real(point1, knd), &
                              real(point2, knd), &
                              real(point3, knd))
  end function
  
  function Plane_Init_3xv3_64(point1, point2, point3) result(res)
    !three points in righthand order
    use ArrayUtilities, only: cross_product
    type(Plane) :: res
    real(real64),intent(in) :: point1(3), point2(3), point3(3)
    
    res = Plane_Init_3xv3_knd(real(point1, knd), &
                              real(point2, knd), &
                              real(point3, knd))
  end function
  
  logical function Plane_Inside(self,x,y,z,eps) result(ins)
    class(Plane),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps

    if (self%GL) then
      ins = self%a*x+self%b*y+self%c*z+self%d >= -eps*sqrt(self%a**2+self%b**2+self%c**2)
    else
      ins = self%a*x+self%b*y+self%c*z+self%d <= eps*sqrt(self%a**2+self%b**2+self%c**2)
    endif  
  end function

  subroutine Plane_Closest(self,xnear,ynear,znear,x,y,z)
    class(Plane),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) ::x,y,z
    real(knd) :: t

    if (abs(self%a)>tiny(1._knd).or. &
        abs(self%b)>tiny(1._knd).or. &
        abs(self%c)>tiny(1._knd)) then
      t = -(self%a*x+self%b*y+self%c*z+self%d)/sqrt(self%a**2+self%b**2+self%c**2)
    else
      t = 0
    endif
    xnear = x+self%a*t
    ynear = y+self%b*t
    znear = z+self%c*t
  end subroutine Plane_Closest



  
  
  
  
  
  
  
  function ConvexPolyhedron_Init(Planes) result(res)
    !Point and an outward normal
    type(ConvexPolyhedron) :: res
    type(Plane),intent(in) :: Planes(:)
    !limitation of gfortran 4.8 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=44672
    allocate(res%Planes(size(Planes)), source=Planes)
    res%nplanes = size(Planes)

  end function

  
  function ConvexPolyhedron_FromTopPoints(Points) result(res)
    !points along the upper flat boundary in the left-hand order (outward normal upwards)
    type(ConvexPolyhedron) :: res
    real(knd), intent(in) :: Points(:,:)
    integer :: i, n
    
    n = size(Points,2)
    
    allocate(res%Planes(n+1))
    res%nplanes = n+1
    
    do i = 1, n-1
      res%Planes(i) = Plane(Points(:,i+1), Points(:,i), [Points(1:2,i), Points(3,i)-(gzmax-gzmin)])
    end do
    res%Planes(n) = Plane(Points(:,1), Points(:,n), [Points(1:2,n), Points(3,n)-(gzmax-gzmin)])
    
    res%Planes(n+1) = Plane(Points(:,1), Points(:,2), Points(:,3))
    
  end function
  
  logical function ConvexPolyhedron_Inside(self,x,y,z,eps) result(ins)
    class(ConvexPolyhedron),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    integer :: i

    if (self%nplanes>0) then
      ins = .true.
      do i = 1,self%nplanes

       ins = self%Planes(i)%Inside(x,y,z,eps)

       if (.not. ins) exit
      enddo
    else
      ins = .false.
    endif
  end function
  
  subroutine ConvexPolyhedron_Closest(self,xnear,ynear,znear,x,y,z)
    class(ConvexPolyhedron),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    real(knd) :: dists(self%nplanes),xP(self%nplanes),yP(self%nplanes),zP(self%nplanes)
    real(knd) :: minv
    real(knd) :: ailine,biline,ciline
    real(knd) :: x0iline,y0iline,z0iline,xln,yln,zln,p
    integer   :: inearest,inearest2,inearest3
    integer   :: i
    real(knd) :: xg(3), ag(3,3), bg(3)
    real(knd) :: rcond


   !Vzdalenosti od rovin, pokud je nejbl. bod roviny uvnitr jine, nebo puv. bod na vnitrni strane -> vzd. *-1
   !Pokud je nejbl. body +- eps. (norm. vektor!,dxmin/100) uvnitr a vne polyh -> hotovo
   !Jinak 2. nejbl. rovina v abs. hodnote -> prusecnice a nejbl bod na ni
   !Nejbl. bod na prusecnici. Pokud +-eps. uvnitr, (najit vekt. v rovine  kolme na prusecnic)-:hotovo
   !Jinak iterativne najit bod na prusecnici uvnitr

    do i = 1,self%nplanes
     call self%Planes(i)%Closest(xP(i),yP(i),zP(i),x,y,z)
     dists(i) = sqrt((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
     if (self%Planes(i)%Inside(x,y,z)) dists(i) = -ABS(dists(i))
    enddo
    !find nearest plane with
    inearest = 0
    minv = huge(minv)
    do i = 1,self%nplanes
     !>=0 caused problems near the corners when on the extended wall, it accepted negative zero
     if (dists(i)>=tiny(dxmin).and.dists(i)<minv) then
       inearest = i
       minv = dists(i)
     endif
    enddo

    if (inearest==0) then
      write(*,*) "no nearest"
      write(*,*) "dists", dists
      write(*,*) "line ",__LINE__," in file ",__FILE__
      call error_stop
    endif

    if (self%Inside(xP(inearest),yP(inearest),zP(inearest), &
                    MIN(dxmin/10000._knd,dymin/10000._knd,dzmin/10000._knd) &
                   )) then
     xnear = xP(inearest)
     ynear = yP(inearest)
     znear = zP(inearest)
     return
    endif

    dists = abs(dists)

    inearest2 = 0
    minv = huge(minv)
    do i = 1,self%nplanes
     if (i/=inearest.and.dists(i)<minv) then
       inearest2 = i
       minv = dists(i)
     endif
    enddo

    inearest3 = 0
    minv = huge(minv)
    do i = 1,self%nplanes
     if (i/=inearest.and.i/=inearest2.and.dists(i)<minv) then
       inearest3 = i
       minv = dists(i)
     endif
    enddo

    ailine = self%Planes(inearest)%b*self%Planes(inearest2)%c-self%Planes(inearest)%c*self%Planes(inearest2)%b
    biline = self%Planes(inearest)%c*self%Planes(inearest2)%a-self%Planes(inearest)%a*self%Planes(inearest2)%c
    ciline = self%Planes(inearest)%a*self%Planes(inearest2)%b-self%Planes(inearest)%b*self%Planes(inearest2)%a

    if ( abs(ailine)<=epsilon(ailine).and. &
         abs(biline)<=epsilon(biline).and. &
         abs(ciline)<=epsilon(ciline) )   then
     write(*,*) "cross product 0"
     write(*,*) "numbers of planes",inearest,inearest2
     write(*,*) "x,y,z",x,y,z
     write(*,*) "pl1", self%Planes(inearest)%a, self%Planes(inearest)%b, &
                       self%Planes(inearest)%c, self%Planes(inearest)%d
     write(*,*) "pl2", self%Planes(inearest2)%a, self%Planes(inearest2)%b, &
                       self%Planes(inearest2)%c, self%Planes(inearest2)%d

     call self%Planes(inearest)%Closest(xln,yln,zln,x,y,z)
     call self%Planes(inearest2)%Closest(xnear,ynear,znear,x,y,z)

     if ((xln-x)**2+(yln-y)**2+(zln-z)**2<(xln-x)**2+(yln-y)**2+(zln-z)**2) then
       xnear = xln
       ynear = yln
       znear = zln
     end if

     return

    endif

    p = sqrt(ailine**2+biline**2+ciline**2)
    ailine = ailine/p
    biline = biline/p
    ciline = ciline/p

    if (abs(ciline)>=0.1_knd) then

       ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,0._knd, &
                           self%Planes(inearest)%b,self%Planes(inearest2)%b,0._knd, &
                           self%Planes(inearest)%c,self%Planes(inearest2)%c,1._knd /), &
                 shape = (/ 3,3 /))

       bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,(zW(Wnz+1)+zW(0))/2._knd /)

       call solve3x3(ag, bg, xg, rcond)

       if (rcond>1e6) write(*,*) "Warning, ill-conditioned matrix in ConvexPolyhedron_Closest!"

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    elseif (abs(biline)>=0.1_knd) then

       ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,0._knd, &
                           self%Planes(inearest)%b,self%Planes(inearest2)%b,1._knd, &
                           self%Planes(inearest)%c,self%Planes(inearest2)%c,0._knd /), &
                 shape = (/ 3,3 /))

       bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,(yV(Vny+1)+yV(0))/2._knd /)

       call solve3x3(ag, bg, xg, rcond)

       if (rcond>1e6) write(*,*) "Warning, ill-conditioned matrix in ConvexPolyhedron_Closest!"

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    else

       ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,1._knd, &
                           self%Planes(inearest)%b,self%Planes(inearest2)%b,0._knd, &
                           self%Planes(inearest)%c,self%Planes(inearest2)%c,0._knd /), &
                  shape = (/ 3,3 /))

       bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,(xU(Unx+1)+xU(0))/2._knd /)

       if (rcond>1e6) write(*,*) "Warning, ill-conditioned matrix in ConvexPolyhedron_Closest!"

       call solve3x3(ag, bg, xg, rcond)

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    endif

    call Closest(Line(x0iline,y0iline,z0iline,ailine,biline,ciline),xln,yln,zln,x,y,z)

    if (self%Inside(xln,yln,zln,min(dxmin/1000._knd,dymin/10000._knd,dzmin/10000._knd))) then
      xnear = xln
      ynear = yln
      znear = zln
      return
    endif



    ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,self%Planes(inearest3)%a, &
                        self%Planes(inearest)%b,self%Planes(inearest2)%b,self%Planes(inearest3)%b, &
                        self%Planes(inearest)%c,self%Planes(inearest2)%c,self%Planes(inearest3)%c /), &
               shape = (/ 3,3 /))
    bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,-self%Planes(inearest3)%d /)

    if (rcond>1e6) write(*,*) "Warning, ill-conditioned matrix in ConvexPolyhedron_Closest!"

    call solve3x3(ag, bg, xg, rcond)

    xnear = xg(1)
    ynear = xg(2)
    znear = xg(3)

  contains

    subroutine solve3x3(a, B, x, det)
      ! Cramers rule from https://groups.google.com/forum/#!topic/comp.lang.fortran/LkZYhHlFIz0
      ! given a(3,3) and B(3), return x(3) such that a*x=B.
      ! If det==0 a is singular, and x is meaningless.
      real(knd), intent(in)  :: a(3,3), b(3)
      real(knd), intent(out) :: x(3), det
      real(knd) :: ainv(3,3)

      ainv(1,1) = a(2,2)*a(3,3) - a(2,3)*a(3,2)
      ainv(1,2) = a(3,2)*a(1,3) - a(3,3)*a(1,2)
      ainv(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
      ainv(2,1) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
      ainv(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
      ainv(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
      ainv(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
      ainv(3,2) = a(3,1)*a(1,2) - a(3,2)*a(1,1)
      ainv(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

      det = a(1,1)*ainv(1,1) + a(1,2)*ainv(2,1) + a(1,3)*ainv(3,1)

      if (abs(det) < 1e-8) then
        det = 1e8
        return
      end if

      det = 1 / det
      x(1) = det * (ainv(1,1)*b(1) + ainv(1,2)*b(2) + ainv(1,3)*b(3))
      x(2) = det * (ainv(2,1)*b(1) + ainv(2,2)*b(2) + ainv(2,3)*b(3))
      x(3) = det * (ainv(3,1)*b(1) + ainv(3,2)*b(2) + ainv(3,3)*b(3))
    end subroutine

  end subroutine ConvexPolyhedron_Closest

  subroutine ConvexPolyhedron_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(ConvexPolyhedron),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: dists(self%nplanes),xP(self%nplanes),yP(self%nplanes),zP(self%nplanes),minv
    integer :: inearest,i

    dists = huge(minv)
    do i = 1,self%nplanes
     call self%Planes(i)%Closest(xP(i),yP(i),zP(i),x,y,z)
     dists(i) = sqrt((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
    enddo

    inearest = 0
    minv = huge(minv)
    do i = 1,self%nplanes
     if (dists(i)>=0.and.dists(i)<minv) then
       inearest = i
       minv = dists(i)
     endif
    enddo

    xnear = xP(inearest)
    ynear = yP(inearest)
    znear = zP(inearest)
  end subroutine


  
  
  
  
  
  
  
  function Polyhedron_Init(filename) result (res)
    type(Polyhedron) :: res
    character(*) :: filename
    
    call res%ReadOff(filename)
    call res%InitBbox
    
  end function
  
  
  subroutine Polyhedron_ReadOff(self,filename)
    !reads geometry from an .off file
    use iso_c_binding, only: c_ptr,c_associated
    class(Polyhedron),intent(out) :: self
    character(*),intent(in) :: filename

    call cgal_polyhedron_read(self%cgalptr, filename)
  
    if (.not.c_associated(self%cgalptr)) then
      write(*,*) "Error reading polyhedron from ",filename
      call error_stop
    end if
  end subroutine
  
  subroutine Polyhedron_InitBbox(self)
    class(Polyhedron),intent(inout) :: self

    associate(b=>self%bbox)
      call cgal_polyhedron_bbox(self%cgalptr, b%xmin, b%ymin, b%zmin, b%xmax, b%ymax, b%zmax)

      call self%InitRefPoints

      b%xmin = max(b%xmin, xU(-2))
      b%ymin = max(b%ymin, yV(-2))
      b%zmin = max(b%zmin, zW(-2))
      b%xmax = min(b%xmax, xU(Prnx+2))
      b%ymax = min(b%ymax, yV(Prny+2))
      b%zmax = min(b%zmax, zW(Prnz+2))

    end associate
  end subroutine
  
  subroutine Polyhedron_InitRefPoints(self)
    class(Polyhedron),intent(inout) :: self
    real(knd) :: p(3)
    integer :: i
    
    allocate(self%ref(self%n_ref_points))
    
    do i = 1, self%n_ref_points
      call random_number(p)
      associate(b=>self%bbox)
        self%ref(i) = r3(b%xmin+(b%xmax-b%xmin)*p(1), &
                         b%ymin+(b%ymax-b%ymin)*p(2), &
                         b%zmax + (0.1)*(b%zmax-b%zmin))
      end associate
    end do
  end subroutine
    
  logical function Polyhedron_Inside(self,x,y,z,eps) result(ins)
    class(Polyhedron),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    real(knd) :: xn, yn, zn

    if (self%in_bbox(x,y,z,eps)) then
         
      ins = &
          cgal_polyhedron_inside(self%cgalptr, &
                                  x,y,z, &
                                  x,y,gzmax+(gzmax-gzmin)/10)

      if (eps>0) then                            
        if (.not.ins) then
          call self%Closest(xn, yn, zn, x, y, z)
          ins = norm2([xn-x, yn-y, zn-z])<eps
        end if
        
        if (.not.ins) ins = &
            cgal_polyhedron_inside(self%cgalptr, &
                                    x-eps,y,z, &
                                    x,y,gzmax+(gzmax-gzmin)/10)
        if (.not.ins) ins = &
            cgal_polyhedron_inside(self%cgalptr, &
                                    x+eps,y,z, &
                                    x,y,gzmax+(gzmax-gzmin)/10)
        if (.not.ins) ins = &
            cgal_polyhedron_inside(self%cgalptr, &
                                    x,y-eps,z, &
                                    x,y,gzmax+(gzmax-gzmin)/10)
        if (.not.ins) ins = &
            cgal_polyhedron_inside(self%cgalptr, &
                                    x,y-eps,z, &
                                    x,y,gzmax+(gzmax-gzmin)/10)
        if (.not.ins) ins = &
            cgal_polyhedron_inside(self%cgalptr, &
                                    x,y,z-eps, &
                                    x,y,gzmax+(gzmax-gzmin)/10)
        if (.not.ins) ins = &
            cgal_polyhedron_inside(self%cgalptr, &
                                    x,y,z+eps, &
                                    x,y,gzmax+(gzmax-gzmin)/10)
                 
      end if
      
      if (ins) then
        ins = .true.
        return
      end if
          
      
      ins = .false.
      
    else
      ins = .false.
    end if
  
  end function Polyhedron_Inside

  subroutine Polyhedron_Closest(self,xnear,ynear,znear,x,y,z)
    class(Polyhedron),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call cgal_polyhedron_closest(self%cgalptr, &
                                 x,y,z, &
                                 xnear, ynear, znear)
  
  end subroutine
  
  logical function Polyhedron_IntersectsRay(self,r) result(intersects)
    class(Polyhedron),intent(in) :: self
    class(Ray),intent(in) :: r

      intersects = &
           cgal_polyhedron_intersects_ray(self%cgalptr, &
                                          r%xc, r%yc, r%zc, &
                                          r%a,  r%b,  r%c )
  end function



  
  
  
  
  
  
  
  
  
  

  logical function Sphere_Inside(self,x,y,z,eps) result(ins)
    class(Sphere),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps

    ins = (self%xc-x)**2 + (self%yc-y)**2 + (self%zc-z)**2 <= (self%r+eps)**2
  end function

  subroutine Sphere_Closest(self,xnear,ynear,znear,x,y,z)
    class(Sphere),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: t,a,b,c

    a = x - self%xc
    b = y - self%yc
    c = z - self%zc
    t = self%r / sqrt(a**2+b**2+c**2)
    xnear = a*t + self%xc
    ynear = b*t + self%yc
    znear = c*t + self%zc
  end subroutine
  
  logical function Sphere_IntersectsRay(self,r) result(res)
     class(Sphere),intent(in) :: self
     class(Ray),intent(in) :: r
     real(knd) :: rc(3),rv(3) !transformed ray center and vector
     real(knd) :: a,b,c,D,t1,t2
     
     rc = [r%xc - self%xc, r%yc - self%yc, r%zc - self%zc]
     rv = [r%a, r%b, r%c]
     
     a = dot_product(rv,rv)
     b = 2 * dot_product(rc,rv)
     c = dot_product(rc,rc) - self%r**2
     
     D = b**2 - 4*a*c
     
     if (D<0) then
       res = .false.
     else
       t1 = (-b - sqrt(D)) / (2*a)
       t2 = (-b + sqrt(D)) / (2*a)
       res = (t1>=0 .or. t2>=0)
     end if
  end function
  
  
  
  
  
  
  
  
  

  function Ellipsoid_Init(xc,yc,zc,a,b,c) result(res)
    type(Ellipsoid) :: res
    real(knd),intent(in) :: xc, yc, zc, a, b, c
    
    res%xc = xc
    res%yc = yc
    res%zc = zc
    res%a  = a
    res%b  = b
    res%c  = c
    
    res%bbox%xmin = xc - a
    res%bbox%xmax = xc + a
    res%bbox%ymin = yc - b
    res%bbox%ymax = yc + b
    res%bbox%zmin = zc - c
    res%bbox%zmax = zc + c    

  end function
  
  logical function Ellipsoid_Inside(self,x,y,z,eps) result(ins)
    class(Ellipsoid),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    
    if (self%in_bbox(x,y,z,eps)) then
      ins =  ((self%xc-x)**2/(self%a)**2 + &
            (self%yc-y)**2/(self%b)**2 + &
            (self%zc-z)**2/(self%c)**2 <= 1._knd+sqrt(eps))
    else
        ins = .false.
    end if
  end function

  subroutine Ellipsoid_Closest(self,xnear,ynear,znear,x,y,z)
    class(Ellipsoid),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: t,a,b,c !auxiliary ray parameters

    ! NOT EXACT!
    a = (x - self%xc)/self%a
    b = (y - self%yc)/self%b
    c = (z - self%zc)/self%c
    t = 1._knd / sqrt(a**2+b**2+c**2)
    xnear = self%a * a*t + self%xc
    ynear = self%b * b*t + self%yc
    znear = self%c * c*t + self%zc
  end subroutine
  
  function Ellipsoid_OutwardNormal(self,x,y,z) result(res)
    real(knd) :: res(3)
    class(Ellipsoid),intent(in) :: self
    real(knd),intent(in) :: x,y,z

    res(1) = 2 * (x - self%xc) / self%a**2
    res(2) = 2 * (y - self%yc) / self%b**2
    res(3) = 2 * (z - self%zc) / self%c**2

    if (norm2(res)>0) res = res / norm2(res)
  end function
  
  logical function Ellipsoid_IntersectsRay(self,r) result(res)
     class(Ellipsoid),intent(in) :: self
     class(Ray),intent(in) :: r
     type(Ray) :: r2
     !FIXME: constructors contain irrelevant properties and no keywords due to problem in ifort as of version 14
     type(Sphere),parameter :: unit_sphere = Sphere( Bbox(0,0,0,0,0,0), 0._knd,  0._knd,&
                                                      0._knd,  1._knd)

     r2 = Ray((r%xc - self%xc)/self%a, (r%yc - self%yc)/self%b, (r%zc - self%zc)/self%c, &
                  r%a/self%a, r%b/self%b, r%c/self%c )

     res = unit_sphere%IntersectsRay(r2)
     
  end function
  


  function CylJacket_Init(point, vec, radius) result(res)
    type(CylJacket) :: res
    real(knd),intent(in) :: point(3), vec(3), radius
    
    res%xc = point(1)
    res%yc = point(2)
    res%zc = point(3)
    res%a = vec(1)
    res%b = vec(2)
    res%c = vec(3)
    res%r = radius
  end function



  function Cylinder_Init(point1, point2, radius, infinite) result(res)
    type(Cylinder) :: res
    real(knd),intent(in) :: point1(3), point2(3), radius
    logical,intent(in),optional :: infinite
    logical :: finite
    
    res%jacket = CylJacket(point1, point2-point1, radius)
    
    finite = .true.
    if (present(infinite)) then
      if (.not.infinite) then
        finite = .false.
      end if
    end if
    
    if (finite) then
      res%plane1 = Plane(point1, point1-point2)
      res%plane2 = Plane(point2, point2-point1)
    end if

  end function
  
  logical function CylJacket_Inside(self,x,y,z,eps) result(ins)
    class(CylJacket),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps

    ins = LineDist(x,y,z,self%xc,self%yc,self%zc,self%a,self%b,self%c) <= self%r+eps
  end function

  subroutine CylJacket_Closest(self,xnear,ynear,znear,x,y,z)
    class(CylJacket),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: t,xl,yl,zl,a,b,c

    call Closest(Line(self%xc,self%yc,self%zc,self%a,self%b,self%c),xl,yl,zl,x,y,z)

    a = x-xl
    b = y-yl
    c = z-zl
    t = self%r/sqrt(a**2+b**2+c**2)

    xnear = a*t+xl
    ynear = b*t+yl
    znear = c*t+zl
  end subroutine

  
  
  
  
  logical function Cylinder_Inside(self,x,y,z,eps) result(ins)
    class(Cylinder),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps

    ins = .true.

    if (.not.self%Jacket%Inside(x,y,z,eps)) ins = .false.

    if (ins.and.allocated(self%Plane1)) then
           if (.not.self%Plane1%Inside(x,y,z,eps)) ins = .false.
    endif
    if (ins.and.allocated(self%Plane2)) then
          if (.not.self%Plane2%Inside(x,y,z,eps)) ins = .false.
    endif
  end function

  subroutine Cylinder_Closest(self,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
   class(Cylinder),intent(in) :: self
   real(knd),intent(out) :: xnear,ynear,znear
   real(knd),intent(in) :: x,y,z
   real(knd) :: xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

   !!!Only for Planes perpendicular to jacket!!!!

   if (allocated(self%Plane1)) then
       if (allocated(self%Plane2)) then
          call self%Jacket%Closest(xJ,yJ,zJ,x,y,z)
          call self%Plane1%Closest(xP1,yP1,zP1,x,y,z)
          call self%Plane2%Closest(xP2,yP2,zP2,x,y,z)
          if (self%Jacket%Inside(x,y,z)) then
             if (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
             else
                xnear = xP2
                ynear = yP2
                znear = zP2
             endif
          elseif (self%Plane1%Inside(x,y,z).and.self%Plane2%Inside(x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          elseif (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
           call self%Jacket%Closest(xnear,ynear,znear,xP1,yP1,zP1)
          else
           call self%Jacket%Closest(xnear,ynear,znear,xP2,yP2,zP2)
          endif
       else
          call self%Jacket%Closest(xJ,yJ,zJ,x,y,z)
          call self%Plane1%Closest(xP1,yP1,zP1,x,y,z)
          if (self%Jacket%Inside(x,y,z)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
          elseif (self%Plane1%Inside(x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          else
           call self%Jacket%Closest(xnear,ynear,znear,xP1,yP1,zP1)
         endif
       endif

    else

      call self%Jacket%Closest(xnear,ynear,znear,x,y,z)

    endif
  end subroutine Cylinder_Closest

  subroutine Cylinder_ClosestOut(self,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
    class(Cylinder),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

    if (allocated(self%Plane1)) then
      call self%Plane1%Closest(xP1,yP1,zP1,x,y,z)
    else
     xP1 = sqrt(huge(xP1))/10
     yP1 = sqrt(huge(yP1))/10
     zP1 = sqrt(huge(zP1))/10
    endif
    if (allocated(self%Plane2)) then
      call self%Plane2%Closest(xP2,yP2,zP2,x,y,z)
    else
     xP2 = sqrt(huge(xP2))/10
     yP2 = sqrt(huge(yP2))/10
     zP2 = sqrt(huge(zP2))/10
    endif
    call self%Jacket%Closest(xJ,yJ,zJ,x,y,z)

    if (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2).and.&
         PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xJ,yJ,zJ)) then
                 xnear = xP1
                 ynear = yP1
                 znear = zP1
    elseif (PointDist(x,y,z,xP2,yP2,zP2)<=PointDist(x,y,z,xP1,yP1,zP1).and.&
             PointDist(x,y,z,xP2,yP2,zP2)<=PointDist(x,y,z,xJ,yJ,zJ)) then
                 xnear = xP2
                 ynear = yP2
                 znear = zP2
    elseif (PointDist(x,y,z,xJ,yJ,zJ)<=PointDist(x,y,z,xP2,yP2,zP2).and.&
             PointDist(x,y,z,xJ,yJ,zJ)<=PointDist(x,y,z,xP1,yP1,zP1)) then
                 xnear = xJ
                 ynear = yJ
                 znear = zJ
    endif
  end subroutine Cylinder_ClosestOut


  
  
  
  
  
  
  
  subroutine Terrain_GridCoords(x2,y2,xi,yj,comp)
    real(knd),intent(in) :: x2,y2
    integer,intent(out) :: xi,yj,comp
    real(knd) :: x,y,distPr,distU,distV
    integer :: xPri,yPrj,xUi,yVj,i

    x = x2
    y = y2

    if (gridtype==GRID_UNIFORM) then

        xPri = min(max(nint( (x - xU(0))/dxmin + 0.5_knd ),1),Prnx+1)

        yPrj = min(max(nint( (y - yV(0))/dymin + 0.5_knd ),1),Prny+1)

        xUi = min(max(nint( (x-xU(0))/dxmin ),0), Unx+1)

        yVj = min(max(nint( (y-yV(0))/dymin ),0), Vny+1)
    else


        xPri = Prnx+1
        do i=0,Prnx+1
         if (xU(i)>=x) then
                       xPri = i
                       exit
                      endif
        enddo

        yPrj = Prny+1
        do i=0,Prny+1
         if (yV(i)>=y) then
                       yPrj = i
                       exit
                      endif
        enddo

        xUi = Prnx+1
        do i = 0,Prnx+1
         if (xPr(i+1)>=x) then
                       xUi = i
                       exit
                      endif
        enddo

        yVj = Prny+1
        do i = 0,Prny+1
         if (yPr(i+1)>=y) then
                       yVj = i
                       exit
                      endif
        enddo

    endif

    distPr = hypot(x-xPr(xPri), y-yPr(yPrj))
    distU = hypot(x-xU(xUi), y-yPr(yPrj))
    distV = hypot(x-xPr(xPri), y-yV(yVj))

    if (distU<distPr.and.distU<distV) then
     xi = xUi
     yj = yPrj
     comp = 1
    elseif (distV<distPr.and.distV<distU) then
     xi = xPri
     yj = yVj
     comp = 2
    else
     xi = xPri
     yj = yPrj
     comp = 3
    endif
  end subroutine Terrain_GridCoords
  
  
  pure function Terrain_LocalPlane(self, xi, yj, comp) result(res)
    type(Plane) :: res
    class(Terrain), intent(in) :: self
    integer, intent(in) :: xi, yj, comp
    real(knd) :: a, b, zloc
    
    if (comp==1) then  !Construct a tangent plane using first diferences

      zloc = self%UPoints(xi,yj)%elev
      
      if (xi<=0) then
        a = (self%PrPoints(1,yj)%elev - self%PrPoints(0,yj)%elev) / dxmin
      else if (xi>=Prnx+1) then
        a = (self%PrPoints(Prnx+1,yj)%elev - self%PrPoints(Prnx,yj)%elev) / dxmin
      else
        a = (self%PrPoints(xi+1,yj)%elev - self%PrPoints(xi,yj)%elev) / dxmin
      end if
      
      if (yj<=0) then
        b = (self%UPoints(xi,1)%elev - self%UPoints(xi,0)%elev) / dymin
      else if (yj>=Uny+1) then
        b = (self%UPoints(xi,Uny+1)%elev - self%UPoints(xi,Uny)%elev) / dymin
      else
        b = (self%UPoints(xi,yj+1)%elev - self%UPoints(xi,yj-1)%elev) / (2*dymin)
      end if

      res%a = a
      res%b = b
      res%c = -1._knd
      res%d = - a*xU(xi) - b*yPr(yj) + zloc

    elseif (comp==2) then

      zloc = self%VPoints(xi,yj)%elev
      
      if (xi<=0) then
        a = (self%VPoints(1,yj)%elev - self%VPoints(0,yj)%elev) / dxmin
      else if (xi>=Vnx+1) then
        a = (self%VPoints(Vnx+1,yj)%elev - self%VPoints(Vnx,yj)%elev) / dxmin
      else
        a = (self%VPoints(xi+1,yj)%elev - self%VPoints(xi-1,yj)%elev) / (2*dxmin)
      end if
      
      if (yj<=0) then
        b = (self%PrPoints(xi,1)%elev - self%PrPoints(xi,0)%elev) / dymin
      else if (yj>=Prny+1) then
        b = (self%PrPoints(xi,Prny+1)%elev - self%PrPoints(xi,Prny)%elev) / dymin
      else
        b = (self%PrPoints(xi,yj+1)%elev - self%PrPoints(xi,yj)%elev) / dymin
      end if

      res%a = a
      res%b = b
      res%c = -1._knd
      res%d = - a*xPr(xi) - b*yV(yj) + zloc

    elseif (comp==3) then

      zloc = self%PrPoints(xi,yj)%elev
      
      if (xi<=0) then
        a = (self%UPoints(1,yj)%elev - self%UPoints(0,yj)%elev) / dxmin
      else if (xi>=Unx+1) then
        a = (self%UPoints(Unx+1,yj)%elev - self%UPoints(Unx,yj)%elev) / dxmin
      else
        a = (self%UPoints(xi,yj)%elev - self%UPoints(xi-1,yj)%elev) / dxmin
      end if
      
      
      if (yj<=0) then
        b = (self%VPoints(xi,1)%elev - self%VPoints(xi,0)%elev) / dymin
      else if (yj>=Vny+1) then
        b = (self%VPoints(xi,Vny+1)%elev - self%VPoints(xi,Vny)%elev) / dymin
      else
        b = (self%VPoints(xi,yj)%elev - self%VPoints(xi,yj-1)%elev) / dymin
      end if

      res%a = a
      res%b = b
      res%c = -1._knd
      res%d = - a*xPr(xi) - b*yPr(yj) + zloc

    endif

    res%gl = .true.
  end function Terrain_LocalPlane
  
  logical function Terrain_Inside(self,x,y,z,eps)
    class(Terrain),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    integer :: xi, yj, comp
    type(Plane) :: Pl
    
    call self%GridCoords(x,y,xi,yj,comp)

    Pl = self%LocalPlane(xi,yj,comp)
    
    Terrain_Inside = Pl%Inside(x, y, z, eps)
  end function

  subroutine Terrain_Closest(self,xnear,ynear,znear,x,y,z)
    class(Terrain),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    type(Plane) :: Pl
    integer     :: xi,yj,comp

    xnear = x
    ynear = y
    znear = z

    call self%GridCoords(x,y,xi,yj,comp)

    Pl = self%LocalPlane(xi,yj,comp)
    
    call Pl%Closest(xnear,ynear,znear,x,y,z)
    
  end subroutine Terrain_Closest
  
  
  
  function Terrain_Init_function(fz, z0) result(res)
    type(Terrain) :: res
    interface
      real(knd) function fz(x,y)
        import
        real(knd), intent(in) :: x,y
      end function
    end interface
    real(knd), optional, intent(in) :: z0
    real(knd) :: z0_l
    logical :: rough
     integer :: i,j
   
    if (present(z0)) then
      z0_l = z0
      rough = .true.
    else
      if (z0B > 0) then
        rough = .true.
        z0_l = z0B
      else
        rough = .false.
        z0_l = 0
      end if
    end if
    
    allocate(res%PrPoints(0:Prnx+1,0:Prny+1))
    allocate(res%UPoints(0:Unx+1,0:Uny+1))
    allocate(res%VPoints(0:Vnx+1,0:Vny+1))
    
    do j = 0, Prny+1
      do i = 0, Prnx+1
        res%PrPoints(i,j)%elev = fz(xPr(i), yPr(j))
        res%PrPoints(i,j)%rough = .true.
        res%PrPoints(i,j)%z0 = z0_l
      end do
    end do
    
    do j = 0, Uny+1
      do i = 0, Unx+1
        res%UPoints(i,j)%elev = fz(xU(i), yPr(j))
        res%UPoints(i,j)%rough = rough
        res%UPoints(i,j)%z0 = z0_l
      end do
    end do
    
    do j = 0, Vny+1
      do i = 0, Vnx+1
        res%VPoints(i,j)%elev = fz(xPr(i), yV(j))
        res%VPoints(i,j)%rough = rough
        res%VPoints(i,j)%z0 = z0_l
      end do
    end do
    
    res%bbox%xmin = xPr(0)
    res%bbox%xmax = xPr(Prnx+1)
    res%bbox%ymin = yPr(0)
    res%bbox%ymax = yPr(Prny+1)
    res%bbox%zmin = -huge(1.)
    res%bbox%zmax = max(maxval(res%PrPoints%elev), &
                        maxval(res%UPoints%elev), &
                        maxval(res%VPoints%elev))
    
  end function Terrain_Init_function



  
  function Terrain_Init_uniform_DEM(elevation, z0_map, displacement_map, z0) result(res)
    use ElevationModels, only: map
    type(Terrain) :: res
    class(map), intent(in) :: elevation
    class(map), optional, intent(in) :: z0_map
    class(map), optional, intent(in) :: displacement_map
    real(knd), optional, intent(in) :: z0
    real(knd) :: z0_l
    logical :: rough
    integer :: i,j
    
    if (present(z0_map)) then
      rough = .true.
      z0_l = -1
    else if (present(z0)) then
      z0_l = z0
      rough = .true.
    else
      if (z0B > 0) then
        rough = .true.
        z0_l = z0B
      else
        rough = .false.
        z0_l = 0
      end if
    end if
    
    allocate(res%PrPoints(0:Prnx+1,0:Prny+1))
    allocate(res%UPoints(0:Unx+1,0:Uny+1))
    allocate(res%VPoints(0:Vnx+1,0:Vny+1))
    
    do j = 0, Prny+1
      do i = 0, Prnx+1
        res%PrPoints(i,j)%elev = elevation%value(xPr(i), yPr(j))
        
        if (present(displacement_map)) &
          res%PrPoints(i,j)%elev = res%PrPoints(i,j)%elev + &
                                     displacement_map%value(xPr(i), yPr(j))
        
        res%PrPoints(i,j)%rough = rough

        if (present(z0_map)) then
          res%PrPoints(i,j)%z0 = z0_map%value(xPr(i), yPr(j))
        else
          res%PrPoints(i,j)%z0 = z0_l
        end if
      end do
    end do
    
    do j = 0, Uny+1
      do i = 0, Unx+1
        res%UPoints(i,j)%elev = elevation%value(xU(i), yPr(j))
        
        if (present(displacement_map)) &
          res%UPoints(i,j)%elev = res%UPoints(i,j)%elev + &
                                     displacement_map%value(xU(i), yPr(j))
                                     
        res%UPoints(i,j)%rough = rough

        if (present(z0_map)) then
          res%UPoints(i,j)%z0 = z0_map%value(xU(i), yPr(j))
        else
          res%UPoints(i,j)%z0 = z0_l
        end if
      end do
    end do
    
    do j = 0, Vny+1
      do i = 0, Vnx+1
        res%VPoints(i,j)%elev = elevation%value(xPr(i), yV(j))
        
        if (present(displacement_map)) &
          res%VPoints(i,j)%elev = res%VPoints(i,j)%elev + &
                                     displacement_map%value(xPr(i), yV(j))
                                     
        res%VPoints(i,j)%rough = rough

        if (present(z0_map)) then
          res%VPoints(i,j)%z0 = z0_map%value(xPr(i), yV(j))
        else
          res%VPoints(i,j)%z0 = z0_l
        end if
      end do
    end do
    
    res%bbox%xmin = xPr(0)
    res%bbox%xmax = xPr(Prnx+1)
    res%bbox%ymin = yPr(0)
    res%bbox%ymax = yPr(Prny+1)
    res%bbox%zmin = -huge(1.)
    res%bbox%zmax = max(maxval(res%PrPoints%elev), &
                        maxval(res%UPoints%elev), &
                        maxval(res%VPoints%elev))
                        
  end function Terrain_Init_uniform_DEM




  
  logical function Translation_in_bbox(self,x,y,z,eps) result(in)
    class(Translation),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    
    in = self%original%in_bbox(x - self%shift%x, & 
                               y - self%shift%y, &
                               z - self%shift%z, &
                               eps)
  end function


  logical function Translation_Inside(self,x,y,z,eps) result(ins)
    class(Translation),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps
    
    ins = self%original%InsideEps(x - self%shift%x, & 
                                  y - self%shift%y, &
                                  z - self%shift%z, &
                                  eps)
    
  end function
  
  
  subroutine Translation_Closest(self,xnear,ynear,znear,x,y,z)
    class(Translation),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%original%Closest(xnear, &
                               ynear, &
                               znear, &
                               x - self%shift%x, & 
                               y - self%shift%y, &
                               z - self%shift%z)
    xnear = xnear + self%shift%x
    ynear = ynear + self%shift%y
    znear = znear + self%shift%z
  end subroutine
  
  subroutine Translation_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(Translation),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%original%ClosestOut(xnear, &
                                  ynear, &
                                  znear, &
                                  x - self%shift%x, & 
                                  y - self%shift%y, &
                                  z - self%shift%z)
    xnear = xnear + self%shift%x
    ynear = ynear + self%shift%y
    znear = znear + self%shift%z
  end subroutine
  
  
  logical function Translation_IntersectsRay(self,r) result(res)
    class(Translation),intent(in) :: self
    class(Ray),intent(in) :: r
    
    res = self%original%IntersectsRay(Ray(r%xc - self%shift%x, &
                               r%yc - self%shift%y, &
                               r%zc - self%shift%z, &
                               r%a, &
                               r%b, &
                               r%c))
  end function
  
  elemental function Translation_Init_3r(original,sx,sy,sz) result(res)
    type(Translation) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: sx,sy,sz
   
    allocate(res%original, source=original)
    res%shift = [sx,sy,sz]
  end function
  
  elemental function Translation_Init_3r3(original,shift) result(res)
    type(Translation) :: res
    class(GeometricShape),intent(in) :: original
    type(r3),intent(in) :: shift
   
    allocate(res%original, source=original)
    res%shift = shift
  end function
  
  function Translation_Init_v3(original,shift) result(res)
    type(Translation) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: shift(3)
   
    allocate(res%original, source=original)
    res%shift = shift
  end function
  
  
  
  logical function Scaling_in_bbox(self,x,y,z,eps) result(in)
    class(Scaling),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    
    in = self%original%in_bbox(x / self%factor%x, & 
                               y / self%factor%y, &
                               z / self%factor%z, &
                               eps / (self%factor%x*self%factor%y*self%factor%z)**(1._knd/3))
  end function


  logical function Scaling_Inside(self,x,y,z,eps) result(ins)
    class(Scaling),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps

    ins = self%original%InsideEps(x / self%factor%x, & 
                                  y / self%factor%y, &
                                  z / self%factor%z, &
                                  eps / (self%factor%x*self%factor%y*self%factor%z)**(1._knd/3))
    
  end function
  
  
  subroutine Scaling_Closest(self,xnear,ynear,znear,x,y,z)
    class(Scaling),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%original%Closest(xnear, &
                               ynear, &
                               znear, &
                               x / self%factor%x, & 
                               y / self%factor%y, &
                               z / self%factor%z)
    xnear = xnear * self%factor%x
    ynear = ynear * self%factor%y
    znear = znear * self%factor%z
  end subroutine
  
  subroutine Scaling_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(Scaling),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%original%ClosestOut(xnear, &
                                  ynear, &
                                  znear, &
                                  x / self%factor%x, & 
                                  y / self%factor%y, &
                                  z / self%factor%z)
    xnear = xnear * self%factor%x
    ynear = ynear * self%factor%y
    znear = znear * self%factor%z
  end subroutine
  
  
  logical function Scaling_IntersectsRay(self,r) result(res)
    class(Scaling),intent(in) :: self
    class(Ray),intent(in) :: r
    
    res = self%original%IntersectsRay(Ray(r%xc / self%factor%x, &
                               r%yc / self%factor%y, &
                               r%zc / self%factor%z, &
                               r%a, &
                               r%b, &
                               r%c))
  end function
  
  elemental function Scaling_Init_r(original,scalar_factor) result(res)
    type(Scaling) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: scalar_factor
   
    allocate(res%original, source=original)
    res%factor = scalar_factor
  end function
  
  elemental function Scaling_Init_r3(original,factor) result(res)
    type(Scaling) :: res
    class(GeometricShape),intent(in) :: original
    type(r3),intent(in) :: factor
   
    allocate(res%original, source=original)
    res%factor = factor
  end function
  
  function Scaling_Init_v3(original,factor) result(res)
    type(Scaling) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: factor(3)
   
    allocate(res%original, source=original)
    res%factor = factor
  end function
  

  
  
  logical function LinearTransform_in_bbox(self,x,y,z,eps) result(in)
    class(LinearTransform),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    real(knd) :: xyz(3)
    
    xyz = matmul(self%inv_matrix, [x,y,z])
    
    in = self%original%in_bbox(xyz(1), & 
                               xyz(2), &
                               xyz(3), &
                               eps)
  end function


  logical function LinearTransform_Inside(self,x,y,z,eps) result(ins)
    class(LinearTransform),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps
    real(knd) :: xyz(3)
    
    xyz = matmul(self%inv_matrix, [x,y,z])
    
    ins = self%original%InsideEps(xyz(1), & 
                                  xyz(2), &
                                  xyz(3), &
                                  eps)
    
  end function
  
  
  subroutine LinearTransform_Closest(self,xnear,ynear,znear,x,y,z)
    class(LinearTransform),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: xyz(3),xyznear(3)

    xyz = matmul(self%inv_matrix, [x,y,z])
    
    call self%original%Closest(xyznear(1), &
                               xyznear(2), &
                               xyznear(3), &
                               xyz(1), & 
                               xyz(2), &
                               xyz(3))
                               
    xyznear = matmul(self%matrix, xyznear)
                        
    xnear = xyznear(1)
    ynear = xyznear(2)
    znear = xyznear(3)
  end subroutine
  
  subroutine LinearTransform_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(LinearTransform),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: xyz(3),xyznear(3)

    xyz = matmul(self%inv_matrix, [x,y,z])
    
    call self%original%ClosestOut(xyznear(1), &
                               xyznear(2), &
                               xyznear(3), &
                               xyz(1), & 
                               xyz(2), &
                               xyz(3))
                               
    xyznear = matmul(self%matrix, xyznear)
                        
    xnear = xyznear(1)
    ynear = xyznear(2)
    znear = xyznear(3)
  end subroutine
  
  
  logical function LinearTransform_IntersectsRay(self,r) result(res)
    class(LinearTransform),intent(in) :: self
    class(Ray),intent(in) :: r
    real(knd) :: xyz(3),abc(3)
    
    xyz = matmul(self%inv_matrix, [r%xc, r%yc, r%zc])
    abc = matmul(self%inv_matrix, [r%a, r%b, r%c])
    
    res = self%original%IntersectsRay(Ray(xyz, abc))
  end function
  
  elemental function LinearTransform_Init_scale_r(original,scalar_factor) result(res)
    type(LinearTransform) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: scalar_factor
   
    allocate(res%original, source=original)
    
    res%matrix = unit_matrix_3 * scalar_factor
    res%inv_matrix = unit_matrix_3 / scalar_factor
  end function
  
  function LinearTransform_Init_scale_v3(original,factor) result(res)
    type(LinearTransform) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: factor(3)
    integer :: i
   
    allocate(res%original, source=original)
    
    res%matrix = 0
    res%inv_matrix = 0

    forall(i=1:3)  res%matrix(i,i) = factor(i)
    forall(i=1:3)  res%inv_matrix(i,i) = 1._knd/factor(i)
  end function
  
  elemental function LinearTransform_Init_scale_r3(original,factor) result(res)
    type(LinearTransform) :: res
    class(GeometricShape),intent(in) :: original
    type(r3),intent(in) :: factor
    real(knd) :: v(3)
    integer :: i
   
    allocate(res%original, source=original)
    
    v = factor
    
    res%matrix = 0
    res%inv_matrix = 0

    forall(i=1:3)  res%matrix(i,i) = v(i)
    forall(i=1:3)  res%inv_matrix(i,i) = 1._knd/v(i)
  end function
  
  pure function rotation_matrix_x(phi) result(res)
    real(knd) :: res(3,3)
    real(knd),intent(in) :: phi
    real(knd) :: c,s
    
    c = cos(phi)
    s = sin(phi)
    res(1,1) = 1
    res(2,1) = 0
    res(3,1) = 0
    res(1,2) = 0
    res(2,2) = c
    res(3,2) = s
    res(1,3) = 0
    res(2,3) = -s
    res(3,3) = c
  end function
  
  pure function rotation_matrix_y(phi) result(res)
    real(knd) :: res(3,3)
    real(knd),intent(in) :: phi
    real(knd) :: c,s
    
    c = cos(phi)
    s = sin(phi)
    res(1,1) = c
    res(2,1) = 0
    res(3,1) = -s
    res(1,2) = 0
    res(2,2) = 1
    res(3,2) = 0
    res(1,3) = s
    res(2,3) = 0
    res(3,3) = c
  end function
  
  pure function rotation_matrix_z(phi) result(res)
    real(knd) :: res(3,3)
    real(knd),intent(in) :: phi
    real(knd) :: c,s
    
    c = cos(phi)
    s = sin(phi)
    res(1,1) = c
    res(2,1) = s
    res(3,1) = 0
    res(1,2) = -s
    res(2,2) = c
    res(3,2) = 0
    res(1,3) = 0
    res(2,3) = 0
    res(3,3) = 1
  end function
  
  
  elemental function LinearTransform_Init_rot(original,axis,phi) result(res)
    type(LinearTransform) :: res
    class(GeometricShape),intent(in) :: original
    integer,intent(in) :: axis
    real(knd),intent(in) :: phi
  
    allocate(res%original, source=original)

    if (axis==1) then
      res%matrix = rotation_matrix_x(phi)
      res%inv_matrix = rotation_matrix_x(-phi)
    else if (axis==2) then
      res%matrix = rotation_matrix_y(phi)
      res%inv_matrix = rotation_matrix_y(-phi)
    else if (axis==3) then
      res%matrix = rotation_matrix_z(phi)
      res%inv_matrix = rotation_matrix_z(-phi)
    end if
  end function


  
  
  logical function Union_in_bbox(self,x,y,z,eps) result(in)
    class(Union),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    integer :: i
    
    in = any([ ( self%items(i)%in_bbox(x,y,z,eps), i=1,self%size ) ])

  end function


  logical function Union_Inside(self,x,y,z,eps) result(ins)
    class(Union),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps
    integer :: i
    
    ins = any([ ( self%items(i)%InsideEps(x,y,z,eps), i=1,self%size ) ])
    
  end function
  
  
  subroutine Union_Closest(self,xnear,ynear,znear,x,y,z)
    class(Union),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: xs(self%size),ys(self%size),zs(self%size)
    integer :: i

    do i=1,self%size
      call self%items(i)%Closest(xs(i),ys(i),zs(i),x,y,z)
    end do
    
    associate (j => minloc(hypot(xs,hypot(ys,zs)), dim=1))
      xnear = xs(j)
      ynear = ys(j)
      znear = zs(j)
    end associate
  end subroutine
  
  subroutine Union_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(Union),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: xs(self%size),ys(self%size),zs(self%size)
    integer :: i
    !FIXME: this will work assuming we are close enough to the boundary
    do i=1,self%size
      call self%items(i)%Closest(xs(i),ys(i),zs(i),x,y,z)
    end do
    
    associate (j => minloc(hypot(xs,hypot(ys,zs)), dim=1))
      xnear = xs(j)
      ynear = ys(j)
      znear = zs(j)
    end associate
  end subroutine  
  
  logical function Union_IntersectsRay(self,r) result(res)
    class(Union),intent(in) :: self
    class(Ray),intent(in) :: r
    integer :: i
    
    res = any([ ( self%items(i)%IntersectsRay(Ray(r%xc,r%yc,r%zc,r%a,r%b,r%c)), &
                  i=1,self%size ) ])
  end function
  
  function Union_Init(items) result(res)
     type(Union) :: res
     class(GeometricShape),intent(in) :: items(:)
    
     allocate(res%items(size(items)), source=items)
     res%size = size(res%items)
  end function
  
  
  function Union_Init_File(filename) result(res)
    type(Union) :: res    
    character(*),intent(in) :: filename
    integer :: l
    
    l = len_trim(filename)
    if (filename(l-4:l)=='.obst') then
      res = Union_Init_Obst(filename)
    else if (filename(l-4:l)=='.geom') then
      res = Union_Init_Geom(filename)
    else if (filename(l-4:l)=='.ltop') then
      res = Union_Init_TopPoints(filename,.false.)
    else if (filename(l-4:l)=='.rtop') then
      res = Union_Init_TopPoints(filename,.true.)
    end if
  end function
  
  
  function Union_Init_TopPoints(filename, right) result(res)
    use Strings, only: upcase
    type(Union) :: res    
    character(*),intent(in) :: filename
    logical, intent(in) :: right
    type(top_points_container), allocatable :: top_points(:)
    integer :: i
    
    top_points = TopPoints(filename, right)
    
    allocate(ConvexPolyhedron :: res%items(size(top_points)))
    
    select type(items => res%items)
      type is (ConvexPolyhedron)
        do i = 1, size(top_points)
          items(i) = ConvexPolyhedron_FromTopPoints(top_points(i)%points)
        end do
    end select
    
  end function Union_Init_TopPoints
  
  
  function Union_Init_Obst(filename) result(res)
    use Strings, only: upcase
    type(Union) :: res    
    character(*),intent(in) :: filename
    integer :: unit,io
    character(180) :: line
    type(ConvexPolyhedron) :: poly
    type(ConvexPolyhedron),allocatable :: items(:)

    open(newunit=unit,file=filename,action='read',status='old',iostat=io)
    
    allocate(items(0))

    if (io==0) then
      do
        read(unit,'(a)',iostat=io) line

        if (io/=0) exit

        line = adjustl(line)

        if (len_trim(line)>0) then

          if (upcase(line(1:10))=='POLYHEDRON') then
          
            call ReadPolyhedron(poly,line(11:))
            
           call add_element(items,poly)
!             items = [items, poly]
            
          end if

        end if

      end do

      close(unit)

    else

      write(*,*) "Could not open file ",filename
      call error_stop

    end if
    
    call move_alloc(items,res%items)
    
    res%size = size(res%items)
    
    contains
    
      subroutine ReadPolyhedron(poly,restline)
        type(ConvexPolyhedron),intent(out) :: poly
        character(*),intent(in)  :: restline
        integer :: nPlanes,i,io

        read(restline,*,iostat=io) nPlanes

        if (io/=0) then
          write(*,*) "Expected number of planes in polyhedron, received '",trim(restline),"' instead."
          call error_stop
        end if

        allocate(poly%Planes(nPlanes))

        poly%nplanes = nPlanes

        do i=1,nPlanes
          call ReadPlane(poly%Planes(i))
        end do

      end subroutine ReadPolyhedron


      subroutine ReadPlane(Pl)
        use Strings
        type(Plane),intent(out) :: Pl
        character(180) :: line
        integer :: io

        read(unit,'(a)',iostat=io) line
        if (io/=0) then
          write(*,*) "Error reading the line with the plane definition."
          call error_stop
        end if

        if (count_multispaces(line) == 4) then
          read(line,*,iostat=io) Pl%a,Pl%b,Pl%c,Pl%d,Pl%gl
        else
          io = 999
        end if
        if (io/=0) then
          write(*,*) "Error parsing the line with the plane definition."
          call error_stop
        end if
      end subroutine ReadPlane
      
      subroutine add_element(a,e)
        type(ConvexPolyhedron),allocatable,intent(inout) :: a(:)
        type(ConvexPolyhedron),intent(in) :: e
        type(ConvexPolyhedron),allocatable :: tmp(:)

        if (.not.allocated(a)) then
          a = [e]
        else
          call move_alloc(a,tmp)
          allocate(a(size(tmp)+1))
          a(1:size(tmp)) = tmp
          a(size(tmp)+1) = e
        end if
      end subroutine
  end function Union_Init_Obst
  
  
  function Union_Init_Geom(filename) result(res)
    use Strings, only: upcase
    type(Union) :: res    
    character(*),intent(in) :: filename
    integer :: unit, io
    integer :: ipoint, ipolyhedron
    integer :: npoints, npolyhedra
    type(ConvexPolyhedron),allocatable :: items(:)
    real(knd),allocatable :: points(:,:)
    character(180) :: line

    open(newunit=unit,file=filename,action='read',status='old',iostat=io)
    
    if (io==0) then
      
      call read_points
      
      call read_polyhedra

      close(unit)

    else

      write(*,*) "Could not open file ",filename
      call error_stop

    end if
    
    call move_alloc(items,res%items)
    
    res%size = size(res%items)
    
    contains
    
      subroutine empty_lines(line)
        character(*), intent(out) :: line

        do
          read (unit, '(a)', iostat=io) line
          if (io/=0) then
            write (*,*) "Error reading from file ", filename
            call error_stop
          end if

          if (len_trim(line)>0) exit
        end do
      end subroutine
    
      subroutine read_points
      
        call read_npoints
        
        allocate(points(3,0:npoints-1))
        
        do ipoint = 0, npoints-1
          call read_point(points(:,ipoint))
        end do
      end subroutine
      
      subroutine read_npoints
        call empty_lines(line)
        call read_npoints_str(line)
        call read_npoints_num(line)
      end subroutine
      
      subroutine read_npoints_str(line)
        character(*), intent(inout) :: line
        line = adjustl(line)
        if (.not.upcase(line(1:7))=='POINTS:') then
          write (*,*) "Error reading points header in file ", filename
          write (*,*) "Expected 'POINTS:', found: ", trim(line)
          call error_stop
        end if
        line = line(8:)
      end subroutine
      
      subroutine read_npoints_num(line)
        character(*), intent(inout) :: line
        
        line = adjustl(line)
        read(line,*,iostat=io) npoints
        if (io/=0) call npoints_error
      end subroutine
      
      subroutine npoints_error
        write (*,*) "Error reading points header in file ", filename
        call error_stop
      end subroutine

      subroutine read_point(p)
        real(knd), intent(out) :: p(3)

        read (unit, '(a)', iostat=io) line
        if (io/=0) call point_error
        
        read (line, *, iostat=io) p
        if (io/=0) call point_error
      end subroutine

      subroutine point_error
        write (*,*) "Error reading point ",ipoint,"in file ", filename
        call error_stop
      end subroutine

      subroutine read_polyhedra
      
        call read_npolyhedra
        
        allocate(items(npolyhedra))
        
        do ipolyhedron = 1, npolyhedra
          call read_polyhedron(items(ipolyhedron))
        end do
      end subroutine
      
      subroutine read_npolyhedra
        do
          read (unit, '(a)', iostat=io) line
          if (io/=0) call npolyhedra_error
          if (len_trim(line)>0) exit
        end do
        call read_npolyhedra_str(line)
        call read_npolyhedra_num(line)
      end subroutine
      
      subroutine read_npolyhedra_str(line)
        character(*), intent(inout) :: line
        integer :: ind

        line = adjustl(line)
        if (.not.upcase(line(1:12))=='POLYHEDRONS:' .and. &
            .not.upcase(line(1:10))=='POLYHEDRA:') then
          write (*,*) "Error reading polyhedra header in file ", filename
          write (*,*) "Expected 'POLYHEDRA:', found: ", trim(line)
          call error_stop
        end if
        ind = index(line(1:12),':')
        line = line(ind+1:)
      end subroutine
      
      subroutine read_npolyhedra_num(line)
        character(*), intent(inout) :: line
        
        line = adjustl(line)
        read(line,*,iostat=io) npolyhedra
        if (io/=0) call npolyhedra_error
      end subroutine
      
      subroutine npolyhedra_error
        write (*,*) "Error reading polyhedra header in file ", filename
        call error_stop
      end subroutine

      subroutine read_polyhedron(poly)
        type(ConvexPolyhedron), intent(out) :: poly
        type(Plane), allocatable :: planes(:)

        call empty_lines(line)
        
        do
          call read_plane(planes, io)
          if (io/=0) exit
        end do
        
        if (.not.allocated(planes)) call polyhedron_error
        
        poly = ConvexPolyhedron(planes)
      end subroutine

      subroutine polyhedron_error
        write (*,*) "Error reading polyhedron ",ipolyhedron,"in file ", filename
        call error_stop
      end subroutine
      
      subroutine read_plane(planes, iostat)
        type(Plane), allocatable, intent(inout) :: planes(:)
        integer, intent(out) :: iostat
        integer :: ind(3)
        
        read (unit, '(a)', iostat=iostat) line
        if (iostat/=0) return
        read(line, *, iostat=iostat) ind
        if (iostat/=0) return
        call add_element( planes, Plane(points(:,ind(1)), &
                                        points(:,ind(2)), &
                                        points(:,ind(3))) ) 
      end subroutine

      subroutine add_element(a,e)
        type(Plane),allocatable,intent(inout) :: a(:)
        type(Plane),intent(in) :: e
        type(Plane),allocatable :: tmp(:)

        if (.not.allocated(a)) then
          a = [e]
        else
          call move_alloc(a,tmp)
          allocate(a(size(tmp)+1))
          a(1:size(tmp)) = tmp
          a(size(tmp)+1) = e
        end if
      end subroutine
  end function Union_Init_Geom
  
  
  
  subroutine ArrayFromObst(res,filename)
    use Strings, only: upcase
    class(GeometricShape),allocatable,intent(out) :: res(:)
    character(*),intent(in) :: filename
    integer :: unit,io
    character(180) :: line
    type(ConvexPolyhedron) :: poly
    type(ConvexPolyhedron),allocatable :: items(:)

    open(newunit=unit,file=filename,action='read',status='old',iostat=io)
    
    allocate(items(0))

    if (io==0) then
      do
        read(unit,'(a)',iostat=io) line

        if (io/=0) exit

        line = adjustl(line)

        if (len_trim(line)>0) then

          if (upcase(line(1:10))=='POLYHEDRON') then
          
            call ReadPolyhedron(poly,line(11:))
            
            call add_element(items,poly)

           end if
        end if

      end do

      close(unit)

    else

      write(*,*) "Could not open file ",filename
      call error_stop

    end if
    
    call move_alloc(items,res)
    
    contains
    
      subroutine ReadPolyhedron(poly,restline)
        type(ConvexPolyhedron),intent(out) :: poly
        character(*),intent(in)  :: restline
        integer :: nPlanes,i,io

        read(restline,*,iostat=io) nPlanes

        if (io/=0) then
          write(*,*) "Expected number of planes in polyhedron, received '",trim(restline),"' instead."
          call error_stop
        end if

        allocate(poly%Planes(nPlanes))

        poly%nplanes = nPlanes

        do i=1,nPlanes
          call ReadPlane(poly%Planes(i))
        end do

      end subroutine ReadPolyhedron


      subroutine ReadPlane(Pl)
        use Strings
        type(Plane),intent(out) :: Pl
        character(180) :: line
        integer :: io

        read(unit,'(a)',iostat=io) line
        if (io/=0) then
          write(*,*) "Error reading the line with the plane definition."
          call error_stop
        end if

        if (count_multispaces(line) == 4) then
          read(line,*,iostat=io) Pl%a,Pl%b,Pl%c,Pl%d,Pl%gl
        else
          io = 999
        end if
        if (io/=0) then
          write(*,*) "Error parsing the line with the plane definition."
          call error_stop
        end if
      end subroutine ReadPlane
      
      subroutine add_element(a,e)
        type(ConvexPolyhedron),allocatable,intent(inout) :: a(:)
        type(ConvexPolyhedron),intent(in) :: e
        type(ConvexPolyhedron),allocatable :: tmp(:)

        if (.not.allocated(a)) then
          a = [e]
        else
          call move_alloc(a,tmp)
          allocate(a(size(tmp)+1))
          a(1:size(tmp)) = tmp
          a(size(tmp)+1) = e
        end if
      end subroutine
  end subroutine
  
  
  function TopPoints(filename, right) result(res)
    type(top_points_container), allocatable :: res(:)
    character(*),intent(in) :: filename
    logical, intent(in) :: right
    logical :: ex, body_opened
    real(knd), allocatable :: points(:,:)
    integer :: u, stat

    allocate(res(0))
    
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
        res = [res, top_points_container(points)]
      else
        res = [res, top_points_container(points(:,size(points,2):1:-1))]
      end if

      deallocate(points)
      body_opened = .false.
    end subroutine

  end function TopPoints
  
  
  function GeometricShape_FromFile(filename) result(res)
    use Strings
    class(GeometricShape), allocatable :: res
    character(*), intent(in) :: filename
    character(5) :: suffix
    logical :: ex

    inquire(file=filename,exist=ex)

    if (.not.ex) &
      call error_stop("Error, obstacle file '"//trim(filename)//"' does not exist.")
    
    suffix = filename(index(filename,'.',back=.true.):)

    if (suffix=='.obst'.or.suffix=='.geom') then
      allocate(res, source = Union(filename))
    else if (suffix=='.off') then
      allocate(res, source = Polyhedron(filename))
    else if (suffix=='.ltop'.or.suffix=='.rtop') then
      allocate(res, source = Union(filename))
    else
      write(*,*) "Unknown file format for geometric shape '"//suffix//"'."
      call error_stop
    end if
  end function

end module GeometricShapes







module Body_class
  use Kinds, only: knd
  use GeometricShapes, only: GeometricShape, Ray
  use Parameters
  
  implicit none

  private

  public Body,Inside


  type,abstract :: Body
     integer :: numofbody
     class(GeometricShape),allocatable :: GeometricShape
  contains
     procedure :: Inside => CInside  !Hack around yet inidentified problem in GCC.
     procedure :: Closest
     procedure :: ClosestOut
     procedure :: ClosestOnLineOut
     procedure :: OutwardNormal
     procedure :: IntersectsRay
  end type
  interface Inside
    module procedure CInside
  end interface

contains


  logical function CInside(self,x,y,z,eps)
    class(Body),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in),optional :: eps
    real(knd) :: x2,y2,z2
    real(knd) :: lx,ly,lz

    if (.not.allocated(self%GeometricShape)) then
      CInside = .false.
    else
      x2 = x
      y2 = y
      z2 = z

      lx = xU(Prnx) - xU(0)
      ly = yV(Prny) - yV(0)
      lz = zW(Prnz) - zW(0)

      if (Btype(Ea)==BC_PERIODIC.and.x2>xU(Prnx+1)) x2 = x2-lx
      if (Btype(No)==BC_PERIODIC.and.y2>yV(Prny+1)) y2 = y2-ly
      if (Btype(To)==BC_PERIODIC.and.z2>zW(Prnz+1)) z2 = z2-lz

      if (Btype(We)==BC_PERIODIC.and.x2<xU(0)) x2 = x2+lx
      if (Btype(So)==BC_PERIODIC.and.y2<yV(0)) y2 = y2+ly
      if (Btype(Bo)==BC_PERIODIC.and.z2<zW(0)) z2 = z2+lz
      CInside = self%GeometricShape%Inside(x2,y2,z2,eps)

    end if

  end function CInside

  

  subroutine Closest(self,xnear,ynear,znear,x,y,z)
    class(Body),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%GeometricShape%Closest(xnear,ynear,znear,x,y,z)

  end subroutine Closest

  

  subroutine ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(Body),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%GeometricShape%ClosestOut(xnear,ynear,znear,x,y,z)

  end subroutine ClosestOut



  real(knd) function ClosestOnLineOut(self,x,y,z,x2,y2,z2) result(t)
    !Find t, such that x+(x2-x)*t lies on the boundary of the SB
    ! x lies inside SB
    class(Body),intent(in) :: self
    real(knd),intent(in) :: x,y,z,x2,y2,z2
    real(knd) :: t1,t2
    integer :: i
    

    t1 = 0
    t2 = 1
    !First, find a point lying outside. We should have the right direction.
    if (self%Inside(x2,y2,z2,0._knd)) then
     t1 = 1
     do
      t2 = t2*1.1_knd
      if (.not. self%Inside(x+(x2-x)*t2,y+(y2-y)*t2,z+(z2-z)*t2,0._knd)) exit
      if (t2>=2) then
        t = 10
        return
      end if
     enddo
    endif
    
    t = (t1+t2)/2._knd

    do i = 1,20         !The bisection method with maximum 20 iterations (should be well enough)
     if (self%Inside(x+(x2-x)*t,y+(y2-y)*t,z+(z2-z)*t,0._knd)) then
      t1 = t
     else
      t2 = t
     endif
     t = (t1+t2)/2._knd

     if (abs(t1-t2)<MIN(dxmin/1000._knd,dymin/1000._knd,dzmin/1000._knd))   exit
    enddo

  end function ClosestOnLineOut
  

  function OutwardNormal(self,x,y,z) result(res)
    real(knd) :: res(3)
    class(Body),intent(in) :: self
    real(knd),intent(in) :: x,y,z

    res = self%GeometricShape%OutwardNormal(x,y,z)
 
  end function

  logical function IntersectsRay(self,r) result(intersects)
    class(Body),intent(in) :: self
    class(Ray),intent(in) :: r

    intersects = self%GeometricShape%IntersectsRay(r)
 
  end function

end module Body_class






