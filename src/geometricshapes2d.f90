module GeometricShapes2D
  use Parameters

  implicit none

  public :: GeometricShape2D, Circle

  private

  type,abstract :: GeometricShape2D
  contains
    procedure(GeometricShape2D_Inside), deferred :: Inside 
  end type

  type, extends(GeometricShape2D) :: Circle
    real(knd) :: xc, yc, r
  contains
    procedure :: Inside => Circle_Inside
  end type

  abstract interface
    logical function GeometricShape2D_Inside(self,x,y)
      import
      class(GeometricShape2D), intent(in) :: self
      real(knd), intent(in) :: x,y
    end function
  end interface


contains

  logical function Circle_Inside(self,x,y) result(ins)
    class(Circle), intent(in) :: self
    real(knd), intent(in) :: x,y

    ins = (self%xc-x)**2 + (self%yc-y)**2 <= (self%r)**2
  end function
  
end module GeometricShapes2D

