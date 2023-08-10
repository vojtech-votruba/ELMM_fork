subroutine  CustomSolidBodies
   use Parameters
   use r3_type
   use GeometricShapes
   use SolidBodies
   use Body_class

   implicit none

   type(ConvexPolyhedron) :: poly
   integer :: i, j
   real(knd) :: base(3,4), points(3,4), period(2)

   obstacles_bbox(Bo) = 0
   obstacles_bbox(To) = 1

   base(:,1) = [ 2.5_knd,-0.5_knd, 1.0_knd]
   base(:,2) = [ 2.5_knd, 0.5_knd, 1.0_knd]
   base(:,3) = [ 0.5_knd, 0.5_knd, 1.0_knd]
   base(:,4) = [ 0.5_knd,-0.5_knd, 1.0_knd]

   period = [3.0_knd, 2.0_knd]
   
   points = base
   
   do j = -2, 3
     do i = -2, 5
       points(1,:) = base(1,:) + i * period(1)
       points(2,:) = base(2,:) + j * period(2)
       call AddSolidBody(SolidBody(ConvexPolyhedron_FromTopPoints(points)))
     end do
   end do
   
   
   
   
end subroutine  CustomSolidBodies
