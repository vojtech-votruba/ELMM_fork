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
   obstacles_bbox(To) = 3

   base(:,4) = [-0.5_knd, 2.5_knd, 1.0_knd]
   base(:,3) = [ 0.5_knd, 2.5_knd, 1.0_knd]
   base(:,2) = [ 0.5_knd, 0.5_knd, 1.0_knd]
   base(:,1) = [-0.5_knd, 0.5_knd, 1.0_knd]

   period = [2.0_knd, 3.0_knd]
   
   points = base
   
!    poly = ConvexPolyhedron_FromTopPoints(points)
!    print *, poly%Inside(0._knd, 1._knd, 0.5_knd)
!    poly = ConvexPolyhedron_FromTopPoints(points(:,4:1:-1))
!    print *, poly%Inside(0._knd, 1._knd, 0.5_knd)
   
   do j = -3, 2
     do i = -3, 12
       points(1,:) = base(1,:) + i * period(1)
       points(2,:) = base(2,:) + j * period(2)

       if (i==0.and.j==-1) then
         points(3,:) = base(3,:) * 3
       else
         points(3,:) = base(3,:)
       end if

       call AddSolidBody(SolidBody(ConvexPolyhedron_FromTopPoints(points)))
     end do
   end do
   
   
   
   
end subroutine  CustomSolidBodies
