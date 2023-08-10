function CustomSurfaceTemperature(x,y,z,t) result(res)
   use Kinds
   real(knd) :: res
   real(knd), intent(in) :: x, y, z
   real(tim), intent(in) :: t

   res =  265 - t / (4 * 3600) ! cooling 0.25 K / hour

end function

subroutine  CustomSolidBodies
end subroutine  CustomSolidBodies
