module ieee_arithmetic
  !poor man's replacement for the intrinsic module
  !do not use if the compiler supports the F2003 version
  !make sure not to use fast math optimizations
  use iso_fortran_env
  
  !common extension isnan may actually fail with optimizations above
!   intrinsic isnan
  
  interface ieee_is_nan
    module procedure ieee_is_nan_real32
    module procedure ieee_is_nan_real64
  end interface

contains
  logical elemental function ieee_is_nan_real32(x) result(res)
    real(real32), intent(in) :: x

    res = x /= x
  end function

  logical elemental function ieee_is_nan_real64(x) result(res)
    real(real64), intent(in) :: x

    res = x /= x
  end function

end module
