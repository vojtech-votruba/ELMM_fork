module Endianness

  use iso_fortran_env
  
  implicit none

  private
  public :: GetEndianness, BigEnd, SwapB, littleendian

  logical, save :: littleendian = .true.

  interface BigEnd
    module procedure BigEnd16
    module procedure BigEnd32
    module procedure BigEnd32_a1
    module procedure BigEnd32_a2
    module procedure BigEnd32_a3
    module procedure BigEnd32_a4
    module procedure BigEnd64
    module procedure BigEnd64_a1
    module procedure BigEnd64_a2
    module procedure BigEnd64_a3
    module procedure BigEnd64_a4
  end interface

  interface SwapB
    module procedure SwapB32
    module procedure SwapB64
  end interface

  contains

    subroutine GetEndianness
      character(4) :: bytes !may not work on some processors

      bytes = transfer(1_int32,bytes)
      if (ichar(bytes(4:4))==1) then
        littleendian=.false.
      else
        littleendian=.true.
      endif
    end subroutine GetEndianness

   
    elemental function BigEnd16(x) result(res)
      integer(int16) :: res
      integer(int16),intent(in)::x
      character(2) :: bytes
      
      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes)
        res = ichar(bytes(2:2),int16)
        res = ior( ishft(ichar(bytes(1:1),int16),8), res )
      endif
    end function
    
    function BigEnd32(x) result(res)
      real(real32),intent(in) :: x
      real(real32) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a1(x) result(res)
      real(real32),intent(in) :: x(:)
      real(real32) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a2(x) result(res)
      real(real32),intent(in) :: x(:,:)
      real(real32) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a3(x) result(res)
      real(real32),intent(in) :: x(:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a4(x) result(res)
      real(real32),intent(in) :: x(:,:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64(x) result(res)
      real(real64),intent(in) :: x
      real(real64) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a1(x) result(res)
      real(real64),intent(in) :: x(:)
      real(real64) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a2(x) result(res)
      real(real64),intent(in) :: x(:,:)
      real(real64) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a3(x) result(res)
      real(real64),intent(in) :: x(:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a4(x) result(res)
      real(real64),intent(in) :: x(:,:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    elemental function SwapB32(x) result(res)
      real(real32) :: res
      real(real32),intent(in) :: x
      character(4) :: bytes
      integer(int32) :: t
      real(real32) :: rbytes, rt
      !nicer looking TRANSFER is problematic to optimize by ifort
      ! gfortran would be OK with that
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(4:4),int32)

      t = ior( ishftc(ichar(bytes(3:3),int32),8),  t )

      t = ior( ishftc(ichar(bytes(2:2),int32),16), t )

      t = ior( ishftc(ichar(bytes(1:1),int32),24), t )

      res = rt
        
    end function
    
    elemental function SwapB64(x) result(res)
      real(real64) :: res
      real(real64),intent(in) :: x
      character(8) :: bytes
      integer(int64) :: t
      real(real64) :: rbytes, rt
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(8:8),int64)

      t = ior( ishftc(ichar(bytes(7:7),int64),8),  t )

      t = ior( ishftc(ichar(bytes(6:6),int64),16), t )

      t = ior( ishftc(ichar(bytes(5:5),int64),24), t )

      t = ior( ishftc(ichar(bytes(4:4),int64),32), t )

      t = ior( ishftc(ichar(bytes(3:3),int64),40), t )

      t = ior( ishftc(ichar(bytes(2:2),int64),48), t )

      t = ior( ishftc(ichar(bytes(1:1),int64),56), t )

      res = rt

    end function

end module Endianness

