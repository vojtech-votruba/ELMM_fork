module FreeUnit
  implicit none

  contains

  !newunit()
  !Finds a not opened logical unit.
  !Not needed in Fortran 2008.

  subroutine newunit(unit)
     integer,intent(out) :: unit
     logical :: isOpen

     integer, parameter :: MIN_UNIT_NUMBER = 10
     integer, parameter :: MAX_UNIT_NUMBER = 99

     do unit = MIN_UNIT_NUMBER, MAX_UNIT_NUMBER
        inquire(unit = unit, opened = isOpen)
        if (.not. isOpen) then
           return
        end if
     end do
     
     unit = -1
  end subroutine newunit

end module FreeUnit
