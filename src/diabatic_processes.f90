module DiabaticProcesses
  use Parameters

  implicit none

  logical :: enable_diabatic_processes = .false.

contains

  subroutine diabatic_heat_terms(Temperature)
    use Radiation
    real(knd), contiguous, intent(inout) :: Temperature(-2:,-2:,-2:)

    !individual types of diabatic heating terms are called from here
    if (enable_radiation_profile) then
      call apply_radiation_heat_profile(Temperature)
    end if
  end subroutine

end module DiabaticProcesses
