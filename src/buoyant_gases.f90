module BuoyantGases
  use Parameters
  use PhysicalProperties
  
  implicit none
  
  ! At the moment buoyant gases should not be combined with air moisture.

  ! Stores information necessary for computation of density
  ! of mixture of gases based on their concantration stored in Scalar.
  ! The concentrations shall be mass per volume kg . m^-3.
  
  ! The conversion to mass fraction w = m_g / (m_g + m_d) shall
  ! assume pressure and temperature that would give the reference value of density rho_air_ref.

  ! The array does not have to cover all scalars simulated. The simulation shall safely proceed
  ! assuming the other scalars are neutral.
  
  ! The array may cover more scalars than actually simulated.
  
  logical :: enable_buoyant_scalars = .false.
  
  type buoyant_scalar
    real(knd) :: m_mol = m_mol_air
    real(knd) :: R = Rd_air_ref  ! R_gas_universal / m_mol
    real(knd) :: eps = 0         ! m_mol_air_ref/m_mol - 1
  end type
  
  ! Index of the array correspends to the number of the scalar.
  type(buoyant_scalar), allocatable :: buoyant_scalars(:)
  
  interface buoyant_scalar
    module procedure buoyant_scalar_init
  end interface
contains

   
  elemental function buoyant_scalar_init(m_mol) result(res)
    type(buoyant_scalar) :: res
    real(knd), intent(in) :: m_mol
    
    res%m_mol = m_mol
    res%R = R_gas_universal / m_mol
    res%eps = m_mol_air / m_mol - 1
  end function
  
  elemental function mass_fraction(sc, rho_g)
    real(knd) mass_fraction
    type(buoyant_scalar), intent(in) :: sc
    real(knd), intent(in) :: rho_g
    ! rho_g is the mass density of the scalar in kg / m^3
    ! mass fraction is rho_g / (rho_g + rho_d)
    
    real(knd) :: p_r, p_g, rho_d
    
    !TODO: avoid re-computing this constant value
    p_r = rho_air_ref * Rd_air_ref * temperature_ref
    
    p_g = rho_g * sc%R * temperature_ref
    
    rho_d = (p_r - p_g) / (Rd_air_ref * temperature_ref)
    
    mass_fraction = rho_g / (rho_g + rho_d)
    
    mass_fraction = min(1._knd, max(0._knd, mass_fraction))
  end function
   
end module

