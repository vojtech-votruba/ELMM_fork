module WaterThermodynamics
  use Parameters
  use PhysicalProperties
  
  implicit none
  
  ! array Temperature becomes the linearized liquid water potential temperature (Betts 1973)
  ! theta_l = theta - (Lv_water/Cp_air * theta/temperature) * q_l

  ! a diagnostic quantity, does not need to be treated the same way as Temperature and Moisture
  real(knd), allocatable, dimension(:,:,:)   :: LiquidWater    !Total specific humidity q_t.

contains

  subroutine compute_liquid_water_content(Temperature, Moisture, Pr)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Temperature, Moisture
    real(knd), contiguous, intent(in) :: Pr(-1:,-1:,-1:)
    real(knd) :: p
    integer :: i, j, k
    
    do k = 1, Wnz+1
      do j = 1, Prny
        do i = 1, Prnx
          p = reference_pressure_z(k) + Pr(i,j,k) * rho_air_ref
          LiquidWater(i,j,k) = &
            liquid_water_condensed(Moisture(i,j,k), Temperature(i,j,k), p)
          ! now the thermodynamic temperature can be computed if needed
        end do
      end do
    end do
  end subroutine

  ! diagnostic function to compute the liquid water contant from the 
  ! total water content q_t, liquid temperature q_l and pressure p
  function liquid_water_condensed(q_t, theta_l, p) result(q_l)
    real(knd) :: q_l
    real(knd), intent(in) :: q_t, theta_l, p
    real(knd) :: q_s
    
    q_s = saturation_humidity(theta_l, p, q_t)
    
    q_l = max(q_t - q_s, 0._knd)
  end function

  function saturation_humidity(theta_l, p, q_t) result(q_s)
    real(knd) :: q_s ! kg/kg
    real(knd), intent(in) :: theta_l, p, q_t ! K, Pa, kg/kg
    real(knd) :: T_l, e_s, q_sl, beta_l
    real(knd), parameter :: reference_pressure = 100000
    
    ! liquid temperature
    T_l = (p / reference_pressure)**(Rd_air_ref / Cp_air_ref) * theta_l
    
    ! water vapour saturation pressure
    e_s = water_v_saturation_p(T_l)
    
    ! saturation specific humidity at T_l
    q_sl = 0.622 * e_s / (p - 0.378 * e_s)
    
    ! saturation specific humidity at temperature T (Clausius-Clapeiron)
    beta_l = 0.622 * (Lv_water_ref / (Rd_air_ref * T_l)) * (Lv_water_ref / (Cp_air_ref * T_l))
    q_s = q_sl * (1 + beta_l * q_t) / (1 + beta_l *q_sl)      
    
  end function


  function water_v_saturation_p(T_l) result(e_s)
    real(knd) :: e_s ! in Pascals!
    real(knd), intent(in) :: T_l
    real(knd), parameter :: e_s0 = 610.78 ! Pa
    real(knd), parameter :: T_trip = 273.16 ! K
    real(knd), parameter :: a = 17.27
    real(knd), parameter :: b = 35.86 ! K
    
    e_s = e_s0 * exp(a * (T_l - T_trip) / (T_l - b))
  end function
  
  
  
end module WaterThermodynamics

