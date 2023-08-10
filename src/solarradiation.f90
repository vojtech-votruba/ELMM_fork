module SolarRadiation
  use Parameters
  use Pressure, only: pressure_solution
  
  implicit none
  
  real(knd),parameter :: SB_sigma = 5.6704E-8_knd

  real(knd) :: sun_azimuth, sun_elevation !in degrees
  real(knd) :: vector_to_sun(3) !unit vector pointing to sun in grid coordinates
  
  integer :: svf_nrays = 50
  
  real(knd),allocatable :: svf_vecs(:,:)
  
contains
  
  pure function solar_direct_flux() result(res)
    real(knd) :: res
    res = 800._knd
  end function

  pure function solar_diffuse_flux() result(res)
    real(knd) :: res
    res = 200._knd
  end function

  subroutine InitSolarRadiation
    real(knd) :: horiz_component, horiz_angle
    real(knd) :: azimuth, z
    integer :: i
      
    sun_azimuth = 237
    sun_elevation = 53
    
    vector_to_sun(3) = sin(sun_elevation * pi / 180)
    horiz_component = cos(sun_elevation * pi / 180)
    horiz_angle = (x_axis_azimuth - sun_azimuth) * pi / 180
    
    vector_to_sun(1) = horiz_component * cos(horiz_angle)
    vector_to_sun(2) = horiz_component * sin(horiz_angle)
    
    allocate(svf_vecs(3,svf_nrays))
    
    do i=1,svf_nrays
      !second algorithm at http://mathworld.wolfram.com/SpherePointPicking.html
      call random_number(azimuth)
      call random_number(z)

      azimuth = azimuth * pi * 2
      horiz_component = sqrt(1 - z**2)

      svf_vecs(1,i) = horiz_component * cos(azimuth)
      svf_vecs(2,i) = horiz_component * sin(azimuth)
      svf_vecs(3,i) = z
    end do
    
  end subroutine
  
 
  
  
  pure function out_lw_radiation(emissivity, T) result(res)
    real(knd) :: res
    real(knd),intent(in) :: emissivity, T
    
    res = SB_sigma * emissivity * T**4
  end function
  
  pure function in_lw_radiation() result(res)
    real(knd) :: res
    !http://www.the-cryosphere-discuss.net/2/487/2008/tcd-2-487-2008.pdf
    res = clear_sky_lw_radiation() * cloud_factor_lw()
  end function
  
  pure function clear_sky_lw_radiation() result(res)
    real(knd) :: res
    real(knd) :: wv_press, atmosphere_emissivity
    
    wv_press = moisture_ref / (moisture_ref + 0.622) * &
                               pressure_solution%bottom_pressure
    
    atmosphere_emissivity = 0.23 + &
               0.44 * (wv_press / temperature_ref)**(0.125_knd)
    
    res = SB_sigma * atmosphere_emissivity * temperature_ref**4
  end function
  
  pure function cloud_factor_lw() result(res)
    real(knd) :: res
    real(knd) :: cloud_fraction
    cloud_fraction = 0
    res = 1 + 0.22 * cloud_fraction**2
  end function
  
end module SolarRadiation
