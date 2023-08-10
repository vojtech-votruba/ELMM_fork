module Kinds

  use iso_fortran_env
  use iso_c_binding, only: c_int, c_float, c_double, c_size_t

  implicit none

!for reference when iso_fortran_env is not available:
!   integer, parameter :: int8   = selected_int_kind(1)
!   integer, parameter :: int32  = selected_int_kind(9)
!   integer, parameter :: int64  = selected_int_kind(10)
!   integer, parameter :: real32 = selected_real_kind(p = 6, r = 37)
!   integer, parameter :: real64 = selected_real_kind(p = 15, r = 200)

  integer, parameter :: dbl = real64, sng = real32

#ifdef DPREC
  integer, parameter :: knd = DBL                                       !knd is default real kind for the whole program, choosing double
#else
  integer, parameter :: knd = SNG                                       !knd is default real kind for the whole program, choosing single
#endif

  integer, parameter :: tim = knd                                       !Kind for time variables, can be double for very small timesteps.
                                                                       !It may affect performance

  integer, parameter :: sint = kind(1)                                  !To save memory a smaller type can be used for some integer
  integer, parameter :: slog = sint                                     ! and logical arrays. Note the same KIND value is guaranteed
                                                                       ! the default intrinsic types.
                                                                       !This can have some negative effect on speed however
end module Kinds



module Parameters

  use Kinds
  use Stop_procedures

  implicit none

  save

  real(knd), parameter :: pi = acos(-1.0_knd)
  real(knd), parameter :: Karman = 0.41_knd
  real(knd), parameter :: BoltzC = 1.3806503e-23_knd

  integer :: Unx, Uny, Unz     !dimensions of grid for velocity component U

  integer :: Vnx, Vny, Vnz     !dimensions of grid for velocity component V

  integer :: Wnx, Wny, Wnz     !dimensions of grid for velocity component W

  integer :: Prnx, Prny, Prnz  !dimensions of grid for pressure
  
  !relevant for distributed parallelism
  integer :: gPrnx, gPrny, gPrnz !global grid dimensions
  integer :: gUnx, gUny, gUnz
  integer :: gVnx, gVny, gVnz
  integer :: gWnx, gWny, gWnz
  
  integer :: offset_to_global_x = 0, offset_to_global_y = 0, offset_to_global_z = 0
  integer :: gPrns(3), offsets_to_global(3) = 0


  !the domain's (global) extents
  real(knd) :: gxmin, gxmax, gymin, gymax, gzmin, gzmax

  !the image's extents
  real(knd) :: im_xmin, im_xmax, im_ymin, im_ymax, im_zmin, im_zmax

  real(knd), allocatable :: xU(:), xPr(:), yV(:), yPr(:), zW(:), zPr(:)          !coordinates of grid points
  real(knd), allocatable :: dxU(:), dxPr(:), dyV(:), dyPr(:), dzW(:), dzPr(:)    !dxPr(i)=xU(i)-xU(i-1), dxU(i)=xPr(i+1)-xPr(i)

  real(knd), allocatable :: gxU(:), gxPr(:), gyV(:), gyPr(:), gzW(:), gzPr(:)          !coordinates of grid points






  
  type time_step_control
    integer   :: max_number_of_time_steps = 2**30 !maximum number of time steps
    
    integer :: check_period = 100 !how often check for a change in the input file 
                                  !or if the time limit has been exceeded
                                  
    real(knd) :: clock_time_limit = huge(1._tim)

    logical   :: variable_time_steps = .true.

    logical   :: enable_U_scaling = .false. !enables that the constant time step to be computed from 
                                           ! U_time_scaling and from CFL
    logical   :: enable_CFL_check = .false. !enables check of maximum CFL 
                                           ! U_time_scaling and from CFL
    real(tim) :: dt_constant = 0

    real(tim) :: dt

    real(tim) :: dt_max           !minimal time step for diagnosing diverging simulation
    
    real(tim) :: dt_min           !maximum time step, larger value will not be used
    
    real(knd) :: U_scaling(3) = 0 !velocity vector used to compute constant_U

    real(knd) :: U_max(3) = 0     !velocity used to compute minimal_dt

    real(knd) :: U_min(3) = 0     !velocity used to compute maximal_dt

    real(knd) :: CFL = 0.3        !CFL used to compute dt, or computed dt if CFL is fixed

    real(knd) :: CFL_max = 0      !CFL which will result in diagnosing the simulation as diverging

    real(tim) :: effective_time   !the time including the partial time-steps during Runge-Kutta stages

    real(tim) :: time             !time of the start of the time step

    real(tim) :: start_time, end_time

  end type

  type(time_step_control) :: time_stepping


  real(knd) :: dxmin, dymin, dzmin  !minimum grid spacing, dimensions of the domain


  real(knd) :: molecular_viscosity = 1._knd / 70000
  real(knd) :: molecular_diffusivity = 1._knd / 70000 / 0.7 !default Prandtl number 0.7


  real(knd) :: pr_gradient_x = 0, pr_gradient_y = 0, pr_gradient_z = 0
  real(knd), allocatable :: pr_gradient_profile_x(:), pr_gradient_profile_y(:)

  logical :: enable_pr_gradient_x_uniform = .false.
  logical :: enable_pr_gradient_y_uniform = .false.
  logical :: enable_pr_gradient_z_uniform = .false.
  
  logical :: enable_pr_gradient_x_profile = .false.
  logical :: enable_pr_gradient_y_profile = .false.


  real(knd) :: grav_acc = 9.81, Coriolis_parameter = 0

  real(knd) :: ShearInletTypeParameter, Uinlet, Uinlet_vec(3)

  real(knd) :: z0W, z0E, z0S, z0N, z0B, z0T

  real(knd) :: windangle = 0

  real(knd) :: totalscalsource

  integer :: scalsourcetype

  real(knd) :: epsCN, epsPoisson, eps, debugparam

  real(tim) :: timefram1, timefram2, timeavg1, timeavg2
  
  integer :: discretization_order = 2

  integer :: advection_method
  integer :: frames
  integer :: steady
  !task type can enable some hard-coded special cases, generally should be left 0
  integer :: task_type = 0
  integer :: averaging

  logical :: enable_ibm_mass_sources = .true.

  logical :: explicit_diffusion = .true.
  logical :: explicit_scalar_diffusion = .true.

  logical :: enable_buoyancy = .false. !1 if enabled, zero otherwise
  logical :: enable_moisture = .false.
  logical :: enable_liquid   = .false. !enable condensation of water vapor
  integer :: num_of_scalars  = 0

  logical :: enable_radiation = .false.

  integer :: partdistrib, computedeposition, computegravsettling
  integer :: maxCNiter, maxPOISSONiter, endstep

  integer :: inlettype, gridtype, profiletype


  real(knd), allocatable :: Uin(:,:), Vin(:,:), Win(:,:), Uoutb(:,:)

  real(knd), allocatable :: TempIn(:,:), MoistIn(:,:)
  
  real(knd), allocatable :: reference_pressure_z(:)







  real(knd) :: x_axis_azimuth = 90!true geographic heading of the x axis in degrees

  real(knd), allocatable :: Viscosity(:,:,:), TDiff(:,:,:)  !(turbulent) viscosity, (turbulent) thermal diffusivity

  integer(sint), allocatable, dimension(:,:,:) :: Utype, Vtype, Wtype, Prtype !number of solid body inside which the point is or 0
  integer, allocatable, dimension(:,:) :: Unull, Vnull, Wnull                !indexes of points to be nulled every timestep

  integer   :: nUnull, nVnull, nWnull  !second dimension of arrays above (number of points)
  
  integer :: n_free_im, n_free_im_U, n_free_im_V, n_free_im_W
  integer :: n_free_domain, n_free_domain_U, n_free_domain_V, n_free_domain_W
  integer :: n_full_im, n_full_im_U, n_full_im_V, n_full_im_W
  integer :: n_full_domain, n_full_domain_U, n_full_domain_V, n_full_domain_W

  logical   :: xgridfromfile, ygridfromfile, zgridfromfile
  integer   :: initcondsfromfile

  integer, parameter :: We = 1, Ea = 2, So = 3, No = 4, Bo = 5, To = 6

  real(knd) :: temperature_ref = 295
  real(knd) :: moisture_ref = 0.001 !TODO: compute from relative humidity

  integer, dimension(6) :: Btype      !boundary condition types for velocity, see below for values
  integer, dimension(6) :: TempBtype  !boundary condition types for temperature
  integer, dimension(6) :: MoistBtype !boundary condition types for temperature
  integer, dimension(6) :: ScalBtype  !boundary condition types for scalars
  
  integer, dimension(6) :: PoissonBtype !boundary conditions for the pressure solver


  real(knd), dimension(3,6)   :: sideU     !velocities on boundaries in case of Dirichlet BC
  real(knd), dimension(6)     :: sideTemp  !temperatures or temperature fluxes on boundaries
  real(knd), dimension(6)     :: sideMoist !moistures or moisture fluxes on boundaries
  real(knd), dimension(6)     :: sideScal  !scalars or scalar fluxes on boundaries

  real(knd), allocatable      :: BsideTArr(:,:), BsideTFLArr(:,:)
  real(knd), allocatable      :: BsideMArr(:,:), BsideMFLArr(:,:)


  integer, parameter :: ScalarTypeTemperature = 1, &
                        ScalarTypeMoisture = 2, &
                        ScalarTypePassive = 3

  integer, parameter :: BC_NOSLIP = 1, BC_FREESLIP = 2, BC_PERIODIC = 3, BC_DIRICHLET = 4, BC_NEUMANN = 5, BC_CONSTFLUX = 6, &  !boundary condition types
                        BC_TURBULENT_INLET = 7, BC_INLET_FROM_FILE = 8, &
                        BC_AUTOMATIC_FLUX = 9, &
                        BC_MPI_BOUNDS_MIN = 1000, BC_MPI_BOUNDS_MAX = 1010, BC_MPI_BOUNDARY = 1000, BC_MPI_PERIODIC = 1001, &
                        BC_DOMAIN_BOUNDS_MIN = 2000, BC_DOMAIN_BOUNDS_MAX = 2010, &
                        BC_DOMAIN_NESTED = 2001, BC_DOMAIN_COPY = 2002, BC_DOMAIN_NESTED_TURBULENT = 2003

  !set by user
  logical :: enable_fixed_flow_rate = .false.
  !set from boundary conditions automatically
  logical :: flow_rate_x_fixed = .false., flow_rate_y_fixed = .false., flow_rate_z_fixed = .false.
  real(knd) :: flow_rate_x, flow_rate_y, flow_rate_z

  !inlet types
  integer, parameter :: ZeroInletType = 0, ConstantInletType = 1, ShearInletType = 2, &
                        ParabolicInletType = 3, TurbulentInletType = 4, FromFileInletType = 5, &
                        GeostrophicInletType = 6
  !inlet profile types
  integer, parameter :: PROFILE_CONSTANT = 1, PROFILE_LOGARITHMIC = 2, PROFILE_POWER_LAW = 3, PROFILE_FROM_FILE = 4
  integer, parameter :: GRID_UNIFORM = 0, GRID_VARIABLE_Z = 1, GRID_GENERAL = 2

  !scalar source types
  integer, parameter :: PointSource = 1, LineSource = 2, AreaSource = 3, VolumeSource = 4

  integer(c_int), bind(C, name="debuglevel") :: debuglevel = 0 !amount of information to write out
  
  logical :: init_phase = .true., run_phase = .false.
  
  character(:), allocatable :: output_dir, shared_output_dir, input_dir
  character(2048) :: scratch_dir = ""

  integer(dbl) :: timer_rate
  real(dbl) :: time_steps_time = 0
  real(dbl) :: poisson_solver_time = 0

end module Parameters


module PhysicalProperties
  use Parameters

  implicit none
  
  real(knd), parameter :: R_gas_universal = 8.31446261815324_knd !J.kg^-1.mol^-1

  real(knd), parameter :: m_mol_air = 28.9644e-3_knd !kg.mol^-1  NIST 400 ppm CO2
  
  real(knd), parameter :: Rd_air_ref = R_gas_universal / m_mol_air !J.kg^-1.K^-1

  real(knd), parameter :: rho_air_ref = 1.196_knd !kg.m^-3

  real(knd), parameter :: Cp_air_ref = 1005 !J.kg^-1.K^-1

  real(knd), parameter :: Cv_air_ref = 718 !J.kg^-1.K^-1

  real(knd), parameter :: Lv_water_ref = 2442000 !J.kg^-1

end module


module RK3
  use Kinds

  implicit none

  real(knd), parameter :: RK_alpha(3) = [ 4._knd/15._knd, 1._knd/15._knd,  1._knd/6._knd  ]
  real(knd), parameter :: RK_beta(3)  = [ 8._knd/15._knd, 5._knd/12._knd,  3._knd/4._knd  ]
  real(knd), parameter :: RK_rho(3)   = [ 0._knd,       -17._knd/60._knd, -5._knd/12._knd ]
  integer, parameter   :: RK_stages = 3

end module RK3
