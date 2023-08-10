module TurbInlet

  use Parameters
  use ArrayUtilities
  use rng_par_zig
#ifdef PAR  
  use custom_par, only: PAR_COMM_NULL
#endif  
  !$ use omp_lib
  
  implicit none

  private
  public :: turbulence_generator, turbulence_generator_nesting, &
            default_turbulence_generator, &
            GetInletFromFile

  type turbulence_generator
    integer :: direction

    real(knd) :: T_lag
    real(knd) :: L_y
    real(knd) :: L_z

    real(knd) :: Ustar_surf_inlet, &
                 stress_gradient_inlet, &
                 U_ref_inlet, &
                 z_ref_inlet, &
                 z0_inlet, &
                 power_exponent_inlet

    integer :: filtny, filtnz, bigNy, bigNz

    real(knd), allocatable, dimension(:) :: expsy, expsz
    real(knd) :: compat
    real(knd), allocatable, dimension(:,:,:) :: Ru, Rv, Rw !arrays of randoms
    real(knd), allocatable, dimension(:,:,:) :: Psiu, Psiv, Psiw
    real(knd), allocatable, dimension(:,:,:) :: bfilt !filter coefficients (ii,jj,kk,kz)


    real(knd),allocatable,dimension(:)     :: Ustar_inlet !friction velocity profile at inlet
    real(knd),allocatable,dimension(:,:)   :: Uinavg, Vinavg, Winavg !mean values of U,V,W at inflow
    real(knd),allocatable,dimension(:,:,:) :: transform_tensor
    real(knd),dimension(1:3,1:3) :: relative_stress
#ifdef PAR    
    integer :: comm = PAR_COMM_NULL
#endif    
  contains
    procedure :: init => turbulence_generator_init
    procedure :: time_step => turbulence_generator_time_step
    procedure, private :: init_turbulence_profiles => turbulence_generator_init_turbulence_profiles
    procedure, private :: init_mean_profiles => turbulence_generator_init_mean_profiles
    procedure, private :: bound_Uin => bound_Uin
    procedure, private :: get_mean_profile_from_file => turbulence_generator_get_mean_profile_from_file
    procedure, private :: get_turbulence_profile_from_file => turbulence_generator_get_turbulence_profile_from_file
  end type

  

  type, extends(turbulence_generator) :: turbulence_generator_nesting
   real(knd), allocatable, dimension(:,:) :: sgs_tke
  contains
    procedure, private :: init_turbulence_profiles => nesting_init_turbulence_profiles
    procedure, private :: init_mean_profiles => nesting_init_mean_profiles
    procedure, public :: update_turbulence_profiles => nesting_update_turbulence_profiles
  end type


  type TInlet
    real(knd),allocatable,dimension(:,:) :: U, V, W, temperature
  end type

  type(turbulence_generator) :: default_turbulence_generator

contains

  subroutine turbulence_generator_init(g)
#ifdef PAR
    use custom_par, only: iim, jim, kim, nxims, nyims, nzims, par_co_sum, par_co_min, &
                          comm_plane_xz,comm_plane_yz
    use exchange_par
#endif
    class(turbulence_generator), intent(inout) :: g
    integer :: i, j, k
    integer :: jlo, jup, klo, kup, ny
    real(knd) :: bysum, bzsum
    integer :: tid
    
    ! The synthetic turbulence generation method by Xie and Castro, 2008, https://dx.doi.org/10.1007/s10494-008-9151-5

#ifdef PAR
    if (g%direction==2) then
      if (.not. (jim==1 .or. jim==nyims)) return
      if (nyims > 1 .and. jim==nyims .and. Btype(No) /= TurbulentInletType) return
      g%comm = comm_plane_xz
    else
      !should not happen for 2D decomposition in Y and Z
      if (.not. (iim==1 .or. iim==nxims)) return
      if (nxims > 1 .and. iim==nxims .and. Btype(Ea) /= TurbulentInletType) return
      g%comm = comm_plane_yz
    end if
#endif

    call g%init_turbulence_profiles

    call g%init_mean_profiles

#ifdef PAR
    g%filtny = min(max(nint(g%L_y / dymin),1), ceiling(1._knd * gPrny / 3),Prny)
    g%filtnz = min(max(nint(g%L_z / dzmin),1), ceiling(1._knd * gPrnz / 3),Prnz)

    g%filtny = par_co_min(g%filtny, g%comm)
    g%filtnz = par_co_min(g%filtnz, g%comm)
#else
    g%filtny = min(max(nint(g%L_y / dymin),1), ceiling(1._knd * Prny / 3))
    g%filtnz = min(max(nint(g%L_z / dzmin),1), ceiling(1._knd * Prnz / 3))
#endif

    if (g%direction==2) then
      ny = Prnx
    else
       ny = Prny
    end if
    
    g%bigNy = 2 * g%filtny
    g%bigNz = 2 * g%filtnz

    jlo = -g%bigNy + 1
    jup = ny + g%bigNy

    klo = -g%bigNz + 1
    kup = Prnz + g%bigNz
      

        

    allocate(g%Ru(jlo:jup,klo:kup,2))
    allocate(g%Rv(jlo:jup,klo:kup,2))
    allocate(g%Rw(jlo:jup,klo:kup,2))
    allocate(g%Psiu(1:ny,1:Prnz,1:2))
    allocate(g%Psiv(1:ny,1:Prnz,1:2))
    allocate(g%Psiw(1:ny,1:Prnz,1:2))

    !$omp parallel private(j,k,tid)
    tid = 0
    !$ tid = omp_get_thread_num()
    
    !$omp do collapse(2)
    do k = klo, kup
     do j = jlo, jup
       call rng_norm(g%Ru(j,k,1), tid)
       call rng_norm(g%Rv(j,k,1), tid)
       call rng_norm(g%Rw(j,k,1), tid)
     end do
    end do

    !$omp do collapse(2)
    do k = klo, kup
     do j = jlo, jup
       call rng_norm(g%Ru(j,k,2), tid)
       call rng_norm(g%Rv(j,k,2), tid)
       call rng_norm(g%Rw(j,k,2), tid)
     end do
    end do

    !$omp end parallel


    if ((g%direction==2 .and. ((Btype(We)==BC_PERIODIC).or.(Btype(Ea)==BC_PERIODIC))) .or. &
        (g%direction==1 .and. ((Btype(So)==BC_PERIODIC).or.(Btype(No)==BC_PERIODIC)))) then
      !$omp parallel workshare
      forall(k = klo:kup, j = jlo:0)
        g%Ru(j,k,1:2) = g%Ru(j+ny,k,1:2)
        g%Rv(j,k,1:2) = g%Rv(j+ny,k,1:2)
        g%Rw(j,k,1:2) = g%Rw(j+ny,k,1:2)
      end forall
      forall(k = klo:kup,j = ny+1:jup)
        g%Ru(j,k,1:2) = g%Ru(j-ny,k,1:2)
        g%Rv(j,k,1:2) = g%Rv(j-ny,k,1:2)
        g%Rw(j,k,1:2) = g%Rw(j-ny,k,1:2)
      end forall
      !$omp end  parallel workshare
    end if

    if  ((Btype(Bo)==BC_PERIODIC).or.(Btype(To)==BC_PERIODIC)) then
      !$omp parallel workshare
      forall(k = klo:0, j = jlo:jup)
        g%Ru(j,k,1:2) = g%Ru(j,k+Prnz,1:2)
        g%Rv(j,k,1:2) = g%Rv(j,k+Prnz,1:2)
        g%Rw(j,k,1:2) = g%Rw(j,k+Prnz,1:2)
      end forall
      forall(k = Prnz+1:kup,j = jlo:jup)
        g%Ru(j,k,1:2) = g%Ru(j,k-Prnz,1:2)
        g%Rv(j,k,1:2) = g%Rv(j,k-Prnz,1:2)
        g%Rw(j,k,1:2) = g%Rw(j,k-Prnz,1:2)
      end forall
      !$omp end  parallel workshare
    end if

#ifdef PAR
    if (g%direction==2) then
      do i=1,2
        call par_exchange_boundaries_xz(g%Ru(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
        call par_exchange_boundaries_xz(g%Rv(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
        call par_exchange_boundaries_xz(g%Rw(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
      end do
    else
      do i=1,2
        call par_exchange_boundaries_yz(g%Ru(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
        call par_exchange_boundaries_yz(g%Rv(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
        call par_exchange_boundaries_yz(g%Rw(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
      end do
    end if
#endif


    allocate(g%bfilt(-g%bigNy:g%bigNy, -g%bigNz:g%bigNz,1))
    allocate(g%expsy(-g%bigNy:g%bigNy), g%expsz(-g%bigNz:g%bigNz))

    bysum = 0
    do i = -g%bigNy, g%bigNy
     g%expsy(i) = exp(-pi * abs(i) / (g%bigNy))
     bysum = bysum + g%expsy(i)**2
    end do
    bysum = sqrt(bysum)
    g%expsy = g%expsy / bysum


    bzsum = 0
    do i = -g%bigNz, g%bigNz
     g%expsz(i) = exp(-pi * abs(i) / (g%bigNz))
     bzsum = bzsum + g%expsz(i)**2
    end do
    bzsum = sqrt(bzsum)
    g%expsz = g%expsz / bzsum

    !$omp parallel private(j,k)
    !$omp do collapse(2)
    do k = -g%bigNz, g%bigNz
      do j = -g%bigNy, g%bigNy
         g%bfilt(j,k,1) = g%expsy(j) * g%expsz(k)
      end do
    end do
    !$omp end  do
    !$omp workshare
    forall(k = 1:Prnz, j = 1:ny)
         g%Psiu(j,k,1) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Ru(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,1))
         g%Psiv(j,k,1) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Rv(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,1))
         g%Psiw(j,k,1) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Rw(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,1))
    end forall
    !$omp end  workshare nowait
    

    if (allocated(g%Vinavg).and.g%direction==2) then
      !$omp workshare
      g%compat = sum(g%Vinavg(1:Vnx,1:Vnz))
      !$omp end  workshare
    else if (allocated(g%Uinavg)) then
      !$omp workshare
      g%compat = sum(g%Uinavg(1:Uny,1:Unz))
      !$omp end  workshare
    else
      !$omp single
      g%compat = 0
      !$omp end single
    end if
    !$omp end  parallel
#ifdef PAR
    g%compat = par_co_sum(g%compat, g%comm)
#endif

  end subroutine


  subroutine turbulence_generator_time_step(g, Uin, Vin, Win, dt)
#ifdef PAR
    use custom_par, only: iim, jim, kim, nxims, nyims, nzims, par_co_sum, par_co_min
    use exchange_par
#endif
    class(turbulence_generator), intent(inout) :: g
    real(knd), intent(out) :: Uin(-2:,-2:), Vin(-2:,-2:), Win(-2:,-2:)
    real(knd), intent(in) :: dt
    integer :: i, j, k
    integer :: jlo, jup, klo, kup, ny
    real(knd) :: Ui, Vi, Wi, p
    integer :: tid
    
#ifdef PAR
    if (g%direction==2) then
      if (.not. (jim==1 .or. jim==nyims)) return
      if (nyims > 1 .and. jim==nyims .and. Btype(No) /= TurbulentInletType) return
    else
      !should not happen for 2D decomposition in Y and Z
      if (.not. (iim==1 .or. iim==nxims)) return
      if (nxims > 1 .and. iim==nxims .and. Btype(Ea) /= TurbulentInletType) return
    end if
#endif

    if (g%direction==2) then
      ny = Prnx
    else
      ny = Prny
    end if

    jlo = -g%bigNy + 1
    jup = ny + g%bigNy
    
    klo = -g%bigNz + 1
    kup = Prnz + g%bigNz


    !$omp parallel do private(j,k,Ui,Vi,Wi)  
    do k = 1, Prnz
      do j = 1, ny
         g%Psiu(j,k,2) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Ru(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,2))
         g%Psiv(j,k,2) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Rv(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,2))
         g%Psiw(j,k,2) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Rw(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,2))
      end  do
    end  do
    !$omp end  parallel do

    call multiply(g%Psiu(:,:,1), exp(-pi * dt / (2._knd * g%T_lag)))
    call add_multiplied(g%Psiu(:,:,1), g%Psiu(:,:,2), sqrt(1 - exp(-pi * dt / (g%T_lag))))

    call multiply(g%Psiv(:,:,1),exp(-pi * dt / (2._knd * g%T_lag)))
    call add_multiplied(g%Psiv(:,:,1), g%Psiv(:,:,2), sqrt(1 - exp(-pi * dt / (g%T_lag))))

    call multiply(g%Psiw(:,:,1),exp(-pi * dt / (2._knd * g%T_lag)))
    call add_multiplied(g%Psiw(:,:,1), g%Psiw(:,:,2), sqrt(1 - exp(-pi * dt / (g%T_lag))))
    
    call set(Uin, 0._knd)
    call set(Vin, 0._knd)
    call set(Win, 0._knd)

    if (allocated(g%Uinavg) .and. allocated(g%Vinavg) .and. allocated(g%Winavg)) then
      !$omp parallel do private(j,k,Ui,Vi,Wi)
      do k = 1, Prnz
       do j = 1, ny
        Ui = g%Psiu(j,k,1)
        Vi = g%Psiv(j,k,1)
        Wi = g%Psiw(j,k,1)

        Uin(j,k) = g%Uinavg(j,k) + g%transform_tensor(1,j,k) * Ui   !a12,a13,a23 = 0

        Vin(j,k) = g%Vinavg(j,k) + g%transform_tensor(2,j,k) * Ui &
                                 + g%transform_tensor(3,j,k) * Vi

        Win(j,k) = g%Winavg(j,k) + g%transform_tensor(4,j,k) * Ui &
                                 + g%transform_tensor(5,j,k) * Vi &
                                 + g%transform_tensor(6,j,k) * Wi
       end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(j,k,Ui,Vi,Wi)
      do k = 1, Prnz
       do j = 1, Prny
        Ui = g%Psiu(j,k,1)
        Vi = g%Psiv(j,k,1)
        Wi = g%Psiw(j,k,1)

        Uin(j,k) = g%transform_tensor(1,j,k) * Ui   !a12,a13,a23 = 0

        Vin(j,k) =  g%transform_tensor(2,j,k) * Ui &
                  + g%transform_tensor(3,j,k) * Vi

        Win(j,k) =  g%transform_tensor(4,j,k) * Ui &
                  + g%transform_tensor(5,j,k) * Vi &
                  + g%transform_tensor(6,j,k) * Wi
       end do
      end do
      !$omp end parallel do
    end if

    if (g%direction==2) then
      !$omp parallel workshare
      p = sum(Vin(1:Vnx,1:Vnz))
      !$omp end parallel workshare      
    else
      !$omp parallel workshare
      p = sum(Uin(1:Uny,1:Unz))
      !$omp end parallel workshare
    end if


#ifdef PAR
    p = par_co_sum(p, g%comm)
#endif
    p = p - g%compat                  !To ensure the g%compatibility condition.
    
    if (g%direction==2) then
      p = p / (gVnx * gVnz)      
      call add(Vin, p)
    else
      p = p / (gUny * gUnz)
      call add(Uin, p)
    end if


    call g%bound_Uin(1, Uin)

    call g%bound_Uin(2, Vin)

    call g%bound_Uin(3, Win)

    !$omp parallel private(j,k,tid)
    tid = 0
    !$ tid = omp_get_thread_num()
    
    !$omp do collapse(2) 
    do k = klo, kup
     do j = jlo, jup
       call rng_norm(g%Ru(j,k,2), tid)
       call rng_norm(g%Rv(j,k,2), tid)
       call rng_norm(g%Rw(j,k,2), tid)
     end do
    end do
    !$omp end do

    if ((g%direction==2 .and. ((Btype(We)==BC_PERIODIC).or.(Btype(Ea)==BC_PERIODIC))) .or. &
        (g%direction==1 .and. ((Btype(So)==BC_PERIODIC).or.(Btype(No)==BC_PERIODIC)))) then

      !$omp do collapse(2)
      do k = klo, kup
        do j = jlo, 0
          g%Ru(j,k,1:2) = g%Ru(j+ny,k,1:2)
          g%Rv(j,k,1:2) = g%Rv(j+ny,k,1:2)
          g%Rw(j,k,1:2) = g%Rw(j+ny,k,1:2)
        end do
      end do
      !$omp end do nowait
      !$omp do collapse(2)
      do k = klo, kup
        do j = ny+1, jup
          g%Ru(j,k,1:2) = g%Ru(j-ny,k,1:2)
          g%Rv(j,k,1:2) = g%Rv(j-ny,k,1:2)
          g%Rw(j,k,1:2) = g%Rw(j-ny,k,1:2)
        end do
      end do
      !$omp end do
    end if
    
    
    if  ((Btype(Bo)==BC_PERIODIC) .or. (Btype(To)==BC_PERIODIC)) then
      !$omp do collapse(2)
      do k = klo, 0
        do j = jlo, jup
          g%Ru(j,k,1:2) = g%Ru(j,k+Prnz,1:2)
          g%Rv(j,k,1:2) = g%Rv(j,k+Prnz,1:2)
          g%Rw(j,k,1:2) = g%Rw(j,k+Prnz,1:2)
        end do
      end do
      !$omp end do nowait
      !$omp do collapse(2)
      do k = Prnz+1, kup
        do j = jlo, jup
          g%Ru(j,k,1:2) = g%Ru(j,k-Prnz,1:2)
          g%Rv(j,k,1:2) = g%Rv(j,k-Prnz,1:2)
          g%Rw(j,k,1:2) = g%Rw(j,k-Prnz,1:2)
        end do
      end do
      !$omp end do
    end if
    !$omp end parallel

#ifdef PAR
    if (g%direction==2) then
      do i=1,2
        call par_exchange_boundaries_xz(g%Ru(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
        call par_exchange_boundaries_xz(g%Rv(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
        call par_exchange_boundaries_xz(g%Rw(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
      end do
    else
      do i=1,2
        call par_exchange_boundaries_yz(g%Ru(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
        call par_exchange_boundaries_yz(g%Rv(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
        call par_exchange_boundaries_yz(g%Rw(:,:,i), ny, Prnz, Btype, &
                                        jlo, klo, g%bigNy, g%bigNz)
      end do
    end if
#endif

  end subroutine turbulence_generator_time_step












  subroutine turbulence_generator_init_turbulence_profiles(g)
    class(turbulence_generator), intent(inout) :: g
    integer :: j,k

    allocate(g%Ustar_inlet(1:Prnz))
    if (g%direction==2) then
      allocate(g%transform_tensor(1:6,1:Prnx,1:Prnz))
    else
      allocate(g%transform_tensor(1:6,1:Prny,1:Prnz))
    end if
    
    if (profiletype==PROFILE_FROM_FILE) then
      call g%get_turbulence_profile_from_file("inlet_stress_profile.txt")
      return
    end if
    
    if ((profiletype==PROFILE_LOGARITHMIC .and. g%Ustar_surf_inlet<=0) .or. &
        (profiletype==PROFILE_POWER_LAW .and. g%Ustar_surf_inlet<=0)) then
      g%Ustar_surf_inlet = abs(Karman * g%U_ref_inlet / log(g%z0_inlet / g%z_ref_inlet))
    end if


    do k = 1, Prnz
      g%Ustar_inlet(k) = g%Ustar_surf_inlet * sqrt(max(1 + g%stress_gradient_inlet * zPr(k),1E-5_knd))
    end do


    do k = 1, Prnz
        ! tt1 = a11,tt2 = a21,tt3 = a22, tt4 = a31, tt5 = a32, tt6 = a33
        g%transform_tensor(1,:,k) = sqrt((g%Ustar_inlet(k)**2) * g%relative_stress(1,1))
        g%transform_tensor(2,:,k) = (g%Ustar_inlet(k)**2) * g%relative_stress(2,1) / g%transform_tensor(1,:,k)
        g%transform_tensor(3,:,k) = sqrt((g%Ustar_inlet(k)**2) * g%relative_stress(2,2) - g%transform_tensor(2,:,k)**2)
        g%transform_tensor(4,:,k) = (g%Ustar_inlet(k)**2) * g%relative_stress(3,1) / g%transform_tensor(1,:,k)
        g%transform_tensor(5,:,k) = ((g%Ustar_inlet(k)**2) * g%relative_stress(3,2)-&
                                 g%transform_tensor(2,:,k) * g%transform_tensor(4,:,k)) / g%transform_tensor(3,:,k)
        g%transform_tensor(6,:,k) = sqrt((g%Ustar_inlet(k)**2) * g%relative_stress(3,3)-&
                                 g%transform_tensor(4,:,k)**2 - g%transform_tensor(5,:,k)**2)

    end do
  end subroutine

  subroutine turbulence_generator_init_mean_profiles(g)
    use ieee_exceptions
    class(turbulence_generator), intent(inout) :: g
    real(knd) :: Ustar_prof, utmp
    integer :: k, maxj, maxk
    logical :: fix_direction
    
    CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .TRUE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .TRUE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .TRUE.)
    
    if (g%direction==2) then
      allocate(g%Uinavg(-2:Unx+3,-2:Unz+3), g%Vinavg(-2:Vnx+3,-2:Vnz+3), g%Winavg(-2:Wnx+3,-2:Wnz+3))
      g%Vinavg = 0
      g%Winavg = 0
    else
      allocate(g%Uinavg(-2:Uny+3,-2:Unz+3), g%Vinavg(-2:Vny+3,-2:Vnz+3), g%Winavg(-2:Wny+3,-2:Wnz+3))
      g%Vinavg = 0
      g%Winavg = 0
    end if

    if  (profiletype==PROFILE_CONSTANT) then

      do k = 1, Prnz

        g%Uinavg(:,k) = Uinlet_vec(1)
        g%Vinavg(:,k) = Uinlet_vec(2)
        g%Winavg(:,k) = Uinlet_vec(3)

      end do

     fix_direction = .false.

    else if (profiletype==PROFILE_LOGARITHMIC) then

      if (g%U_ref_inlet/=0.and.g%z_ref_inlet>0) then

        Ustar_prof = g%U_ref_inlet * Karman / log(g%z_ref_inlet / g%z0_inlet)

        do k = 1, Prnz
          utmp = (Ustar_prof / Karman) * log(zPr(k) / g%z0_inlet)
          if (sign(1._knd,Ustar_prof) * utmp<abs(Ustar_prof) / 2) &
            utmp = (Ustar_prof / 2) * zPr(k) / (g%z0_inlet * 1.22)
          g%Uinavg(:,k) = utmp
        end do

        fix_direction = .true.
        
     else

        utmp = (g%Ustar_inlet(1) / Karman) * log(zPr(1) / g%z0_inlet)
        g%Uinavg(:,1) = utmp
        do k = 2, Prnz
          utmp = (g%Ustar_inlet(k) / Karman) * log(zPr(k) / zPr(k-1)) + utmp
          if (sign(1._knd,g%Ustar_inlet(1)) * utmp < abs(g%Ustar_inlet(1)) / 2) &
            utmp = (g%Ustar_inlet(1) / 2) * zPr(k) / (g%z0_inlet * 1.22)
          g%Uinavg(:,k) = utmp
        end do

        fix_direction = .true.
        
      end  if

    else if (profiletype==PROFILE_POWER_LAW) then

      do k = 1, Prnz
        g%Uinavg(:,k) = g%U_ref_inlet * (zPr(k) / g%z_ref_inlet)**g%power_exponent_inlet
      end do

      fix_direction = .false.
      
    else if (profiletype==PROFILE_FROM_FILE) then
    
      call g%get_mean_profile_from_file("inlet_mean_profile.txt")
      
      fix_direction = .false.

    else

      g%Uinavg = 0

    end if

    if (fix_direction) then
      maxj = min(ubound(g%Uinavg,1),ubound(g%Vinavg,1),ubound(g%Winavg,1))
      maxk = min(ubound(g%Uinavg,2),ubound(g%Vinavg,2),ubound(g%Winavg,2))
      if (norm2(Uinlet_vec(2:3))>0) then
        g%Vinavg(-2:maxj,-2:maxk) = (Uinlet_vec(2)/Uinlet) * g%Uinavg(-2:maxj,-2:maxk)
        g%Winavg(-2:maxj,-2:maxk) = (Uinlet_vec(3)/Uinlet) * g%Uinavg(-2:maxj,-2:maxk)
        g%Uinavg(-2:maxj,-2:maxk) = (Uinlet_vec(1)/Uinlet) * g%Uinavg(-2:maxj,-2:maxk)
      else if (windangle/=0) then
        g%Vinavg(-2:maxj,-2:maxk) = &
            g%Uinavg(-2:maxj,-2:maxk) * sin(windangle / 180._knd * pi)
        g%Uinavg(:,:) = g%Uinavg * cos(windangle / 180._knd * pi)
      end if
    end if

    if (g%direction==2) then
      where(Utype(:,0,:) > 0) g%Uinavg = 0
      where(Vtype(:,0,:) > 0) g%Vinavg = 0
      where(Wtype(:,0,:) > 0) g%Winavg = 0
    else
      where(Utype(0,:,:) > 0) g%Uinavg = 0
      where(Vtype(0,:,:) > 0) g%Vinavg = 0
      where(Wtype(0,:,:) > 0) g%Winavg = 0
    end if

  end subroutine turbulence_generator_init_mean_profiles
  
  
  subroutine turbulence_generator_get_mean_profile_from_file(g, fname)
    use ieee_exceptions
    use Interpolation
    class(turbulence_generator), intent(inout) :: g
    character(*), intent(in) :: fname
    character(1024) :: line
    integer :: nh, n, i, j, k
    integer :: io, io1, io2, unit
    integer :: components
    real(knd) :: r4(4)
    real(knd), allocatable :: z(:), up(:), vp(:), wp(:), cu(:,:), cv(:,:), cw(:,:)
    
    
    CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .TRUE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .TRUE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .TRUE.)
    
    open(newunit=unit,file=fname,status="old",action="read",iostat = io)
    if (io/=0) then
      call error_stop("Cannot open profile of the mean inlet profile in: "//fname)
    end if

    n = 0
    nh = 0
    components = 0
    do
      read(unit,'(a)',iostat=io) line
      if (io/=0) exit
      
      if (scan(line(1:1),"#!")>0.and.n==0) then
        nh = nh + 1
        cycle
      end if
      if (components==1) then
        read(line,*,iostat=io) r4(1:2)
        if (io/=0) exit
      else if (components==2) then
        read(line,*,iostat=io) r4(1:3)
        if (io/=0) exit
      else if (components==3) then
        read(line,*,iostat=io) r4
        if (io/=0) exit
      else if (components == 0) then
        read(line,*,iostat=io) r4
        read(line,*,iostat=io2) r4(1:3)
        read(line,*,iostat=io1) r4(1:2)
        if (io==3) then
          components = 3
        else if (io2==0) then
          components = 2
        else if (io1==0) then
          components = 1
        else
          exit
        end if
      end if

      n = n + 1
    end do
    
    rewind(unit)
    
    if (n>0) then
      allocate(z(n), up(n), vp(n), wp(n))
      allocate(cu(0:4,n), cv(0:4,n), cw(0:4,n))
      do i = 1, nh
        read(unit,'(a)')
      end do
      do i = 1, n
        read(unit,'(a)') line
        if (components==3) then
          read(line, *) z(i), up(i), vp(i), wp(i)
        else if (components==2) then
          read(line, *) z(i), up(i), vp(i)
          wp(i) = 0        
        else if (components==1) then
          read(line, *) z(i), up(i)
          vp(i) =0; wp(i) = 0        
        end if
      end do
    else
      write(*,*) "Error, mean inlet profile file '",trim(fname),"'empty."
      call error_stop
    end if

    if (n > 1) then
      call linear_interpolation(z, up, cu)
      call linear_interpolation(z, vp, cv)
      call linear_interpolation(z, wp, cw)

      j = 1
      do k = 1, Prnz
        g%Uinavg(:,k) = linear_interpolation_eval(zPr(k), z, cu, j)
        g%Vinavg(:,k) = linear_interpolation_eval(zPr(k), z, cv, j)
        g%Winavg(:,k) = linear_interpolation_eval(zW(k),  z, cw, j)        
      end do

    else
      g%Uinavg = up(1)
      g%Vinavg = vp(1)
      g%Winavg = wp(1)
    end if

  end subroutine turbulence_generator_get_mean_profile_from_file
  

  subroutine turbulence_generator_get_turbulence_profile_from_file(g, fname)
    use Interpolation
    use Strings
    class(turbulence_generator), intent(inout) :: g
    character(*), intent(in) :: fname
    character(1024) :: line
    integer :: nh, n, i, j, k, comp
    integer :: io, io4, unit
    integer :: components
    real(knd) :: r7(7)
    real(knd), allocatable :: z(:), stress(:,:), c_stress(:,:,:)
    real(knd) :: uu, vv, ww, uv, uw, vw, tmp
    
    open(newunit=unit,file=fname,status="old",action="read",iostat = io)
    if (io/=0) then
      call error_stop("Cannot open profile of the mean inlet profile in: "//fname)
    end if

    n = 0
    nh = 0
    components = 0
    do
      read(unit,'(a)',iostat=io) line
      if (io/=0) exit
      
      if (scan(line(1:1),"#!")>0.and.n==0) then
        nh = nh + 1
        cycle
      end if
      if (components==4) then
        read(line,*,iostat=io) r7(1:5)
        if (io/=0) exit
      else if (components==6) then
        read(line,*,iostat=io) r7
        if (io/=0) exit
      else if (components == 0) then
        read(line,*,iostat=io) r7
        read(line,*,iostat=io4) r7(1:5)
        if (io==0) then
          components = 6
        else if (io4==0) then
          components = 4
        else
          exit
        end if
      end if

      n = n + 1
    end do
    
    rewind(unit)
    
    if (n>0) then
      allocate(z(n), stress(n,6))  !order uu, vv, ww, uv, uw, vw
      allocate(c_stress(0:4,n,6))
      do i = 1, nh
        read(unit,'(a)')
      end do
      do i = 1, n
        read(unit,'(a)') line
        if (components==4) then
          read(line, *) z(i), stress(i,[1,2,3,5])
          stress(i,4) = 0
          stress(i,6) = 0
        else
          read(line, *) z(i), stress(i,:)
        end if
      end do
    else
      write(*,*) "Error, mean inlet profile file '",trim(fname),"'empty."
      call error_stop
    end if

    if (n > 1) then
      do comp = 1, 6
        call linear_interpolation(z, stress(:,comp), c_stress(:,:,comp))
      end do

      j = 1
      do k = 1, Prnz
        !a small inconsistency in z due to the staggered grid
        uu = linear_interpolation_eval(zPr(k), z, c_stress(:,:,1), j)
        vv = linear_interpolation_eval(zPr(k), z, c_stress(:,:,2), j)
        ww = linear_interpolation_eval(zW(k),  z, c_stress(:,:,3), j)
        uv = linear_interpolation_eval(zPr(k), z, c_stress(:,:,4), j)
        uw = linear_interpolation_eval(zPr(k), z, c_stress(:,:,5), j)
        vw = linear_interpolation_eval(zPr(k), z, c_stress(:,:,6), j)

        ! different numbering! See Xie and Castro, 2008
        ! tt1 = a11,tt2 = a21,tt3 = a22, tt4 = a31, tt5 = a32, tt6 = a33
        if (uu <= 0) then
          write(*,*) "interpolated stresses at z =",zPr(k)
          write(*,*) "uu,vv,ww,uv,uw,vw:",uu, vv, ww, uv, uw, vw
          call error_stop("Error in the interpolated inlet stress profile.&
                         & Stress component uu not positive at&
                         & k = " // itoa(k))
        end if
        if (vv <= 0) then
          write(*,*) "interpolated stresses at z =",zPr(k)
          write(*,*) "uu,vv,ww,uv,uw,vw:",uu, vv, ww, uv, uw, vw
          call error_stop("Error in the interpolated inlet stress profile.&
                         & Stress component vv not positive at&
                         & k = " // itoa(k))
        end if
        if (ww <= 0) then
          write(*,*) "interpolated stresses at z =",zPr(k)
          write(*,*) "uu,vv,ww,uv,uw,vw:",uu, vv, ww, uv, uw, vw
          call error_stop("Error in the interpolated inlet stress profile.&
                         & Stress component ww not positive at&
                         & k = " // itoa(k))
        end if
        g%transform_tensor(1,1,k) = sqrt(uu)
        g%transform_tensor(2,1,k) = uv / g%transform_tensor(1,1,k)
        
        tmp = vv - g%transform_tensor(2,1,k)**2
        if (tmp <= 0) then
          write(*,*) "interpolated stresses at z =",zPr(k)
          write(*,*) "uu,vv,ww,uv,uw,vw:",uu, vv, ww, uv, uw, vw
          call error_stop("Error in the interpolated inlet stress profile.&
                         & Argument of sqrt() for tensor component 3 (vv) not positive at&
                         & k = " // itoa(k))
        end if
        g%transform_tensor(3,1,k) = sqrt(tmp)

        g%transform_tensor(4,1,k) = uw / g%transform_tensor(1,1,k)
        
        g%transform_tensor(5,1,k) = (vw - &
                                     g%transform_tensor(2,1,k) * g%transform_tensor(4,1,k)) / &
                                    g%transform_tensor(3,1,k)

        tmp = ww - g%transform_tensor(4,1,k)**2 - g%transform_tensor(5,1,k)**2
        if  (tmp < 0) then
          write(*,*) "interpolated stresses at z =",zPr(k)
          write(*,*) "uu,vv,ww,uv,uw,vw:",uu, vv, ww, uv, uw, vw
          call error_stop("Error in the interpolated inlet stress profile.&
                         & Argument of sqrt() for tensor component 6 (ww) negative at&
                         & k = " // itoa(k))
        end if
        g%transform_tensor(6,1,k) = sqrt(tmp)
        
        do comp = 1, 6
          g%transform_tensor(comp,:,k) = g%transform_tensor(comp,1,k)
        end do
      end do


    else
      uu = stress(1,1)
      vv = stress(1,2)
      ww = stress(1,3)
      uv = stress(1,4)
      uw = stress(1,5)
      vw = stress(1,6)
      
      g%transform_tensor(1,:,:) = sqrt(uu)
      g%transform_tensor(2,:,:) = uv / g%transform_tensor(1,:,:)
      g%transform_tensor(3,:,:) = sqrt(vv - g%transform_tensor(2,:,:)**2)
      g%transform_tensor(4,:,:) = uw / g%transform_tensor(1,:,:)
      g%transform_tensor(5,:,:) = (vw - &
                                      g%transform_tensor(2,:,:) * &
                                      g%transform_tensor(4,:,:)) &
                                  / g%transform_tensor(3,:,:)
      g%transform_tensor(6,:,:) = sqrt(ww - &
                                        g%transform_tensor(4,:,:)**2 - &
                                        g%transform_tensor(5,:,:)**2)
    end if

  end subroutine turbulence_generator_get_turbulence_profile_from_file
  

  subroutine nesting_init_turbulence_profiles(g)
    class(turbulence_generator_nesting), intent(inout) :: g
    integer :: j,k
    real(knd) :: stress

    g%relative_stress = 0
    g%relative_stress(1,1) = 1
    g%relative_stress(2,2) = 1
    g%relative_stress(3,3) = 1

    allocate(g%transform_tensor(1:6,1:Prny,1:Prnz))
    g%transform_tensor = 0

    do k = 1, Prnz
     do j = 1, Prny  ! tt1 = a11,tt2 = a21,tt3 = a22, tt4 = a31, tt5 = a32, tt6 = a33
        stress = sqrt(2*g%sgs_tke(j,k)/3)
        
        g%transform_tensor(1,j,k) = stress
        g%transform_tensor(3,j,k) = stress
        g%transform_tensor(6,j,k) = stress
     end do
    end do
  end subroutine

  subroutine nesting_update_turbulence_profiles(g)
    class(turbulence_generator_nesting), intent(inout) :: g
    integer :: j,k
    real(knd) :: stress

    do k = 1, Prnz
     do j = 1, Prny  ! tt1 = a11,tt2 = a21,tt3 = a22, tt4 = a31, tt5 = a32, tt6 = a33
        stress = sqrt(2*g%sgs_tke(j,k)/3)

        g%transform_tensor(1,j,k) = stress
        g%transform_tensor(3,j,k) = stress
        g%transform_tensor(6,j,k) = stress
     end do
    end do
  end subroutine

  subroutine nesting_init_mean_profiles(g)
    class(turbulence_generator_nesting), intent(inout) :: g
    !leave the avg profiles unallocated
  end subroutine nesting_init_mean_profiles


  subroutine GetInletFromFile(t)
    real(TIM),intent(in):: t
    integer,save:: called = 0
    integer :: Prny2, Prnz2, Vny2, Wnz2
    real(knd) :: dx2
    character(12):: fname
    integer,save:: inletfnum

    type(TInlet),pointer,save:: In1=>null(),In2=>null(),Inp=>null() !pointer to inlets to avoid transfers, time(In1)<time(In2)
    real(TIM),save:: t1,t2 !time if file1, file 2
    real(TIM) :: tp
    integer,save:: io

    real(knd) :: c1,c2

    if (called==0) then
       open(102,file="inletframeinfo.unf",form='unformatted',status='old',action='read',iostat = io)
       if (io/=0) then
        write(*,*) 'Error while opening file inletframeinfo.unf'
       end if
       read(102) Prny2, Prnz2  !for check of consistency of grids before use
       read(102) Vny2
       read(102) Wnz2
       read(102) dx2

       allocate(In1,In2)
       allocate(In1%U(Uny,Unz),In1%V(Vny,Vnz),In1%W(Wny,Wnz))
       allocate(In2%U(Uny,Unz),In2%V(Vny,Vnz),In2%W(Wny,Wnz))
       if (enable_buoyancy) allocate(In1%temperature(Prny,Prnz),In2%temperature(Prny,Prnz))

       if ((Prny/=Prny2).or.(Prnz/=Prnz2).or.(Vny/=Vny2).or.(Wnz/=Wnz2).or.((dx2-dxPr(0))/dx2>0.1)) then
        call error_stop("Mismatch of computational grid and inlet file.")
       end if
       called = 1
       inletfnum = 1

       fname(1:5)="frame"
       write(fname(6:8),"(I3.3)") inletfnum
       fname(9:12)=".unf"

       open(11,file = fname,form='unformatted',status='old',action='read',iostat = io)
       if (io/=0) then
        write(*,*) "Error while opening file ", fname
        call error_stop
       end if
       read(102) t1
       call ReadBC_INLET_FROM_FILE(11,In1)
       close(11)

       inletfnum = inletfnum+1

       fname(1:5)="frame"
       write(fname(6:8),"(I3.3)") inletfnum
       fname(9:12)=".unf"

        open(11,file = fname,form='unformatted',status='old',action='read',iostat = io)
       if (io/=0) then
        write(*,*) "Error while opening file ", fname
        call error_stop
       end if
       read(102) t2
       call ReadBC_INLET_FROM_FILE(11,In2)
       close(11)
    end if

    io = 0

    do while (t2<t.and.io==0)
      read(102,iostat = io) tp
      if (io==0) then
       t1 = t2
       t2 = tp

       inletfnum = inletfnum+1
       fname(1:5)="frame"
       write(fname(6:8),"(I3.3)") inletfnum
       fname(9:12)=".unf"

       open(11,file = fname,form='unformatted',status='old',action='read',iostat = io)
       if (io/=0) then
        write(*,*) "Error while opening file ", fname
        call error_stop
       end if
       Inp=>In1
       In1=>In2
       In2=>Inp
       call ReadBC_INLET_FROM_FILE(11,In2)
       close(11)
      end if
    end do


    if (t<=t1) then
     c1 = 1
     c2 = 0
    else if (t>=t2) then
     c1 = 0
     c2 = 1
    else
     c1=(t2-t)/(t2-t1)
     c2=(t-t1)/(t2-t1)
    end if
    write(*,*) "t",t1,t2,t
    write(*,*) "c",c1,c2

    Uin(1:Uny,1:Unz) = c1*In1%U+c2*In2%U
    Vin(1:Vny,1:Vnz) = c1*In1%V+c2*In2%V
    Win(1:Wny,1:Wnz) = c1*In1%W+c2*In2%W
    if (enable_buoyancy) TempIn(1:Prny,1:Prnz) = c1*In1%temperature+c2*In2%temperature

  end subroutine GetInletFromFile

  subroutine  ReadBC_INLET_FROM_FILE(unitnum,In)
    integer,intent(in):: unitnum
    type(TInlet),intent(inout)::In

    read(unitnum) In%U
    read(unitnum) In%V
    read(unitnum) In%W
    if (enable_buoyancy) then
         read(unitnum) In%temperature
    end if
  end subroutine ReadBC_INLET_FROM_FILE


  subroutine bound_Uin(g, component,Uin)
#ifdef PAR
    use exchange_par, only: par_exchange_boundaries_xz, par_exchange_boundaries_yz
#endif
    class(turbulence_generator), intent(in) :: g
    integer,  intent(in)    :: component
    real(knd),intent(inout) :: Uin(-2:,-2:)
    integer :: ny, nz, side_m, side_p
    
    if (g%direction==2) then
      if (component==1) then
        ny = Unx
        nz = Unz
      else if (component==2) then
        ny = Vnx
        nz = Vnz
      else
        ny = Wnx
        nz = Wnz
      end if
      
      side_m = We
      side_p = Ea
    else
      if (component==1) then
        ny = Uny
        nz = Unz
      else if (component==2) then
        ny = Vny
        nz = Vnz
      else
        ny = Wny
        nz = Wnz
      end if
      
      side_m = So
      side_p = No
    end if
   
#ifdef PAR
    if (g%direction==2) then
      call par_exchange_boundaries_xz(Uin, ny, nz, Btype, &
                                      -2, -2, 2, 2)
    else
      call par_exchange_boundaries_yz(Uin, ny, nz, Btype, &
                                      -2, -2, 2, 2)
    end if
#endif

    if (Btype(Bo)==BC_DIRICHLET) then
      if (component==3) then
        Uin(1:ny,0) = sideU(component,Bo)
        Uin(1:ny,-1) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,1))
      else
        Uin(1:ny,0) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,1))
        Uin(1:ny,-1) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,2))
      end  if
    else if (Btype(Bo)==BC_PERIODIC) then
        Uin(1:ny,-1:0) = Uin(1:ny,nz-1:nz)
    else if (Btype(Bo)==BC_NOSLIP.or.(component==3.and.Btype(Bo)==BC_FREESLIP)) then
      if (component==3) then
        Uin(1:ny,0) = 0
        Uin(1:ny,-1)=-Uin(1:ny,1)
      else
        Uin(1:ny,-1:0)=-Uin(1:ny,2:1:-1)
      end if
    else if (Btype(Bo)==BC_NEUMANN.or.(component/=3.and.Btype(Bo)==BC_FREESLIP)) then
        Uin(1:ny,0) = Uin(1:ny,1)
        Uin(1:ny,-1) = Uin(1:ny,1)
    else if (Btype(Bo)==BC_PERIODIC) then
        Uin(1:ny,-1:0) = Uin(1:ny,nz-1:nz)
    end if

    if (Btype(To)==BC_DIRICHLET) then
      if (component==3) then
        Uin(1:ny,nz+1) = sideU(component,To)
        Uin(1:ny,nz+2) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz))
      else
        Uin(1:ny,nz+1) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz))
        Uin(1:ny,nz+2) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz-1))
      end  if
    else if (Btype(To)==BC_PERIODIC) then
        Uin(1:ny,nz+1:nz+2) = Uin(1:ny,nz-1:nz)
    else if (Btype(To)==BC_NOSLIP.or.(component==3.and.(Btype(To)==BC_FREESLIP))) then
      if (component==3) then
        Uin(1:ny,nz+1) = 0
        Uin(1:ny,nz+2)=-Uin(1:ny,nz)
      else
        Uin(1:ny,nz+1:nz+2)=-Uin(1:ny,nz:nz-1:-1)
      end if
    else if (Btype(To)==BC_NEUMANN.or.(component/=3.and.(Btype(To)==BC_FREESLIP))) then
        Uin(1:ny,nz+1) = Uin(1:ny,nz)
        Uin(1:ny,nz+2) = Uin(1:ny,nz)
    else if (Btype(To)==BC_PERIODIC) then
        Uin(1:ny,nz+1:nz+2) = Uin(1:ny,1:2)
    end if

    if (Btype(side_m)==BC_DIRICHLET) then
      if (component==2) then
        Uin(0,-1:nz+2) = sideU(component,side_m)
        Uin(-1,-1:nz+2) = sideU(component,side_m)+(sideU(component,side_m)-Uin(1,-1:nz+2))
      else
        Uin(0,-1:nz+2) = sideU(component,side_m)+(sideU(component,side_m)-Uin(1,-1:nz+2))
        Uin(-1,-1:nz+2) = sideU(component,side_m)+(sideU(component,side_m)-Uin(2,-1:nz+2))
      end if
    else if (Btype(side_m)==BC_NOSLIP.or.(component==2.and.Btype(side_m)==BC_FREESLIP)) then
      if (component==2) then
        Uin(0,-1:nz+2) = 0
        Uin(-1,-1:nz+2)=-Uin(1,-1:nz+2)
      else
        Uin(0,-1:nz+2)=-Uin(1,-1:nz+2)
        Uin(-1,-1:nz+2)=-Uin(2,-1:nz+2)
      end if
    else if (Btype(side_m)==BC_NEUMANN.or.(component/=2.and.Btype(side_m)==BC_FREESLIP)) then
        Uin(0,-1:nz+2) = Uin(1,-1:nz+2)
        Uin(-1,-1:nz+2) = Uin(1,-1:nz+2)
    else if (Btype(side_m)==BC_PERIODIC) then  !Periodic BC
        Uin(0,-1:nz+2) = Uin(ny,-1:nz+2)
        Uin(-1,-1:nz+2) = Uin(ny-1,-1:nz+2)
    end if

    if (Btype(side_p)==BC_DIRICHLET) then
      if (component==2) then
        Uin(ny+1,-1:nz+2) = sideU(component,side_p)
        Uin(ny+2,-1:nz+2) = sideU(component,side_p)+(sideU(component,side_p)-Uin(ny,-1:nz+2))
      else
        Uin(ny+1,-1:nz+2) = sideU(component,side_p)+(sideU(component,side_p)-Uin(ny,-1:nz+2))
        Uin(ny+2,-1:nz+2) = sideU(component,side_p)+(sideU(component,side_p)-Uin(ny-1,-1:nz+2))
      end if
    else if (Btype(side_p)==BC_NOSLIP.or.(component==2.and.Btype(side_m)==BC_FREESLIP)) then
      if (component==2) then
        Uin(ny+1,-1:nz+2) = 0
        Uin(ny+2,-1:nz+2)=-Uin(ny,-1:nz+2)
      else
        Uin(ny+1:ny+2,-1:nz+2)=-Uin(ny:ny-1:-1,-1:nz+2)
      end if
    else if (Btype(side_p)==BC_NEUMANN.or.(component/=2.and.Btype(side_p)==BC_FREESLIP)) then
        Uin(ny+1,-1:nz+2) = Uin(ny,-1:nz+2)
        Uin(ny+2,-1:nz+2) = Uin(ny,-1:nz+2)
    else if (Btype(side_p)==BC_PERIODIC) then  !Periodic BC
        Uin(ny+1:ny+2,-1:nz+2) = Uin(1:2,-1:nz+2)
    end if


  end  subroutine bound_Uin





end  module TurbInlet
