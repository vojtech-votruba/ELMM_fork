module VerticalProfiles

  use Parameters
  use Strings
  use Directories
  
  implicit none

  type Profiles
    real(knd), dimension(:), allocatable :: u, v, uu, vv, ww, &
                                            uusgs, vvsgs, wwsgs, &
                                            tkesgs, dissip, &
                                            temp, tt, &
                                            tempfl, tempflsgs, tempflvisc, &
                                            moist, mm, &
                                            moistfl, moistflsgs, moistflvisc, &
                                            uw, uwsgs, &
                                            vw, vwsgs, &
                                            uz_visc, vz_visc

    real(knd), dimension(:,:), allocatable :: scal, scalfl, &  !which scalar, height
                                              scalflsgs, ss
                                           
  contains
    procedure :: Allocate => Profiles_Allocate
    procedure :: Save => Profiles_Save
    procedure :: Add => Profiles_Add
  end type
  
  type, abstract :: DerivedProfiles
    real(knd) :: start, end
    type(Profiles) :: profs
    character(1024) :: subdir
  contains
    procedure(DerivedProfiles_TimeStep), deferred :: TimeStep
    procedure(DerivedProfiles_Save), deferred :: Save
  end type
  
  type, extends(DerivedProfiles) :: InstantaneousProfiles
    real(knd) :: interval
    integer :: nth = 1
  contains
    procedure :: TimeStep => InstantaneousProfiles_TimeStep
    procedure :: Save => InstantaneousProfiles_Save
  end type
  
  type, extends(DerivedProfiles) :: TimeAveragedProfiles
    logical :: saved = .false.
  contains
    procedure :: TimeStep => TimeAveragedProfiles_TimeStep
    procedure :: Save => TimeAveragedProfiles_Save
    procedure :: ForceSave => TimeAveragedProfiles_ForceSave
  end type
  
  
  type ProfileSwitches
    real(knd) :: average_start = 0, average_end = -1
    real(knd) :: instant_start = 0, instant_end = -1, instant_interval = huge(1.0_knd)
    real(knd) :: running_start = 0, running_end = -1, running_interval = huge(1.0_knd)
  end type
  
  type(Profiles) :: current_profiles
  
  logical :: enable_profiles = .false.
  
  
  real(knd), dimension(:), allocatable :: n_all_Pr, n_free_Pr, n_free_U, n_free_V, n_free_W, &
                                          n_free_PrW, n_free_UW, n_free_VW, n_free_UW_sgs, n_free_VW_sgs
                                            
  real(knd) :: n_free_PrW_surf
  
  abstract interface
    subroutine DerivedProfiles_TimeStep(self, profs, time, dt)
      import
      class(DerivedProfiles), intent(inout) :: self
      type(Profiles), intent(in) :: profs
      real(tim), intent(in) :: time, dt
    end subroutine
    subroutine DerivedProfiles_Save(self)
      import
      class(DerivedProfiles), intent(inout) :: self
    end subroutine
  end interface
  
  interface InstantaneousProfiles
    module procedure InstantaneousProfiles_Init
  end interface

  interface TimeAveragedProfiles
    module procedure TimeAveragedProfiles_Init
  end interface


contains
  
  function InstantaneousProfiles_Init(start, end, interval, subdir) result(res)
    type(InstantaneousProfiles) :: res
    real(knd), intent(in) ::  start, end, interval
    character(*), intent(in) :: subdir
    
    res%start = start
    res%end = end
    res%interval = interval
    res%subdir = subdir
    
    call res%profs%Allocate
    
    call create_directory(output_dir // trim(subdir))
  end function


  subroutine InstantaneousProfiles_TimeStep(self, profs, time, dt)
    class(InstantaneousProfiles), intent(inout) :: self
    type(Profiles), intent(in) :: profs
    real(tim), intent(in) :: time, dt

   if ( (time-dt < self%start + (self%nth-1)*self%interval) .and. &
         (time >= self%start + (self%nth-1)*self%interval ) .and. &
         (time <= self%end)) then
         
      self%profs = profs
      
      call self%Save
      self%nth = self%nth + 1
    end if

  end subroutine
  
  
  subroutine InstantaneousProfiles_Save(self)
    class(InstantaneousProfiles), intent(inout) :: self

    call self%profs%save(output_dir // trim(self%subdir), self%nth)
  end subroutine
  
  
  
  
  
  
  function TimeAveragedProfiles_Init(start, end, subdir) result(res)
    type(TimeAveragedProfiles) :: res
    real(knd), intent(in) ::  start, end
    character(*), intent(in) :: subdir
    
    res%start = start
    res%end = end
    res%subdir = subdir
    
    call res%profs%Allocate
    
    call create_directory(output_dir // trim(subdir))
  end function


  subroutine TimeAveragedProfiles_TimeStep(self, profs, time, dt)
    class(TimeAveragedProfiles), intent(inout) :: self
    type(Profiles), intent(in) :: profs
    real(tim), intent(in) :: time, dt

    if ( (time >= self%start) .and. (time-dt <= self%end) ) then
      
      call self%profs%add(profs, &
                          min(dt, &
                              time - self%start, &
                              self%end - (time-dt), &
                              self%end - self%start) &
                            / (self%end - self%start))
                            
    else if (time > self%end .and. .not. self%saved) then
      call self%Save
    end if
  end subroutine

  
  subroutine TimeAveragedProfiles_Save(self)
    class(TimeAveragedProfiles), intent(inout) :: self

    call self%profs%save(output_dir // trim(self%subdir))
    self%saved = .true.
  end subroutine
  
  subroutine TimeAveragedProfiles_ForceSave(self)
    class(TimeAveragedProfiles), intent(inout) :: self

    if (.not.self%saved) call self%Save
  end subroutine
  
  
  
  
  
  
  
  
  
  
  subroutine Profiles_Allocate(self)
    class(Profiles), intent(out) :: self
  
    allocate(self%u(0:Unz+1), self%v(0:Vnz+1))
    allocate(self%uu(1:Unz), self%vv(1:Vnz), self%ww(0:Prnz))
    allocate(self%uusgs(1:Unz), self%vvsgs(1:Vnz), self%wwsgs(0:Prnz))
    allocate(self%uw(0:Prnz), self%uwsgs(0:Prnz))
    allocate(self%vw(0:Prnz), self%vwsgs(0:Prnz))
    allocate(self%uz_visc(0:Prnz), self%vz_visc(0:Prnz))
    
    allocate(self%tkesgs(1:Prnz), self%dissip(1:Prnz))

    allocate(self%temp(0:Prnz+1), self%tt(1:Prnz))
    allocate(self%tempfl(0:Prnz), self%tempflsgs(0:Prnz), self%tempflvisc(0:Prnz))

    allocate(self%moist(1:Prnz), self%mm(1:Prnz))
    allocate(self%moistfl(0:Prnz), self%moistflsgs(0:Prnz), self%moistflvisc(0:Prnz))

    allocate(self%scal(num_of_scalars,1:Prnz), self%scalfl(num_of_scalars,0:Prnz))
    allocate(self%scalflsgs(num_of_scalars,0:Prnz), self%ss(num_of_scalars,1:Prnz))

    self%u = 0
    self%v = 0
    self%uu = 0
    self%vv = 0
    self%ww = 0
    self%uusgs = 0
    self%vvsgs = 0
    self%wwsgs = 0
    
    self%tkesgs = 0
    self%dissip = 0
    
    self%uw = 0
    self%vw = 0
    self%uwsgs = 0
    self%vwsgs = 0
    self%uz_visc = 0
    self%vz_visc = 0

    self%temp = 0
    self%tempfl = 0
    self%tempflsgs = 0
    self%tempflvisc = 0
    self%tt = 0

    self%moist = 0
    self%moistfl = 0
    self%moistflsgs = 0
    self%moistflvisc = 0
    self%mm = 0

    if (num_of_scalars>0) then
      self%scal = 0
      self%scalfl = 0
      self%scalflsgs = 0
      self%ss = 0
    end if
     
  end subroutine


  
  
  subroutine Profiles_Save(p, dir, nth)
    use FreeUnit
    use custom_par, only: iim, jim, kim, nzims

    class(Profiles), intent(inout) :: p
    character(*), intent(in) :: dir
    integer, intent(in), optional :: nth
    integer :: i, j, k, unit, flux_start
    real(knd), allocatable :: ttsgs(:), mmsgs(:), sssgs(:,:)
    character(:), allocatable :: nth_char
    
    if (present(nth)) then
      nth_char = "-" // itoa(nth)
    else
      nth_char =''
    end if

    call newunit(unit)

    call ReduceProfile(p%u,      n_free_U)
    call ReduceProfile(p%v,      n_free_V)
    
#ifdef PAR
      if (iim==1.and.jim==1) then
        call profile_exchange_bc(p%u)
        call profile_exchange_bc(p%v)
      end if
#endif      
    
    call ReduceProfile(p%uu,     n_free_U(1:Unz))
    call ReduceProfile(p%vv,     n_free_V(1:Vnz))
    call ReduceProfile(p%ww,     n_free_W(0:Prnz))
    call ReduceProfile(p%uusgs,  n_free_U(1:Unz))
    call ReduceProfile(p%vvsgs,  n_free_V(1:Vnz))
    call ReduceProfile(p%wwsgs,  n_free_W(0:Prnz))
    call ReduceProfile(p%tkesgs, n_free_Pr(1:Prnz))
    call ReduceProfile(p%dissip, n_free_Pr(1:Prnz))
    call ReduceProfile(p%uw,     n_free_UW)
    call ReduceProfile(p%vw,     n_free_VW)
    call ReduceProfile(p%uwsgs,  n_free_UW_sgs)
    call ReduceProfile(p%vwsgs,  n_free_VW_sgs)
    call ReduceProfile(p%uz_visc,  n_free_UW_sgs)
    call ReduceProfile(p%vz_visc,  n_free_VW_sgs)
    
    if (enable_buoyancy) then
      call ReduceProfile(p%temp(1:Prnz),      n_free_Pr(1:Prnz))
#ifdef PAR
      if (iim==1.and.jim==1) then
        call profile_exchange_bc(p%temp)
      end if
#endif      
      call ReduceProfile(p%tempfl,    n_all_Pr(0:Prnz))
      call ReduceProfile(p%tempflsgs, [n_free_PrW_surf, n_free_PrW(1:Prnz)])
      call ReduceProfile(p%tempflvisc,n_free_PrW(0:Prnz))
      call ReduceProfile(p%tt,        n_free_Pr(1:Prnz))
    end if
    
    if (enable_moisture) then
      call ReduceProfile(p%moist,      n_free_Pr(1:Prnz))
      call ReduceProfile(p%moistfl,    n_all_Pr(0:Prnz))
      call ReduceProfile(p%moistflsgs, n_free_PrW(0:Prnz))
      call ReduceProfile(p%moistflvisc,n_free_PrW(0:Prnz))
      call ReduceProfile(p%mm,         n_free_Pr(1:Prnz))
    end if
      
    do i = 1, num_of_scalars
      call ReduceProfile(p%scal(i,:),      n_free_Pr(1:Prnz))
      call ReduceProfile(p%scalfl(i,:),    n_all_Pr(0:Prnz))
      call ReduceProfile(p%scalflsgs(i,:), n_free_PrW(0:Prnz))
      call ReduceProfile(p%ss(i,:),        n_free_Pr(1:Prnz))
    end do
      
    if (iim==1.and.jim==1) then
      if (kim > 1) then
        flux_start = 1
      else
        flux_start = 0
      end if

      p%uu = p%uu - p%u(1:Unz)**2
      p%vv = p%vv - p%v(1:Vnz)**2
      
      p%uusgs = p%uusgs + 2 * p%tkesgs / 3
      p%vvsgs = p%vvsgs + 2 * p%tkesgs / 3
      p%wwsgs(1:Prnz-1) = p%wwsgs(1:Prnz-1) + 2 * ((p%tkesgs(1:Prnz-1)+p%tkesgs(2:Prnz))/2) / 3

      open(unit, file = dir // "profu" // nth_char // ".txt")
      do k = 1,Unz !Unz == Vnz
        write(unit,*) zPr(k), p%u(k), p%v(k)
      end do
      close(unit)

      open(unit, file = dir // "profuu" // nth_char // ".txt")
      do k = 1,Unz
        write(unit,*) zPr(k), p%uu(k), p%uusgs(k)
      end do
      close(unit)

      open(unit, file = dir // "profvv" // nth_char // ".txt")
      do k = 1,Vnz
        write(unit,*) zPr(k), p%vv(k), p%vvsgs(k)
      end do
      close(unit)

      open(unit, file = dir // "profww" // nth_char // ".txt")
      do k = 1,Wnz
        write(unit,*) zW(k), p%ww(k), p%wwsgs(k)
      end do
      close(unit)

      open(unit, file = dir // "proftke" // nth_char // ".txt")
      do k = 1,Prnz
        write(unit,*) zPr(k), (p%uu(k)+p%vv(k)+(p%ww(k-1)+p%ww(k))/2) / 2, &
                             p%tkesgs(k)
      end do
      close(unit)

      open(unit, file = dir // "profdissip" // nth_char // ".txt")
      do k = 1,Prnz
        write(unit,*) zPr(k), p%dissip(k)
      end do
      close(unit)

      open(unit, file = dir // "profuw" // nth_char // ".txt")
      do k = flux_start,Prnz
        write(unit,*) zW(k), p%uw(k), p%uwsgs(k), p%uz_visc(k)
      end do
      close(unit)

      open(unit, file = dir // "profvw" // nth_char // ".txt")
      do k = flux_start,Prnz
        write(unit,*) zW(k), p%vw(k), p%vwsgs(k), p%vz_visc(k)
      end do
      close(unit)

      if (enable_buoyancy) then

        open(unit, file = dir // "proftemp" // nth_char // ".txt")
        do k = 1,Prnz
        write(unit,*) zPr(k), p%temp(k)
        end do
        close(unit)

        open(unit, file = dir // "proftempfl" // nth_char // ".txt")
        do k = flux_start,Prnz
        write(unit,*) zW(k), p%tempfl(k), p%tempflsgs(k), p%tempflvisc(k)
        end do
        close(unit)

        p%tt = p%tt - p%temp(1:Prnz)**2
        !Niewstadt, Mason, Moeng, Schumann, 1993 - LES of CBL: A Comparison of Four Computer Codes 
        ttsgs = ((p%tempflsgs(0:Prnz-1)+p%tempflsgs(1:Prnz))/2)**2 / &
                                (0.67_knd**4 * (p%tkesgs(1:Prnz)))
  
        open(unit, file = dir // "proftt" // nth_char // ".txt")
        do k = 1,Prnz
        write(unit,*) zPr(k), p%tt(k), ttsgs(k)
        end do
        close(unit)
      end if

      if (enable_buoyancy) then
        !TODO: This is only valid without moisture. Otherwise we need the virtual temperature (flux).
        block
          real(knd) :: num, denom
          real(knd) :: du_dz, dv_dz, dth_dz
          real(knd) :: Ri, Rf, Prt, Km, Kh
          real(knd) :: shear, N_BV, u_star, theta_star
          real(knd) :: L_MO, L_Oz, L_Ko, l_mix

          open(unit, file = dir // "profother" // nth_char // ".txt")
          write(unit,'(*(a,2x))') "#", "z", "Ri", "Rf", "shear", "N_BV", "u_star", "theta_star", &
                                  "L_MO", "L_Oz", "L_Ko", "Prt", "Km", "Kh", "l_mix"
          do k = 1,Prnz
            if (kim==1 .and. k==1) then
              du_dz = ( p%u(2) - 0 ) / (dzmin)
              dv_dz = ( p%v(2) - 0 ) / (dzmin)
              dth_dz = (p%temp(2)-p%temp(1)) / (dzmin)
            else if (kim==nzims .and. k==Prnz) then
              du_dz = ( p%u(k) - p%u(k-1) ) / (dzmin)
              dv_dz = ( p%v(k) - p%v(k-1) ) / (dzmin)
              dth_dz = ( p%temp(k) - p%temp(k-1)) / (dzmin)
            else
              du_dz = ( p%u(k+1) - p%u(k-1) ) / (2*dzmin)
              dv_dz = ( p%v(k+1) - p%v(k-1) ) / (2*dzmin)
              dth_dz = ( p%temp(k+1) - p%temp(k-1) ) / (2*dzmin)
            end if
            
            shear = hypot(du_dz, dv_dz)
            
            num = (grav_acc/temperature_ref) * dth_dz
            denom = shear**2
            
            if (abs(denom)>1E-5_knd*abs(num) .and. abs(denom)>epsilon(1._knd)) then
              Ri = num / denom
            else
              Ri = 100000._knd * sign(1.0_knd, num) * sign(1.0_knd, denom)
            end if

            
            num = (grav_acc/temperature_ref) * (p%tempfl(k)+p%tempflsgs(k))
            denom = (p%uw(k)+p%uwsgs(k)) * du_dz + (p%vw(k)+p%vwsgs(k)) * dv_dz

            if (abs(denom)>1E-5_knd*abs(num) .and. abs(num)>epsilon(1._knd)) then
              Rf = num / denom
            else
              Rf = 100000._knd * sign(1.0_knd, num) * sign(1.0_knd, denom)
            end if
            

            Km = (hypot(p%uw(k)+p%uwsgs(k), p%vw(k)+p%vwsgs(k))) / shear
            Kh = - (p%tempfl(k)+p%tempflsgs(k)) / dth_dz
            
            Prt = Km / Kh
            
            !allowing NaN for dth_dz<0
            N_BV = sqrt((grav_acc/temperature_ref) * dth_dz)
            
            u_star = sqrt(hypot(p%uw(k)+p%uwsgs(k), p%vw(k)+p%vwsgs(k)))
            
            theta_star = - (p%tempfl(k)+p%tempflsgs(k)) / u_star
            
            !local MO
            L_MO = u_star**2 * temperature_ref / (Karman * grav_acc * theta_star)
            
            L_Oz = sqrt(p%dissip(k) / N_BV**3)
            
            !or eta, the Kolmogorov microscale
            L_Ko = (molecular_viscosity**3 / p%dissip(k))**(1/4._knd)
            
            l_mix = sqrt(Km / shear)
                        
            write(unit,'(*(g0,2x))') zPr(k), Ri, Rf, shear, N_BV, u_star, theta_star, &
                                     L_MO, L_Oz, L_Ko, Prt, Km, Kh, l_mix
          end do
          close(unit)
        end block
      else
        block
          real(knd) :: du_dz, dv_dz
          real(knd) :: Km
          real(knd) :: shear, u_star
          real(knd) :: L_Ko, l_mix
          
          open(unit, file = dir // "profother" // nth_char // ".txt")
          write(unit,'(*(a,2x))') "#", "z", "shear", "u_star", "L_Ko", "Km", "l_mix"
          do k = 1,Prnz
            if (kim==1 .and. k==1) then
              du_dz = ( p%u(2) - 0 ) / (dzmin)
              dv_dz = ( p%v(2) - 0 ) / (dzmin)
            else if (kim==nzims .and. k==Prnz) then
              du_dz = ( p%u(k) - p%u(k-1) ) / (dzmin)
              dv_dz = ( p%v(k) - p%v(k-1) ) / (dzmin)
            else
              du_dz = ( p%u(k+1) - p%u(k-1) ) / (2*dzmin)
              dv_dz = ( p%v(k+1) - p%v(k-1) ) / (2*dzmin)
            end if
            
            shear = hypot(du_dz, dv_dz)
            
            Km = (hypot(p%uw(k)+p%uwsgs(k), p%vw(k)+p%vwsgs(k))) / shear
            
            u_star = sqrt(hypot(p%uw(k)+p%uwsgs(k), p%vw(k)+p%vwsgs(k)))
                        
            !or eta, the Kolmogorov microscale
            L_Ko = (molecular_viscosity**3 / p%dissip(k))**(1/4._knd)
            
            l_mix = sqrt(Km / shear)
                        
            write(unit,'(*(g0,2x))') zPr(k), shear, u_star, L_Ko, Km, l_mix
          end do
          close(unit)
        end block
      end if !enable_buoyancy == 1


      if (enable_moisture) then

         open(unit, file = dir // "profmoist" // nth_char // ".txt")
         do k = 1,Prnz
          write(unit,*) zPr(k), p%moist(k)
         end do
         close(unit)

         open(unit, file = dir // "profmoistfl" // nth_char // ".txt")
         do k = flux_start,Prnz
          write(unit,*) zW(k), p%moistfl(k), p%moistflsgs(k), p%moistflvisc(k)
         end do
         close(unit)

         p%mm = p%mm - p%moist**2
         !see p%ttsgs above
         mmsgs = ((p%moistflsgs(0:Prnz-1)+p%moistflsgs(1:Prnz))/2)**2 / &
                        (0.67_knd**4 * (p%tkesgs(1:Prnz)))
   
         open(unit, file = dir // "profmm" // nth_char // ".txt")
         do k = 1,Prnz
          write(unit,*) zPr(k), p%mm(k), mmsgs(k)
         end do
         close(unit)

      end if !enable_moisture == 1


      if (num_of_scalars>0) then

         open(unit, file = dir // "profscal" // nth_char // ".txt")
         do k = 1,Prnz
          write(unit,*) zPr(k), (p%scal(i,k), i = 1,num_of_scalars)
         end do
         close(unit)

         open(unit, file = dir // "profscalfl" // nth_char // ".txt")
         do k = flux_start,Prnz
          write(unit,*) zW(k), (p%scalfl(i,k), p%scalflsgs(i,k), i = 1,num_of_scalars)
         end do
         close(unit)

         p%ss = p%ss - p%scal**2
         
         allocate(sssgs(num_of_scalars, 1:Prnz))
         !see p%ttsgs above
         do i = 1, num_of_scalars
           sssgs(i,:) = ((p%scalfl(i,0:Prnz-1)+p%scalfl(i,1:Prnz))/2)**2 / &
                               (0.67_knd**4 * (p%tkesgs(1:Prnz)))
         end do
         
         open(unit, file = dir // "profss" // nth_char // ".txt")
         do k = 1,Prnz
          write(unit,*) zPr(k), (p%ss(i,k), sssgs(i,k), i = 1,num_of_scalars)
         end do
         close(unit)

      end if !num_of_scalars>0
     
    end if

#ifdef PAR    
  contains
  
    subroutine profile_exchange_bc(prof)
      use custom_par
      real(knd), intent(inout) :: prof(0:)      
      !set halo cells in profiles for gradient computation
      !TODO: consider fourth order gradients?
      integer :: reqs(4)
      integer :: ie
      
      reqs = MPI_REQUEST_NULL
      
      if (kim<nzims) then
        call MPI_ISend(prof(Prnz), 1, PAR_KND, kim, 1301, comm_row_z, reqs(1), ie)
        call MPI_IRecv(prof(Prnz+1), 1, PAR_KND, kim, 1301, comm_row_z, reqs(2), ie)
      end if
        
      if (kim>1) then
        call MPI_ISend(prof(1), 1, PAR_KND, kim-2, 1301, comm_row_z, reqs(3), ie)
        call MPI_IRecv(prof(0), 1, PAR_KND, kim-2, 1301, comm_row_z, reqs(4), ie)
      end if
        
      call MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE, ie)
    end subroutine
#endif
  end subroutine Profiles_Save
  
  
  subroutine Profiles_Add(self, current, weight)
    class(Profiles), intent(inout) :: self
    type(Profiles), intent(in) :: current
    real(knd), intent(in) :: weight
  
    self%u = self%u + current%u * weight
    self%v = self%v + current%v * weight
    self%uw = self%uw + current%uw * weight
    self%uwsgs = self%uwsgs + current%uwsgs * weight
    self%vw = self%vw + current%vw * weight
    self%vwsgs = self%vwsgs + current%vwsgs * weight
    self%uz_visc = self%uz_visc + current%uz_visc * weight
    self%vz_visc = self%vz_visc + current%vz_visc * weight
    
    self%uu = self%uu + current%uu * weight
    self%vv = self%vv + current%vv * weight
    self%ww = self%ww + current%ww * weight
    self%uusgs = self%uusgs + current%uusgs * weight
    self%vvsgs = self%vvsgs + current%vvsgs * weight
    self%wwsgs = self%wwsgs + current%wwsgs * weight
    self%tkesgs = self%tkesgs + current%tkesgs * weight
    self%dissip = self%dissip + current%dissip * weight

    if (enable_buoyancy) then
      self%temp = self%temp + current%temp * weight
      self%tempfl = self%tempfl + current%tempfl * weight
      self%tempflsgs = self%tempflsgs + current%tempflsgs * weight
      self%tempflvisc = self%tempflvisc + current%tempflvisc * weight
      self%tt = self%tt + current%tt * weight
    end if

    if (enable_buoyancy) then
      self%moist = self%moist + current%moist * weight
      self%moistfl = self%moistfl + current%moistfl * weight
      self%moistflsgs = self%moistflsgs + current%moistflsgs * weight
      self%moistflvisc = self%moistflvisc + current%moistflvisc * weight
      self%mm = self%mm + current%mm * weight
    end if

    if (num_of_scalars>0) then
      self%scal = self%scal + current%scal * weight
      self%scalfl = self%scalfl + current%scalfl * weight
      self%scalflsgs = self%scalflsgs + current%scalflsgs * weight
      self%ss = self%ss + current%ss * weight
    end if
  end subroutine Profiles_Add
  
  
  
  
  subroutine ReduceProfile(prof, n_free)
#ifdef PAR
    use custom_par
#endif
    real(knd), intent(inout) :: prof(:)
    real(knd), intent(in)    :: n_free(:)
  
#ifdef PAR
    call par_sum_to_master_horizontal(prof)
    if (iim==1.and.jim==1) prof = prof / n_free
#else
    prof = prof / n_free
#endif
  end subroutine
  
end module VerticalProfiles
