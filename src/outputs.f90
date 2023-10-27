module Outputs
  use Parameters
  use Boundaries
  use Scalars
  use Wallmodels, only: GroundDeposition, &
                        GroundUstar, GroundUstarUVW, &
                        GroundTFlux, GroundTFluxUVW, &
                        GroundMFlux, GroundMFluxUVW, &
                        wallmodeltype
  use ImmersedBoundary
  use Turbinlet, only: turb_gen => default_turbulence_generator
  use Endianness
  use FreeUnit
  use VTKFrames, only: SaveVTKFrames, FinalizeVTKFrames
#ifdef PAR
  use Frames_ParallelIO, only: SaveFrames_ParallelIO, FinalizeFrames_ParallelIO
#endif
  use SurfaceFrames, only: SaveSurfaceFrames, FinalizeSurfaceFrames
  use StaggeredFrames, only: SaveStaggeredFrames, FinalizeStaggeredFrames
  use Output_helpers
  use VerticalProfiles
  

  implicit none
 

  private
  public store, display, probes, scalar_probes,  &
         OutTStep, Output, CreateOutputDirectories, AllocateOutputs, ReadProbes,  &
         ProfileSwitches, current_profiles, profiles_config, enable_profiles
         
  type(ProfileSwitches) :: profiles_config
  
  type(TimeAveragedProfiles), allocatable :: average_profiles
  type(InstantaneousProfiles), allocatable :: instant_profiles
  type(TimeAveragedProfiles), allocatable :: running_average_profiles(:)

  real(knd), allocatable :: U_avg(:,:,:),V_avg(:,:,:),W_avg(:,:,:) !<u>
  real(knd), allocatable :: UU_prime(:,:,:),VV_prime(:,:,:),WW_prime(:,:,:) !<uu>, <u'u'> must be computed before saving
  real(knd), allocatable :: UV_prime(:,:,:),UW_prime(:,:,:),VW_prime(:,:,:) !<uv>, <u'v'> must be computed before saving
  real(knd), allocatable :: UU_prime_sgs(:,:,:),VV_prime_sgs(:,:,:),WW_prime_sgs(:,:,:) !<uu>, <u'u'> must be computed before saving
  real(knd), allocatable :: UV_prime_sgs(:,:,:),UW_prime_sgs(:,:,:),VW_prime_sgs(:,:,:) !<uv>, <u'v'> must be computed before saving
  real(knd), allocatable :: TKE_prime_sgs(:,:,:) !subgrid TKE approximation
  
  real(knd), allocatable :: Pr_avg(:,:,:) !<p>
  real(knd), allocatable :: Temperature_avg(:,:,:) !<theta>
  real(knd), allocatable :: Moisture_avg(:,:,:) !<q>
  real(knd), allocatable :: Scalar_avg(:,:,:,:) !<c>

  real(knd), allocatable :: Scalar_variance(:,:,:,:)

  real(knd), allocatable :: Scalar_max(:,:,:,:)

  real(knd), allocatable :: Scalar_intermitency(:,:,:,:)

  real(knd), allocatable :: Scalar_fl_U_avg(:,:,:,:) !<cu>, <c'u'> must be computed before saving
  real(knd), allocatable :: Scalar_fl_V_avg(:,:,:,:)
  real(knd), allocatable :: Scalar_fl_W_avg(:,:,:,:)

  real(knd), allocatable :: Scalar_fl_U_sgs(:,:,:,:)
  real(knd), allocatable :: Scalar_fl_V_sgs(:,:,:,:)
  real(knd), allocatable :: Scalar_fl_W_sgs(:,:,:,:)

  real(TIM), allocatable, dimension(:) :: times                                !times of the timesteps

  real(knd), allocatable, dimension(:) :: delta_time, tke, dissip

  real(knd), allocatable, dimension(:) :: pr_gradient_x_time, pr_gradient_y_time, pr_gradient_z_time

  real(knd), allocatable, dimension(:,:) :: ustar, tstar, mstar                !first index differentiates flux from friction number
                                                                             !second index is time
  real(knd), allocatable, dimension(:,:) :: U_time,V_time,W_time,temp_time,moist_time  !position, time

  real(knd), allocatable, dimension(:,:,:) :: scalp_time                        !which scalar, position, time
  real(knd), allocatable, dimension(:,:) :: scalsum_time                        !which scalar, time

  real(knd), allocatable, dimension(:,:,:) :: momentum_fluxes_time, momentum_fluxes_sgs_time   !component, position, time
                                                       !components: 1,1; 1,2; 1,3; 2,2; 2,3; 3,3 for 1..6
  real(knd), allocatable, dimension(:,:,:,:) :: scalar_fluxes_time !component, scalar, position, time
                                                                 !components x,y,z

  type TProbe
    integer :: Ui,Uj,Uk,Vi,Vj,Vk,Wi,Wj,Wk    !grid coordinates of probes in the U,V,W grids
    integer :: i,j,k                         !grid coordinates of probes in the scalar grid
    real(knd) :: x,y,z                       !physical coordinates of probes
    logical :: inside
    integer :: number
  end type TProbe

  !for flow variables including scalar ones (temperature, moisture)
  type(TProbe), allocatable, dimension(:), save :: probes

  !for fluxes
  type(TProbe), allocatable, dimension(:), save :: flux_probes

  !for passive scalars
  type(TProbe), allocatable, dimension(:), save :: scalar_probes
  
  integer :: time_series_max_length = 10000 !how often save the time series
  integer :: time_series_step = 0

  type TOutputSwitches
    integer :: U = 1
    integer :: U_interp = 0
    integer :: V = 1
    integer :: V_interp = 0
    integer :: W = 1
    integer :: W_interp = 0

    integer :: out = 1
    integer :: avg = 1

    integer :: scalars = 1
    integer :: scalars_avg = 1
    integer :: scalars_variance = 0
    integer :: scalars_max = 0
    integer :: scalars_intermitency = 0
    
    real(knd) :: scalars_intermitency_threshold = epsilon(1._knd)

    integer :: deposition = 0


    integer :: out_U = 1
    integer :: out_vorticity = 0
    integer :: out_Pr = 1
    integer :: out_Prtype = 0
    integer :: out_lambda2 = 0
    integer :: out_temperature = 1
    integer :: out_moisture = 1
    integer :: out_divergence = 0
    integer :: out_viscosity = 0

    integer :: avg_U = 1  !1..only in avg.vtk, 2..only in separate Xavg.vtk, 3..both
    integer :: avg_vorticity = 0
    integer :: avg_Pr = 1
    integer :: avg_Prtype = 1
    integer :: avg_temperature = 1
    integer :: avg_moisture = 1

    integer :: avg_UU_prime = 0 !1..only in  avg.vtk, 2..only in separate Xavg.vtk, 3..both
    integer :: avg_flux_scalar = 0
    integer :: avg_flux_scalar_sgs = 0

    integer :: delta_time = 0
    integer :: tke = 0
    integer :: dissip = 0
    integer :: scalsum_time = 1
    integer :: scaltotsum_time = 0
    integer :: ustar = 1
    integer :: tstar = 1
    integer :: mstar = 1
    
    integer :: probes_fluxes = 0
  end type TOutputSwitches

  type(TOutputSwitches), save :: store

  type TDisplaySwitches
    integer :: delta = 0
    integer :: ustar = 0
    integer :: tstar = 0
    integer :: mstar = 0
  end type

  type(TDisplaySwitches), save :: display

  !line feed
  character, parameter :: lf = achar(10)

contains


  subroutine ReadProbes(ps,nps,pfile)
    type(TProbe), allocatable, intent(out) :: ps(:)
    integer, intent(out)     :: nps
    character(*), intent(in) ::pfile
    integer :: i,io,unit
    real(knd) :: tmp3(3)

    open(newunit=unit, file=pfile, status="old", action="read",iostat=io)
    if (io/=0) then
      call error_stop("Error: File "//trim(pfile)//" could not be opened!")
    else
      nps = 0
      do
        read(unit,*,iostat=io) tmp3
        if (io/=0) exit
        nps = nps + 1
      end do

      allocate(ps(nps))

      rewind(unit)
      do i=1,nps
        read(unit,*) ps(i)%x,ps(i)%y,ps(i)%z
      end do
      close(unit)
    end if
  end subroutine
  

  subroutine CreateOutputDirectories
    use custom_par
    
    integer :: j, k, u, io
  
#if defined(_WIN32) || defined(_WIN64)
    call system("mkdir "//output_dir)
#else
    do j = 1, nims
      call par_sync_all()
      
      if (myim==j) then
        do k = 1, 10
          call system("mkdir -p "//output_dir)
          open(newunit=u, file=output_dir//"test", status="replace", iostat=io)
          if (io==0) then
            close(u, status="delete")
            exit
          end if
          call sleep(1)
          
          if (k==10) call error_stop("Error, unable to create "//output_dir)
        end do
      end if
      
      call par_sync_all()
    end do
#endif
  end subroutine


  subroutine AllocateOutputs
    use custom_par
    use Strings
    use Wallmodels
    
    integer :: k

    call GetEndianness

    if (store%avg_U==0.and.store%avg_UU_prime>1) store%avg_U = 1

    call par_sync_out("  ...allocating arrays for 3D statistics.")

    if (averaging==1) then
      if (store%avg_U>0) then
        allocate(U_avg(-2:Unx+3,-2:Uny+3,-2:Unz+3))
        allocate(V_avg(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
        allocate(W_avg(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
        U_avg = 0
        V_avg = 0
        W_avg = 0
      end if

      if (store%avg_UU_prime>0) then
        allocate(UU_prime(-2:Unx+3,-2:Uny+3,-2:Unz+3))
        allocate(VV_prime(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
        allocate(WW_prime(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
        allocate(UV_prime(1:Prnx,1:Prny,1:Prnz))
        allocate(UW_prime(1:Prnx,1:Prny,1:Prnz))
        allocate(VW_prime(1:Prnx,1:Prny,1:Prnz))
        UU_prime = 0
        VV_prime = 0
        WW_prime = 0
        UV_prime = 0
        UW_prime = 0
        VW_prime = 0
        allocate(TKE_prime_sgs(1:Prnx,1:Prny,1:Prnz))
        allocate(UU_prime_sgs(-2:Unx+3,-2:Uny+3,-2:Unz+3))
        allocate(VV_prime_sgs(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
        allocate(WW_prime_sgs(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
        allocate(UV_prime_sgs(1:Prnx,1:Prny,1:Prnz))
        allocate(UW_prime_sgs(1:Prnx,1:Prny,1:Prnz))
        allocate(VW_prime_sgs(1:Prnx,1:Prny,1:Prnz))
        TKE_prime_sgs = 0
        UU_prime_sgs = 0
        VV_prime_sgs = 0
        WW_prime_sgs = 0
        UV_prime_sgs = 0
        UW_prime_sgs = 0
        VW_prime_sgs = 0
      end if

      if (store%avg_Pr==1) then
        allocate(Pr_avg(1:Prnx,1:Prny,1:Prnz))
        Pr_avg = 0
      end if

      if (enable_buoyancy.and.store%avg_temperature==1) then
        allocate(Temperature_avg(-2:Prnx+3,-2:Prny+3,-2:Prnz+3))
        Temperature_avg = 0
      end if

      if (enable_moisture.and.store%avg_moisture==1) then
        allocate(Moisture_avg(-2:Prnx+3,-2:Prny+3,-2:Prnz+3))
        Moisture_avg = 0
      end if
    end if

    if (num_of_scalars>0.and.store%scalars_avg==1) then
      allocate(Scalar_avg(-2:Prnx+3,-2:Prny+3,-2:Prnz+3,num_of_scalars))
      Scalar_avg = 0
    else
      allocate(Scalar_avg(0,0,0,0))
    end if
    
    if (num_of_scalars>0.and.store%scalars_variance==1) then
      allocate(Scalar_variance(-2:Prnx+3,-2:Prny+3,-2:Prnz+3,num_of_scalars))
      Scalar_variance = 0
    else
      allocate(Scalar_variance(0,0,0,0))
    end if

    if (num_of_scalars>0.and.store%scalars_max==1) then
      allocate(Scalar_max(-2:Prnx+3,-2:Prny+3,-2:Prnz+3,num_of_scalars))
      Scalar_max = 0
    else
      allocate(Scalar_max(0,0,0,0))
    end if

    if (num_of_scalars>0.and.store%scalars_intermitency==1) then
      allocate(Scalar_intermitency(-2:Prnx+3,-2:Prny+3,-2:Prnz+3,num_of_scalars))
      Scalar_intermitency = 0
    else
      allocate(Scalar_intermitency(0,0,0,0))
    end if

    if (num_of_scalars>0.and.store%avg_flux_scalar==1) then
      allocate(Scalar_fl_U_avg(Unx,Uny,Unz,num_of_scalars))
      allocate(Scalar_fl_V_avg(Vnx,Vny,Vnz,num_of_scalars))
      allocate(Scalar_fl_W_avg(Wnx,Wny,Wnz,num_of_scalars))
      Scalar_fl_U_avg = 0
      Scalar_fl_V_avg = 0
      Scalar_fl_W_avg = 0
      if (store%avg_flux_scalar_sgs==1) then
        allocate(Scalar_fl_U_sgs(Unx,Uny,Unz,num_of_scalars))
        allocate(Scalar_fl_V_sgs(Vnx,Vny,Vnz,num_of_scalars))
        allocate(Scalar_fl_W_sgs(Wnx,Wny,Wnz,num_of_scalars))
        Scalar_fl_U_sgs = 0
        Scalar_fl_V_sgs = 0
        Scalar_fl_W_sgs = 0
      end if
    end if

    
    call par_sync_out("  ...preparing probes.")

    do k = 1,size(probes)
      associate(p => probes(k))
        p%number = k
        
        p%inside = InDomain(p%x,p%y,p%z)
        
        call GridCoords_interp(p%i,p%j,p%k,p%x,p%y,p%z)

        p%i = max(p%i,1)
        p%j = max(p%j,1)
        p%k = max(p%k,1)
        p%i = min(p%i,Prnx)
        p%j = min(p%j,Prny)
        p%k = min(p%k,Prnz)

        call GridCoords_interp_U(p%Ui,p%Uj,p%Uk,p%x,p%y,p%z)
        call GridCoords_interp_V(p%Vi,p%Vj,p%Vk,p%x,p%y,p%z)
        call GridCoords_interp_W(p%Wi,p%Wj,p%Wk,p%x,p%y,p%z)

        if (Utype(p%Ui, p%Uj, p%Uk)>0) then
          do
            p%Uk = p%Uk + 1
            if (Utype(p%Ui, p%Uj, p%Uk)<=0 .or. p%Uk>=Unz) exit
          end do
        end if

        if (Vtype(p%Vi, p%Vj, p%Vk)>0) then
          do
            p%Vk = p%Vk + 1
            if (Vtype(p%Vi, p%Vj, p%Vk)<=0 .or. p%Vk>=Vnz) exit
          end do
        end if

        if (Wtype(p%Wi, p%Wj, p%Wk)>0) then
          do
            p%Wk = p%Wk + 1
            if (Wtype(p%Wi, p%Wj, p%Wk)<=0 .or. p%Wk>=Wnz) exit
          end do
        end if

      end associate
    end do
    
    do k = 1,size(scalar_probes)
      associate(p => scalar_probes(k))
        p%number = k
        
        p%inside = InDomain(p%x,p%y,p%z)
        
        call GridCoords(p%i,p%j,p%k,p%x,p%y,p%z)

        p%i = max(p%i,1)
        p%j = max(p%j,1)
        p%k = max(p%k,1)
        p%i = min(p%i,Prnx)
        p%j = min(p%j,Prny)
        p%k = min(p%k,Prnz)

        if (Prtype(p%i, p%j, p%k)>0) then
          do
            p%k = p%k + 1
            if (Prtype(p%i, p%j, p%k)<=0 .or. p%k>=Prnz) exit
          end do
        end if
      end associate
    end do


    probes = pack(probes, probes%inside)

    scalar_probes = pack(scalar_probes, scalar_probes%inside)

      
    call par_sync_out("  ...allocating probe time series.")

    if (size(probes)>0) then

      allocate(U_time(size(probes),1:time_series_max_length), &
                V_time(size(probes),1:time_series_max_length), &
                W_time(size(probes),1:time_series_max_length))
      U_time = huge(1.0_knd)
      V_time = huge(1.0_knd)
      W_time = huge(1.0_knd)

      if (store%probes_fluxes==1) then
        allocate(momentum_fluxes_time(6,size(probes),1:time_series_max_length))
        momentum_fluxes_time = huge(1.0_knd)
        allocate(momentum_fluxes_sgs_time(6,size(probes),1:time_series_max_length))
        momentum_fluxes_sgs_time = huge(1.0_knd)
        allocate(scalar_fluxes_time(3,1:num_of_scalars,size(probes),1:time_series_max_length))
        scalar_fluxes_time = huge(1.0_knd)
      else
        allocate(momentum_fluxes_time(0,0,0))
        allocate(momentum_fluxes_sgs_time(0,0,0))
        allocate(scalar_fluxes_time(0,0,0,0))
      end if

      if (enable_buoyancy) then
        allocate(temp_time(size(probes),1:time_series_max_length))
        temp_time = huge(1.0_knd)
      end if

      if (enable_moisture) then
        allocate(moist_time(size(probes),1:time_series_max_length))
        moist_time = huge(1.0_knd)
      end if

    end if

    call par_sync_out("  ...allocating global time series.")

    allocate(times(1:time_series_max_length))
    times = huge(1.0_knd)   

    if (store%delta_time==1) then
      allocate(delta_time(1:time_series_max_length))
      delta_time = huge(1.0_knd)
    end if

    if (store%tke==1) then
      allocate(tke(0:time_series_max_length))
      tke = huge(1.0_knd)
    end if

    if (store%dissip==1) then
      allocate(dissip(1:time_series_max_length))
      dissip = huge(1.0_knd)
      dissip(1) = 0
    end if

    if (wallmodeltype>0.and.(display%ustar==1.or.store%ustar==1)) then
      allocate(ustar(2,1:time_series_max_length))
      ustar = huge(1.0)
    end if

    if (wallmodeltype>0.and.enable_buoyancy.and.(display%tstar==1.or.store%tstar==1)) then
      allocate(tstar(2,1:time_series_max_length))
      tstar = huge(1.0)
    end if

    if (wallmodeltype>0.and.enable_buoyancy.and.(display%mstar==1.or.store%mstar==1)) then
      allocate(mstar(2,1:time_series_max_length))
      mstar = huge(1.0)
    end if

    if (num_of_scalars>0) then
      if (store%scalsum_time==1.or.store%scaltotsum_time==1) then
        allocate(scalsum_time(1:num_of_scalars,1:time_series_max_length))
        scalsum_time = huge(1.0_knd)
      end if

      if (size(scalar_probes)>0) then
        allocate(scalp_time(1:num_of_scalars,1:size(scalar_probes),1:time_series_max_length))
        scalp_time = huge(1.0_knd)
      end if
    end if

    if (flow_rate_x_fixed) then
      allocate(pr_gradient_x_time(1:time_series_max_length))
      pr_gradient_x_time = huge(1.0_knd)
    end if

    if (flow_rate_y_fixed) then
      allocate(pr_gradient_y_time(1:time_series_max_length))
      pr_gradient_y_time = huge(1.0_knd)
    end if

    if (flow_rate_z_fixed) then
      allocate(pr_gradient_z_time(1:time_series_max_length))
      pr_gradient_z_time = huge(1.0_knd)
    end if

    if (enable_profiles) then

      call par_sync_out("  ...preparing profiles.")

      if (.not.allocated(times)) allocate(times(1:time_series_max_length))
    
      call current_profiles%allocate
      
      profiles_config%average_end = min(profiles_config%average_end, time_stepping%end_time)
      profiles_config%running_end = min(profiles_config%running_end, time_stepping%end_time)
      profiles_config%instant_end = min(profiles_config%instant_end, time_stepping%end_time)
      
      profiles_config%average_start = max(profiles_config%average_start, time_stepping%start_time)
      profiles_config%running_start = max(profiles_config%running_start, time_stepping%start_time)
      profiles_config%instant_start = max(profiles_config%instant_start, time_stepping%start_time)
      
      if (profiles_config%average_end > profiles_config%average_start) then
        average_profiles = TimeAveragedProfiles(profiles_config%average_start, &
                                                profiles_config%average_end, &
                                                "average_profiles/")
      end if
    
      if (profiles_config%running_end > profiles_config%running_start .and. &
          profiles_config%running_interval > 0) then
          
        allocate(running_average_profiles( &
                    ceiling( (profiles_config%running_end - profiles_config%running_start) / &
                            profiles_config%running_interval)))
        do k = 1, size(running_average_profiles)
            running_average_profiles(k) = &
              TimeAveragedProfiles(profiles_config%running_start + (k-1)*profiles_config%running_interval, &
                                  profiles_config%running_start + k*profiles_config%running_interval, &
                                  "average_profiles-"//itoa(k)//"/")
        end do
        
      end if
    
      if (profiles_config%instant_end > profiles_config%instant_start) then
        instant_profiles = &
          InstantaneousProfiles(profiles_config%instant_start, &
                                profiles_config%instant_end, &
                                profiles_config%instant_interval, &
                                "instant_profiles/")
      end if
    
      allocate(n_free_U(0:Unz+1))
      allocate(n_free_V(0:Vnz+1))
      allocate(n_free_W(0:Wnz+1))
      allocate(n_free_Pr(0:Prnz+1))
      allocate(n_all_Pr(0:Prnz+1))
      allocate(n_free_PrW(0:Prnz))
      allocate(n_free_UW(0:Prnz))
      allocate(n_free_VW(0:Prnz))
      allocate(n_free_UW_sgs(0:Prnz))
      allocate(n_free_VW_sgs(0:Prnz))

      do k = 0, Unz+1
        n_free_U(k) = count(Utype(1:Unx,1:Uny,k) <= 0)
      end do
      do k = 0, Vnz+1
        n_free_V(k) = count(Vtype(1:Vnx,1:Vny,k) <= 0)
      end do
      do k = 0, Wnz+1
        n_free_W(k) = count(Wtype(1:Wnx,1:Wny,k) <= 0)
      end do
      do k = 0, Prnz+1
        n_free_Pr(k) = count(Prtype(1:Prnx,1:Prny,k) <= 0)
      end do
      do k = 0, Prnz+1
        n_all_Pr(k) = gPrnx * gPrny
      end do
      do k = 0, Prnz
        n_free_PrW(k) = count(Prtype(1:Prnx,1:Prny,k+1) <= 0 .or. Prtype(1:Prnx,1:Prny,k) <= 0)
      end do
      do k = 0, Prnz
        n_free_UW(k) = count((Utype(1:Unx,1:Uny,k+1)<=0.or.Utype(1:Unx,1:Uny,k)<=0) .and. &
                          (Wtype(2:Unx+1,1:Uny,k)<=0.or.Wtype(1:Unx,1:Uny,k)<=0))
      end do
      do k = 0, Prnz
        n_free_VW(k) = count((Vtype(1:Vnx,1:Vny,k+1)<=0.or.Vtype(1:Vnx,1:Vny,k)<=0) .and. &
                          (Wtype(1:Vnx,2:Vny+1,k)<=0.or.Wtype(1:Vnx,1:Vny,k)<=0))
      end do
      do k = 0, Prnz
        n_free_UW_sgs(k) = count((Utype(1:Unx,1:Uny,k+1)<=0.or.Utype(1:Unx,1:Uny,k)<=0)) 
      end do
      do k = 0, Prnz
        n_free_VW_sgs(k) = count((Vtype(1:Vnx,1:Vny,k+1)<=0.or.Vtype(1:Vnx,1:Vny,k)<=0))
      end do
      
      if (kim==1) then
        n_free_PrW_surf = 0
        do k = 1, size(WMPoints)
          if (WMPoints(k)%zk==1.and.Prtype(WMPoints(k)%xi,WMPoints(k)%yj,1)<=0) then
            n_free_PrW_surf = n_free_PrW_surf + 1
          end if
        end do
      endif

#ifdef PAR
      n_free_U = par_co_sum(n_free_U, comm = comm_plane_xy)
      n_free_V = par_co_sum(n_free_V, comm = comm_plane_xy)
      n_free_W = par_co_sum(n_free_W, comm = comm_plane_xy)
      n_free_Pr = par_co_sum(n_free_Pr, comm = comm_plane_xy)
      n_free_PrW = par_co_sum(n_free_PrW, comm = comm_plane_xy)
      n_free_PrW_surf = par_co_sum(n_free_PrW_surf, comm = comm_plane_xy)
      n_free_UW = par_co_sum(n_free_UW, comm = comm_plane_xy)
      n_free_VW = par_co_sum(n_free_VW, comm = comm_plane_xy)
      n_free_UW_sgs = par_co_sum(n_free_UW_sgs, comm = comm_plane_xy)
      n_free_VW_sgs = par_co_sum(n_free_VW_sgs, comm = comm_plane_xy)
#endif

      
    else
    
      if (Btype(To)==BC_AUTOMATIC_FLUX) then
        allocate(current_profiles%uw(0:Prnz), current_profiles%uwsgs(0:Prnz))
        allocate(current_profiles%vw(0:Prnz), current_profiles%vwsgs(0:Prnz))
      end if
      
      if (TempBtype(To)==BC_AUTOMATIC_FLUX) then
        allocate(current_profiles%tempfl(0:Prnz), current_profiles%tempflsgs(0:Prnz))
      else
        !to avoid the necessity of an allocatable dummy argument
        allocate(current_profiles%tempfl(0))
      end if

      if (MoistBtype(To)==BC_AUTOMATIC_FLUX) then
        allocate(current_profiles%moistfl(0:Prnz), current_profiles%moistflsgs(0:Prnz))
      else
        !to avoid the necessity of an allocatable dummy argument
        allocate(current_profiles%moistfl(0))
      end if

    end if

  end subroutine AllocateOutputs





  subroutine InitTimeSeries
    character(5) :: prob
    integer :: k
    
    call create(output_dir//"times.unf")
    
    do k = 1,size(probes)

      write(prob,"(i0)") probes(k)%number

      call create(output_dir//"Utimep"//trim(prob)//".unf")

      call create(output_dir//"Vtimep"//trim(prob)//".unf")

      call create(output_dir//"Wtimep"//trim(prob)//".unf")
      
      if (store%probes_fluxes==1) then
        call create(output_dir//"stresstimep"//trim(prob)//".unf")
        
        call create(output_dir//"sgstresstimep"//trim(prob)//".unf")

        call create(output_dir//"scalfltimep"//trim(prob)//".unf")
      end if

      if (enable_buoyancy) then
        call create(output_dir//"temptimep"//trim(prob)//".unf")
      end if

      if (enable_moisture) then
        call create(output_dir//"moisttimep"//trim(prob)//".unf")
      end if

    end do

    do k = 1,size(scalar_probes)

      write(prob,"(i0)") scalar_probes(k)%number

      if (num_of_scalars>0) then
        call create(output_dir//"scaltimep"//trim(prob)//".unf")
      end if

    end do

    if (store%delta_time==1) then
      call create(output_dir//"delta_time.unf")
    end if

    if (store%tke==1) then
      call create(output_dir//"tke.unf")
    end if

    if (store%tke==1.and.store%dissip==1) then
      call create(output_dir//"dissip.unf")
    end if

    if (wallmodeltype>0.and.store%ustar==1) then
      call create(output_dir//"Retau.unf")
    end if

    if (wallmodeltype>0.and.enable_buoyancy.and.allocated(tstar).and.store%tstar==1) then
      call create(output_dir//"tflux.unf")
    end if

    if (wallmodeltype>0.and.enable_buoyancy.and.allocated(mstar).and.store%tstar==1) then
      call create(output_dir//"mflux.unf")
    end if


    if (flow_rate_x_fixed) then
      call create(output_dir//"pr_gradient_x.unf")
    end if

    if (flow_rate_y_fixed) then
      call create(output_dir//"pr_gradient_y.unf")
    end if

    if (flow_rate_z_fixed) then
      call create(output_dir//"pr_gradient_z.unf")
    end if


    if (num_of_scalars>0.and.store%scalsum_time==1) then
      call create(output_dir//"scalsumtime.unf")
    end if

    if (num_of_scalars>0.and.store%scaltotsum_time==1) then
      call create(output_dir//"scaltotsumtime.unf")
    end if
    
  contains
    subroutine create(fname)
      character(*) :: fname
      integer :: unit, io
      open(newunit=unit,file=fname,status="replace",iostat=io)
      if (io/=0) call error_stop("Error, unable to create "//fname)
      close(unit)
    end subroutine  
  end subroutine InitTimeSeries







  subroutine OutTStep(U,V,W,Pr,Temperature,Moisture,Scalar,dt,delta)
    use Wallmodels, only: ComputeViscsWM
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)   :: U,V,W
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(in)      :: Pr
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)   :: Temperature
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)   :: Moisture
    real(knd), dimension(-2:,-2:,-2:,:), contiguous, intent(in) :: Scalar
    real(knd), intent(in) :: dt
    real(knd), intent(in) :: delta

    integer :: step !just a shorter name
    
    integer :: l,i,j,k
    real(knd) :: S, S2, ground_ustar_Pr, ground_ustar_UVW
    real(knd) :: time_weight
    real(knd) :: fl_L, fl_R
    logical, save :: called = .false.
    
#ifdef CUSTOM_OUTPUT
    interface
      subroutine CustomTimeStepOutput
      end subroutine
    end interface   
#endif

    if (.not.called) then
      call InitTimeSeries
      called = .true.
    end if
    
    time_series_step = time_series_step+1
    step = time_series_step

    ! We compute the subgrid stresses from the eddy viscosity even at the walls
    ! We should not set also TDiff
    if (wallmodeltype>0.and.(store%probes_fluxes==1.or.enable_profiles)) then
      call ComputeViscsWM(U,V,W,Pr,Temperature,Moisture)
    end if

    times(step) = time_stepping%time

    if (store%scalsum_time==1.or.store%scaltotsum_time==1) then
      do l = 1,num_of_scalars
         S = 0
         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
            if (Prtype(i,j,k)<=0) S = S + Scalar(i,j,k,l)*dxPr(i)*dyPr(j)*dzPr(k)
           end do
          end do
         end do
         scalsum_time(l,step) = S
      end do
    end if



    do k = 1,size(probes)
      associate (p=> probes(k))
        U_time(k,step)=Trilinint((p%x-xU(p%Ui))/(xU(p%Ui+1)-xU(p%Ui)), &
                             (p%y-yPr(p%Uj))/(yPr(p%Uj+1)-yPr(p%Uj)), &
                             (p%z-zPr(p%Uk))/(zPr(p%Uk+1)-zPr(p%Uk)), &
                             U(p%Ui,p%Uj,p%Uk),U(p%Ui+1,p%Uj,p%Uk), &
                             U(p%Ui,p%Uj+1,p%Uk),U(p%Ui,p%Uj,p%Uk+1), &
                             U(p%Ui+1,p%Uj+1,p%Uk),U(p%Ui+1,p%Uj,p%Uk+1), &
                             U(p%Ui,p%Uj+1,p%Uk+1),U(p%Ui+1,p%Uj+1,p%Uk+1))

        V_time(k,step)=Trilinint((p%x-xPr(p%Vi))/(xPr(p%Vi+1)-xPr(p%Vi)), &
                             (p%y-yV(p%Vj))/(yV(p%Vj+1)-yV(p%Vj)), &
                             (p%z-zPr(p%Vk))/(zPr(p%Vk+1)-zPr(p%Vk)), &
                             V(p%Vi,p%Vj,p%Vk),V(p%Vi+1,p%Vj,p%Vk), &
                             V(p%Vi,p%Vj+1,p%Vk),V(p%Vi,p%Vj,p%Vk+1), &
                             V(p%Vi+1,p%Vj+1,p%Vk),V(p%Vi+1,p%Vj,p%Vk+1), &
                             V(p%Vi,p%Vj+1,p%Vk+1),V(p%Vi+1,p%Vj+1,p%Vk+1))

        W_time(k,step)=Trilinint((p%x-xPr(p%Wi))/(xPr(p%Wi+1)-xPr(p%Wi)), &
                             (p%y-yPr(p%Wj))/(yPr(p%Wj+1)-yPr(p%Wj)), &
                             (p%z-zW(p%Wk))/(zW(p%Wk+1)-zW(p%Wk)), &
                             W(p%Wi,p%Wj,p%Wk),W(p%Wi+1,p%Wj,p%Wk), &
                             W(p%Wi,p%Wj+1,p%Wk),W(p%Wi,p%Wj,p%Wk+1), &
                             W(p%Wi+1,p%Wj+1,p%Wk),W(p%Wi+1,p%Wj,p%Wk+1), &
                             W(p%Wi,p%Wj+1,p%Wk+1),W(p%Wi+1,p%Wj+1,p%Wk+1))


        if (enable_buoyancy) then
          temp_time(k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)), &
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)), &
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)), &
                             Temperature(p%i,p%j,p%k), &
                             Temperature(p%i+1,p%j,p%k), &
                             Temperature(p%i,p%j+1,p%k), &
                             Temperature(p%i,p%j,p%k+1), &
                             Temperature(p%i+1,p%j+1,p%k), &
                             Temperature(p%i+1,p%j,p%k+1), &
                             Temperature(p%i,p%j+1,p%k+1), &
                             Temperature(p%i+1,p%j+1,p%k+1))
        end if

        if (enable_moisture) then
          moist_time(k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)), &
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)), &
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)), &
                             Moisture(p%i,p%j,p%k), &
                             Moisture(p%i+1,p%j,p%k), &
                             Moisture(p%i,p%j+1,p%k), &
                             Moisture(p%i,p%j,p%k+1), &
                             Moisture(p%i+1,p%j+1,p%k), &
                             Moisture(p%i+1,p%j,p%k+1), &
                             Moisture(p%i,p%j+1,p%k+1), &
                             Moisture(p%i+1,p%j+1,p%k+1))
        end if
        
        if (store%probes_fluxes==1) then 
          !NOTE: valid only for central schemes
          momentum_fluxes_time(1,k,step) = (U_time(k,step))**2
          momentum_fluxes_time(2,k,step) = U_time(k,step) * V_time(k,step)
          momentum_fluxes_time(3,k,step) = U_time(k,step) * W_time(k,step)
          momentum_fluxes_time(4,k,step) = (V_time(k,step))**2
          momentum_fluxes_time(5,k,step) = V_time(k,step) * W_time(k,step)
          momentum_fluxes_time(6,k,step) = (W_time(k,step))**2
          
          !NOTE: subgrid part evaluated in the center of the cell
          momentum_fluxes_sgs_time(1,k,step) = (U(p%i, p%j, p%k) - U(p%i-1, p%j, p%k)) / &
                                             dxPr(p%i)
          momentum_fluxes_sgs_time(4,k,step) = (V(p%i, p%j, p%k) - V(p%i, p%j-1, p%k)) / &
                                             dyPr(p%j)
          momentum_fluxes_sgs_time(6,k,step) = (W(p%i, p%j, p%k) - W(p%i, p%j, p%k-1)) / &
                                             dzPr(p%k)

          fl_L = ( V(p%i+1, p%j, p%k)+V(p%i+1, p%j-1, p%k) - &
                   V(p%i-1, p%j, p%k)+V(p%i-1, p%j-1, p%k) ) / &
                 (2*(xPr(p%i+1)-xPr(p%i-1)))
          fl_R = ( U(p%i, p%j+1, p%k)+U(p%i-1, p%j+1, p%k) - &
                   U(p%i, p%j-1, p%k)+U(p%i-1, p%j-1, p%k) ) / &
                 (2*(yPr(p%j+1)-yPr(p%j-1)))
          momentum_fluxes_sgs_time(2,k,step) = (fl_L + fl_R) / 2
                 
          fl_L = ( W(p%i+1, p%j, p%k)+W(p%i+1, p%j, p%k-1) - &
                   W(p%i-1, p%j, p%k)+W(p%i-1, p%j, p%k-1) ) / &
                 (2*(xPr(p%i+1)-xPr(p%i-1)))
          fl_R = ( U(p%i, p%j, p%k+1)+U(p%i-1, p%j, p%k+1) - &
                   U(p%i, p%j, p%k-1)+U(p%i-1, p%j, p%k-1) ) / &
                 (2*(zPr(p%k+1)-zPr(p%k-1)))
          momentum_fluxes_sgs_time(3,k,step) = (fl_L + fl_R) / 2
                                  
          fl_L = ( W(p%i, p%j+1, p%k)+W(p%i, p%j+1, p%k-1) - &
                   W(p%i, p%j-1, p%k)+W(p%i, p%j-1, p%k-1) ) / &
                 (2*(yPr(p%j+1)-yPr(p%j-1)))
          fl_R = ( V(p%i, p%j, p%k+1)+V(p%i, p%j-1, p%k+1) - &
                   V(p%i, p%j, p%k-1)+V(p%i, p%j-1, p%k-1) ) / &
                 (2*(zPr(p%k+1)-zPr(p%k-1)))
          momentum_fluxes_sgs_time(5,k,step) = (fl_L + fl_R) / 2
                 
          momentum_fluxes_sgs_time(:,k,step) = Viscosity(p%i, p%j, p%k) * &
                                             momentum_fluxes_sgs_time(:,k,step)
        end if
      end associate
    end do

    do k = 1,size(scalar_probes)    
      associate (p=> scalar_probes(k))
        do l = 1,num_of_scalars
          scalp_time(l,k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)), &
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)), &
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)), &
                             Scalar(p%i,p%j,p%k,l), &
                             Scalar(p%i+1,p%j,p%k,l), &
                             Scalar(p%i,p%j+1,p%k,l), &
                             Scalar(p%i,p%j,p%k+1,l), &
                             Scalar(p%i+1,p%j+1,p%k,l), &
                             Scalar(p%i+1,p%j,p%k+1,l), &
                             Scalar(p%i,p%j+1,p%k+1,l), &
                             Scalar(p%i+1,p%j+1,p%k+1,l))
        end do
      end associate
    end do

    if (store%tke==1) then
      tke(step)=TotKE(U,V,W)
    end if

    if (store%tke==1.and.store%dissip==1.and.step>1) then
      dissip(step)=(tke(step-1)-tke(step))/(times(step)-times(step-1))
    end if

    if (store%delta_time==1.and.dt>0) then
      delta_time(step)=delta
    end if

    if (flow_rate_x_fixed) then
      pr_gradient_x_time(step) = pr_gradient_x + pr_gradient_x_dynamic
    end if

    if (flow_rate_y_fixed) then
      pr_gradient_y_time(step) = pr_gradient_y + pr_gradient_y_dynamic
    end if

    if (flow_rate_z_fixed) then
      pr_gradient_z_time(step) = pr_gradient_z + pr_gradient_z_dynamic
    end if



    if (display%delta==1) then
      if (master) write(*,*) "delta: ",delta
    end if




    if ((averaging == 1) .and. &
        ((time_stepping%time >= timeavg1) .and. &
         (time_stepping%time - time_stepping%dt < timeavg2))) then

      time_weight = min(time_stepping%dt, &
                        time_stepping%time - timeavg1, &
                        timeavg2 - (time_stepping%time - time_stepping%dt), &
                        timeavg2 - timeavg1) / &
                    (timeavg2 - timeavg1)

      if (store%avg_U>0) then
         U_avg = U_avg + U * time_weight
         V_avg = V_avg + V * time_weight
         W_avg = W_avg + W * time_weight
      end if

      if (store%avg_UU_prime>0) then
        UU_prime = UU_prime + U**2 * time_weight
        VV_prime = VV_prime + V**2 * time_weight
        WW_prime = WW_prime + W**2 * time_weight
        UV_prime = UV_prime + (U(0:Prnx-1,1:Prny,1:Prnz)+U(1:Prnx,1:Prny,1:Prnz)) * &
                              (V(1:Prnx,0:Prny-1,1:Prnz)+V(1:Prnx,1:Prny,1:Prnz)) * &
                              time_weight / 4
        UW_prime = UW_prime + (U(0:Prnx-1,1:Prny,1:Prnz)+U(1:Prnx,1:Prny,1:Prnz)) * &
                              (W(1:Prnx,1:Prny,0:Prnz-1)+W(1:Prnx,1:Prny,1:Prnz)) * &
                              time_weight / 4
        VW_prime = VW_prime + (V(1:Prnx,0:Prny-1,1:Prnz)+V(1:Prnx,1:Prny,1:Prnz)) * &
                              (W(1:Prnx,1:Prny,0:Prnz-1)+W(1:Prnx,1:Prny,1:Prnz)) * &
                              time_weight / 4
        call AddSubgridStresses(TKE_prime_sgs, UU_prime_sgs, VV_prime_sgs, WW_prime_sgs, &
                                UV_prime_sgs, UW_prime_sgs, VW_prime_sgs, &
                                U, V, W, &
                                time_weight)
      end if

      if (store%avg_Pr==1) then
        Pr_avg = Pr_avg + Pr(1:Prnx,1:Prny,1:Prnz) * time_weight
      end if

      if (enable_buoyancy.and.store%avg_temperature==1) then
        Temperature_avg = Temperature_avg + temperature * time_weight
      end if

      if (enable_moisture.and.store%avg_moisture==1) then
        Moisture_avg = Moisture_avg + moisture * time_weight
      end if

      if (num_of_scalars>0.and.store%scalars_avg==1) then
        Scalar_avg = Scalar_avg + Scalar * time_weight
      end if

      if (num_of_scalars>0.and.store%scalars_variance==1) then
        Scalar_variance = Scalar_variance + Scalar**2 * time_weight
      end if

      if (num_of_scalars>0.and.store%scalars_intermitency==1) then
        where (Scalar>=store%scalars_intermitency_threshold)  &
          Scalar_intermitency = Scalar_intermitency + time_weight
      end if
    end if

    if ((averaging == 1) .and. &
        ((time_stepping%time >= timeavg1) .and. &
         (time_stepping%time - time_stepping%dt < timeavg2))) then

      if (num_of_scalars > 0 .and. store%avg_flux_scalar == 1) then
        do i=1,num_of_scalars
          call AddScalarAdvVector(Scalar_fl_U_avg(:,:,:,i), &
                               Scalar_fl_V_avg(:,:,:,i), &
                               Scalar_fl_W_avg(:,:,:,i), &
                               Scalar(:,:,:,i), &
                               U,V,W,time_weight, &
                               scalar_fluxes_time(:,i,:,step), &
                               scalar_probes%i, scalar_probes%j, scalar_probes%k)
          if (store%avg_flux_scalar_sgs==1) then
            call AddScalarDiffVector(Scalar_fl_U_sgs(:,:,:,i), &
                                  Scalar_fl_V_sgs(:,:,:,i), &
                                  Scalar_fl_W_sgs(:,:,:,i), &
                                  Scalar(:,:,:,i),time_weight, &
                                  scalar_fluxes_time(:,i,:,step), &
                                  scalar_probes%i, scalar_probes%j, scalar_probes%k)
          else
            call AddScalarDiffVector(Scalar_fl_U_avg(:,:,:,i), &
                                  Scalar_fl_V_avg(:,:,:,i), &
                                  Scalar_fl_W_avg(:,:,:,i), &
                                  Scalar(:,:,:,i),time_weight, &
                                  scalar_fluxes_time(:,i,:,step), &
                                  scalar_probes%i, scalar_probes%j, scalar_probes%k)
          end if
        end do
      end if
      
    end if


    if (wallmodeltype>0.and.(display%ustar==1.or.store%ustar==1)) then

      ground_ustar_Pr = GroundUstar()

      if (molecular_viscosity > 0) then
        S2 = ground_ustar_Pr / molecular_viscosity
      else
        S2 = 0
      end if

      if (display%ustar==1) then
        if (allocated(turb_gen%Ustar_inlet)) then
         if (master) write(*,'(*(a10,es15.3))') "ustar:",ground_ustar_Pr,"Re_tau:",S2,"u*inlet",turb_gen%Ustar_inlet(1)
        else
         if (master) write(*,'(*(a10,es15.3))') "ustar:",ground_ustar_Pr,"Re_tau:",S2
        end if
      end if

      ground_ustar_UVW = GroundUstarUVW()

      if (molecular_viscosity > 0) then
        S2 = ground_ustar_UVW / molecular_viscosity
      else
        S2 = 0
      end if

      if (display%ustar==1) then
        if (allocated(turb_gen%Ustar_inlet)) then
         if (master) write(*,'(*(a10,es15.3))') "ustarUVW:",ground_ustar_UVW,"Re_tau:",S2,"u*inlet",turb_gen%Ustar_inlet(1)
        else
         if (master) write(*,'(*(a10,es15.3))') "ustarUVW:",ground_ustar_UVW,"Re_tau:",S2
        end if
      end if
      if (store%ustar==1) then
        ustar(:,step)=[ S2 , ground_ustar_UVW ]
      end if
      
    end if


    if (wallmodeltype>0.and.enable_buoyancy.and.(display%tstar==1.or.store%tstar==1)) then
     S2 = GroundTFlux()
     S = - S2 /  ground_ustar_Pr

     if (display%tstar==1) then
       if (master) write(*,'(*(a10,es15.3))') "Tstar:", S, "Tflux:", S2
     end if

     if (store%tstar==1) then
       tstar(:,step) = [ S2,S ]
     end if

     S2 = GroundTFluxUVW()
     S = - S2 / ground_ustar_UVW 

     if (display%tstar==1) then
       if (master) write(*,'(*(a10,es15.3))') "TstarUVW:" ,S, "Tflux:", S2
     end if

    end if

    if (wallmodeltype>0.and.enable_moisture.and.(display%mstar==1.or.store%mstar==1)) then
     S2 = GroundMFlux()
     S = - S2 /  ground_ustar_Pr

     if (display%mstar==1) then
       if (master) write(*,'(*(a10,es15.3))') "Mstar:", S, "Mflux:", S2
     end if

     if (store%mstar==1) then
       mstar(:,step) = [ S2,S ]
     end if

     S2 = GroundMFluxUVW()
     S = - S2 / ground_ustar_UVW 

     if (display%mstar==1) then
       if (master) write(*,'(*(a10,es15.3))') "MstarUVW:" ,S, "Mflux:", S2
     end if

    end if

    if ((enable_profiles) .or. &
       Btype(To)==BC_AUTOMATIC_FLUX) then
      call StressProfiles(U,V,W)
    end if

    if (enable_profiles) then

      call FluxSGSProfiles(W,Temperature,Moisture,Scalar)

    else if (TempBtype(To)==BC_AUTOMATIC_FLUX .or. &
             MoistBtype(To)==BC_AUTOMATIC_FLUX .or. &
             ScalBType(To)==BC_AUTOMATIC_FLUX) then

      if (TempBtype(To)==BC_AUTOMATIC_FLUX.and.enable_buoyancy) &
        call TemperatureFluxSGSProfile(W,Temperature)

      if (MoistBtype(To)==BC_AUTOMATIC_FLUX.and.enable_moisture) &
        call MoistureFluxSGSProfile(W,Moisture)

      if (ScalBtype(To)==BC_AUTOMATIC_FLUX.and.num_of_scalars>0) &
        call ScalarFluxSGSProfile(W,Scalar)

    end if

    if (enable_profiles) then
      call BLProfiles(U,V,W,Temperature,Moisture,Scalar)
      
      if (allocated(average_profiles)) &
        call average_profiles%TimeStep(current_profiles, time_stepping%time, time_stepping%dt)

      if (allocated(instant_profiles)) &
        call instant_profiles%TimeStep(current_profiles, time_stepping%time, time_stepping%dt)
      
      if (allocated(running_average_profiles)) then
        do i = 1, size(running_average_profiles)
          call running_average_profiles(i)%TimeStep(current_profiles, time_stepping%time, time_stepping%dt)
        end do
      end if
    end if


    call SaveVTKFrames(time_stepping%time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    
#ifdef PAR
    call SaveFrames_ParallelIO(time_stepping%time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
#endif

    call SaveSurfaceFrames(time_stepping%time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)

    call SaveStaggeredFrames(time_stepping%time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    
    if (time_series_step==time_series_max_length) then
      call OutputTimeSeries
      time_series_step = 0
    end if
    
#ifdef CUSTOM_OUTPUT
    call CustomTimeStepOutput
#endif

  end subroutine OutTstep


  subroutine AddSubgridStresses(TKE, UU, VV, WW, UV, UW, VW, &
                                U, V, W, weight)
    use Filters, only: filtertype, filter_ratios
    real(knd), dimension( 1:, 1:, 1:), contiguous, intent(inout) :: TKE
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: UU, VV, WW
    real(knd), dimension( 1:, 1:, 1:), contiguous, intent(inout) :: UV, UW, VW
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U, V, W
    real(knd), intent(in) :: weight
    real(knd) :: Ax, Ay, Az, Ax2, Ay2, Az2
    integer :: i, j, k
    real(knd) :: width
    real(knd), parameter :: Ck = 0.1 !model constant for the sgs_tke eddy viscosity sgs model

    width = filter_ratios(filtertype) * (dxmin*dymin*dzmin)**(1._knd/3._knd)
    
    Ax = weight/(2*dxmin)
    Ay = weight/(2*dymin)
    Az = weight/(2*dzmin)

    Ax2 = weight/(8*dxmin)
    Ay2 = weight/(8*dymin)
    Az2 = weight/(8*dzmin)

    !NOTE:neglecting part caused by molecular viscosity
    !$omp parallel private(i,j,k)
    !$omp do
    do k=1,Unz
      do j=1,Uny
        do i=1,Unx
          UU(i,j,k) =  UU(i,j,k) + &
                        Ax * (  (Viscosity(i+1,j,k) - molecular_viscosity) &
                                 * (U(i+1,j,k)-U(i,j,k)) &
                              + (Viscosity(i,j,k) - molecular_viscosity) &
                                 * (U(i,j,k)-U(i-1,j,k)))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k=1,Vnz
      do j=1,Vny
        do i=1,Vnx
          VV(i,j,k) = VV(i,j,k) + &
                       Ay * (  (Viscosity(i,j+1,k) - molecular_viscosity) &
                                * (V(i,j+1,k)-V(i,j,k)) &
                             + (Viscosity(i,j,k) - molecular_viscosity) &
                                * (V(i,j,k)-V(i,j-1,k)))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k=1,Wnz
      do j=1,Wny
        do i=1,Wnx
          WW(i,j,k) = WW(i,j,k) + &
                       Az * (  (Viscosity(i,j,k+1) - molecular_viscosity) &
                                * (W(i,j,k+1)-W(i,j,k)) &
                             + (Viscosity(i,j,k) - molecular_viscosity) &
                                * (W(i,j,k)-W(i,j,k-1)))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k=1,Prnz
      do j=1,Prny
        do i=1,Prnx
          UV(i,j,k) = UV(i,j,k) + (Viscosity(i,j,k) - molecular_viscosity) * &
                          (Ay2 * (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k)) + &
                           Ax2 * (V(i+1,j,k)-V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k)))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k=1,Prnz
      do j=1,Prny
        do i=1,Prnx
          UW(i,j,k) = UW(i,j,k) + (Viscosity(i,j,k) - molecular_viscosity) * &
                          (Az2 * (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) + &
                           Ax2 * (W(i+1,j,k)-W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1)))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k=1,Prnz
      do j=1,Prny
        do i=1,Prnx
          VW(i,j,k) = VW(i,j,k) + (Viscosity(i,j,k) - molecular_viscosity) * &
                          (Az2 * (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) + &
                           Ay2 * (W(i,j+1,k)-W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1)))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k=1,Prnz
      do j=1,Prny
        do i=1,Prnx
          TKE(i,j,k) = TKE(i,j,k) + weight * ((Viscosity(i,j,k) - molecular_viscosity) / (Ck * width))**2
        end do
      end do
    end do
    !$omp end do nowait
    !$omp end parallel
  end subroutine AddSubgridStresses


  subroutine OutputTimeSeries
    character(5) :: prob
    integer :: k,unit

    call newunit(unit)
    
    open(unit,file=output_dir//"times.unf", &
         access="stream",status="old",position="append")
    write(unit) times(1:time_series_step)
    close(unit)

    do k = 1,size(probes)

      write(prob,"(i0)") probes(k)%number

      open(unit,file=output_dir//"Utimep"//trim(prob)//".unf", &
           access="stream",status="old",position="append")
      write(unit) U_time(k,1:time_series_step)
      close(unit)

      open(unit,file=output_dir//"Vtimep"//trim(prob)//".unf", &
           access="stream",status="old",position="append")
      write(unit) V_time(k,1:time_series_step)
      close(unit)

      open(unit,file=output_dir//"Wtimep"//trim(prob)//".unf", &
           access="stream",status="old",position="append")
      write(unit) W_time(k,1:time_series_step)
      close(unit)
      
      if (store%probes_fluxes==1) then
        open(unit,file=output_dir//"stresstimep"//trim(prob)//".unf", &
             access="stream",status="old",position="append")
        write(unit) momentum_fluxes_time(:,k,1:time_series_step)
        close(unit)
        
        open(unit,file=output_dir//"sgstresstimep"//trim(prob)//".unf", &
             access="stream",status="old",position="append")
        write(unit) momentum_fluxes_time(:,k,1:time_series_step)
        close(unit)

        open(unit,file=output_dir//"scalfltimep"//trim(prob)//".unf", &
             access="stream",status="old",position="append")
        write(unit) scalar_fluxes_time(:,:,k,1:time_series_step)
        close(unit)
      end if

      if (enable_buoyancy) then
        open(unit,file=output_dir//"temptimep"//trim(prob)//".unf", &
             access="stream",status="old",position="append")
        write(unit) temp_time(k,1:time_series_step)
        close(unit)
      end if

      if (enable_moisture) then
        open(unit,file=output_dir//"moisttimep"//trim(prob)//".unf", &
             access="stream",status="old",position="append")
        write(unit) moist_time(k,1:time_series_step)
        close(unit)
      end if

    end do

    do k = 1,size(scalar_probes)

      write(prob,"(i0)") scalar_probes(k)%number

      if (num_of_scalars>0) then
        open(unit,file=output_dir//"scaltimep"//trim(prob)//".unf", &
             access="stream",status="old",position="append")
        write(unit) scalp_time(:,k,1:time_series_step)
        close(unit)
      end if

    end do

    if (store%delta_time==1) then
      open(unit,file=output_dir//"delta_time.unf", &
           access="stream",status="old",position="append")
      write(unit) delta_time(1:time_series_step)
      close(unit)
    end if

    if (store%tke==1) then
      open(unit,file=output_dir//"tke.unf", &
           access="stream",status="old",position="append")
      write(unit) tke(1:time_series_step)
      close(unit)
    end if

    if (store%tke==1.and.store%dissip==1) then
      open(unit,file=output_dir//"dissip.unf", &
           access="stream",status="old",position="append")
      write(unit) dissip(1:time_series_step)
      close(unit)
    end if

    if (wallmodeltype>0.and.store%ustar==1) then
      open(unit,file=output_dir//"Retau.unf", &
           access="stream",status="old",position="append")
      write(unit) ustar(:,1:time_series_step)
      close(unit)
    end if

    if (wallmodeltype>0.and.enable_buoyancy.and.allocated(tstar).and.store%tstar==1) then
      open(unit,file=output_dir//"tflux.unf", &
           access="stream",status="old",position="append")
      write(unit) tstar(:,1:time_series_step)
      close(unit)
    end if

    if (wallmodeltype>0.and.enable_moisture.and.allocated(mstar).and.store%mstar==1) then
      open(unit,file=output_dir//"mflux.unf", &
           access="stream",status="old",position="append")
      write(unit) mstar(:,1:time_series_step)
      close(unit)
    end if

    if (flow_rate_x_fixed) then
      open(unit,file=output_dir//"pr_gradient_x.unf", &
           access="stream",status="old",position="append")
      write(unit) pr_gradient_x_time(1:time_series_step)
      close(unit)
    end if

    if (flow_rate_y_fixed) then
      open(unit,file=output_dir//"pr_gradient_y.unf", &
           access="stream",status="old",position="append")
      write(unit) pr_gradient_y_time(1:time_series_step)
      close(unit)
    end if

    if (flow_rate_z_fixed) then
      open(unit,file=output_dir//"pr_gradient_z.unf", &
           access="stream",status="old",position="append")
      write(unit) pr_gradient_z_time(1:time_series_step)
      close(unit)
    end if

    if (num_of_scalars>0.and.store%scalsum_time==1) then
      open(unit,file=output_dir//"scalsumtime.unf", &
           access="stream",status="old",position="append")
      write(unit) scalsum_time(:,1:time_series_step)
      close(unit)
    end if

    if (num_of_scalars>0.and.store%scaltotsum_time==1) then
      open(unit,file=output_dir//"scaltotsumtime.unf", &
           access="stream",status="old",position="append")
      write(unit) sum(scalsum_time(:,1:time_series_step), dim=1)
      close(unit)
    end if
    
  end subroutine OutputTimeSeries
  
  
  

  subroutine OutputProfiles
    integer :: i


    if (enable_profiles) then

      if (allocated(average_profiles)) call average_profiles%ForceSave
      
      if (allocated(running_average_profiles)) then
        do i = 1, size(running_average_profiles)
          call running_average_profiles(i)%ForceSave
        end do
      end if

    end if

  end subroutine OutputProfiles
  
  






  subroutine OutputOut(U,V,W,Pr,Temperature,Moisture)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(in) :: Pr
    real(knd), contiguous, intent(in) :: Temperature(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: Moisture(-2:,-2:,-2:)
    character(70) :: str
    real(real32), allocatable :: tmp(:,:,:,:)
    integer :: i,j,k,unit
    real(knd), parameter :: C1 = 9._knd / 8, C3 = 1._knd / (8*3)

    if (store%out==1) then

       call newunit(unit);

       open(unit,file=output_dir//"out.vtk", &
         access='stream',status='replace',form="unformatted",action="write")

       write(unit) "# vtk DataFile Version 2.0", lf
       write(unit) "CLMM output file", lf
       write(unit) "BINARY", lf
       write(unit) "DATASET RECTILINEAR_GRID", lf
       str="DIMENSIONS"
       write(str(12:),*) Prnx,Prny,Prnz
       write(unit) str, lf
       str="X_COORDINATES"
       write(str(15:),'(i5,2x,a)') Prnx,"float"
       write(unit) str, lf
       write(unit) BigEnd(real(xPr(1:Prnx), real32)), lf
       str="Y_COORDINATES"
       write(str(15:),'(i5,2x,a)') Prny,"float"
       write(unit) str, lf
       write(unit) BigEnd(real(yPr(1:Prny), real32)), lf
       str="Z_COORDINATES"
       write(str(15:),'(i5,2x,a)') Prnz,"float"
       write(unit) str, lf
       write(unit) BigEnd(real(zPr(1:Prnz), real32)), lf
       str="POINT_DATA"
       write(str(12:),*) Prnx*Prny*Prnz
       write(unit) str, lf

       if (store%out_Pr==1) then
         write(unit) "SCALARS p float", lf
         write(unit) "LOOKUP_TABLE default", lf

         write(unit) BigEnd(real(Pr(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (enable_buoyancy.and.store%out_temperature==1) then
         write(unit) "SCALARS temperature float", lf
         write(unit) "LOOKUP_TABLE default", lf

         write(unit) BigEnd(real(Temperature(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (enable_moisture.and.store%out_moisture==1) then
         write(unit) "SCALARS moisture float", lf
         write(unit) "LOOKUP_TABLE default", lf

         write(unit) BigEnd(real(Moisture(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (store%out_Prtype==1) then
         write(unit) "SCALARS ptype float", lf
         write(unit) "LOOKUP_TABLE default", lf

         write(unit) BigEnd(real(Prtype(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (store%out_divergence==1) then
         write(unit) "SCALARS div float", lf
         write(unit) "LOOKUP_TABLE default", lf

         if (gridtype==GRID_UNIFORM) then
           if (discretization_order == 4) then
             write(unit) BigEnd(real( &
                     ( C1*(U(1:Prnx,1:Prny,1:Prnz)-U(0:Prnx-1,1:Prny,1:Prnz)) - &
                       C3*(U(2:Prnx+1,1:Prny,1:Prnz)-U(-1:Prnx-2,1:Prny,1:Prnz)) ) / dxmin &
                   + ( C1*(V(1:Prnx,1:Prny,1:Prnz)-V(1:Prnx,0:Prny-1,1:Prnz)) - &
                       C3*(V(1:Prnx,2:Prny+1,1:Prnz)-V(1:Prnx,-1:Prny-2,1:Prnz)) ) / dymin &
                   + ( C1*(W(1:Prnx,1:Prny,1:Prnz)-W(1:Prnx,1:Prny,0:Prnz-1)) - &
                       C3*(W(1:Prnx,1:Prny,2:Prnz+1)-W(1:Prnx,1:Prny,-1:Prnz-2)) ) / dzmin &
                                , real32))
           else
             write(unit) BigEnd(real((U(1:Prnx,1:Prny,1:Prnz) - &
                                      U(0:Prnx-1,1:Prny,1:Prnz))/dxmin + &
                                     (V(1:Prnx,1:Prny,1:Prnz) - &
                                      V(1:Prnx,0:Prny-1,1:Prnz))/dymin + &
                                     (W(1:Prnx,1:Prny,1:Prnz) - &
                                      W(1:Prnx,1:Prny,0:Prnz-1))/dzmin &
                                , real32))
           end if
         else
           do k = 1,Prnz
            do j = 1,Prny
             do i = 1,Prnx
               write(unit) BigEnd(real((U(i,j,k)-U(i-1,j,k))/(dxPr(i)) + &
                                       (V(i,j,k)-V(i,j-1,k))/(dyPr(j)) + &
                                       (W(i,j,k)-W(i,j,k-1))/(dzPr(k)), real32))
             end do
            end do
           end do
         end if

         write(unit) lf
       end if

       if (store%out_lambda2==1) then
         allocate(tmp(1:Prnx,1:Prny,1:Prnz,1))

         write(unit) "SCALARS lambda2 float", lf
         write(unit) "LOOKUP_TABLE default", lf

         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
             tmp(i,j,k,1) = BigEnd(real(Lambda2(i,j,k,U,V,W), real32))
           end do
          end do
         end do

         write(unit) tmp

         write(unit) lf
         deallocate(tmp)
       end if

       if (store%out_viscosity==1) then
         write(unit) "SCALARS visc float", lf
         write(unit) "LOOKUP_TABLE default", lf

         write(unit) BigEnd(real(Viscosity(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (store%out_U==1.or.store%out_vorticity==1) allocate(tmp(1:3,1:Prnx,1:Prny,1:Prnz))

       if (store%out_U==1) then
         write(unit) "VECTORS u float", lf

         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
             tmp(1:3,i,j,k) = [ BigEnd(real((U(i,j,k)+U(i-1,j,k))/2._knd, real32)), &
                                BigEnd(real((V(i,j,k)+V(i,j-1,k))/2._knd, real32)), &
                                BigEnd(real((W(i,j,k)+W(i,j,k-1))/2._knd, real32)) ]
           end do
          end do
         end do

         write(unit) tmp

         write(unit) lf
       end if

       if (store%out_vorticity==1) then
         write(unit) "VECTORS vort float", lf

         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
             tmp(1:3,i,j,k) = [ BigEnd(real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin) &
                                          -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin), real32)), &
                                BigEnd(real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin) &
                                          -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin), real32)), &
                                BigEnd(real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin) &
                                          -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin), real32)) ]
           end do
          end do
         end do

         write(unit) tmp

         write(unit) lf
       end if
       close(unit)
    end if  !store%out
  end subroutine OutputOut








  subroutine OutputScalars(Scalar)
    real(knd), dimension(-2:,-2:,-2:,:), contiguous, intent(in) :: Scalar
    character(70) :: str
    real(knd), dimension(:,:,:), allocatable :: depos
    character(8) ::  scalname
    integer :: l,unit

    scalname = "scalar"

    if (num_of_scalars>0) then
      if (store%scalars==1) then

          call newunit(unit)

          open(unit,file=output_dir//"scalars.vtk", &
            access='stream',status='replace',form="unformatted",action="write")

          write(unit) "# vtk DataFile Version 2.0", lf
          write(unit) "CLMM output file", lf
          write(unit) "BINARY", lf
          write(unit) "DATASET RECTILINEAR_GRID", lf
          str="DIMENSIONS"
          write(str(12:),*) Prnx,Prny,Prnz
          write(unit) str, lf
          str="X_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prnx,"float"
          write(unit) str, lf
          write(unit) BigEnd(real(xPr(1:Prnx), real32)), lf
          str="Y_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prny,"float"
          write(unit) str, lf
          write(unit) BigEnd(real(yPr(1:Prny), real32)), lf
          str="Z_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prnz,"float"
          write(unit) str, lf
          write(unit) BigEnd(real(zPr(1:Prnz), real32)), lf
          str="POINT_DATA"
          write(str(12:),*) Prnx*Prny*Prnz
          write(unit) str, lf

          do l = 1,num_of_scalars
            write(scalname(7:8),"(I2.2)") l
            write(unit) "SCALARS ", scalname , " float", lf
            write(unit) "LOOKUP_TABLE default", lf

            write(unit) BigEnd(real(Scalar(1:Prnx,1:Prny,1:Prnz,l), real32))

            write(unit) lf
          end do
          close(unit)
      end if !store%scalars

      if (computedeposition>0.and.store%deposition==1) then

          allocate(depos(1:Prnx,1:Prny,num_of_scalars))
          depos = GroundDeposition()

          call newunit(unit)

          open(unit,file=output_dir//"deposition.vtk", &
            access='stream',status='replace',form="unformatted",action="write")

          write(unit) "CLMM output file", lf
          write(unit) "BINARY", lf
          write(unit) "DATASET RECTILINEAR_GRID", lf
          str="DIMENSIONS"
          write(str(12:),*) Prnx,Prny,1
          write(unit) str, lf
          str="X_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prnx,"float"
          write(unit) str, lf
          write(unit) BigEnd(real(xPr(1:Prnx), real32)), lf
          str="Y_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prny,"float"
          write(unit) str, lf
          write(unit) BigEnd(real(yPr(1:Prny), real32))
          str="Z_COORDINATES"
          write(str(15:),'(i5,2x,a)') 1,"float"
          write(unit) str, lf
          write(unit) BigEnd(real(zW(0), real32)), lf
          str="POINT_DATA"
          write(str(12:),*) Prnx*Prny
          write(unit) str, lf

          do l = 1,num_of_scalars
            write(scalname(7:8),"(I2.2)") l
            write(unit) "SCALARS ", scalname , " float", lf
            write(unit) "LOOKUP_TABLE default", lf

            write(unit) BigEnd(real(depos(1:Prnx,1:Prny,l), real32))

            write(unit) lf
          end do

          close(unit)

          deallocate(depos)

      end if  !store%deposition
    end if  !num_of_scalars
  end subroutine OutputScalars










  subroutine OutputAvg(U, V, W, &
                       Pr, Temperature, Moisture, &
                       UU, VV, WW, UV, UW, VW, &
                       TKE_sgs, UU_sgs, VV_sgs, WW_sgs, UV_sgs, UW_sgs, VW_sgs)
    real(knd), dimension(:,:,:), allocatable, intent(inout) :: U, V, W
    real(knd), dimension(:,:,:), allocatable, intent(inout) :: UU, VV, WW
    real(knd), dimension(:,:,:), allocatable, intent(inout) :: UV, UW, VW
    real(knd), dimension(:,:,:), allocatable, intent(inout) :: TKE_sgs
    real(knd), dimension(:,:,:), allocatable, intent(inout) :: UU_sgs, VV_sgs, WW_sgs
    real(knd), dimension(:,:,:), allocatable, intent(inout) :: UV_sgs, UW_sgs, VW_sgs
    real(knd), dimension(:,:,:), allocatable, intent(inout) :: Pr
    real(knd), dimension(:,:,:), allocatable, intent(inout) :: Temperature
    real(knd), dimension(:,:,:), allocatable, intent(inout) :: Moisture
    character(70) :: str
    integer :: i,j,k,unit
    real(real32), allocatable :: tmp(:,:,:,:), sc_tmp(:,:,:)
    real(knd) :: time_factor

    if (averaging==1 .and. time_stepping%time > timeavg1) then
      if (time_stepping%time < timeavg2) then
        time_factor = (timeavg2 - timeavg1) / (time_stepping%time - timeavg1)
        
        U = U * time_factor
        V = V * time_factor
        W = W * time_factor
        
        if (store%avg_Pr==1) then
          Pr = Pr * time_factor
        end if

        if (enable_buoyancy.and.store%avg_temperature==1) then
          Temperature = Temperature * time_factor
        end if

        if (enable_moisture.and.store%avg_moisture==1) then
          Moisture = Moisture * time_factor
        end if                
        
        if (store%avg_UU_prime>0) then
          UU = UU * time_factor
          VV = VV * time_factor
          WW = WW * time_factor
          UV = UV * time_factor
          UW = UW * time_factor
          VW = VW * time_factor
          UU_sgs = UU_sgs * time_factor
          VV_sgs = VV_sgs * time_factor
          WW_sgs = WW_sgs * time_factor
          UV_sgs = UV_sgs * time_factor
          UW_sgs = UW_sgs * time_factor
          VW_sgs = VW_sgs * time_factor
          TKE_sgs = TKE_sgs * time_factor
        end if
        
      end if
    
      if (store%avg_UU_prime>0) then
        UU = UU - U**2
        VV = VV - V**2
        WW = WW - W**2
        UV = UV - (U(0:Prnx-1,1:Prny,1:Prnz)+U(1:Prnx,1:Prny,1:Prnz)) * &
                  (V(1:Prnx,0:Prny-1,1:Prnz)+V(1:Prnx,1:Prny,1:Prnz)) / 4
        UW = UW - (U(0:Prnx-1,1:Prny,1:Prnz)+U(1:Prnx,1:Prny,1:Prnz)) * &
                  (W(1:Prnx,1:Prny,0:Prnz-1)+W(1:Prnx,1:Prny,1:Prnz)) / 4
        VW = VW - (V(1:Prnx,0:Prny-1,1:Prnz)+V(1:Prnx,1:Prny,1:Prnz)) * &
                  (W(1:Prnx,1:Prny,0:Prnz-1)+W(1:Prnx,1:Prny,1:Prnz)) / 4
      end if

      if (store%avg==1) then
        allocate(tmp(1:3,1:Prnx,1:Prny,1:Prnz))

        call newunit(unit)

        open(unit,file=output_dir//"avg.vtk", &
          access='stream',status='replace',form="unformatted",action="write")

        write(unit) "# vtk DataFile Version 2.0", lf
        write(unit) "CLMM output file", lf
        write(unit) "BINARY", lf
        write(unit) "DATASET RECTILINEAR_GRID", lf
        str="DIMENSIONS"
        write(str(12:),*) Prnx,Prny,Prnz
        write(unit) str, lf
        str="X_COORDINATES"
        write(str(15:),'(i5,2x,a)') Prnx,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(xPr(1:Prnx), real32)), lf
        str="Y_COORDINATES"
        write(str(15:),'(i5,2x,a)') Prny,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(yPr(1:Prny), real32)), lf
        str="Z_COORDINATES"
        write(str(15:),'(i5,2x,a)') Prnz,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(zPr(1:Prnz), real32)), lf
        str="POINT_DATA"
        write(str(12:),*) Prnx*Prny*Prnz
        write(unit) str, lf

        if (store%avg_Pr==1) then
          write(unit) "SCALARS p float", lf
          write(unit) "LOOKUP_TABLE default", lf

          write(unit) BigEnd(real(Pr(1:Prnx,1:Prny,1:Prnz), real32))

          write(unit) lf
        end if

        if (store%avg_Prtype==1) then
          write(unit) "SCALARS ptype float", lf
          write(unit) "LOOKUP_TABLE default", lf

          write(unit) BigEnd(real(Prtype(1:Prnx,1:Prny,1:Prnz), real32))

          write(unit) lf
        end if

        if (enable_buoyancy.and.store%avg_temperature==1) then
          write(unit) "SCALARS temperature float", lf
          write(unit) "LOOKUP_TABLE default", lf

          write(unit) BigEnd(real(Temperature(1:Prnx,1:Prny,1:Prnz), real32))

          write(unit) lf
        end if

        if (enable_moisture.and.store%avg_moisture==1) then
          write(unit) "SCALARS moisture float", lf
          write(unit) "LOOKUP_TABLE default", lf

          write(unit) BigEnd(real(Moisture(1:Prnx,1:Prny,1:Prnz), real32))

          write(unit) lf
        end if

        if (btest(store%avg_U,0)) then
          write(unit) "VECTORS u float", lf

          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              tmp(:,i,j,k) = [ BigEnd(real((U(i,j,k)+U(i-1,j,k))/2._knd, real32)), &
                               BigEnd(real((V(i,j,k)+V(i,j-1,k))/2._knd, real32)), &
                               BigEnd(real((W(i,j,k)+W(i,j,k-1))/2._knd, real32)) ]
            end do
           end do
          end do

          write(unit) tmp

          write(unit) lf
        end if

        if (btest(store%avg_UU_prime,0)) then
        
          allocate(sc_tmp(1:Prnx, 1:Prny, 1:Prnz))
          
          write(unit) "SCALARS uu_prime float", lf
          write(unit) "LOOKUP_TABLE default", lf
          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              sc_tmp(i,j,k) = BigEnd(real((UU(i,j,k)+UU(i-1,j,k))/2._knd, real32))
            end do
           end do
          end do
          write(unit) sc_tmp, lf
          
          write(unit) "SCALARS vv_prime float", lf
          write(unit) "LOOKUP_TABLE default", lf
          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              sc_tmp(i,j,k) = BigEnd(real((VV(i,j,k)+VV(i,j-1,k))/2._knd, real32))
            end do
           end do
          end do
          write(unit) sc_tmp, lf
          
          write(unit) "SCALARS ww_prime float", lf
          write(unit) "LOOKUP_TABLE default", lf
          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              sc_tmp(i,j,k) = BigEnd(real((WW(i,j,k)+WW(i,j,k-1))/2._knd, real32))
            end do
           end do
          end do
          write(unit) sc_tmp, lf
          
          write(unit) "SCALARS uu_prime_sgs float", lf
          write(unit) "LOOKUP_TABLE default", lf
          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              sc_tmp(i,j,k) = BigEnd(real((UU_sgs(i,j,k)+UU_sgs(i-1,j,k))/2._knd, real32))
            end do
           end do
          end do
          write(unit) sc_tmp, lf
          
          write(unit) "SCALARS vv_prime_sgs float", lf
          write(unit) "LOOKUP_TABLE default", lf
          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              sc_tmp(i,j,k) = BigEnd(real((VV_sgs(i,j,k)+VV_sgs(i,j-1,k))/2._knd, real32))
            end do
           end do
          end do
          write(unit) sc_tmp, lf
          
          write(unit) "SCALARS ww_prime_sgs float", lf
          write(unit) "LOOKUP_TABLE default", lf
          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              sc_tmp(i,j,k) = BigEnd(real((WW_sgs(i,j,k)+WW_sgs(i,j,k-1))/2._knd, real32))
            end do
           end do
          end do
          write(unit) sc_tmp, lf
          
          write(unit) "SCALARS tke_sgs float", lf
          write(unit) "LOOKUP_TABLE default", lf
          write(unit) BigEnd(real(TKE_sgs, real32)), lf
          
          write(unit) "SCALARS uv_prime float", lf
          write(unit) "LOOKUP_TABLE default", lf
          write(unit) BigEnd(real(UV, real32)), lf
          
          write(unit) "SCALARS uw_prime float", lf
          write(unit) "LOOKUP_TABLE default", lf
          write(unit) BigEnd(real(UW, real32)), lf
          
          write(unit) "SCALARS vw_prime float", lf
          write(unit) "LOOKUP_TABLE default", lf
          write(unit) BigEnd(real(VW, real32)), lf
          
          write(unit) "SCALARS uv_prime_sgs float", lf
          write(unit) "LOOKUP_TABLE default", lf
          write(unit) BigEnd(real(UV_sgs, real32)), lf
          
          write(unit) "SCALARS uw_prime_sgs float", lf
          write(unit) "LOOKUP_TABLE default", lf
          write(unit) BigEnd(real(UW_sgs, real32)), lf
          
          write(unit) "SCALARS vw_prime_sgs float", lf
          write(unit) "LOOKUP_TABLE default", lf
          write(unit) BigEnd(real(VW_sgs, real32)), lf
        end if

        if (store%avg_vorticity==1) then
          write(unit) "VECTORS vort float", lf

          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              tmp(:,i,j,k) = [ BigEnd(real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin) &
                                         -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin), real32)), &
                               BigEnd(real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin) &
                                         -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin), real32)), &
                               BigEnd(real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin) &
                                         -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin), real32)) ]
            end do
           end do
          end do

          write(unit) tmp

          write(unit) lf
        end if

        close(unit)

      end if !store%avg

      if (btest(store%avg_U,1)) &
        call OutputUVW(U, V, W, &
                       output_dir//"Uavg.vtk", &
                       output_dir//"Vavg.vtk", &
                       output_dir//"Wavg.vtk", &
                       .true.)

      if (btest(store%avg_UU_prime,1)) then
        call OutputUVW(UU, VV, WW, &
                       output_dir//"Urms.vtk", &
                       output_dir//"Vrms.vtk", &
                       output_dir//"Wrms.vtk", &
                       .true.)
        call OutputUVW(UU_sgs, VV_sgs, WW_sgs, &
                       output_dir//"Urms_sgs.vtk", &
                       output_dir//"Vrms_sgs.vtk", &
                       output_dir//"Wrms_sgs.vtk", &
                       .true.)
      end if

    end if !averaging

  end subroutine OutputAvg

  subroutine OutputScalarStats(S_avg, S_var, S_max, S_int)
    real(knd), dimension(-2:,-2:,-2:,:), contiguous, intent(inout) :: S_avg, S_var
    real(knd), dimension(-2:,-2:,-2:,:), contiguous, intent(in) :: S_max, S_int
    character(70) :: str
    integer :: unit
    real(knd) :: time_factor

    if (averaging==1.and. &
         (store%scalars_avg==1.or. &
          store%scalars_max==1.or. &
          store%scalars_intermitency==1) &
        .and. num_of_scalars>0 &
        .and. time_stepping%time > timeavg1) then
       
      if (time_stepping%time < timeavg2) then
        time_factor = (timeavg2 - timeavg1) / (time_stepping%time - timeavg1)
        
        S_avg = S_avg * time_factor
        if (store%scalars_variance==1) then
          S_var = S_var * time_factor
        end if
      end if
      
      if (store%scalars_variance==1) S_var = S_var - S_avg**2
    
      open(newunit=unit,file=output_dir//"scalars_avg.vtk", &
        access='stream',status='replace',form="unformatted",action="write")

      write(unit) "# vtk DataFile Version 2.0", lf
      write(unit) "CLMM output file", lf
      write(unit) "BINARY", lf
      write(unit) "DATASET RECTILINEAR_GRID", lf
      str="DIMENSIONS"
      write(str(12:),*) Prnx,Prny,Prnz
      write(unit) str, lf
      str="X_COORDINATES"
      write(str(15:),'(i5,2x,a)') Prnx,"float"
      write(unit) str, lf
      write(unit) BigEnd(real(xPr(1:Prnx), real32)), lf
      str="Y_COORDINATES"
      write(str(15:),'(i5,2x,a)') Prny,"float"
      write(unit) str, lf
      write(unit) BigEnd(real(yPr(1:Prny), real32)), lf
      str="Z_COORDINATES"
      write(str(15:),'(i5,2x,a)') Prnz,"float"
      write(unit) str, lf
      write(unit) BigEnd(real(zPr(1:Prnz), real32)), lf
      str="POINT_DATA"
      write(str(12:),*) Prnx*Prny*Prnz
      write(unit) str, lf

      if (store%scalars_avg==1) then
        call aux(S_avg,'_avg')
      end if
      
      if (store%scalars_variance==1) then
        call aux(S_var,'_variance')
      end if
      
      if (store%scalars_max==1) then
        call aux(S_max,'_max')
      end if
      
      if (store%scalars_intermitency==1) then
        call aux(S_int,'_intermitency')
      end if
      
      close(unit)

    end if
    
  contains
  
    subroutine aux(S,suff)
      real(knd), dimension(-2:,-2:,-2:,:), contiguous, intent(in) :: S
      character(*), intent(in) :: suff
      integer :: l
      character(8) ::  scalname="scalar00"

      do l = 1,num_of_scalars
        write(scalname(7:8),"(I2.2)") l
        write(unit) "SCALARS ", scalname//suff, " float", lf
        write(unit) "LOOKUP_TABLE default", lf

        write(unit) BigEnd(real(S(1:Prnx,1:Prny,1:Prnz,l), real32))

        write(unit) lf
      end do
    end subroutine
  end subroutine OutputScalarStats




  subroutine OutputUVW(U,V,W,fnameU,fnameV,fnameW,avg_mode_arg)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
    character(*), intent(in) :: fnameU,fnameV,fnameW
    logical,optional, intent(in) :: avg_mode_arg
    character(70) :: str
    integer :: i,unit
    logical :: avg_mode

    if (present(avg_mode_arg)) then
      avg_mode = avg_mode_arg
    else
      avg_mode = .false.
    end if

    if (store%U==1.or.avg_mode) then

        call newunit(unit)

        open(unit,file=fnameU, &
          access='stream',status='replace',form="unformatted",action="write")

        write(unit) "# vtk DataFile Version 2.0", lf
        write(unit) "CLMM output file", lf
        write(unit) "BINARY", lf
        write(unit) "DATASET RECTILINEAR_GRID", lf
        str="DIMENSIONS"
        write(str(12:),*) Unx,Uny,Unz
        write(unit) str, lf
        str="X_COORDINATES"
        write(str(15:),'(i5,2x,a)') Unx,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(xU(1:Unx), real32)), lf
        str="Y_COORDINATES"
        write(str(15:),'(i5,2x,a)') Uny,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(yPr(1:Uny), real32)), lf
        str="Z_COORDINATES"
        write(str(15:),'(i5,2x,a)') Unz,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(zPr(1:Unz), real32)), lf
        str="POINT_DATA"
        write(str(12:),*) Unx*Uny*Unz
        write(unit) str, lf


        write(unit) "SCALARS U float", lf
        write(unit) "LOOKUP_TABLE default", lf

        write(unit) BigEnd(real(U(1:Unx,1:Uny,1:Unz), real32))

        write(unit) lf

        if (store%U_interp==1) then
          write(unit) "SCALARS Utype float", lf
          write(unit) "LOOKUP_TABLE default", lf

          write(unit) BigEnd(real(Utype(1:Unx,1:Uny,1:Unz), real32))

          write(unit) lf
        end if

        close(unit)
    end if

    if (store%V==1.or.avg_mode) then

        call newunit(unit)

        open(unit,file=fnameV, &
          access='stream',status='replace',form="unformatted",action="write")

        write(unit) "# vtk DataFile Version 2.0", lf
        write(unit) "CLMM output file", lf
        write(unit) "BINARY", lf
        write(unit) "DATASET RECTILINEAR_GRID", lf
        str="DIMENSIONS"
        write(str(12:),*) Vnx,Vny,Vnz
        write(unit) str, lf
        str="X_COORDINATES"
        write(str(15:),'(i5,2x,a)') Vnx,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(xPr(1:Vnx), real32)), lf
        str="Y_COORDINATES"
        write(str(15:),'(i5,2x,a)') Vny,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(yV(1:Vny), real32)), lf
        str="Z_COORDINATES"
        write(str(15:),'(i5,2x,a)') Vnz,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(zPr(1:Vnz), real32)), lf
        str="POINT_DATA"
        write(str(12:),*) Vnx*Vny*Vnz
        write(unit) str, lf


        write(unit) "SCALARS V float", lf
        write(unit) "LOOKUP_TABLE default", lf

        write(unit) BigEnd(real(V(1:Vnx,1:Vny,1:Vnz), real32))

        write(unit) lf

        if (store%V_interp==1) then
          write(unit) "SCALARS Vtype float", lf
          write(unit) "LOOKUP_TABLE default", lf

          write(unit) BigEnd(real(Vtype(1:Vnx,1:Vny,1:Vnz), real32))

          write(unit) lf
        end if

        close(unit)
    end if !store%V

    if (store%W==1.or.avg_mode) then

        call newunit(unit)

        open(unit,file=fnameW, &
          access='stream',status='replace',form="unformatted",action="write")

        write(unit) "# vtk DataFile Version 2.0", lf
        write(unit) "CLMM output file", lf
        write(unit) "BINARY", lf
        write(unit) "DATASET RECTILINEAR_GRID", lf
        str="DIMENSIONS"
        write(str(12:),*) Wnx,Wny,Wnz
        write(unit) str, lf
        str="X_COORDINATES"
        write(str(15:),'(i5,2x,a)') Wnx,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(xPr(1:Wnx), real32)), lf
        str="Y_COORDINATES"
        write(str(15:),'(i5,2x,a)') Wny,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(yPr(1:Wny), real32)), lf
        str="Z_COORDINATES"
        write(str(15:),'(i5,2x,a)') Wnz,"float"
        write(unit) str, lf
        write(unit) BigEnd(real(zW(1:Wnz), real32)), lf
        str="POINT_DATA"
        write(str(12:),*) Wnx*Wny*Wnz
        write(unit) str, lf


        write(unit) "SCALARS W float", lf
        write(unit) "LOOKUP_TABLE default", lf

        write(unit) BigEnd(real(W(1:Wnx,1:Wny,1:Wnz), real32))

        write(unit) lf

        if (store%W_interp==1) then
          write(unit) "SCALARS Wtype float", lf
          write(unit) "LOOKUP_TABLE default", lf

          write(unit) BigEnd(real(Wtype(1:Wnx,1:Wny,1:Wnz), real32))

          write(unit) lf
        end if

        close(unit)
    end if !store%W

    if (store%U_interp/=0 .and. .not.avg_mode) then

      call newunit(unit)

      open(unit,file=output_dir//"Uinterp.txt")
      do i = 1,size(UIBPoints)
        write(unit,*) "xi,yj,zk",UIBPoints(i)%xi,UIBPoints(i)%yj,UIBPoints(i)%zk
        write(unit,*) "interp",UIBPoints(i)%interp
        write(unit,*) "xi",UIBPoints(i)%IntPoints%xi
        write(unit,*) "yj",UIBPoints(i)%IntPoints%yj
        write(unit,*) "zk",UIBPoints(i)%IntPoints%zk
        write(unit,*) "coefs",UIBPoints(i)%IntPoints%coef
      end do
      close(unit)
    end if

    if (store%V_interp/=0 .and. .not.avg_mode) then

      call newunit(unit)

      open(unit,file=output_dir//"Vinterp.txt")
      do i = 1,size(VIBPoints)
        write(unit,*) "xi,yj,zk",VIBPoints(i)%xi,VIBPoints(i)%yj,VIBPoints(i)%zk
        write(unit,*) "interp",VIBPoints(i)%interp
        write(unit,*) "xi",VIBPoints(i)%IntPoints%xi
        write(unit,*) "yj",VIBPoints(i)%IntPoints%yj
        write(unit,*) "zk",VIBPoints(i)%IntPoints%zk
        write(unit,*) "coefs",VIBPoints(i)%IntPoints%coef
      end do
      close(unit)
    end if

    if (store%W_interp/=0) then

      call newunit(unit)

      open(unit,file=output_dir//"Winterp.txt")
      do i = 1,size(WIBPoints)
        write(unit,*) "xi,yj,zk",WIBPoints(i)%xi,WIBPoints(i)%yj,WIBPoints(i)%zk
        write(unit,*) "interp",WIBPoints(i)%interp
        write(unit,*) "xi",WIBPoints(i)%IntPoints%xi
        write(unit,*) "yj",WIBPoints(i)%IntPoints%yj
        write(unit,*) "zk",WIBPoints(i)%IntPoints%zk
        write(unit,*) "coefs",WIBPoints(i)%IntPoints%coef
      end do
      close(unit)
    end if

    if (store%U_interp/=0.and.store%V_interp/=0.and.store%W_interp/=0) then

      call newunit(unit)

      open(unit,file=output_dir//"Scinterp.txt")
      do i = 1,size(ScalFlIBPoints)
        write(unit,*) "xi,yj,zk",ScalFlIBPoints(i)%xi,ScalFlIBPoints(i)%yj,ScalFlIBPoints(i)%zk
        write(unit,*) "interp",ScalFlIBPoints(i)%interp
        write(unit,*) "xi",ScalFlIBPoints(i)%IntPoints%xi
        write(unit,*) "yj",ScalFlIBPoints(i)%IntPoints%yj
        write(unit,*) "zk",ScalFlIBPoints(i)%IntPoints%zk
        write(unit,*) "coefs",ScalFlIBPoints(i)%IntPoints%coef
      end do
      close(unit)
    end if

  end subroutine OutputUVW


  subroutine OutputAvgFluxes
    real(knd), allocatable :: Scalar_fl_U_adv(:,:,:,:)
    real(knd), allocatable :: Scalar_fl_V_adv(:,:,:,:)
    real(knd), allocatable :: Scalar_fl_W_adv(:,:,:,:)
    real(knd), allocatable :: Scalar_fl_U_turb(:,:,:,:)
    real(knd), allocatable :: Scalar_fl_V_turb(:,:,:,:)
    real(knd), allocatable :: Scalar_fl_W_turb(:,:,:,:)
    integer :: i
    real(knd) :: time_factor
   
    if (num_of_scalars>0 .and. store%avg_flux_scalar==1 .and. time_stepping%time > timeavg1) then

      if (time_stepping%time < timeavg2) then
        time_factor = (timeavg2 - timeavg1) / (time_stepping%time - timeavg1)
        
        Scalar_fl_U_avg = Scalar_fl_U_avg * time_factor
        Scalar_fl_V_avg = Scalar_fl_V_avg * time_factor
        Scalar_fl_W_avg = Scalar_fl_W_avg * time_factor
        if (store%avg_flux_scalar_sgs==1) then
          Scalar_fl_U_sgs = Scalar_fl_U_sgs * time_factor
          Scalar_fl_V_sgs = Scalar_fl_V_sgs * time_factor
          Scalar_fl_W_sgs = Scalar_fl_W_sgs * time_factor
        end if
      end if
      
      
      !FIXME: just to allocate, mold= does not work correctly as of gcc 4.8
!       allocate(Scalar_fl_U_adv,mold=Scalar_fl_U_avg)
      Scalar_fl_U_adv = Scalar_fl_U_avg
      Scalar_fl_V_adv = Scalar_fl_V_avg
      Scalar_fl_W_adv = Scalar_fl_W_avg
      Scalar_fl_U_adv = 0
      Scalar_fl_V_adv = 0
      Scalar_fl_W_adv = 0

      do i=1,num_of_scalars
        call AddScalarAdvVector(Scalar_fl_U_adv(:,:,:,i), &
                             Scalar_fl_V_adv(:,:,:,i), &
                             Scalar_fl_W_adv(:,:,:,i), &
                             Scalar_avg(:,:,:,i), &
                             U_avg,V_avg,W_avg, &
                             1._knd, &
                             scalar_fluxes_time(:,i,:,1), &
                             scalar_probes%i, scalar_probes%j, scalar_probes%k)
      end do

      if (store%avg_flux_scalar_sgs==1) then
        Scalar_fl_U_avg = Scalar_fl_U_avg + Scalar_fl_U_sgs
        Scalar_fl_V_avg = Scalar_fl_V_avg + Scalar_fl_V_sgs
        Scalar_fl_W_avg = Scalar_fl_W_avg + Scalar_fl_W_sgs
      end if
      
      Scalar_fl_U_turb = Scalar_fl_U_avg - Scalar_fl_U_adv
      Scalar_fl_V_turb = Scalar_fl_V_avg - Scalar_fl_V_adv
      Scalar_fl_W_turb = Scalar_fl_W_avg - Scalar_fl_W_adv
      

      call SaveScalarVTKFluxes(Scalar_fl_U_avg, &
                               Scalar_fl_U_adv, &
                               Scalar_fl_U_turb, &
                               Scalar_fl_U_sgs, &
                               output_dir//"scalflu.vtk",xU,yPr,zPr)
      call SaveScalarVTKFluxes(Scalar_fl_V_avg, &
                               Scalar_fl_V_adv, &
                               Scalar_fl_V_turb, &
                               Scalar_fl_V_sgs, &
                               output_dir//"scalflv.vtk",xPr,yV,zPr)
      call SaveScalarVTKFluxes(Scalar_fl_W_avg, &
                               Scalar_fl_W_adv, &
                               Scalar_fl_W_turb, &
                               Scalar_fl_W_sgs, &
                               output_dir//"scalflw.vtk",xPr,yPr,zW)

    end if

  end subroutine

  subroutine SaveScalarVTKFluxes(Avg,Adv,Turb,Sgs,file_name,x,y,z)
    real(knd), dimension(:,:,:,:), contiguous, intent(in) :: Avg, Adv, Turb
    real(knd), dimension(:,:,:,:), allocatable, intent(in) :: Sgs
    character(*), intent(in) :: file_name
    real(knd), allocatable, intent(in) :: x(:), y(:), z(:)
    integer :: nx,ny,nz
    character(70) :: str
    character(13) ::  scalname
    integer :: l, unit

    nx = size(Avg,1)
    ny = size(Avg,2)
    nz = size(Avg,3)

    call newunit(unit)

    open(unit,file=trim(file_name), &
      access='stream',status='replace',form="unformatted",action="write")

    write(unit) "# vtk DataFile Version 2.0", lf
    write(unit) "CLMM output file", lf
    write(unit) "BINARY", lf
    write(unit) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),*) nx,ny,nz
    write(unit) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(x(1:nx), real32)), lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(y(1:ny), real32)), lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(z(1:nz), real32)), lf
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit) str, lf

    scalname = "scalfl-avg-"

    do l = 1,num_of_scalars
      write(scalname(12:13),"(I2.2)") l
      write(unit) "SCALARS ", scalname , " float", lf
      write(unit) "LOOKUP_TABLE default", lf

      write(unit) BigEnd(real(Avg(:,:,:,l), real32))

      write(unit) lf
    end do

    scalname = "scalfl-adv-"

    do l = 1,num_of_scalars
      write(scalname(12:13),"(I2.2)") l
      write(unit) "SCALARS ", scalname , " float", lf
      write(unit) "LOOKUP_TABLE default", lf

      write(unit) BigEnd(real(Adv(:,:,:,l), real32))

      write(unit) lf
    end do

    scalname = "scalfl-trb-"

    do l = 1,num_of_scalars
      write(scalname(12:13),"(I2.2)") l
      write(unit) "SCALARS ", scalname , " float", lf
      write(unit) "LOOKUP_TABLE default", lf

      write(unit) BigEnd(real(Turb(:,:,:,l), real32))

      write(unit) lf
    end do
    
    if (store%avg_flux_scalar_sgs==1) then
      scalname = "scalfl-sgs-"

      do l = 1,num_of_scalars
        write(scalname(12:13),"(I2.2)") l
        write(unit) "SCALARS ", scalname , " float", lf
        write(unit) "LOOKUP_TABLE default", lf

        write(unit) BigEnd(real(Sgs(:,:,:,l), real32))

        write(unit) lf
      end do
    end if
    
    close(unit)
  end subroutine SaveScalarVTKFluxes



  subroutine Output(U,V,W,Pr,Temperature,Moisture,Scalar)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U,V,W
    real(knd), contiguous, intent(inout) :: Pr(-1:,-1:,-1:)
    real(knd), contiguous, intent(in) :: Temperature(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: Moisture(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: Scalar(-2:,-2:,-2:,:)
    
#ifdef CUSTOM_OUTPUT
    interface
      subroutine CustomOutput
      end subroutine
    end interface
#endif

    call BoundUVW(U, V, W)

    call OutputOut(U,V,W,Pr,Temperature,Moisture)

    call OutputScalars(Scalar)

    call OutputUVW(U,V,W,output_dir//"U.vtk",output_dir//"V.vtk",output_dir//"W.vtk")
    
    call OutputTimeSeries

    if (averaging==1 .and. time_stepping%time>=timeavg1) then

      call OutputAvg(U_avg, V_avg, W_avg, &
                     Pr_avg, Temperature_avg, Moisture_avg, &
                     UU_prime, VV_prime, WW_prime, &
                     UV_prime, UW_prime, VW_prime, &
                     TKE_prime_sgs, UU_prime_sgs, VV_prime_sgs, WW_prime_sgs, &
                     UV_prime_sgs, UW_prime_sgs, VW_prime_sgs)
      
      call OutputScalarStats(Scalar_avg,Scalar_variance,Scalar_max,Scalar_intermitency)

      call OutputAvgFluxes

    end if

    call OutputProfiles

    call FinalizeVTKFrames
    
#ifdef PAR
    call FinalizeFrames_ParallelIO
#endif

    call FinalizeSurfaceFrames

    call FinalizeStaggeredFrames
    
#ifdef CUSTOM_OUTPUT
    call CustomOutput
#endif

    if (len_trim(scratch_dir)>0) then
      call CopyFromScratch(output_dir)
    end if

  end subroutine Output





  subroutine StressProfiles(U,V,W)
    use Filters, only: filtertype, filter_ratios
    use Subgrid, only: TKEDissipation
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
    integer   :: i,j,k
    real(knd) :: S, width
    real(knd), parameter :: Ck = 0.1 !model constant for the sgs_tke eddy viscosity sgs model

    width = filter_ratios(filtertype) * (dxmin*dymin*dzmin)**(1._knd/3._knd)
    
    !$omp parallel private(i,j,k,S)
    !$omp do
    do k = 1,Prnz
      S = 0
      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,k)<=0) then
           S = S + ( (Viscosity(i,j,k) - molecular_viscosity) / (Ck * width) )**2
         end if
       end do
      end do
      current_profiles%tkesgs(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 1,Prnz
      S = 0
      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,k)<=0) then
           S = S + TKEDissipation(i, j, k, U, V, W)
         end if
       end do
      end do
      current_profiles%dissip(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
      S = 0
      do j = 1,Uny
       do i = 1,Unx
         if ((Utype(i,j,k+1)<=0.or.Utype(i,j,k)<=0).and.(Wtype(i+1,j,k)<=0.or.Wtype(i,j,k)<=0)) then
           S = S + (U(i,j,k+1) + U(i,j,k)) * (W(i+1,j,k) + W(i,j,k)) / 4
         end if
       end do
      end do
      current_profiles%uw(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
      S = 0
      do j = 1,Vny
       do i = 1,Vnx
         if ((Vtype(i,j,k+1)<=0.or.Vtype(i,j,k)<=0).and.(Wtype(i,j+1,k)<=0.or.Wtype(i,j,k)<=0)) then
           S = S + (V(i,j,k+1) + V(i,j,k)) * (W(i,j+1,k) + W(i,j,k)) / 4
         end if
       end do
      end do
      current_profiles%vw(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
      S = 0
      do j = 1,Uny
       do i = 1,Unx
         if (Utype(i,j,k+1)<=0.or.Utype(i,j,k)<=0) then
           S = S - ( (Viscosity(i+1,j,k+1) + &
                      Viscosity(i+1,j,k) + &
                      Viscosity(i,j,k+1) + &
                      Viscosity(i,j,k)) / 4 &
                     - molecular_viscosity) * &
                   (U(i,j,k+1)-U(i,j,k)) / dzmin
         end if
       end do
      end do
      current_profiles%uwsgs(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
      S = 0
      do j = 1,Vny
       do i = 1,Vnx
         if (Vtype(i,j,k+1)<=0.or.Vtype(i,j,k)<=0) then
           S = S - ( (Viscosity(i,j+1,k+1) + &
                      Viscosity(i,j+1,k) + &
                      Viscosity(i,j,k+1) + &
                      Viscosity(i,j,k)) / 4 &
                     - molecular_viscosity) * &
                   (V(i,j,k+1)-V(i,j,k)) / dzmin
         end if
       end do
      end do
      current_profiles%vwsgs(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
      S = 0
      do j = 1,Uny
       do i = 1,Unx
         if ((Utype(i,j,k+1)<=0.or.Utype(i,j,k)<=0)) then
           S = S - molecular_viscosity * (U(i,j,k+1) - U(i,j,k)) / dzmin
         end if
       end do
      end do
      current_profiles%uz_visc(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
      S = 0
      do j = 1,Vny
       do i = 1,Vnx
         if ((Vtype(i,j,k+1)<=0.or.Vtype(i,j,k)<=0)) then
           S = S - molecular_viscosity * (V(i,j,k+1) - V(i,j,k)) / dzmin
         end if
       end do
      end do
      current_profiles%vz_visc(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 1,Unz
      S = 0
      do j = 1,Uny
       do i = 1,Unx
         if (Utype(i,j,k)<=0) then
           S = S - ((Viscosity(i,j,k+1)+Viscosity(i,j,k)) / 2 - molecular_viscosity) &
                   * (U(i+1,j,k)-U(i-1,j,k)) / dxmin / 2
         end if
       end do
      end do
      current_profiles%uusgs(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 1,Vnz
      S = 0
      do j = 1,Vny
       do i = 1,Vnx
         if (Vtype(i,j,k)<=0) then
           S = S - ((Viscosity(i,j+1,k)+Viscosity(i,j,k)) / 2 - molecular_viscosity) &
                   * (V(i,j+1,k)-V(i,j-1,k)) / dymin / 2
         end if
       end do
      end do
      current_profiles%vvsgs(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Wnz
      S = 0
      do j = 1,Wny
       do i = 1,Wnx
         if (Wtype(i,j,k)<=0) then
           S = S - ((Viscosity(i,j,k+1)+Viscosity(i,j,k)) / 2 - molecular_viscosity) &
                   * (W(i,j,k+1)-W(i,j,k-1)) / dzmin / 2
         end if
       end do
      end do
      current_profiles%wwsgs(k) = S
    end do
    !$omp end do nowait
    !$omp end parallel
    
  end subroutine StressProfiles
  
  
  subroutine FluxSGSProfiles(W,Temperature,Moisture,Scalar)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Temperature
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Moisture
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(in) :: Scalar

    
    if (enable_buoyancy) call TemperatureFluxSGSProfile(W,Temperature)
    if (enable_moisture) call MoistureFluxSGSProfile(W,Moisture)
    if (num_of_scalars>0) call ScalarFluxSGSProfile(W,Scalar)
  end subroutine


  subroutine TemperatureFluxSGSProfile(W,Temperature)
    use Wallmodels

    use custom_par, only: kim

    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Temperature
    real(knd) :: S
    integer   :: i,j,k
    !current_profiles%tempfl is computed directly during advection step

    !$omp parallel 
    !$omp do private(i,j,k,S)
    do k = 1,Prnz
      S = 0
      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
           S = S - ( ( (TDiff(i,j,k+1)+TDiff(i,j,k)) / 2 - molecular_diffusivity) &
                    * (Temperature(i,j,k+1)-Temperature(i,j,k))) / dzmin
         end if
       end do
      end do
      current_profiles%tempflsgs(k) = S
    end do  

    !$omp do private(i,j,k,S)
    do k = 0,Prnz
      S = 0
      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
           S = S - ( molecular_diffusivity &
                    * (Temperature(i,j,k+1)-Temperature(i,j,k))) / dzmin
         end if
       end do
      end do
      current_profiles%tempflvisc(k) = S
    end do

    if (kim==1) then
      !$omp single
      S = 0
      !$omp end single
      !$omp do private(i) reduction(+:S)
      do i = 1, size(WMPoints)
        if (WMPoints(i)%zk==1.and.Prtype(WMPoints(i)%xi,WMPoints(i)%yj,1)<=0) then
          S = S + WMPoints(i)%temperature_flux
        end if
      end do
      
      !$omp single
      if (abs(S)<abs(current_profiles%tempflvisc(0))) then
        current_profiles%tempflvisc(0) = 0
      else
        S = S - current_profiles%tempflvisc(0)
      end if
      current_profiles%tempflsgs(0) = S
      !$omp end single
    endif
    !$omp end parallel
   
  end subroutine

  subroutine MoistureFluxSGSProfile(W,Moisture)
    use Wallmodels

    use custom_par, only: kim

    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Moisture
    real(knd) :: S
    integer   :: i,j,k
    !current_profiles%moistfl is computed directly during advection step

    !$omp parallel 
    !$omp do private(i,j,k,S)
    do k = 0,Prnz
      S = 0
      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
           S = S - ( ( (TDiff(i,j,k+1)+TDiff(i,j,k)) / 2 - molecular_diffusivity) &
                    * (Moisture(i,j,k+1)-Moisture(i,j,k))) / dzmin
         end if
       end do
      end do
      current_profiles%moistflsgs(k) = S
    end do

    !$omp do private(i,j,k,S)
    do k = 0,Prnz
      S = 0
      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
           S = S - ( molecular_diffusivity &
                    * (Moisture(i,j,k+1)-Moisture(i,j,k))) / dzmin
         end if
       end do
      end do
      current_profiles%moistflvisc(k) = S
    end do

    if (kim==1) then
      !$omp single
      S = 0
      !$omp end single
      !$omp do private(i) reduction(+:S)
      do i = 1, size(WMPoints)
        if (WMPoints(i)%zk==1.and.Prtype(WMPoints(i)%xi,WMPoints(i)%yj,1)<=0) then
          S = S + WMPoints(i)%moisture_flux
        end if
      end do

      !$omp single      
      if (abs(S)<abs(current_profiles%moistflvisc(0))) then
        current_profiles%moistflvisc(0) = 0
      else
        S = S - current_profiles%moistflvisc(0)
      end if
      current_profiles%moistflsgs(0) = S
      !$omp end single
    endif
    !$omp end parallel
   
  end subroutine
      
  subroutine ScalarFluxSGSProfile(W,Scalar)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: W
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(in) :: Scalar
    real(knd) :: S
    integer   :: i,j,k,l

    !$omp parallel private(i,j,k,S)
    do l = 1,num_of_scalars
      !$omp do
      do k = 0,Prnz
        S = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
            S = S + 0.5_knd * (Scalar(i,j,k+1,l)+Scalar(i,j,k,l)) * (W(i,j,k))
           end if
         end do
        end do
        current_profiles%scalfl(l,k) = S
      end do
      !$omp end do nowait
      !$omp do
      do k = 0,Prnz
        S = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
             S = S - ( ( (TDiff(i,j,k+1)+TDiff(i,j,k)) / 2 - molecular_diffusivity) &
                      * (Scalar(i,j,k+1,l)-Scalar(i,j,k,l))) / dzmin
           end if
         end do
        end do
        current_profiles%scalflsgs(l,k) = S
      end do
      !$omp end do
    end do
    !$omp end parallel
  end subroutine

  


  subroutine BLProfiles(U,V,W,Temperature,Moisture,Scalar)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Temperature
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Moisture
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(in) :: Scalar
    real(knd) :: S
    integer   :: i,j,k,l

    !$omp parallel private(i,j,k,S)
    !$omp do
    do k = 0,Unz+1
      S = 0
      do j = 1,Uny
       do i = 1,Unx
         if (Utype(i,j,k)<=0) then
           S = S + U(i,j,k)
         end if
       end do
      end do
      current_profiles%u(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 1,Vnz+1
      S = 0
      do j = 1,Vny
       do i = 1,Vnx
         if (Vtype(i,j,k)<=0) then
           S = S + V(i,j,k)
         end if
       end do
      end do
      current_profiles%v(k) = S
    end do
    !$omp end do nowait

    if (size(Temperature)>0) then
      !$omp do
      do k = 1,Prnz
        S = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k)<=0) then
             S = S + Temperature(i,j,k)
           end if
         end do
        end do
        current_profiles%temp(k) = S
      end do
      !$omp end do
    end if

    if (size(Moisture)>0) then
      !$omp do
      do k = 1,Prnz
        S = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k)<=0) then
             S = S + Moisture(i,j,k)
           end if
         end do
        end do
        current_profiles%moist(k) = S
      end do
      !$omp end do
    end if

    if (size(Scalar)>0) then
      do l = 1,num_of_scalars
        !$omp do
        do k = 1,Prnz
          S = 0
          do j = 1,Prny
           do i = 1,Prnx
             if (Prtype(i,j,k)<=0) then
               S = S + Scalar(i,j,k,l)
             end if
           end do
          end do
          current_profiles%scal(l,k) = S
        end do
        !$omp end do
      end do
    end if

    !$omp do
    do k = 1,Unz
      S = 0
      do j = 1,Uny
       do i = 1,Unx
         if (Utype(i,j,k)<=0) then
           S = S + (U(i,j,k))**2
         end if
       end do
      end do
      current_profiles%uu(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 1,Vnz
      S = 0
      do j = 1,Vny
       do i = 1,Vnx
         if (Vtype(i,j,k)<=0) then
          S = S + (V(i,j,k))**2
         end if
       end do
      end do
      current_profiles%vv(k) = S
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Wnz
      S = 0
      do j = 1,Wny
       do i = 1,Wnx
         if (Wtype(i,j,k)<=0) then
           S = S + (W(i,j,k))**2
         end if
       end do
      end do
      current_profiles%ww(k) = S
    end do
    !$omp end do

    if (enable_buoyancy) then
      !$omp do
      do k = 1,Prnz
        S = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k)<=0) then
             S = S + (Temperature(i,j,k))**2
           end if
         end do
        end do
        current_profiles%tt(k) = S
      end do
      !$omp end do
    end if ! size(Temperature)


    if (enable_moisture) then
      !$omp do
      do k = 1,Prnz
        S = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k)<=0) then
             S = S + (Moisture(i,j,k))**2
           end if
         end do
        end do
        current_profiles%mm(k) = S
      end do
      !$omp end do nowait

    end if ! size(Moisture)



    do l = 1,num_of_scalars
      !$omp do
      do k = 1,Prnz
        S = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k)<=0) then
            S = S + (Scalar(i,j,k,l))**2
           end if
         end do
        end do
        current_profiles%ss(l,k) = S
      end do
      !$omp end do
    end do
    !$omp end parallel

  end subroutine BLProfiles




  subroutine OutputU2(U,V,W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
    integer :: unit
    character(70) :: str

    call newunit(unit)

    open(unit,file=output_dir//"U2.vtk")
    write(unit) "# vtk DataFile Version 2.0", lf
    write(unit) "CLMM output file", lf
    write(unit) "BINARY", lf
    write(unit) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),*) Unx,Uny,Unz
    write(unit) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') Unx,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(xU(1:Unx), real32)), lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') Uny,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(yPr(1:Uny), real32)), lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') Unz,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(zPr(1:Unz), real32)), lf
    str="POINT_DATA"
    write(str(12:),*) Unx*Uny*Unz
    write(unit) str, lf


    write(unit) "SCALARS U float", lf
    write(unit) "LOOKUP_TABLE default", lf

    write(unit) BigEnd(real(U(1:Unx,1:Uny,1:Unz), real32)), lf

    write(unit) lf
    close(unit)


    open(unit,file=output_dir//"V2.vtk")
    write(unit) "# vtk DataFile Version 2.0", lf
    write(unit) "CLMM output file", lf
    write(unit) "BINARY", lf
    write(unit) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),*) Vnx,Vny,Vnz
    write(unit) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') Vnx,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(xPr(1:Vnx), real32)), lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') Vny,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(yV(1:Vny), real32)), lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') Vnz,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(zPr(1:Vnz), real32)), lf
    str="POINT_DATA"
    write(str(12:),*) Vnx*Vny*Vnz
    write(unit) str, lf


    write(unit) "SCALARS V float", lf
    write(unit) "LOOKUP_TABLE default", lf

    write(unit) BigEnd(real(V(1:Vnx,1:Vny,1:Vnz), real32)), lf

    write(unit) lf
    close(unit)


    open(unit,file=output_dir//"W2.vtk")
    write(unit) "# vtk DataFile Version 2.0", lf
    write(unit) "CLMM output file", lf
    write(unit) "BINARY", lf
    write(unit) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),*) Wnx,Wny,Wnz
    write(unit) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') Wnx,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(xPr(1:Wnx), real32)), lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') Wny,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(yPr(1:Wny), real32)), lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') Wnz,"float"
    write(unit) str, lf
    write(unit) BigEnd(real(zW(1:Wnz), real32)), lf
    str="POINT_DATA"
    write(str(12:),*) Wnx*Wny*Wnz
    write(unit) str, lf


    write(unit) "SCALARS W float", lf
    write(unit) "LOOKUP_TABLE default", lf

    write(unit) BigEnd(real(W(1:Wnx,1:Wny,1:Wnz), real32)), lf

    write(unit) lf
    close(unit)
  end subroutine OutputU2


  subroutine OUTINLET(U,V,W,Temperature)
    use Parameters, t_s => time_stepping
    !for output of 2d data for use as an inilet condition later
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Temperature
    integer, save ::fnum
    integer, save :: called = 0

    if ((t_s%time>=timefram1).and.(t_s%time<=timefram2+(timefram2-timefram1)/(frames-1))&
        .and.(t_s%time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
     if (called==0) then
      open(101,file=output_dir//"inletframeinfo.unf",form='unformatted',status='replace',action='write')
      write(101) Prny,Prnz  !for check of consistency of grids before use
      write(101) Vny
      write(101) Wnz
      write(101) dxPr(0)
      called = 1
      fnum = 0
     end if
     fnum = fnum+1
     write(101) t_s%time-timefram1
     call OUTINLETFrame(U,V,W,Temperature,fnum)
    elseif (t_s%time>timefram2+(timefram2-timefram1)/(frames-1).and.called==1) then
      close(101)
      called = 2
    end if

  end subroutine OUTINLET


  subroutine OUTINLETFrame(U,V,W,Temperature,n)
    real(knd), intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(knd), dimension(-2:,-2:,-2:), intent(in)   :: Temperature
    integer :: n
    character(12) :: fname
    integer :: mini,maxi,minj,maxj,mink,maxk,unit

    call GridCoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._knd,(yV(Prny+1)+yV(0))/2._knd,(zW(Prnz+1)+zW(0))/2._knd)
    maxi = mini
    minj = 1
    maxj = Prny
    mink = 1
    maxk = Prnz


    fname(1:5)="frame"
    write(fname(6:8),"(I3.3)") n
    fname(9:12)=".unf"
    if (master) write(*,*) "Saving frame:",fname(1:6),"   time:", time_stepping%time

    call newunit(unit)

    open(unit,file = fname,form='unformatted',access='sequential',status='replace',action='write')


    write(unit) U(mini,1:Uny,1:Unz)
    write(unit) V(mini,1:Vny,1:Vnz)
    write(unit) W(mini,1:Wny,1:Wnz)
    if (enable_buoyancy) then
         write(unit) Temperature(mini,1:Prny,1:Prnz)
    end if
    close(unit)

  end subroutine OUTINLETFrame


  subroutine CopyFromScratch(out_dir)
! #ifdef PAR
!     use custom_par
!     integer :: ie
! #else
!     integer, parameter :: nims = 1, myim = 1
! #endif
    character(*), intent(in) :: out_dir
!     character(len(out_dir)) :: final_dir
!     integer :: i, l
!     integer, parameter :: par = 8 !number of paralelly copying images
! 
!     do i = 1, ceiling(real(nims)/par)
!       if (myim>(i-1)*par.and.myim<=i*par) then
!       
!         l = len_trim(scratch_dir)
! 
!         if (l>0) then
! 
!           final_dir = out_dir(l+1:)
! 
!           !FIXME: change to execute_command_line when more widely supported
! #if defined(_WIN32) || defined(_WIN64)
!           call system("mkdir "//trim(final_dir))
!           call system("Robocopy "//out_dir//"* "//trim(final_dir)//" /E")
! #else
!           !Will try 10 times
!           call system('n=0; until mkdir -p '//trim(final_dir)// &
!                                 ' || [ "$n" -ge "9" ]; do n=$(($n+1)); done')
!           call system('n=0; until cp -r '//out_dir//'* '//trim(final_dir)// &
!                                 ' || [ "$n" -ge "9" ]; do n=$(($n+1)); done')
! #endif
!         end if
!       end if
! #ifdef PAR
!       call MPI_Barrier(domain_comm, ie)
! #endif
!     end do
  end subroutine



  pure real(knd) function TriLinInt(a,b,c,vel000,vel100,vel010,vel001,vel110,vel101,vel011,vel111)
    real(knd), intent(in) :: a,b,c,vel000,vel100,vel010,vel001,vel110,vel101,vel011,vel111

    TriLinInt=   (1-a)*(1-b)*(1-c)*vel000+&
                 a*(1-b)*(1-c)*vel100+&
                 (1-a)*b*(1-c)*vel010+&
                 (1-a)*(1-b)*c*vel001+&
                 a*b*(1-c)*vel110+&
                 a*(1-b)*c*vel101+&
                 (1-a)*b*c*vel011+&
                 a*b*c*vel111

  end function TriLinInt



end module Outputs
