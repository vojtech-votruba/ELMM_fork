program CLMM

  use Parameters
  use Initial
  use TimeSteps,  only: TMarchRK3
  use Dynamics,   only: Dynamics_Deallocate
  use Scalars,    only: Scalars_Deallocate
  use Outputs,    only: AllocateOutputs, OutTStep, Output
  use Endianness, only: GetEndianness
  use custom_par


  implicit none




  real(knd), allocatable, dimension(:,:,:)   :: U           !Velocity component in x direction -- horizontal
  real(knd), allocatable, dimension(:,:,:)   :: V           !Velocity component in y direction -- horizontal
  real(knd), allocatable, dimension(:,:,:)   :: W           !Velocity component in z direction -- vertical
  real(knd), allocatable, dimension(:,:,:)   :: Pr          !Pressure (kinematic pressure with possible subgrid TKE part)
  real(knd), allocatable, dimension(:,:,:)   :: Temperature !If buoyancy enabled then liquid potential temperature, othrewise potential temperature.
  real(knd), allocatable, dimension(:,:,:)   :: Moisture    !Total specific humidity q_t.
  real(knd), allocatable, dimension(:,:,:,:) :: Scalar      !last index is a number of scalar (because of paging)

  real(knd) :: delta = 0

  real(tim) :: dt = 0 !time_step

  integer   :: time_step

  real(knd) :: time_step_time

  integer(dbl) :: time_steps_timer_count_start, &
                  time_steps_timer_count_1, &
                  time_steps_timer_count_2, &
                  timer_max_count
  
  integer(dbl) :: timer_count_time_limit

  logical :: error_exit = .false.


  call par_init

  call GetEndianness


  call system_clock(count = time_steps_timer_count_start)  
  
  
  call par_sync_out("Reading parameters...")
  call ReadConfiguration


  call par_sync_out("Setting up boundary conditions...")
  call InitBoundaryConditions


  call par_sync_out("Allocating arrays...")
  call AllocateGlobals

  call par_sync_out("Preparing output data and files...")
  call AllocateOutputs

  call par_sync_out("Setting up initial conditions...")
  call InitialConditions(U, V, W, Pr, &
                         Temperature ,Moisture, Scalar, &
                         time_stepping%dt)


  call InitTimer

  time_stepping%time = time_stepping%start_time
  time_step = 0

  call OutTStep(U, V, W, Pr, &
                Temperature, Moisture, Scalar, &
                time_stepping%dt, delta)
  
  init_phase = .false.
  run_phase = .true.

  if (time_stepping%end_time > time_stepping%start_time) then

    call par_sync_out("Computing...")

    time_step = 1

    do

      if (time_step > time_stepping%max_number_of_time_steps) exit

      call system_clock(count = time_steps_timer_count_1)

#ifdef PAR
      if (master) then
        if (enable_multiple_domains) then
          write (*,'(a,i0,a,i12,a,f12.6)') "domain: ", domain_index, &
                                           "   tstep: ", time_step, &
                                           "   time: ",  time_stepping%time
        else
          write (*,'(a,i12,a,f12.6)') "tstep: ", time_step, "   time: ", time_stepping%time
        end if
      end if
#else
      write (*,'(a,i12,a,f12.6)') "tstep: ", time_step, "   time: ", time_stepping%time
#endif



      call TMarchRK3(U, V, W, Pr, &
                     Temperature, Moisture, Scalar, &
                     delta)



      if (time_stepping%variable_time_steps) then
        time_stepping%time = time_stepping%time + time_stepping%dt
        if (abs(time_stepping%time - time_stepping%end_time)<5*spacing(time_stepping%end_time)) then
          time_stepping%time = time_stepping%end_time
        end if
      else
        time_stepping%time = time_stepping%start_time + time_step * time_stepping%dt
        if (abs(time_stepping%time - time_stepping%end_time)<time_stepping%dt/10) then
          time_stepping%time = time_stepping%end_time
        end if
      end if



      call OutTStep(U, V, W, Pr, &
                    Temperature, Moisture, Scalar, &
                    time_stepping%dt, delta)


      call system_clock(count = time_steps_timer_count_2)
      if (master) then
        time_step_time = real(time_steps_timer_count_2 - time_steps_timer_count_1, knd) / &
          real(timer_rate, knd)
        time_steps_time = time_steps_time + time_step_time
      end if

      if ((steady==1) .and. (delta<eps)) then
        if (master) write (*,*) "Steady state reached."
        exit
      endif

      if ((steady==0) .and. (time_stepping%time>=time_stepping%end_time)) then
        if (master) write (*,*) "Time limit reached."
        exit
      endif

      if (time_step>=3 .and. &
          time_stepping%variable_time_steps .and. &
          time_stepping%dt < time_stepping%dt_min) then
        if (master) write (*,*) "Solution diverged. 1"
        error_exit = .true.
        exit
      endif

      if (time_step>=3 .and. &
          .not.time_stepping%variable_time_steps .and. &
          time_stepping%enable_CFL_check .and. &
          time_stepping%CFL > time_stepping%CFL_max) then
        if (master) write (*,*) "Solution diverged. 2"
        error_exit = .true.
        exit
      endif
      
      !don't check that often to avoid synchronizing unnecesarilly
      if (mod(time_step,time_stepping%check_period)==0) then
        if (master) then
          error_exit = time_steps_timer_count_2 - time_steps_timer_count_start > timer_count_time_limit
          if (error_exit) write(*,*) "Maximum clock time exceeded."
        end if
#ifdef PAR        
        call par_co_broadcast(error_exit, master_im)
#endif        
        if (error_exit) exit
      end if

      time_step = time_step + 1
      
    enddo

  endif

  call par_sync_all
  if (master .and. time_steps_time>0) &
    write(*,'(1x,a,f0.2)') "Total wall clock time for time steps: ", time_steps_time
  if (master .and. poisson_solver_time>0) &
    write(*,'(1x,a,f0.2)') "Wall clock time for Poisson solver: ", poisson_solver_time
#ifdef PAR
  block
    use domains_bc_par
    if (master .and. enable_multiple_domains .and. time_communicating_domains>0) &
      write(*,'(1x,a,i0,a,f0.2)') "Wall clock time waiting and communicating with other domains on domain ", &
        domain_index,": ",time_communicating_domains
  end block
#endif
  call par_sync_all

  call Dynamics_Deallocate
  call Scalars_Deallocate

  call par_sync_out("Saving results...")


  call Output(U,V,W,Pr,Temperature,Moisture,Scalar)


  call par_sync_out("saved.")

  call DeallocateGlobals

#ifdef PAR
  if (error_exit .and. enable_multiple_domains .and. number_of_domains>1) then
    call error_stop("Error, one of the domains stopped.")
  else
    call par_finalize
  end if
#endif


contains


  subroutine AllocateGlobals
    use WaterThermodynamics
    allocate(U(-2:Unx+3,-2:Uny+3,-2:Unz+3))
    allocate(V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
    allocate(W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
    allocate(Pr(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))

    U = 0
    V = 0
    W = 0
    Pr = 0


    if (enable_buoyancy) then
      allocate(Temperature(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
      Temperature = 0
    else
      allocate(Temperature(0,0,0))
    endif

    if (enable_moisture) then
      allocate(Moisture(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
      Moisture = 0      
    else
      allocate(Moisture(0,0,0))
    endif
    
    if (enable_liquid) then
      allocate(LiquidWater(1:Prnx,1:Prny,1:Wnz+1))
      LiquidWater = 0      
    else
      allocate(LiquidWater(0,0,0))
    endif


    if (num_of_scalars>0) then
      allocate(Scalar(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,num_of_scalars))
      Scalar = 0
    else
      allocate(Scalar(0,0,0,0))
    endif


    allocate(Viscosity(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))

    if (enable_buoyancy .or. &
        enable_moisture .or. &
        num_of_scalars>0)      then

      allocate(TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
      TDiff = 0

    else
      allocate(TDiff(0,0,0))
    endif

  end subroutine AllocateGlobals


  subroutine DeallocateGlobals
    !To assist memory leak checking because some compilers do not
    !deallocate main program and module variables automatically
    !(in F2008 they have the SAVE attribute).

    deallocate(U, V, W, Pr)

    deallocate(Temperature, Moisture, Scalar)

    deallocate(Viscosity, TDiff)

    deallocate(xU, yV, zW, dxU, dyV, dzW)
    deallocate(xPr, yPr, zPr, dxPr, dyPr, dzPr)

    deallocate(Utype, Vtype, Wtype, Prtype)

    deallocate(Uin, Vin, Win)

    if (allocated(TempIn)) deallocate(TempIn)
    if (allocated(MoistIn)) deallocate(MoistIn)

    if (allocated(BsideTArr)) deallocate(BsideTArr)
    if (allocated(BsideTFlArr)) deallocate(BsideTFlArr)

    if (allocated(BsideMArr)) deallocate(BsideMArr)
    if (allocated(BsideMFlArr)) deallocate(BsideMFlArr)


  end subroutine
  
  
  subroutine InitTimer
    call system_clock(count_rate = timer_rate)
    call system_clock(count_max = timer_max_count)   
    timer_count_time_limit = int( min(time_stepping%clock_time_limit &
                                        * real(timer_rate, knd),  &
                                      real(timer_max_count, knd) * 0.999_dbl) &
                                , dbl)  
  end subroutine


end program CLMM

