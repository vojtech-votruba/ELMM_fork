module Custom_gabls4
  use Kinds
  use Interpolation
  use Parameters
#ifdef PAR  
  use custom_mpi
  use custom_par
#endif
  
  implicit none
  
  logical :: initialized = .false.
  
  integer :: last = 1
  real(knd) :: surf_temp
  
  real(knd), allocatable :: coefs(:,:), table(:,:)
  
  real(knd), allocatable, dimension(:) :: temp, tempfl, up, vp, uu, vv, ww, uw, vw, tkesgs
  
#ifdef PAR
  real(knd), allocatable, dimension(:) :: u_part, v_part, uu_part, vv_part, ww_part, tkesgs_part
  integer, allocatable, dimension(:) :: ns_z, offs_z
  
  interface
    subroutine MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, &
                          RECVTYPE, ROOT, COMM, IERROR)
       INTEGER    SENDBUF, RECVBUF(*)
       INTEGER    SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, ROOT
       INTEGER    COMM, IERROR
    end subroutine
    subroutine MPI_GATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNTS, &
                           DISPLS, RECVTYPE, ROOT, COMM, IERROR)
      use Kinds
      real(KND)    SENDBUF(*), RECVBUF(*)
      INTEGER    SENDCOUNT, SENDTYPE, RECVCOUNTS(*), DISPLS(*)
      INTEGER    RECVTYPE, ROOT, COMM, IERROR
    end subroutine
  end interface

#endif

  integer, parameter :: h2 = 1, h3 = 2, h7 = 3, h9 = 4, h10 = 5, h15 = 6, &
                        h18 = 7, h23 = 8, h25 = 9, h30 = 10, h33 = 11, &
                        h38 = 12, h42 = 13
  real(knd), parameter :: z(13) = [2._knd, 3.3_knd, 7.03_knd, 8.8_knd, 10._knd, 15.43_knd, &
                                   17.9_knd, 22.79_knd, 25.3_knd, 30.15_knd, 32.7_knd, &
                                   37.51_knd, 41.9_knd]
  integer :: zk(13), zwk(13)
  
  character(:), allocatable :: gabls_series_dir
  
  integer ::  gabls_series_step = 0
  
  integer, parameter :: gabls_series_max_length = 10000
  
  real(knd), dimension(gabls_series_max_length) :: t2, &
                                                   t3, &
                                                   u3, &
                                                   v3, &
                                                   uw3, &
                                                   vw3, &
                                                   Tw3, &
                                                   tke3, &
                                                   t9, &
                                                   u9, &
                                                   v9, &
                                                   u10, &
                                                   v10, &
                                                   uw15, &
                                                   vw15, &
                                                   Tw15, &
                                                   tke15, &
                                                   uw7, &
                                                   vw7, &
                                                   Tw7, &
                                                   tke7, &
                                                   t18, &
                                                   u18, &
                                                   v18, &
                                                   uw23, &
                                                   vw23, &
                                                   Tw23, &
                                                   tke23, &
                                                   t25, &
                                                   u25, &
                                                   v25, &
                                                   uw30, &
                                                   vw30, &
                                                   Tw30, &
                                                   tke30, &
                                                   t33, &
                                                   u33, &
                                                   v33, &
                                                   uw38, &
                                                   vw38, &
                                                   Tw38, &
                                                   tke38, &
                                                   t42, &
                                                   u42, &
                                                   v42, &
                                                   tsurf, &
                                                   hpbl
  
contains

  subroutine read_table(fname)
    character(*), intent(in) :: fname
    integer :: io, u, n, i
    real(knd) :: a
    character(80) :: line
    
    open(newunit=u, file=fname, status="old", action="read")
    read(u,*)
    n = 0
    do
      read(u,'(a)',iostat=io) line
      if (io/=0) exit
      n = n + 1
    end do

    allocate(table(n,2))

    rewind(u)
    read(u,*)
    
    do i = 1, n
      read(u,'(a)',iostat=io) line
      read(line, *) table(i,:)
    end do

    close(u)

  end subroutine
    
  subroutine initialize_surface_temp
  
    call read_table("input/Tg.txt")

    table(:,2) = table(:,2) * (1000.0_knd/651._knd)**0.2854_knd

    allocate(coefs(0:3,size(table,1)))
 
    call cubic_spline(table(:,1), table(:,2), coefs)
    
    surf_temp = cubic_spline_eval(time_stepping%start_time, table(:,1), coefs, last)

    initialized = .true.
  
  end subroutine
  
  pure function interp(z, zk, a) result(res)
    real(knd) :: res
    real(knd), intent(in) :: z
    integer, intent(in) :: zk
    real(knd), intent(in) :: a(:)
    
    
    res = a(zk) + &
            (a(zk+1) - a(zk)) * &
            (z - zPr(zk)) / (zPr(zk+1) - zPr(zk))
  end function

  pure function interpw(z, zk, a) result(res)
    real(knd) :: res
    real(knd), intent(in) :: z
    integer, intent(in) :: zk
    real(knd), intent(in) :: a(:)
    
    
    res = a(zk) + &
            (a(zk+1) - a(zk)) * &
            (z - zW(zk)) / (zW(zk+1) - zW(zk))
  end function
  
  subroutine init_series   
    integer :: i, k
    
#ifdef PAR    
    integer :: ierr
      
    if (master) then
    
      allocate(ns_z(nzims), offs_z(nzims))
      
    else if (iim==1 .and. jim==1) then
    
      allocate(ns_z(0), offs_z(0))
      
    end if
    
    if (iim==1 .and. jim==1) &
      call MPI_Gather(Prnz, 1, MPI_INTEGER, ns_z, 1, MPI_INTEGER, 0, comm_row_z, ierr)
      
    if (master) then
    
      offs_z(1) = 0
      do  k = 2, nzims
        offs_z(k) = offs_z(k-1) + ns_z(k-1)
      end do
      
      allocate(uu(sum(ns_z)))
      allocate(vv(sum(ns_z)))
      allocate(ww(sum(ns_z)))
      allocate(tkesgs(sum(ns_z)))
      
    else if (iim==1 .and. jim==1) then
    
      allocate(uu(0))
      allocate(vv(0))
      allocate(ww(sum(ns_z)))
      allocate(tkesgs(0))
      
    end if
    
    if (master) gabls_series_dir = "output/im-1-1-1/"
#else    
    gabls_series_dir = "output/"
#endif

    if (master) then
      k = 0
      do i = 1, size(z)
        call find_zk(z(i), zk(i))
      end do
      k = 0
      do i = 1, size(z)
        call find_zwk(z(i), zwk(i))
      end do
      
      if (any(zk>Prnz).or.any(zwk>Prnz)) &
        call error_stop( "Bottom image domain too shallow for GABLS time series.")
      
      if (any(zk<1).or.any(zwk<1)) &
        call error_stop( "Vertical resolution too coarse for GABLS time series.")
      
    
      call system("mkdir -p "//gabls_series_dir)
      
      call create("t2.unf")
      call create("t3.unf")
      call create("u3.unf")
      call create("v3.unf")
      call create("uw3.unf")
      call create("vw3.unf")
      call create("Tw3.unf")
      call create("tke3.unf")
      call create("t9.unf")
      call create("u9.unf")
      call create("v9.unf")
      call create("u10.unf")
      call create("v10.unf")
      call create("uw7.unf")
      call create("vw7.unf")
      call create("Tw7.unf")
      call create("tke7.unf")
      call create("uw15.unf")
      call create("vw15.unf")
      call create("Tw15.unf")
      call create("tke15.unf")
      call create("t18.unf")
      call create("u18.unf")
      call create("v18.unf")
      call create("uw23.unf")
      call create("vw23.unf")
      call create("Tw23.unf")
      call create("tke23.unf")
      call create("t25.unf")
      call create("u25.unf")
      call create("v25.unf")
      call create("uw30.unf")
      call create("vw30.unf")
      call create("Tw30.unf")
      call create("tke30.unf")
      call create("t33.unf")
      call create("u33.unf")
      call create("v33.unf")
      call create("uw38.unf")
      call create("vw38.unf")
      call create("Tw38.unf")
      call create("tke38.unf")
      call create("t42.unf")
      call create("u42.unf")
      call create("v42.unf")
      call create("tsurf.unf")
      call create("hpbl.unf")
    end if

    
  contains
  
    subroutine find_zk(zval, zkval)
      real(knd), intent(in) :: zval
      integer, intent(out) :: zkval
      do while (dzmin/2 + k * dzmin < zval .and. k <= Prnz+2)
        k = k +1
      end do
      zkval = k
    end subroutine
    
    subroutine find_zwk(zval, zkval)
      real(knd), intent(in) :: zval
      integer, intent(out) :: zkval
      do while (k*dzmin < zval .and. k <= Prnz+2)
        k = k + 1
      end do
      zkval = k
    end subroutine
    
    subroutine create(fname)
      character(*) :: fname
      integer :: unit
      open(newunit=unit,file=gabls_series_dir//fname,status="replace")
      close(unit)
    end subroutine  
  end subroutine init_series 
  
  subroutine OutputGablsSeries
    use Freeunit
    
    integer :: unit
    
    call newunit(unit)
    
    call save(t2, "t2.unf")
    call save(t3, "t3.unf")
    call save(u3, "u3.unf")
    call save(v3, "v3.unf")
    call save(uw3, "uw3.unf")
    call save(vw3, "vw3.unf")
    call save(Tw3, "Tw3.unf")
    call save(tke3, "tke3.unf")
    call save(t9, "t9.unf")
    call save(u9, "u9.unf")
    call save(v9, "v9.unf")
    call save(u10, "u10.unf")
    call save(v10, "v10.unf")
    call save(uw7, "uw7.unf")
    call save(vw7, "vw7.unf")
    call save(Tw7, "Tw7.unf")
    call save(tke7, "tke7.unf")
    call save(uw15, "uw15.unf")
    call save(vw15, "vw15.unf")
    call save(Tw15, "Tw15.unf")
    call save(tke15, "tke15.unf")
    call save(t18, "t18.unf")
    call save(u18, "u18.unf")
    call save(v18, "v18.unf")
    call save(uw23, "uw23.unf")
    call save(vw23, "vw23.unf")
    call save(Tw23, "Tw23.unf")
    call save(tke23, "tke23.unf")
    call save(t25, "t25.unf")
    call save(u25, "u25.unf")
    call save(v25, "v25.unf")
    call save(uw30, "uw30.unf")
    call save(vw30, "vw30.unf")
    call save(Tw30, "Tw30.unf")
    call save(tke30, "tke30.unf")
    call save(t33, "t33.unf")
    call save(u33, "u33.unf")
    call save(v33, "v33.unf")
    call save(uw38, "uw38.unf")
    call save(vw38, "vw38.unf")
    call save(Tw38, "Tw38.unf")
    call save(tke38, "tke38.unf")
    call save(t42, "t42.unf")
    call save(u42, "u42.unf")
    call save(v42, "v42.unf")
    call save(tsurf, "tsurf.unf")
    call save(hpbl, "hpbl.unf")
  contains
    subroutine save(a, fname)
      real(knd), intent(in), contiguous :: a(:)
      character(*), intent(in) :: fname
      
      open(unit,file=gabls_series_dir//fname, &
         access="stream",status="old",position="append")
      write(unit) a(1:gabls_series_step)
      close(unit)
    end subroutine
  end subroutine OutputGablsSeries

end module Custom_gabls4


subroutine CustomTimeStepOutput
  use Custom_gabls4
  use Parameters
  use Outputs
  
  implicit none

  integer :: step
  integer :: ierr
  
  real(knd) :: tke_k, tke_last, tke_01, z90
  integer :: k, n
  
  gabls_series_step = gabls_series_step + 1

#ifdef PAR
  if (kim==1) then
    temp = current_profiles%temp(1:Prnz)
    tempfl = current_profiles%tempfl(1:Prnz)
    tempfl = tempfl + current_profiles%tempflsgs(1:Prnz)
    uw = current_profiles%uw(1:Prnz)
    vw = current_profiles%vw(1:Prnz)
    uw = uw + current_profiles%uwsgs(1:Prnz)
    vw = vw + current_profiles%vwsgs(1:Prnz)
    
    call par_sum_to_master_horizontal(temp)
    call par_sum_to_master_horizontal(tempfl)
    call par_sum_to_master_horizontal(uw)
    call par_sum_to_master_horizontal(vw)
  end if

  n = gPrnx * gPrny
  
  u_part = current_profiles%u(1:Prnz)
  v_part = current_profiles%v(1:Prnz)
  call par_sum_to_master_horizontal(u_part)
  call par_sum_to_master_horizontal(v_part)
  uu_part = current_profiles%uu(1:Prnz)
  vv_part = current_profiles%vv(1:Prnz)
  ww_part = current_profiles%ww(1:Prnz)

  tkesgs_part = current_profiles%tkesgs(1:Prnz)
  
  call par_sum_to_master_horizontal(uu_part)
  call par_sum_to_master_horizontal(vv_part)
  call par_sum_to_master_horizontal(ww_part)
  call par_sum_to_master_horizontal(tkesgs_part)
  
  if (iim==1.and.jim==1) then
    u_part = u_part / n
    v_part = v_part / n
    uu_part = uu_part / n - u_part**2
    vv_part = vv_part / n - v_part**2
  
    call MPI_Gatherv(uu_part, Prnz, MPI_KND, uu, ns_z, offs_z, MPI_KND, 0, comm_row_z, ierr)
    call MPI_Gatherv(vv_part, Prnz, MPI_KND, vv, ns_z, offs_z, MPI_KND, 0, comm_row_z, ierr)
    call MPI_Gatherv(ww_part, Prnz, MPI_KND, ww, ns_z, offs_z, MPI_KND, 0, comm_row_z, ierr)
    call MPI_Gatherv(tkesgs_part, Prnz, MPI_KND, tkesgs, ns_z, offs_z, MPI_KND, 0, comm_row_z, ierr)
  end if
  
  if (master) then
    temp = temp / n
    tempfl = tempfl / n
    up = u_part
    vp = v_part
    ww = ww / n
    uw = uw / n
    vw = vw / n
    tkesgs = tkesgs / n
  
#else
  n = Prnx * Prny
  temp = current_profiles%temp(1:Prnz) / n
  tempfl = current_profiles%tempfl(1:Prnz) / n
  tempfl = tempfl + current_profiles%tempflsgs(1:Prnz) / n
  up = current_profiles%u(1:Prnz) / n 
  vp = current_profiles%v(1:Prnz) / n 
  uu = current_profiles%uu(1:Prnz)
  vv = current_profiles%vv(1:Prnz)
  ww = current_profiles%ww(1:Prnz)

  uu = uu / n
  vv = vv / n
  ww = ww / n
  uu = uu - up**2
  vv = vv - vp**2
  uw = current_profiles%uw(1:Prnz) / n 
  vw = current_profiles%vw(1:Prnz) / n 
  uw = uw + current_profiles%uwsgs(1:Prnz) / n 
  vw = vw + current_profiles%vwsgs(1:Prnz) / n 
  tkesgs = current_profiles%tkesgs(1:Prnz) / n 
#endif
  
    step = gabls_series_step

    t2(step) = interp(z(h2), zk(h2), temp)

    t3(step) = interp(z(h3), zk(h3), temp)
    u3(step) = interp(z(h3), zk(h3), up)
    v3(step) = interp(z(h3), zk(h3), vp)
    uw3(step) = interpw(z(h3), zwk(h3), uw)
    vw3(step) = interpw(z(h3), zwk(h3), vw)
    Tw3(step) = interpw(z(h3), zwk(h3), tempfl)
    tke3(step) = ( interp(z(h3), zk(h3), uu) + &
                   interp(z(h3), zk(h3), vv) + &
                   interpw(z(h3), zwk(h3), ww)) / 2 + &
                 interp(z(h3), zk(h3), tkesgs)

    t9(step) = interp(z(h9), zk(h9), temp)
    u9(step) = interp(z(h9), zk(h9), up)
    v9(step) = interp(z(h9), zk(h9), vp)

    u10(step) = interp(z(h10), zk(h10), up)
    v10(step) = interp(z(h10), zk(h10), vp)

    uw7(step) = interpw(z(h7), zwk(h7), uw)
    vw7(step) = interpw(z(h7), zwk(h7), vw)
    Tw7(step) = interpw(z(h7), zwk(h7), tempfl)
    tke7(step) = ( interp(z(h7), zk(h7), uu) + &
                    interp(z(h7), zk(h7), vv) + &
                    interpw(z(h7), zwk(h7), ww)) / 2 + &
                  interp(z(h7), zk(h7), tkesgs)

    uw15(step) = interpw(z(h15), zwk(h15), uw)
    vw15(step) = interpw(z(h15), zwk(h15), vw)
    Tw15(step) = interpw(z(h15), zwk(h15), tempfl)
    tke15(step) = ( interp(z(h15), zk(h15), uu) + &
                    interp(z(h15), zk(h15), vv) + &
                    interpw(z(h15), zwk(h15), ww)) / 2 + &
                  interp(z(h15), zk(h15), tkesgs)

    t18(step) = interp(z(h18), zk(h18), temp)
    u18(step) = interp(z(h18), zk(h18), up)
    v18(step) = interp(z(h18), zk(h18), vp)

    uw23(step) = interpw(z(h23), zwk(h23), uw)
    vw23(step) = interpw(z(h23), zwk(h23), vw)
    Tw23(step) = interpw(z(h23), zwk(h23), tempfl)
    tke23(step) = ( interp(z(h23), zk(h23), uu) + &
                    interp(z(h23), zk(h23), vv) + &
                    interpw(z(h23), zwk(h23), ww)) / 2 + &
                  interp(z(h23), zk(h23), tkesgs)

    t25(step) = interp(z(h25), zk(h25), temp)
    u25(step) = interp(z(h25), zk(h25), up)
    v25(step) = interp(z(h25), zk(h25), vp)

    uw30(step) = interpw(z(h30), zwk(h30), uw)
    vw30(step) = interpw(z(h30), zwk(h30), vw)
    Tw30(step) = interpw(z(h30), zwk(h30), tempfl)
    tke30(step) = ( interp(z(h30), zk(h30), uu) + &
                    interp(z(h30), zk(h30), vv) + &
                    interpw(z(h30), zwk(h30), ww)) / 2 + &
                  interp(z(h30), zk(h30), tkesgs)

    t33(step) = interp(z(h33), zk(h33), temp)
    u33(step) = interp(z(h33), zk(h33), up)
    v33(step) = interp(z(h33), zk(h33), vp)

    uw38(step) = interpw(z(h38), zwk(h38), uw)
    vw38(step) = interpw(z(h38), zwk(h38), vw)
    Tw38(step) = (interpw(z(h38), zwk(h38), tempfl))
    tke38(step) = ( interp(z(h38), zk(h38), uu) + &
                    interp(z(h38), zk(h38), vv) + &
                    interpw(z(h38), zwk(h38), ww) ) / 2 + &
                  interp(z(h38), zk(h38), tkesgs)

    t42(step) = interp(z(h42), zk(h42), temp)
    u42(step) = interp(z(h42), zk(h42), up)
    v42(step) = interp(z(h42), zk(h42), vp)
    
    tsurf(step) = surf_temp
    
    hpbl(step) = gzmax
    
    tke_01 = tke3(step) / 10
    tke_last = tke3(step)
    do k = 3, size(tkesgs)-1
      tke_k = ( uu(k) + vv(k) + (ww(k) + ww(k-1)) / 2 ) &
                / 2 + &
              tkesgs(k)
      if (tke_k<tke_01) then
        z90 = dzmin*(k-0.5_knd) + dzmin * (tke_last - tke_01) / (tke_last - tke_k)
        hpbl(step) = z90 / 0.9
        exit
      end if
      
      tke_last = tke_k
    end do

#ifdef PAR
  end if
#endif

  if (gabls_series_step==gabls_series_max_length) then
    if (master) call OutputGablsSeries
    gabls_series_step = 0
  end if

end subroutine CustomTimeStepOutput


subroutine CustomOutput
  use Custom_gabls4
  
  if (master) call OutputGablsSeries
end subroutine









function CustomSurfaceTemperature(x,y,z,t) result(res)
   use Kinds
   use Custom_gabls4, only: surf_temp
   
   implicit none
   
   real(knd) :: res
   real(knd), intent(in) :: x, y, z
   real(tim), intent(in) :: t
   
   res = surf_temp

end function

subroutine CustomTimeStepProcedure
  use Custom_gabls4
  use Parameters
  
  implicit none

  surf_temp = cubic_spline_eval(time_stepping%time, table(:,1), coefs, last)
end subroutine








subroutine CustomConfiguration_first
  use Custom_gabls4

  call initialize_surface_temp
  
end subroutine

subroutine CustomConfiguration_last
  use Kinds
  use Scalars
  use Custom_gabls4
  use Frames_common
  use VTKFrames
  use Sponge, only: enable_top_sponge, top_sponge_bottom, sponge_to_profiles
  
  implicit none
  
  integer :: n, u, i, io
  real(knd), allocatable :: theta(:,:)
  character(80) :: line
  
  call init_series
  
  open(newunit=u, file="input/theta_profile.txt", status="old", action="read")
  n = 0
  do
    read(u,'(a)',iostat=io) line
    if (io/=0) exit
    n = n + 1
  end do

  allocate(theta(n,2))

  rewind(u)
  
  do i = 1, n
    read(u,'(a)',iostat=io) line
    read(line, *) theta(i,:)
  end do

  close(u)
  
  call move_alloc(theta, TemperatureProfileObj%points)
  
!   call AddDomain(TFrameDomain('3d', &
!                3, 1, 0.0_knd, &
!                TFrameTimes(3600._knd, 82800._knd, 12), TFrameFlags(U=1, Pr=1, temperature=1)))
               
  enable_top_sponge = .true.
  top_sponge_bottom = 600

  sponge_to_profiles = .true.
  
end subroutine

subroutine CustomBoundaryConditions
  use Parameters
  use Wallmodels
  use Sponge, only: U_sponge_avg, V_sponge_avg, W_sponge_avg
  
  implicit none
  
  integer :: i, j
  
  WMPoints%z0 = 0.001
  WMPoints%z0H = 0.0001
  
  do j = 1, 3
    do i = 1, 6
      WMPointsUVW(i,j)%points%z0 = 0.001
      WMPointsUVW(i,j)%points%z0H = 0.0001
    end do
  end do

  U_sponge_avg = Uin(1,1:Unz)
  V_sponge_avg = Vin(1,1:Vnz)
  W_sponge_avg = Win(1,1:Wnz)
end subroutine

subroutine  CustomSolidBodies
end subroutine  CustomSolidBodies
