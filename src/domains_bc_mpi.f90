module domains_bc_par

  use custom_par
  use domains_bc_types
  use domains_bc_init
  use Parameters

  implicit none

  private

  public par_init_domain_boundary_conditions, &
         par_exchange_domain_bounds, &
         par_update_domain_bounds, &
         par_update_pr_gradient, &
         par_update_domain_bounds_UVW, &
         par_update_domain_bounds_temperature, &
         par_update_domain_bounds_moisture, &
         par_domain_bound_relaxation, &
         par_receive_initial_conditions, &
         par_send_initial_conditions, &
         par_domain_two_way_nesting_feedback, &
         domain_spatial_ratio, &
         domain_time_step_ratio, &
         time_communicating_domains




  !send buffers should be indexed from 1, there may be more of them from several nested domains
  !they are not accessed by the boundary conditions routines
  type(dom_bc_buffer_copy), allocatable :: domain_bc_send_buffers(:)
  !receive buffers should be indexed with the index of the domain side (West to Top)
  ! to ease finding the right buffer to a given nested boundary 
  type(dom_bc_buffer_copy), allocatable :: domain_bc_recv_buffers_copy(:)

  !receive buffers should be indexed with the index of the domain side (West to Top)
  ! to ease finding the right buffer to a given nested boundary 
  type(dom_bc_buffer_turbulence_generator), allocatable :: domain_bc_recv_buffers(:)
  
  real(dbl), protected :: time_communicating_domains = 0

contains


  subroutine par_init_domain_boundary_conditions()

    call par_init_domain_boundary_conditions_implementation(domain_bc_send_buffers, &
                                                            domain_bc_recv_buffers)
  end subroutine



  subroutine par_exchange_domain_bounds(U, V, W, Temperature, Moisture, Scalar, time, dt, receive, send)
    !if necessary send or receive the new boundary conditions
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V ,W 
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: Temperature, Moisture
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(in) :: Scalar
    real(tim), intent(in) :: time, dt
    logical, intent(in), optional :: send, receive

    logical :: send_l, receive_l

    integer, allocatable :: requests(:)
    integer :: di, i, j, k, scal
    integer :: ie
    integer(int64) :: t1, t2, timer_rate
    
    call system_clock(count=t1, count_rate=timer_rate)

    send_l = .true.
    receive_l = .true.

    if (present(send)) send_l = send
    if (present(receive)) receive_l = receive

    allocate(requests(0))

    if (enable_multiple_domains) then

      if (send_l) then

        if (allocated(domain_bc_send_buffers)) then
          do i = lbound(domain_bc_send_buffers,1), &
                 ubound(domain_bc_send_buffers,1)
            associate(b => domain_bc_send_buffers(i))

              if (b%enabled) then
                !derivatives approximated by the backward difference

                if (b%time<=-tiny(-1.0_tim)) then
                  b%dU_dt = 0
                  b%dV_dt = 0
                  b%dW_dt = 0
                  if (enable_buoyancy) then
                    b%dTemperature_dt = 0
                  end if
                  if (enable_moisture) then
                    b%dMoisture_dt = 0
                  end if
                  if (num_of_scalars > 0) then
                    b%dScalar_dt = 0
                  end if
                else
                  b%dU_dt(:,:,:) = (U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2) - b%U) / dt
                  b%dV_dt(:,:,:) = (V(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2) - b%V) / dt
                  b%dW_dt(:,:,:) = (W(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2) - b%W) / dt
                  if (enable_buoyancy) then
                    b%dTemperature_dt(:,:,:) = (Temperature(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2) - b%Temperature) / dt
                  end if
                  if (enable_moisture) then
                    b%dMoisture_dt(:,:,:) = (Moisture(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2) - b%Moisture) / dt
                  end if
                  if (num_of_scalars > 0) then
                    do scal = 1, num_of_scalars
                      b%dScalar_dt(:,:,:,scal) = (Scalar(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2,scal)-b%Scalar(:,:,:,scal)) / dt
                    end do
                  end if
                end if

                b%U(:,:,:) = U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2)
                b%V(:,:,:) = V(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2)
                b%W(:,:,:) = W(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2)
                if (enable_buoyancy) then
                  b%Temperature(:,:,:) = Temperature(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2)
                end if
                if (enable_moisture) then
                  b%Moisture(:,:,:) = Moisture(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2)
                end if
                if (num_of_scalars > 0) then
                  do scal = 1, num_of_scalars
                    b%Scalar(:,:,:,scal) = Scalar(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2,scal)
                  end do
                end if

                b%time = time

                call send_arrays(b)
              end if

            end associate
          end do
        end if

        
        if (allocated(domain_parent_buffers)) then
          do di = 1, size(domain_parent_buffers)
            if (allocated(domain_parent_buffers(di)%bs)) then
              do k = lbound(domain_parent_buffers(di)%bs,3), ubound(domain_parent_buffers(di)%bs,3)
                do j = lbound(domain_parent_buffers(di)%bs,2), ubound(domain_parent_buffers(di)%bs,2)
                  do i = lbound(domain_parent_buffers(di)%bs,1), ubound(domain_parent_buffers(di)%bs,1)

                    associate(b => domain_parent_buffers(di)%bs(i,j,k))
                      if (b%exchange_pr_gradient_x) then
                        b%pr_gradient_x_dt = (pr_gradient_x - b%pr_gradient_x) / dt
                        b%pr_gradient_x = pr_gradient_x
                        call send_scalar(b, 11, b%pr_gradient_x)
                        call send_scalar(b, 12, b%pr_gradient_x_dt)
                      end if

                      if (b%exchange_pr_gradient_y) then
                        b%pr_gradient_y_dt = (pr_gradient_y - b%pr_gradient_y) / dt
                        b%pr_gradient_y = pr_gradient_y
                        call send_scalar(b, 13, b%pr_gradient_y)
                        call send_scalar(b, 14, b%pr_gradient_y_dt)
                      end if

                      if (b%exchange_pr_gradient_z) then
                        b%pr_gradient_z_dt = (pr_gradient_z - b%pr_gradient_z) / dt
                        b%pr_gradient_z = pr_gradient_z
                        call send_scalar(b, 15, b%pr_gradient_z)
                        call send_scalar(b, 16, b%pr_gradient_z_dt)
                      end if
                    end associate

                  end do
                end do
              end do
            end if
          end do
        end if

      end if

      if (receive_l) then

        if (allocated(domain_bc_recv_buffers_copy)) then
          do i = lbound(domain_bc_recv_buffers_copy,1), &
                 ubound(domain_bc_recv_buffers_copy,1)
            associate(b => domain_bc_recv_buffers_copy(i))
    
              if (b%enabled) then
                call recv_arrays_copy(b)
                b%time = time
              end if

            end associate
          end do
        end if

        if (parent_domain>0) then
          associate (b => domain_child_buffer)
            if (b%time_step<=1) then
              if (b%exchange_pr_gradient_x) then
                call recv_scalar(b, 11, b%pr_gradient_x)
                call recv_scalar(b, 12, b%pr_gradient_x_dt)
              end if

              if (b%exchange_pr_gradient_y) then
                call recv_scalar(b, 13, b%pr_gradient_y)
                call recv_scalar(b, 14, b%pr_gradient_y_dt)
              end if

              if (b%exchange_pr_gradient_z) then
                call recv_scalar(b, 15, b%pr_gradient_z)
                call recv_scalar(b, 16, b%pr_gradient_z_dt)
              end if

              b%time = time
            end if
            if (b%time_step==b%time_step_ratio) then
              b%time_step=1
            else
              b%time_step = b%time_step + 1
            end if
          end associate
        end if


        if (allocated(domain_bc_recv_buffers)) then
          do i = lbound(domain_bc_recv_buffers,1), &
                 ubound(domain_bc_recv_buffers,1)
            associate(b => domain_bc_recv_buffers(i))
    
              if (b%enabled) then

                b%interpolate = .false.

                !for 0 and 1, initialization and the first time-step
                if (b%time_step<=1) then
                  call recv_arrays(b)

                  b%time = time
                  b%interpolate = .true.
                end if
                if (b%time_step==b%time_step_ratio) then
                  b%time_step=1
                else
                  b%time_step = b%time_step + 1
                end if
              end if

            end associate
          end do
        end if


        call MPI_Waitall(size(requests), requests, MPI_STATUSES_IGNORE, ie)
        if (ie/=0) call error_stop("Error, MPI_Waitall in par_exchange_domain_bounds returns", ie)

        if (receive_l) then
          if (allocated(domain_bc_recv_buffers)) then
            do i = lbound(domain_bc_recv_buffers,1), &
                   ubound(domain_bc_recv_buffers,1)
              associate(b => domain_bc_recv_buffers(i))
      
                if (b%enabled.and.b%interpolate) then
                  call par_interpolate_buffers(b)
                end if

                if (b%turb_generator_enabled) then
                  if (b%interpolate) call b%compute_sgs_tke
                  call b%turb_generator%time_step(b%U_turb, b%V_turb, b%V_turb, dt)

                end if

              end associate
            end do
          end if
       end if

      end if
      
    end if

    call system_clock(count=t2)
    time_communicating_domains = time_communicating_domains + real(t2-t1, dbl)/real(timer_rate,dbl)
    
  contains

    subroutine send_arrays(b)
      type(dom_bc_buffer_copy), intent(inout) :: b
      
      call send_array(b%U, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 1)
      call send_array(b%V, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 2)
      call send_array(b%W, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 3)

      call send_array(b%dU_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 1)
      call send_array(b%dV_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 2)
      call send_array(b%dW_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 3)

      if (enable_buoyancy) then
        call send_array(b%Temperature, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 4)
        call send_array(b%dTemperature_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 4)
      end if

      if (enable_moisture) then
        call send_array(b%Moisture, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 5)
        call send_array(b%dMoisture_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 5)
      end if

      if (num_of_scalars > 0) then
        call send_array4(b%Scalar, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 6)
        call send_array4(b%dScalar_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 6)
      end if

    end subroutine

    subroutine recv_arrays_copy(b)
      type(dom_bc_buffer_copy), intent(inout) :: b
      real(knd) :: avg
    
      call recv_array(b%U, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv_array(b%V, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv_array(b%W, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)    

      call recv_array(b%dU_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv_array(b%dV_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv_array(b%dW_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)

      if (enable_buoyancy) then
        call recv_array(b%Temperature, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
        call recv_array(b%dTemperature_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
      end if

      if (enable_moisture) then
        call recv_array(b%Moisture, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
        call recv_array(b%dMoisture_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
      end if

      if (num_of_scalars > 0) then
        call recv_array4(b%Scalar, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 6)
        call recv_array4(b%dScalar_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 6)
      end if
    end subroutine

    subroutine recv_arrays(b)
      class(dom_bc_buffer_refined), intent(inout) :: b
      real(knd) :: avg
    
      call recv_array(b%r_U, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv_array(b%r_V, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv_array(b%r_W, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)    

      call recv_array(b%r_dU_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv_array(b%r_dV_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv_array(b%r_dW_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)

      if (enable_buoyancy) then
        call recv_array(b%r_Temperature, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
        call recv_array(b%r_dTemperature_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
      end if

      if (enable_moisture) then
        call recv_array(b%r_Moisture, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
        call recv_array(b%r_dMoisture_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
      end if

      if (num_of_scalars > 0) then
        call recv_array4(b%r_Scalar, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 6)
        call recv_array4(b%r_dScalar_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 6)
      end if

    end subroutine

    subroutine send_array(a, rank, comm, tag)
      real(knd), intent(in), contiguous :: a(:,:,:)
      integer, intent(in) :: rank, comm, tag
      integer :: request, ie

      !ISsend because we do want synchronisation here. Otherwise a fast domain
      !could run too quickly, buffer too many requests and fill up the memory
      call MPI_ISsend(a, size(a), MPI_KND, rank, mod(tag, MPI_TAG_MAX), comm, request, ie)
      if (ie/=0) call error_stop("error sending MPI message.")

      requests = [requests, request]
    end subroutine

    subroutine send_array4(a, rank, comm, tag)
      real(knd), intent(in), contiguous :: a(:,:,:,:)
      integer, intent(in) :: rank, comm, tag
      integer :: request, ie

      !ISsend because we do want synchronisation here. Otherwise a fast domain
      !could run too quickly, buffer too many requests and fill up the memory
      call MPI_ISsend(a, size(a), MPI_KND, rank, tag, comm, request, ie)
      if (ie/=0) call error_stop("error sending MPI message.")

      requests = [requests, request]
    end subroutine

    subroutine recv_array(a, rank, comm, tag)
      real(knd), intent(inout), contiguous :: a(:,:,:)
      integer, intent(in) :: rank, comm, tag
      integer :: request, ie

      call MPI_IRecv(a, size(a), MPI_KND, rank, mod(tag, MPI_TAG_MAX), comm, request, ie)
      if (ie/=0) call error_stop("error receiving MPI message.")

      requests = [requests, request]
    end subroutine

    subroutine recv_array4(a, rank, comm, tag)
      real(knd), intent(inout), contiguous :: a(:,:,:,:)
      integer, intent(in) :: rank, comm, tag
      integer :: request, ie

      call MPI_IRecv(a, size(a), MPI_KND, rank, tag, comm, request, ie)
      if (ie/=0) call error_stop("error receiving MPI message.")

      requests = [requests, request]
    end subroutine

    subroutine send_scalar(b, tag_base, a)
      class(dom_buffer_pr_gradient), intent(inout) :: b
      integer, intent(in) :: tag_base
      real(knd), intent(in) :: a
      integer :: request, ie
      
      call MPI_ISsend(a, 1, MPI_KND, &
                     b%remote_rank, &
                     mod(tag_base, MPI_TAG_MAX), &
                     b%comm, &
                     request, ie)
      if (ie/=0) call error_stop("error sending MPI message.")

      requests = [requests, request]
    end subroutine

    subroutine recv_scalar(b, tag_base, a)
      class(dom_buffer_pr_gradient), intent(inout) :: b
      integer, intent(in) :: tag_base
      real(knd), intent(inout) :: a
      integer :: request, ie
      
      call MPI_IRecv(a, 1, MPI_KND, &
                      b%remote_rank, &
                      mod(tag_base, MPI_TAG_MAX), &
                      b%comm, &
                      request, ie)
      if (ie/=0) call error_stop("error sending MPI message.")

      requests = [requests, request]
    end subroutine

  end subroutine par_exchange_domain_bounds




  subroutine par_interpolate_buffers(b)
    class(dom_bc_buffer_refined), intent(inout) :: b
    integer :: scal

    if (any(b%interp_order==0)) then

      call interpolate_spectral(b)
      call interpolate_U_trilinear(b%r_dU_dt, b%dU_dt)
      call interpolate_V_trilinear(b%r_dV_dt, b%dV_dt)
      call interpolate_W_trilinear(b%r_dW_dt, b%dW_dt)

      !HACK
      if (b%direction/=To .and. kim==1 .and. Btype(Bo)==BC_NOSLIP) then
        block
          integer :: lo, k
          lo = b%spatial_ratio / 2 + mod(b%spatial_ratio,2)

          do k = lo-1,  1, -1
            b%U(:,:,k) = 2* b%U(:,:,k+1) - b%U(:,:,k+2)
          end do
          do k = b%Uk1,  0
            b%U(:,:,k) = - b%U(:,:,1-k)
          end do
        end block
      end if

    else !(all(b%interp_order==2)) then
      call interpolate_U_trilinear(b%r_U, b%U, .true.)
      call interpolate_U_trilinear(b%r_dU_dt, b%dU_dt)

      call interpolate_V_trilinear(b%r_V, b%V)
      call interpolate_V_trilinear(b%r_dV_dt, b%dV_dt)

      call interpolate_W_trilinear(b%r_W, b%W)
      call interpolate_W_trilinear(b%r_dW_dt, b%dW_dt)

      if (enable_buoyancy) then
        call interpolate_scalar_trilinear(b%r_Temperature, b%Temperature)
        call interpolate_scalar_trilinear(b%r_dTemperature_dt, b%dTemperature_dt)
      end if

      if (enable_moisture) then
        call interpolate_scalar_trilinear(b%r_Moisture, b%Moisture)
        call interpolate_scalar_trilinear(b%r_dMoisture_dt, b%dMoisture_dt)
      end if

      do scal = 1, num_of_scalars
        call interpolate_scalar_trilinear(b%r_Scalar(:,:,:,1), b%Scalar(:,:,:,1))
        call interpolate_scalar_trilinear(b%r_dScalar_dt(:,:,:,1), b%dScalar_dt(:,:,:,1))
      end do

!     else
!       call interpolate_U_spline(b%r_U, b%U)
!       call interpolate_U_spline(b%r_dU_dt, b%dU_dt)
! 
!       call interpolate_V_spline(b%r_V, b%V)
!       call interpolate_V_spline(b%r_dV_dt, b%dV_dt)
! 
!       call interpolate_W_spline(b%r_W, b%W)
!       call interpolate_W_spline(b%r_dW_dt, b%dW_dt)
    end if

  contains

!     NOTE: the following implementation of bspline interpolation is working, but only in double precision.
!     It also should be made faster.
!
!     subroutine interpolate_U_spline(in, out)
!       use bspline_oo_module
!       real(knd), intent(in), allocatable :: in(:,:,:)
!       real(knd), intent(inout), allocatable :: out(:,:,:)
!       integer :: i, j, k
!       type(bspline_3d) :: spl
!       integer :: iflag, idx, idy, idz
!       real(real64) :: val
! 
!       idx = 0; idy = 0; idz = 0
!       
!       call spl%initialize(b%r_xU(b%r_Ui1:b%r_Ui2), &
!                           b%r_y(b%r_Uj1:b%r_Uj2), &
!                           b%r_z(b%r_Uk1:b%r_Uk2), &
!                           in, &
!                           b%interp_order(1), b%interp_order(2), b%interp_order(3), &
!                           iflag)
!       if (iflag/=1) then
!         write(*,*) "Error in spline initialization, iflag:",iflag
!         stop
!       end if
!       
!       
!       do k = b%Uk1, b%Uk2
!         do j = b%Uj1, b%Uj2
!           do i = b%Ui1, b%Ui2
!             call spl%evaluate(xU(i), yPr(j), zPr(k), &
!                               idx, idy, idz, val, iflag)
!             if (iflag/=0) then
!               write(*,*) "Error in spline interpolation of U at", &
!                           xU(i), yPr(j), zPr(k)
!               write(*,*) "flag:",iflag
!               write(*,*) "xs:",b%r_xU(b%r_Ui1:b%r_Ui2)
!               write(*,*) "ys:",b%r_y(b%r_Uj1:b%r_Uj2)
!               write(*,*) "zs:",b%r_z(b%r_Uk1:b%r_Uk2)
!               call error_stop()
!             end if
! 
!             out(i,j,k) = val
!           end do
!         end do
!       end do
! 
!       call spl%destroy()
!     end subroutine
! 
!     subroutine interpolate_V_spline(in, out)
!       use bspline_oo_module
!       real(knd), intent(in), allocatable :: in(:,:,:)
!       real(knd), intent(inout), allocatable :: out(:,:,:)
!       integer :: i, j, k
!       type(bspline_3d) :: spl
!       integer :: iflag, idx, idy, idz
!       real(real64) :: val
! 
!       idx = 0; idy = 0; idz = 0
!       
!       call spl%initialize(b%r_x(b%r_Vi1:b%r_Vi2), &
!                           b%r_yV(b%r_Vj1:b%r_Vj2), &
!                           b%r_z(b%r_Vk1:b%r_Vk2), &
!                           in, &
!                           b%interp_order(1), b%interp_order(2), b%interp_order(3), &
!                           iflag)
!       if (iflag/=1) then
!         write(*,*) "Error in spline initialization, iflag:",iflag
!         stop
!       end if
!       
!       
!       do k = b%Vk1, b%Vk2
!         do j = b%Vj1, b%Vj2
!           do i = b%Vi1, b%Vi2
!             call spl%evaluate(xPr(i), yV(j), zPr(k), &
!                               idx, idy, idz, val, iflag)
!             if (iflag/=0) then
!               write(*,*) "Error in spline interpolation of V at", &
!                           xPr(i), yV(j), zPr(k)
!               write(*,*) "flag:",iflag
!               write(*,*) "xs:",b%r_x(b%r_Vi1:b%r_Vi2)
!               write(*,*) "ys:",b%r_yV(b%r_Vj1:b%r_Vj2)
!               write(*,*) "zs:",b%r_z(b%r_Vk1:b%r_Vk2)
!               call error_stop()
!             end if
! 
!             out(i,j,k) = val
!           end do
!         end do
!       end do
! 
!       call spl%destroy()
!     end subroutine
! 
!     subroutine interpolate_W_spline(in, out)
!       use bspline_oo_module
!       real(knd), intent(in), allocatable :: in(:,:,:)
!       real(knd), intent(inout), allocatable :: out(:,:,:)
!       integer :: i, j, k
!       type(bspline_3d) :: spl
!       integer :: iflag, idx, idy, idz
!       real(real64) :: val
! 
!       idx = 0; idy = 0; idz = 0
!       
!       call spl%initialize(b%r_x(b%r_Wi1:b%r_Wi2), &
!                           b%r_y(b%r_Wj1:b%r_Wj2), &
!                           b%r_zW(b%r_Wk1:b%r_Wk2), &
!                           b%r_W, &
!                           b%interp_order(1), b%interp_order(2), b%interp_order(3), &
!                           iflag)
!       if (iflag/=1) then
!         write(*,*) "Error in spline initialization, iflag:",iflag
!         stop
!       end if
!       
!       
!       do k = b%Wk1, b%Wk2
!         do j = b%Wj1, b%Wj2
!           do i = b%Wi1, b%Wi2
!             call spl%evaluate(xPr(i), yPr(j), zW(k), &
!                               idx, idy, idz, val, iflag)
!             if (iflag/=0) then
!               write(*,*) "Error in spline interpolation of W at", &
!                           xPr(i), yPr(j), zW(k)
!               write(*,*) "flag:",iflag
!               write(*,*) "xs:",b%r_x(b%r_Wi1:b%r_Wi2)
!               write(*,*) "ys:",b%r_y(b%r_Wj1:b%r_Wj2)
!               write(*,*) "zs:",b%r_zW(b%r_Wk1:b%r_Wk2)
!               call error_stop()
!             end if
! 
!             out(i,j,k) = val
!           end do
!         end do
!       end do
! 
!       call spl%destroy()
!     end subroutine

    subroutine interpolate_U_trilinear(in, out, treat_bottom)
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout), allocatable :: out(:,:,:)
      logical, intent(in), optional :: treat_bottom
      integer :: xi, yj, zk
      integer :: i, j, k
      integer :: lo
      
      !still a HACK
      !it makes the mean flow velocity to be compatible with the lowes cell of the parent
      if (.false..and.present(treat_bottom) .and. b%direction/=To .and. kim==1 .and. Btype(Bo)==BC_NOSLIP) then

        lo = b%spatial_ratio + 1

        !$omp parallel do private(i, j, k, xi, yj, zk)
        do k = b%Uk1, b%Uk2
          do j = b%Uj1, b%Uj2
            do i = b%Ui1, b%Ui2
              call U_r_index(xU(i), yPr(j), zPr(k), xi, yj, zk)

              out(i,j,k) = &
                TriLinInt((xU(i)   - b%r_xU(xi)) / b%r_dx, &
                          (yPr(j)  - b%r_y(yj) ) / b%r_dy, &
                          (zPr(k)  - b%r_z(zk) ) / b%r_dz, &
                          in(xi  , yj  , zk  ), &
                          in(xi+1, yj  , zk  ), &
                          in(xi  , yj+1, zk  ), &
                          in(xi  , yj  , zk+1), &
                          in(xi+1, yj+1, zk  ), &
                          in(xi+1, yj  , zk+1), &
                          in(xi  , yj+1, zk+1), &
                          in(xi+1, yj+1, zk+1))
            end do
          end do
        end do

        lo = b%spatial_ratio/2 + mod(b%spatial_ratio,2)

        block
          real(knd) :: ustar, h, Uout, us(b%spatial_ratio)

          do j = b%Uj1, b%Uj2,b%spatial_ratio
            do i = b%Ui1, b%Ui2, b%spatial_ratio
              call U_r_index(xU(i), yPr(j), zPr(1), xi, yj, zk)
!this interpolates to the centre of a Pr cell. 
!move to a U cell?
              Uout = TriLinInt((xU(i)   - b%r_xU(xi)) / b%r_dx, &
                          (yPr(j)  - b%r_y(yj) ) / b%r_dy, &
                          1._knd, &
                          in(xi  , yj  , zk  ), &
                          in(xi+1, yj  , zk  ), &
                          in(xi  , yj+1, zk  ), &
                          in(xi  , yj  , zk+1), &
                          in(xi+1, yj+1, zk  ), &
                          in(xi+1, yj  , zk+1), &
                          in(xi  , yj+1, zk+1), &
                          in(xi+1, yj+1, zk+1))

              h = (zPr(b%spatial_ratio)+zW(0))/2
              ustar = Uout * 0.41 /log(h/z0B)

              do k = 1, lo-1
                us(k) = ustar/0.41 * log(zPr(k)/z0B)
              end do

              do k = 1, lo-1
                associate (o => out(i:min(i+b%spatial_ratio-1,b%Ui2),j:min(j+b%spatial_ratio-1,b%Uj2),k))
                  o = o * us(k) / (sum(o)/size(o))
                end associate
              end do

            end do
          end do
          do k = b%Uk1,  0
            out(:,:,k) = - out(:,:,1-k)
          end do
        end block

      else

        !$omp parallel do private(i, j, k, xi, yj, zk)
        do k = b%Uk1, b%Uk2
          do j = b%Uj1, b%Uj2
            do i = b%Ui1, b%Ui2
              call U_r_index(xU(i), yPr(j), zPr(k), xi, yj, zk)

              out(i,j,k) = &
                TriLinInt((xU(i)   - b%r_xU(xi)) / b%r_dx, &
                          (yPr(j)  - b%r_y(yj) ) / b%r_dy, &
                          (zPr(k)  - b%r_z(zk) ) / b%r_dz, &
                          in(xi  , yj  , zk  ), &
                          in(xi+1, yj  , zk  ), &
                          in(xi  , yj+1, zk  ), &
                          in(xi  , yj  , zk+1), &
                          in(xi+1, yj+1, zk  ), &
                          in(xi+1, yj  , zk+1), &
                          in(xi  , yj+1, zk+1), &
                          in(xi+1, yj+1, zk+1))
            end do
          end do
        end do

      end if

    end subroutine


    pure subroutine U_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_xU,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_xU(lx))/b%r_dx ), 0) + lx, ubound(b%r_xU,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

    subroutine interpolate_V_trilinear(in, out)
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout), allocatable :: out(:,:,:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = b%Vk1, b%Vk2
        do j = b%Vj1, b%Vj2
          do i = b%Vi1, b%Vi2
            call V_r_index(xPr(i), yV(j), zPr(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xPr(i) - b%r_x(xi) ) / b%r_dx, &
                        (yV(j)  - b%r_yV(yj)) / b%r_dy, &
                        (zPr(k) - b%r_z(zk) ) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine

    pure subroutine V_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_yV,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_yV(ly))/b%r_dy ), 0) + ly, ubound(b%r_yV,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

    subroutine interpolate_W_trilinear(in, out)
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout), allocatable :: out(:,:,:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = b%Wk1, b%Wk2
        do j = b%Wj1, b%Wj2
          do i = b%Wi1, b%Wi2
            call W_r_index(xPr(i), yPr(j), zW(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xPr(i) - b%r_x(xi) ) / b%r_dx, &
                        (yV(j)  - b%r_y(yj) ) / b%r_dy, &
                        (zPr(k) - b%r_zW(zk)) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine

    pure subroutine W_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_zW,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_zW(lz))/b%r_dz ), 0) + lz, ubound(b%r_zW,1)-1)
    end subroutine

    subroutine interpolate_scalar_trilinear(in, out)
      real(knd), intent(in) :: in(b%r_Pri1:,b%r_Prj1:,b%r_Prk1:)
      real(knd), intent(inout) :: out(b%Pri1:,b%Prj1:,b%Prk1:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = b%Prk1, b%Prk2
        do j = b%Prj1, b%Prj2
          do i = b%Pri1, b%Pri2
            call Pr_r_index(xPr(i), yPr(j), zPr(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xPr(i) - b%r_x(xi)) / b%r_dx, &
                        (yPr(j) - b%r_y(yj)) / b%r_dy, &
                        (zPr(k) - b%r_z(zk)) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine


    pure subroutine Pr_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

  end subroutine



  pure real(knd) function TriLinInt(a, b, c, &
                                   val000, val100, val010, val001, val110, val101, val011, val111)
    real(knd), intent(in) :: a, b, c
    real(knd), intent(in) :: val000, val100, val010, val001, val110, val101, val011, val111

    TriLinInt =  (1-a) * (1-b) * (1-c) * val000 + &
                 a     * (1-b) * (1-c) * val100 + &
                 (1-a) * b     * (1-c) * val010 + &
                 (1-a) * (1-b) * c     * val001 + &
                 a     * b     * (1-c) * val110 + &
                 a     * (1-b) * c     * val101 + &
                 (1-a) * b     * c     * val011 + &
                 a     * b     * c     * val111
  end function TriLinInt



  subroutine interpolate_spectral(b)
    use fftw3
    use ArrayUtilities
    class(dom_bc_buffer_refined), intent(inout) :: b
    integer :: r
    integer :: i, j, k
    real(knd) :: slope

    slope = ( avg(b%r_U(:,:,b%r_Uk2)) - avg(b%r_U(:,:,b%r_Uk1)) ) / &
            ( b%r_z(b%r_Uk2) - b%r_z(b%r_Uk1) )

    do k = b%r_Uk1, b%r_Uk2
      do j = b%r_Uj1, b%r_Uj2
        do i = b%r_Ui1, b%r_Ui2
          b%r_U(i,j,k) = b%r_U(i,j,k) - slope * (b%r_z(k)-b%r_z(b%r_Uk1))
        end do
      end do
    end do

    call do_transforms(b%interpolation%forward_U, &
                       b%interpolation%backward_U, &
                       b%r_U, &
                       b%interpolation%trans_U, &
                       b%interpolation%fft_r_U, &
                       b%interpolation%fft_U)

    call do_transforms(b%interpolation%forward_V, &
                       b%interpolation%backward_V, &
                       b%r_V, &
                       b%interpolation%trans_V, &
                       b%interpolation%fft_r_V, &
                       b%interpolation%fft_V)

    call do_transforms(b%interpolation%forward_W, &
                       b%interpolation%backward_W, &
                       b%r_W, &
                       b%interpolation%trans_W, &
                       b%interpolation%fft_r_W, &
                       b%interpolation%fft_W)

    r = b%spatial_ratio

    b%U = b%interpolation%trans_U(2*r-1 : 2*r-1 + b%Ui2-b%Ui1, &
                                  2*r+1 : 2*r+1 + b%Uj2-b%Uj1, &
                                  2*r+1 : 2*r+1 + b%Uk2-b%Uk1)

    do k = b%Uk1, b%Uk2
      do j = b%Uj1, b%Uj2
        do i = b%Ui1, b%Ui2
          b%U(i,j,k) = b%U(i,j,k) + slope * (zPr(k)-b%r_z(b%Uk1))
        end do
      end do
    end do

    b%V = b%interpolation%trans_V(2*r+1 : 2*r+1 + b%Vi2-b%Vi1, &
                                  2*r-1 : 2*r-1 + b%Vj2-b%Vj1, &
                                  2*r+1 : 2*r+1 + b%Vk2-b%Vk1)

    b%W = b%interpolation%trans_W(2*r+1 : 2*r+1 + b%Wi2-b%Wi1, &
                                  2*r+1 : 2*r+1 + b%Wj2-b%Wj1, &
                                  2*r-1 : 2*r-1 + b%Wk2-b%Wk1)

  contains

    subroutine do_transforms(forw, back, in, out, in_cf, out_cf)
      type(c_ptr), intent(in) :: forw, back
      real(knd), dimension(:,:,:), contiguous    :: in, out
      complex(knd), dimension(:,:,:), contiguous :: in_cf, out_cf
      integer :: nxi, nyi, nzi
      integer :: nxo, nyo, nzo

#ifdef DPREC
      call fftw_execute_dft_r2c(forw, in, in_cf)
#else
      call fftwf_execute_dft_r2c(forw, in, in_cf)
#endif

      nxi = size(in,1)
      nyi = size(in,2)
      nzi = size(in,3)
      nxo = size(out,1)
      nyo = size(out,2)
      nzo = size(out,3)

      out_cf = 0

      out_cf(1:nxi/2+1,1:nyi/2,1:nzi/2)         = in_cf(1:nxi/2+1,1:nyi/2,1:nzi/2)
      out_cf(1:nxi/2+1,nyo-nyi/2+1:,1:nzi/2)     = in_cf(1:nxi/2+1,nyi-nyi/2+1:,1:nzi/2)
      out_cf(1:nxi/2+1,1:nyi/2,nzo-nzi/2+1:)     = in_cf(1:nxi/2+1,1:nyi/2,nzi-nzi/2+1:)
      out_cf(1:nxi/2+1,nyo-nyi/2+1:,nzo-nzi/2+1:) = in_cf(1:nxi/2+1,nyi-nyi/2+1:,nzi-nzi/2+1:)
      if (mod(nxi,2)==0) out_cf(nxi/2+1,:,:) = out_cf(nxi/2+1,:,:)/2
      
#ifdef DPREC
      call fftw_execute_dft_c2r(back, out_cf, out)
#else
      call fftwf_execute_dft_c2r(back, out_cf, out)
#endif
      
      out = out / ((nxi)*(nyi)*(nzi))

    end subroutine

  end subroutine






  subroutine par_update_domain_bounds_UVW(U, V, W, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    real(knd), intent(in) :: eff_time
    integer :: bi
    real(knd) :: S, t_diff
    integer :: i, j, k

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do bi = lbound(domain_bc_recv_buffers_copy,1), &
                ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(bi))
            if (b%enabled) then

              where (Utype(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)<=0) &
                U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = b%U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)                 

              where (Vtype(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2)<=0) &
                V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = b%V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2)

              where (Wtype(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2)<=0) &
                W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = b%W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2)


              if (eff_time > b%time) then
                t_diff = eff_time - b%time

                where (Utype(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)<=0) &
                  U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) + &
                                                         b%dU_dt(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) * t_diff
                  

                where (Vtype(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2)<=0) &
                  V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) + &
                                                         b%dV_dt(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) * t_diff

                where (Wtype(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2)<=0) &
                  W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) + &
                                                         b%dW_dt(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) * t_diff
              end if

            end if
          end associate
        end do
      end if

      if (allocated(domain_bc_recv_buffers)) then
        do bi = lbound(domain_bc_recv_buffers,1), &
                ubound(domain_bc_recv_buffers,1)
          associate(b => domain_bc_recv_buffers(bi))
            if (b%enabled) then

              where (Utype(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)<=0)
                U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = b%U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)
              else where
                U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = 0
              end where

              where (Vtype(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2)<=0)
                V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = b%V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2)
              else where
                V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = 0
              end where

              where (Wtype(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2)<=0)
                W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = b%W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2)
              else where
                W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = 0
              end where


              if (eff_time > b%time) then
                t_diff = eff_time - b%time

                where (Utype(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)<=0)
                  U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) + &
                                                         b%dU_dt(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) * t_diff
                else where
                  U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = 0
                end where
                  

                where (Vtype(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2)<=0)
                  V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) + &
                                                         b%dV_dt(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) * t_diff
                else where
                  V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = 0
                end where

                where (Wtype(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2)<=0)
                  W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) + &
                                                         b%dW_dt(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) * t_diff
                else where
                  W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = 0
                end where
              end if

            end if
          end associate
        end do
      end if

      if (allocated(domain_bc_recv_buffers)) then
        do bi = lbound(domain_bc_recv_buffers,1), &
                ubound(domain_bc_recv_buffers,1)
          associate(b => domain_bc_recv_buffers(bi))
            if (b%enabled.and.b%turb_generator_enabled) then
              select case (b%direction)
                case (We, Ea)
                  do k = b%bUk1, b%bUk2
                    do j = b%bUj1, b%bUj2
                      do i = b%bUi1, b%bUi2
                        if (Utype(i,j,k)<=0) U(i,j,k) = U(i,j,k) + b%U_turb(j,k)
                      end do
                    end do
                  end do
                  do k = b%bVk1, b%bVk2
                    do j = b%bVj1, b%bVj2
                      do i = b%bVi1, b%bVi2
                        if (Vtype(i,j,k)<=0) V(i,j,k) = V(i,j,k) + b%V_turb(j,k)
                      end do
                    end do
                  end do
                  do k = b%bWk1, b%bWk2
                    do j = b%bWj1, b%bWj2
                      do i = b%bWi1, b%bWi2
                        if (Wtype(i,j,k)<=0) W(i,j,k) = W(i,j,k) + b%W_turb(j,k)
                      end do
                    end do
                  end do
                case (So, No)
                  do k = b%bUk1, b%bUk2
                    do j = b%bUj1, b%bUj2
                      do i = b%bUi1, b%bUi2
                        if (Utype(i,j,k)<=0) U(i,j,k) = U(i,j,k) + b%U_turb(i,k)
                      end do
                    end do
                  end do
                  do k = b%bVk1, b%bVk2
                    do j = b%bVj1, b%bVj2
                      do i = b%bVi1, b%bVi2
                        if (Vtype(i,j,k)<=0) V(i,j,k) = V(i,j,k) + b%V_turb(i,k)
                      end do
                    end do
                  end do
                  do k = b%bWk1, b%bWk2
                    do j = b%bWj1, b%bWj2
                      do i = b%bWi1, b%bWi2
                        if (Wtype(i,j,k)<=0) W(i,j,k) = W(i,j,k) + b%W_turb(i,k)
                      end do
                    end do
                  end do
                case default
                  continue
              end select
            end if
          end associate
        end do
      end if

    end if

  end subroutine



  subroutine par_update_domain_bounds_temperature(Temperature, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: Temperature
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: i

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do i = lbound(domain_bc_recv_buffers_copy,1), &
               ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(i))

              where (Prtype(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)<=0) &
                Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                  b%Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)

              t_diff = eff_time - b%time
              if (t_diff > 0) then
                where (Prtype(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)<=0) &
                  Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                    Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) + &
                    b%dTemperature_dt(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) * t_diff
              end if

          end associate
        end do
      end if

    end if

  end subroutine

  subroutine par_update_domain_bounds_moisture(Moisture, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: Moisture
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: i

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do i = lbound(domain_bc_recv_buffers_copy,1), &
               ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(i))

              where (Prtype(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)<=0) &
                Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                  b%Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)

              t_diff = eff_time - b%time
              if (t_diff > 0) then
                where (Prtype(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)<=0) &
                  Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                    Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) + &
                    b%dMoisture_dt(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) * t_diff
              end if

          end associate
        end do
      end if

    end if

  end subroutine

  subroutine par_update_domain_bounds_scalar(Scalar, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(inout) :: Scalar
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: i, scal

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do i = lbound(domain_bc_recv_buffers_copy,1), &
               ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(i))

            t_diff = eff_time - b%time

            do scal = 1, num_of_scalars
              where (Prtype(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)<=0) &
                Scalar(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2,scal) = &
                b%Scalar(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2,scal)

              if (t_diff > 0) then
                where (Prtype(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)<=0) &
                  Scalar(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2,scal) = &
                    Scalar(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2,scal) + &
                    b%dScalar_dt(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2,scal) * t_diff
              end if
            end do

          end associate
        end do
      end if

    end if

  end subroutine

  subroutine par_update_pr_gradient(eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff

    if (parent_domain>0) then
      associate(b => domain_child_buffer)
        t_diff = eff_time - b%time

        if (b%exchange_pr_gradient_x) then
          pr_gradient_x = b%pr_gradient_x + b%pr_gradient_x_dt * t_diff
        end if
        if (b%exchange_pr_gradient_y) then
          pr_gradient_y = b%pr_gradient_y + b%pr_gradient_y_dt * t_diff
        end if
        if (b%exchange_pr_gradient_z) then
          pr_gradient_z = b%pr_gradient_z + b%pr_gradient_z_dt * t_diff
        end if
      end associate
    end if

  end subroutine

  subroutine par_update_domain_bounds(U, V, W, Temperature, Moisture, Scalar, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V ,W 
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: Temperature, Moisture
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(inout) :: Scalar
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: bi

    call par_update_pr_gradient(eff_time)

    call par_update_domain_bounds_UVW(U, V, W, eff_time)

    call par_update_domain_bounds_temperature(Temperature, eff_time)

    call par_update_domain_bounds_moisture(Moisture, eff_time)

    call par_update_domain_bounds_scalar(Scalar, eff_time)

  end subroutine


  subroutine par_domain_bound_relaxation(U, V, W, Temperature, Moisture, Scalar, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V ,W 
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: Temperature, Moisture
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(inout) :: Scalar
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: bi, width
    real(knd), allocatable :: ca(:), cb(:)

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do bi = lbound(domain_bc_recv_buffers_copy,1), &
                ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(bi))

            if (b%relaxation) then

              t_diff = eff_time - b%time

              select case(b%direction)
                case (We)
                  call relax_domain_copy_We(b)
                case (Ea)
                  call relax_domain_copy_Ea(b)
                case (So)
                  call relax_domain_copy_So(b)
                case (No)
                  call relax_domain_copy_No(b)
                case (Bo)
                case (To)
                  call relax_domain_copy_To(b)
              end select
            end if
          end associate
        end do
      end if

     if (allocated(domain_bc_recv_buffers)) then
        do bi = lbound(domain_bc_recv_buffers,1), &
                ubound(domain_bc_recv_buffers,1)
          associate(b => domain_bc_recv_buffers(bi))

            if (b%relaxation) then

              width = b%spatial_ratio * 2

              cb = cb_table(width) * b%relax_factor
              where(cb>1) cb = 1
              ca = 1 - cb

              t_diff = eff_time - b%time

              select case(b%direction)
                case (We)
                  call relax_domain_We(b)
                case (Ea)
                  call relax_domain_Ea(b)
                case (So)
                  call relax_domain_So(b)
                case (No)
                  call relax_domain_No(b)
                case (Bo)
                  call relax_domain_Bo(b)
                case (To)
                  call relax_domain_To(b)
              end select
            end if
          end associate
        end do
      end if

    end if
  contains

    subroutine relax_domain_copy_We(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(1,1:Uny,1:Unz) = ca(1) * U(1,1:Uny,1:Unz) + &
              cb(1) * (b%U(1,1:Uny,1:Unz) + t_diff * b%dU_dt(1,1:Uny,1:Unz))
      U(2,1:Uny,1:Unz) = ca(2) * U(2,1:Uny,1:Unz) + &
              cb(2) * (b%U(2,1:Uny,1:Unz) + t_diff * b%dU_dt(2,1:Uny,1:Unz))

      V(1,1:Vny,1:Vnz) = ca(1) * V(1,1:Vny,1:Vnz) + &
              cb(1) * (b%V(1,1:Vny,1:Vnz) + t_diff * b%dV_dt(1,1:Vny,1:Vnz))
      V(2,1:Vny,1:Vnz) = ca(2) * V(2,1:Vny,1:Vnz) + &
              cb(2) * (b%V(2,1:Vny,1:Vnz) + t_diff * b%dV_dt(2,1:Vny,1:Vnz))

      W(1,1:Wny,1:Wnz) = ca(1) * W(1,1:Wny,1:Wnz) + &
              cb(1) * (b%W(1,1:Wny,1:Wnz) + t_diff * b%dW_dt(1,1:Wny,1:Wnz))
      W(2,1:Wny,1:Wnz) = ca(2) * W(2,1:Wny,1:Wnz) + &
              cb(2) * (b%W(2,1:Wny,1:Wnz) + t_diff * b%dW_dt(2,1:Wny,1:Wnz))
    end subroutine

    subroutine relax_domain_copy_Ea(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(Unx,1:Uny,1:Unz) = ca(1) * U(Unx,1:Uny,1:Unz) + &
              cb(1) * (b%U(Unx,1:Uny,1:Unz) + t_diff * b%dU_dt(Unx,1:Uny,1:Unz))
      U(Unx-1,1:Uny,1:Unz) = ca(2) * U(Unx-1,1:Uny,1:Unz) + &
              cb(2) * (b%U(Unx-1,1:Uny,1:Unz) + t_diff * b%dU_dt(Unx-1,1:Uny,1:Unz))

      V(Vnx,1:Vny,1:Vnz) = ca(1) * V(Vnx,1:Vny,1:Vnz) + &
              cb(1) * (b%V(Vnx,1:Vny,1:Vnz) + t_diff * b%dV_dt(Vnx,1:Vny,1:Vnz))
      V(Vnx-1,1:Vny,1:Vnz) = ca(2) * V(Vnx-1,1:Vny,1:Vnz) + &
              cb(2) * (b%V(Vnx-1,1:Vny,1:Vnz) + t_diff * b%dV_dt(Vnx-1,1:Vny,1:Vnz))

      W(Wnx,1:Wny,1:Wnz) = ca(1) * W(Wnx,1:Wny,1:Wnz) + &
              cb(1) * (b%W(Wnx,1:Wny,1:Wnz) + t_diff * b%dW_dt(Wnx,1:Wny,1:Wnz))
      W(Wnx-1,1:Wny,1:Wnz) = ca(2) * W(Wnx-1,1:Wny,1:Wnz) + &
              cb(2) * (b%W(Wnx-1,1:Wny,1:Wnz) + t_diff * b%dW_dt(Wnx-1,1:Wny,1:Wnz))
    end subroutine

    subroutine relax_domain_copy_So(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(1:Unx,1,1:Unz) = ca(1) * U(1:Unx,1,1:Unz) + &
              cb(1) * (b%U(1:Unx,1,1:Unz) + t_diff * b%dU_dt(1:Unx,1,1:Unz))
      U(1:Unx,2,1:Unz) = ca(2) * U(1:Unx,2,1:Unz) + &
              cb(2) * (b%U(1:Unx,2,1:Unz) + t_diff * b%dU_dt(1:Unx,2,1:Unz))

      V(1:Vnx,1,1:Vnz) = ca(1) * V(1:Vnx,1,1:Vnz) + &
              cb(1) * (b%V(1:Vnx,1,1:Vnz) + t_diff * b%dV_dt(1:Vnx,1,1:Vnz))
      V(1:Vnx,2,1:Vnz) = ca(2) * V(1:Vnx,2,1:Vnz) + &
              cb(2) * (b%V(1:Vnx,2,1:Vnz) + t_diff * b%dV_dt(1:Vnx,2,1:Vnz))

      W(1:Wnx,1,1:Wnz) = ca(1) * W(1:Wnx,1,1:Wnz) + &
              cb(1) * (b%W(1:Wnx,1,1:Wnz) + t_diff * b%dW_dt(1:Wnx,1,1:Wnz))
      W(1:Wnx,2,1:Wnz) = ca(2) * W(1:Wnx,2,1:Wnz) + &
              cb(2) * (b%W(1:Wnx,2,1:Wnz) + t_diff * b%dW_dt(1:Wnx,2,1:Wnz))
    end subroutine

    subroutine relax_domain_copy_No(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(1:Unx,Uny,1:Unz) = ca(1) * U(1:Unx,Uny,1:Unz) + &
              cb(1) * (b%U(1:Unx,Uny,1:Unz) + t_diff * b%dU_dt(1:Unx,Uny,1:Unz))
      U(1:Unx,Uny-1,1:Unz) = ca(2) * U(1:Unx,Uny-1,1:Unz) + &
              cb(2) * (b%U(1:Unx,Uny-1,1:Unz) + t_diff * b%dU_dt(1:Unx,Uny-1,1:Unz))

      V(1:Vnx,Vny,1:Vnz) = ca(1) * V(1:Vnx,Vny,1:Vnz) + &
              cb(1) * (b%V(1:Vnx,Vny,1:Vnz) + t_diff * b%dV_dt(1:Vnx,Vny,1:Vnz))
      V(1:Vnx,Vny-1,1:Vnz) = ca(2) * V(1:Vnx,Vny-1,1:Vnz) + &
              cb(2) * (b%V(1:Vnx,Vny-1,1:Vnz) + t_diff * b%dV_dt(1:Vnx,Vny-1,1:Vnz))

      W(1:Wnx,Wny,1:Wnz) = ca(1) * W(1:Wnx,Wny,1:Wnz) + &
              cb(1) * (b%W(1:Wnx,Wny,1:Wnz) + t_diff * b%dW_dt(1:Wnx,Wny,1:Wnz))
      W(1:Wnx,Wny-1,1:Wnz) = ca(2) * W(1:Wnx,Wny-1,1:Wnz) + &
              cb(2) * (b%W(1:Wnx,Wny-1,1:Wnz) + t_diff * b%dW_dt(1:Wnx,Wny-1,1:Wnz))
    end subroutine

    subroutine relax_domain_copy_To(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(1:Unx,1:Uny,Unz) = ca(1) * U(1:Unx,1:Uny,Unz) + &
              cb(1) * (b%U(1:Unx,1:Uny,Unz) + t_diff * b%dU_dt(1:Unx,1:Uny,Unz))
      U(1:Unx,1:Uny,Unz-1) = ca(2) * U(1:Unx,1:Uny,Unz-1) + &
              cb(2) * (b%U(1:Unx,1:Uny,Unz-1) + t_diff * b%dU_dt(1:Unx,1:Uny,Unz-1))

      V(1:Vnx,1:Vny,Vnz) = ca(1) * V(1:Vnx,1:Vny,Vnz) + &
              cb(1) * (b%V(1:Vnx,1:Vny,Vnz) + t_diff * b%dV_dt(1:Vnx,1:Vny,Vnz))
      V(1:Vnx,1:Vny,Vnz-1) = ca(2) * V(1:Vnx,1:Vny,Vnz-1) + &
              cb(2) * (b%V(1:Vnx,1:Vny,Vnz-1) + t_diff * b%dV_dt(1:Vnx,1:Vny,Vnz-1))

      W(1:Wnx,1:Wny,Wnz) = ca(1) * W(1:Wnx,1:Wny,Wnz) + &
              cb(1) * (b%W(1:Wnx,1:Wny,Wnz) + t_diff * b%dW_dt(1:Wnx,1:Wny,Wnz))
      W(1:Wnx,1:Wny,Wnz-1) = ca(2) * W(1:Wnx,1:Wny,Wnz-1) + &
              cb(2) * (b%W(1:Wnx,1:Wny,Wnz-1) + t_diff * b%dW_dt(1:Wnx,1:Wny,Wnz-1))
    end subroutine



    subroutine relax_domain_We(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i, scal

      width = b%spatial_ratio * 2

      !$omp parallel private(i)
      !$omp do
      do i = 1, width
        where (Utype(i,1:Uny,1:Unz)<=0) &
          U(i,1:Uny,1:Unz) = ca(i) * U(i,1:Uny,1:Unz) + &
                cb(i) * (b%U(i,1:Uny,1:Unz) + t_diff * b%dU_dt(i,1:Uny,1:Unz))

        where (Vtype(i,1:Vny,1:Vnz)<=0) &
          V(i,1:Vny,1:Vnz) = ca(i) * V(i,1:Vny,1:Vnz) + &
                cb(i) * (b%V(i,1:Vny,1:Vnz) + t_diff * b%dV_dt(i,1:Vny,1:Vnz))

        where (Wtype(i,1:Wny,1:Wnz)<=0) &
          W(i,1:Wny,1:Wnz) = ca(i) * W(i,1:Wny,1:Wnz) + &
                cb(i) * (b%W(i,1:Wny,1:Wnz) + t_diff * b%dW_dt(i,1:Wny,1:Wnz))
      end do
      !$omp end do nowait
      
      !$omp do collapse(2)
      do scal = 1, num_of_scalars
        do i = 1, width
          where (Prtype(i,1:Prny,1:Prnz)<=0) &
            Scalar(i,1:Prny,1:Prnz,scal) = ca(i) * Scalar(i,1:Prny,1:Prnz,scal) + &
                  cb(i) * (b%Scalar(i,1:Prny,1:Prnz,scal) + t_diff * b%dScalar_dt(i,1:Prny,1:Prnz,scal))
        end do
      end do
      !$omp end parallel
    end subroutine

    subroutine relax_domain_Ea(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i, wi, scal

      width = b%spatial_ratio * 2

      !$omp parallel private(i, wi)
      !$omp do
      do wi = 1, width
        i = Unx - wi + 1
        where (Utype(i,1:Uny,1:Unz)<=0) &
          U(i,1:Uny,1:Unz) = ca(wi) * U(i,1:Uny,1:Unz) + &
                             cb(wi) * (b%U(i,1:Uny,1:Unz) + &
                                       t_diff * b%dU_dt(i,1:Uny,1:Unz))

        i = Vnx - wi + 1
        where (Vtype(i,1:Vny,1:Vnz)<=0) &
          V(i,1:Vny,1:Vnz) = ca(wi) * V(i,1:Vny,1:Vnz) + &
                             cb(wi) * (b%V(i,1:Vny,1:Vnz) + &
                                       t_diff * b%dV_dt(i,1:Vny,1:Vnz))

        i = Wnx - wi + 1
        where (Wtype(i,1:Wny,1:Wnz)<=0) &
          W(i,1:Wny,1:Wnz) = ca(wi) * W(i,1:Wny,1:Wnz) + &
                             cb(wi) * (b%W(i,1:Wny,1:Wnz) + &
                                       t_diff * b%dW_dt(i,1:Wny,1:Wnz))
      end do
      !$omp end do nowait
      
      !$omp do collapse(2)
      do scal = 1, num_of_scalars
        do wi = 1, width
          i = Prnx - wi + 1
          where (Prtype(i,1:Prny,1:Prnz)<=0) &
            Scalar(i,1:Prny,1:Prnz,scal) = ca(wi) * Scalar(i,1:Prny,1:Prnz,scal) + &
                  cb(wi) * (b%Scalar(i,1:Prny,1:Prnz,scal) + t_diff * b%dScalar_dt(i,1:Prny,1:Prnz,scal))
        end do
      end do
      !$omp end parallel
    end subroutine

    subroutine relax_domain_So(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i, scal

      width = b%spatial_ratio * 2

      !$omp parallel private(i)
      !$omp do
      do i = 1, width
        where (Utype(1:Unx,i,1:Unz)<=0) &
          U(1:Unx,i,1:Unz) = ca(i) * U(1:Unx,i,1:Unz) + &
                cb(i) * (b%U(1:Unx,i,1:Unz) + t_diff * b%dU_dt(1:Unx,i,1:Unz))


        where (Vtype(1:Vnx,i,1:Vnz)<=0) &
          V(1:Vnx,i,1:Vnz) = ca(i) * V(1:Vnx,i,1:Vnz) + &
                cb(i) * (b%V(1:Vnx,i,1:Vnz) + t_diff * b%dV_dt(1:Vnx,i,1:Vnz))


        where (Wtype(1:Wnx,i,1:Wnz)<=0) &
          W(1:Wnx,i,1:Wnz) = ca(i) * W(1:Wnx,i,1:Wnz) + &
                cb(i) * (b%W(1:Wnx,i,1:Wnz) + t_diff * b%dW_dt(1:Wnx,i,1:Wnz))
      end do
      !$omp end do nowait
      
      !$omp do collapse(2)
      do scal = 1, num_of_scalars
        do i = 1, width
          where (Prtype(1:Prnx,i,1:Prnz)<=0) &
            Scalar(1:Prnx,i,1:Prnz,scal) = ca(i) * Scalar(1:Prnx,i,1:Prnz,scal) + &
                  cb(i) * (b%Scalar(1:Prnx,i,1:Prnz,scal) + t_diff * b%dScalar_dt(1:Prnx,i,1:Prnz,scal))
        end do
      end do
      !$omp end parallel
    end subroutine

    subroutine relax_domain_No(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i, wi, scal

      width = b%spatial_ratio * 2

      !$omp parallel private(i, wi)
      !$omp do
      do wi = 1, width
        i = Uny - wi + 1
        where (Utype(1:Unx,i,1:Unz)<=0) &
          U(1:Unx,i,1:Unz) = ca(wi) * U(1:Unx,i,1:Unz) + &
                             cb(wi) * (b%U(1:Unx,i,1:Unz) + &
                                       t_diff * b%dU_dt(1:Unx,i,1:Unz))

        i = Vny - wi + 1
        where (Vtype(1:Vnx,i,1:Vnz)<=0) &
          V(1:Vnx,i,1:Vnz) = ca(wi) * V(1:Vnx,i,1:Vnz) + &
                             cb(wi) * (b%V(1:Vnx,i,1:Vnz) + &
                                       t_diff * b%dV_dt(1:Vnx,i,1:Vnz))

        i = Wny - wi + 1
        where (Wtype(1:Wnx,i,1:Wnz)<=0) &
          W(1:Wnx,i,1:Wnz) = ca(wi) * W(1:Wnx,i,1:Wnz) + &
                             cb(wi) * (b%W(1:Wnx,i,1:Wnz) + &
                                       t_diff * b%dW_dt(1:Wnx,i,1:Wnz))
      end do
      !$omp end do nowait
      
      !$omp do collapse(2)
      do scal = 1, num_of_scalars
        do wi = 1, width
          i = Prny - wi + 1
          where (Prtype(1:Prnx,i,1:Prnz)<=0) &
            Scalar(1:Prnx,i,1:Prnz,scal) = ca(wi) * Scalar(1:Prnx,i,1:Prnz,scal) + &
                  cb(wi) * (b%Scalar(1:Prnx,i,1:Prnz,scal) + t_diff * b%dScalar_dt(1:Prnx,i,1:Prnz,scal))
        end do
      end do
      !$omp end parallel
    end subroutine

    subroutine relax_domain_Bo(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i, scal

      width = b%spatial_ratio * 2

      !$omp parallel private(i)
      !$omp do
      do i = 1, width
        where (Utype(1:Unx,1:Uny,i)<=0) &
          U(1:Unx,1:Uny,i) = ca(i) * U(1:Unx,1:Uny,i) + &
                cb(i) * (b%U(1:Unx,1:Uny,i) + t_diff * b%dU_dt(1:Unx,1:Uny,i))


        where (Vtype(1:Vnx,1:Vny,i)<=0) &
          V(1:Vnx,1:Vny,i) = ca(i) * V(1:Vnx,1:Vny,i) + &
                cb(i) * (b%V(1:Vnx,1:Vny,i) + t_diff * b%dV_dt(1:Vnx,1:Vny,i))


        where (Wtype(1:Wnx,1:Wny,i)<=0) &
          W(1:Wnx,1:Wny,i) = ca(i) * W(1:Wnx,1:Wny,i) + &
                cb(i) * (b%W(1:Wnx,1:Wny,i) + t_diff * b%dW_dt(1:Wnx,1:Wny,i))
      end do
      !$omp end do nowait
      
      !$omp do collapse(2)
      do scal = 1, num_of_scalars
        do i = 1, width
          where (Prtype(1:Prnx,1:Prny,i)<=0) &
            Scalar(1:Prnx,1:Prny,i,scal) = ca(i) * Scalar(1:Prnx,1:Prny,i,scal) + &
                  cb(i) * (b%Scalar(1:Prnx,1:Prny,i,scal) + t_diff * b%dScalar_dt(1:Prnx,1:Prny,i,scal))
        end do
      end do
      !$omp end parallel
    end subroutine

    subroutine relax_domain_To(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i, wi, scal

      width = b%spatial_ratio * 2

      !$omp parallel private(i, wi)
      !$omp do
      do wi = 1, width
        i = Unz - wi + 1
        where (Utype(1:Unx,1:Uny,i)<=0) &
          U(1:Unx,1:Uny,i) = ca(wi) * U(1:Unx,1:Uny,i) + &
                             cb(wi) * (b%U(1:Unx,1:Uny,i) + &
                                       t_diff * b%dU_dt(1:Unx,1:Uny,i))

        i = Vnz - wi + 1
        where (Vtype(1:Vnx,1:Vny,i)<=0) &
          V(1:Vnx,1:Vny,i) = ca(wi) * V(1:Vnx,1:Vny,i) + &
                             cb(wi) * (b%V(1:Vnx,1:Vny,i) + &
                                       t_diff * b%dV_dt(1:Vnx,1:Vny,i))

        i = Wnz - wi + 1
        where (Wtype(1:Wnx,1:Wny,i)<=0) &
          W(1:Wnx,1:Wny,i) = ca(wi) * W(1:Wnx,1:Wny,i) + &
                             cb(wi) * (b%W(1:Wnx,1:Wny,i) + &
                                       t_diff * b%dW_dt(1:Wnx,1:Wny,i))
      end do
      !$omp end do nowait
      
      !$omp do collapse(2)
      do scal = 1, num_of_scalars
        do wi = 1, width
          i = Prnz - wi + 1
          where (Prtype(1:Prnx,1:Prny,i)<=0) &
            Scalar(1:Prnx,1:Prny,i,scal) = ca(wi) * Scalar(1:Prnx,1:Prny,i,scal) + &
                  cb(wi) * (b%Scalar(1:Prnx,1:Prny,i,scal) + t_diff * b%dScalar_dt(1:Prnx,1:Prny,i,scal))
        end do
      end do
      !$omp end parallel
    end subroutine

    pure function ca_table(width) result(res)
      integer, intent(in) :: width
      real(knd) :: res(width)
      integer :: i
     
      res = 1 - cb_table(width)
    end function

    pure function cb_table(width) result(res)
      integer, intent(in) :: width
      real(knd) :: res(width)
      integer :: i
     
      do i = 1, width
        res(i) = (1 - DampF((width - i + 0.5_knd)/width)) / 10
      end do
   end function

    pure function DampF(x) result(res)
      real(knd) :: res
      real(knd), intent(in) :: x

      if (x<=0) then
        res = 1
      else if (x>=1) then
        res = 0
      else
        res = (1 - 0.04_knd*x**2) * &
                ( 1 - (1 - exp(10._knd*x**2)) / (1 - exp(10._knd)) )
      end if
    end function

  end subroutine par_domain_bound_relaxation




  subroutine par_receive_initial_conditions(receive, U, V, W, Pr, Temperature, Moisture, Scalar)
    real(knd), intent(inout) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:) ,W(-2:,-2:,-2:), Pr(-1:,-1:,-1:)
    real(knd), intent(inout) :: Temperature(-2:,-2:,-2:), Moisture(-2:,-2:,-2:), Scalar(-2:,-2:,-2:,1:)
    logical, intent(in) :: receive
    integer :: err
 
    if (parent_domain>0) then
      call MPI_Send(receive, 1, MPI_LOGICAL, &
                    domain_child_buffer%remote_rank, 3001, domain_child_buffer%comm, err)

      if (receive) then
        call receive_and_interpolate(domain_child_buffer, U, V, W, Pr, Temperature, Moisture, Scalar)
      end if
    end if

  end subroutine

  subroutine receive_and_interpolate(b, U, V, W, Pr, Temperature, Moisture, Scalar)
    type(dom_child_buffer), intent(inout) :: b
    real(knd), intent(inout) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:) ,W(-2:,-2:,-2:), Pr(-1:,-1:,-1:)
    real(knd), intent(inout) :: Temperature(-2:,-2:,-2:), Moisture(-2:,-2:,-2:), Scalar(-2:,-2:,-2:,1:)
    real(knd), allocatable :: tmp(:,:,:)
    integer :: n
    integer :: err

    !to avoid unpredictable allocation of large arrays during passing to MPI_Recv
    !avoiding MPI derived types for each staggered grid at least for now
    allocate(tmp(b%r_i1:b%r_i2,b%r_j1:b%r_j2,b%r_k1:b%r_k2))

    n = size(tmp)

    call MPI_Recv(tmp, n, MPI_KND, &
                  b%remote_rank, 3002, b%comm, MPI_STATUS_IGNORE, err)
    call interpolate_U_trilinear(tmp, U)

    call MPI_Recv(tmp, n, MPI_KND, &
                  b%remote_rank, 3003, b%comm, MPI_STATUS_IGNORE, err)
    call interpolate_V_trilinear(tmp, V)

    call MPI_Recv(tmp, n, MPI_KND, &
                  b%remote_rank, 3004, b%comm, MPI_STATUS_IGNORE, err)
    call interpolate_W_trilinear(tmp, W)

!     call MPI_Recv(tmp, n, MPI_KND, &
!                   b%remote_rank, 3005, b%comm, MPI_STATUS_IGNORE, err)
!     call interpolate_trilinear(tmp, Pr, 1)

    if (enable_buoyancy) then
      call MPI_Recv(tmp, n, MPI_KND, &
                    b%remote_rank, 3006, b%comm, MPI_STATUS_IGNORE, err)
      call interpolate_trilinear(tmp, Pr, -1)
    end if

    if (enable_moisture) then
      call MPI_Recv(tmp, n, MPI_KND, &
                    b%remote_rank, 3007, b%comm, MPI_STATUS_IGNORE, err)
      call interpolate_trilinear(tmp, Pr, -1)
    end if

  contains

    subroutine interpolate_trilinear(in, out, lb)
      integer, intent(in) :: lb
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout) :: out(lb:,lb:,lb:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = 1, Prnz
        do j = 1, Prny
          do i = 1, Prnx
            call r_index(xPr(i), yPr(j), zPr(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xPr(i)  - b%r_z(xi)) / b%r_dx, &
                        (yPr(j)  - b%r_y(yj)) / b%r_dy, &
                        (zPr(k)  - b%r_z(zk)) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine

    pure subroutine r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

    subroutine interpolate_U_trilinear(in, out)
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout) :: out(-2:,-2:,-2:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            call U_r_index(xU(i), yPr(j), zPr(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xU(i)   - b%r_xU(xi)) / b%r_dx, &
                        (yPr(j)  - b%r_y(yj) ) / b%r_dy, &
                        (zPr(k)  - b%r_z(zk) ) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine

    pure subroutine U_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_xU,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_xU(lx))/b%r_dx ), 0) + lx, ubound(b%r_xU,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

    subroutine interpolate_V_trilinear(in, out)
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout) :: out(-2:,-2:,-2:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            call V_r_index(xPr(i), yV(j), zPr(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xPr(i) - b%r_x(xi) ) / b%r_dx, &
                        (yV(j)  - b%r_yV(yj)) / b%r_dy, &
                        (zPr(k) - b%r_z(zk) ) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine

    pure subroutine V_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_yV,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_yV(ly))/b%r_dy ), 0) + ly, ubound(b%r_yV,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

    subroutine interpolate_W_trilinear(in, out)
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout) :: out(-2:,-2:,-2:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = 1, Wnz
        do j = 1, Wny
          do i = 1, Wnx
            call W_r_index(xPr(i), yPr(j), zW(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xPr(i) - b%r_x(xi) ) / b%r_dx, &
                        (yV(j)  - b%r_y(yj) ) / b%r_dy, &
                        (zPr(k) - b%r_zW(zk)) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine

    pure subroutine W_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_zW,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_zW(lz))/b%r_dz ), 0) + lz, ubound(b%r_zW,1)-1)
    end subroutine


  end subroutine receive_and_interpolate
  
  
  
  subroutine par_send_initial_conditions(U, V, W, Pr, Temperature, Moisture, Scalar)
    real(knd), intent(in) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:) ,W(-2:,-2:,-2:), Pr(-1:,-1:,-1:)
    real(knd), intent(in) :: Temperature(-2:,-2:,-2:), Moisture(-2:,-2:,-2:), Scalar(-2:,-2:,-2:,:)

    real(knd), allocatable :: tmp(:,:,:)
    integer :: di, i, j, k, n
    integer :: err
    logical :: send

    if (allocated(domain_parent_buffers)) then
      do di = 1, size(domain_parent_buffers)
        if (allocated(domain_parent_buffers(di)%bs)) then
          do k = lbound(domain_parent_buffers(di)%bs,3), ubound(domain_parent_buffers(di)%bs,3)
            do j = lbound(domain_parent_buffers(di)%bs,2), ubound(domain_parent_buffers(di)%bs,2)
              do i = lbound(domain_parent_buffers(di)%bs,1), ubound(domain_parent_buffers(di)%bs,1)

                associate(b => domain_parent_buffers(di)%bs(i,j,k))

                  call MPI_Recv(send, 1, MPI_LOGICAL, &
                                b%remote_rank, 3001, b%comm, MPI_STATUS_IGNORE, err)

                  if (send) then
                    allocate(tmp(b%i1:b%i2,b%j1:b%j2,b%k1:b%k2))
                    n = size(tmp)

                    tmp = U(b%i1:b%i2,b%j1:b%j2,b%k1:b%k2)
                    call MPI_Send(tmp, n, MPI_KND, &
                                  b%remote_rank, 3002, b%comm, err)

                    tmp = V(b%i1:b%i2,b%j1:b%j2,b%k1:b%k2)
                    call MPI_Send(tmp, n, MPI_KND, &
                                  b%remote_rank, 3003, b%comm, err)

                    tmp = W(b%i1:b%i2,b%j1:b%j2,b%k1:b%k2)
                    call MPI_Send(tmp, n, MPI_KND, &
                                  b%remote_rank, 3004, b%comm, err)

!                     tmp = Pr(b%i1:b%i2,b%j1:b%j2,b%k1:b%k2)
!                     call MPI_Send(tmp, n, MPI_KND, &
!                                   b%remote_rank, 3005, b%comm, err)

                    if (enable_buoyancy) then
                      tmp = Temperature(b%i1:b%i2,b%j1:b%j2,b%k1:b%k2)
                      call MPI_Send(tmp, n, MPI_KND, &
                                    b%remote_rank, 3005, b%comm, err)
                    end if

                    if (enable_moisture) then
                      tmp = Moisture(b%i1:b%i2,b%j1:b%j2,b%k1:b%k2)
                      call MPI_Send(tmp, n, MPI_KND, &
                                    b%remote_rank, 3005, b%comm, err)
                    end if

                    deallocate(tmp)
                  end if
                end associate

              end do
            end do
          end do
        end if
      end do
    end if
  end subroutine par_send_initial_conditions
  
  
  
  
  subroutine par_domain_two_way_nesting_feedback(U, V, W, Temperature, Moisture, Scalar, &
                                                time, dt)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V ,W 
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: Temperature, Moisture
    real(knd), dimension(-2:,-2:,-2:,1:), contiguous, intent(inout) :: Scalar
    real(tim), intent(in) :: time, dt
    real(knd), allocatable :: tmp(:,:,:), tmp4(:,:,:,:)
    integer :: di, i, j, k, n, scal
    integer :: m(6), tmp_lb(3)
    integer :: err
    integer(int64) :: t1, t2, timer_rate
    
    call system_clock(count=t1, count_rate=timer_rate)
                                            
    if (parent_domain>0) then
      associate(b=>domain_child_buffer)
        if (b%is_two_way_nested .and. b%time_step==b%time_step_ratio) then
          allocate(tmp(b%r_i1:b%r_i2,b%r_j1:b%r_j2,b%r_k1:b%r_k2))

          n = size(tmp)

          tmp_lb = lbound(tmp)

          call filter(tmp, tmp_lb, U, [-2, -2, -2])
          call MPI_Send(tmp, n, MPI_KND, &
                        b%remote_rank, 4002, b%comm, err)

          call filter(tmp, tmp_lb, V, [-2, -2, -2])
          call MPI_Send(tmp, n, MPI_KND, &
                        b%remote_rank, 4003, b%comm, err)

          call filter(tmp, tmp_lb, W, [-2, -2, -2])
          call MPI_Send(tmp, n, MPI_KND, &
                        b%remote_rank, 4004, b%comm, err)

          if (enable_buoyancy) then
            call filter(tmp, tmp_lb, Temperature, [-1, -1, -1])
            call MPI_Send(tmp, n, MPI_KND, &
                          b%remote_rank, 4005, b%comm, err)
          end if

          if (enable_moisture) then
            call filter(tmp, tmp_lb, Moisture, [-1, -1, -1])
            call MPI_Send(tmp, n, MPI_KND, &
                          b%remote_rank, 4006, b%comm, err)
          end if

          if (num_of_scalars > 0) then
            allocate(tmp4(b%r_i1:b%r_i2,b%r_j1:b%r_j2,b%r_k1:b%r_k2,1:num_of_scalars))
            do scal = 1, num_of_scalars
              call filter(tmp4(:,:,:,scal), tmp_lb, Scalar(:,:,:,scal), [-1, -1, -1])
            end do
            call MPI_Send(tmp4, n*num_of_scalars, MPI_KND, &
                          b%remote_rank, 4007, b%comm, err)
            deallocate(tmp4)
          end if

          deallocate(tmp)
        end if
      end associate
    end if
    
    if (allocated(domain_parent_buffers)) then
      do di = 1, size(domain_parent_buffers)
        if (allocated(domain_parent_buffers(di)%bs)) then
          do k = lbound(domain_parent_buffers(di)%bs,3), ubound(domain_parent_buffers(di)%bs,3)
            do j = lbound(domain_parent_buffers(di)%bs,2), ubound(domain_parent_buffers(di)%bs,2)
              do i = lbound(domain_parent_buffers(di)%bs,1), ubound(domain_parent_buffers(di)%bs,1)

                associate(b => domain_parent_buffers(di)%bs(i,j,k))

                  if (b%is_two_way_nested) then
                    allocate(tmp(b%i1:b%i2,b%j1:b%j2,b%k1:b%k2))
                    n = size(tmp)

                    m = 1
                    
                    where (is_boundary_domain_boundary) m = 3
                    if (b%i1==0      .and. BType(We)<BC_DOMAIN_NESTED) m(We) = 1
                    if (b%i2==Prnx+1 .and. BType(Ea)<BC_DOMAIN_NESTED) m(Ea) = 1
                    if (b%j1==0      .and. BType(So)<BC_DOMAIN_NESTED) m(So) = 1
                    if (b%j2==Prny+1 .and. BType(No)<BC_DOMAIN_NESTED) m(No) = 1
                    if (b%k1==0      .and. BType(Bo)<BC_DOMAIN_NESTED) m(Bo) = 1
                    if (b%k2==Prnz+1 .and. BType(To)<BC_DOMAIN_NESTED) m(To) = 1

                    m(Ea) = 5
                    
                    call MPI_Recv(tmp, n, MPI_KND, &
                                  b%remote_rank, 4002, b%comm, MPI_STATUS_IGNORE, err)
                    where (Utype(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))<=0) &
                      U(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6)) = &
                        tmp(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))

                    call MPI_Recv(tmp, n, MPI_KND, &
                                  b%remote_rank, 4003, b%comm, MPI_STATUS_IGNORE, err)
                    where (Vtype(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))<=0) &
                      V(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6)) = &
                      tmp(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))

                    call MPI_Recv(tmp, n, MPI_KND, &
                                  b%remote_rank, 4004, b%comm, MPI_STATUS_IGNORE, err)
                    where (Wtype(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))<=0) &
                      W(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6)) = &
                      tmp(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))

                    if (enable_buoyancy) then
                      call MPI_Recv(tmp, n, MPI_KND, &
                                    b%remote_rank, 4005, b%comm, MPI_STATUS_IGNORE, err)
                      where (Prtype(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))<=0) &
                        Temperature(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6)) = &
                        tmp(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))
                    end if

                    if (enable_moisture) then
                      call MPI_Recv(tmp, n, MPI_KND, &
                                    b%remote_rank, 4006, b%comm, MPI_STATUS_IGNORE, err)
                      where (Prtype(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))<=0) &
                        Moisture(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6)) = &
                        tmp(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))
                    end if

                    if (num_of_scalars > 0) then
                      allocate(tmp4(b%i1:b%i2,b%j1:b%j2,b%k1:b%k2,1:num_of_scalars))
                      call MPI_Recv(tmp4, n*num_of_scalars, MPI_KND, &
                                    b%remote_rank, 4007, b%comm, MPI_STATUS_IGNORE, err)
                      do scal = 1, num_of_scalars
                        where (Prtype(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6))<=0) &
                          Scalar(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6),scal) = &
                          tmp4(b%i1+m(1):b%i2-m(2),b%j1+m(3):b%j2-m(4),b%k1+m(5):b%k2-m(6),scal)
                      end do
                      deallocate(tmp4)
                    end if

                    deallocate(tmp)
                  end if
                end associate

              end do
            end do
          end do
        end if
      end do
    end if
    
    call system_clock(count=t2)
    time_communicating_domains = time_communicating_domains + real(t2-t1, dbl)/real(timer_rate,dbl)
    
  contains
  
    subroutine filter(t, t_lb, a, a_lb)
      integer, intent(in) :: t_lb(3), a_lb(3)
      real(knd), contiguous, intent(inout) :: t(t_lb(1):,t_lb(2):,t_lb(3):)
      real(knd), contiguous, intent(in) :: a(a_lb(1):,a_lb(2):,a_lb(3):)
      integer :: i, j, k
      integer :: ii, jj, kk
      integer :: s
      
      associate(b=>domain_child_buffer)
        s = b%spatial_ratio
        
        do k = b%r_k1+1, b%r_k2-1
          kk = (k - (b%r_k1+1))*s
          do j = b%r_j1+1, b%r_j2-1 
            jj = (j - (b%r_j1+1))*s
            do i = b%r_i1+1, b%r_i2-1
              ii = (i - (b%r_i1+1))*s
              t(i,j,k) = sum(a(ii + 1 : ii + s, &
                               jj + 1 : jj + s, &
                               kk + 1 : kk + s)) &
                          / s**3
            end do
          end do
        end do
      end associate
    end subroutine
    
  end subroutine par_domain_two_way_nesting_feedback





end module domains_bc_par
