module domains_bc_types
!   use mpi, only: MPI_PROC_NULL, MPI_COMM_NULL, MPI_STATUS_SIZE
  
  use iso_c_binding, only: c_ptr
  use custom_mpi
  use Parameters
  use TurbInlet
 
  implicit none
  
  type dom_buffer_base
    !MPI communicator for remote communication
    integer :: comm = MPI_COMM_NULL
    integer :: remote_rank = MPI_PROC_NULL
    integer :: remote_image(3) = 0

    !whether the child is two way nested in the parent
    logical :: is_two_way_nested = .false.

    !time at which the data is valid
    real(tim) :: time = -tiny(1.0_tim)

    !the ratio of grid resolutions between parent and child domains
    integer :: spatial_ratio = 1
    !The buffers are transferred every `time_step_ratio` time steps
    integer :: time_step_ratio = 1
    !count in a cycle depending on time_step_ratio for child domains
    integer :: time_step = 0
  end type

  type, extends(dom_buffer_base) :: dom_buffer_pr_gradient
    !whether child gets the pressure gradient from parent
    !makes sence when parent changes pressure-gradient dynamically
    logical :: exchange_pr_gradient_x = .false.
    logical :: exchange_pr_gradient_y = .false.
    logical :: exchange_pr_gradient_z = .false.

    real(knd) :: pr_gradient_x = 0
    real(knd) :: pr_gradient_y = 0
    real(knd) :: pr_gradient_z = 0
    real(knd) :: pr_gradient_x_dt = 0
    real(knd) :: pr_gradient_y_dt = 0
    real(knd) :: pr_gradient_z_dt = 0
  end type

  !exchanges necessary information with the child domain
  !one instance per child image
  type, extends(dom_buffer_pr_gradient) :: dom_parent_buffer

    !child image extents (0..n+1)
    integer :: i1, i2, j1, j2, k1, k2
  end type


  !exchanges necessary information with the parent domain
  !one instance only, there is only one parent image
  type, extends(dom_buffer_pr_gradient) :: dom_child_buffer
    !receive buffer grid, fields are interpolated from this grid to the finer grid
    integer :: r_i1, r_i2
    integer :: r_j1, r_j2
    integer :: r_k1, r_k2

    real(knd) :: r_dx, r_dy, r_dz
    real(knd) :: r_x0, r_y0, r_z0

    real(knd), allocatable :: r_xU(:), r_yV(:), r_zW(:)
    real(knd), allocatable :: r_x(:), r_y(:), r_z(:)
  end type

  type dom_parent_buffers_container
    integer :: remote_domain
    integer :: iim1 = 0, iim2 = -1, jim1 = 0, jim2 = -1, kim1 = 0, kim2 = -1
    type(dom_parent_buffer), allocatable :: bs(:,:,:)
  end type

  type(dom_child_buffer), allocatable :: domain_child_buffer

  !index is the number of child domain
  type(dom_parent_buffers_container), allocatable :: domain_parent_buffers(:)





  type, extends(dom_buffer_base) :: dom_bc_buffer_copy
    logical :: enabled = .false.

    logical :: rescale_compatibility = .false.
 
    logical :: relaxation = .true.
    real(knd) :: relax_factor = 1
    integer :: relaxation_width = 2

    integer :: remote_domain
    !This is the index of the cell boundary corresponding to the domain boundary 
    ! or to the nested domain boundary.
    integer :: position
    !The direction of the buffer from position, 1 to 6 (We to To).
    integer :: direction = 0

    !grid coordinates of the buffer (from 1 to 2)
    integer :: Ui1, Ui2, Vi1, Vi2, Wi1, Wi2, Pri1, Pri2
    integer :: Uj1, Uj2, Vj1, Vj2, Wj1, Wj2, Prj1, Prj2
    integer :: Uk1, Uk2, Vk1, Vk2, Wk1, Wk2, Prk1, Prk2

    !grid coordinates of the boundary region (from 1 to 2)
    integer :: bUi1, bUi2, bVi1, bVi2, bWi1, bWi2, bPri1, bPri2
    integer :: bUj1, bUj2, bVj1, bVj2, bWj1, bWj2, bPrj1, bPrj2
    integer :: bUk1, bUk2, bVk1, bVk2, bWk1, bWk2, bPrk1, bPrk2

    !the most simple type, just a copy, no interpolation
    !the grid resolution must be identical
    !the boundary position must always exactly coincide with the grid cell boundaries
    real(knd), allocatable, dimension(:,:,:) :: U, V, W, Pr, Temperature, Moisture
    real(knd), allocatable, dimension(:,:,:,:) :: Scalar
    real(knd), allocatable, dimension(:,:,:) :: dU_dt, dV_dt, dW_dt, dPr_dt, &
                                                dTemperature_dt, dMoisture_dt
    real(knd), allocatable, dimension(:,:,:,:) :: dScalar_dt

    !whether the interpolation from the receive buffers is necessary in this tim-step
    logical :: interpolate = .true.
  end type


  type spectral_interpolation
    !FFTW plans
    type(c_ptr) :: forward_U, forward_V, forward_W
    type(c_ptr) :: backward_U, backward_V, backward_W

    complex(knd), dimension(:,:,:), allocatable :: fft_r_U, fft_r_V, fft_r_W
    complex(knd), dimension(:,:,:), allocatable :: fft_U, fft_V, fft_W
    real(knd),    dimension(:,:,:), allocatable :: trans_U, trans_V, trans_W
  end type



  type, extends(dom_bc_buffer_copy) :: dom_bc_buffer_refined
    !order of accuracy of the spatial interpolation in the respective directions (0..spectral, 2..linear, 3..parabolic)
    integer :: interp_order(3) = 2

    !receive buffer grid, fields are interpolated from this grid to the finer grid
    integer :: r_Ui1, r_Ui2, r_Vi1, r_Vi2, r_Wi1, r_Wi2, r_Pri1, r_Pri2
    integer :: r_Uj1, r_Uj2, r_Vj1, r_Vj2, r_Wj1, r_Wj2, r_Prj1, r_Prj2
    integer :: r_Uk1, r_Uk2, r_Vk1, r_Vk2, r_Wk1, r_Wk2, r_Prk1, r_Prk2

    real(knd) :: r_dx, r_dy, r_dz
    real(knd) :: r_x0, r_y0, r_z0

    real(knd), allocatable :: r_xU(:), r_yV(:), r_zW(:)
    real(knd), allocatable :: r_x(:), r_y(:), r_z(:)

    !receive buffers
    real(knd), allocatable, dimension(:,:,:) :: r_U, r_V, r_W, r_Pr, r_Temperature, r_Moisture
    real(knd), allocatable, dimension(:,:,:,:) :: r_Scalar
    real(knd), allocatable, dimension(:,:,:) :: r_dU_dt, r_dV_dt, r_dW_dt, r_dPr_dt, &
                                                r_dTemperature_dt, r_dMoisture_dt
    real(knd), allocatable, dimension(:,:,:,:) :: r_dScalar_dt

    type(spectral_interpolation) :: interpolation   
  end type


  type, extends(dom_bc_buffer_refined) :: dom_bc_buffer_turbulence_generator
    logical :: turb_generator_enabled = .false.
    type(turbulence_generator_nesting), allocatable :: turb_generator
    real(knd), allocatable, dimension(:,:) :: U_turb, V_turb, W_turb
  contains
    procedure :: compute_sgs_tke => dom_bc_buffer_turbulence_generator_compute_sgs_tke
  end type
  
contains

  subroutine dom_bc_buffer_turbulence_generator_compute_sgs_tke(b)
    use Filters
    class(dom_bc_buffer_turbulence_generator), intent(inout) :: b
    real(knd), allocatable, dimension(:,:,:) :: Uf, Vf, Wf
    real(knd) :: uu
    integer :: i, j, k, xi, yj, zk

    select case (b%direction)
      case(We)
        allocate(Uf(b%r_Ui1+2-1:b%r_Ui1+2+1,b%r_Uj1:b%r_Uj2,b%r_Uk1:b%r_Uk2))
        allocate(Vf(b%r_Vi1+2-1:b%r_Vi1+2+1,b%r_Vj1:b%r_Vj2,b%r_Vk1:b%r_Vk2))
        allocate(Wf(b%r_Wi1+2-1:b%r_Wi1+2+1,b%r_Wj1:b%r_Wj2,b%r_Wk1:b%r_Wk2))

        call FilterTopHatSimple(Uf, b%r_U(b%r_Ui1+2-1:b%r_Ui1+2+1,b%r_Uj1:b%r_Uj2,b%r_Uk1:b%r_Uk2))
        call FilterTopHatSimple(Vf, b%r_V(b%r_Vi1+2-1:b%r_Vi1+2+1,b%r_Vj1:b%r_Vj2,b%r_Vk1:b%r_Vk2))
        call FilterTopHatSimple(Wf, b%r_W(b%r_Wi1+2-1:b%r_Wi1+2+1,b%r_Wj1:b%r_Wj2,b%r_Wk1:b%r_Wk2))

        Uf(b%r_Ui1+2,:,:) = (b%r_U(b%r_Ui1+2,:,:) - Uf(b%r_Ui1+2,:,:))**2
        Vf(b%r_Ui1+2,:,:) = (b%r_V(b%r_Vi1+2,:,:) - Vf(b%r_Vi1+2,:,:))**2
        Wf(b%r_Ui1+2,:,:) = (b%r_W(b%r_Wi1+2,:,:) - Wf(b%r_Wi1+2,:,:))**2

        do k = 1, Prnz
          do j = 1, Prny
            b%turb_generator%sgs_tke = 0

            call U_r_index(im_xmin, yPr(j), zPr(k), xi, yj, zk)
            uu = BiLinInt((yPr(j)  - b%r_y(yj) ) / b%r_dy, &
                          (zPr(k)  - b%r_z(zk) ) / b%r_dz, &
                          Uf(b%r_Ui1+2  , yj  , zk  ), &
                          Uf(b%r_Ui1+2  , yj+1, zk  ), &
                          Uf(b%r_Ui1+2  , yj  , zk+1), &
                          Uf(b%r_Ui1+2  , yj+1, zk+1))
            uu = max(uu,0._knd)
            b%turb_generator%sgs_tke = b%turb_generator%sgs_tke + uu

            

            call V_r_index(im_xmin, yV(j), zPr(k), xi, yj, zk)
            uu = BiLinInt((yV(j)   - b%r_y(yj) ) / b%r_dy, &
                          (zPr(k)  - b%r_z(zk) ) / b%r_dz, &
                          Vf(b%r_Vi1+2  , yj  , zk  ), &
                          Vf(b%r_Vi1+2  , yj+1, zk  ), &
                          Vf(b%r_Vi1+2  , yj  , zk+1), &
                          Vf(b%r_Vi1+2  , yj+1, zk+1))
            uu = max(uu,0._knd)
            b%turb_generator%sgs_tke = b%turb_generator%sgs_tke + uu

            call W_r_index(im_xmin, yPr(j), zW(k), xi, yj, zk)
            uu = BiLinInt((yPr(j)  - b%r_y(yj) ) / b%r_dy, &
                          (zW(k)   - b%r_z(zk) ) / b%r_dz, &
                          Wf(b%r_Wi1+2  , yj  , zk  ), &
                          Wf(b%r_Wi1+2  , yj+1, zk  ), &
                          Wf(b%r_Wi1+2  , yj  , zk+1), &
                          Wf(b%r_Wi1+2  , yj+1, zk+1))
            uu = max(uu,0._knd)
            b%turb_generator%sgs_tke = b%turb_generator%sgs_tke + uu

          end do
        end do

        call b%turb_generator%update_turbulence_profiles
            
    end select

  contains

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

    pure real(knd) function BiLinInt(a, b, &
                                     val00, val10, val01, val11)
      real(knd), intent(in) :: a, b
      real(knd), intent(in) :: val00, val10, val01, val11

      BiLinInt =  (1-a) * (1-b) * val00 + &
                   a     * (1-b) * val10 + &
                   (1-a) * b     * val01 + &
                   a     * b     * val11

    end function BiLinInt

  end subroutine dom_bc_buffer_turbulence_generator_compute_sgs_tke


end module domains_bc_types

