module custom_mpi
  use mpi
  use Kinds

  integer, parameter :: MPI_real32 = MPI_REAL4, MPI_real64 = MPI_REAL8

#ifdef DPREC
  integer, parameter :: MPI_knd = MPI_real64
#else
  integer, parameter :: MPI_knd = MPI_real32
#endif

  real(real32), pointer, contiguous :: MPI_IN_PLACE_real32(:)
  real(real64), pointer, contiguous :: MPI_IN_PLACE_real64(:)

  integer, parameter :: MPI_TAG_MAX = 32768

end module

module custom_par
  use Kinds
  use custom_mpi
  use stop_procedures
  
  implicit none

  private :: MPI_knd, MPI_real32, MPI_real64, MPI_IN_PLACE_real32, MPI_IN_PLACE_real64  
  
  integer, parameter :: PAR_COMM_NULL = MPI_COMM_NULL
  integer, parameter :: PAR_DATATYPE_NULL = MPI_DATATYPE_NULL
  integer, parameter :: PAR_KND = MPI_KND
  integer, parameter :: PAR_REAL32 = MPI_real32
  integer, parameter :: PAR_REAL64 = MPI_real64
  integer, protected :: PAR_TRIPLET = MPI_DATATYPE_NULL
  integer, protected :: PAR_TRIPLET32 = MPI_DATATYPE_NULL
  integer, protected :: PAR_TRIPLET64 = MPI_DATATYPE_NULL
  
  logical :: enable_multiple_domains = .false.
  
  !TODO: put current domain properties into a derived type
  
  !number of the domain
  integer :: domain_index = 1
  
  !number of the domain
  integer :: number_of_domains = 1

  integer :: parent_domain = 0

  integer :: parent_image(3) = [0, 0, 0]
  
  integer, allocatable :: child_domains(:)

  !child domains which intersect this image
  integer, allocatable :: image_child_domains(:)
  
  !whether given domain boundary is receiving and using boundary conditions from parent
  logical :: is_domain_boundary_nested(6) = .true.

  logical :: is_boundary_domain_boundary(6) = .true.

  !is this domain double nested in its parent?
  logical :: is_this_domain_two_way_nested = .false.

  !whether given domain boundary is generating additional synthetic turbulence
  logical :: has_domain_boundary_turbulence_generator(6) = .false.

  !whether given domain boundary has a relaxation zone
  logical :: has_domain_boundary_relaxation(6) = .true.
  integer, target :: domain_boundary_relaxation_width(6) = 2 !in parent domain units
  real(knd), target :: domain_boundary_relaxation_factor(6) = 1  

  !whether receive initial conditions from the parent
  logical :: receive_initial_conditions_from_parent = .true.

  !ratio of grid sizes between the outer domain and this domain
  integer :: domain_spatial_ratio = 1

  !ratio of time-step size between the outer domain and this domain
  integer :: domain_time_step_ratio = 1

  
  type domain_proc_grid
    integer, allocatable :: arr(:,:,:)
  end type

  !numbers of images in individual domains
  integer, allocatable :: domain_nims(:), &
                          domain_nxims(:), &
                          domain_nyims(:), &
                          domain_nzims(:), &
                          domain_comms(:), &
                          domain_groups(:), &
                          domain_comms_union(:,:), &
                          domain_groups_union(:,:)

  !image numbers and rank numbers (in world_comm)
  type(domain_proc_grid), allocatable :: domain_images_grid(:), &
                                         domain_ranks_grid(:)

  !is this domain double nested in its parent?
  logical, allocatable :: domain_is_domain_two_way_nested(:)

  !boundary, domain index
  logical, allocatable :: domain_is_boundary_nested(:,:)

  type domain_computational_grids
    real(knd) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(knd) :: dx, dy, dz
    
    real(knd), allocatable :: xPr(:), yPr(:), zPr(:), xU(:), yV(:), zW(:)
    
    !numbers of cells
    integer   :: nx, ny, nz
    !spatial ratio relative to the parent of this domain (if applicable)
    integer   :: spatial_ratio
    !offsets to global grid numbering in the domain offsets_x(nxims), offsets_y(nyims), offsets_z(nzims)
    integer, allocatable :: offsets_x(:), offsets_y(:), offsets_z(:)
    !numbers of cells nxs(nxims), nys(nyims), nzs(nzims)
    integer, allocatable :: nxs(:), nys(:), nzs(:)
    !extents images in the domain xmins(nxims)...
    real(knd), allocatable :: xmins(:), xmaxs(:), ymins(:), ymaxs(:), zmins(:), zmaxs(:)
    !whether the image plus side of the image is an internal or periodic boundary
    logical, allocatable :: internal_or_periodic_Ea(:), &
                            internal_or_periodic_No(:), &
                            internal_or_periodic_To(:)
  end type


  !the actual computational grids of each domain
  type(domain_computational_grids), allocatable :: domain_grids(:)

  
  !at which number starts numbering of this domain?
  integer :: first_domain_rank_in_world = 0
  integer :: first_domain_im_in_world = 1
  
  integer :: world_comm_size
  
  
  integer :: nims, npx = -1, npy = -1, npz = -1
  integer :: npxyz(3) = -1, pxyz(3)
  integer :: nxims, nyims, nzims
  integer :: myim, myrank, iim, jim, kim
  integer :: w_im, e_im, s_im, n_im, b_im, t_im
  integer :: w_rank, e_rank, s_rank, n_rank, b_rank, t_rank
  integer :: neigh_ims(6), neigh_ranks(6)

  integer :: my_world_im, my_world_rank
  
  !domain master image numbers in the domain communicator
  integer, parameter :: master_im = 1
  integer, parameter :: master_rank = 0
  
  integer, parameter :: world_comm = MPI_COMM_WORLD
  integer :: domain_comm = MPI_COMM_WORLD
  integer :: poisfft_comm = MPI_COMM_NULL
  integer :: cart_comm = MPI_COMM_NULL
  integer :: domain_masters_comm

  integer :: world_group = MPI_GROUP_NULL
  integer :: domain_masters_group
  
  !MPI communicators which include the inner or the outer domain
  integer :: inner_comm = MPI_COMM_NULL, outer_comm = MPI_COMM_NULL
  
  integer :: comm_plane_yz = MPI_COMM_NULL, comm_plane_xz = MPI_COMM_NULL, comm_plane_xy = MPI_COMM_NULL
  integer :: comm_row_x = MPI_COMM_NULL, comm_row_y = MPI_COMM_NULL, comm_row_z = MPI_COMM_NULL
  integer :: my_plane_yz_im, my_plane_xz_im, my_plane_xy_im, my_row_x_im, my_row_y_im, my_row_z_im
  
  
  integer :: cart_comm_dim = -1
  integer, allocatable :: images_grid(:,:,:), ranks_grid(:,:,:)
  
  interface par_co_reduce
    module procedure par_co_reduce_32
    module procedure par_co_reduce_64
    module procedure par_co_reduce_32_1d
    module procedure par_co_reduce_64_1d
    module procedure par_co_reduce_logical
    module procedure par_co_reduce_int
    module procedure par_co_reduce_int_1d
  end interface
  
  interface par_co_min
    module procedure par_co_min_32
    module procedure par_co_min_64
    module procedure par_co_min_int
    module procedure par_co_min_32_comm
    module procedure par_co_min_64_comm
    module procedure par_co_min_int_comm
  end interface
  
  interface par_co_max
    module procedure par_co_max_32
    module procedure par_co_max_64
  end interface
  
  interface par_co_sum
    module procedure par_co_sum_32
    module procedure par_co_sum_32_1d
    module procedure par_co_sum_32_comm
    module procedure par_co_sum_32_comm_1d
    module procedure par_co_sum_64
    module procedure par_co_sum_64_1d
    module procedure par_co_sum_64_comm
    module procedure par_co_sum_64_comm_1d
    module procedure par_co_sum_int
    module procedure par_co_sum_int_1d
    module procedure par_co_sum_int_comm
  end interface

  interface par_co_sum_plane_xy
    module procedure par_co_sum_plane_xy_32
    module procedure par_co_sum_plane_xy_64
  end interface

  interface par_co_sum_plane_yz
    module procedure par_co_sum_plane_yz_32
    module procedure par_co_sum_plane_yz_64
  end interface

  interface par_co_sum_plane_xz
    module procedure par_co_sum_plane_xz_32
    module procedure par_co_sum_plane_xz_64
  end interface

  interface par_sum_to_master_horizontal
    module procedure par_sum_to_master_horizontal_32_1d
    module procedure par_sum_to_master_horizontal_64_1d
  end interface

  interface par_co_broadcast
    module procedure par_co_broadcast_real32
    module procedure par_co_broadcast_real64
    module procedure par_co_broadcast_int
    module procedure par_co_broadcast_logical
  end interface

  interface par_broadcast_from_last_x
    module procedure par_broadcast_from_last_x_real32
    module procedure par_broadcast_from_last_x_real64
  end interface

  interface par_broadcast_from_last_y
    module procedure par_broadcast_from_last_y_real32
    module procedure par_broadcast_from_last_y_real64
  end interface

  interface par_broadcast_from_last_z
    module procedure par_broadcast_from_last_z_real32
    module procedure par_broadcast_from_last_z_real64
  end interface
  
  interface par_send
    module procedure par_send_real32
    module procedure par_send_real64
    module procedure par_send_int
    module procedure par_send_logical
  end interface

  interface par_recv
    module procedure par_recv_real32
    module procedure par_recv_real64
    module procedure par_recv_int
    module procedure par_recv_logical
  end interface

contains

  subroutine par_sync_all()
    integer :: ie
    call MPI_Barrier(domain_comm, ie)
    if (ie/=0) call error_stop("Error in MPI Barrier.")
  end subroutine

  subroutine par_sync_out(str)
    character(*), intent(in) :: str
    call par_sync_all()
    if (master) write(*,*) str
    call par_sync_all()
  end subroutine

 
  integer function par_this_image(comm) result(res)
    integer, intent(in), optional :: comm
    integer :: ie

    if (present(comm)) then
      call MPI_Comm_rank(comm, res, ie)
    else
      call MPI_Comm_rank(domain_comm, res, ie)
    end if
    res = res + 1
    if (ie/=0) call error_stop("MPI_Comm_rank ERROR")
  end function
  

  integer function par_num_images() result(res)
    integer :: ie
    call MPI_Comm_size(domain_comm, res, ie)  
    if (ie/=0) call error_stop("MPI_Comm_size ERROR")
  end function

  
  integer function par_image_index(sub) result(res)
    integer, intent(in) :: sub(3)
    integer :: ie
    
    if (cart_comm_dim==-1) then
      call MPI_Cartdim_get(cart_comm, cart_comm_dim, ie)
      if (ie/=0) call error_stop("MPI_Cartdim_get")
    end if
    call MPI_Cart_rank(cart_comm, sub(3:1:-1)-1, res, ie)  
    if (ie/=0) call error_stop("MPI_Cart_rank")
    res = res + 1
  end function
  

  subroutine par_get_image_coords()
    integer :: ie
    integer :: i, j, k
    
    call MPI_Cart_coords(cart_comm, myrank, 3, pxyz, ie)
    if (ie/=0) call error_stop("MPI_Cart_coords")
           
    pxyz = pxyz(3:1:-1)
           
    iim = pxyz(1) + 1
    jim = pxyz(2) + 1
    kim = pxyz(3) + 1
    
    write(*,*) myim, "coords:",iim,jim,kim
    
    if (iim>1) then
      w_im = par_image_index([iim-1, jim, kim])
    else
      w_im = par_image_index([nxims, jim, kim])
    end if
    if (iim<nxims) then
      e_im = par_image_index([iim+1, jim, kim])
    else
      e_im = par_image_index([1, jim, kim])
    end if
    
    if (jim>1) then
      s_im = par_image_index([iim, jim-1, kim])
    else
      s_im = par_image_index([iim, nyims, kim])
    end if
    if (jim<nyims) then
      n_im = par_image_index([iim, jim+1, kim])
    else
      n_im = par_image_index([iim, 1, kim])
    end if

    if (kim>1) then
      b_im = par_image_index([iim, jim, kim-1])
    else
      b_im = par_image_index([iim, jim, nzims])
    end if
    if (kim<nzims) then
      t_im = par_image_index([iim, jim, kim+1])
    else
      t_im = par_image_index([iim, jim, 1])
    end if
    
    w_rank = w_im - 1
    e_rank = e_im - 1
    s_rank = s_im - 1
    n_rank = n_im - 1
    b_rank = b_im - 1
    t_rank = t_im - 1
    
    neigh_ranks = [w_rank, e_rank, s_rank, n_rank, b_rank, t_rank]
    
    neigh_ims = [w_im, e_im, s_im, n_im, b_im, t_im]
    
    allocate(images_grid(1:nxims, 1:nyims, 1:nzims))
    do k = 1, nzims
      do j = 1, nyims
        do i = 1, nxims
          images_grid(i,j,k) = par_image_index([i, j, k])
        end do
      end do
    end do
    
    ranks_grid = images_grid - 1
     
  end subroutine
  
  
  
  subroutine par_init
    use Kinds
    use iso_c_binding

    integer :: ie
    integer :: required, provided

    required = MPI_THREAD_SERIALIZED
  
    call MPI_Init_thread(required, provided, ie)
    if (ie/=0) call error_stop("Error in MPI_Init")

    if (provided<required) then
      write(*,*) "------------------------------"
      write(*,*) "Error, the provided MPI threading support smaller than required!"
      write(*,*) "required:", required
      write(*,*) "provided:", provided
      write(*,*) "Trying to continue anyway, but a crash is likely and the results will be questionable."
      write(*,*) "------------------------------"

      call par_sync_all()
    end if

    call c_f_pointer(my_loc(MPI_IN_PLACE), MPI_IN_PLACE_real32, [1])
    call c_f_pointer(my_loc(MPI_IN_PLACE), MPI_IN_PLACE_real64, [1])
    
    call MPI_Errhandler_set(world_comm, MPI_ERRORS_ARE_FATAL, ie)
    if (ie/=0) call error_stop("Error in MPI_Errhandler_set")
    
    
    call MPI_Comm_size(world_comm, world_comm_size, ie)
    if (ie/=0) call error_stop("Error calling MPI_Comm_size")
    
    my_world_im = par_this_image(world_comm)
    my_world_rank = my_world_im - 1
    
    call MPI_Type_contiguous(3, MPI_KND, PAR_TRIPLET, ie)
    if (ie/=0) call error_stop("Error in MPI_Type_contiguous")
    call MPI_Type_commit(PAR_TRIPLET, ie)
    if (ie/=0) call error_stop("Error in MPI_Type_commit")

    call MPI_Type_contiguous(3, MPI_REAL32, PAR_TRIPLET32, ie)
    if (ie/=0) call error_stop("Error in MPI_Type_contiguous")
    call MPI_Type_commit(PAR_TRIPLET32, ie)
    if (ie/=0) call error_stop("Error in MPI_Type_commit")

    call MPI_Type_contiguous(3, MPI_REAL64, PAR_TRIPLET64, ie)
    if (ie/=0) call error_stop("Error in MPI_Type_contiguous")
    call MPI_Type_commit(PAR_TRIPLET64, ie)
    if (ie/=0) call error_stop("Error in MPI_Type_commit")

  contains

    type(c_ptr) function my_loc(t)
      integer, intent(in), target :: t
      
      !This is not standard conforming.
      !We assume the address here is the same as the address
      ! that MPI_IN_PLACE uses when passed as an argument.
      my_loc = c_loc(t)
    end function
      
  end subroutine par_init
  
  
  
  subroutine par_init_domains
    integer :: ie
    !domain of a (world-ordered) image
    integer, allocatable :: ims_domain(:)
    integer, allocatable :: check_n(:)
    integer :: dom, dom2, i

    
    interface
      subroutine MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, &
                 RECVTYPE, COMM, IERROR)
        integer ::  SENDBUF, RECVBUF(*)
        integer :: SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM
        integer :: IERROR
      end subroutine
    end interface
    
    if (number_of_domains < domain_index) then
      write(*,*) "Error, domain index", domain_index, " larger than the number of domains", number_of_domains
      call error_stop()
    end if
    
    if (.not. allocated(child_domains)) allocate(child_domains(0))

    if (parent_domain==0) is_domain_boundary_nested = .false.
    
    allocate(check_n(world_comm_size))
    
    call MPI_AllGather(number_of_domains, 1, MPI_INTEGER, &
                       check_n, 1, MPI_INTEGER, &
                       world_comm, ie)
                       
    if (any(check_n /= number_of_domains)) then
      write(*,*) "Error, number of domains is not defined equally for all images."
      call error_stop()
    end if

  
    allocate(domain_nims(number_of_domains))
    allocate(domain_nxims(number_of_domains))
    allocate(domain_nyims(number_of_domains))
    allocate(domain_nzims(number_of_domains))
    
    domain_nims = 0
    domain_nxims = 0
    domain_nyims = 0
    domain_nzims = 0
    
    domain_nims(domain_index) = 1


    allocate(domain_comms(number_of_domains))
    allocate(domain_groups(number_of_domains))
    allocate(domain_comms_union(number_of_domains-1, 2:number_of_domains))
    allocate(domain_groups_union(number_of_domains-1, 2:number_of_domains))

    domain_comms = MPI_COMM_NULL
    domain_groups = MPI_GROUP_NULL
    domain_comms_union = MPI_COMM_NULL
    domain_groups_union = MPI_GROUP_NULL
    
    
    call MPI_Allreduce(MPI_IN_PLACE, domain_nims, number_of_domains, &
                       MPI_INTEGER, MPI_SUM, world_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Allreduce for domain_nims.")
                       
    if (enable_multiple_domains .and. domain_nims(domain_index) /= product(npxyz)) then
      write(*,*) "Error, npxyz must be specified and equal to the number of MPI processes for each domain."
      write(*,*) domain_nims(domain_index), " /= ", npxyz(1) * npxyz(2) * npxyz(3), " for domain ", domain_index
      call error_stop()
    end if
    
    if (any(domain_nims <= 0)) then
      dom = minloc(domain_nims, 1)
      write(*,*) "Error, ", domain_nims(dom), " images in domain", dom
      call error_stop()
    end if

    call MPI_Comm_group(world_comm, world_group, ie)
    if (ie/=0) call error_stop("Error calling MPI_Comm_group.")

    domain_comm = world_comm
    allocate(ims_domain(world_comm_size))
    ims_domain = 0
    ims_domain(par_this_image(world_comm)) = domain_index

    call MPI_Allreduce(MPI_IN_PLACE, ims_domain, world_comm_size, &
                       MPI_INTEGER, MPI_SUM, world_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Allreduce for ims_domain.")

    do i = 2, world_comm_size
      if (ims_domain(i-1)>ims_domain(i)) &
        call error_stop("Error, MPI ranks must be ordered to individual domains in an increasing order.")
    end do

    first_domain_rank_in_world = sum(domain_nims(1:domain_index-1))
    first_domain_im_in_world = first_domain_rank_in_world + 1

    do dom = 1, number_of_domains
      call MPI_Group_incl(world_group, domain_nims(dom), &
                          [( sum(domain_nims(1:dom-1)) + i - 1, i = 1,  domain_nims(dom) )], &
                          domain_groups(dom), &
                          ie)
      if (ie/=0) call error_stop("Error calling MPI_Group_incl.")

      call MPI_Comm_create(world_comm, domain_groups(dom), domain_comms(dom), ie)
      if (ie/=0) call error_stop("Error calling MPI_Comm_create.")
    end do

    domain_comm = domain_comms(domain_index)  
    
    do dom = 1, number_of_domains - 1
      do dom2 = dom + 1, number_of_domains
        call MPI_Group_union(min(domain_groups(dom), domain_groups(dom2)), &
                             max(domain_groups(dom), domain_groups(dom2)), &
                             domain_groups_union(dom, dom2), ie)
        if (ie/=0) call error_stop("Error calling MPI_Group_union.")

        call MPI_Comm_create(world_comm, domain_groups_union(dom, dom2), domain_comms_union(dom, dom2), ie)
        if (ie/=0) call error_stop("Error calling MPI_Comm_create.")
      end do
    end do
    
    allocate(domain_is_domain_two_way_nested(number_of_domains))
    domain_is_domain_two_way_nested = .false.

    domain_is_domain_two_way_nested(domain_index) = is_this_domain_two_way_nested

    !instead of scatter due to the 3D topology
    call MPI_Allreduce(MPI_IN_PLACE, domain_is_domain_two_way_nested, &
                       size(domain_is_domain_two_way_nested), MPI_LOGICAL, &
                       MPI_LOR, world_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Allreduce for domain_is_domain_two_way_nested.")

    
    

    allocate(domain_is_boundary_nested(6,number_of_domains))
    domain_is_boundary_nested = .false.

    domain_is_boundary_nested(:,domain_index) = is_domain_boundary_nested
    !instead of scatter due to the 3D topology
    call MPI_Allreduce(MPI_IN_PLACE, domain_is_boundary_nested, &
                       size(domain_is_boundary_nested), MPI_LOGICAL, &
                       MPI_LOR, world_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Allreduce for domain_is_boundary_nested.")

  end subroutine par_init_domains

  


  subroutine par_init_domain_grids
    use Parameters
    integer :: dom
    integer :: ie

    domain_nxims(domain_index) = nxims
    domain_nyims(domain_index) = nyims
    domain_nzims(domain_index) = nzims

    is_boundary_domain_boundary = [iim==1, iim==nxims, jim==1, jim==nyims, kim==1, kim==nzims]

    !broadcast the process domain dimensions from one image of each domain
    do dom = 1, number_of_domains
      call MPI_Bcast(domain_nxims(dom), 1, &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_nyims(dom), 1, &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_nzims(dom), 1, &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
    end do

    allocate(domain_images_grid(number_of_domains))
    allocate(domain_ranks_grid(number_of_domains))


    do dom = 1, number_of_domains
      allocate(domain_images_grid(dom)%arr(domain_nxims(dom), &
                                           domain_nyims(dom), &
                                           domain_nzims(dom)))
      domain_images_grid(dom)%arr = 0
      allocate(domain_ranks_grid(dom)%arr(domain_nxims(dom), &
                                          domain_nyims(dom), &
                                          domain_nzims(dom)))
      domain_ranks_grid(dom)%arr = 0
    end do
    domain_images_grid(domain_index)%arr(iim, jim, kim) = my_world_im

    !instead of scatter due to the 3D topology
    call MPI_Allreduce(MPI_IN_PLACE, domain_images_grid(domain_index)%arr, &
                       size(domain_images_grid(domain_index)%arr), MPI_INTEGER, MPI_SUM, domain_comm, ie)

    !broadcast all domain grids from one image of each domain
    do dom = 1, number_of_domains
      call MPI_Bcast(domain_images_grid(dom)%arr, domain_nims(dom), &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
      domain_ranks_grid(dom)%arr = domain_images_grid(dom)%arr - 1
    end do 


  
    allocate(domain_grids(number_of_domains))

    do dom = 1, number_of_domains
      allocate(domain_grids(dom)%nxs(domain_nxims(dom)))
      domain_grids(dom)%nxs = 0

      allocate(domain_grids(dom)%nys(domain_nyims(dom)))
      domain_grids(dom)%nys = 0

      allocate(domain_grids(dom)%nzs(domain_nzims(dom)))
      domain_grids(dom)%nzs = 0

      allocate(domain_grids(dom)%offsets_x(domain_nxims(dom)))
      domain_grids(dom)%offsets_x = 0

      allocate(domain_grids(dom)%offsets_y(domain_nyims(dom)))
      domain_grids(dom)%offsets_y = 0

      allocate(domain_grids(dom)%offsets_z(domain_nzims(dom)))
      domain_grids(dom)%offsets_z = 0

      allocate(domain_grids(dom)%xmins(domain_nxims(dom)))
      domain_grids(dom)%xmins = 0
      allocate(domain_grids(dom)%xmaxs(domain_nxims(dom)))
      domain_grids(dom)%xmaxs = 0

      allocate(domain_grids(dom)%ymins(domain_nyims(dom)))
      domain_grids(dom)%ymins = 0
      allocate(domain_grids(dom)%ymaxs(domain_nyims(dom)))
      domain_grids(dom)%ymaxs = 0

      allocate(domain_grids(dom)%zmins(domain_nzims(dom)))
      domain_grids(dom)%zmins = 0
      allocate(domain_grids(dom)%zmaxs(domain_nzims(dom)))
      domain_grids(dom)%zmaxs = 0

    end do

    domain_grids(domain_index)%xmin = gxmin
    domain_grids(domain_index)%xmax = gxmax
    domain_grids(domain_index)%ymin = gymin
    domain_grids(domain_index)%ymax = gymax
    domain_grids(domain_index)%zmin = gzmin
    domain_grids(domain_index)%zmax = gzmax

    domain_grids(domain_index)%dx = dxmin
    domain_grids(domain_index)%dy = dymin
    domain_grids(domain_index)%dz = dzmin

    domain_grids(domain_index)%nx = gPrnx
    domain_grids(domain_index)%ny = gPrny
    domain_grids(domain_index)%nz = gPrnz
    
    domain_grids(domain_index)%spatial_ratio = domain_spatial_ratio

    do dom = 1, number_of_domains
      call MPI_Bcast(domain_grids(dom)%xmin, 1, &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%xmax, 1, &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%ymin, 1, &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%ymax, 1, &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%zmin, 1, &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%zmax, 1, &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)

      call MPI_Bcast(domain_grids(dom)%dx, 1, &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%dy, 1, &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%dz, 1, &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)

      call MPI_Bcast(domain_grids(dom)%nx, 1, &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%ny, 1, &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%nz, 1, &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)

      call MPI_Bcast(domain_grids(dom)%spatial_ratio, 1, &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)

    end do

    domain_grids(domain_index)%nxs(iim) = Prnx
    domain_grids(domain_index)%nys(jim) = Prny
    domain_grids(domain_index)%nzs(kim) = Prnz
    
    domain_grids(domain_index)%offsets_x(iim) = offset_to_global_x
    domain_grids(domain_index)%offsets_y(jim) = offset_to_global_y
    domain_grids(domain_index)%offsets_z(kim) = offset_to_global_z
    
    domain_grids(domain_index)%xmins(iim) = im_xmin
    domain_grids(domain_index)%xmaxs(iim) = im_xmax
    domain_grids(domain_index)%ymins(jim) = im_ymin
    domain_grids(domain_index)%ymaxs(jim) = im_ymax
    domain_grids(domain_index)%zmins(kim) = im_zmin
    domain_grids(domain_index)%zmaxs(kim) = im_zmax

    
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%nxs, &
                       size(domain_grids(domain_index)%nxs), MPI_INTEGER, MPI_SUM, comm_row_x, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%nys, &
                       size(domain_grids(domain_index)%nys), MPI_INTEGER, MPI_SUM, comm_row_y, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%nzs, &
                       size(domain_grids(domain_index)%nzs), MPI_INTEGER, MPI_SUM, comm_row_z, ie)

    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%offsets_x, &
                       size(domain_grids(domain_index)%offsets_x), MPI_INTEGER, MPI_SUM, comm_row_x, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%offsets_y, &
                       size(domain_grids(domain_index)%offsets_y), MPI_INTEGER, MPI_SUM, comm_row_y, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%offsets_z, &
                       size(domain_grids(domain_index)%offsets_z), MPI_INTEGER, MPI_SUM, comm_row_z, ie)

    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%xmins, &
                       size(domain_grids(domain_index)%xmins), MPI_KND, MPI_SUM, comm_row_x, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%xmaxs, &
                       size(domain_grids(domain_index)%xmaxs), MPI_KND, MPI_SUM, comm_row_x, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%ymins, &
                       size(domain_grids(domain_index)%ymins), MPI_KND, MPI_SUM, comm_row_y, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%ymaxs, &
                       size(domain_grids(domain_index)%ymaxs), MPI_KND, MPI_SUM, comm_row_y, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%zmins, &
                       size(domain_grids(domain_index)%zmins), MPI_KND, MPI_SUM, comm_row_z, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%zmaxs, &
                       size(domain_grids(domain_index)%zmaxs), MPI_KND, MPI_SUM, comm_row_z, ie)

    do dom = 1, number_of_domains
      call MPI_Bcast(domain_grids(dom)%nxs, size(domain_grids(dom)%nxs), &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%nys, size(domain_grids(dom)%nys), &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%nzs, size(domain_grids(dom)%nzs), &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
    end do

    do dom = 1, number_of_domains
      call MPI_Bcast(domain_grids(dom)%offsets_x, size(domain_grids(dom)%offsets_x), &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%offsets_y, size(domain_grids(dom)%offsets_y), &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%offsets_z, size(domain_grids(dom)%offsets_z), &
                     MPI_INTEGER, sum(domain_nims(1:dom-1)), world_comm, ie)
    end do

    do dom = 1, number_of_domains
      call MPI_Bcast(domain_grids(dom)%xmins, size(domain_grids(dom)%xmins), &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%xmaxs, size(domain_grids(dom)%xmaxs), &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%ymins, size(domain_grids(dom)%ymins), &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%ymaxs, size(domain_grids(dom)%ymaxs), &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%zmins, size(domain_grids(dom)%zmins), &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%zmaxs, size(domain_grids(dom)%zmaxs), &
                     MPI_KND, sum(domain_nims(1:dom-1)), world_comm, ie)
    end do

    do dom = 1, number_of_domains
      allocate(domain_grids(dom)%internal_or_periodic_Ea(domain_nxims(dom)))
      allocate(domain_grids(dom)%internal_or_periodic_No(domain_nyims(dom)))
      allocate(domain_grids(dom)%internal_or_periodic_To(domain_nzims(dom)))
      domain_grids(dom)%internal_or_periodic_Ea = .false.
      domain_grids(dom)%internal_or_periodic_No = .false.
      domain_grids(dom)%internal_or_periodic_To = .false.
    end do

    domain_grids(domain_index)%internal_or_periodic_Ea(iim) = (iim<nxims) .or. Btype(Ea)==BC_PERIODIC
    domain_grids(domain_index)%internal_or_periodic_No(jim) = (jim<nyims) .or. Btype(No)==BC_PERIODIC
    domain_grids(domain_index)%internal_or_periodic_To(kim) = (kim<nzims) .or. Btype(To)==BC_PERIODIC


    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%internal_or_periodic_Ea, &
                       size(domain_grids(domain_index)%internal_or_periodic_Ea), MPI_LOGICAL, MPI_LOR, comm_row_x, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%internal_or_periodic_No, &
                       size(domain_grids(domain_index)%internal_or_periodic_No), MPI_LOGICAL, MPI_LOR, comm_row_y, ie)
    call MPI_Allreduce(MPI_IN_PLACE, domain_grids(domain_index)%internal_or_periodic_To, &
                       size(domain_grids(domain_index)%internal_or_periodic_To), MPI_LOGICAL, MPI_LOR, comm_row_z, ie)

    do dom = 1, number_of_domains
      call MPI_Bcast(domain_grids(dom)%internal_or_periodic_Ea, size(domain_grids(dom)%internal_or_periodic_Ea), &
                     MPI_LOGICAL, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%internal_or_periodic_No, size(domain_grids(dom)%internal_or_periodic_No), &
                     MPI_LOGICAL, sum(domain_nims(1:dom-1)), world_comm, ie)
      call MPI_Bcast(domain_grids(dom)%internal_or_periodic_To, size(domain_grids(dom)%internal_or_periodic_To), &
                     MPI_LOGICAL, sum(domain_nims(1:dom-1)), world_comm, ie)
    end do
  end subroutine par_init_domain_grids



  
  subroutine par_init_grid
    use PoisFFT
    use iso_c_binding
    use Parameters

    use strings, only: itoa
    
    integer(c_intptr_t) :: ng(3), nxyz(3), off(3), nxyz2(3), nsxyz2(3)
    integer :: ie
    integer :: pos
    character(80) :: str_dir

    call par_init_domains

    nims = par_num_images()

    myim = par_this_image()

    myrank = myim - 1

    master = (myim == master_im)
    
    if (any(npxyz<1)) then
      call try_2d_decompose_process_grid
    end if

    nxims = npxyz(1)
    nyims = npxyz(2)
    nzims = npxyz(3)
    
    if (product(npxyz)/=nims) then
      if (master) then
        if (enable_multiple_domains) then
          write(*,'(*(g0))') " Error on domain ", domain_index
        else
          write(*,'(*(g0))') " Error:"
        end if
        write(*,'(*(g0))') " Could not decompose the processes to Nx x Ny x Nz process grid."
        write(*,'(*(g0))') " Last grid tried: ", npxyz(1), " x ", npxyz(2), " x ", npxyz(3), ' = ', product(npxyz)
        write(*,'(*(g0))') " Number of processes: ", nims
        write(*,'(*(g0))') " Supply the requested process grid as 'ELMM npxyz=1,Ny,Nz'."
      end if
      call error_stop("Invalid process grid.")
    end if

    !The number of images does not have to be 1 necesarilly, but it is 
    ! for the optimal preformance of the FFT in Poisson solver PoisFFT.
    if (npxyz(1)/=1) call error_stop("The number of images in x direction must be one.")
    
    !Create the 3D Cartesian communicator "cart_comm" for use in CLMM
    call MPI_Cart_create(domain_comm, 3, int(npxyz(3:1:-1)), &
                         [.false.,.false.,.false.], &
                         .true., cart_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Cart_create.")
    
    domain_comm = cart_comm
    
    myim = par_this_image()
    
    myrank = myim - 1
    !very unlikely to be changed
    master = (myrank==0)
    
    !creates a 2D! Cartesian communicator "poisfft_comm"
    !2D because of the PFFT library
    !the dimensions in the poisfft_comm are z,y
    call PoisFFT_InitMPIGrid(domain_comm, npxyz(3:2:-1), poisfft_comm, ie)
    
    call par_init_sub_comms
    
    call par_get_image_coords
    
    gPrnx = Prnx
    gPrny = Prny
    gPrnz = Prnz
    
    gPrns = [gPrnx, gPrny, gPrnz]
    
    ng = gPrns

    call PoisFFT_LocalGridSize(3,ng,cart_comm,nxyz,off,nxyz2,nsxyz2)
    if (any(nxyz/=nxyz2).or.any(off/=nsxyz2)) call error_stop("Process grid not supported by the FFT library.")
    
    call MPI_Barrier(domain_comm, ie)
    if (ie/=0) call error_stop("Error in MPI Barrier.")
    
    write(*,*) iim,jim,kim, "nxyz:", nxyz
    
    call MPI_Barrier(domain_comm, ie)
    if (ie/=0) call error_stop("Error in MPI Barrier.")
    
    if (par_co_any(any(nxyz<=0))) then
      call error_stop("Zero or negative grid size on one or more images. Try different process grid.")
    end if
    
    Prnx = int(nxyz(1))
    Prny = int(nxyz(2))
    Prnz = int(nxyz(3))
    
    offset_to_global_x = int(off(1))
    offset_to_global_y = int(off(2))
    offset_to_global_z = int(off(3))
    
    offsets_to_global = int(off)

    im_xmin = gxmin + offset_to_global_x * dxmin
    im_xmax = im_xmin + Prnx * dxmin

    im_ymin = gymin + offset_to_global_y * dymin
    im_ymax = im_ymin + Prny * dymin

    if (gridtype==GRID_VARIABLE_Z) then
      im_zmin = gzW(offset_to_global_z)
      im_zmax = gzW(offset_to_global_z + Prnz)
    else
      im_zmin = gzmin + offset_to_global_z * dzmin
      im_zmax = im_zmin + Prnz * dzmin
    end if


    call par_init_domain_grids
   
    
    output_dir = output_dir // "im-" // itoa(iim) // &
                                 "-" // itoa(jim) // &
                                 "-" // itoa(kim) // "/"

    input_dir = input_dir // "im-" // itoa(iim) // &
                                           "-" // itoa(jim) // &
                                           "-" // itoa(kim) // "/"

  contains
  
    subroutine try_2d_decompose_process_grid
      integer :: i
      npxyz(1) = 1
      npxyz(2) = nint(sqrt(real(nims)))
      npxyz(3) = nims / npxyz(2)
      if (master)    write (*,*) "Trying to decompose in",npxyz,"process grid."
      
      if (product(npxyz) /= nims) then
        do i = 1, min(4, npxyz(2) - 1)
          npxyz(2) = npxyz(2) - 1
          npxyz(3) = nims / npxyz(2)
          if (master)    write (*,*) "Trying to decompose in",npxyz,"process grid."          
          if (product(npxyz) == nims) exit
        end do
      end if
    end subroutine try_2d_decompose_process_grid
    
  end subroutine par_init_grid

  
  subroutine par_init_boundaries
    use Parameters
    
    call helper(Btype)
    call helper(TempBtype)
    call helper(MoistBtype)
    call helper(ScalBtype)
    
  contains
  
    subroutine helper(Bt)
      integer, intent(inout) :: Bt(:)
      
      if (nxims>1) then
        if (iim>1) then
          Bt(We) = BC_MPI_BOUNDARY
        end if
        if (iim<nxims) then
          Bt(Ea) = BC_MPI_BOUNDARY
        end if
        if (iim==1.and.Bt(We)==BC_PERIODIC) Bt(We) = BC_MPI_PERIODIC
        if (iim==nxims.and.Bt(Ea)==BC_PERIODIC) Bt(Ea) = BC_MPI_PERIODIC
      end if
    
      if (nyims>1) then
        if (jim>1) then
          Bt(So) = BC_MPI_BOUNDARY
        end if
        if (jim<nyims) then
          Bt(No) = BC_MPI_BOUNDARY
        end if
        if (jim==1.and.Bt(So)==BC_PERIODIC) Bt(So) = BC_MPI_PERIODIC
        if (jim==nyims.and.Bt(No)==BC_PERIODIC) Bt(No) = BC_MPI_PERIODIC
      end if
    
      if (nzims>1) then
        if (kim>1) then
          Bt(Bo) = BC_MPI_BOUNDARY
        end if
        if (kim<nzims) then
          Bt(To) = BC_MPI_BOUNDARY
        end if
        if (kim==1.and.Bt(Bo)==BC_PERIODIC) Bt(Bo) = BC_MPI_PERIODIC
        if (kim==nzims.and.Bt(To)==BC_PERIODIC) Bt(To) = BC_MPI_PERIODIC
      end if
    end subroutine
    
  end subroutine par_init_boundaries
  
  subroutine par_finalize
    integer :: ie
    call MPI_Barrier(domain_comm, ie)
    if (ie/=0) call error_stop("Error when waiting before finalizing MPI.")
    call MPI_Finalize(ie)
    if (ie/=0) call error_stop("Error finalizing MPI.")
  end subroutine
  
  
  
  subroutine par_init_sub_comms
    integer :: ie
  
    !note the order of the dimensions is reversed
    call MPI_Cart_sub(cart_comm, [.true., .true., .false.], comm_plane_yz, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
    call MPI_Cart_sub(cart_comm, [.true., .false., .true.], comm_plane_xz, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
    call MPI_Cart_sub(cart_comm, [.false., .true., .true.], comm_plane_xy, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
  
    call MPI_Cart_sub(cart_comm, [.false., .false., .true.], comm_row_x, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
    call MPI_Cart_sub(cart_comm, [.false., .true., .false.], comm_row_y, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
    call MPI_Cart_sub(cart_comm, [.true., .false., .false.], comm_row_z, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
    
    my_plane_yz_im = par_this_image(comm_plane_yz)
    my_plane_xz_im = par_this_image(comm_plane_xz)
    my_plane_xy_im = par_this_image(comm_plane_xy)

    my_row_x_im = par_this_image(comm_row_x)
    my_row_y_im = par_this_image(comm_row_y)
    my_row_z_im = par_this_image(comm_row_z)

  end subroutine
  
  
  
  
  
  
  subroutine par_recv_real32(a, src)
    real(real32), intent(out) :: a
    integer, intent(in) :: src
    integer :: ie
    
    call MPI_Recv(a, 1, MPI_real32, src-1, 111, domain_comm, MPI_STATUS_IGNORE, ie)
  end subroutine

  subroutine par_recv_real64(a, src)
    real(real64), intent(out) :: a
    integer, intent(in) :: src
    integer :: ie
    
    call MPI_Recv(a, 1, MPI_real64, src-1, 111, domain_comm, MPI_STATUS_IGNORE, ie)
  end subroutine

  subroutine par_recv_int(a, src)
    integer, intent(out) :: a
    integer, intent(in) :: src
    integer :: ie
    
    call MPI_Recv(a, 1, MPI_INTEGER, src-1, 111, domain_comm, MPI_STATUS_IGNORE, ie)
  end subroutine

  subroutine par_recv_logical(a, src)
    logical, intent(out) :: a
    integer, intent(in) :: src
    integer :: ie
    
    call MPI_Recv(a, 1, MPI_LOGICAL, src-1, 111, domain_comm, MPI_STATUS_IGNORE, ie)
  end subroutine


  subroutine par_send_real32(a, dest)
    real(real32), intent(out) :: a
    integer, intent(in) :: dest
    integer :: ie
    
    call MPI_Send(a, 1, MPI_real32, dest-1, 111, domain_comm, ie)
  end subroutine

  subroutine par_send_real64(a, dest)
    real(real64), intent(out) :: a
    integer, intent(in) :: dest
    integer :: ie
    
    call MPI_Send(a, 1, MPI_real64, dest-1, 111, domain_comm, ie)
  end subroutine

  subroutine par_send_int(a, dest)
    integer, intent(out) :: a
    integer, intent(in) :: dest
    integer :: ie
    
    call MPI_Send(a, 1, MPI_INTEGER, dest-1, 111, domain_comm, ie)
  end subroutine

  subroutine par_send_logical(a, dest)
    logical, intent(out) :: a
    integer, intent(in) :: dest
    integer :: ie
    
    call MPI_Send(a, 1, MPI_LOGICAL, dest-1, 111, domain_comm, ie)
  end subroutine


  
  
  
  
  
  function par_co_reduce_32(x,op,comm) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer, intent(in) :: op, comm
    integer :: ie
    
    interface
      subroutine MPI_Allreduce(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        real(real32) :: SENDBUF, RECVBUF
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
      end subroutine
    end interface
    
    call MPI_Allreduce(x, res, &
                       count=1, datatype=MPI_KND, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_32.")
  end function

  function par_co_reduce_64(x,op,comm) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer, intent(in) :: op, comm
    integer :: ie
    
    interface
      subroutine MPI_Allreduce(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        real(real64) :: SENDBUF, RECVBUF
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
      end subroutine
    end interface
    
    call MPI_Allreduce(x, res, &
                       count=1, datatype=MPI_KND, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_64.")
  end function

  function par_co_reduce_32_1d(x,op,comm) result(res)
    real(real32),intent(in) :: x(:)
    real(real32) :: res(size(x))
    integer, intent(in) :: op, comm
    integer :: ie
    
    interface
      subroutine MPI_Allreduce(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
        real(real32) :: SENDBUF(count), RECVBUF(count)
      end subroutine
    end interface
    
    call MPI_Allreduce(x, res, &
                       count=size(x), datatype=MPI_KND, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_32.")
  end function

  function par_co_reduce_64_1d(x,op,comm) result(res)
    real(real64),intent(in) :: x(:)
    real(real64) :: res(size(x))
    integer, intent(in) :: op, comm
    integer :: ie
    
    interface
      subroutine MPI_Allreduce(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
        real(real64) :: SENDBUF(count), RECVBUF(count)
      end subroutine
    end interface
    
    call MPI_Allreduce(x, res, &
                       count=size(x), datatype=MPI_KND, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_64.")
  end function

  function par_co_reduce_logical(x,op,comm) result(res)
    logical :: res
    logical,intent(in) :: x
    integer, intent(in) :: op, comm
    integer :: ie
    
    interface
      subroutine MPI_Allreduce(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        logical :: SENDBUF, RECVBUF
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
      end subroutine
    end interface
    
    call MPI_Allreduce(x, res, &
                       count=1, datatype=MPI_LOGICAL, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_logical.")
  end function

  function par_co_reduce_int(x,op,comm) result(res)
    integer :: res
    integer,intent(in) :: x
    integer, intent(in) :: op, comm
    integer :: ie
    
    interface
      subroutine MPI_Allreduce(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        integer :: SENDBUF, RECVBUF
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
      end subroutine
    end interface
    
    call MPI_Allreduce(x, res, &
                       count=1, datatype=MPI_INTEGER, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_logical.")
  end function

  function par_co_reduce_int_1d(x,op,comm) result(res)
    integer,intent(in) :: x(:)
    integer :: res(size(x))
    integer, intent(in) :: op, comm
    integer :: ie
    
    interface
      subroutine MPI_Allreduce(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
        integer :: SENDBUF(count), RECVBUF(count)
      end subroutine
    end interface
    
    call MPI_Allreduce(x, res, &
                       count=size(x), datatype=MPI_INTEGER, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_logical.")
  end function

  function par_co_min_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_MIN, domain_comm)
  end function

  function par_co_min_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_MIN, domain_comm)
  end function

  function par_co_min_int(x) result(res)
    integer :: res
    integer,intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_MIN, domain_comm)
  end function

  function par_co_min_32_comm(x, comm) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer, intent(in) :: comm
    
    res = par_co_reduce(x, MPI_MIN, comm)
  end function

  function par_co_min_64_comm(x, comm) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer, intent(in) :: comm
    integer :: ie
    
    res = par_co_reduce(x, MPI_MIN, comm)
  end function

  function par_co_min_int_comm(x, comm) result(res)
    integer :: res
    integer,intent(in) :: x
    integer, intent(in) :: comm
    integer :: ie
    
    res = par_co_reduce(x, MPI_MIN, comm)
  end function

  function par_co_max_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_MAX, domain_comm)
  end function

  function par_co_max_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_MAX, domain_comm)
  end function

  function par_co_sum_int(x) result(res)
    integer :: res
    integer,intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, domain_comm)
  end function

  function par_co_sum_int_1d(x) result(res)
    integer,intent(in) :: x(:)
    integer :: res(size(x))
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, domain_comm)
  end function

  function par_co_sum_int_comm(x, comm) result(res)
    integer :: res
    integer,intent(in) :: x
    integer, intent(in) :: comm
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_sum_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, domain_comm)
  end function

  function par_co_sum_64_1d(x) result(res)
    real(real64),intent(in) :: x(:)
    real(real64) :: res(size(x))
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, domain_comm)
  end function

  function par_co_sum_32_1d(x) result(res)
    real(real32),intent(in) :: x(:)
    real(real32) :: res(size(x))
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, domain_comm)
  end function

  function par_co_sum_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, domain_comm)
  end function

  function par_co_sum_32_comm(x, comm) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer, intent(in) :: comm
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_sum_64_comm(x, comm) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer, intent(in) :: comm
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_sum_32_comm_1d(x, comm) result(res)
    real(real32),intent(in) :: x(:)
    real(real32) :: res(size(x))
    integer, intent(in) :: comm
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_sum_64_comm_1d(x, comm) result(res)
    real(real64),intent(in) :: x(:)
    real(real64) :: res(size(x))
    integer, intent(in) :: comm
    integer :: ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_any(x) result(res)
    logical :: res
    logical,intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_LOR, domain_comm)
  end function


  function par_co_all(x) result(res)
    logical :: res
    logical,intent(in) :: x
    integer :: ie
    
    res = par_co_reduce(x, MPI_LAND, domain_comm)
  end function

  
  function par_co_sum_plane_xy_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x

    res =  par_co_sum_32_comm(x, comm_plane_xy)
  end function

  function par_co_sum_plane_xy_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x

    res =  par_co_sum_64_comm(x, comm_plane_xy)
  end function

  function par_co_sum_plane_yz_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x

    res =  par_co_sum_32_comm(x, comm_plane_yz)
  end function

  function par_co_sum_plane_yz_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x

    res =  par_co_sum_64_comm(x, comm_plane_yz)
  end function

  function par_co_sum_plane_xz_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x

    res =  par_co_sum_32_comm(x, comm_plane_xz)
  end function

  function par_co_sum_plane_xz_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x

    res =  par_co_sum_64_comm(x, comm_plane_xz)
  end function


  subroutine par_sum_to_master_horizontal_32_1d(x)
    real(real32), intent(inout) :: x(:)
    integer :: ie
    interface
      subroutine MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR)
        import
        INTEGER    COUNT, DATATYPE, OP, ROOT, COMM, IERROR
        real(real32)  SENDBUF(*), RECVBUF(count)
      end subroutine
    end interface
    !sums the values across horizontal plane of images to the master image (iim==1 and jim==1)
    if (iim==1.and.jim==1) then
      call MPI_Reduce(MPI_IN_PLACE_real32, x, size(x), MPI_real32, MPI_SUM, 0, comm_plane_xy, ie)
    else
      call MPI_Reduce(x, x, size(x), MPI_real32, MPI_SUM, 0, comm_plane_xy, ie)
    end if
  end subroutine

  subroutine par_sum_to_master_horizontal_64_1d(x)
    real(real64), intent(inout) :: x(:)
    integer :: ie
    interface
      subroutine MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR)
        import
        INTEGER    COUNT, DATATYPE, OP, ROOT, COMM, IERROR
        real(real64)  SENDBUF(*), RECVBUF(count)
      end subroutine
    end interface
    !sums the values across horizontal plane of images to the master image (iim==1 and jim==1)
    if (iim==1.and.jim==1) then
      call MPI_Reduce(MPI_IN_PLACE_real64, x, size(x), MPI_real64, MPI_SUM, 0, comm_plane_xy, ie)
    else
      call MPI_Reduce(x, x, size(x), MPI_real64, MPI_SUM, 0, comm_plane_xy, ie)
    end if
  end subroutine



  subroutine par_co_broadcast_real32(x, source_image)
    real(real32), intent(inout) :: x
    integer, intent(in) :: source_image
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real32) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real32, source_image-1, domain_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

  subroutine par_co_broadcast_real64(x, source_image)
    real(real64), intent(inout) :: x
    integer, intent(in) :: source_image
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real64) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real64, source_image-1, domain_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

  subroutine par_co_broadcast_int(x, source_image)
    integer, intent(inout) :: x
    integer, intent(in) :: source_image
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        integer ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_INTEGER, source_image-1, domain_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

  subroutine par_co_broadcast_logical(x, source_image)
    logical, intent(inout) :: x
    integer, intent(in) :: source_image
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        logical ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_LOGICAL, source_image-1, domain_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine



  subroutine par_broadcast_from_last_x_real32(x)
    real(real32), intent(inout) :: x
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real32) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real32, nxims-1, comm_row_x, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

  subroutine par_broadcast_from_last_x_real64(x)
    real(real64), intent(inout) :: x
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real64) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real64, nxims-1, comm_row_x, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

  subroutine par_broadcast_from_last_y_real32(x)
    real(real32), intent(inout) :: x
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real32) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real32, nyims-1, comm_row_y, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

  subroutine par_broadcast_from_last_y_real64(x)
    real(real64), intent(inout) :: x
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real64) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real64, nyims-1, comm_row_y, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

  subroutine par_broadcast_from_last_z_real32(x)
    real(real32), intent(inout) :: x
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real32) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real32, nzims-1, comm_row_z, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

  subroutine par_broadcast_from_last_z_real64(x)
    real(real64), intent(inout) :: x
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real64) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real64, nzims-1, comm_row_z, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

end module custom_par
