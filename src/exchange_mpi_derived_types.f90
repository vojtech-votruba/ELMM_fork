module exchange_mpi_derived_types
  use Kinds
  use Parameters, only: We, Ea, So, No, Bo, To, BC_MPI_PERIODIC
  use custom_mpi, only: MPI_ORDER_FORTRAN, MPI_KND, MPI_DATATYPE_NULL
  
  implicit none
  
  integer :: recv_mpi_types(6,7) = MPI_DATATYPE_NULL, send_mpi_types(6,7) = MPI_DATATYPE_NULL
  !6 directions
  !second index:
  ! 1 U, 2 V, 3 W
  ! 4 scalars starting from -1 (Temperature, Scalar)
  ! 5 scalars starting from -0 (Pr)
  ! 6 Q

contains

  subroutine init_mpi_derived_types(s_types, r_types, nx, ny, nz, lb, width)
    integer, intent(out) :: s_types(6), r_types(6)
    integer, intent(in) :: nx, ny, nz, lb, width
    integer :: i, ierr
    
    interface
      subroutine MPI_TYPE_CREATE_SUBARRAY(NDIMS, ARRAY_OF_SIZES, ARRAY_OF_SUBSIZES, &
                                          ARRAY_OF_STARTS, ORDER, OLDTYPE, NEWTYPE, IERROR)
        integer :: NDIMS, ARRAY_OF_SIZES(*), ARRAY_OF_SUBSIZES(*)
        integer :: ARRAY_OF_STARTS(*), ORDER, OLDTYPE, NEWTYPE, IERROR

      end subroutine
      subroutine MPI_TYPE_COMMIT(DATATYPE, IERROR)
        integer :: DATATYPE, IERROR
      end subroutine
    end interface
  
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [width, ny, nz], &
                                  off(1, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(We), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [width, ny, nz], &
                                  off(nx+1-width, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(Ea), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [nx+2*width, width, nz], &
                                  off(1-width, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(So), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [nx+2*width, width, nz], &
                                  off(1-width, ny+1-width, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(No), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [nx+2*width, ny+2*width, width], &
                                  off(1-width, 1-width, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(Bo), &
                                  ierr)
 
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [nx+2*width, ny+2*width, width], &
                                  off(1-width, 1-width, nz+1-width), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(To), &
                                  ierr)
                                  
 
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [width, ny, nz], &
                                  off(1-width, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(We), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [width, ny, nz], &
                                  off(nx+1, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(Ea), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [nx+2*width, width, nz], &
                                  off(1-width, 1-width, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(So), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [nx+2*width, width, nz], &
                                  off(1-width, ny+1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(No), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [nx+2*width, ny+2*width, width], &
                                  off(1-width, 1-width, 1-width), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(Bo), &
                                  ierr)
 
    call MPI_Type_create_subarray(3, &
                                  [nx+(1-lb)*2, ny+(1-lb)*2, nz+(1-lb)*2], &
                                  [nx+2*width, ny+2*width, width], &
                                  off(1-width, 1-width, nz+1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(To), &
                                  ierr)
 
   do i = 1, 6
     call MPI_Type_commit(s_types(i), ierr)
     call MPI_Type_commit(r_types(i), ierr)
   end do
 
  contains
    
    function off(i, j, k)
      integer :: off(3)
      integer, intent(in) :: i, j, k
      
      off(1) = i - lb
      off(2) = j - lb
      off(3) = k - lb
    end function
  
  end subroutine init_mpi_derived_types
  

  

  subroutine init_mpi_derived_types_Pr(s_types, r_types, nx, ny, nz, size_x, size_y, size_z)
    integer, intent(out) :: s_types(6), r_types(6)
    !the allocated size depends on the boundary conditions (Unx+3, Vny+3, Wnz+3)
    integer, intent(in) :: size_x, size_y, size_z
    integer, intent(in) :: nx, ny, nz
    integer :: i, ierr
    
    interface
      subroutine MPI_TYPE_CREATE_SUBARRAY(NDIMS, ARRAY_OF_SIZES, ARRAY_OF_SUBSIZES, &
                                          ARRAY_OF_STARTS, ORDER, OLDTYPE, NEWTYPE, IERROR)
        integer :: NDIMS, ARRAY_OF_SIZES(*), ARRAY_OF_SUBSIZES(*)
        integer :: ARRAY_OF_STARTS(*), ORDER, OLDTYPE, NEWTYPE, IERROR

      end subroutine
      subroutine MPI_TYPE_COMMIT(DATATYPE, IERROR)
        integer :: DATATYPE, IERROR
      end subroutine
    end interface
  
    call MPI_Type_create_subarray(3, &
                                  [size_x, size_y, size_z], &
                                  [2, ny, nz], &
                                  off(1, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(We), &
                                  ierr)                               
                                  
    call MPI_Type_create_subarray(3, &
                                  [size_x, size_y, size_z], &
                                  [1, ny, nz], &
                                  off(nx, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(Ea), &
                                  ierr)                               
                                  
    call MPI_Type_create_subarray(3, &
                                  [size_x, size_y, size_z], &
                                  [nx, 2, nz], &
                                  off(1, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(So), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [size_x, size_y, size_z], &
                                  [nx, 1, nz], &
                                  off(1, ny, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(No), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [size_x, size_y, size_z], &
                                  [nx, ny, 2], &
                                  off(1, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(Bo), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [size_x, size_y, size_z], &
                                  [nx, ny, 1], &
                                  off(1, 1, nz), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(To), &
                                  ierr)
                                  

                                  
    call MPI_Type_create_subarray(3, &
                                  [size_x, size_y, size_z], &
                                  [1, ny, nz], &
                                  off(0, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(We), &
                                  ierr)

    call MPI_Type_create_subarray(3, &
                                  [size_x, size_y, size_z], &
                                  [nx, 1, nz], &
                                  off(1, 0, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(So), &
                                  ierr)                               
                                  
    call MPI_Type_create_subarray(3, &
                                  [size_x, size_y, size_z], &
                                  [nx, ny, 1], &
                                  off(1, 1, 0), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(Bo), &
                                  ierr)                               
                                  
    if (size_x >= nx+3) &
      call MPI_Type_create_subarray(3, &
                                    [size_x, size_y, size_z], &
                                    [2, ny, nz], &
                                    off(nx+1, 1, 1), &
                                    MPI_ORDER_FORTRAN, &
                                    MPI_KND, &
                                    r_types(Ea), &
                                    ierr)
                                  
    if (size_y >= ny+3) &
      call MPI_Type_create_subarray(3, &
                                    [size_x, size_y, size_z], &
                                    [nx, 2, nz], &
                                    off(1, ny+1, 1), &
                                    MPI_ORDER_FORTRAN, &
                                    MPI_KND, &
                                    r_types(No), &
                                    ierr)
                                  
    if (size_z >= nz+3) &
      call MPI_Type_create_subarray(3, &
                                    [size_x, size_y, size_z], &
                                    [nx, ny, 2], &
                                    off(1, 1, nz+1), &
                                    MPI_ORDER_FORTRAN, &
                                    MPI_KND, &
                                    r_types(To), &
                                    ierr)
 
   do i = 1, 6
     if (s_types(i)/=MPI_DATATYPE_NULL) call MPI_Type_commit(s_types(i), ierr)
     if (r_types(i)/=MPI_DATATYPE_NULL) call MPI_Type_commit(r_types(i), ierr)
   end do
 
  contains
    
    function off(i, j, k)
      integer :: off(3)
      integer, intent(in) :: i, j, k
      
      off(1) = i - 0
      off(2) = j - 0
      off(3) = k - 0
    end function
  
  end subroutine init_mpi_derived_types_Pr
  



  subroutine init_mpi_derived_types_Q(s_types, r_types, nx, ny, nz)
    integer, intent(out) :: s_types(6), r_types(6)
    integer, intent(in) :: nx, ny, nz
    integer :: i, ierr
    
    interface
      subroutine MPI_TYPE_CREATE_SUBARRAY(NDIMS, ARRAY_OF_SIZES, ARRAY_OF_SUBSIZES, &
                                          ARRAY_OF_STARTS, ORDER, OLDTYPE, NEWTYPE, IERROR)
        integer :: NDIMS, ARRAY_OF_SIZES(*), ARRAY_OF_SUBSIZES(*)
        integer :: ARRAY_OF_STARTS(*), ORDER, OLDTYPE, NEWTYPE, IERROR

      end subroutine
      subroutine MPI_TYPE_COMMIT(DATATYPE, IERROR)
        integer :: DATATYPE, IERROR
      end subroutine
    end interface
  
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [1, ny, nz], &
                                  off(0, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(We), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [1, ny, nz], &
                                  off(nx+1, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(Ea), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [nx+2, 1, nz], &
                                  off(0, 0, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(So), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [nx+2, 1, nz], &
                                  off(0, ny+1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(No), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [nx+2, ny+2, 1], &
                                  off(0, 0, 0), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(Bo), &
                                  ierr)
 
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [nx+2, ny+2, 1], &
                                  off(0, 0, nz+1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  s_types(To), &
                                  ierr)
                                  
 
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [1, ny, nz], &
                                  off(1, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(We), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [1, ny, nz], &
                                  off(nx, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(Ea), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [nx+2, 1, nz], &
                                  off(0, 1, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(So), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [nx+2, 1, nz], &
                                  off(0, ny, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(No), &
                                  ierr)
                                  
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [nx+2, ny+2, 1], &
                                  off(0, 0, 1), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(Bo), &
                                  ierr)
 
    call MPI_Type_create_subarray(3, &
                                  [nx+2, ny+2, nz+2], &
                                  [nx+2, ny+2, 1], &
                                  off(0, 0, nz), &
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  r_types(To), &
                                  ierr)
 
   do i = 1, 6
     call MPI_Type_commit(s_types(i), ierr)
     call MPI_Type_commit(r_types(i), ierr)
   end do
 
  contains
    
    function off(i, j, k)
      integer :: off(3)
      integer, intent(in) :: i, j, k
      
      off(1) = i
      off(2) = j
      off(3) = k
    end function
  
  end subroutine init_mpi_derived_types_Q
  
  
  
end module exchange_mpi_derived_types

