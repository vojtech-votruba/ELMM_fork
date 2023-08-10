#ifdef PAR
module Frames_ParallelIO
  use iso_c_binding, only: c_ptr
  use Kinds
  use Parameters
  use Pthreads
  use Endianness
  use Frames_common
  use custom_par, only: PAR_COMM_NULL, PAR_DATATYPE_NULL
  use VTKFrames
  
  implicit none
  
  private
  
  public TFrameDomain_ParallelIO, AddDomain, &
         InitFrames_ParallelIO, SaveFrames_ParallelIO, FinalizeFrames_ParallelIO

  type, extends(TFrameDomain) :: TFrameDomain_ParallelIO
    ! suffix changed in the constructor to .eaf = ELMM array frame
    character(4)  :: header_suffix = ".efh" ! efh = ELMM frame header
    
    logical :: master = .false.
    
    integer :: glob_scalar_storage_size = 0, glob_vector_storage_size = 0
    
    integer :: number_of_arrays = 0
    
    integer :: comm = PAR_COMM_NULL
    integer :: par_datatype_scalar = PAR_DATATYPE_NULL, par_datatype_vector = PAR_DATATYPE_NULL
    integer :: par_filetype_scalar = PAR_DATATYPE_NULL, par_filetype_vector = PAR_DATATYPE_NULL
    
    integer :: ngxyz(3)
    integer :: lo_im_ijk(3) = [1,1,1]
    integer :: lo_ijk(3) = [1,1,1]
  contains
    procedure :: ParallelDatatypes => TFrameDomain_ParallelIO_ParallelDatatypes
    procedure :: SaveHeader => TFrameDomain_ParallelIO_SaveHeader
    procedure :: DoSave => TFrameDomain_ParallelIO_DoSave    
    procedure :: SetRest => TFrameDomain_ParallelIO_SetRest
  end type
  
  interface TFrameDomain_ParallelIO
    module procedure TFrameDomain_ParallelIO_Init
  end interface
  
  interface AddDomain
    module procedure AddFrameDomain_ParallelIO
  end interface
  
  type(TFrameDomain_ParallelIO), allocatable, target :: FrameDomains_ParallelIO(:)
  
  character, parameter :: lf = achar(10)

  integer :: frame_unit = 3999000
  
contains
  
  
  subroutine InitFrames_ParallelIO
    integer :: i
  
    if (allocated(FrameDomains_ParallelIO)) then
      do i=1, size(FrameDomains_ParallelIO)
        call FrameDomains_ParallelIO(i)%SetRest(num_of_scalars)
      end do
      
      FrameDomains_ParallelIO = pack(FrameDomains_ParallelIO, FrameDomains_ParallelIO%active)
    end if
  
  end subroutine
  

  

  subroutine add_element_fd(a,e)
    type(TFrameDomain_ParallelIO),allocatable,intent(inout) :: a(:)
    type(TFrameDomain_ParallelIO),intent(in) :: e
    type(TFrameDomain_ParallelIO),allocatable :: tmp(:)

    if (.not.allocated(a)) then
      a = [e]
    else
      call move_alloc(a,tmp)
      allocate(a(size(tmp)+1))
      a(1:size(tmp)) = tmp
      a(size(tmp)+1) = e
    end if
  end subroutine
  
  
    
  subroutine AddFrameDomain_ParallelIO(D)
    type(TFrameDomain_ParallelIO),intent(in) :: D
    
    if (allocated(FrameDomains_ParallelIO)) then
      if (any(FrameDomains_ParallelIO%base_name==D%base_name)) &
        call error_stop("Error, duplicate frame basename '"//trim(D%base_name)//"'. &
                        &Check for duplicate parallel frame labels.")
    end if
    
    call add_element_fd(FrameDomains_ParallelIO, D)

  end subroutine
   
  
  
  subroutine SaveFrames_ParallelIO(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    real(knd),intent(in) :: time
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),contiguous,intent(in) :: Pr(-1:,-1:,-1:), &
                                       Temperature(-1:,-1:,-1:), Viscosity(-1:,-1:,-1:), &
                                       Moisture(-1:,-1:,-1:), Scalar(-1:,-1:,-1:,1:)
    integer :: i
    
    if (allocated(FrameDomains_ParallelIO)) then
      do i=1,size(FrameDomains_ParallelIO)
        call FrameDomains_ParallelIO(i)%Save(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
      end do
    end if
  end subroutine

  subroutine FinalizeFrames_ParallelIO
    integer :: i
    if (allocated(FrameDomains_ParallelIO)) then
      do i=1,size(FrameDomains_ParallelIO)
        call FrameDomains_ParallelIO(i)%Finalize
      end do
      deallocate(FrameDomains_ParallelIO)
    end if
  end subroutine  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  function TFrameDomain_ParallelIO_Init(label,dimension,direction,position, &
                                        time_params,frame_flags,bbox) result(D)
    type(TFrameDomain_ParallelIO) :: D
    character(*) :: label
    integer,intent(in) :: dimension,direction
    real(knd),intent(in) :: position
    type(TFrameTimes),intent(in) :: time_params
    type(TFrameFlags),intent(in) :: frame_flags
    type(bounding_box), intent(in), optional :: bbox

    D%TFrameDomain = TFrameDomain(label,dimension,direction,position,time_params,frame_flags,bbox)
    
    D%suffix = ".eaf"
    D%base_name = trim(shared_output_dir)//"frame-"//label

  end function
  
  

  subroutine TFrameDomain_ParallelIO_SetRest(D, num_of_scalars)
    use mpi
    use Boundaries, only: GridCoords
    use custom_par, only: iim, jim, kim, domain_comm, &
                          comm_plane_xy, comm_plane_yz, comm_plane_xz
    class(TFrameDomain_ParallelIO),intent(inout) :: D
    integer, intent(in) :: num_of_scalars
    integer :: ie, rank
    
    D%active = D%InDomain()
    
    if (D%dimension == 3) then
      if (allocated(D%bbox)) then
        if (D%active) then
          call MPI_Comm_split(domain_comm, 1, 0, D%comm, ie)
          call MPI_Comm_rank(D%comm, rank, ie)
          D%master = rank==0
        else
          call MPI_Comm_split(domain_comm, MPI_UNDEFINED, 0, D%comm, ie)
        end if
      else
        D%comm = domain_comm
        D%master = master
      end if
    else
      if (D%direction == 1) then
        D%comm = comm_plane_yz
        D%master = (jim == 1 .and. kim == 1)
      else if (D%direction == 2) then
        D%comm = comm_plane_xz
        D%master = (iim == 1 .and. kim == 1)
      else if (D%direction == 3) then
        D%comm = comm_plane_xy
        D%master = (iim == 1 .and. jim == 1)
      end if
    end if
    
    if (.not.D%active) return

    call D%TFrameDomain%SetRest(num_of_scalars)
    
    call D%ParallelDatatypes
    
    if (D%master) call D%SaveHeader
  end subroutine TFrameDomain_ParallelIO_SetRest
    
    
  
  subroutine TFrameDomain_ParallelIO_ParallelDatatypes(D)
    use custom_par
    use iso_c_binding
    
    class(TFrameDomain_ParallelIO), intent(inout), target :: D
    integer :: plane_type = PAR_DATATYPE_NULL, vec_plane_type = PAR_DATATYPE_NULL
    integer :: plane_size, vec_plane_size, scalar_size, vector_size
    integer, dimension(3) :: nxyz, off
    integer :: mini, maxi, minj, maxj, mink, maxk
    integer :: ie
    integer :: nps
    integer, dimension(:), allocatable :: is, js, ks, nxs, nys, nzs
    integer :: ngx, ngy, ngz
    integer :: lo_i, lo_j, lo_k
    
    mini = D%minPri
    maxi = D%maxPri
    minj = D%minPrj
    maxj = D%maxPrj
    mink = D%minPrk
    maxk = D%maxPrk
    
    ! To avoid the MPI integer size limitation.
    ! Not easily achievable on present hardware anyway.
    plane_size = (maxi - mini + 1) * (maxj - minj + 1)
    call MPI_Type_contiguous(plane_size, PAR_REAL32, plane_type, ie)
    call MPI_Type_commit(plane_type, ie)
    
    vec_plane_size = (maxi - mini + 1) * (maxj - minj + 1)
    call MPI_Type_contiguous(vec_plane_size, PAR_TRIPLET32, vec_plane_type, ie)
    call MPI_Type_commit(vec_plane_type, ie)
    
    scalar_size = maxk - mink + 1
    call MPI_Type_contiguous(scalar_size, plane_type, D%par_datatype_scalar, ie)
    call MPI_Type_commit(D%par_datatype_scalar, ie)
    
    vector_size = maxk - mink + 1
    call MPI_Type_contiguous(vector_size, vec_plane_type, D%par_datatype_vector, ie)
    call MPI_Type_commit(D%par_datatype_vector, ie)
    
    if (allocated(D%bbox)) then
      nxyz = [maxi-mini+1, maxj-minj+1, maxk-mink+1]
    
    
      call MPI_Comm_size(D%comm, nps, ie)
      allocate(is(nps))
      call MPI_AllGather(iim, 1, MPI_INTEGER, is, 1 ,MPI_INTEGER, D%comm, ie)
      lo_i = minval(is)
!       hi_i = maxval(is)
      deallocate(is)
      
      allocate(js(nps))
      call MPI_AllGather(jim, 1, MPI_INTEGER, js, 1 ,MPI_INTEGER, D%comm, ie)
      lo_j = minval(js)
!       hi_j = maxval(js)
      deallocate(js)
      
      allocate(ks(nps))
      call MPI_AllGather(kim, 1, MPI_INTEGER, ks, 1 ,MPI_INTEGER, D%comm, ie)
      lo_k = minval(ks)
!       hi_k = maxval(ks)
      deallocate(ks)
      
      D%lo_im_ijk = [lo_i, lo_j, lo_k]
      
      allocate(nxs(nps))
      if (jim==lo_j.and.kim==lo_k) then
        call MPI_AllGather(nxyz(1), 1, MPI_INTEGER, nxs, 1 ,MPI_INTEGER, D%comm, ie)
      else
        call MPI_AllGather(0, 1, MPI_INTEGER, nxs, 1 ,MPI_INTEGER, D%comm, ie)
      end if
      nxs = pack(nxs, nxs>0)
      D%ngxyz(1) = sum(nxs)
      off(1) = sum(nxs(1:iim-lo_i))
      D%lo_ijk(1) = mini + offsets_to_global(1)
      call MPI_AllReduce(MPI_IN_PLACE, D%lo_ijk(1), 1, MPI_INTEGER, MPI_MIN, D%comm, ie)
      deallocate(nxs)
      
      allocate(nys(nps))
      if (iim==lo_i.and.kim==lo_k) then
        call MPI_AllGather(nxyz(2), 1, MPI_INTEGER, nys, 1 ,MPI_INTEGER, D%comm, ie)
      else
        call MPI_AllGather(0, 1, MPI_INTEGER, nys, 1 ,MPI_INTEGER, D%comm, ie)
      end if
      nys = pack(nys, nys>0)
      D%ngxyz(2) = sum(nys)
      off(2) = sum(nys(1:jim-lo_j))
      D%lo_ijk(2) = minj + offsets_to_global(2)
      call MPI_AllReduce(MPI_IN_PLACE, D%lo_ijk(2), 1, MPI_INTEGER, MPI_MIN, D%comm, ie)
      deallocate(nys)
      
      allocate(nzs(nps))
      if (jim==lo_j.and.iim==lo_i) then
        call MPI_AllGather(nxyz(3), 1, MPI_INTEGER, nzs, 1 ,MPI_INTEGER, D%comm, ie)
      else
        call MPI_AllGather(0, 1, MPI_INTEGER, nzs, 1 ,MPI_INTEGER, D%comm, ie)
      end if
      nzs = pack(nzs, nzs>0)
      D%ngxyz(3) = sum(nzs)
      off(3) = sum(nzs(1:kim-lo_k))
      D%lo_ijk(3) = mink + offsets_to_global(3)
      call MPI_AllReduce(MPI_IN_PLACE, D%lo_ijk(3), 1, MPI_INTEGER, MPI_MIN, D%comm, ie)
      deallocate(nzs)

      
    else
      D%ngxyz = gPrns
      nxyz = [Prnx, Prny, Prnz]      
      off = offsets_to_global
      D%lo_ijk = [1,1,1]
      if (D%dimension==2) then
        if (D%direction==1) then
          D%ngxyz(1) = 1
          nxyz(1) = 1
          off(1) = 0
          D%lo_ijk(1) = mini + offsets_to_global(1)
        else if (D%direction==2) then
          D%ngxyz(2) = 1
          nxyz(2) = 1
          off(2) = 0
          D%lo_ijk(2) = minj + offsets_to_global(2)
        else if (D%direction==3) then
          D%ngxyz(3) = 1
          nxyz(3) = 1
          off(3) = 0
          D%lo_ijk(3) = mink + offsets_to_global(3)
        end if
      end if
    end if
    
    call MPI_Type_create_subarray(3, D%ngxyz, nxyz, off, &
        MPI_ORDER_FORTRAN, PAR_REAL32, D%par_filetype_scalar, ie)
    call MPI_Type_commit(D%par_filetype_scalar, ie)

    call MPI_Type_create_subarray(3, D%ngxyz, nxyz, off, &
        MPI_ORDER_FORTRAN, PAR_TRIPLET32, D%par_filetype_vector, ie)
    call MPI_Type_commit(D%par_filetype_vector, ie)

    D%glob_scalar_storage_size = product(D%ngxyz) * &
                                 storage_size(1_real32) / CHARACTER_STORAGE_SIZE
    D%glob_vector_storage_size = 3 * D%glob_scalar_storage_size
    
  end subroutine TFrameDomain_ParallelIO_ParallelDatatypes
  
  
  
  subroutine TFrameDomain_ParallelIO_SaveHeader(D)
    ! To be executed only by D%master.
    ! To be used to decode individual frames and create vtk files.
    use custom_par
    use iso_c_binding
    
    class(TFrameDomain_ParallelIO), intent(inout), target :: D
    
    integer(int32) :: nxs(3)
    real(knd), allocatable :: xs(:), ys(:), zs(:)
    
    integer(int32), parameter :: scalar_flag = 1001, vector_flag = 1002
    
    character(2512) :: file_name
    character(8) ::  scalname
    integer :: i, j, k, sc

    scalname = "scalar00"
    file_name = ""

    file_name = trim(D%base_name)//trim(D%header_suffix)
    
    open(D%unit,file = file_name, &
      access='stream',status='replace',form="unformatted",action="write")
      
    write(D%unit) 1_int32                                         ! endianness
    write(D%unit) int(storage_size(1_real32) / &
                        CHARACTER_STORAGE_SIZE, &
                      int32)                                      ! size of real used

    xs = [(gxmin + dxmin/2 + (i-1)*dxmin, i = D%lo_ijk(1), D%ngxyz(1)+D%lo_ijk(1)-1)]
    ys = [(gymin + dymin/2 + (j-1)*dymin, j = D%lo_ijk(2), D%ngxyz(2)+D%lo_ijk(2)-1)]
    zs = [(gzPr(k), k = D%lo_ijk(3), D%ngxyz(3)+D%lo_ijk(3)-1)]
    
    if (D%dimension==2) then
      if (D%direction==1) then
        xs = [xPr(D%minPri)]        
      else if (D%direction==2) then
        ys = [yPr(D%minPrj)]
      else
        zs = [zPr(D%minPrk)]
      end if
    end if
    
    D%number_of_arrays = 0
    
    write(D%unit) int(D%ngxyz, int32)
    write(D%unit) real(xs,real32)
    write(D%unit) real(ys,real32)
    write(D%unit) real(zs,real32)

    if (D%flags%Pr==1) then
      write(D%unit) scalar_flag
      write(D%unit) "p", lf
      D%number_of_arrays = D%number_of_arrays + 1
    end if

    if (D%flags%lambda2==1) then
      write(D%unit) scalar_flag
      write(D%unit) "lambda2", lf
      D%number_of_arrays = D%number_of_arrays + 1
    end if

    if (D%flags%scalars==1) then
      scalname(1:6)="scalar"
      do sc = 1,num_of_scalars
        write(scalname(7:8),"(I2.2)") sc
        write(D%unit) scalar_flag
        write(D%unit) scalname , lf
        D%number_of_arrays = D%number_of_arrays + 1
      end do
    elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
      write(D%unit) scalar_flag
      write(D%unit) "scalar", lf
      D%number_of_arrays = D%number_of_arrays + 1
    end if

    if (enable_buoyancy.and.D%flags%temperature==1) then
      write(D%unit) scalar_flag
      write(D%unit) "temperature", lf
      D%number_of_arrays = D%number_of_arrays + 1
    end if

    if (enable_moisture.and.D%flags%moisture==1) then
      write(D%unit) scalar_flag
      write(D%unit) "moisture", lf
      D%number_of_arrays = D%number_of_arrays + 1
    end if

    if (enable_liquid.and.D%flags%liquid_water==1) then
      write(D%unit) scalar_flag
      write(D%unit) "liquid_water", lf
      D%number_of_arrays = D%number_of_arrays + 1
    end if

    if (enable_buoyancy.and.D%flags%temperature_flux==1) then
      write(D%unit) scalar_flag
      write(D%unit) "temperature_flux"
      D%number_of_arrays = D%number_of_arrays + 1
    end if

    if (enable_buoyancy.and.D%flags%moisture_flux==1) then
      write(D%unit) scalar_flag
      write(D%unit) "moisture_flux", lf
      D%number_of_arrays = D%number_of_arrays + 1
    end if

    if (D%flags%scalar_flux==1) then
      if (D%flags%scalars==1) then
        write(D%unit) scalar_flag
        scalname(1:6)="scalar_flux"
        do sc = 1,num_of_scalars
          write(scalname(7:8),"(I2.2)") sc
          write(D%unit) scalname, lf
          D%number_of_arrays = D%number_of_arrays + 1
        end do
      elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
        write(D%unit) scalar_flag
        write(D%unit) "scalar_flux", lf
        D%number_of_arrays = D%number_of_arrays + 1
      end if
    end if

    if (D%flags%U==1) then
      write(D%unit) vector_flag
      write(D%unit) "u", lf
      D%number_of_arrays = D%number_of_arrays + 1
    end if

    if (D%flags%vorticity==1) then
      write(D%unit) vector_flag
      write(D%unit) "vorticity", lf
      D%number_of_arrays = D%number_of_arrays + 1
    end if
    close(D%unit)
  end subroutine  TFrameDomain_ParallelIO_SaveHeader
  

  
  recursive subroutine SaveBuffers_Parallel(Dptr) bind(C)
    use iso_c_binding, only: c_f_pointer
    use custom_par
    type(c_ptr),value :: Dptr
    type(TFrameDomain_ParallelIO), pointer :: D
    character(2512) :: file_name
    character(70) :: str
    character(8) ::  scalname
    integer :: sc, ie
    integer(MPI_OFFSET_KIND) :: pos
    integer :: view
    integer, allocatable :: requests(:)
    integer :: ireq

    scalname = "scalar00"
    file_name = ""
    
    call c_f_pointer(Dptr, D)

    write(file_name,'(a,i0,a)') trim(D%base_name)//"-",D%frame_number,trim(D%suffix)

    ! Open the file only if it does not already exist
    call MPI_File_open(D%comm, file_name, IOR(IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY), MPI_MODE_EXCL), MPI_INFO_NULL, D%unit, ie)
    
    ! Deal with an already existing file 
    if (ie/=0) then
      call MPI_Barrier(D%comm, ie)
      
      if (master) call MPI_File_delete(file_name,MPI_INFO_NULL, ie)
      if (ie/=0) call error_stop("Error calling MPI_file_delete."//"file_name: '" &
                                 //file_name//"' Error code:", ie)
      
      call MPI_Barrier(D%comm, ie)
      
      call MPI_File_open(D%comm,file_name, IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY), MPI_INFO_NULL, D%unit, ie)
      if (ie/=0) call error_stop("Error calling MPI_file_open, second try after delete."//"file_name: '" &
                                 //file_name//"' Error code:" , ie)
    end if
    
#if NONBLOCKING_MPI_IO    
    allocate(requests(D%number_of_arrays))
    requests = MPI_REQUEST_NULL
    ireq = 0
#endif    
    
    pos = 0
    view = 0

    if (D%flags%Pr==1) then
      call save_scalar(D%Pr)
    end if

    if (D%flags%lambda2==1) then
      call save_scalar(D%lambda2)
    end if

    if (D%flags%scalars==1) then
      scalname(1:6)="scalar"
      do sc = 1,num_of_scalars
        call save_scalar(D%Scalar(:,:,:,sc))
      end do
    elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
      call save_scalar(D%Scalar(:,:,:,1))
    end if

    if (enable_buoyancy.and.D%flags%temperature==1) then
      call save_scalar(D%Temperature)
    end if

    if (enable_moisture.and.D%flags%moisture==1) then
      call save_scalar(D%Moisture)
    end if

    if (enable_liquid.and.D%flags%liquid_water==1) then
      call save_scalar(D%LiquidWater)
    end if

    if (enable_buoyancy.and.D%flags%temperature_flux==1) then
      call save_scalar(D%TemperatureFl)
    end if

    if (enable_buoyancy.and.D%flags%moisture_flux==1) then
      call save_scalar(D%MoistureFl)
    end if

    if (D%flags%scalar_flux==1) then
      if (D%flags%scalars==1) then
        scalname(1:6)="scalar_flux"
        do sc = 1,num_of_scalars
          call save_scalar(D%ScalarFl(:,:,:,sc))
        end do
      elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
        call save_scalar(D%ScalarFl(:,:,:,1))
      end if
    end if

    if (D%flags%U==1) then
      call save_vector(D%U)
    end if

    if (D%flags%vorticity==1) then
      call save_vector(D%vorticity)
    end if
    
#ifdef NONBLOCKING_MPI_IO    
    call MPI_Waitall(size(requests), requests, MPI_STATUSES_IGNORE, ie)
#endif
    
    call MPI_File_close(D%unit, ie)
  contains
    subroutine save_scalar(sc)
      real(real32), intent(in), contiguous :: sc(:,:,:)
      if (view/=1) then
        call MPI_File_set_view(D%unit, pos, D%par_datatype_scalar, &
                               D%par_filetype_scalar, "native", MPI_INFO_NULL, ie)
        view = 1
      end if
      
#ifdef NONBLOCKING_MPI_IO    
      ireq = ireq + 1
      call MPI_File_iwrite_all(D%unit, sc, 1, D%par_datatype_scalar, requests(ireq), ie)
#else
      call MPI_File_write_all(D%unit, sc, 1, D%par_datatype_scalar, MPI_STATUS_IGNORE, ie)
#endif

      pos = pos + D%glob_scalar_storage_size
    end subroutine
    
    subroutine save_vector(vec)
      real(real32), intent(in), contiguous :: vec(:,:,:,:)
      
      if (view/=2) then
        call MPI_File_set_view(D%unit, pos, D%par_datatype_vector, &
                               D%par_filetype_vector, "native", MPI_INFO_NULL, ie)
        view = 2
      end if

#ifdef NONBLOCKING_MPI_IO    
      ireq = ireq + 1
      call MPI_File_iwrite_all(D%unit, vec, 1, D%par_datatype_vector, requests(ireq), ie)
#else
      call MPI_File_write_all(D%unit, vec, 1, D%par_datatype_vector, MPI_STATUS_IGNORE, ie)
#endif

      pos = pos + D%glob_vector_storage_size
    end subroutine
  end subroutine SaveBuffers_Parallel
  

  subroutine TFrameDomain_ParallelIO_DoSave(D)
    use iso_c_binding, only: c_loc, c_funloc
    use stop_procedures, only: error_stop
    class(TFrameDomain_ParallelIO),target,asynchronous,intent(inout) :: D
    integer :: err

!     err = 1
!     
!     select type (Dnp => D)
!       type is (TFrameDomain_ParallelIO)
!         call pthread_create_opaque(Dnp%threadptr, &
!                                    c_funloc(SaveBuffers), &
!                                    c_loc(Dnp), err)
!       class default
!         call error_stop("Error: wrong type of D in TFrameDomain_DoSave.")
!     end select
! 
!     if (err==0) then
!       D%in_progress = .true.
!     else
!       write (*,*) "Error in creating frame thread. Will run again synchronously. Code:",err
      select type (Dnp => D)
        type is (TFrameDomain_ParallelIO)
          call SaveBuffers_Parallel(c_loc(Dnp))
        class default
          call error_stop("Error: wrong type of D in TFrameDomain_DoSave.")
      end select
!     end if

  end subroutine TFrameDomain_ParallelIO_DoSave
  

    

end module Frames_ParallelIO
#endif
