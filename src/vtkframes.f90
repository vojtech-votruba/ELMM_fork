module VTKFrames
  use iso_c_binding, only: c_ptr
  use Kinds
  use Parameters
  use Pthreads
  use Endianness
  use Frames_common
  
  implicit none
  
  private
  
  public TFrameFlags, TFrameDomain, AddDomain, SaveVTKFrames, FinalizeVTKFrames, InitVTKFrames, SaveBuffers
  
  type TFrameFlags
    integer :: U = 1
    integer :: vorticity = 0
    integer :: Pr = 0
    integer :: lambda2 = 0
    integer :: scalars = 1
    integer :: sumscalars = 0
    integer :: temperature = 1
    integer :: moisture = 1
    integer :: liquid_water = 1
    integer :: temperature_flux = 0
    integer :: moisture_flux = 0
    integer :: scalar_flux = 0
  end type
  
  type, extends(TFrameBase) :: TFrameDomain
    integer   :: dimension,direction
    real(knd) :: position
    
    type(TFrameFlags) :: flags
    
    !big endian copies of coordinates
    real(real32),allocatable,dimension(:) :: xPr_be, yPr_be, zPr_be
    
    real(real32), allocatable :: Pr(:,:,:)
    real(real32), allocatable :: U(:,:,:,:)
    real(real32), allocatable :: Temperature(:,:,:), Moisture(:,:,:), LiquidWater(:,:,:), Scalar(:,:,:,:)
    real(real32), allocatable :: TemperatureFl(:,:,:), MoistureFl(:,:,:), ScalarFl(:,:,:,:)
    real(real32), allocatable :: Lambda2(:,:,:), Vorticity(:,:,:,:)
    
    character(4)  :: suffix = ".vtk"
  contains
    procedure :: InDomain => TFrameDomain_InDomain
    procedure :: Fill => TFrameDomain_Fill
    procedure :: DoSave => TFrameDomain_DoSave
    procedure :: SetRest => TFrameDomain_SetRest
    procedure :: Finalize => TFrameDomain_Finalize
  end type TFrameDomain
  
  interface TFrameDomain
    module procedure TFrameDomain_Init
  end interface
  
  interface AddDomain
    module procedure AddFrameDomain
  end interface
  
  type(TFrameDomain), allocatable, target :: FrameDomains(:)
  
  character, parameter :: lf = achar(10)

  integer :: frame_unit = 1999000
  
contains

  
  

  subroutine InitVTKFrames
    integer :: i
  
    if (allocated(FrameDomains)) then
      FrameDomains = pack(FrameDomains, FrameDomains%InDomain())
      
      do i=1, size(FrameDomains)
        call FrameDomains(i)%SetRest(num_of_scalars)
      end do
    end if
   
  end subroutine

  

  subroutine add_element_fd(a,e)
    type(TFrameDomain),allocatable,intent(inout) :: a(:)
    type(TFrameDomain),intent(in) :: e
    type(TFrameDomain),allocatable :: tmp(:)

    if (.not.allocated(a)) then
      a = [e]
    else
      call move_alloc(a,tmp)
      allocate(a(size(tmp)+1))
      a(1:size(tmp)) = tmp
      a(size(tmp)+1) = e
    end if
  end subroutine
  
  
    
  subroutine AddFrameDomain(D)
    type(TFrameDomain),intent(in) :: D
    
    if (allocated(FrameDomains)) then
      if (any(FrameDomains%base_name==D%base_name)) &
        call error_stop("Error, duplicate frame basename '"//trim(D%base_name)//"'. &
                        &Check for duplicate frame labels.")
    end if

    call add_element_fd(FrameDomains, D)

  end subroutine
   
  
  
  subroutine SaveVTKFrames(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    real(knd),intent(in) :: time
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),contiguous,intent(in) :: Pr(-1:,-1:,-1:), &
                                       Temperature(-2:,-2:,-2:), Viscosity(-1:,-1:,-1:), &
                                       Moisture(-2:,-2:,-2:), Scalar(-2:,-2:,-2:,1:)
    integer :: i
    
    if (allocated(FrameDomains)) then
      do i=1,size(FrameDomains)
        call FrameDomains(i)%Save(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
      end do
    end if
  end subroutine

  subroutine FinalizeVTKFrames
    integer :: i
    if (allocated(FrameDomains)) then
      do i=1,size(FrameDomains)
        call FrameDomains(i)%Finalize
      end do
      deallocate(FrameDomains)
    end if
  end subroutine  
  
  
  
  
  
  
  
  function TFrameDomain_Init(label,dimension,direction,position,time_params,frame_flags,bbox) result(D)
    type(TFrameDomain) :: D
    character(*) :: label
    integer,intent(in) :: dimension,direction
    real(knd),intent(in) :: position
    type(TFrameTimes),intent(in) :: time_params
    type(TFrameFlags),intent(in) :: frame_flags
    type(bounding_box), intent(in), optional :: bbox

    D%base_name = trim(output_dir)//"frame-"//label

    D%frame_times = time_params

    D%flags = frame_flags

    D%unit = frame_unit
    frame_unit = frame_unit + 1

    D%dimension = dimension
    D%direction = direction
    D%position = position
    
    if (present(bbox)) D%bbox = bbox

  end function
  
  
 
  

  subroutine TFrameDomain_SetRest(D, num_of_scalars)
    use Boundaries, only: GridCoords
    class(TFrameDomain),intent(inout) :: D
    integer, intent(in) :: num_of_scalars
    integer :: mini,maxi,minj,maxj,mink,maxk
    integer :: bb_mini,bb_maxi,bb_minj,bb_maxj,bb_mink,bb_maxk

    if (D%dimension==3) then

      D%minPri = 1
      D%maxPri = Prnx
      D%minPrj = 1
      D%maxPrj = Prny
      D%minPrk = 1
      D%maxPrk = Prnz
      
      if (allocated(D%bbox)) then
        call GridCoords(bb_mini, bb_minj, bb_mink, &
                        D%bbox%xmin, D%bbox%ymin, D%bbox%zmin)
        call GridCoords(bb_maxi, bb_maxj, bb_maxk, &
                        D%bbox%xmax, D%bbox%ymax, D%bbox%zmax)
                        
        D%minPri = max(D%minPri, bb_mini)
        D%minPrj = max(D%minPrj, bb_minj)
        D%minPrk = max(D%minPrk, bb_mink)
        D%maxPri = min(D%maxPri, bb_maxi)
        D%maxPrj = min(D%maxPrj, bb_maxj)
        D%maxPrk = min(D%maxPrk, bb_maxk)
      end if

    else

      if (D%direction==1) then

        call GridCoords(D%minPri, D%minPrj, D%minPrk, D%position, &
                        (yV(Prny+1)+yV(0))/2._knd, (zW(Prnz+1)+zW(0))/2._knd )

        D%maxPri = D%minPri
        D%minPrj = 1
        D%maxPrj = Prny
        D%minPrk = 1
        D%maxPrk = Prnz

      elseif (D%direction==2) then

        call GridCoords(D%minPri, D%minPrj, D%minPrk, (xU(Prnx+1)+xU(0))/2._knd, &
                        D%position, (zW(Prnz+1)+zW(0))/2._knd )

        D%maxPrj = D%minPrj
        D%minPri = 1
        D%maxPri = Prnx
        D%minPrk = 1
        D%maxPrk = Prnz

      else

        call GridCoords(D%minPri, D%minPrj, D%minPrk, (xU(Prnx+1)+xU(0))/2._knd, &
                        (yV(Prny+1)+yV(0))/2._knd, D%position )

        D%maxPrk = D%minPrk
        D%minPri = 1
        D%maxPri = Prnx
        D%minPrj = 1
        D%maxPrj = Prny

      end if

    end if
    
    D%sizePr = max( (D%maxPri - D%minPri + 1) * (D%maxPrj - D%minPrj + 1) * (D%maxPrk - D%minPrk + 1) , 0 )
    
    mini = D%minPri
    minj = D%minPrj
    mink = D%minPrk
    maxi = D%maxPri
    maxj = D%maxPrj
    maxk = D%maxPrk

    D%xPr = xPr(mini:maxi)
    D%yPr = yPr(minj:maxj)
    D%zPr = zPr(mink:maxk)
    
    D%xPr_be = BigEnd(real(D%xPr, real32))
    D%yPr_be = BigEnd(real(D%yPr, real32))
    D%zPr_be = BigEnd(real(D%zPr, real32))
    
    if (D%flags%Pr==1) then
      allocate(D%Pr(mini:maxi,minj:maxj,mink:maxk))
    end if

    if (D%flags%lambda2==1) then
      allocate(D%lambda2(mini:maxi,minj:maxj,mink:maxk))
    end if

    if (D%flags%scalars==1) then
      allocate(D%Scalar(mini:maxi,minj:maxj,mink:maxk,num_of_scalars))
    elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
      allocate(D%Scalar(mini:maxi,minj:maxj,mink:maxk,1))
    end if

    if (enable_buoyancy.and.D%flags%temperature==1) then
      allocate(D%Temperature(mini:maxi,minj:maxj,mink:maxk))
    end if

    if (enable_moisture.and.D%flags%moisture==1) then
      allocate(D%Moisture(mini:maxi,minj:maxj,mink:maxk))
    end if

    if (enable_liquid.and.D%flags%liquid_water==1) then
      allocate(D%LiquidWater(mini:maxi,minj:maxj,mink:maxk))
    end if

    if (enable_buoyancy.and.D%flags%temperature_flux==1) then
      allocate(D%TemperatureFl(mini:maxi,minj:maxj,mink:maxk))
    end if

    if (enable_moisture.and.D%flags%moisture_flux==1) then
      allocate(D%MoistureFl(mini:maxi,minj:maxj,mink:maxk))
    end if

    if (D%flags%scalar_flux==1) then
      if (D%flags%scalars==1) then
        allocate(D%ScalarFl(mini:maxi,minj:maxj,mink:maxk,num_of_scalars))
      elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
        allocate(D%ScalarFl(mini:maxi,minj:maxj,mink:maxk,1))
      end if
    end if

    if (D%flags%U==1) then
      allocate(D%U(3,mini:maxi,minj:maxj,mink:maxk))
    end if

    if (D%flags%vorticity==1) then
      allocate(D%vorticity(3,mini:maxi,minj:maxj,mink:maxk))
    end if
  end subroutine TFrameDomain_SetRest
  
  
  elemental logical function TFrameDomain_InDomain(D) result(res)
    class(TFrameDomain), intent(in) :: D
    
    if (D%dimension==3) then
      if (allocated(D%bbox)) then
        res = (D%bbox%xmin < im_xmax) .and. &
              (D%bbox%ymin < im_ymax) .and. &
              (D%bbox%zmin < im_zmax) .and. &
              (D%bbox%xmax > im_xmin) .and. &
              (D%bbox%ymax > im_ymin) .and. &
              (D%bbox%zmax > im_zmin)
      else
        res = .true.
      end if
    else if (D%dimension==2) then
      if (D%direction==1) then
        res = (D%position>=xU(0) .and. &
               (D%position<xU(Prnx) .or. &
                (D%position==xU(Prnx).and.Btype(Ea)/=BC_MPI_BOUNDARY)))
      else if (D%direction==2) then   
        res = (D%position>=yV(0) .and. &
               (D%position<yV(Prny) .or. &
                (D%position==yV(Prny).and.Btype(No)/=BC_MPI_BOUNDARY)))
      else if (D%direction==3) then   
        res = (D%position>=zW(0) .and. &
               (D%position<zW(Prnz) .or. &
                (D%position==zW(Prnz).and.Btype(To)/=BC_MPI_BOUNDARY)))
      else
        res = .false.
      end if
    else
      res = .false.
    end if
  end function
  
  
  
  
  

  subroutine TFrameDomain_Fill(D, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    use WaterThermodynamics, only: LiquidWater
    use Output_helpers, only: Lambda2, ScalarVerticalFlux
    !Fill the output buffers for asynchronous output
    class(TFrameDomain),intent(inout) :: D
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),contiguous,intent(in) :: Pr(-1:,-1:,-1:), Viscosity(-1:,-1:,-1:), &
                                       Temperature(-2:,-2:,-2:), Moisture(-2:,-2:,-2:), &
                                       Scalar(-2:,-2:,-2:,1:)
    integer :: i,j,k,l
    integer :: mini,maxi,minj,maxj,mink,maxk

    mini = D%minPri
    minj = D%minPrj
    mink = D%minPrk
    maxi = D%maxPri
    maxj = D%maxPrj
    maxk = D%maxPrk

    if (D%flags%Pr==1) then
      D%Pr = real(Pr(mini:maxi,minj:maxj,mink:maxk), real32)
    end if

    if (D%flags%lambda2==1) then
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
          if (Prtype(i,j,k)<=0) then
            D%lambda2(i,j,k) = real(Lambda2(i,j,k,U,V,W), real32)
          else
            D%lambda2(i,j,k) = 0
          end if
        end do
       end do
      end do

    end if

    if (D%flags%scalars==1) then
      do l = 1,num_of_scalars
       do k = mink,maxk
        do j = minj,maxj
         do i = mini,maxi
           if (Prtype(i,j,k)<=0) then
             D%Scalar(i,j,k,l) = real(Scalar(i,j,k,l), real32)
           else
             D%Scalar(i,j,k,l) = 0
           end if
         end do
        end do
       end do
      end do
    elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
          if (Prtype(i,j,k)<=0) then
            D%Scalar(i,j,k,1) = real(sum(Scalar(i,j,k,:)), real32)
          else
            D%Scalar(i,j,k,1) = 0
          end if
        end do
       end do
      end do
    end if

    if (enable_buoyancy.and.D%flags%temperature==1) then
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
          if (Prtype(i,j,k)<=0) then
            D%Temperature(i,j,k) = real(Temperature(i,j,k), real32)
          else
            D%Temperature(i,j,k) = real(temperature_ref, real32)
          end if
        end do
       end do
      end do
    end if

    if (enable_moisture.and.D%flags%moisture==1) then
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
          if (Prtype(i,j,k)<=0) then
            D%Moisture(i,j,k) = real(Moisture(i,j,k), real32)
          else
            D%Moisture(i,j,k) = real(moisture_ref, real32)
          end if
        end do
       end do
      end do
    end if

    if (enable_liquid.and.D%flags%liquid_water==1) then
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
          if (Prtype(i,j,k)<=0) then
            D%LiquidWater(i,j,k) = real(LiquidWater(i,j,k), real32)
          else
            D%LiquidWater(i,j,k) = real(0, real32)
          end if
        end do
       end do
      end do
    end if

    if (enable_buoyancy.and.D%flags%temperature_flux==1) then
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
          if (Prtype(i,j,k)<=0) then
            D%TemperatureFl(i,j,k) = real(ScalarVerticalFlux(i,j,k,Temperature,W), real32)
          else
            D%TemperatureFl(i,j,k) = 0
          end if
        end do
       end do
      end do
    end if

    if (enable_moisture.and.D%flags%moisture_flux==1) then
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
          if (Prtype(i,j,k)<=0) then
            D%MoistureFl(i,j,k) = real(ScalarVerticalFlux(i,j,k,Moisture,W), real32)
          else
            D%MoistureFl(i,j,k) = 0
          end if
        end do
       end do
      end do
    end if

    if (D%flags%scalar_flux==1) then
      if (D%flags%scalars==1) then
        do l = 1,num_of_scalars
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
              if (Prtype(i,j,k)<=0) then
                D%ScalarFl(i,j,k,l) = real( ScalarVerticalFlux(i,j,k,Scalar(:,:,:,l),W) , real32)
              else
                D%ScalarFl(i,j,k,l) = 0
              end if
            end do
           end do
          end do
        end do
      elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
        do k = mink,maxk
         do j = minj,maxj
          do i = mini,maxi
            if (Prtype(i,j,k)<=0) then
              D%ScalarFl(i,j,k,1) = real(ScalarVerticalFlux(i,j,k,sum(Scalar,4),W) , real32)
            else
              D%ScalarFl(i,j,k,1) = 0
            end if
          end do
         end do
        end do
      end if
    end if

    if (D%flags%U==1) then
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
          if (Prtype(i,j,k)<=0) then
            D%U(1,i,j,k) = real( (U(i,j,k)+U(i-1,j,k))/2, real32)
            D%U(2,i,j,k) = real( (V(i,j,k)+V(i,j-1,k))/2, real32)
            D%U(3,i,j,k) = real( (W(i,j,k)+W(i,j,k-1))/2, real32)
          else
            D%U(:,i,j,k) = 0
          end if
        end do
       end do
      end do
    end if

    if (D%flags%vorticity==1) then
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
          if (Prtype(i,j,k)<=0) then
            D%vorticity(1,i,j,k) = real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin) &
                                     - (V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin), real32)
            D%vorticity(2,i,j,k) = real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin) &
                                     - (W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin), real32)
            D%vorticity(3,i,j,k) = real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin) &
                                     - (U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin), real32)
          else
            D%vorticity(:,i,j,k) = 0
          end if
        end do
       end do
      end do
    end if
  end subroutine TFrameDomain_Fill
  
  
  recursive subroutine SaveBuffers(Dptr) bind(C)
    use iso_c_binding, only: c_f_pointer
    type(c_ptr),value :: Dptr
    type(TFrameDomain),pointer :: D
    character(2512) :: file_name
    character(70) :: str
    character(8) ::  scalname
    integer :: mini, maxi, minj, maxj, mink, maxk
    integer :: sc

    scalname = "scalar00"
    file_name = ""
    
    call c_f_pointer(Dptr, D)

    write(file_name,'(a,i0,a)') trim(D%base_name)//"-",D%frame_number,trim(D%suffix)

    if (littleendian) then
      if (allocated(D%Pr))            D%Pr = SwapB(D%Pr)
      if (allocated(D%U))             D%U = SwapB(D%U)
      if (allocated(D%Temperature))   D%Temperature = SwapB(D%Temperature)
      if (allocated(D%Moisture))      D%Moisture = SwapB(D%Moisture)
      if (allocated(D%LiquidWater))   D%LiquidWater = SwapB(D%LiquidWater)
      if (allocated(D%Scalar))        D%Scalar = SwapB(D%Scalar)
      if (allocated(D%TemperatureFl)) D%TemperatureFl = SwapB(D%TemperatureFl)
      if (allocated(D%MoistureFl))    D%MoistureFl = SwapB(D%MoistureFl)
      if (allocated(D%ScalarFl))      D%ScalarFl = SwapB(D%ScalarFl)
      if (allocated(D%Lambda2))       D%Lambda2 = SwapB(D%Lambda2)
      if (allocated(D%Vorticity))     D%Vorticity = SwapB(D%Vorticity)
    end if

    mini = D%minPri
    minj = D%minPrj
    mink = D%minPrk
    maxi = D%maxPri
    maxj = D%maxPrj
    maxk = D%maxPrk

    open(D%unit,file = file_name, &
      access='stream',status='replace',form="unformatted",action="write")

    write(D%unit) "# vtk DataFile Version 2.0", lf
    write(D%unit) "CLMM output file", lf
    write(D%unit) "BINARY", lf
    write(D%unit) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),*) maxi-mini+1,maxj-minj+1,maxk-mink+1
    write(D%unit) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') maxi-mini+1,"float"
    write(D%unit) str, lf
    write(D%unit) D%xPr_be, lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') maxj-minj+1,"float"
    write(D%unit) str, lf
    write(D%unit) D%yPr_be, lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') maxk-mink+1,"float"
    write(D%unit) str, lf
    write(D%unit) D%zPr_be, lf
    str="POINT_DATA"
    write(str(12:),*) D%sizePr
    write(D%unit) str, lf

    if (D%flags%Pr==1) then
      write(D%unit) "SCALARS p float", lf
      write(D%unit) "LOOKUP_TABLE default", lf
      write(D%unit) D%Pr, lf
    end if

    if (D%flags%lambda2==1) then
      write(D%unit) "SCALARS lambda2 float", lf
      write(D%unit) "LOOKUP_TABLE default", lf
      write(D%unit) D%lambda2, lf
    end if

    if (D%flags%scalars==1) then
      scalname(1:6)="scalar"
      do sc = 1,num_of_scalars
       write(scalname(7:8),"(I2.2)") sc
       write(D%unit) "SCALARS ", scalname , " float", lf
       write(D%unit) "LOOKUP_TABLE default", lf
       write(D%unit) D%Scalar(:,:,:,sc), lf
      end do
    elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
      write(D%unit) "SCALARS ", "scalar" , " float", lf
      write(D%unit) "LOOKUP_TABLE default", lf
      write(D%unit) D%Scalar(:,:,:,1), lf
    end if

    if (enable_buoyancy.and.D%flags%temperature==1) then
      write(D%unit) "SCALARS temperature float", lf
      write(D%unit) "LOOKUP_TABLE default", lf
      write(D%unit) D%Temperature, lf
    end if

    if (enable_moisture.and.D%flags%moisture==1) then
      write(D%unit) "SCALARS moisture float", lf
      write(D%unit) "LOOKUP_TABLE default", lf
      write(D%unit) D%Moisture, lf
    end if

    if (enable_liquid.and.D%flags%liquid_water==1) then
      write(D%unit) "SCALARS liquid_water float", lf
      write(D%unit) "LOOKUP_TABLE default", lf
      write(D%unit) D%LiquidWater, lf
    end if

    if (enable_buoyancy.and.D%flags%temperature_flux==1) then
      write(D%unit) "SCALARS temperature_flux float", lf
      write(D%unit) "LOOKUP_TABLE default", lf
      write(D%unit) D%TemperatureFl, lf
    end if

    if (enable_buoyancy.and.D%flags%moisture_flux==1) then
      write(D%unit) "SCALARS moisture_flux float", lf
      write(D%unit) "LOOKUP_TABLE default", lf
      write(D%unit) D%MoistureFl, lf
    end if

    if (D%flags%scalar_flux==1) then
      if (D%flags%scalars==1) then
        scalname(1:6)="scalar_flux"
        do sc = 1,num_of_scalars
          write(scalname(7:8),"(I2.2)") sc
          write(D%unit) "SCALARS ", scalname , " float", lf
          write(D%unit) "LOOKUP_TABLE default", lf
          write(D%unit) D%ScalarFl(:,:,:,sc), lf
        end do
      elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
        write(D%unit) "SCALARS scalar_flux float", lf
        write(D%unit) "LOOKUP_TABLE default", lf
        write(D%unit) D%ScalarFl(:,:,:,1), lf
      end if
    end if

    if (D%flags%U==1) then
      write(D%unit) "VECTORS u float", lf
      write(D%unit) D%U
      write(D%unit) lf
    end if

    if (D%flags%vorticity==1) then
      write(D%unit) "VECTORS vorticity float", lf
      write(D%unit) D%vorticity, lf
    end if
    close(D%unit)

  end subroutine SaveBuffers
  
  
  subroutine TFrameDomain_DoSave(D)
    use iso_c_binding, only: c_loc, c_funloc
    use stop_procedures, only: error_stop
    class(TFrameDomain),target,asynchronous,intent(inout) :: D
    integer :: err

!     err = 1
!     
!     select type (Dnp => D)
!       type is (TFrameDomain)
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
        type is (TFrameDomain)
          call SaveBuffers(c_loc(Dnp))
        class default
          call error_stop("Error: wrong type of D in TFrameDomain_DoSave.")
      end select
!     end if

  end subroutine TFrameDomain_DoSave


  
  subroutine TFrameDomain_Finalize(D)
    class(TFrameDomain),intent(inout) :: D

    if (D%in_progress) call D%Wait

    call D%SaveTimes

  end subroutine

  
  
end module VTKFrames
