module SurfaceFrames
  use Parameters
  use Frames_common
  use VTKFrames
  
  implicit none
  
  !first version simply assumes that only the images with kim==1 are concerned
  
  private
  
  public AddSurfaceFrameDomain, TSurfaceFrameDomain, SaveSurfaceFrames, FinalizeSurfaceFrames, InitSurfaceFrames

  type, extends(TFrameDomain) :: TSurfaceFrameDomain
    integer :: number_of_body = 1
    integer, allocatable :: k(:,:)
  contains
    procedure :: Fill => TSurfaceFrameDomain_Fill
    procedure :: DoSave => TSurfaceFrameDomain_DoSave
    procedure :: SetRest => TSurfaceFrameDomain_SetRest
  end type
  
  interface AddDomain
    module procedure AddSurfaceFrameDomain
  end interface
  
  interface TSurfaceFrameDomain
    module procedure TSurfaceFrameDomain_Init
  end interface
  
  type(TSurfaceFrameDomain), allocatable :: SurfaceFrameDomains(:)
  
  integer :: frame_unit = 2999000
  
contains

  elemental logical function InDomain(D) result(res)
    type(TSurfaceFrameDomain), intent(in) :: D
    !assume kim == 1
    res = (offset_to_global_z == 0)
  end function
  
  subroutine InitSurfaceFrames
    integer :: i
  
    if (allocated(SurfaceFrameDomains)) then
      SurfaceFrameDomains = pack(SurfaceFrameDomains, InDomain(SurfaceFrameDomains))
      
      do i=1, size(SurfaceFrameDomains)
        call SurfaceFrameDomains(i)%SetRest(num_of_scalars)
      end do
    end if
  
  end subroutine


  subroutine add_element_fd(a,e)
    type(TSurfaceFrameDomain),allocatable,intent(inout) :: a(:)
    type(TSurfaceFrameDomain),intent(in) :: e
    type(TSurfaceFrameDomain),allocatable :: tmp(:)

    if (.not.allocated(a)) then
      a = [e]
    else
      call move_alloc(a,tmp)
      allocate(a(size(tmp)+1))
      a(1:size(tmp)) = tmp
      a(size(tmp)+1) = e
    end if
  end subroutine
    
  subroutine AddSurfaceFrameDomain(D)
    type(TSurfaceFrameDomain),intent(in) :: D

    call add_element_fd(SurfaceFrameDomains, D)

  end subroutine

  subroutine TSurfaceFrameDomain_SetRest(D, num_of_scalars)
    use Boundaries, only: GridCoords
    use Endianness, only: BigEnd
    class(TSurfaceFrameDomain),intent(inout) :: D
    integer, intent(in) :: num_of_scalars
    integer :: mini,maxi,minj,maxj
    integer :: i, j, k


    D%maxPrk = 1
    D%maxPrk = 1
    D%minPri = 1
    D%maxPri = Prnx
    D%minPrj = 1
    D%maxPrj = Prny
    D%minPrk = 1
    D%maxPrk = 1

        
    
    D%sizePr = max( (D%maxPri - D%minPri + 1) * (D%maxPrj - D%minPrj + 1) , 0 )
    
    mini = D%minPri
    minj = D%minPrj

    maxi = D%maxPri
    maxj = D%maxPrj

    D%xPr = xPr(mini:maxi)
    D%yPr = yPr(minj:maxj)
    D%zPr = zPr(1:1)
    
    D%xPr_be = BigEnd(real(D%xPr, real32))
    D%yPr_be = BigEnd(real(D%yPr, real32))
    D%zPr_be = BigEnd(real(D%zPr, real32))
    
    if (D%flags%Pr==1) then
      allocate(D%Pr(mini:maxi,minj:maxj,1))
    end if

    if (D%flags%lambda2==1) then
      allocate(D%lambda2(mini:maxi,minj:maxj,1))
    end if

    if (D%flags%scalars==1) then
      allocate(D%Scalar(mini:maxi,minj:maxj,1,num_of_scalars))
    elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
      allocate(D%Scalar(mini:maxi,minj:maxj,1,1))
    end if

    if (enable_buoyancy.and.D%flags%temperature==1) then
      allocate(D%Temperature(mini:maxi,minj:maxj,1))
    end if

    if (enable_moisture.and.D%flags%moisture==1) then
      allocate(D%Moisture(mini:maxi,minj:maxj,1))
    end if

    if (enable_buoyancy.and.D%flags%temperature_flux==1) then
      allocate(D%TemperatureFl(mini:maxi,minj:maxj,1))
    end if

    if (enable_moisture.and.D%flags%moisture_flux==1) then
      allocate(D%MoistureFl(mini:maxi,minj:maxj,1))
    end if

    if (D%flags%scalar_flux==1) then
      if (D%flags%scalars==1) then
        allocate(D%ScalarFl(mini:maxi,minj:maxj,1,num_of_scalars))
      elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
        allocate(D%ScalarFl(mini:maxi,minj:maxj,1,1))
      end if
    end if

    if (D%flags%U==1) then
      allocate(D%U(3,mini:maxi,minj:maxj,1))
    end if

    if (D%flags%vorticity==1) then
      allocate(D%vorticity(3,mini:maxi,minj:maxj,1))
    end if
    
    allocate(D%k(mini:maxi,minj:maxj))
    
    do j = minj ,maxj
      do i = mini, maxi
        k = 1
        do
          if ((Prtype(i,j,k)>0.and.Prtype(i,j,k)/=D%number_of_body) .or. k>Prnz) then
            D%k(i,j) = -1
            exit
          else if (Prtype(i,j,k)<=0) then
            D%k(i,j) = k
            exit
          else
            k = k + 1
          end if
        end do
      end do
    end do
          
  end subroutine TSurfaceFrameDomain_SetRest
  
  
  
  
  subroutine SaveSurfaceFrames(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    real(knd),intent(in) :: time
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),contiguous,intent(in) :: Pr(-1:,-1:,-1:), &
                                       Temperature(-1:,-1:,-1:), Viscosity(-1:,-1:,-1:), &
                                       Moisture(-1:,-1:,-1:), Scalar(-1:,-1:,-1:,1:)
    integer :: i
    
    if (allocated(SurfaceFrameDomains)) then
      do i=1,size(SurfaceFrameDomains)
        call SurfaceFrameDomains(i)%Save(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
      end do
    end if
  end subroutine
  

  subroutine FinalizeSurfaceFrames
    integer :: i
    if (allocated(SurfaceFrameDomains)) then
      do i=1,size(SurfaceFrameDomains)
        call SurfaceFrameDomains(i)%Finalize
      end do
      deallocate(SurfaceFrameDomains)
    end if
  end subroutine  
  
  
  
  
  
  
  function TSurfaceFrameDomain_Init(label,dimension,direction,position,time_params,frame_flags) result(D)
    type(TSurfaceFrameDomain) :: D
    character(*) :: label
    integer,intent(in) :: dimension,direction
    real(knd),intent(in) :: position
    type(TFrameTimes),intent(in) :: time_params
    type(TFrameFlags),intent(in) :: frame_flags

    D%TFrameDomain = TFrameDomain(label,dimension,direction,position,time_params,frame_flags)

  end function
  
  
  
  

  subroutine TSurfaceFrameDomain_Fill(D, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    use Output_helpers, only: Lambda2, ScalarVerticalFlux
    !Fill the output buffers for asynchronous output
    class(TSurfaceFrameDomain),intent(inout) :: D
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),contiguous,intent(in) :: Pr(-1:,-1:,-1:), Viscosity(-1:,-1:,-1:), &
                                       Temperature(-1:,-1:,-1:), Moisture(-1:,-1:,-1:), &
                                       Scalar(-1:,-1:,-1:,1:)
    integer :: i,j,l
    integer :: mini,maxi,minj,maxj

    mini = D%minPri
    minj = D%minPrj

    maxi = D%maxPri
    maxj = D%maxPrj

    if (D%flags%Pr==1) then
      do j = minj, maxj
        do i = mini, maxi
          if (D%k(i,j)>0) then
            D%Pr(i,j,1) = real(Pr(i,j,D%k(i,j)), real32)
          else
            D%Pr(i,j,1) = -huge(1._real32)
          end if
        end do
      end do
    end if

    if (D%flags%lambda2==1) then
      do j = minj, maxj
        do i = mini, maxi
          if (D%k(i,j)>0) then
            D%lambda2(i,j,1) = real(Lambda2(i,j,D%k(i,j),U,V,W), real32)
          else
            D%lambda2(i,j,1) = -huge(1._real32)
          end if
        end do
      end do
    end if

    if (D%flags%scalars==1) then
      do l = 1,num_of_scalars
        do j = minj, maxj
          do i = mini, maxi
            if (D%k(i,j)>0) then
              D%Scalar(i,j,1,l) = real(Scalar(i,j,D%k(i,j),l), real32)
            else
              D%Scalar(i,j,1,l) = -huge(1._real32)
            end if
          end do
        end do
      end do
    elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
      do j = minj, maxj
        do i = mini, maxi
          if (D%k(i,j)>0) then
            D%Scalar(i,j,1,1) = real(sum(Scalar(i,j,D%k(i,j),:)), real32)
          else
            D%Scalar(i,j,1,1) = -huge(1._real32)
          end if
        end do
      end do
    end if

    if (enable_buoyancy.and.D%flags%temperature==1) then
      do j = minj, maxj
        do i = mini, maxi
          if (D%k(i,j)>0) then
            D%Temperature(i,j,1) = real(Temperature(i,j,D%k(i,j)), real32)
          else
            D%Temperature(i,j,1) = -huge(1._real32)
          end if
        end do
      end do
    end if

    if (enable_moisture.and.D%flags%moisture==1) then
      do j = minj, maxj
        do i = mini, maxi
          if (D%k(i,j)>0) then
            D%Moisture(i,j,1) = real(Moisture(i,j,D%k(i,j)), real32)
          else
            D%Moisture(i,j,1) = -huge(1._real32)
          end if
        end do
      end do
    end if

    if (enable_buoyancy.and.D%flags%temperature_flux==1) then
      do j = minj, maxj
        do i = mini, maxi
          if (D%k(i,j)>0) then
            D%TemperatureFl(i,j,1) = real(ScalarVerticalFlux(i,j,D%k(i,j),Temperature,W), real32)
          else
            D%TemperatureFl(i,j,1) = -huge(1._real32)
          end if
        end do
      end do
    end if

    if (enable_moisture.and.D%flags%moisture_flux==1) then
      do j = minj, maxj
        do i = mini, maxi
          if (D%k(i,j)>0) then
            D%MoistureFl(i,j,1) = real(ScalarVerticalFlux(i,j,D%k(i,j),Moisture,W), real32)
          else
            D%MoistureFl(i,j,1) = -huge(1._real32)
          end if
        end do
      end do
    end if

    if (D%flags%scalar_flux==1) then
      if (D%flags%scalars==1) then
        do l = 1,num_of_scalars
          do j = minj, maxj
            do i = mini, maxi
              if (D%k(i,j)>0) then
                D%ScalarFl(i,j,1,l) = real( ScalarVerticalFlux(i,j,D%k(i,j),Scalar(:,:,:,l),W) , real32)
              else
                D%ScalarFl(i,j,1,l) = -huge(1._real32)
              end if
            end do
          end do
        end do
      elseif (D%flags%sumscalars==1.and.num_of_scalars>0) then
        do j = minj, maxj
          do i = mini, maxi
            if (D%k(i,j)>0) then
              D%ScalarFl(i,j,1,1) = real(ScalarVerticalFlux(i,j,D%k(i,j),sum(Scalar,4),W) , real32)
            else
              D%ScalarFl(i,j,1,1) = -huge(1._real32)
            end if
          end do
        end do
      end if
    end if

    if (D%flags%U==1) then
      do j = minj, maxj
        do i = mini, maxi
          if (D%k(i,j)>0) then
            D%U(1,i,j,1) = real( (U(i,j,D%k(i,j))+U(i-1,j,D%k(i,j)))/2, real32)
            D%U(2,i,j,1) = real( (V(i,j,D%k(i,j))+V(i,j-1,D%k(i,j)))/2, real32)
            D%U(3,i,j,1) = real( (W(i,j,D%k(i,j))+W(i,j,D%k(i,j)-1))/2, real32)
          else
            D%U(:,i,j,1) = -huge(1._real32)
          end if
        end do
      end do
    end if

    if (D%flags%vorticity==1) then
      do j = minj, maxj
        do i = mini, maxi
          if (D%k(i,j)>0) then
            D%vorticity(1,i,j,1) = real((W(i,j+1,D%k(i,j))-W(i,j-1,D%k(i,j))+ &
                                         W(i,j+1,D%k(i,j)-1)-W(i,j-1,D%k(i,j)-1)) &
                                        /(4*dxmin) &
                                      - (V(i,j,D%k(i,j)+1)-V(i,j,D%k(i,j)-1)+ &
                                         V(i,j-1,D%k(i,j)+1)-V(i,j-1,D%k(i,j)-1)) &
                                        /(4*dymin), real32)
            D%vorticity(2,i,j,1) = real((U(i,j,D%k(i,j)+1)-U(i,j,D%k(i,j)-1)+&
                                         U(i-1,j,D%k(i,j)+1)-U(i-1,j,D%k(i,j)-1)) &
                                        /(4*dxmin) &
                                      - (W(i+1,j,D%k(i,j))-W(i-1,j,D%k(i,j))+&
                                         W(i+1,j,D%k(i,j)-1)-W(i-1,j,D%k(i,j)-1)) &
                                        /(4*dymin), real32)
            D%vorticity(3,i,j,1) = real((V(i+1,j,D%k(i,j))-V(i-1,j,D%k(i,j))+ &
                                         V(i+1,j-1,D%k(i,j))-V(i-1,j-1,D%k(i,j))) &
                                        /(4*dxmin) &
                                      - (U(i,j+1,D%k(i,j))-U(i,j-1,D%k(i,j))+ &
                                         U(i-1,j+1,D%k(i,j))-U(i-1,j-1,D%k(i,j))) &
                                        /(4*dymin), real32)
          else
            D%vorticity(:,i,j,1) = -huge(1._real32)
          end if
        end do
      end do
    end if
  end subroutine TSurfaceFrameDomain_Fill
  
  
  subroutine TSurfaceFrameDomain_DoSave(D)
    use iso_c_binding, only: c_loc, c_funloc
    use stop_procedures, only: error_stop
    class(TSurfaceFrameDomain),target,asynchronous,intent(inout) :: D
    integer :: err

    err = 1
    
    associate (Dnp => D%TFrameDomain)
      call pthread_create_opaque(Dnp%threadptr, &
                                 c_funloc(SaveBuffers), &
                                 c_loc(Dnp), err)
    end associate

    if (err==0) then
      D%in_progress = .true.
    else
      write (*,*) "Error in creating frame thread. Will run again synchronously. Code:",err
      associate (Dnp => D%TFrameDomain)
        call SaveBuffers(c_loc(Dnp))
      end associate
    end if

  end subroutine TSurfaceFrameDomain_DoSave


end module SurfaceFrames
