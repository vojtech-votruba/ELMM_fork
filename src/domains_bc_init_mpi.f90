module domains_bc_init

  use domains_bc_types
  use custom_par

  implicit  none

  private

  public par_init_domain_boundary_conditions_implementation

contains

  subroutine par_init_domain_boundary_conditions_implementation(domain_bc_send_buffers, &
                                                                domain_bc_recv_buffers)
    type(dom_bc_buffer_copy), allocatable :: domain_bc_send_buffers(:)
    type(dom_bc_buffer_turbulence_generator), allocatable :: domain_bc_recv_buffers(:)
    integer :: xi, dii, djj, dkk
    integer :: num
    integer :: child_domain, i_child_domain, side
    integer :: n_send_buffers
    integer :: err
    integer, allocatable :: requests(:)
    integer :: request

    logical, allocatable :: child_domain_in_this_image(:)

    real(knd), parameter :: eps = 100*epsilon(1._knd)
 

    if (enable_multiple_domains) then

      allocate(requests(0))



      ! 1 Find out which child domains are actually immersed in this image.

      allocate(child_domain_in_this_image(size(child_domains)))

      do i_child_domain = 1, size(child_domains)
        child_domain_in_this_image(i_child_domain) = (im_xmin < domain_grids(i_child_domain)%xmax - eps .and. &
                                                      im_xmax > domain_grids(i_child_domain)%xmin + eps .and. &
                                                      im_ymin < domain_grids(i_child_domain)%ymax - eps .and. &
                                                      im_ymax > domain_grids(i_child_domain)%ymin + eps .and. &
                                                      im_zmin < domain_grids(i_child_domain)%zmax - eps .and. &
                                                      im_zmax > domain_grids(i_child_domain)%zmin + eps)
      end do

      image_child_domains = pack(child_domains, child_domain_in_this_image)


      ! 2 Find out which child domain images are immersed in this image.
      !   They always form 3D rectangular blocks.

      allocate(domain_parent_buffers(size(image_child_domains)))

      do i_child_domain = 1, size(image_child_domains)
        associate(bs => domain_parent_buffers(i_child_domain))

          child_domain = image_child_domains(i_child_domain)
          bs%remote_domain = child_domain

          do dkk = 1, domain_nzims(child_domain)
            do djj = 1, domain_nyims(child_domain)
              do dii = 1, domain_nxims(child_domain)

                !if child domain immersed in this image
                if (im_xmin < domain_grids(child_domain)%xmaxs(dii) - eps .and. &
                    im_xmax > domain_grids(child_domain)%xmins(dii) + eps .and. &
                    im_ymin < domain_grids(child_domain)%ymaxs(djj) - eps .and. &
                    im_ymax > domain_grids(child_domain)%ymins(djj) + eps .and. &
                    im_zmin < domain_grids(child_domain)%zmaxs(dkk) - eps .and. &
                    im_zmax > domain_grids(child_domain)%zmins(dkk) + eps) then

                  bs%iim2 = dii
                  bs%jim2 = djj
                  bs%kim2 = dkk                

                  if (bs%iim1==0) then
                      bs%iim1 = dii
                      bs%jim1 = djj
                      bs%kim1 = dkk
                  end if
                end if

              end do
            end do
          end do

        end associate
      end do


      
      ! 3 Find out the parent image.

      !   Only one parent image shall exist.
      if (parent_domain > 0) then
outer:  do dkk = 1, domain_nzims(parent_domain)
          do djj = 1, domain_nyims(parent_domain)
            do dii = 1, domain_nxims(parent_domain)
              if (im_xmin < domain_grids(parent_domain)%xmaxs(dii) - eps .and. &
                  im_xmax > domain_grids(parent_domain)%xmins(dii) + eps .and. &
                  im_ymin < domain_grids(parent_domain)%ymaxs(djj) - eps .and. &
                  im_ymax > domain_grids(parent_domain)%ymins(djj) + eps .and. &
                  im_zmin < domain_grids(parent_domain)%zmaxs(dkk) - eps .and. &
                  im_zmax > domain_grids(parent_domain)%zmins(dkk) + eps) then
                parent_image = [dii, djj, dkk]
                exit outer
              end if
            end do
          end do
        end do outer
        
        if (any(parent_image<1)) then
          write(*,'(a,i0,a)') "Error, no parent image for domain ",domain_index,". Is the domain properly nested?"
          call error_stop
        end if

      end if
      
      ! 4 Set up the domain parrent buffers for each child image immersed in this image.

      !send buffers should be counted from 1, there may be more of them from several nested domains
      n_send_buffers = 0
      do i_child_domain = 1, size(image_child_domains)
        associate(bs => domain_parent_buffers(i_child_domain))

          child_domain = image_child_domains(i_child_domain)

          ! 4b Create the domain parent buffers.

          allocate(bs%bs( &
                     bs%iim1:bs%iim2, &
                     bs%jim1:bs%jim2, &
                     bs%kim1:bs%kim2))

          do dkk = bs%kim1, bs%kim2
            do djj = bs%jim1, bs%jim2
              do dii = bs%iim1, bs%iim2
                call create_domain_parent_buffer(bs%bs(dii,djj,dkk), &
                                                 child_domain, [dii, djj, dkk])
              end do
            end do
          end do


          ! 4b Compute the number of the parent (send) boundary buffers.

          do side = We, Ea
            if (domain_is_boundary_nested(side,child_domain)) &
              n_send_buffers = n_send_buffers + &
                                 (bs%jim2-bs%jim1+1) * (bs%kim2-bs%kim1+1)
          end do 
          do side = So, No
            if (domain_is_boundary_nested(side,child_domain)) &
              n_send_buffers = n_send_buffers + &
                                 (bs%iim2-bs%iim1+1) * (bs%kim2-bs%kim1+1)
          end do 
          do side = Bo, To
            if (domain_is_boundary_nested(side,child_domain)) &
              n_send_buffers = n_send_buffers + &
                                 (bs%iim2-bs%iim1+1) * (bs%jim2-bs%jim1+1)
          end do 

        end associate
      end do

      
      ! 5 Set up the parent (send) buffers for the domain nesting boundaries
          
      allocate(domain_bc_send_buffers(n_send_buffers))

      num = 0
      do i_child_domain = 1, size(image_child_domains)
        associate(bs => domain_parent_buffers(i_child_domain))

          child_domain = image_child_domains(i_child_domain)

          if (domain_is_boundary_nested(We,child_domain) .and. bs%iim1==1) then   
            do dkk = bs%kim1, bs%kim2
              do djj = bs%jim1, bs%jim2
                num = num + 1

                call create_boundary_parent_buffer(num=num, &
                       dir=We, domain=child_domain, im=[1,djj,dkk])
             end do
            end do
          end if

          if (domain_is_boundary_nested(Ea,child_domain) .and. bs%iim2==domain_nxims(child_domain)) then   
            do dkk = bs%kim1, bs%kim2
              do djj = bs%jim1, bs%jim2
                num = num + 1

                call create_boundary_parent_buffer(num=num, &
                       dir=Ea, domain=child_domain, im=[domain_nxims(child_domain),djj,dkk])
             end do
            end do
          end if

          if (domain_is_boundary_nested(So,child_domain) .and. bs%jim1==1) then   
            do dkk = bs%kim1, bs%kim2
              do dii = bs%iim1, bs%iim2
                num = num + 1

                call create_boundary_parent_buffer(num=num, &
                       dir=So, domain=child_domain, im=[dii,1,dkk])
             end do
            end do
          end if

          if (domain_is_boundary_nested(No,child_domain) .and. bs%jim2==domain_nyims(child_domain)) then   
            do dkk = bs%kim1, bs%kim2
              do dii = bs%iim1, bs%iim2
                num = num + 1

                call create_boundary_parent_buffer(num=num, &
                       dir=No, domain=child_domain, im=[dii,domain_nyims(child_domain),dkk])
             end do
            end do
          end if

         if (domain_is_boundary_nested(Bo,child_domain) .and. bs%kim1==1) then   
             do djj = bs%jim1, bs%jim2
              do dii = bs%iim1, bs%iim2
                num = num + 1

                call create_boundary_parent_buffer(num=num, &
                       dir=Bo, domain=child_domain, im=[dii,djj,1])

             end do
            end do
          end if

          if (domain_is_boundary_nested(To,child_domain) .and. bs%kim2==domain_nzims(child_domain)) then   
            do djj = bs%jim1, bs%jim2
              do dii = bs%iim1, bs%iim2
                num = num + 1

                call create_boundary_parent_buffer(num=num, &
                       dir=To, domain=child_domain, im=[dii,djj,domain_nzims(child_domain)])

             end do
            end do
          end if

        end associate
      end do


      ! 6 Set up the domain child buffer and the boundary child buffers.
      !   Ther shall be only one parent image.

      if (parent_domain>0) then
        allocate(domain_bc_recv_buffers(We:To))

        allocate(domain_child_buffer)

        call create_domain_child_buffer(domain_child_buffer)

        do side = We, To
          if (is_domain_boundary_nested(side).and.is_boundary_domain_boundary(side)) then
            Btype(side) = BC_DOMAIN_COPY
            TempBtype(side) = BC_DOMAIN_COPY
            MoistBtype(side) = BC_DOMAIN_COPY
            !HACK!
            if (side/=We) ScalBtype(side) = BC_DOMAIN_COPY

            call create_boundary_child_buffer(dir=side, domain=parent_domain)
          endif
       end do

      end if


      ! 7 Wait until all synchronization messages are sent or received.

      call MPI_Waitall(size(requests), requests, MPI_STATUSES_IGNORE, err)


      ! 8 Set that the uniform pressure force should be used when its value will be being received from the parent.

      if (parent_domain>0) then
        if (domain_child_buffer%exchange_pr_gradient_x) &
          enable_pr_gradient_x_uniform = .true.

        if (domain_child_buffer%exchange_pr_gradient_y) &
          enable_pr_gradient_y_uniform = .true.

        if (domain_child_buffer%exchange_pr_gradient_z) &
          enable_pr_gradient_z_uniform = .true.
      end if

    end if !enable_multiple_domains


contains

    function x_coord(x) result(res)
      integer :: res
      real(knd), intent(in) :: x
      res = min( max(nint( (x - im_xmin)/dxmin ),0) , Unx+1)
    end function

    function y_coord(y) result(res)
      integer :: res
      real(knd), intent(in) :: y
      res = min( max(nint( (y - im_ymin)/dymin ),0) , Vny+1)
    end function

    function z_coord(z) result(res)
      integer :: res
      real(knd), intent(in) :: z
      res = min( max(nint( (z - im_zmin)/dzmin ),0) , Wnz+1)
    end function

    subroutine create_boundary_parent_buffer(num, dir, domain, im)
      integer :: Ui1, Ui2, Vi1, Vi2, Wi1, Wi2, Pri1, Pri2
      integer :: Uj1, Uj2, Vj1, Vj2, Wj1, Wj2, Prj1, Prj2
      integer :: Uk1, Uk2, Vk1, Vk2, Wk1, Wk2, Prk1, Prk2
      integer, intent(in) :: num, dir, domain, im(3)
      real(knd) :: cxmax, cxmin, cymax, cymin, czmax, czmin
      integer :: cxi1, cxi2, cyj1, cyj2, czk1, czk2
      integer :: pos

      associate (b=>domain_bc_send_buffers(num))

        b%comm = world_comm
        b%remote_rank = domain_ranks_grid(domain)%arr(im(1),im(2),im(3))
        b%remote_domain = domain
        
        b%direction = dir
        
        b%spatial_ratio = domain_grids(domain)%spatial_ratio

        b%enabled = .true.

        !get the child image extent
        cxmin = domain_grids(domain)%xmins(im(1))
        cxmax = domain_grids(domain)%xmaxs(im(1))
        cymin = domain_grids(domain)%ymins(im(2))
        cymax = domain_grids(domain)%ymaxs(im(2))
        czmin = domain_grids(domain)%zmins(im(3))
        czmax = domain_grids(domain)%zmaxs(im(3))

        !find the child grid indexes
        cxi1 = x_coord(cxmin)
        cxi2 = x_coord(cxmax)
        cyj1 = y_coord(cymin)
        cyj2 = y_coord(cymax)
        czk1 = z_coord(czmin)
        czk2 = z_coord(czmax)

        select case (dir)
          case (We)
            pos = cxi1
          case (Ea)
            pos = cxi2
          case (So)
            pos = cyj1
          case (No)
            pos = cyj2
          case (Bo)
            pos = czk1
          case (To)
            pos = czk2
        end select

        if (dir==We .or. dir==Ea) then

          if (dir==We) then
            Ui1 = pos-2
            Ui2 = pos+2

            Vi1 = pos-2
            Vi2 = pos+3

            Wi1 = pos-2
            Wi2 = pos+3

            Pri1 = pos-1
            Pri2 = pos+3
          else if (dir==Ea) then
            Ui1 = pos-2
            Ui2 = pos+2

            Vi1 = pos-2
            Vi2 = pos+3

            Wi1 = pos-2
            Wi2 = pos+3

            Pri1 = pos-2
            Pri2 = pos+2
          end if

          Uj1 = cyj1-2
          Uj2 = cyj2+3
          Uk1 = czk1-2
          Uk2 = czk2+3

          Vj1 = cyj1-2
          if (domain_grids(domain)%internal_or_periodic_No(im(2))) then
            Vj2 = cyj2+3
          else
            Vj2 = cyj2+2
          end if
          Vk1 = czk1-2
          Vk2 = czk2+3

          Wj1 = cyj1-2
          Wj2 = cyj2+3
          Wk1 = czk1-2
          if (domain_grids(domain)%internal_or_periodic_To(im(3))) then
            Wk2 = czk2+3
          else
            Wk2 = czk2+2
          end if

          Prj1 = cyj1-1
          Prj2 = cyj2+2
          Prk1 = czk1-1
          Prk2 = czk2+2

        else if (dir==So .or. dir==No) then

          if (dir==So) then
            Uj1 = pos-2
            Uj2 = pos+3

            Vj1 = pos-2
            Vj2 = pos+2

            Wj1 = pos-2
            Wj2 = pos+3

            Prj1 = pos-1
            Prj2 = pos+3
          else if (dir==No) then
            Uj1 = pos-2
            Uj2 = pos+3

            Vj1 = pos-2
            Vj2 = pos+2

            Wj1 = pos-2
            Wj2 = pos+3

            Prj1 = pos-2
            Prj2 = pos+2
          end if

          Ui1 = cxi1-2
          if (domain_grids(domain)%internal_or_periodic_Ea(im(1))) then
            Ui2 = cxi2+3
          else
            Ui2 = cxi2+2
          end if
          Uk1 = czk1-2
          Uk2 = czk2+3

          Vi1 = cxi1-2
          Vi2 = cxi2+3
          Vk1 = czk1-2
          Vk2 = czk2+3

          Wi1 = cxi1-2
          Wi2 = cxi2+3
          Wk1 = czk1-2
          if (domain_grids(domain)%internal_or_periodic_To(im(3))) then
            Wk2 = czk2+3
          else
            Wk2 = czk2+2
          end if

          Pri1 = cxi1-1
          Pri2 = cxi2+2
          Prk1 = czk1-1
          Prk2 = czk2+2

        else

          if (dir==Bo) then
            Uk1 = pos-2
            Uk2 = pos+3

            Vk1 = pos-2
            Vk2 = pos+3

            Wk1 = pos-2
            Wk2 = pos+2

            Prk1 = pos-1
            Prk2 = pos+3
          else if (dir==To) then
            Uk1 = pos-2
            Uk2 = pos+3

            Vk1 = pos-2
            Vk2 = pos+3

            Wk1 = pos-2
            Wk2 = pos+2

            Prk1 = pos-2
            Prk2 = pos+2
          end if

          Ui1 = cxi1-2
          if (domain_grids(domain)%internal_or_periodic_Ea(im(1))) then
            Ui2 = cxi2+3
          else
            Ui2 = cxi2+2
          end if
          Uj1 = cyj1-2
          Uj2 = cyj2+3

          Vi1 = cxi1-2
          Vi2 = cxi2+3
          Vj1 = cyj1-2
          if (domain_grids(domain)%internal_or_periodic_No(im(2))) then
            Vj2 = cyj2+3
          else
            Vj2 = cyj2+2
          end if

          Wi1 = cxi1-2
          Wi2 = cxi2+3
          Wj1 = cyj1-2
          Wj2 = cyj2+3

          Pri1 = cxi1-1
          Pri2 = cxi2+2
          Prj1 = cyj1-1
          Prj2 = cyj2+2

        end if

        b%position = pos

        b%Ui1 = Ui1 
        b%Ui2 = Ui2 
        b%Uj1 = Uj1 
        b%Uj2 = Uj2 
        b%Uk1 = Uk1 
        b%Uk2 = Uk2 

        b%Vi1 = Vi1 
        b%Vi2 = Vi2 
        b%Vj1 = Vj1 
        b%Vj2 = Vj2 
        b%Vk1 = Vk1 
        b%Vk2 = Vk2 

        b%Wi1 = Wi1 
        b%Wi2 = Wi2 
        b%Wj1 = Wj1 
        b%Wj2 = Wj2 
        b%Wk1 = Wk1 
        b%Wk2 = Wk2 

        b%Pri1 = Pri1
        b%Pri2 = Pri2
        b%Prj1 = Prj1
        b%Prj2 = Prj2
        b%Prk1 = Prk1
        b%Prk2 = Prk2

        allocate(b%U(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
        allocate(b%V(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
        allocate(b%W(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

        allocate(b%dU_dt(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
        allocate(b%dV_dt(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
        allocate(b%dW_dt(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

        if (enable_buoyancy) then
          allocate(b%Temperature(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          allocate(b%dTemperature_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
        end if
        if (enable_moisture) then
          allocate(b%Moisture(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          allocate(b%dMoisture_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
        end if

        if (num_of_scalars > 0) then
          allocate(b%Scalar(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2,num_of_scalars))
          allocate(b%dScalar_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2,num_of_scalars))
        end if

      end associate

    end subroutine create_boundary_parent_buffer


    subroutine create_boundary_child_buffer(dir, domain)
      use fftw3
      integer :: Ui1, Ui2, Vi1, Vi2, Wi1, Wi2, Pri1, Pri2
      integer :: Uj1, Uj2, Vj1, Vj2, Wj1, Wj2, Prj1, Prj2
      integer :: Uk1, Uk2, Vk1, Vk2, Wk1, Wk2, Prk1, Prk2
      integer, intent(in) :: dir, domain
      integer :: i, width
      integer :: nx, ny, nz

      associate(b => domain_bc_recv_buffers(dir))

        b%comm = world_comm
        b%remote_rank = domain_ranks_grid(domain)%arr(parent_image(1), &
                                                      parent_image(2), &
                                                      parent_image(3))       
        b%remote_domain = domain
        b%direction = dir
        b%enabled = .true.

        b%relaxation = has_domain_boundary_relaxation(dir)
        b%relaxation_width = domain_boundary_relaxation_width(dir)
        b%relax_factor = domain_boundary_relaxation_factor(dir)

!spectral interpolation is disabled by default. Big problems with mean profile distortion.
!         if (b%direction==We) b%interp_order = 0

        b%spatial_ratio = domain_spatial_ratio
        b%time_step_ratio = domain_time_step_ratio

        width = b%spatial_ratio * b%relaxation_width

        select case  (dir)
          case (We)
            b%position = 0

            Ui1 = -2
            Ui2 = width
            Uj1 = -2
            Uj2 = Uny+3
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = -2
            Vi2 = width
            Vj1 = -2
            Vj2 = Vny+3
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = -2
            Wi2 = width
            Wj1 = -2
            Wj2 = Wny+3
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = -1
            Pri2 = width
            Prj1 = -1
            Prj2 = Prny+2
            Prk1 = -1
            Prk2 = Prnz+2

            b%r_Ui1 = -2
            b%r_Ui2 = +2
            b%r_Uj1 = -2
            b%r_Uj2 = Uny/b%spatial_ratio+3
            b%r_Uk1 = -2
            b%r_Uk2 = Unz/b%spatial_ratio+3
            
            b%r_Vi1 = -2
            b%r_Vi2 = +3
            b%r_Vj1 = -2
            b%r_Vj2 = Vny/b%spatial_ratio+3
            b%r_Vk1 = -2
            b%r_Vk2 = Vnz/b%spatial_ratio+3
            
            b%r_Wi1 = -2
            b%r_Wi2 = +3
            b%r_Wj1 = -2
            b%r_Wj2 = Wny/b%spatial_ratio+3
            b%r_Wk1 = -2
            b%r_Wk2 = Wnz/b%spatial_ratio+3
            
            b%r_Pri1 = -1
            b%r_Pri2 = +3
            b%r_Prj1 = -1
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = -1
            b%r_Prk2 = Prnz/b%spatial_ratio+2

          case (Ea)
            b%position = Prnx

            Ui1 = Unx-width+1
            Ui2 = Prnx+2
            Uj1 = -2
            Uj2 = Uny+3
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = Vnx-width+1
            Vi2 = Prnx+3
            Vj1 = -2
            Vj2 = Vny+3
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = Wnx-width+1
            Wi2 = Prnx+3
            Wj1 = -2
            Wj2 = Wny+3
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = Prnx-width+1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = Prny+2
            Prk1 = -1
            Prk2 = Prnz+2

            b%r_Ui1 = Prnx/b%spatial_ratio-2
            b%r_Ui2 = Prnx/b%spatial_ratio+2
            b%r_Uj1 = -2
            b%r_Uj2 = Uny/b%spatial_ratio+3
            b%r_Uk1 = -2
            b%r_Uk2 = Unz/b%spatial_ratio+3
            
            b%r_Vi1 = Prnx/b%spatial_ratio-2
            b%r_Vi2 = Prnx/b%spatial_ratio+3
            b%r_Vj1 = -2
            b%r_Vj2 = Vny/b%spatial_ratio+3
            b%r_Vk1 = -2
            b%r_Vk2 = Vnz/b%spatial_ratio+3
            
            b%r_Wi1 = Prnx/b%spatial_ratio-2
            b%r_Wi2 = Prnx/b%spatial_ratio+3
            b%r_Wj1 = -2
            b%r_Wj2 = Wny/b%spatial_ratio+3
            b%r_Wk1 = -2
            b%r_Wk2 = Wnz/b%spatial_ratio+3
            
            b%r_Pri1 = Prnx/b%spatial_ratio-2
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = -1
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = -1
            b%r_Prk2 = Prnz/b%spatial_ratio+2

          case (So)
            b%position = 0

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = -2
            Uj2 = width
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = -2
            Vj2 = width
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = -2
            Wj2 = width
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = width
            Prk1 = -1
            Prk2 = Prnz+2

            b%r_Ui1 = -2
            b%r_Ui2 = Unx/b%spatial_ratio+3
            b%r_Uj1 = -2
            b%r_Uj2 = +3
            b%r_Uk1 = -2
            b%r_Uk2 = Unz/b%spatial_ratio+3
            
            b%r_Vi1 = -2
            b%r_Vi2 = Vnx/b%spatial_ratio+3
            b%r_Vj1 = -2
            b%r_Vj2 = +2
            b%r_Vk1 = -2
            b%r_Vk2 = Vnz/b%spatial_ratio+3
            
            b%r_Wi1 = -2
            b%r_Wi2 = Wnx/b%spatial_ratio+3
            b%r_Wj1 = -2
            b%r_Wj2 = +3
            b%r_Wk1 = -2
            b%r_Wk2 = Wnz/b%spatial_ratio+3
            
            b%r_Pri1 = -1
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = -1
            b%r_Prj2 = +3
            b%r_Prk1 = -1
            b%r_Prk2 = Prnz/b%spatial_ratio+2

          case (No)
            b%position = Prny

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = Uny-width+1
            Uj2 = Prny+3
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = Vny-width+1
            Vj2 = Prny+2
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = Wny-width+1
            Wj2 = Prny+3
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = Prny-width+1
            Prj2 = Prny+2
            Prk1 = -1
            Prk2 = Prnz+2

            b%r_Ui1 = -2
            b%r_Ui2 = Unx/b%spatial_ratio+3
            b%r_Uj1 = Prny/b%spatial_ratio-2
            b%r_Uj2 = Prny/b%spatial_ratio+3
            b%r_Uk1 = -2
            b%r_Uk2 = Unz/b%spatial_ratio+3
            
            b%r_Vi1 = -2
            b%r_Vi2 = Vnx/b%spatial_ratio+3
            b%r_Vj1 = Prny/b%spatial_ratio-2
            b%r_Vj2 = Prny/b%spatial_ratio+2
            b%r_Vk1 = -2
            b%r_Vk2 = Vnz/b%spatial_ratio+3
            
            b%r_Wi1 = -2
            b%r_Wi2 = Wnx/b%spatial_ratio+3
            b%r_Wj1 = Prny/b%spatial_ratio-2
            b%r_Wj2 = Prny/b%spatial_ratio+3
            b%r_Wk1 = -2
            b%r_Wk2 = Wnz/b%spatial_ratio+3
            
            b%r_Pri1 = -1
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = Prny/b%spatial_ratio-2
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = -1
            b%r_Prk2 = Prnz/b%spatial_ratio+2

          case (Bo)
            b%position = 0

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = -2
            Uj2 = Uny+3
            Uk1 = -2
            Uk2 = width

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = -2
            Vj2 = Vny+3
            Vk1 = -2
            Vk2 = width

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = -2
            Wj2 = Wny+3
            Wk1 = -2
            Wk2 = width

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = Prny+2
            Prk1 = -1
            Prk2 = width

            b%r_Ui1 = -2
            b%r_Ui2 = Unx/b%spatial_ratio+3
            b%r_Uj1 = -2
            b%r_Uj2 = Uny/b%spatial_ratio+3
            b%r_Uk1 = -2
            b%r_Uk2 = +3
            
            b%r_Vi1 = -2
            b%r_Vi2 = Vnx/b%spatial_ratio+3
            b%r_Vj1 = -2
            b%r_Vj2 = Vny/b%spatial_ratio+3
            b%r_Vk1 = -2
            b%r_Vk2 = +3
            
            b%r_Wi1 = -2
            b%r_Wi2 = Wnx/b%spatial_ratio+3
            b%r_Wj1 = -2
            b%r_Wj2 = Wny/b%spatial_ratio+3
            b%r_Wk1 = -2
            b%r_Wk2 = +2
            
            b%r_Pri1 = -1
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = -1
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = -1
            b%r_Prk2 = +3

          case (To)
            b%position = Prnz

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = -2
            Uj2 = Uny+3
            Uk1 = Unz-width+1
            Uk2 = Prnz+3

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = -2
            Vj2 = Vny+3
            Vk1 = Vnz-width+1
            Vk2 = Prnz+3

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = -2
            Wj2 = Wny+3
            Wk1 = Wnz-width+1
            Wk2 = Prnz+2

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = Prny+2
            Prk1 = Prnz-width+1
            Prk2 = Prnz+2

            b%r_Ui1 = -2
            b%r_Ui2 = Unx/b%spatial_ratio+3
            b%r_Uj1 = -2
            b%r_Uj2 = Uny/b%spatial_ratio+3
            b%r_Uk1 = Prnz/b%spatial_ratio-2
            b%r_Uk2 = Prnz/b%spatial_ratio+3
            
            b%r_Vi1 = -2
            b%r_Vi2 = Vnx/b%spatial_ratio+3
            b%r_Vj1 = -2
            b%r_Vj2 = Vny/b%spatial_ratio+3
            b%r_Vk1 = Prnz/b%spatial_ratio-2
            b%r_Vk2 = Prnz/b%spatial_ratio+3
            
            b%r_Wi1 = -2
            b%r_Wi2 = Wnx/b%spatial_ratio+3
            b%r_Wj1 = -2
            b%r_Wj2 = Wny/b%spatial_ratio+3
            b%r_Wk1 = Prnz/b%spatial_ratio-2
            b%r_Wk2 = Prnz/b%spatial_ratio+2
            
            b%r_Pri1 = -1
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = -1
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = Prnz/b%spatial_ratio-2
            b%r_Prk2 = Prnz/b%spatial_ratio+2

        end select    

        b%Ui1 = Ui1 
        b%Ui2 = Ui2 
        b%Uj1 = Uj1 
        b%Uj2 = Uj2 
        b%Uk1 = Uk1 
        b%Uk2 = Uk2 

        b%Vi1 = Vi1 
        b%Vi2 = Vi2 
        b%Vj1 = Vj1 
        b%Vj2 = Vj2 
        b%Vk1 = Vk1 
        b%Vk2 = Vk2 

        b%Wi1 = Wi1 
        b%Wi2 = Wi2 
        b%Wj1 = Wj1 
        b%Wj2 = Wj2 
        b%Wk1 = Wk1 
        b%Wk2 = Wk2 

        b%Pri1 = Pri1
        b%Pri2 = Pri2
        b%Prj1 = Prj1
        b%Prj2 = Prj2
        b%Prk1 = Prk1
        b%Prk2 = Prk2


        b%bUi1 = Ui1 
        b%bUi2 = Ui2 
        b%bUj1 = Uj1 
        b%bUj2 = Uj2 
        b%bUk1 = Uk1 
        b%bUk2 = Uk2 

        b%bVi1 = Vi1 
        b%bVi2 = Vi2 
        b%bVj1 = Vj1 
        b%bVj2 = Vj2 
        b%bVk1 = Vk1 
        b%bVk2 = Vk2 

        b%bWi1 = Wi1 
        b%bWi2 = Wi2 
        b%bWj1 = Wj1 
        b%bWj2 = Wj2 
        b%bWk1 = Wk1 
        b%bWk2 = Wk2 

        b%bPri1 = Pri1
        b%bPri2 = Pri2
        b%bPrj1 = Prj1
        b%bPrj2 = Prj2
        b%bPrk1 = Prk1
        b%bPrk2 = Prk2

        select case  (dir)
          case (We)
            b%bUi2 = Ui2 - width 
            b%bVi2 = Vi2 - width 
            b%bWi2 = Wi2 - width 
            b%bPri2 = Pri2 - width 
          case (Ea)
            b%bUi1 = Ui1 + width 
            b%bVi1 = Vi1 + width 
            b%bWi1 = Wi1 + width 
            b%bPri1 = Pri1 + width 
          case (So)
            b%bUj2 = Uj2 - width 
            b%bVj2 = Vj2 - width 
            b%bWj2 = Wj2 - width 
            b%bPrj2 = Prj2 - width 
          case (No)
            b%bUj1 = Uj1 + width 
            b%bVj1 = Vj1 + width 
            b%bWj1 = Wj1 + width 
            b%bPrj1 = Prj1 + width 
          case (Bo)
            b%bUk2 = Uk2 - width 
            b%bVk2 = Vk2 - width 
            b%bWk2 = Wk2 - width 
            b%bPrk2 = Prk2 - width 
          case (To)
            b%bUk1 = Uk1 + width 
            b%bVk1 = Vk1 + width 
            b%bWk1 = Wk1 + width 
            b%bPrk1 = Prk1 + width 
        end select

         
        b%r_x0 = im_xmin
        b%r_y0 = im_ymin
        b%r_z0 = im_zmin
        b%r_dx = domain_grids(domain)%dx
        b%r_dy = domain_grids(domain)%dy
        b%r_dz = domain_grids(domain)%dz
        
        allocate(b%r_xU(b%r_Ui1:b%r_Ui2))
        allocate(b%r_yV(b%r_Vj1:b%r_Vj2))
        allocate(b%r_zW(b%r_Wk1:b%r_Wk2))
        allocate(b%r_x(b%r_Vi1:b%r_Vi2))
        allocate(b%r_y(b%r_Uj1:b%r_Uj2))
        allocate(b%r_z(b%r_Uk1:b%r_Uk2))
        b%r_xU = [(b%r_x0 + i*b%r_dx , i = b%r_Ui1,b%r_Ui2)]
        b%r_yV = [(b%r_y0 + i*b%r_dy , i = b%r_Vj1,b%r_Vj2)]
        b%r_zW = [(b%r_z0 + i*b%r_dz , i = b%r_Wk1,b%r_Wk2)]
        b%r_x = [(b%r_x0 + (i-0.5_knd)*b%r_dx , i = b%r_Vi1,b%r_Vi2)]
        b%r_y = [(b%r_y0 + (i-0.5_knd)*b%r_dy , i = b%r_Uj1,b%r_Uj2)]
        b%r_z = [(b%r_z0 + (i-0.5_knd)*b%r_dz , i = b%r_Vk1,b%r_Vk2)]


        allocate(b%U(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
        allocate(b%V(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
        allocate(b%W(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

        allocate(b%dU_dt(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
        allocate(b%dV_dt(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
        allocate(b%dW_dt(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

        if (enable_buoyancy) then
          allocate(b%Temperature(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          allocate(b%dTemperature_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
        end if
        if (enable_moisture) then
          allocate(b%Moisture(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          allocate(b%dMoisture_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
        end if

        if (num_of_scalars > 0) then
          allocate(b%Scalar(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2,num_of_scalars))
          allocate(b%dScalar_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2,num_of_scalars))
        end if

        allocate(b%r_U(b%r_Ui1:b%r_Ui2,b%r_Uj1:b%r_Uj2,b%r_Uk1:b%r_Uk2))
        allocate(b%r_V(b%r_Vi1:b%r_Vi2,b%r_Vj1:b%r_Vj2,b%r_Vk1:b%r_Vk2))
        allocate(b%r_W(b%r_Wi1:b%r_Wi2,b%r_Wj1:b%r_Wj2,b%r_Wk1:b%r_Wk2))

        allocate(b%r_dU_dt(b%r_Ui1:b%r_Ui2,b%r_Uj1:b%r_Uj2,b%r_Uk1:b%r_Uk2))
        allocate(b%r_dV_dt(b%r_Vi1:b%r_Vi2,b%r_Vj1:b%r_Vj2,b%r_Vk1:b%r_Vk2))
        allocate(b%r_dW_dt(b%r_Wi1:b%r_Wi2,b%r_Wj1:b%r_Wj2,b%r_Wk1:b%r_Wk2))

        if (enable_buoyancy) then
          allocate(b%r_Temperature(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2))
          allocate(b%r_dTemperature_dt(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2))
        end if
        if (enable_moisture) then
          allocate(b%r_Moisture(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2))
          allocate(b%r_dMoisture_dt(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2))
        end if

        if (num_of_scalars > 0) then
          allocate(b%r_Scalar(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2,num_of_scalars))
          allocate(b%r_dScalar_dt(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2,num_of_scalars))
        end if


        if (any(b%interp_order==0)) then !.and.b%direction==We) then
          if (b%direction==We) b%relaxation = .false.
          allocate(b%interpolation%fft_r_U(size(b%r_U,1)/2+1,size(b%r_U,2),size(b%r_U,3)))
          allocate(b%interpolation%fft_r_V(size(b%r_V,1)/2+1,size(b%r_V,2),size(b%r_V,3)))
          allocate(b%interpolation%fft_r_W(size(b%r_W,1)/2+1,size(b%r_W,2),size(b%r_W,3)))
          
!           call pfft_plan_dft_r2c(int([size(b%r_U,3),size(b%r_U,2),size(b%r_U,1)],c_intptr_t), &
!                                  b%r_U, b%interpolation%fft_r_U, 
          b%interpolation%forward_U = fftw_plan_gen(size(b%r_U,3), size(b%r_U,2), size(b%r_U,1), &
                b%r_U, b%interpolation%fft_r_U, FFTW_UNALIGNED)
          b%interpolation%forward_V = fftw_plan_gen(size(b%r_V,3), size(b%r_V,2), size(b%r_V,1), &
                b%r_V, b%interpolation%fft_r_V, FFTW_UNALIGNED)
          b%interpolation%forward_W = fftw_plan_gen(size(b%r_W,3), size(b%r_W,2), size(b%r_W,1), &
                b%r_W, b%interpolation%fft_r_W, FFTW_UNALIGNED)

          nx = size(b%r_U,1)*b%spatial_ratio
          ny = size(b%r_U,2)*b%spatial_ratio
          nz = size(b%r_U,3)*b%spatial_ratio
          allocate(b%interpolation%fft_U(nx/2+1,ny,nz))
          allocate(b%interpolation%trans_U(nx,ny,nz))
          b%interpolation%backward_U = fftw_plan_gen(nz, ny, nx, &
                b%interpolation%fft_U, b%interpolation%trans_U, FFTW_UNALIGNED)

          nx = size(b%r_V,1)*b%spatial_ratio
          ny = size(b%r_V,2)*b%spatial_ratio
          nz = size(b%r_V,3)*b%spatial_ratio
          allocate(b%interpolation%fft_V(nx/2+1,ny,nz))
          allocate(b%interpolation%trans_V(nx,ny,nz))
          b%interpolation%backward_V = fftw_plan_gen(nz, ny, nx, &
                b%interpolation%fft_V, b%interpolation%trans_V, FFTW_UNALIGNED)

          nx = size(b%r_W,1)*b%spatial_ratio
          ny = size(b%r_W,2)*b%spatial_ratio
          nz = size(b%r_W,3)*b%spatial_ratio
          allocate(b%interpolation%fft_W(nx/2+1,ny,nz))
          allocate(b%interpolation%trans_W(nx,ny,nz))
          b%interpolation%backward_W = fftw_plan_gen(nz, ny, nx, &
                b%interpolation%fft_W, b%interpolation%trans_W, FFTW_UNALIGNED)

        end if



        if (has_domain_boundary_turbulence_generator(dir)) then
          b%turb_generator_enabled = .true.
          b%relaxation = .false.

          allocate(b%turb_generator)

          select case (dir)
            case (We:Ea)
              allocate(b%U_turb(Uj1:Uj2,Uk1:Uk2))
              allocate(b%V_turb(Vj1:Vj2,Vk1:Vk2))
              allocate(b%W_turb(Wj1:Wj2,Wk1:Wk2))
              allocate(b%turb_generator%sgs_tke(1:Prny, 1:Prnz))
            case (So:No)
              allocate(b%U_turb(Ui1:Ui2,Uk1:Uk2))
              allocate(b%V_turb(Vi1:Vi2,Vk1:Vk2))
              allocate(b%W_turb(Wi1:Wi2,Wk1:Wk2))
              allocate(b%turb_generator%sgs_tke(1:Prnx, 1:Prnz))
            case (Bo:To)
              allocate(b%U_turb(Ui1:Ui2,Uj1:Uj2))
              allocate(b%V_turb(Vi1:Vi2,Vj1:Vj2))
              allocate(b%W_turb(Wi1:Wi2,Wj1:Wj2))
              allocate(b%turb_generator%sgs_tke(1:Prnx, 1:Prny))
          end select

          b%U_turb = 0
          b%V_turb = 0
          b%W_turb = 0

          b%turb_generator%sgs_tke = 0.01
          b%turb_generator%L_y = dymin * b%spatial_ratio / 2
          b%turb_generator%L_z = dzmin * b%spatial_ratio / 2
          b%turb_generator%T_lag = time_stepping%dt_constant * b%time_step_ratio * 2
          call b%turb_generator%init()
        end if

      end associate

    end subroutine create_boundary_child_buffer



    subroutine create_domain_parent_buffer(b, child_domain, im)
      use Strings, only: itoa
      type(dom_parent_buffer), intent(out) :: b
      integer, intent(in) :: child_domain
      integer, intent(in) :: im(3)
      real(knd) :: cxmax, cxmin, cymax, cymin, czmax, czmin
      integer :: cxi1, cxi2, cyj1, cyj2, czk1, czk2

      b%comm = world_comm
      b%remote_image = im
      b%remote_rank = domain_ranks_grid(child_domain)%arr(im(1),im(2),im(3))

      b%spatial_ratio = domain_grids(child_domain)%spatial_ratio

      b%is_two_way_nested = domain_is_domain_two_way_nested(child_domain)
      
      !Get the child image extent.
      cxmin = domain_grids(child_domain)%xmins(im(1))
      cxmax = domain_grids(child_domain)%xmaxs(im(1))
      cymin = domain_grids(child_domain)%ymins(im(2))
      cymax = domain_grids(child_domain)%ymaxs(im(2))
      czmin = domain_grids(child_domain)%zmins(im(3))
      czmax = domain_grids(child_domain)%zmaxs(im(3))

      !Find the child grid indexes.
      cxi1 = x_coord(cxmin)
      cxi2 = x_coord(cxmax)
      cyj1 = y_coord(cymin)
      cyj2 = y_coord(cymax)
      czk1 = z_coord(czmin)
      czk2 = z_coord(czmax)
      
      !Check that grid cell dimensions correspond to the specified spatial refinement ratio.
      if (abs(dxmin/domain_grids(child_domain)%dx - b%spatial_ratio) > 100*epsilon(1.0_real32)) then
        write(*,*) "Error, inconsistent specification of the grid refinement ratio and cell sizes of domain ",child_domain,"."
        write(*,*) "spatial_ratio specified:",b%spatial_ratio
        write(*,*) "parent dx:",dxmin
        write(*,*) "child dx:",domain_grids(child_domain)%dx
        write(*,*) "ratio:",dxmin/domain_grids(child_domain)%dx
        call error_stop()
      end if
      if (abs(dymin/domain_grids(child_domain)%dy - b%spatial_ratio) > 100*epsilon(1.0_real32)) then
        write(*,*) "Error, inconsistent specification of the grid refinement ratio and cell sizes of domain ",child_domain,"."
        write(*,*) "spatial_ratio specified:",b%spatial_ratio
        write(*,*) "parent dy:",dymin
        write(*,*) "child dy:",domain_grids(child_domain)%dy
        call error_stop()
      end if
      if (abs(dzmin/domain_grids(child_domain)%dz - b%spatial_ratio) > 100*epsilon(1.0_real32)) then
        write(*,*) "Error, inconsistent specification of the grid refinement ratio and cell sizes of domain ",child_domain,"."
        write(*,*) "spatial_ratio specified:",b%spatial_ratio
        write(*,*) "parent dz:",dzmin
        write(*,*) "child dz:",domain_grids(child_domain)%dz
        call error_stop()
      end if
      
      
      !Check that child image does not exceed this (parent) image.
      !It may exceed it if no neighbouring boundary is nested and two-way nesting is disabled.
      if ((cxmin < im_xmin - dxmin/100) .and. &
          (b%is_two_way_nested  .or. &
           any(domain_is_boundary_nested([We,So,No,Bo,To],child_domain)))) then
        write(*,*) "Child image",im,"has lower x bound", cxmin, &
                   " which is lower than the lower x bound", im_xmin , &
                   "of parent image",iim,jim,kim
        call error_stop()
      end if
      if ((cxmax > im_xmax + dxmin/100) .and. &
          (b%is_two_way_nested  .or. &
           any(domain_is_boundary_nested([Ea,So,No,Bo,To],child_domain)))) then
        write(*,*) "Child image",im,"has upper x bound", cxmax, &
                   " which is higher than the upper x bound", im_xmax , &
                   "of parent image",iim,jim,kim
        call error_stop()
      end if
      if ((cymin < im_ymin - dymin/100) .and. &
          (b%is_two_way_nested  .or. &
           any(domain_is_boundary_nested([We,Ea,So,Bo,To],child_domain)))) then
        write(*,*) "Child image",im,"has lower y bound", cymin, &
                   " which is lower than the lower y bound", im_ymin , &
                   "of parent image",iim,jim,kim
        call error_stop()
      end if
      if ((cymax > im_ymax + dymin/100) .and. &
          (b%is_two_way_nested  .or. &
           any(domain_is_boundary_nested([We,Ea,No,Bo,To],child_domain)))) then
        write(*,*) "Child image",im,"has upper y bound", cymax, &
                   " which is higher than the upper y bound", im_ymax , &
                   "of parent image",iim,jim,kim
        call error_stop()
      end if
      if ((czmin < im_zmin - dzmin/100) .and. &
          (b%is_two_way_nested  .or. &
           any(domain_is_boundary_nested([We,Ea,So,No,Bo],child_domain)))) then
        write(*,*) "Child image",im,"has lower z bound", czmin, &
                   " which is lower than the lower z bound", im_zmin , &
                   "of parent image",iim,jim,kim
        call error_stop()
      end if
      if ((czmax > im_zmax + dzmin/100) .and. &
          (b%is_two_way_nested  .or. &
           any(domain_is_boundary_nested([We,Ea,So,No,To],child_domain)))) then
        write(*,*) "Child image",im,"has upper z bound", czmax, &
                   " which is higher than the upper z bound", im_zmax , &
                   "of parent image",iim,jim,kim
        call error_stop()
      end if
      
      
      
      !Check that the grids are aligned.
      !Disregard if the boundary is outside of the parent domain (checked above).
      if ((cxmin>im_xmin-dxmin/10 .and. cxmin<im_xmax+dxmin/10) .and. &
          .not.(abs(cxi1*dxmin + im_xmin - cxmin) < dxmin/100)) &
        call error_stop("The west boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned &
                       &or the image boundary divides a parent-domain cell")
      if ((cxmax>im_xmin-dxmin/10 .and. cxmax<im_xmax+dxmin/10) .and. &
          .not.(abs(cxi2*dxmin + im_xmin - cxmax) < dxmin/100)) &
        call error_stop("The east boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned &
                       &or the image boundary divides a parent-domain cell")
      if ((cymin>im_ymin-dymin/10 .and. cymin<im_ymax+dymin/10) .and. &
          .not.(abs(cyj1*dymin + im_ymin - cymin) < dymin/100)) &
        call error_stop("The south boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned &
                       &or the image boundary divides a parent-domain cell")
      if ((cymax>im_ymin-dymin/10 .and. cymax<im_ymax+dymin/10) .and. &
          .not.(abs(cyj2*dymin + im_ymin - cymax) < dymin/100)) &
        call error_stop("The north boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned &
                       &or the image boundary divides a parent-domain cell")
      if ((czmin>im_zmin-dzmin/10 .and. czmin<im_zmax+dzmin/10) .and. &
          .not.(abs(czk1*dzmin + im_zmin - czmin) < dzmin/100)) &
        call error_stop("The bottom boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned &
                       &or the image boundary divides a parent-domain cell")
      if ((czmax>im_zmin-dzmin/10 .and. czmax<im_zmax+dzmin/10) .and. &
          .not.(abs(czk2*dzmin + im_zmin - czmax) < dzmin/100)) &
        call error_stop("The top boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned &
                       &or the image boundary divides a parent-domain cell")


      

      !corresponds to index 0 on the child grid
      b%i1 = cxi1
      !we need index Prnx+1 on the child grid
      b%i2 = cxi2 + 1

      b%j1 = cyj1
      b%j2 = cyj2 + 1

      b%k1 = czk1
      b%k2 = czk2 + 1




      b%exchange_pr_gradient_x = enable_fixed_flow_rate .and. flow_rate_x_fixed
      b%exchange_pr_gradient_y = enable_fixed_flow_rate .and. flow_rate_y_fixed
      b%exchange_pr_gradient_z = enable_fixed_flow_rate .and. flow_rate_z_fixed

      call MPI_ISend(b%exchange_pr_gradient_x, 1, MPI_LOGICAL, &
                     b%remote_rank, 1, b%comm, &
                     request, err)
      requests = [requests, request]

      call MPI_ISend(b%exchange_pr_gradient_y, 1, MPI_LOGICAL, &
                     b%remote_rank, 2, b%comm, &
                     request, err)
      requests = [requests, request]

      call MPI_ISend(b%exchange_pr_gradient_z, 1, MPI_LOGICAL, &
                     b%remote_rank, 2, b%comm, &
                     request, err)
      requests = [requests, request]
    end subroutine



    subroutine create_domain_child_buffer(b)
      type(dom_child_buffer), intent(out) :: b
      integer :: i, j, k

      b%comm = world_comm
      b%remote_image = parent_image
      b%remote_rank = domain_ranks_grid(parent_domain)%arr(parent_image(1), &
                                                           parent_image(2), &
                                                           parent_image(3))
      b%spatial_ratio = domain_spatial_ratio
      b%time_step_ratio = domain_time_step_ratio

      b%is_two_way_nested = is_this_domain_two_way_nested

      b%r_i1 = 0
      b%r_i2 = Prnx/b%spatial_ratio + 1
      b%r_j1 = 0
      b%r_j2 = Prny/b%spatial_ratio + 1
      b%r_k1 = 0
      b%r_k2 = Prnz/b%spatial_ratio + 1

      b%r_x0 = im_xmin
      b%r_y0 = im_ymin
      b%r_z0 = im_zmin
      b%r_dx = domain_grids(parent_domain)%dx
      b%r_dy = domain_grids(parent_domain)%dy
      b%r_dz = domain_grids(parent_domain)%dz
      
      allocate(b%r_xU(b%r_i1:b%r_i2))
      allocate(b%r_yV(b%r_j1:b%r_j2))
      allocate(b%r_zW(b%r_k1:b%r_k2))
      allocate(b%r_x(b%r_i1:b%r_i2))
      allocate(b%r_y(b%r_j1:b%r_j2))
      allocate(b%r_z(b%r_k1:b%r_k2))
      b%r_xU(:) = [(b%r_x0 + i*b%r_dx , i = b%r_i1,b%r_i2)]
      b%r_yV(:) = [(b%r_y0 + i*b%r_dy , i = b%r_j1,b%r_j2)]
      b%r_zW(:) = [(b%r_z0 + i*b%r_dz , i = b%r_k1,b%r_k2)]
      b%r_x(:) = [(b%r_x0 + (i-0.5_knd)*b%r_dx , i = b%r_i1,b%r_i2)]
      b%r_y(:) = [(b%r_y0 + (i-0.5_knd)*b%r_dy , i = b%r_j1,b%r_j2)]
      b%r_z(:) = [(b%r_z0 + (i-0.5_knd)*b%r_dz , i = b%r_k1,b%r_k2)]

      call MPI_IRecv(b%exchange_pr_gradient_x, 1, MPI_LOGICAL, &
                     b%remote_rank, 1, b%comm, &
                     request, err)
      requests = [requests, request]

      call MPI_IRecv(b%exchange_pr_gradient_y, 1, MPI_LOGICAL, &
                     b%remote_rank, 2, b%comm, &
                     request, err)
      requests = [requests, request]

      call MPI_IRecv(b%exchange_pr_gradient_z, 1, MPI_LOGICAL, &
                     b%remote_rank, 2, b%comm, &
                     request, err)
      requests = [requests, request]
    end subroutine

  end subroutine

end module domains_bc_init
