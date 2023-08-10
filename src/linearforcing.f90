module LinearForcing
  use Parameters
#ifdef PAR
  use custom_par
#endif

  implicit none
  
  private
  
  public init_linear_forcing, linear_forcing
  
  logical :: enable_linear_forcing = .false.
  logical :: enable_linear_forcing_tke0_epsilon0 = .false.
  logical :: enable_forcing_variable_means = .false.
  logical :: enable_filter = .false.
  
  real(knd) :: linear_forcing_constant = 0
  
  real(knd) :: forcing_tke0 = 0
  real(knd) :: forcing_epsilon0 = 0
  
  real(knd) :: Umean(3) = 0
  
  real(knd), allocatable, dimension(:,:,:) :: Uf, Vf, Wf
  
contains

  subroutine init_linear_forcing(enable, enable_tke0_epsilon0, variable_means, filter, &
                                 forcing_constant, tke0, epsilon0, mean_velocity)
    logical, intent(in) :: enable, enable_tke0_epsilon0, variable_means, filter
    real(knd), intent(in) :: forcing_constant, tke0, epsilon0
    real(knd), intent(in) :: mean_velocity(3)
    
    enable_linear_forcing = enable
    enable_linear_forcing_tke0_epsilon0 = enable_tke0_epsilon0
    enable_forcing_variable_means = variable_means
    enable_filter = filter
    
    forcing_tke0 = tke0
    forcing_epsilon0 = epsilon0
    linear_forcing_constant = forcing_constant
    Umean = mean_velocity
  
    if (enable_linear_forcing) then
      if (enable_linear_forcing_tke0_epsilon0) then
        linear_forcing_constant = forcing_epsilon0 / (2 * forcing_tke0)
      else
        linear_forcing_constant = forcing_constant
      end if
      
      if (.not.enable_forcing_variable_means)  Umean = mean_velocity
    end if
    
    if (enable_filter) then
      allocate(Uf(-2:Unx+3,-2:Uny+3,-2:Unz+3))
      allocate(Vf(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
      allocate(Wf(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
    end if

  end subroutine
  
  
  subroutine  linear_forcing(U2, V2, W2, U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2, V2, W2

    if (enable_linear_forcing) then
      if (enable_linear_forcing_tke0_epsilon0) then
        if (enable_filter) then
          call linear_forcing_tke0_epsilon0_filter(U2, V2, W2, U, V, W)
        else
          call linear_forcing_tke0_epsilon0(U2, V2, W2, U, V, W)
        end if
      else
        call linear_forcing_simple(U2, V2, W2, U, V, W)
      end if
    end if
  end subroutine

  subroutine linear_forcing_simple(U2, V2, W2, U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2, V2, W2
    integer :: i, j, k
    
    if (enable_forcing_variable_means)  call mean_UVW(Umean, U, V, W)
  
    !$omp parallel
    !$omp do private(i,j,k) collapse(3)
    do k = 1, Unz
      do j = 1, Uny
        do i = 1, Unx
          U2(i,j,k) = U2(i,j,k) + (U(i,j,k) - Umean(1)) * linear_forcing_constant
        end do
      end do
    end do
    !$omp do private(i,j,k) collapse(3)
    do k = 1, Vnz
      do j = 1, Vny
        do i = 1, Vnx
          V2(i,j,k) = V2(i,j,k) + (V(i,j,k) - Umean(2)) * linear_forcing_constant
        end do
      end do
    end do
    !$omp do private(i,j,k) collapse(3)
    do k = 1, Wnz
      do j = 1, Wny
        do i = 1, Wnx
          W2(i,j,k) = W2(i,j,k) + (W(i,j,k) - Umean(3)) * linear_forcing_constant
        end do
      end do
    end do
    !$omp end parallel
  end subroutine

  
  subroutine linear_forcing_tke0_epsilon0(U2, V2, W2, U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2, V2, W2
    integer :: i, j, k
    real(knd) :: tke, c
  
    if (enable_forcing_variable_means)  call mean_UVW(Umean, U, V, W)
  
    call mean_TKE(tke, U, V, W, Umean)

    tke = max(0.01_knd * forcing_tke0, min(100._knd * forcing_tke0, tke))

    c = forcing_tke0 / tke
  
    !$omp parallel
    !$omp do private(i,j,k) collapse(3)
    do k = 1, Unz
      do j = 1, Uny
        do i = 1, Unx
          U2(i,j,k) = U2(i,j,k) + (U(i,j,k) - Umean(1)) * linear_forcing_constant * c
        end do
      end do
    end do
  
    !$omp do private(i,j,k) collapse(3)
    do k = 1, Vnz
      do j = 1, Vny
        do i = 1, Vnx
          V2(i,j,k) = V2(i,j,k) + (V(i,j,k) - Umean(2)) * linear_forcing_constant * c
        end do
      end do
    end do
  
    !$omp do private(i,j,k) collapse(3)
    do k = 1, Wnz
      do j = 1, Wny
        do i = 1, Wnx
          W2(i,j,k) = W2(i,j,k) + (W(i,j,k) - Umean(3)) * linear_forcing_constant * c
        end do
      end do
    end do
    !$omp end parallel
  end subroutine

  
  subroutine linear_forcing_tke0_epsilon0_filter(U2, V2, W2, U, V, W)
    use Filters, only: Filter3ord
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2, V2, W2
    integer :: i, j, k
    real(knd) :: tke, c
  
    if (enable_forcing_variable_means)  call mean_UVW(Umean, U, V, W)
  
    call mean_TKE(tke, U, V, W, Umean)
    
    tke = max(0.01_knd * forcing_tke0, min(100._knd * forcing_tke0, tke))

    c = forcing_tke0 / tke
  
    Uf = U
    Vf = V
    Wf = W
    
    call Filter3ord(Uf, Utype, 1)
    call Filter3ord(Uf, Utype, 2)
    call Filter3ord(Uf, Utype, 3)
    call Filter3ord(Vf, Vtype, 1)
    call Filter3ord(Vf, Vtype, 2)
    call Filter3ord(Vf, Vtype, 3)
    call Filter3ord(Wf, Wtype, 1)
    call Filter3ord(Wf, Wtype, 2)
    call Filter3ord(Wf, Wtype, 3)

    !$omp parallel
    !$omp do private(i,j,k) collapse(3)
    do k = 1, Unz
      do j = 1, Uny
        do i = 1, Unx
          U2(i,j,k) = U2(i,j,k) + (Uf(i,j,k) - Umean(1)) * linear_forcing_constant * c
        end do
      end do
    end do
  
    !$omp do private(i,j,k) collapse(3)
    do k = 1, Vnz
      do j = 1, Vny
        do i = 1, Vnx
          V2(i,j,k) = V2(i,j,k) + (Vf(i,j,k) - Umean(2)) * linear_forcing_constant * c
        end do
      end do
    end do
  
    !$omp do private(i,j,k) collapse(3)
    do k = 1, Wnz
      do j = 1, Wny
        do i = 1, Wnx
          W2(i,j,k) = W2(i,j,k) + (Wf(i,j,k) - Umean(3)) * linear_forcing_constant * c
        end do
      end do
    end do
    !$omp end parallel
  end subroutine
  

  subroutine mean_UVW(Umean, U, V, W)
    real(knd), intent(out) :: Umean(3)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
    integer :: i, j, k
    
    Umean = 0
    
    if (n_full_im_U==0.and. &
        n_full_im_V==0.and. &
        n_full_im_W==0) then
      !$omp parallel reduction(+:Umean)
      !$omp do collapse(3)
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            Umean(1) = Umean(1) + U(i,j,k)
          end do
        end do
      end do
      !$omp do collapse(3)
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            Umean(2) = Umean(2) + V(i,j,k)
          end do
        end do
      end do
      !$omp do collapse(3)
      do k = 1, Wnz
        do j = 1, Wny
          do i = 1, Wnx
            Umean(3) = Umean(3) + W(i,j,k)
          end do
        end do
      end do
      !$omp end parallel
    else
      !$omp parallel reduction(+:Umean)
      !$omp do collapse(3)
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            if (Utype(i,j,k)<=0) Umean(1) = Umean(1) + U(i,j,k)
          end do
        end do
      end do
      !$omp do collapse(3)
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            if (Vtype(i,j,k)<=0) Umean(2) = Umean(2) + V(i,j,k)
          end do
        end do
      end do
      !$omp do collapse(3)
      do k = 1, Wnz
        do j = 1, Wny
          do i = 1, Wnx
            if (Wtype(i,j,k)<=0) Umean(3) = Umean(3) + W(i,j,k)
          end do
        end do
      end do
      !$omp end parallel
    end if  
#ifdef PAR
    Umean = par_co_sum(Umean)
#endif

    Umean(1) = Umean(1) / n_free_domain_U
    Umean(2) = Umean(2) / n_free_domain_V
    Umean(3) = Umean(3) / n_free_domain_W
  end subroutine

  subroutine mean_TKE(tke, U, V, W, Umean)
    real(knd), intent(out) :: tke
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
    real(knd), intent(in) :: Umean(3)
    integer :: i, j, k
    
    tke = 0
    
    if (n_full_im_U==0.and. &
        n_full_im_V==0.and. &
        n_full_im_W==0) then
        
      !$omp parallel do private(i,j,k) reduction(+:tke) collapse(3)
      do k = 1, Prnz
        do j = 1, Prny
          do i = 1, Prnx
            tke = tke + (U(i,j,k) + U(i-1,j,k) - 2*Umean(1))**2 &
                      + (V(i,j,k) + V(i,j-1,k) - 2*Umean(2))**2 &
                      + (W(i,j,k) + W(i,j,k-1) - 2*Umean(3))**2
          end do
        end do
      end do
      !$omp end parallel do
      
      tke = tke / 8
    else
      !$omp parallel do private(i,j,k) reduction(+:tke) collapse(3)
      do k = 1, Prnz
        do j = 1, Prny
          do i = 1, Prnx
            if (Prtype(i,j,k)<=0) &
              tke = tke + (U(i,j,k) + U(i-1,j,k) - 2*Umean(1))**2 &
                        + (V(i,j,k) + V(i,j-1,k) - 2*Umean(2))**2 &
                        + (W(i,j,k) + W(i,j,k-1) - 2*Umean(3))**2
          end do
        end do
      end do
      !$omp end parallel do
     
      tke = tke / 8
    end if
       
#ifdef PAR
    tke = par_co_sum(tke)
#endif

    tke = tke / n_free_domain
  end subroutine

end module

