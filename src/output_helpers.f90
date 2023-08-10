module Output_helpers

  use Kinds
  use Parameters

contains


  function TotKE(U,V,W) result(res)
#ifdef PAR
  use custom_par
#endif
    real(knd) :: res
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd) :: Um,Vm,Wm
    integer :: i,j,k
    real(knd) :: lx,ly,lz

    lx = gxmax - gxmin
    ly = gymax - gymin
    lz = gzmax - gzmin

    res = 0

    !$omp parallel do private(i,j,k) reduction(+:res)
    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       res = res + (U(i-1,j,k) + U(i,j,k))**2 + &
                   (V(i,j-1,k) + V(i,j,k))**2 + &
                   (W(i,j,k-1) + W(i,j,k))**2
      end do
     end do
    end do
    !$omp end parallel do
    
#ifdef PAR
    res = par_co_sum(res)
#endif    

    !2**2 for the average above, 2 from the definition of kinetic energy
    res = res / 8
    !average
    res = res / (gPrnx*gPrny*gPrnz)
    !integrate

    res = res*lx*ly*lz
  end function TotKE

  pure real(knd) function VorticityMag(i,j,k,U,V,W) result(res)
    integer,intent(in) :: i,j,k
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W

    res = sum(Vorticity(i,j,k,U,V,W)**2)
    res = Sqrt(res)
  end function VorticityMag


  pure real(knd) function Lambda2(i,j,k,U,V,W) result(res)
    integer,intent(in) :: i,j,k
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W

    res = ((U(i,j,k)-U(i-1,j,k))/dxmin)**2
    res = res + ((V(i,j,k)-V(i,j-1,k))/dymin)**2
    res = res + ((W(i,j,k)-W(i,j,k-1))/dzmin)**2
    res = res + ((V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4*dxmin))**2
    res = res + ((W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4*dxmin))**2
    res = res + ((U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4*dymin))**2
    res = res + ((W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4*dymin))**2
    res = res + ((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4*dzmin))**2
    res = res + ((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4*dzmin))**2
    res = -Sqrt(res)
    res = VorticityMag(i,j,k,U,V,W) + res
  end function Lambda2
  
  pure real(knd) function ScalarVerticalFlux(i,j,k,Scal,W)
    !not accurate!  
    integer, intent(in)   :: i,j,k
    real(knd), contiguous, intent(in) :: Scal(-1:,-1:,-1:), W(-2:,-2:,-2:)

    ScalarVerticalFlux = Scal(i,j,k) * (W(i,j,k)+W(i,j,k-1))/2 + &
                       TDiff(i,j,k) * (Scal(i,j,k+1)-Scal(i,j,k-1)) / (zW(k+1)-zW(k-1))
  end function ScalarVerticalFlux

  pure function Vorticity(i,j,k,U,V,W)
    real(knd), dimension(3) :: Vorticity
    integer, intent(in)     :: i,j,k
    real(knd),dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W

    Vorticity = [         (W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                              -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin), &
                           (U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                             -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin), &
                           (V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                             -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin) ]
  end function Vorticity

end module



module Directories
  use Parameters, only: output_dir
  
  implicit none
  
contains
  
  subroutine create_directory(dir)
    character(*), intent(in) :: dir
    integer :: k, io, u
    
#if defined(_WIN32) || defined(_WIN64)
   call system("mkdir "//dir)
#else
   do k = 1, 1000
     call system("mkdir -p "//dir)
     open(newunit=u, file=dir//"test", status="replace", iostat=io)
     if (io==0) then
       close(u, status="delete")
       exit
     end if
     call sleep(1)
   end do
#endif
  end subroutine
end module

