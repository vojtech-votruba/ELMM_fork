module Limiters

 use Parameters, only : knd
 
 implicit none

  private
  public FluxLimiter, Limiter, limiter_type, limiter_parameter
  
  integer, parameter :: minmodlim = 1, extminmodlim = 2, gammalim = 3, &
                        vanalbadalim = 4, vanleerlim = 5, superbeelim = 6
  integer   :: limiter_type
  real(knd) :: limiter_parameter

 contains

  elemental real(knd) function FluxLimiter(r)
    real(knd),intent(in):: r
    real(knd) :: K

    K=(1+2._knd*r)/3._knd  !assuming kappa=1/3 scheme (Koren, 1993; Hundsdorfer, 1994)
    if (limiter_parameter>0) then
     FluxLimiter=max(0._knd,min(2._knd*r,min(limiter_parameter,K)))
    elseif (limiter_parameter<0) then
     FluxLimiter=K
    else
     FluxLimiter=0
    endif
  endfunction FluxLimiter


  elemental real(knd) function Limiter(a,b,c)
    real(knd),intent(in):: a,b
    real(knd),intent(in),optional:: c
    real(knd) :: R
    real(knd),parameter:: epsil=0.00001_knd

    if (limiter_type==minmodlim) then
     Limiter=MinMod(a,b)
    elseif (limiter_type==extminmodlim) then
     if (.not.present(c)) then
       Limiter=MinMod(limiter_parameter*a,limiter_parameter*b,(a+b)/2._knd)
     else
       Limiter=MinMod(limiter_parameter*a,limiter_parameter*b,c)
     endif
    elseif (limiter_type==gammalim) then
     if (abs(b)>epsil.and.abs(a)>epsil) then
      R=a/b
      Limiter=b*GammaLimiter(R)
     else
      Limiter=MinMod(a,b)
     endif
    elseif (limiter_type==vanalbadalim) then
     if (abs(b)>epsil) then
      R=a/b
      Limiter=b*VanAlbada(R)
     else
      Limiter=MinMod(a,b)
     endif
    elseif (limiter_type==vanleerlim) then
     if (abs(b)>epsil) then
      R=a/b
      Limiter=b*VanLeer(R)
     else
      Limiter=MinMod(a,b)
     endif
    elseif (limiter_type==superbeelim) then
     if (abs(b)>epsil) then
      R=a/b
      Limiter=b*SuperBee(R)
     else
      Limiter=MinMod(a,b)
     endif
    elseif (limiter_type==0) then
     Limiter=0
    else
     Limiter=(a+b)/2._knd
    endif
  endfunction Limiter

  elemental real(knd) function Heaviside(x)
    real(knd),intent(in):: x
    Heaviside=(1._knd+SIGN(1._knd,x))/2._knd
  endfunction Heaviside

  elemental real(knd) function GammaLimiter(R)
    real(knd),intent(in):: R
    GammaLimiter =   ((1-limiter_parameter)/limiter_parameter) &
                   * R                                         &
                   * (Heaviside(R) - Heaviside(r-(limiter_parameter/(1-limiter_parameter)))) &
                   + Heaviside(r-(limiter_parameter/(1-limiter_parameter)))
  endfunction GammaLimiter

  elemental real(knd) function VanAlbada(R)
    real(knd),intent(in):: R
    VanAlbada=(R+R*R)/(1+R*R)
  endfunction VanAlbada

  elemental real(knd) function VanLeer(R)
    real(knd),intent(in):: R
    VanLeer=(R+abs(R))/(1+abs(R))
  endfunction VanLeer

  elemental real(knd) function SuperBee(R)
    real(knd),intent(in):: R
    SuperBee=MAX(0._knd,MAX(MIN(2*R,1._knd),MIN(R,2._knd)))
  endfunction SuperBee

  elemental real(knd) function MinMod(a,b,c)
    real(knd),intent(in):: a,b
    real(knd),intent(in),optional:: c

    if (.not. present(c)) then

     MinMod=(SIGN(1._knd,a)+SIGN(1._knd,b))*MIN(ABS(a),ABS(b))/2._knd

    else

     if ((a>0).and.(b>0).and.(c>0)) then
         MinMod=MIN(a,b,c)
     elseif ((a<0).and.(b<0).and.(c<0)) then
         MinMod=MAX(a,b,c)
     else
         MinMod=0
     endif
    endif
  endfunction MinMod

end module Limiters

