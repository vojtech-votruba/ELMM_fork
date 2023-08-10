module interpolation
  use Kinds

  implicit none

contains

  !*****************************************************************************
  !
  !  cubic_spline() defines an interpolatory cubic spline.
  !
  !  Discussion:
  !
  !    A tridiagonal linear system for the unknown slopes S(I) of
  !    F at x(I), I=1,..., N, is generated and then solved by Gauss
  !    elimination, with S(I) ending up in C(2,I), for all I.
  !
  !  Modified:
  !
  !    20 March 2015 Vladimir Fuka (modernization and adaptation of ver. 14 February 2007) 
  !
  !  Author:
  !
  !    Carl de Boor
  !
  !  Reference:
  !
  !    Carl de Boor,
  !    A Practical Guide to Splines,
  !    Springer, 2001,
  !    ISBN: 0387953663,
  !    LC: QA1.A647.v27.
  !
  subroutine cubic_spline(x, y, c, dydx_1, dydx_n, d2ydx2_1, d2ydx2_n)
    real(knd), intent(in) :: x(:), y(:)
    real(knd), intent(out) :: c(:,:)
    real(knd), intent(in), optional :: dydx_1, dydx_n, d2ydx2_1, d2ydx2_n
    integer :: bc1
    integer :: bcn
    real(knd) :: divdf1
    real(knd) :: divdf3
    real(knd) :: dx
    real(knd) :: g
    integer :: n, i

    n = size(x)

    c = 0
    
    c(1,:) = y
    
    if (present(dydx_1)) then
      bc1 = 1
      c(2,1) = dydx_1
    else if (present(d2ydx2_1)) then
      bc1 = 2
      c(2,1) = d2ydx2_1
    else
      bc1 = 0
    end if
    
    if (present(dydx_n)) then
      bcn = 1
      c(2,n) = dydx_n
    else if (present(d2ydx2_n)) then
      bcn = 2
      c(2,n) = d2ydx2_n
    else
      bcn = 0
    end if
    
  !
  !  C(3,*) and C(4,*) are used initially for temporary storage.
  !
  !  Store first differences of the x sequence in C(3,*).
  !
  !  Store first divided difference of data in C(4,*).
  !
    do i = 2, n
      c(3,i) = x(i) - x(i-1)
    end do

    do i = 2, n 
      c(4,i) = ( c(1,i) - c(1,i-1) ) / ( x(i) - x(i-1) )
    end do
  !
  !  Construct the first equation from the boundary condition
  !  at the left endpoint, of the form:
  !
  !    C(4,1) * S(1) + C(3,1) * S(2) = C(2,1)
  !
  !  bc1 = 0: Not-a-knot
  !
    if ( bc1 == 0 ) then

      if ( n <= 2 ) then
        c(4,1) = 1
        c(3,1) = 1
        c(2,1) = 2 * c(4,2)
        call n2
        return
      end if

      c(4,1) = c(3,3)
      c(3,1) = c(3,2) + c(3,3)
      c(2,1) = ( ( c(3,2) + 2 * c(3,1) ) * c(4,2) * c(3,3) &
        + c(3,2)**2 * c(4,3) ) / c(3,1)
  !
  !  bc1 = 1: derivative specified.
  !
    else if ( bc1 == 1 ) then

      c(4,1) = 1
      c(3,1) = 0

      if ( n == 2 ) then
        call n2
        return
      end if
  !
  !  Second derivative prescribed at left end.
  !
    else

      c(4,1) = 2
      c(3,1) = 1
      c(2,1) = 3 * c(4,2) - c(3,2) / 2 * c(2,1)

      if ( n == 2 ) then
        call n2
        return
      end if

    end if
  !
  !  If there are interior knots, generate the corresponding
  !  equations and carry out the forward pass of Gauss elimination,
  !  after which the I-th equation reads:
  !
  !    C(4,I) * S(I) + C(3,I) * S(I+1) = C(2,I).
  !
    do i = 2, n-1
      g = -c(3,i+1) / c(4,i-1)
      c(2,i) = g * c(2,i-1) + 3 * ( c(3,i) * c(4,i+1) + c(3,i+1) * c(4,i) )
      c(4,i) = g * c(3,i-1) + 2 * ( c(3,i) + c(3,i+1))
    end do

  !
  !  Construct the last equation from the second boundary condition, of
  !  the form
  !
  !    -G * C(4,N-1) * S(N-1) + C(4,N) * S(N) = C(2,N)
  !
  !  If slope is prescribed at right end, one can go directly to
  !  back-substitution, since the C array happens to be set up just
  !  right for it at this point.
  !
    if ( bcn == 1 ) then

      continue

    else if ( bcn > 1 ) then

      c(2,n) = 3 * c(4,n) + c(3,n) / 2 * c(2,n)
      c(4,n) = 2
      g = -1 / c(4,n-1)
      c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
      c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)

    else if ( n /= 3 .or. bc1 /= 0 ) then

      !
      !  Not-a-knot and 3 <= N, and either 3 < N or also not-a-knot
      !  at left end point.
      !
      g = c(3,n-1) + c(3,n)
      c(2,n) = ( ( c(3,n) + 2 * g ) * c(4,n) * c(3,n-1) + c(3,n)**2 &
        * ( c(1,n-1) - c(1,n-2) ) / c(3,n-1) ) / g
      g = - g / c(4,n-1)
      c(4,n) = c(3,n-1)
      c(4,n) = c(4,n) + g * c(3,n-1)
      c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)

    else

      !
      !  N = 3 and not-a-knot also at left.
      !
      c(2,n) = 2 * c(4,n)
      c(4,n) = 1
      g = -1 / c(4,n-1)
      c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
      c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)

    end if

    call back_solve

  contains

    subroutine n2
      !
      !  N = 2.
      !
      
      if ( bcn == 2  ) then

        c(2,n) = 3 * c(4,n) + c(3,n) / 2 * c(2,n)
        c(4,n) = 2
        g = -1 / c(4,n-1)
        c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
        c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
     
      else if ( bcn == 0 .and. bc1 /= 0 ) then

        c(2,n) = 2 * c(4,n)
        c(4,n) = 1
        g = -1 / c(4,n-1)
        c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
        c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)

      else if ( bcn == 0 .and. bc1 == 0 ) then

        c(2,n) = c(4,n)

      end if

      call back_solve
    end subroutine

    subroutine back_solve
      !
      !  Back solve the upper triangular system 
      !
      !    C(4,I) * S(I) + C(3,I) * S(I+1) = B(I)
      !
      !  for the slopes C(2,I), given that S(N) is already known.
      !
      do i = n-1, 1, -1
        c(2,i) = ( c(2,i) - c(3,i) * c(2,i+1) ) / c(4,i)
      end do
      !
      !  Generate cubic coefficients in each interval, that is, the
      !  derivatives at its left endpoint, from value and slope at its
      !  endpoints.
      !
      do i = 2, n
        dx = c(3,i)
        divdf1 = ( c(1,i) - c(1,i-1) ) / dx
        divdf3 = c(2,i-1) + c(2,i) - 2 * divdf1
        c(3,i-1) = ( divdf1 - c(2,i-1) - divdf3 ) / dx
        c(4,i-1) = divdf3 / dx**2
      end do
    end subroutine

  end subroutine cubic_spline


  function cubic_spline_eval(x, xs, c, i) result(res)
    real(knd) :: res
    real(knd), intent(in) :: x
    real(knd), intent(in) :: xs(:)
    real(knd), intent(in) :: c(0:,:)
    integer, intent(inout) :: i
    real(knd) :: h
    integer :: n

    n = size(xs)
    i = min(max(1, i), n-1)

    if (x<xs(1)) then
      res = c(0, 1)
      return
    end if

    if (x>xs(n)) then
      res = c(0, n)
      return
    end if

    !if (x < xs(i))
    do while (x < xs(i))
      i = i - 1
    end do

    !if (x > xs(i+1))
    do while (x > xs(i+1))
      i = i + 1
    end do

    h = x - xs(i)
    res = c(0,i) + h * (c(1,i) + h * (c(2,i) + h * c(3,i)))

  end function cubic_spline_eval
  

  
  
  subroutine linear_interpolation(x, y, c)
    real(knd), intent(in) :: x(:), y(:)
    real(knd), intent(out) :: c(0:,:)
    integer :: i
    
    c(0, :) = y
    do i = 1, size(x)-1
      c(1, i) = (y(i+1)-y(i)) / (x(i+1)-x(i))
    end do
    c(1,size(x)) = 0
  end subroutine
  
  function linear_interpolation_eval(x, xs, c, i) result(res)
    real(knd) :: res
    real(knd), intent(in) :: x
    real(knd), intent(in) :: xs(:)
    real(knd), intent(in) :: c(0:,:)
    integer, intent(inout) :: i
    real(knd) :: h
    integer :: n

    n = size(xs)
    i = min(max(1, i), n-1)

    if (x<xs(1)) then
      res = c(0, 1)
      return
    end if

    if (x>xs(n)) then
      res = c(0, n)
      return
    end if

    !if (x < xs(i))
    do while (x < xs(i))
      i = i - 1
    end do

    !if (x > xs(i+1))
    do while (x > xs(i+1))
      i = i + 1
    end do

    h = x - xs(i)
    res = c(0,i) + h * c(1,i)

  end function

end module
