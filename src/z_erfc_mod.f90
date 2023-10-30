module z_erfc

contains


subroutine zfun_erfc(zIn, zOut)
!  calls wofz_f90 and multiplies the result by i * sqrt(pi)

implicit none
integer, parameter :: DBL = selected_real_kind(p=13,r=200)
complex(kind=DBL), intent(in) :: zIn
complex(kind=DBL), intent(out) :: zOut
logical :: flag
real(kind=DBL), parameter :: &
    zero=0.0d0,&
    four=4.0d0,&
    one=1.0d0
!cant use this unless f2003 or later
!real(kind = DBL), parameter :: sqrtpi = sqrt(four * atan(one))
real(kind=DBL), parameter :: sqrtpi = 1.77245385091d0

call wofz_f90(zIn, zOut, flag)
zOut = zOut * cmplx(zero, one) * sqrtpi

return
end subroutine zfun_erfc




subroutine wofz_f90(zIn,zOut,flag)

!  algorithm 680, collected algorithms from acm.
!  this work published in transactions on mathematical software,
!  vol. 16, no. 1, pp. 47.
!
!  given a complex number z = (xi,yi), this subroutine computes
!  the value of the faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
!  where erfc is the complex complementary error-function and i
!  means sqrt(-1).
!  the accuracy of the algorithm for z in the 1st and 2nd quadrant
!  is 14 significant digits; in the 3rd and 4th it is 13 significant
!  digits outside a circular region with radius 0.126 around a zero
!  of the function.
!  all real variables in the program are double precision.
!
!  the code contains a few compiler-dependent parameters :
!     rmaxreal = the maximum value of rmaxreal equals the root of
!                rmax = the largest number which can still be
!                implemented on the computer in double precision
!                floating-point arithmetic
!     rmaxexp  = ln(rmax) - ln(2)
!     rmaxgoni = the largest possible argument of a double precision
!                goniometric function (dcos, dsin, ...)
!  the reason why these parameters are needed as they are defined will
!  be explained in the code by means of comments
!
!
!  parameter list
!     xi     = real      part of z
!     yi     = imaginary part of z
!     u      = real      part of w(z)
!     v      = imaginary part of w(z)
!     flag   = an error flag indicating whether overflow will
!              occur or not; type logical;
!              the values of this variable have the following
!              meaning :
!              flag=.false. : no error condition
!              flag=.true.  : overflow will occur, the routine
!                             becomes inactive
!  xi, yi      are the input-parameters
!  u, v, flag  are the output-parameters
!
!  furthermore the parameter factor equals 2/sqrt(pi)
!
!  the routine is not underflow-protected but any variable can be
!  put to 0 upon underflow;
!
!  reference - gpm poppe, cmj wijers; more efficient computation of
!  the complex error-function, acm trans. math. software (1990)
!  vol. 16, no. 1, pp. 47.


implicit none
integer, parameter :: DBL = selected_real_kind(p=13,r=200)
complex(kind=DBL), intent(in) :: zIn
complex(kind=DBL), intent(out) :: zOut
logical, intent(out) :: flag

real(kind=DBL) :: c, daux, h, h2=0, qlambda=0, qrho,&
    rx, ry, sx,&
    sy, tx, ty, u, u1, u2, v
real(kind=DBL) :: v1, v2, w1, x, xabs, xabsq, xaux, xi,&
    xquad, xsum, y, yabs, yi, yquad, ysum
integer :: i, j, kapn, n, np1, nu
logical :: a, b
real(kind=DBL), parameter :: &
    factor=1.12837916709551257388d0, &
    rmaxreal=0.5d+154,&
    rmaxexp=708.503061461606d0,&
    rmaxgoni=3.53711887601422d+15,&
    zero=0.0d0,&
    one=1.0d0

xi = real(zIn)
yi = aimag(zIn)
flag = .false.
xabs = dabs(xi)
yabs = dabs(yi)
x = xabs/6.3
y = yabs/4.4
!print*, 'I am in the right routine'
! The following if-statement protects
! qrho = (x**2 + y**2) against overflow

overflowLabel1: &
if ( (xabs.gt.rmaxreal) .or. (yabs.gt.rmaxreal) ) then

   flag = .true.

else

    qrho = x**2 + y**2

    xabsq = xabs**2
    xquad = xabsq - yabs**2
    yquad = 2*xabs*yabs

    a = qrho.lt.0.085264d0

    qrhoBreakdown: &
    if ( a ) then
        !print*, 'I am in the power series'

        !  if (qrho.lt.0.085264d0) then the faddeeva-function is evaluated
        !  using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
        !  n is the minimum number of terms needed to obtain the required
        !  accuracy

        qrho = (1-0.85*y)*dsqrt(qrho)
        n = idnint(6+72*qrho)
        j = 2*n + 1
        xsum = one/j
        ysum = zero

        do i = n , 1 , -1
           j = j - 2
           xaux = (xsum*xquad-ysum*yquad)/i
           ysum = (xsum*yquad+ysum*xquad)/i
           xsum = xaux + one/j
        enddo

        u1 = -factor*(xsum*yabs+ysum*xabs) + 1.0
        v1 = factor*(xsum*xabs-ysum*yabs)
        daux = dexp(-xquad)
        u2 = daux*dcos(yquad)
        v2 = -daux*dsin(yquad)

        u = u1*u2 - v1*v2
        v = u1*v2 + v1*u2

    else

    !  if (qrho.gt.1.o) then w(z) is evaluated using the laplace
    !  continued fraction
    !  nu is the minimum number of terms needed to obtain the required
    !  accuracy
    !
    !  if ((qrho.gt.0.085264d0).and.(qrho.lt.1.0)) then w(z) is evaluated
    !  by a truncated taylor expansion, where the laplace continued fraction
    !  is used to calculate the derivatives of w(z)
    !  kapn is the minimum number of terms in the taylor expansion needed
    !  to obtain the required accuracy
    !  nu is the minimum number of terms of the continued fraction needed
    !  to calculate the derivatives with the required accuracy

        if ( qrho.gt.one ) then

            h       = zero
            kapn    = 0
            qrho    = dsqrt(qrho)
            nu      = idint(3+(1442/(26*qrho+77)))
            !nu      = idint(3+(1442/(26*qrho+77)))/10
            !print*, 'in qrho gt one'

        else

            qrho    = (one-y)*dsqrt(one-qrho)
            h       = 1.88*qrho
            h2      = 2*h
            kapn    = idnint(7+34*qrho)
            nu      = idnint(16+26*qrho)
            !print*, 'in qrho gt one else'

        endif

        b = (h.gt.zero)

        if ( b ) qlambda = h2**kapn

        rx = zero
        ry = zero
        sx = zero
        sy = zero

        do n = nu, 0, -1

           np1  = n + 1
           tx   = yabs + h + np1*rx
           ty   = xabs - np1*ry
           c    = 0.5d0/(tx**2+ty**2)
           rx   = c*tx
           ry   = c*ty

           if ( (b) .and. (n.le.kapn) ) then

              tx = qlambda + sx
              sx = rx*tx - ry*sy
              sy = ry*tx + rx*sy
              qlambda = qlambda/h2

           endif

        enddo

        if ( h.eq.zero) then
           u = factor*rx
           v = factor*ry
        else
           u = factor*sx
           v = factor*sy
        endif

        if ( yabs.eq.zero) u = dexp(-xabs**2)

    endif qrhoBreakdown

    ! Evaluation of w(z) in the other quadrants

    yiLTzero: &
    if ( yi.lt.zero) then

        aLabel2: &
        if ( a ) then

            u2 = 2*u2
            v2 = 2*v2

        else

            xquad = -xquad

            ! The following if-statement protects 2*exp(-z**2)
            ! against overflow

            overflowLabel: &
            if ( (yquad.gt.rmaxgoni) .or. (xquad.gt.rmaxexp) ) then

                flag = .true.
                zOut = cmplx(u,v,DBL)
                return

            else

                w1 = 2*dexp(xquad)
                u2 = w1*dcos(yquad)
                v2 = -w1*dsin(yquad)

            endif overflowLabel

        endif aLabel2

        u = u2 - u
        v = v2 - v

        if ( xi.gt.zero) v = -v

        else
            if ( xi.lt.zero) v = -v

        endif yiLTzero

    zOut = cmplx(u,v,DBL)
    return

endif overflowLabel1

end subroutine wofz_f90

end module z_erfc
