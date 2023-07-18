!module METS2AORSA_MYRA

!  private
!  public::  GET_WMAT_MYRA, intplt1d, WINTERP1D_3

!contains

  subroutine GET_SUM_MYRA(k_uper, sum, W, ZSPEC, ASPEC, BMAG, &
       & lmax, ENORM, UPARMIN, UPARMAX, &
       & NUPAR, NUPER, UPER, UPAR, DFDUPER, DFDUPAR,  &
       & ealphak, ebetak, ebk, nkdim1, nkdim2, mkdim1, mkdim2,   &
       & nkx1, nkx2, nky1, nky2, &
       & uxx, uxy, uxz, &
       & uyx, uyy, uyz, &
       & uzx, uzy, uzz, &
       & nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i_global, j_global, &
       & lmaxdim, ndist, nzeta)


    implicit none

    integer:: i_global, j_global, ier, nxdim, nydim, k, lmaxdim, ndist
    integer, intent(IN):: NUPAR, NUPER, lmax
    integer:: nkx1, nkx2, nky1, nky2, j_upar, k_uper
    integer:: nkdim1, nkdim2, mkdim1, mkdim2
    integer:: NHARM, IHARM, M, N, i, nzeta

    real::  uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz
    real::  xkphi
    real::  xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
    real::  xkperpn, xkperpni, xkrhon, xketan, xkbn, beta
    real, intent(IN):: W, ZSPEC, ASPEC, BMAG
    real, intent(IN):: ENORM, UPARMIN, UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER, DFDUPAR
    real:: W2, WCW, RRP, RRM, WC, WCI
    real:: MUT0, SQMUT0, SQMUT0I, KPARA1, NPARA1
    real:: ISQ2, SQ2, NWCW, DFACTPAR, DFACTPER, U0, U
    real:: UPAR0, dfdupar0, dfduper0, du, dui, p
    real:: time, t1, tmsec, second1, dummy, WI, uperpk
    real:: dzeta, dzetai, zetamax, zetamin, zeta0, zetamax1, zetamax2
    real:: A1, A2, A3
    real:: temp1, temp2, temp3

    real, dimension(:),     allocatable :: zetai
    real, dimension(:,:),   allocatable :: Jni
    real, dimension(:,:,:), allocatable :: Jn
    real, dimension(:,:),   allocatable :: cosbeta
    real, dimension(:,:),   allocatable :: sinbeta
    real, dimension(:,:),   allocatable :: NPARA_sav

    complex::epsx, epsy, epsz
    complex::ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
    &        ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
    &           ebk(nkdim1 : nkdim2, mkdim1 : mkdim2)
    complex::cexpkx, cexpky, zi, zeta

    complex::xx(nkdim1 : nkdim2, 1 : nxdim),   &
     &      yy(mkdim1 : mkdim2, 1 : nydim)

    complex::cexpn, cexpnp1, cexpnm1, cexp11

    complex::cexp1, cexp2, cexp0
    complex::sum1_11, sum1_31
    complex::sum2_1, sum2_2, sum2_3, sum
    complex::b(100)


    real, parameter:: EOVERAMU = 9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
    
    real, parameter:: MPC2 = 938271998.38
    real, parameter:: C = 2.99792458e8
    real, parameter:: PI = 3.141592653597932384


    allocate(zetai(nzeta + 1) )
    allocate(Jni(-lmaxdim : lmaxdim, nzeta + 1) )
    allocate( Jn(-lmaxdim : lmaxdim, nkdim1 : nkdim2, mkdim1 : mkdim2))
    allocate(cosbeta(nkdim1 : nkdim2, mkdim1 : mkdim2) )
    allocate(sinbeta(nkdim1 : nkdim2, mkdim1 : mkdim2) )
    allocate(NPARA_sav(nkdim1 : nkdim2, mkdim1 : mkdim2) )

    zi = cmplx(0., 1.)

    uperpk = uper(k_uper)

    W2 = W * W
    WI = 1.0 / W

    WCW = BMAG * ZSPEC * EOVERMH / ASPEC / W
    WC = WCW * W
    WCI = 1.0 /WC

    MUT0 = 0.5 * MPC2 * ASPEC / ENORM
    SQMUT0 = SQRT(MUT0)
    SQMUT0I = 1.0 / SQMUT0

    ISQ2 = SQRT(0.5)
    SQ2 = SQRT(2.0)
    NHARM = lmax

    du = (upar(nupar) - upar(1)) / (nupar - 1)
    dui = 1.0 / du

    if(nzeta .eq. 1)then

    ! -------------------------------------------------------- !
    ! ---Don't interpolate: precalculate all Bessel functions-- !
    ! -------------------------------------------------------- !

       do n = nkx1, nkx2
          do m = nky1, nky2

             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             xkbn   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi
        
             NPARA_sav(n, m) = xkbn * C / W

             xkperpn = sqrt(xkrhon**2 + xketan**2)
             if(xkperpn .eq. 0.0)xkperpn = 1.0e-08
        
             cosbeta(n, m) = xkrhon / xkperpn
             sinbeta(n, m) = xketan / xkperpn

             zeta = xkperpn * uper(k_uper) * c * sqmut0i / wc

             call besjc(zeta, nharm + 2, b, ier)
             if(ier .ne. 0) write(6, *) "ier = ", ier

             do IHARM = 0, NHARM + 1
                Jn(iharm,  n, m) = b(iharm + 1)
                Jn(-iharm, n, m) = (-1.0)**iharm * Jn(iharm, n, m)
             end do

          end do
       end do

    else

!      t1 = second1(dummy)

       ! -------------------------------------- !
       ! ---Interpolate; calculate zeta mesh -- !
       ! -------------------------------------- !

       zetamax = 0.0
       zetamin = 0.0

       do n = nkx1, nkx2
          do m = nky1, nky2

             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             xkbn   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi

             xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08
        
             zeta0 = xkperpn * uperpk * c * sqmut0i * wci
        
             if (zeta0 .gt. zetamax) zetamax = zeta0
             if (zeta0 .lt. zetamin) zetamin = zeta0
          end do
       end do


        if(zetamax .eq. zetamin)then
          zetamax =  1.0e-06
          zetamin = -1.0e-06
       end if

       dzeta = (zetamax - zetamin) / (nzeta - 1)
       dzetai = 1.0 / dzeta

       ! ------------------------------------------------- !
       ! ---Pre-calculate Bessel functions on zeta mesh -- !
       ! ------------------------------------------------- !
        
       do i = 1, nzeta + 1
          zetai(i) = zetamin + (i - 1) * dzeta
          zeta = cmplx(zetai(i), 0.0)
        
          call besjc(zeta, nharm + 2, b, ier)
!          if(ier .ne. 0) write(6, *) "ier = ", ier
        
          do iharm = 0, NHARM + 1
             Jni(iharm,  i) = b(iharm + 1)
             Jni(-iharm, i) = (-1.0)**iharm * b(iharm + 1)
          end do
       end do





       ! --------------------------------- !
       ! ---Interpolate Bessel functions-- !
       ! --------------------------------- !

       do n = nkx1, nkx2
          do m = nky1, nky2

             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             xkbn   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi
        
             NPARA_sav(n, m) = xkbn * C * WI

             xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08
             xkperpni = 1.0 / xkperpn
        
             cosbeta(n, m) = xkrhon * xkperpni
             sinbeta(n, m) = xketan * xkperpni
        
             zeta0 = xkperpn * uperpk * c * sqmut0i * wci
        
             i = int((zeta0 - zetamin) * dzetai) + 1
             p = (zeta0 - zetai(i)) * dzetai
             A1 = 0.5 * P * (P - 1.)
             A2 = 1. - P * P
             A3 = 0.5 * P * (P + 1.)
        
             do iharm = -NHARM - 1, NHARM + 1
        
                Jn(iharm, n, m) = Jni(iharm, i) + p * (Jni(iharm, i + 1) - Jni(iharm, i))
                if(i .ne. 1 )then
                   Jn(iharm, n, m) = A1 * Jni(iharm, i - 1) + A2 * Jni(iharm, i) + A3 * Jni(iharm, i + 1)
                end if
        
             end do
        
          end do
       end do


    end if



    ! ------------------------ !
    ! ---Sum over harmonics--- !
    ! ------------------------ !

    sum = 0.0

    do IHARM = -NHARM, NHARM

       NWCW = real(IHARM) * WCW


       sum1_11 = 0.0
       sum1_31 = 0.0

       sum2_1 = 0.0
       sum2_2 = 0.0
       sum2_3 = 0.0

      ! --------------------------- !
      ! ---Sum over Fourier modes-- !
      ! --------------------------- !

       do n = nkx1, nkx2

          do m = nky1, nky2

             NPARA1 = NPARA_sav(n, m)
        
             cexpkx = xx(n, i_global)
             cexpky = yy(m, j_global)

             cexp11 = cosbeta(n, m) + zi * sinbeta(n, m)
             cexpn  = cexp11 ** iharm
                
             cexpnp1 = cexpn * cexp11
             cexpnm1 = cexpn / cexp11

             cexp1 = cexpkx * cexpky * cexpnp1
             cexp2 = cexpkx * cexpky * cexpnm1
             cexp0 = cexpkx * cexpky * cexpn
                
             epsx = isq2 * (ealphak(n, m) - zi * ebetak(n, m)) * cexp1
             epsy = isq2 * (ealphak(n, m) + zi * ebetak(n, m)) * cexp2
             epsz = ebk(n, m) * cexp0
        
             sum2_1 = sum2_1 + conjg(epsx) * Jn(IHARM + 1, n, m)
             sum2_2 = sum2_2 + conjg(epsy) * Jn(IHARM - 1, n, m)
             sum2_3 = sum2_3 + conjg(epsz) * Jn(IHARM, n, m)

             ! ------------------------ !
             ! -- Resonance relation -- !
             ! ------------------------ !
             RRP = 1.0 - NWCW - NPARA1 * UPARMAX * SQMUT0i
             RRM = 1.0 - NWCW - NPARA1 * UPARMIN * SQMUT0i

             if (RRP * RRM .le. 0) then

                ! --------------- !
                ! -- Resonance -- !
                ! --------------- !

                UPAR0 = SQMUT0 / NPARA1 * (1. - NWCW)
        
                i = int((UPAR0 - UPAR(1)) * dui) + 1
                p = (UPAR0 - UPAR(i)) * dui

                dfduper0 = dfduper(k_uper, NUPAR)
                if (i .ne. NUPAR) then
                   dfduper0 = dfduper(k_uper, i) + (dfduper(k_uper, i+1) - dfduper(k_uper, i)) * p
                end if
                
                U0 = DFDUPER0
        
                if(ndist .eq. 1) then  !--Non-Maxwellian--!

                   DFACTPAR = NPARA1 * UPAR0 * SQMUT0I
                   DFACTPER = NPARA1 * UPER(k_uper) * SQMUT0I

                   dfdupar0 = dfdupar(k_uper, NUPAR)
                   if(i .ne. NUPAR)then
                      dfdupar0 = dfdupar(k_uper, i) + (dfdupar(k_uper, i+1) - dfdupar(k_uper, i)) * p
                   end if
        
                   U0 = (1. - DFACTPAR) * DFDUPER0 + DFACTPER * DFDUPAR0
                
                end if
        
        
                U  = PI * SQMUT0 / abs(NPARA1) * U0


                temp1 = UPER(k_uper) * UPER(k_uper) * U
                temp2 = SQ2 * UPER(k_uper) * UPAR0 * U
                temp3 = 2.0 * UPAR0 * UPAR0 * U
        

                
                sum1_11 = sum1_11 + temp1 * Jn(IHARM + 1, n, m) * epsx  &
     &                            + temp1 * Jn(IHARM - 1, n, m) * epsy  &  
     &                            + temp2 * Jn(IHARM, n, m)     * epsz

                sum1_31 = sum1_31 + temp2 * Jn(IHARM + 1, n, m) * epsx  &
     &                            + temp2 * Jn(IHARM - 1, n, m) * epsy  &
     &                            + temp3 * Jn(IHARM, n, m)     * epsz          
        
             end if

          end do
       end do
     
       sum = sum + sum2_1 * sum1_11 + sum2_2 * sum1_11 + sum2_3 * sum1_31 

    end do



    deallocate(zetai)
    deallocate(Jni)
    deallocate(Jn)
    deallocate(cosbeta)
    deallocate(sinbeta)
    deallocate(NPARA_sav)


  end subroutine GET_SUM_MYRA


!
!*************************************************************************
!

    subroutine WINTERP1D_3(X, LX, FX, X0, FX0)

!   --------------------------
!   1D quadradic interpolation
!   --------------------------

    implicit none

    integer, intent(IN):: LX
    real, dimension(LX), intent(IN):: X,FX
    real, intent(IN):: X0
    real, intent(inout):: FX0

    real:: FXY0,XMIN,XMAX,DX
    real:: P,A1,A2,A3
    integer:: I

    FX0 = 0.
    XMIN = X(1)
    XMAX = X(LX)
    DX = (XMAX - XMIN) / (LX - 1)

    I = int(1. + (X0 - XMIN) / DX)
    P = (X0 - X(I)) / DX

    if(i .ne. 1 .and. i .ne. lx)then
       A1 = 0.5 * P * (P - 1.)
       A2 = 1. - P * P
       A3 = 0.5 * P * (P + 1.)
       FX0 = A1 * FX(I-1) + A2 * FX(I) + A3 * FX(I+1)

    else
       if(i .eq. lx) fx0 = fx(i)
       if(i .eq. 1)  fx0 = fx(i) + (fx(i+1) - fx(i)) * p
    end if


    end subroutine WINTERP1D_3

!
!*************************************************************************
!

       subroutine intplt1d(x, lx, fx, nlen, x0, fx0)

!      -----------------------
!      1D linear interpolation
!      -----------------------

       implicit none

       integer:: lx, i, nlen
       real:: fx(lx), x(lx), fx0, x0, dx
      
       do i = 1, nlen
          if (x0 .ge. x(i) .and. x0 .lt. x(i+1) ) exit
       end do
       
       if(i == nlen) then
          fx0 = fx(nlen)
       else
          dx = x(i+1) - x(i)
          fx0 = fx(i) + (fx(i+1) - fx(i)) * (x0 - x(i)) / dx
       end if

       end subroutine intplt1d

!
!*************************************************************************
!
   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer:: n
double precision x(n), y(n), b(n), c(n), d(n)
integer:: i, j, gap
double precision h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions 
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination 
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline
!
!*************************************************************************
!
  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer:: n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer:: i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline


!end module METS2AORSA_MYRA


