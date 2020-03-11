c
c*******************************************************************************
c



      subroutine z_table(sgn_kprl, zeta, gamma, xnu, z0, z1, z2)

*     ----------------------------------------------------
*     z_table is a wrapper routine for David Smithe's table
*     lookup Z function
*     ----------------------------------------------------

      implicit none

      complex, intent(out)   :: z0, z1, z2
      complex, intent(in)    :: zeta
      real, intent(in)       :: sgn_kprl, gamma, xnu
      complex, dimension(0:4)  :: zfunc

      real alpha_s
      complex zi

      zi = cmplx(0.0, 1.0)

      alpha_s = -2.0 * gamma


      if (sgn_kprl .ge. 0.0) then

         call ztable_zfunc( real(zeta), alpha_s, xnu, zfunc )
         z0 = zfunc(0)
         z1 = - 0.5 * zfunc(1) - gamma / (2.0 * zi) * zfunc(2)
         z2 = 0.5 * zfunc(0) + 0.25 * zfunc(2)
     .      + gamma / (2.0 * zi) * zfunc(3)
     .      - gamma**2 / 4. * zfunc(4)


      else if (sgn_kprl .lt. 0.0) then

         call ztable_zfunc( - real(zeta), - alpha_s, - xnu, zfunc )
         z0 = - zfunc(0)
         z1 = - 0.5 * zfunc(1) + gamma / (2.0 * zi) * zfunc(2)
         z2 = - 0.5 * zfunc(0) - 0.25 * zfunc(2)
     .      + gamma / (2.0 * zi) * zfunc(3)
     .      + gamma**2 / 4. * zfunc(4)


      end if

      return

      end subroutine z_table



c
c*******************************************************************************
c

      subroutine z_exact(sgn_kprl, zeta, gamma, xnu, z0, z1, z2,
     .   dz0, dz1, dz2)

*     ----------------------------------------------------
*     z_exact is a wrapper routine for David Smithe's own
*     version of "z_smithe"
*     ----------------------------------------------------

      implicit none

      complex, intent(out)   :: z0, z1, z2, dz0, dz1, dz2
      complex, intent(in)    :: zeta
      real, intent(in)       :: sgn_kprl, gamma, xnu
      complex, dimension(0:4)  :: zfunc

      real alpha_s
      complex zi

      zi = cmplx(0.0, 1.0)

      alpha_s = -2.0 * gamma

      if ( sgn_kprl .ge. 0.0 ) then

         call ztable_zfunc( real(zeta), alpha_s, xnu, zfunc )
         z0 = zfunc(0)
         z1 = - 0.5 * zfunc(1) - gamma / (2.0 * zi) * zfunc(2)
         z2 = 0.5 * zfunc(0) + 0.25 * zfunc(2)
     .      + gamma / (2.0 * zi) * zfunc(3)
     .      - gamma**2 / 4. * zfunc(4)

         dz0 = zfunc(1)
         dz1 = - 0.5 * zfunc(2) - gamma / (2.0 * zi) * zfunc(3)
c         dz2 = 0.5 * zfunc(1) + 0.25 * zfunc(3)
c     .      + gamma / (2.0 * zi) * zfunc(4)



      else if ( sgn_kprl .lt. 0.0 ) then

         call ztable_zfunc( - real(zeta), - alpha_s, - xnu, zfunc )
         z0 = - zfunc(0)
         z1 = - 0.5 * zfunc(1) + gamma / (2.0 * zi) * zfunc(2)
         z2 = - 0.5 * zfunc(0) - 0.25 * zfunc(2)
     .      + gamma / (2.0 * zi) * zfunc(3)
     .      + gamma**2 / 4. * zfunc(4)

         dz0 =  zfunc(1)
         dz1 =  0.5 * zfunc(2) - gamma / (2.0 * zi) * zfunc(3)
c         dz2 =  0.5 * zfunc(1) + 0.25 * zfunc(3)
c     .      - gamma / (2.0 * zi) * zfunc(4)

      end if

      return

      end subroutine z_exact


c
c***************************************************************************
c


      subroutine z_approx(sgn_kprl, zeta, gamma, z0, z1, z2)

      implicit none

*     ------------------
*     Finite k_parallel:
*     ------------------
      real gamma, fgam, y0, y, sgn_kprl, descrim
      complex zeta, z0, z1, z2, zfunct, zetat, fzeta


      y0 = 1.5
      y = y0


      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gamma .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gamma * y) - 1.)
     .         / (2. * gamma * y)
         endif

         zetat = fgam * zeta
c        zfunct = fzeta(zetat)
         call zfun (zetat, zfunct)
         z0 = fgam * zfunct
         z1 = fgam * (1.0 + fgam * zeta * zfunct)
         z2 =  fgam**2 * zeta * (1.0 + fgam * zeta * zfunct)
      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gamma .gt. 1.0e-05)then
            descrim = 1. - 4. * gamma * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gamma * y) )
     .         / (2. * gamma * y)
         endif


         zetat = - fgam * zeta
c         zfunct = fzeta( zetat)
         call zfun (zetat, zfunct)
         z0 = - fgam * zfunct
         z1 = y / abs(y) * fgam * (1.0 - fgam * zeta * zfunct)
         z2 =  fgam**2 * zeta * (1.0 - fgam * zeta * zfunct)
      end if


      return
      end

c
c******************************************************************************
c

      subroutine z_approx0(sgn_kprl, alpha, xkprl_eff, zeta_eff,
     .   al, bl, cl)

      implicit none

*     -----------------------
*     Small k_parallel limit:
*     -----------------------
      real sgn_kprl, alpha, xkprl_eff
      complex zeta_eff, al, bl, cl, zfunct, zetat, fzeta

      if(sgn_kprl .ge. 0.0)then
c        zfunct = fzeta(zeta_eff)
         call zfun (zeta_eff, zfunct)

         al = 1. /(xkprl_eff * alpha) * zfunct
         bl = 1. /(xkprl_eff * alpha) * (1. + zeta_eff * zfunct)
         cl = zeta_eff /(xkprl_eff * alpha) *(1. + zeta_eff * zfunct)
      end if


      if(sgn_kprl .lt. 0.0)then
c        zfunct = fzeta(-zeta_eff)
         call zfun (-zeta_eff, zfunct)

         al = - 1. /(xkprl_eff * alpha) * zfunct
         bl = - 1. /(xkprl_eff * alpha) * (1. - zeta_eff * zfunct)
         cl = zeta_eff /(xkprl_eff * alpha) * (1. - zeta_eff * zfunct)
            end if


      return
      end


c
c***************************************************************************
c
      subroutine z2pole(zeta, zfunct)

      complex zeta, zfunct, zi, znum
      data pi/3.141592654/
      data zi/(0.0,1.0)/
      data g/1.141592654/
      data sqpi/1.772453851/
      znum = zi*sqpi + g*zeta
      zfunct = znum/(1.0 - zeta*znum)
      return
      end
c
c******************************************************************************
c

      subroutine z_smithe(sgn_kprl, zeta, gamma, z0, z1, z2)

      implicit none

      real dx, xmax, xnueff, gamma, sgn_kprl, xmax0
      integer j, npts, nmaxdim
      real x(10000), zeta_mod
      complex z0, z1, z2, zi, f(10000), zeta
      complex cexpx(10000)

      parameter (nmaxdim = 10000)

      zi = cmplx(0.,1.)
      xnueff = 0.0
      zeta_mod = sqrt(conjg(zeta) * zeta)


      if(zeta_mod .gt. 1.0e+02)then

         z0 = - 1.0 / zeta - 0.5 / zeta**3
     .                          + 3.0 * gamma / (zi * zeta**4)
         z1 = - 0.5 / zeta**2 + gamma / (zi * zeta**3)
     .                                       - 0.75 / zeta**4
         z2 = - 0.5 / zeta - 0.75 / zeta**3
     .                          + 4.5 * gamma / (zi * zeta**4)
         return

      else

         xmax0 = 7.0
         npts = 10000
         xmax = sgn_kprl * xmax0
         dx = xmax / (npts - 1)


*        ------------
*        calculate Z0
*        ------------
         do  j = 1, npts

            x(j) = (j-1) * dx
            cexpx(j) = exp(zi * zeta * x(j)
     .         - x(j)**2 / 4. * (1. + gamma * x(j))**2
     .         - xnueff / 8. * x(j)**3   )

            f(j) = zi * cexpx(j)
         end do

         call dgrat1(x, f, 1, npts, z0, nmaxdim)


*        ------------
*        calculate Z1
*        ------------
         do  j = 1, npts
            f(j) = 0.5 * x(j)* (1.0 + gamma * x(j)) * cexpx(j)
         end do

         call dgrat1(x, f, 1, npts, z1, nmaxdim)


*        ------------
*        calculate Z2
*        ------------
         do  j = 1, npts
            f(j) = zi / 2.0 * (1.0 -  x(j)**2 / 2.
     .           * (1. + gamma * x(j))**2) * cexpx(j)
         end do

         call dgrat1(x, f, 1, npts, z2, nmaxdim)



         return

      end if

      end

c
c*******************************************************************************
c

      subroutine dgrat1(x, f, nx1, nx2, ans, nmaxdim)

      implicit none

      real x(nmaxdim), dx
      complex f(nmaxdim), ans
      integer n, nx1, nx2, nmaxdim

      ans = 0.0

      do n = nx1, nx2 - 1
         dx = x(n+1) - x(n)
         ans = ans + dx * (f(n) + f(n+1)) / 2.0
      end do

      return
      end


c
c*******************************************************************************
c

      SUBROUTINE ZFUN(Z,FU)
C
C     ROUTINE WHICH EVALUATES THE PLASMA DISPERSION FUNCTION.  USES
C     NUMERICAL INTEGRATION (ABSOLUTE VALUE OF THE COMPLEX ARGUMENT, Z,
C     LESS THAN 5) OR ASYMPTOTIC EXPANSION.
C
      COMPLEX Z,FU,TEMP1,TEMP2,Z2,TPIIOD
      DIMENSION C(21),D(21),E(21),F(21),W(21)
      DATA DELTA/5.0E-1/, YI/-1.0E+0/, N/21/, ITEST/0/
      DATA ZERO/0.0E+0/, OSQPI/5.6418958355E-1/, TPI/6.28318530718E+0/
C
      IF(ITEST .EQ. 1) GO TO 1
      ITEST=1
C***     DEFINE WEIGHTS AND STORE CONSTANTS USED FOR INTEGRATION.
C***     WEIGHTS ARE DERIVED FROM A 3 POINT INTEGRATION SCHEME.
      NM3=N-3
      CONOI=DELTA*OSQPI
      W(1)=CONOI*3.75E-1
      W(2)=CONOI*1.1666666667E+0
      W(3)=CONOI*9.5833333333E-1
      DO 20 I=1,3
        NMIP1=N-I+1
   20 W(NMIP1)=W(I)
      DO 30 I=4,NM3
   30 W(I)=CONOI
      NO2=N/2
      NO2P1=NO2+1
      X=ZERO
      Y=YI
      E(NO2P1)=X
      F(NO2P1)=Y
      TEMP1=CMPLX(X,Y)
      TEMP2=CEXP(-TEMP1*TEMP1)
      C(NO2P1)=REAL(TEMP2)*W(NO2P1)
      D(NO2P1)=AIMAG(TEMP2)*W(NO2P1)
      DO 200 I=1,NO2
      X=DELTA*FLOAT(I)
      NPI=NO2P1+I
      NMI=NO2P1-I
      TEMP1=CMPLX(X,Y)
      TEMP2=CEXP(-TEMP1*TEMP1)
      C(NPI)=REAL(TEMP2)*W(NPI)
      C(NMI)=REAL(TEMP2)*W(NMI)
      D(NPI)=AIMAG(TEMP2)*W(NPI)
      D(NMI)=AIMAG(TEMP2)*W(NMI)*(-1.0E+0)
      E(NMI)=-X
      E(NPI)=X
      F(NPI)=Y
      F(NMI)=Y
  200 CONTINUE
      TPIOD=TPI/DELTA
      TPIIOD=CMPLX(ZERO,TPIOD)
C***     BEGIN CALCULATIONS.
    1 G=REAL(Z)
      YY=AIMAG(Z)
      H=ABS(YY)
      ZA=Z*CONJG(Z)
      IF(ZA .GE. 2.5E+1) GO TO 5
      Z2=Z*Z
C***     NUMERICAL INTEGRATION.
C***     F=1/SQRT(PI)*SUM OF...W(I)*EXP(-X(I)**2)/(X(I)-Z)...I=1,N.
C***     INTEGRATION IS ALONG A LINE X(I) IN THE COMPLEX PLANE, WHERE
C***     THE IMAGINARY PART OF X(I)=YI AND THE DIFFERENCE BETWEEN
C***     SUCCESSIVE REAL PARTS OF X(I)=DELTA.  LIMITS OF INTEGRATION
C***     ARE FROM -DELTA*N/2 TO DELTA*N/2.
C***     COMPUTE THE INTEGRAL BY TAKING THE SUM FROM 1 TO N OF THE
C***     CONSTANTS DIVIDED BY X(I)-Z.  USES REAL ARITHMETIC.
      ZR=0.0E+0
      ZI=0.0E+0
      DO 7 I=1,N
      A=E(I)-G
      B=F(I)-H
      DEN=A*A+B*B
      ODEN=1.0E+0/DEN
      ZR=ZR+(A*C(I)+B*D(I))*ODEN
      ZI=ZI+(A*D(I)-B*C(I))*ODEN
    7 CONTINUE
C***     ADD THE CORRECTION TERM.
      FU=CMPLX(ZR,ZI)+(0.0E+0,-3.5449077018E+0)*CEXP(-Z2-TPIOD*
     *   (H-YI)+TPIIOD*G)
      IF(YY .GE. ZERO) GO TO 6
C***     IMAGINARY PART OF ARGUMENT IS NEGATIVE.
      FU=CONJG(FU)+(0.0E+0,3.5449077018E+0)*CEXP(-Z2)
      GO TO 6
C***     MAGNITUDE OF ARGUMENT IS GREATER THAN 5, USE
C***     ASYMPTOTIC EXPANSION.
    5 CALL AEXPAN(Z,G,H,YY,FU)
    6 RETURN
      END
      SUBROUTINE AEXPAN(Z,G,H,YY,FU)
C
C     ROUTINE WHICH COMPUTES THE PLASMA DISPERSION FUNCTION USING
C     ASYMPTOTIC EXPANSION.  IF THE IMAGINARY PART OF THE ARGUMENT, YY,
C     IS EQUAL TO ZERO REAL ARITHMETIC IS USED.
C
      COMPLEX Z,FU,A,Z2,OZ2
      DATA ZERO/0.0E+0/, N/8/
C
      IF(YY .EQ. ZERO) GO TO 1
C***     COMPLEX ARITHMETIC.
      Z2=Z*Z
      OZ2=1.0E+0/Z2
      FU=-1.0E+0/Z
      A=FU
      EN=5.0E-1
      DO 10 I=1,N
      A=EN*A*OZ2
      FU=FU+A
   10 EN=EN+1.0E+0
      IF(YY .GT. ZERO) GO TO 30
      IF(H .GT. SQRT(G*G+1.72E+2)) GO TO 20
      FU=FU+(0.0E+0,3.5449077018E+0)*CEXP(-Z2)
      GO TO 30
C***     ERROR STOP TO AVOID OVERFLOW.
   20 WRITE (51,100) Z
  100 FORMAT(//1X,'*** ERROR STOP IN FRDCNT ROUTINE, ARGUMENT IS',
     * ' TOO SMALL,  ARG =',1P2E14.6)
      STOP
C***     REAL ARITHMETIC.
    1 X2=G*G
      OX2=1.0E+0/X2
      F=-1.0E+0/G
      B=F
      EN=5.0E-1
      DO 11 I=1,N
      B=EN*B*OX2
      F=F+B
   11 EN=EN+1.0E+0
      C=1.7724538509E+0*EXP(-X2)
      FU=CMPLX(F,C)
   30 RETURN
      END
C
C
C TWO POL APPROXIMATION
C
      SUBROUTINE Z2P(Z,F)
      COMPLEX Z,F
      F=CMPLX(.5,.81)/(CMPLX(.51,-.81)-Z)-CMPLX(.5,-.81)
     1/(CMPLX(.51,.81)+Z)
      RETURN
      END

c
c THIS CODE IS SUPERIOR TO PDF AND PDFLL.  NO ERROR AT ZETA=(4.,0.1)
c  disp(z) is the plasma z function
