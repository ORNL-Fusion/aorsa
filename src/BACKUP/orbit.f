
c
c***************************************************************************
c
      subroutine f_new(x, y, dy)
      use size_mod        

      implicit none

      integer fcount, nnodex, nnodey, nxmx, nymx 
c      integer nmodesmax, mmodesmax

      real x, phi, y(1), dy(1), br, bz, bphi, modb, x_extint, sgn_vprl
      real bratio_phi
      real xprime, yprime, dxdphi, dydphi, dr, dz, rt, capr, xwleft
      real dxdl, dydl, dphidl
                  
      real bxn(nmodesmax, mmodesmax), byn(nmodesmax, mmodesmax), 
     .   bzn(nmodesmax, mmodesmax), bmod(nmodesmax, mmodesmax),
     .   bratio(nmodesmax, mmodesmax)

      real sigma, surf2
      real zbxn (nmodesmax, mmodesmax, 3)
      real zbyn (nmodesmax, mmodesmax, 3)
      real zbzn (nmodesmax, mmodesmax, 3)
      real zbmod(nmodesmax, mmodesmax, 3)
      real zbratio(nmodesmax, mmodesmax, 3)
      real xprimea(nmodesmax), yprimea(mmodesmax)

      common/fcom/fcount, bxn, byn, bzn, bmod, bratio, nxmx, nymx, 
     .   dr, dz,
     .   nnodex, nnodey, rt, xwleft, sgn_vprl, modb, bratio_phi, 
     .   dxdphi, dydphi, capr

      common/spline_com/sigma, zbxn, zbyn, zbzn, zbmod, zbratio, 
     .   xprimea, yprimea

      phi = x
      xprime = y(1)
      yprime = y(2)


      br = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bxn, nxmx, zbxn, sigma)

      bz = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   byn, nxmx, zbyn, sigma)

      bphi =surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bzn, nxmx, zbzn, sigma)

      modb =surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bmod, nxmx, zbmod, sigma)
     
      bratio_phi=surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bratio, nxmx, zbratio, sigma)


      x_extint = xprime + xwleft
      capr = x_extint + rt

      
      dxdl = br
      dydl = bz
      dphidl = bphi / capr


      dy(1) = dxdl
      dy(2) = dydl
      dy(3) = dphidl

      fcount = fcount + 1

 1312 format(1p8e12.4)

      return
      end

c
c***************************************************************************
c
      subroutine f(x, y, dy)
      use size_mod        

      implicit none

      integer fcount, nnodex, nnodey, nxmx, nymx 
c      integer nmodesmax, mmodesmax

      real x, phi, y(1), dy(1), br, bz, bphi, modb, x_extint, sgn_vprl
      real bratio_phi
      real xprime, yprime, dxdphi, dydphi, dr, dz, rt, capr, xwleft
      
*-------------------------------------------------------------
*     450 x 450 modes:
*     IMPORTANT!! The following dimensions must exactly match 
*     those in subroutine the main program in eqdsk_setup.f
*-------------------------------------------------------------
c      parameter (nmodesmax = 450)
c      parameter (mmodesmax = 450)
      
      
      real bxn(nmodesmax, mmodesmax), byn(nmodesmax, mmodesmax), 
     .   bzn(nmodesmax, mmodesmax), bmod(nmodesmax, mmodesmax),
     .   bratio(nmodesmax, mmodesmax)

      real sigma, surf2
      real zbxn (nmodesmax, mmodesmax, 3)
      real zbyn (nmodesmax, mmodesmax, 3)
      real zbzn (nmodesmax, mmodesmax, 3)
      real zbmod(nmodesmax, mmodesmax, 3)
      real zbratio(nmodesmax, mmodesmax, 3)
      real xprimea(nmodesmax), yprimea(mmodesmax)

      common/fcom/fcount, bxn, byn, bzn, bmod, bratio, nxmx, nymx, 
     .   dr, dz,
     .   nnodex, nnodey, rt, xwleft, sgn_vprl, modb, bratio_phi, 
     .   dxdphi, dydphi, capr

      common/spline_com/sigma, zbxn, zbyn, zbzn, zbmod, zbratio, 
     .   xprimea, yprimea

      phi = x
      xprime = y(1)
      yprime = y(2)

c      call intplt2(xprime, yprime, br,   nnodex, nnodey, bxn,
c     .   nxmx, nymx, dr, dz)
c      call intplt2(xprime, yprime, bz,   nnodex, nnodey, byn,
c     .   nxmx, nymx, dr, dz)
c      call intplt2(xprime, yprime, bphi, nnodex, nnodey, bzn,
c     .   nxmx, nymx, dr, dz)
c      call intplt2(xprime, yprime, modb, nnodex, nnodey, bmod,
c     .   nxmx, nymx, dr, dz)
c      call intplt2(xprime, yprime, bratio_phi, nnodex, nnodey, bratio,
c     .   nxmx, nymx, dr, dz)


      br = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bxn, nxmx, zbxn, sigma)

      bz = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   byn, nxmx, zbyn, sigma)

      bphi =surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bzn, nxmx, zbzn, sigma)

      modb =surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bmod, nxmx, zbmod, sigma)
     
      bratio_phi=surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bratio, nxmx, zbratio, sigma)



      x_extint = xprime + xwleft
      capr = x_extint + rt

      dxdphi = sgn_vprl * capr * br / bphi
      dydphi = sgn_vprl * capr * bz / bphi


      dy(1) = dxdphi
      dy(2) = dydphi

      fcount = fcount + 1

 1312 format(1p8e12.4)

      return
      end

c
c***************************************************************************
c


      subroutine fdtau(dtau, nxdim, nydim, len_x, modb_x,
     .               n_theta_max, n_psi_max, norb_dim, sinth2_init, 
     .               modb_init, n_theta_,
     .               i_psi, i, j, nphi_enter, nphi_exit, sgn_vprl)

      implicit none

      integer n_theta_max, n_psi_max, norb_dim, n_theta, i_psi
      integer nxdim, nydim

      real dtau(n_theta_max)
      real len_x(norb_dim), modb_x(norb_dim)

      real sinth2_init(n_theta_max, n_psi_max)

      integer n_theta_(n_psi_max), i, j, nphi_enter, nphi_exit

      real mri, b_ip1, b_i, dl, l_ip1, l_i, sgn_vprl, argi, modb_init
      real costh_i, tanth2_i, sinth2_i, costh2_i

*     -------------------------------------------------
*     This subroutine calculates the time (dtau)
*     that a particle stays in the (i, j) spatial cell
*     -------------------------------------------------

      l_i =   len_x(nphi_enter)
      l_ip1 = len_x(nphi_exit)

      b_i =   modb_x(nphi_enter)
      b_ip1 = modb_x(nphi_exit)

      mri = b_ip1 / b_i
      dl = l_ip1 - l_i
      dl = dl * sgn_vprl

      do n_theta = 1, n_theta_(i_psi)

c         thetai = theta_(n_theta, i_psi)
c         costhi = cos(thetai)
c         sinthi = sin(thetai)
c         tanthi = tan(thetai)

c         write(6, 100)modb_init
c  100    format(1p8e12.4)
c         call exit

         dtau(n_theta) = 0.0

         sinth2_i = b_i / modb_init * sinth2_init(n_theta, i_psi)
	 costh2_i = 1.0 - sinth2_i
	 tanth2_i = sinth2_i / costh2_i
	 costh_i = sqrt(costh2_i)
	 
	 if(sinth2_i .gt. 1.0) go to 300

         argi = 1.0 - (mri - 1.0) * tanth2_i


*        -------
*        Passing
*        -------
         if(argi .ge. 0.0) then          ! passing

            dtau(n_theta) = dl / sgn_vprl
     .            * 2.0 / abs(costh_i) / (1.0 + sqrt(argi))
         end if



*        -------
*        Trapped
*        -------
         if(argi .lt. 0.0) then          ! trapped

            dtau(n_theta) = dl / sgn_vprl
     .           * 2.0 * abs(costh_i)  / (mri - 1.0) / sinth2_i
         end if

  300    continue

      end do
      



      return
      end


c
c***************************************************************************
c



      SUBROUTINE EXTINT (NMAX, X, Y, F, H0, M MAX, ERROR)
C     IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      INTEGER N MAX, M MAX
C     REAL*8 X, H0, Y(NMAX)
      REAL   X, H0, Y(NMAX)
C
C          A DEFERRED-LIMIT INTEGRATOR (J.P. BORIS AND N.K. WINSOR)
C
C         THIS SUBROUTINE INTEGRATES UP TO 100 SIMULTANEOUS FIRST ORDER
C     ORDINARY DIFFERENTIAL EQUATIONS FROM X TO X + H0 BY REPEATED EX-
C     TRAPOLATIONS ON A MIDPOINT RULE. UP TO 10 EXTRAPOLATIONS MAY BE
C     REQUESTED BEFORE REDUCTION OF THE INITIAL STEPSIZE IS CARRIED OUT.
C     (REF. R. BULIRSCH AND J. STOER -  NUMERISCHE MATHEMATIK 8,1 (1966)
C
C     NMAX         THE TOTAL NUMBER OF DEPENDENT VARIABLES BEING INTE-
C                  GRATED. THERE WILL BE ONE FIRST ORDER EQUATION FOR
C                  EACH DEPENDENT VARIABLE.
C
C     X            THE INDEPENDENT VARIABLE, X IS TREATED AS REAL AND
C                  MONOTONIC DURING THE INTEGRATION STEP. THE VALUE OF
C                  X ON RETURN FROM ''EXT INT'' CONTAINS THE VALUE WHICH
C                  IS APPROPRIATE TO THE LENGTH OF THE INTEGRATION ACTU-
C                  ALLY PERFORMED. IF CONVERGENCE HAS BEEN OBSERVED IN
C                  M MAX OF FEWER EXTRAPOLATIONS, X (AT EXIT) = X (AT
C                  ENTRY) + H0. ON ENTRY X MUST BE SET TO THE INITIAL
C                  VALUE OF THE INDEPENDENT VARIABLE FOR THE INTEGRATION
C                  STEP BEING CONTEMPLATED.
C
C     Y            THE DEPENDENT VARIABLES. EACH DEPENDENT VARIABLE Y(N)
C                  (FOR N = 1, 2, .., NMAX) IS INTEGRATED FROM X TO X+H0
C                  IF ADEQUATE CONVERGENCE, AS DEFINED BY THE ''ERROR''
C                  SUBROUTINE, OCCURS IN M MAX OR FEWER EXTRAPOLATIONS.
C
C     F            THE DERIVATIVE SUBROUTINE SUPPLIED BY THE USER FOR
C                  HIS PARTICULAR PROBLEM. IT MUST BE OF THE FORM
C                  F (X, Y, DY) WHERE THE ARRAY DY(N) (FOR N = 1,2, ...,
C                  NMAX) IS RETURNED CONTAINING THE NMAX DERIVATIVES
C                  (DY(N)/DX)(X,Y).
C
C     H0           IS THE BASIC STEPSIZE OF THE INTEGRATION. IF CONVER-
C                  GENCE OCCURS WITHIN MMAX EXTRAPOLATIONS, X RETURNS
C                  FROM EXT INT WITH THE VALUE X + H0. THE VALUES IN THE
C                  ARRAY Y ARE THE VALUES OF THE DEPENDENT VARIABLES AT
C                  THIS VALUE OF X. IF CONVERGENCE DOES NOT OCCUR, H0 IS
C                  HALVED AND THE ENTIRE EXTRAPOLATION PROCEDURE IS
C                  REPEATED AND REPEATED AGAIN UNTIL CONVERGENCE OCCURS.
C                  AN ATTEMPT HAS BEEN MADE TO UTILIZE AS MUCH PREVIOUS-
C                  LY COMPUTED INFORMATION AS POSSIBLE WHEN H0 MUST BE
C                  HALVED FOR CONVERGENCE.
C
C     M MAX        CONTAINS THE NUMBER OF TIMES EXTRAPOLATION IS ATTEM-
C                  PTED BEFORE H0 IS HALVED. THIS VALUE WILL VARY WITH
C                  COMPUTER ROUND-OFF ERROR AND WITH THE TYPE AND NUMBER
C                  OF EQUATIONS BEING INTEGRATED. IN ALL CASES, HOWEVER,
C                  ONE SHOULD SPECIFY AN M MAX 'GQ' 2. THE STEPSIZE FOR
C                  FUTURE ITERATIONS IS SELECTED AND RETURNED IN H0.
C                  THIS NEW VALUE IS CHOSEN SO THAT CONVERGENCE WILL BE
C                  OBSERVED AT ABOUT THE 0.66*MMAX-TH EXTRAPOLATION.
C                  AS WRITTEN, MMAX.LE.10.
C
C     ERROR        IS A SUBROUTINE WHICH IS USED TO DETERMINE THE SATIS-
C                  FACTORY CONVERGENCE OF EACH INDIVIDUAL EQUATION BEING
C                  INTEGRATED. AN EXAMPLE IS GIVEN WHICH CORRESPONDS TO
C                  THE ERROR CRITERION USED BY BULIRSCH AND STOER.
C                      THE CALLING SEQUENCE IS
C                  ERROR (M, DY, CONV, FINISH) WHERE M IS THE ORDER OF
C                  EXTRAPOLATION, DY IS THE VECTOR OF INCREMENTS TO THE
C                  DEPENDENT VARIABLES Y, CONV IS A VECTOR OF LOGICALS
C                  WHICH ARE .TRUE. COMPONENTWISE WHEN THE CORRESPONDING
C                  DEPENDENT VARIABLE HAS CONVERGED, AND FINISH IS
C                  RETURNED .TRUE. WHEN THE ENTIRE SYSTEM OF EQUATIONS
C                  HAS SATISFIED THE CONVERGENCE CRITERION.
C
      LOGICAL LATERL, CONV(100), PREVIN, FINISH
C     REAL*8 STEPFC, X0, U, SUM, YM, BETA, H, DEN, SQRT2, YP
      REAL   STEPFC, X0, U, SUM, YM, BETA, H, DEN, SQRT2, YP
      INTEGER J, K, L, M, N, L MAX, K ASIDE, PTS, MM, MMAXP
      INTEGER KMIN
C     REAL*8 HM(11), S(11), P(11), YBAR(100,11), Y0(100)
      REAL   HM(11), S(11), P(11), YBAR(100,11), Y0(100)
C     REAL*8 Y NEW(100), Y OLD(100), DY(100), DY0(100), Y HOLD(7,10,100)
      REAL   Y NEW(100), Y OLD(100), DY(100), DY0(100), Y HOLD(7,10,100)
C
C
C INITIALIZE..100
      MMAXP = MMAX + 1
      SQRT2 =  SQRT(2.0)
      FINISH = .FALSE.
      LATERL = .FALSE.
      X0 = X
      DO 100 N = 1,NMAX
  100   Y0(N) = Y(N)
      LMAX = (MMAX + 1)/2 + 1
      CALL F (X0, Y, DY0)
C
C START A NEW LEVEL..200
  204 X = X0 + H0
      KMIN = 1
      STEPFC = 2.0**(MMAX/3.0+0.5)
      DO 205 N = 1, NMAX
  205   CONV(N) = .FALSE.
      IF (.NOT.LATERL .OR. (MMAX.LT.1)) GO TO 203
C     ELSE BEGIN SHIFTING OLD INFORMATION INTO POSITION..
      DO 202 N = 1, NMAX
        Y(N) = Y HOLD(2,1,N)
        YBAR(N,1) = Y(N)
        DO 202 L = 1, LMAX
          MM2 = MMAX - 2
          MMIN = IABS(2*L-3)
          DO 202 M = MMIN, MM2
  202       Y HOLD(L,M,N) = Y HOLD (L+1, M+2, N)
C
C COMPUTING BETA AND ASSOCIATED QUANTITIES..300
  203 H = H0/2
      HM(1) = H0/2.0
      HM(2) = H0/4.0
      HM(3) = H0/6.0
      BETA = 0.25/(H0*H0)
C
C EXTRAPOLATIONS OF HIGHER ORDER EXPANSIONS..400
      DO 400 MM = 1, MMAXP
C     BEGINNING THE LOOP OVER MMAX EXTRAPOLATIONS IN THIS LEVEL..
        M = MM - 1
        STEPFC = STEPFC/SQRT2
        PREVIN = .FALSE.
        IF (LATERL.AND.(M.LT.MMAX - 1)) PREVIN = .TRUE.
        KASIDE = 2
        IF (2*(M/2).EQ.M) KASIDE = 3
        L = (M + 1)/2 + 1
        IF (M.GT.2) HM(MM) = HM(MM-2)/2
        H = HM(MM)
C       S(MM) = 1 - DEXP(-BETA*H*H)
        S(MM) = 1 -  EXP(-BETA*H*H)
        IF (PREVIN) GO TO 503
C       ELSE GENERATE THE M-TH MIDPOINT INTEGRAL..
          DO 404 N = 1, NMAX
            Y OLD(N) = Y0(N)
  404       Y NEW(N) = Y0(N) + H*DY0(N)
          CALL F (X0+H, YNEW, DY)
          PTS = (H0*1.000001)/H
          DO 405 K = 2, PTS
C         BEGIN THE MIDPOINT INTEGRATION..
            DO 406 N = 1, NMAX
              U = Y OLD(N) + 2*H*DY(N)
              Y OLD(N) = Y NEW (N)
  406         Y NEW (N) = U
            CALL F (X0 + K*H, Y NEW, DY)
            IF ((K.NE.KASIDE).OR.(L.LT.2)) GO TO 405
C           ELSE BEGIN PUTTING INFORMATION ASIDE..
              DO 408 N = 1, NMAX
  408           Y HOLD (L, M, N) = (YNEW(N) + YOLD(N) + H*DY(N))/2.0
              L = L - 1
              KASIDE = 2*KASIDE
C           END OF PUTTING INFORMATION ASIDE
  405     CONTINUE
C         END OF THE MIDPOINT INTEGRATION
C
C NOW ADVANCE THE DEPENDENT VARIABLES..500
  503   IF (M.GT.0) GO TO 504
        DO 505 N = 1, NMAX
          YBAR(N,1) = (YNEW(N) + YOLD(N) + H*DY(N))/2
  505     Y(N) = YBAR(N,1)
        IF (MMAX.EQ.0) GO TO 700
        GO TO 400
C       NOW DETERMINE THE INTERPOLATIONAL POLYNOMIALS..
  504   IF (STEPFC.LT.1.1) KMIN = KMIN + 1
        DEN = 1
        DO 401 K = KMIN, M
          P(K) = ((H/HM(K))**2)
          DO 410 J = KMIN, M
            IF (J.NE.K) P(K) = P(K)*(S(J) - S(MM))/(S(J) - S(K))
  410     CONTINUE
  401     DEN = DEN - P(K)
C       END DETERMINATION OF THE INTERPOLATIONAL POLYNOMIALS
        DO 500 N=1, N MAX
          IF (CONV(N)) GO TO 500
          YP = Y(N)
          YM = Y HOLD (1,M,N)
          IF (.NOT.PREVIN) YM = (YNEW(N) + YOLD(N) + H*DY(N))/2.0
          SUM = 0.0
          IF (M.LT.2) GO TO 501
          DO 502 J = KMIN, M
  502       SUM = SUM + (YBAR(N,J) - YP)*P(J)
  501     DY(N) = 0.0
          IF (DEN.NE.0.0) DY(N) = ((YM - YP) - SUM)/DEN
          Y(N) = YP + DY(N)
          YBAR(N,MM) = YM
  500   CONTINUE
        CALL ERROR (M, DY, CONV, FINISH)
        IF (FINISH) GO TO 700
  400 CONTINUE
C
C PREPARE FOR THE NEXT LEVEL IF NECESSARY..600
      LATERL = .TRUE.
      H0 = H0/2
      GO TO 204
C
C RETURN..700
  700 H0 = H0 * STEPFC
      RETURN
      END

c
c***************************************************************************
c

      SUBROUTINE ERROR (M, DY, CONV, FINISH)
C
C         THIS VERSION OF THE ERROR ROUTINE CORRESPONDS CLOSELY TO THE
C     VERSION EMPLOYED BY BULIRSCHE AND STOER. THE DIFFERENCES, HOWEVER,
C     ARE NOTEWORTHY. THIS VERSION IS CALLED ONLY AFTER ALL UNCONVERGED
C     DEPENDENT VARIABLES HAVE BEEN EXTRAPOLATED AT ORDER M. THUS AN
C     EXCESSIVE NUMBER OF SUBROUTINE REFERENCES ARE MADE UNNECESSARY.
C         NEW VALUES FOR ANY PARTICULAR DEPENDENT VARIABLE ARE NOT CAL-
C     CULATED AFTER CONVERGENCE, AS DEFINED BY THIS SUBROUTINE, HAS BEEN
C     OBSERVED FOR THAT VARIABLE. THUS POSSIBLE INSTABILITIES ASSOCIATED
C     WITH ATTEMPTS AT OVERCONVERGENCE CAN BE ELIMINATED. ONE OF THESE
C     INSTABILITIES, ARISING FROM RESETTING THE MAGNITUDE VECTOR S AT
C     EVERY EXTRAPOLATION AS DOES THE B - S PROGRAM, IS ELIMINATED BY
C     RESETTING S ONLY AFTER THE CORRESPONDING DEPENDENT VARIABLE HAS
C     CONVERGED ACCORDING TO THE CRITERIA CHOSEN.
C     IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      INTEGER M
C     REAL*8 DY(M)
      REAL   DY(M)
      LOGICAL CONV(M), FINISH
      COMMON /ERRCOM/ EPS, S(100), Y(100), NMAX
C     REAL*8 EPS, S, Y
      REAL   EPS, S, Y
      INTEGER NMAX, NTIMES (100)
C
      IF (M.NE.1) GO TO 1
        DO 3 N = 1, NMAX
    3     NTIMES(N) = 0
        NCONV = 0
    1 DO 2 N = 1, NMAX
C       IF (.NOT.(DABS(DY(N))/S(N).LT.EPS).OR. CONV(N))  GO TO 2
        IF (.NOT.( ABS(DY(N))/S(N).LT.EPS).OR. CONV(N))  GO TO 2
        NTIMES(N) = NTIMES(N) + 1
        IF (NTIMES(N).EQ. 1) NCONV = NCONV + 1
        IF (NTIMES(N).EQ. 2) CONV(N) = .TRUE.
C       IF (DABS(Y(N)).GT. S(N)) S(N) = DABS(Y(N))
        IF ( ABS(Y(N)).GT. S(N)) S(N) =  ABS(Y(N))
    2 CONTINUE
      IF (NCONV.EQ.NMAX) FINISH = .TRUE.
      RETURN
      END

c
c***************************************************************************
c


