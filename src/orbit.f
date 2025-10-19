c JCW this file is a bit of a mess with legacy code and possible bounds errors
c*****************************************************************************
c
      subroutine field_line_trace(sgn_vprl_in, i, j, 
     .   nmodesx, nmodesy,
     .   rwleft, rwright, ytop, ybottom, xprime_eq, yprime_eq,
     .   rmaxis, zmaxis, b0, psio, psimag, psi_tor_max,     
     .   bxn_eq, byn_eq, bzn_eq, bmod_eq, capr, y, psi, rho_pol2d, 
     .   qsafety, dx_eq, dy_eq,
     .   bmod_mid, capr_bpol_mid2, capr_bpol_mid, rho_tor2d,
     .   i_psi, dldb_tot12, dldbavg, 
     .   n_prof_flux, rhomax, norb_dim, lbmax,
     .   len_x, phin_x, capr_x, capz_x, nlen_exit)
     
      use size_mod            
      
      implicit none

      external f2, error

      common/fcom/ bxn, byn, bzn, bmod, bratio, 
     .   dx, dy, rt, xwleft, sgn_vprl, modb, bratio_phi, 
     .   dxdphi, dydphi, caprx, fcount, nxdim, nydim, nnodex, nnodey

      common/spline_com/sigma, zbxn, zbyn, zbzn, zbmod, zbratio, 
     .   xprime, yprime

      common/errcom/eps, s_err(100), y_len(100), nmax
      
      real dx_eq, dy_eq
      integer jmid, idum
      integer nmax, mmax, fcount, nxdim, nydim, n_len, n_len_max

      integer i0, j0, i_err

      real s_err, y_len, sgn_vprl, modb, bratio_phi, sgn_vprl_in
      real xphi, yphi, phi, dy_len(100), psi_tor_max
      real len, phi0
      real h0, eps, delta_b, caprx, r_max, drg, dzg, drg32, dzg32
      real x_extint, y_extint, xprimex, yprimex,  
     .   xprimex0, yprimex0
      integer norb_dim, nphi_enter, nlen_exit, nphii
      integer islpsw, islpsw1, ierr, nrho
      real sigma, rmaxis, zmaxis
          
      integer i_write, n_prof_flux, nphi1, nphi2, nuper, nupar
      
      CHARACTER(128) :: netCDF_file
     
      real capr_x(norb_dim), capz_x(norb_dim), phin_x(norb_dim),
     .   modb_x(norb_dim), dlen_x(norb_dim), dl_vprl(norb_dim),
     .   dl_bratio(norb_dim),
     .   dl_vprlh, len_x(norb_dim), length, length_prev
      real capr_x0, capz_x0
      real delta_x, delta_y, delta_z, delta_l, xprime_prev, delta_phi, 
     .   yprime_prev, phi_prev
      real yprime_want, dxdphi, dydphi, dxdphi_prev, dydphi_prev


      integer npts, ndim, ndeg, lxdata, iout, lipwr, ndval, lenws
      integer nr, ntheta, nrmax, nthmax, neqdsk, meqdsk, nk, mk
      integer ma, mr, mz, ipsi, iflag_gammab, isolve

      integer nxmx, nymx, idiag, jdiag, ieq
      integer ndfmax

      integer mkdim1, mkdim2

      integer nkdim1, nkdim2, nkx1, nkx2, nldim, nldim3,
     .   nky1, nky2, iant, jant1, jant2

      integer nxeqdmax, nyeqdmax

      real rmin, rmax, zmin, zmax, psio, ro, zo

      integer nrhomax

      integer n_theta_max, n_u_max, n_psi_max
     
      parameter (n_theta_max = 200)
      parameter (n_u_max = 200)
      parameter (n_psi_max = 200)

      parameter (nxmx = nmodesmax)
      parameter (nymx = mmodesmax)
      parameter (nrhomax = nmodesmax * 2)

      parameter (nkdim1 = - nmodesmax / 2)
      parameter (nkdim2 =   nmodesmax / 2)

      parameter (mkdim1 = - mmodesmax / 2)
      parameter (mkdim2 =   mmodesmax / 2)

      parameter (nldim  = nxmx * nymx)
      parameter (nldim3 = 3 * nldim)
      parameter (ndfmax = nldim3)

      parameter (nrmax = 201)
      parameter (nthmax = 64)

      parameter (nxeqdmax = 257, nyeqdmax = 257)

      parameter (lxdata = nxeqdmax * nyeqdmax)
      parameter (lipwr = 10)
      parameter (ndim = 2)
      parameter (ndeg = 2)


c*** ceez.f arrays:
      real zbxn(nxmx, nymx, 3)
      real zbyn(nxmx, nymx, 3)
      real zbzn(nxmx, nymx, 3)
      real zbmod(nxmx, nymx, 3)
      real zbratio(nxmx, nymx, 3)

c      real temp(2 *(nxmx + nymx) )
c      real surf2, curv2

      integer n_theta_(n_psi_max)

      real dtau(n_theta_max), dtau_tot(n_theta_max), 
     .   dtau_first(n_theta_max)
     
      real dtau_first_12, dtau_tot12
      real dldb_tot12(nxmx, nymx)
     
      real dtau_tot1(n_theta_max), dtau_first1(n_theta_max)
      real dtau_tot2(n_theta_max), dtau_first2(n_theta_max)

      integer n_theta, n_u, n_psi
      integer i_theta, i_u, i_psi
      real vc, vc_mks, sinthi, modbi, argi, vprl, dtau_sum, modbh,
     .   dtau_tot_sum, dldb_tot_sum, dldb_tot1, dldb_tot2

      real rho_eqdsk(nrmax),
     .     theta_eqdsk(nthmax)

      real psisep, psimag, router, z0
      integer nxeqd, nyeqd

      real dzdrhok, dzdthk, drdrhok, drdthk, xjacob

      real rhomax, lbmax, len1, len2

      integer nmodesx, nmodesy, nwdot, lmax, ibessel,
     .    inu, iprint, iexact,
     .    iroot, iequat, igeom,
     .    iqx, iqprof, iez, icurve, izfunc,
     .    nstep, nabs,
     .    isigma, itemp,
     .    nfreqm,  nkzm,
     .    idens,  ibackground, iabsorb,
     .    nzfun, nnodex, nnodey, i, j,
     .    jequat, iflag, liw, lw, nrhs, icenter, nboundary

      integer nnoderho

      integer n, m, nphi

      real ti0, xnuead, xnu1ad, xnu2ad, xant, te0,
     .    delta0, xwall, xnwall, delta,
     .    epszet, amu1, amu2, z1, z2, eta,
     .    b0, rt, ytop, ybottom, xnurf, aplasm, xnlim,
     .    xn0, flat, b1rat, b2rat, curdnx, curdny, curdnz,
     .    xnuabs, xbnch, xleft, xright,
     .    telim, tilim, ti2lim, ti3lim, rhoplasm,
     .    alphan, alphate, alphati,
     .    dfreq,  dkz, reomg1, reomg2, reomg3,
     .    r0, xnudip, adip, efold,
     .    amu3, z3, eta3, xnu3ad,
     .    xdelta, wdelta, xdelt2, wdelt2, zeffcd,
     .    rzoom1, rzoom2, q0, prfin,
     .    alim, grad, qavg0, ymax,
     .    rhonorm, ekappa, xiota0, rholim, psilim, psilim_, yant,
     .    rwleft, rwright, xwleft, xwright, psi_lim,
     .    rwleft_auto, rwright_auto, ytop_auto, ybottom_auto
     

      real xkthrho, wphase, vsound, domgk, xkthdx, omgestar, rhoi10,
     .   v0i, vthi10, vthe, vthi, vphase, rnz, xn2, eta2,
     .   xk0, shearedge, eta1, xn1, xmi1, xmh, xme, qi1, xmi2,
     .   xmi3, t0i, t0i2, t0i3, t0e, q, teedge, clight, xmu0, xlnlam,
     .   omgci10, omgrf, xmax, qe,i3, qi2, pi, eps0, xn3, qi3,
     .   costh, sinth, radius, rnx,  rny, rnphi, ti02, ti03, twopi


      real xprimec(nxmx), caprc(nxmx), xcourse(nxmx), capr(nxmx),
     .   xprime(nxmx), x(nxmx), dx, dxc, xprime_eq(nxmx) 
     
      real rhon(nrhomax), wdoti1avg(nrhomax), wdoti2avg(nrhomax),
     .   wdoteavg(nrhomax), drho, dvol(nrhomax), fvol(nrhomax),
     .   capr_bpol_mid(nrhomax), bmod_midavg(nrhomax), dldbavg(nrhomax)
            
      real rho_tor2d(nxmx, nymx), rho_pol2d(nxmx, nymx)
      real capr_bpol_mid2(nxmx, nymx),  bmod_mid(nxmx, nymx)
      real bxn_eq(nxmx, nymx), byn_eq(nxmx, nymx), 
     .     bzn_eq(nxmx, nymx), bmod_eq(nxmx, nymx)
      real psi(nxmx, nymx), qsafety(nxmx, nymx) 
      real bxn(nxmx, nymx), byn(nxmx, nymx), 
     .   bzn(nxmx, nymx), bmod(nxmx, nymx), 
     .   bratio(nxmx, nymx)
     
      real rhome, rhomi1, rhomi2, rhomi3, prod      

      real pcedotj1, pcedotje, pcedotj2, pcedotjt, pcedotj3
      real pedotj1, pedotje, pedotj2, pedotjt, pedotj3

      real p, pi1, pi2, pit, pi3, pe

      real yprimec(nymx), ycourse(nymx),
     .     yprime(nymx), y(nymx), dy, dyc, yprime_eq(nymx)

      integer icnc
      logical ismine, ismine1, ismine2, ismine3

      integer lld,nrow_local,ncol_local

      character(4):: suffix

      nxdim = nxmx
      nydim = nymx

      r0 = rmaxis
      z0 = zmaxis
      rt = r0
      
      nxdim = nxmx
      nydim = nymx
      
      nkx2 = nmodesx / 2
      nkx1 = - nmodesx / 2 + 1
      nnodex = nmodesx
      nnoderho = nnodex / 2

      nky2 = nmodesy / 2
      nky1 = - nmodesy / 2 + 1
      nnodey = nmodesy
      
      pi = 3.141592654
      twopi = 2.0 * pi
      
*     --------------------------------------------------------
*     Set spline parameters and calculate spline coefficients:
*       sigma = 0.0 for tensor product cubic splines
*       sigma = 50 for bi-linear interpolation
*       documentation recommend sigma=1.0 as standard value
*      --------------------------------------------------------
      sigma = 1.0
      islpsw = 255
      islpsw1 = 3
      
      sgn_vprl = sgn_vprl_in      
      xprime = xprime_eq
      yprime = yprime_eq
      dx = dx_eq
      dy = dy_eq
                    
      capr_x0 = capr(i) 
      capz_x0 = y(j)
      

      bxn = bxn_eq 
      byn = byn_eq
      bzn = bzn_eq
      bmod = bmod_eq
                
      x_extint = capr_x0 - rt
      y_extint = capz_x0
      
      xwleft = rwleft - rt
      xwright = rwright - rt

      xprimex0 = x_extint - xwleft
      yprimex0 = y_extint - ybottom
      phi0 = 0.0
                                             
*     ---------------------
*     Set extint parameters
*     ---------------------         
      h0 = 1.0e-04
      nmax = 3
      
      mmax = 5
      eps = 1.0e-05
      
      do i_err = 1, nmax
         s_err(i_err) = .1
      end do

*     -------------------------------------------
*     Set initial conditions for field line trace
*     -------------------------------------------
      len1 = 0.0
      len2 =  lbmax / 2.
         
      len = 0.0
      y_len(1) = xprimex0
      y_len(2) = yprimex0
      y_len(3) = phi0                 
                    
      call f(len, y_len, dy_len)           

*     ---------------------------------------------
*     Do norb_dim steps of field line trace in len
*     ----------------------------------------------                    
      do n_len = 1, norb_dim        
         fcount = 0
                       
         call extint(nmax, len, y_len, f2, h0, mmax, error)              
            
         xprimex = y_len(1)
         yprimex = y_len(2)
         phi     = y_len(3)

         x_extint = xprimex + xwleft
         y_extint = yprimex + ybottom       
               
*        --------------------------------------------
*        Save arrays of R, Z, phi as funtions of length
*        --------------------------------------------
         len_x(n_len) = len
         phin_x(n_len) = phi
         capr_x(n_len) = x_extint + rt
         capz_x(n_len) = y_extint
         modb_x(n_len) = modb
                         
c        write(6, 2213) n_len, len_x(n_len), 
c     .     capr_x(n_len), capz_x(n_len), phin_x(n_len), 
c     .     fcount, h0
     
         h0 = 2.0e-02                                           

         if (y_extint .ge. ytop .or. y_extint .le. ybottom) go to 201

         if (capr_x(n_len) .le. rwleft .or. 
     .                          capr_x(n_len) .ge. rwright) go to 201

         if (abs(len) .ge. len2) go to 201

      end do       
*     -----------------------
*     End of field line trace
*     -----------------------

  201 continue
  
      if(n_len .gt. norb_dim) n_len = norb_dim
      nlen_exit = n_len
      len_x = len_x * sgn_vprl
      
                            
           
2213  format(1i10, 1p,4e12.4, 1i5, 1p,1e12.4)      
9318  format(a128)
      return
      end

c
c***************************************************************************
c
      subroutine f2(x, y, dy) !integration of d x/d ell
      use size_mod        

      implicit none

      integer fcount, nnodex, nnodey, nxmx, nymx 

      real x, len, y(3), dy(3), br, bz, bphi, modb, x_extint, sgn_vprl
      real bratio_phi
      real xprime, yprime, dxdphi, dydphi, dr, dz, rt, capr, xwleft
      real drdl, dzdl, dphidl
                  
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

      common/fcom/ bxn, byn, bzn, bmod, bratio,
     .   dr, dz, rt, xwleft, sgn_vprl, modb, bratio_phi, 
     .   dxdphi, dydphi, capr, fcount, nxmx, nymx, nnodex, nnodey

      common/spline_com/sigma, zbxn, zbyn, zbzn, zbmod, zbratio, 
     .   xprimea, yprimea

      len = x
      xprime = y(1)
      yprime = y(2)

      br = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bxn, nxmx, zbxn, sigma)

      bz = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   byn, nxmx, zbyn, sigma)

      bphi =surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bzn, nxmx, zbzn, sigma)

      x_extint = xprime + xwleft
      capr = x_extint + rt
      
      drdl = br 
      dzdl = bz
      dphidl = bphi / capr

      dy(1) = drdl * sgn_vprl
      dy(2) = dzdl * sgn_vprl
      dy(3) = dphidl * sgn_vprl

      fcount = fcount + 1

      return
      end
c
c***************************************************************************
c
      subroutine f(x, y, dy) !integration of d x/d phi
      use size_mod        

      use aorsa2din_mod, only: eqtype
      implicit none

      integer fcount, nnodex, nnodey, nxmx, nymx 
c      integer nmodesmax, mmodesmax

      real x, phi, y(3), dy(3), br, bz, bphi, modb, x_extint, sgn_vprl
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

      common/fcom/ bxn, byn, bzn, bmod, bratio, 
     .   dr, dz, rt, xwleft, sgn_vprl, modb, bratio_phi, 
     .   dxdphi, dydphi, capr, fcount, nxmx, nymx, nnodex, nnodey

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

      !Note bxn etc are Normalized to bmod
      br = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bxn, nxmx, zbxn, sigma)

      bz = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   byn, nxmx, zbyn, sigma)

      bphi =surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bzn, nxmx, zbzn, sigma)

      modb = surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bmod, nxmx, zbmod, sigma)
     
      bratio_phi=surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea,
     .   bratio, nxmx, zbratio, sigma)



      x_extint = xprime + xwleft
      capr = x_extint + rt

      if (eqtype /= 'tokamak' .and. bphi <= 1.0e-4) then !find dX/dell if bphi<0.01 modb
         dxdphi = sgn_vprl * br!Use f2 instead
         dydphi = sgn_vprl * bz
      else 
         dxdphi = sgn_vprl * capr * br / bphi
         dydphi = sgn_vprl * capr * bz / bphi
      end if

      dy(1) = dxdphi
      dy(2) = dydphi
      !dy(3) = dphi/dphi=1

      fcount = fcount + 1

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

         dtau(n_theta) = 0.0

         sinth2_i = b_i / modb_init * sinth2_init(n_theta, i_psi)
         costh2_i = 1.0 - sinth2_i
         tanth2_i = sinth2_i / costh2_i
         costh_i = sqrt(costh2_i)
         
         if(sinth2_i > 1.0) cycle

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


      end do

      return
      end


c
c***************************************************************************
c
!check limits here JCW
      SUBROUTINE EXTINT (NMAX, X, Y, F, H0, MMAX, ERROR)
      INTEGER NMAX, MMAX
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
C                  MMAX OF FEWER EXTRAPOLATIONS, X (AT EXIT) = X (AT
C                  ENTRY) + H0. ON ENTRY X MUST BE SET TO THE INITIAL
C                  VALUE OF THE INDEPENDENT VARIABLE FOR THE INTEGRATION
C                  STEP BEING CONTEMPLATED.
C
C     Y            THE DEPENDENT VARIABLES. EACH DEPENDENT VARIABLE Y(N)
C                  (FOR N = 1, 2, .., NMAX) IS INTEGRATED FROM X TO X+H0
C                  IF ADEQUATE CONVERGENCE, AS DEFINED BY THE ''ERROR''
C                  SUBROUTINE, OCCURS IN MMAX OR FEWER EXTRAPOLATIONS.
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
C     MMAX        CONTAINS THE NUMBER OF TIMES EXTRAPOLATION IS ATTEM-
C                  PTED BEFORE H0 IS HALVED. THIS VALUE WILL VARY WITH
C                  COMPUTER ROUND-OFF ERROR AND WITH THE TYPE AND NUMBER
C                  OF EQUATIONS BEING INTEGRATED. IN ALL CASES, HOWEVER,
C                  ONE SHOULD SPECIFY AN MMAX 'GQ' 2. THE STEPSIZE FOR
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
      REAL   STEPFC, X0, U, SUM, YM, BETA, H, DEN, SQRT2, YP
      INTEGER J, K, L, M, N, LMAX, KASIDE, PTS, MM, MMAXP
      INTEGER KMIN
      REAL   HM(11), S(11), P(11), YBAR(100,11), Y0(100)
      REAL   YNEW(100), YOLD(100), DY(100), DY0(100), YHOLD(7,10,100)


! INITIALIZE..100
      MMAXP = MMAX + 1
      SQRT2 =  SQRT(2.0)
      FINISH = .FALSE.
      LATERL = .FALSE.
      X0 = X
      DO N = 1,NMAX
        Y0(N) = Y(N)
      END DO
      LMAX = (MMAX + 1)/2 + 1
      CALL F (X0, Y, DY0)

! START A NEW LEVEL..200

  204 X = X0 + H0
      KMIN = 1
      STEPFC = 2.0**(MMAX/3.0+0.5)
      DO N = 1, NMAX
         CONV(N) = .FALSE.
      END DO

      IF (.NOT.LATERL .OR. (MMAX.LT.1)) GO TO 203
C     ELSE BEGIN SHIFTING OLD INFORMATION INTO POSITION..
      DO N = 1, NMAX
        Y(N) = YHOLD(2,1,N)
        YBAR(N,1) = Y(N)
        DO L = 1, LMAX
          MM2 = MMAX - 2
          MMIN = IABS(2*L-3)
          DO M = MMIN, MM2
             YHOLD(L,M,N) = YHOLD (L+1, M+2, N)
          END DO
        END DO
      END DO
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
        S(MM) = 1 -  EXP(-BETA*H*H)

        IF (PREVIN) GO TO 503
C       ELSE GENERATE THE M-TH MIDPOINT INTEGRAL..
          DO N = 1, NMAX
            YOLD(N) = Y0(N)
            YNEW(N) = Y0(N) + H*DY0(N)
          END DO
          CALL F (X0+H, YNEW, DY)
          PTS = (H0*1.000001)/H
          DO 405 K = 2, PTS
C         BEGIN THE MIDPOINT INTEGRATION..
            DO N = 1, NMAX
              U = YOLD(N) + 2*H*DY(N)
              YOLD(N) = YNEW (N)
              YNEW (N) = U
            END DO
            CALL F (X0 + K*H, YNEW, DY)
            IF ((K.NE.KASIDE).OR.(L.LT.2)) GO TO 405
C           ELSE BEGIN PUTTING INFORMATION ASIDE..
              DO N = 1, NMAX
                YHOLD (L, M, N) = (YNEW(N) + YOLD(N) + H*DY(N))/2.0
              END DO
              L = L - 1
              KASIDE = 2*KASIDE
C           END OF PUTTING INFORMATION ASIDE
  405     CONTINUE
C         END OF THE MIDPOINT INTEGRATION

C
C NOW ADVANCE THE DEPENDENT VARIABLES..500
  503   IF (M.GT.0) GO TO 504
        DO N = 1, NMAX
          YBAR(N,1) = (YNEW(N) + YOLD(N) + H*DY(N))/2
          Y(N) = YBAR(N,1)
        END DO
        IF (MMAX.EQ.0) GO TO 700
        GO TO 400

C       NOW DETERMINE THE INTERPOLATIONAL POLYNOMIALS..
  504   IF (STEPFC.LT.1.1) KMIN = KMIN + 1
        DEN = 1
        DO K = KMIN, M
          P(K) = ((H/HM(K))**2)
          DO J = KMIN, M
            IF (J.NE.K) P(K) = P(K)*(S(J) - S(MM))/(S(J) - S(K))
          END DO
          DEN = DEN - P(K)
        END DO

C       END DETERMINATION OF THE INTERPOLATIONAL POLYNOMIALS
        DO 500 N=1, NMAX
          IF (CONV(N)) GO TO 500
          YP = Y(N)
          YM = YHOLD (1,M,N)
          IF (.NOT.PREVIN) YM = (YNEW(N) + YOLD(N) + H*DY(N))/2.0
          SUM = 0.0
          IF (M.LT.2) GO TO 501
          DO J = KMIN, M
            SUM = SUM + (YBAR(N,J) - YP)*P(J)
          END DO
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

C RETURN..700
  700 H0 = H0 * STEPFC
      RETURN
      END

c
c***************************************************************************
c
c
c***************************************************************************
c
!check limits here JCW
      SUBROUTINE EXTINT_F90 (NMAX, X, Y, F, H0, MMAX, ERROR)
C     IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      INTEGER NMAX, MMAX
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
C                  MMAX OF FEWER EXTRAPOLATIONS, X (AT EXIT) = X (AT
C                  ENTRY) + H0. ON ENTRY X MUST BE SET TO THE INITIAL
C                  VALUE OF THE INDEPENDENT VARIABLE FOR THE INTEGRATION
C                  STEP BEING CONTEMPLATED.
C
C     Y            THE DEPENDENT VARIABLES. EACH DEPENDENT VARIABLE Y(N)
C                  (FOR N = 1, 2, .., NMAX) IS INTEGRATED FROM X TO X+H0
C                  IF ADEQUATE CONVERGENCE, AS DEFINED BY THE ''ERROR''
C                  SUBROUTINE, OCCURS IN MMAX OR FEWER EXTRAPOLATIONS.
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
C     MMAX        CONTAINS THE NUMBER OF TIMES EXTRAPOLATION IS ATTEM-
C                  PTED BEFORE H0 IS HALVED. THIS VALUE WILL VARY WITH
C                  COMPUTER ROUND-OFF ERROR AND WITH THE TYPE AND NUMBER
C                  OF EQUATIONS BEING INTEGRATED. IN ALL CASES, HOWEVER,
C                  ONE SHOULD SPECIFY AN MMAX 'GQ' 2. THE STEPSIZE FOR
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
      INTEGER J, K, L, M, N, LMAX, KASIDE, PTS, MM, MMAXP
      INTEGER KMIN, icount
C     REAL*8 HM(11), S(11), P(11), YBAR(100,11), Y0(100)
      REAL   HM(11), S(11), P(11), YBAR(100,11), Y0(100)
C     REAL*8 YNEW(100), YOLD(100), DY(100), DY0(100), YHOLD(7,10,100)
      REAL   YNEW(100), YOLD(100), DY(100), DY0(100), YHOLD(7,10,100)
C
C
C INITIALIZE..100
      icount = 0
      MMAXP = MMAX + 1
      SQRT2 =  SQRT(2.0)
      FINISH = .FALSE.
      LATERL = .FALSE.
      X0 = X
      DO N = 1,NMAX
         Y0(N) = Y(N)
      END DO
      LMAX = (MMAX + 1)/2 + 1
      CALL F (X0, Y, DY0)

C START A NEW LEVEL..200
  204 X = X0 + H0
      icount = icount + 1
      !write(*,*) 'extint icount',icount
      KMIN = 1
      STEPFC = 2.0**(MMAX/3.0+0.5)
      DO  N = 1, NMAX
         CONV(N) = .FALSE.
      END DO
      
      IF (.NOT.LATERL .OR. (MMAX.LT.1)) GO TO 203

C     ELSE BEGIN SHIFTING OLD INFORMATION INTO POSITION..
      DO N = 1, NMAX
         Y(N) = YHOLD(2,1,N)
         YBAR(N,1) = Y(N)
         DO L = 1, LMAX
            MM2 = MMAX - 2
            MMIN = IABS(2*L-3)
            DO M = MMIN, MM2
               YHOLD(L,M,N) = YHOLD (L+1, M+2, N)
            END DO
         END DO
      END DO
C
C COMPUTING BETA AND ASSOCIATED QUANTITIES..300
  203 H = H0/2.0
      HM(1) = H0/2.0
      HM(2) = H0/4.0
      HM(3) = H0/6.0
      BETA = 0.25/(H0*H0)

C EXTRAPOLATIONS OF HIGHER ORDER EXPANSIONS..400
      DO MM = 1, MMAXP
         
C     BEGINNING THE LOOP OVER MMAX EXTRAPOLATIONS IN THIS LEVEL..
         M = MM - 1
         STEPFC = STEPFC/SQRT2
         PREVIN = .FALSE.
         IF (LATERL.AND.(M.LT.MMAX - 1)) PREVIN = .TRUE.
         KASIDE = 2
         IF (2*(M/2).EQ.M) KASIDE = 3
         L = (M + 1)/2 + 1
         IF (M.GT.2) HM(MM) = HM(MM-2)/2.
         H = HM(MM)
         S(MM) = 1 -  EXP(-BETA*H*H)
         IF (PREVIN) GO TO 503
C     ELSE GENERATE THE M-TH MIDPOINT INTEGRAL..
         DO N = 1, NMAX
            YOLD(N) = Y0(N)
            YNEW(N) = Y0(N) + H*DY0(N)
         END DO
         CALL F (X0+H, YNEW, DY)
         PTS = (H0*1.000001)/H
         DO K = 2, PTS
C     BEGIN THE MIDPOINT INTEGRATION..
            DO N = 1, NMAX
               U = YOLD(N) + 2*H*DY(N)
               YOLD(N) = YNEW (N)
               YNEW (N) = U
            END DO
                  
            CALL F (X0 + K*H, YNEW, DY)
            IF ((K.NE.KASIDE).OR.(L.LT.2)) CYCLE
C     ELSE BEGIN PUTTING INFORMATION ASIDE..
            DO N = 1, NMAX
               YHOLD (L, M, N) = (YNEW(N) + YOLD(N) + H*DY(N))/2.0
            END DO
            L = L - 1
            KASIDE = 2*KASIDE
C           END OF PUTTING INFORMATION ASIDE
         END DO
         
C        END OF THE MIDPOINT INTEGRATION

C     NOW ADVANCE THE DEPENDENT VARIABLES..500
         
 503     IF (M.GT.0) GO TO 504
         DO N = 1, NMAX
            YBAR(N,1) = (YNEW(N) + YOLD(N) + H*DY(N))/2.
            Y(N) = YBAR(N,1)
         END DO
         IF (MMAX.EQ.0) GO TO 700
         CYCLE
        
C     NOW DETERMINE THE INTERPOLATIONAL POLYNOMIALS..
 504     IF (STEPFC.LT.1.1) KMIN = KMIN + 1
         DEN = 1
         DO K = KMIN, M
            P(K) = ((H/HM(K))**2)
            DO J = KMIN, M
               IF (J.NE.K) P(K) = P(K)*(S(J) - S(MM))/(S(J) - S(K))
            END DO
            DEN = DEN - P(K)
         END DO
C     END DETERMINATION OF THE INTERPOLATIONAL POLYNOMIALS

         DO N=1, NMAX
            IF (CONV(N)) CYCLE
            YP = Y(N)
            YM = YHOLD (1,M,N)
            IF (.NOT.PREVIN) YM = (YNEW(N) + YOLD(N) + H*DY(N))/2.0
            SUM = 0.0
            IF (M.LT.2) GO TO 501
            DO  J = KMIN, M
               SUM = SUM + (YBAR(N,J) - YP)*P(J)
            END DO
 501        DY(N) = 0.0
            IF (DEN.NE.0.0) DY(N) = ((YM - YP) - SUM)/DEN
            Y(N) = YP + DY(N)
            YBAR(N,MM) = YM
         END DO
         CALL ERROR (M, DY, CONV, FINISH)
         IF (FINISH) GO TO 700
      END DO
C
C PREPARE FOR THE NEXT LEVEL IF NECESSARY..600
      LATERL = .TRUE.
      H0 = H0/2.0
      GO TO 204
C
C RETURN..700
 700  H0 = H0 * STEPFC
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

!      INTEGER M
!      REAL   DY(M)
!      LOGICAL CONV(M), FINISH
!      COMMON /ERRCOM/ EPS, S(100), Y(100), NMAX
!      REAL   EPS, S, Y
!      INTEGER NMAX, NTIMES (100)
      IMPLICIT NONE
      COMMON /ERRCOM/ EPS, S(100), Y(100), NMAX
      REAL   EPS, S, Y
      INTEGER NMAX

      INTEGER M
      REAL   DY(NMAX)  !JCW bug fix was size M
      LOGICAL CONV(NMAX), FINISH

      INTEGER NTIMES(100), NCONV, N


      IF (M==1) THEN
         NTIMES(1:NMAX) = 0
         NCONV = 0
      END IF
      
!      write(*,*) 'orbit err',nmax,m,dy,conv,s(1:nmax),eps

      DO N = 1, NMAX
        IF (.NOT.( ABS(DY(N))/S(N).LT.EPS).OR. CONV(N))  CYCLE
        NTIMES(N) = NTIMES(N) + 1
        IF (NTIMES(N).EQ. 1) NCONV = NCONV + 1
        IF (NTIMES(N).EQ. 2) CONV(N) = .TRUE.
        IF ( ABS(Y(N)).GT. S(N)) S(N) =  ABS(Y(N))
      END DO
      
      IF (NCONV.EQ.NMAX) FINISH = .TRUE.
      RETURN
      END

c
c***************************************************************************
c
