      subroutine run_rf2x(nmodesx, nmodesy, rwleft, rwright,
     &   ytop, ybottom, myid, nxdim, nydim,
     &   rt, b0, rho,
     &   redotje, redotji,
     &   redotj1, redotj2,
     &   redotj3, redotj4,
     &   redotj5, redotj6,
     &   xjprl,  wdote,
     &   wdoti1, wdoti2,
     &   wdoti3, wdoti4,
     &   wdoti5, wdoti6, rhomax)
      use size_mod 
      
      implicit none

      integer n, ndata, ndim, ndata2

      real  alphan_slo, alphan3, alphati3, alphati2, alphan2, xn3lim,
     &     xn2lim, xnslo, xnslolim, drho_data, rhot, rhomin

      integer nxdim, nydim
      integer nnode_local, nnode_overlap, ftrap, l, lmaxdim
      integer nnodex_loc, nx_overlap, nnodex_eff
      integer nnodey_loc, ny_overlap, nnodey_eff
      integer istart, ifinish, jstart, jfinish,
     &  i_seg, j_seg, iseg1, iseg2, jseg1, jseg2, i0, j0

      integer nxmx, nymx, idiag, jdiag, ieq, n_prof_flux
      integer ndfmax,
     &    ninteg, nd, izoom1, izoom2, jzoom1, jzoom2, ndf, nmaxe, irnc
      integer nrow, ncol, norder, iprofile

      integer mkdim1, mkdim2

      integer nkdim1, nkdim2, nkx1, nkx2, nldim, nldim3,
     &   nky1, nky2, iant

      real eslow, eslowev, qi_slo, eta_slo, z_slo, xmi_slo, amu_slo,
     &    xn_slo

      real flimiter, parabola, dfdx, dfdy

      real tmem, tsys, tio, ttotal, time0, time, cpu, dummy, second1
      real tmin, gflops, gflopsp, ops, teev
      real sqx, signb, gausspsi, dpsiant

      real dxsqx, dxxsqx, dxysqx, dysqx, dyysqx

      real dthetant0, dpsiant0, psiant, psipne, psipte,
     &   psipti1, psipti2, psipti3

c      integer nmodesmax, mmodesmax 
      integer nrhomax, isolve


      integer nnodex_fine, nnodey_fine

      parameter (nnodex_fine = 500)
      parameter (nnodey_fine = 500)
      

*     -------------------------------------------------------------------
*     450 x 450 modes:
*     IMPORTANT: These numer must match EXACTLY those in main program!!!!
*     -------------------------------------------------------------------
c      parameter (nmodesmax = 450)
c      parameter (mmodesmax = 450)
      
      

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

      parameter (lmaxdim = 99)


      real xnavg(nrhomax), xn1avg(nrhomax), xn2avg(nrhomax),
     &     xn3avg(nrhomax), xna_sloavg(nrhomax)
      real xkteavg(nrhomax), xktiavg(nrhomax),  xkti2avg(nrhomax),
     &     xkti3avg(nrhomax)


      complex zi

      real xnuomg, rhomax

      real gaussian_ne, gaussian_te, gaussian_ti1, gaussian_ti2,
     &   gaussian_ti3

      real frho, dxfrho, dxxfrho, dxyfrho, dyfrho, dyyfrho


      integer nmodesx, nmodesy, nwdot, lmax, ibessel,
     &    inu, iprint, iexact,
     &    iroot, iequat, igeom,
     &    iqx, iqprof, iez, icurve, izfunc,
     &    nstep, nabs,
     &    isigma, itemp,
     &    nfreqm,  nkzm,
     &    idens,  ibackground, iabsorb,
     &    nzfun, nnodex, nnodey, i, j,
     &    jequat, iflag, liw, lw, nrhs, icenter, nboundary

      integer nnoderho

      integer m, nphi, iflag_gammab

      real ti0, xnuead, xnu1ad, xnu2ad, xant, te0,
     &    delta0, xwall, xnwall, delta,
     &    epszet, amu1, amu2, z1, z2, eta,
     &    b0, rt, ytop, ybottom, xnurf, aplasm, xnlim, 
     &    xn0, flat, b1rat, b2rat, curdnx, curdny, curdnz,
     &    xnuabs, xbnch, xleft, xright,
     &    telim, tilim, ti2lim, ti3lim, rhoplasm,
     &    alphan, alphate, alphati, alphaq,
     &    dfreq,  dkz, reomg1, reomg2, reomg3,
     &    r0, xnudip, adip, efold,
     &    amu3, z3, eta3, xnu3ad,
     &    xdelta, wdelta, xdelt2, wdelt2, zeffcd,
     &    rzoom1, rzoom2, yzoom1, yzoom2, q0, prfin,
     &    alim, grad, qavg0, ymax,
     &    ekappa, xiota0, rholim, psilim, psimol, yant,
     &    rwleft, rwright, xwleft, xwright, psi_lim, psi1

      real xkphi(nxmx), xktau, xkrho, rant, xkphi0, xnphi
      real rlim

      real xkthrho, wphase, vsound, domgk, xkthdx, omgestar, rhoi10,
     &   v0i, vthi10, vthe, vthi, vphase, rhoi1overl, rnz, xn2, eta2,
     &   xk0, shearedge, eta1, xn1, xmi1, xmh, xme, qi1, xmi2,
     &   xmi3, t0i, t0i2, t0i3, t0e, q, teedge, clight, xmu0, xlnlam,
     &   omgci10, omgrf, xmax, qe,i3, qi2, pi, eps0, xn3, qi3,
     &   costh, sinth, radius, rnx,  rny, rnphi, ti02, ti03, twopi

      real xjantx, xjanty, xjantz, xjant
      real xjx(nxmx, nymx), xjy(nxmx, nymx), xjz(nxmx, nymx)
      real djxdx, djydy
      real xprimec(nxmx), caprc(nxmx), xcourse(nxmx), capr(nxmx),
     &   xprime(nxmx), x(nxmx), dx, dxc


      real, dimension(:),   allocatable :: x_fine
      real, dimension(:),   allocatable :: y_fine
      real, dimension(:),   allocatable :: capr_fine
      real, dimension(:,:), allocatable :: rho_fine
      real dx_fine, dy_fine

      real rhon(nrhomax),   wdoti1avg(nrhomax), wdoti2avg(nrhomax),
     &  wdoti3avg(nrhomax), wdoti4avg(nrhomax), wdoti5avg(nrhomax),
     &  wdoti6avg(nrhomax),
     &  wdoteavg(nrhomax),  wdotavg(nrhomax), drho
     
      real wdoti1_dvol(nrhomax), wdoti2_dvol(nrhomax),
     &     wdoti3_dvol(nrhomax), wdoti4_dvol(nrhomax), 
     &     wdoti5_dvol(nrhomax), wdoti6_dvol(nrhomax),
     &     wdote_dvol(nrhomax)
     
      real redotj1_dvol(nrhomax), redotj2_dvol(nrhomax),
     &     redotj3_dvol(nrhomax), redotj4_dvol(nrhomax),
     &     redotj5_dvol(nrhomax), redotj6_dvol(nrhomax),
     &     redotje_dvol(nrhomax)
     
      real capr_bpol_mid(nrhomax)
      real capr_bpol_mid2(nxmx, nymx)

      real fyavg(nrhomax), qhat, omgte, omgti,
     &     vthe0, vthi0, xnuee, xnuii, xnu7omg

      real redotj1avg(nrhomax), redotj2avg(nrhomax),
     &     redotjeavg(nrhomax), redotj3avg(nrhomax),
     &     redotj4avg(nrhomax),
     &     redotj5avg(nrhomax), redotj6avg(nrhomax),
     &     redotjsavg(nrhomax),
     &     redotjiavg(nrhomax), xjprlavg(nrhomax) 
      real darea(nrhomax), dvol(nrhomax), volume(nrhomax)

      real xjprl_int(nrhomax), fyp_int(nrhomax), redotje_int(nrhomax)
      
      real redotji1_int(nrhomax), redotji2_int(nrhomax),
     &     redotji3_int(nrhomax), redotji4_int(nrhomax), 
     &     redotji5_int(nrhomax), redotji6_int(nrhomax)
     
      real wdote_int(nrhomax),
     &     wdoti1_int(nrhomax), wdoti2_int(nrhomax), 
     &     wdoti3_int(nrhomax), wdoti4_int(nrhomax),
     &     wdoti5_int(nrhomax), wdoti6_int(nrhomax)
     
     

      real fypavg(nrhomax), fypi1avg(nrhomax), fypi2avg(nrhomax),
     &     fypeavg(nrhomax), fypi3avg(nrhomax)
      real vyavg(nrhomax),  vyi1avg(nrhomax),  vyi2avg(nrhomax)
      real xnupiavg(nrhomax), rhomavg(nrhomax)

      real drhodx, drhodxx, drhody, drhodyy, gradrho

      real dxpsi, dxxpsi, dxypsi, dypsi, dyypsi
        

      real rho(nxdim, nydim), theta(nxmx, nymx),
     &     rhohatx(nxmx, nymx), rhohaty(nxmx, nymx),
     &     theta0(nxmx, nymx),
     &     bx(nxmx, nymx), by(nxmx, nymx), bz(nxmx, nymx),
     &     btau(nxmx, nymx), bzeta(nxmx, nymx),
     &     dxdth(nxmx, nymx), dzdth(nxmx, nymx), xntau(nxmx, nymx),
     &     xn(nxmx, nymx), xkte(nxmx, nymx), xkti(nxmx, nymx),
     &     xkti2(nxmx, nymx), xkti3(nxmx, nymx),
     &     xn1a(nxmx, nymx), xnea(nxmx, nymx), xn2a(nxmx, nymx),
     &     xn3a(nxmx, nymx), omgce(nxmx, nymx),
     &     omgci1(nxmx, nymx), omgci2(nxmx, nymx), omgci3(nxmx, nymx),
     &     omgpe2(nxmx, nymx),
     &     omgp12(nxmx, nymx), omgp22(nxmx, nymx), omgp32(nxmx, nymx),
     &     xiota(nxmx, nymx), xlprl(nxmx, nymx),
     &     xlprl_e(nxmx, nymx), rho_tor2d(nxmx, nymx)
     
      real dxbx, dxxbx, dxybx, dybx, dyybx
      real dxby, dxxby, dxyby, dyby, dyyby
      real dxbz, dxxbz, dxybz, dybz, dyybz
     
      real dxbmod, dxxbmod, dxybmod, dybmod, dyybmod
          
          
      real wdote(nxdim, nydim), wdoti1(nxmx, nymx),
     &     wdoti2(nxmx, nymx), wdoti3(nxmx, nymx), wdot(nxmx, nymx)
      real wdoti4(nxdim, nydim), wdoti5(nxmx, nymx),
     &     wdoti6(nxmx, nymx)

      real xjprl(nxdim, nydim), xjtot, fyptot

      real redotj1(nxdim, nydim), redotj2(nxdim, nydim),
     &     redotj3(nxdim, nydim), redotj4(nxdim, nydim),
     &     redotj5(nxdim, nydim), redotj6(nxdim, nydim),
     &     redotje(nxdim, nydim), redotji(nxdim, nydim)

      real redotjs(nxmx, nymx)

      real betan, betan2, betan3, betan_slo, betate, betati, betati2,
     &    betati3, taue

      real yprimec(nxmx), ycourse(nxmx),
     &     yprime(nxmx), y(nxmx), dy, dyc

      integer icnc
      logical ismine, ismine1, ismine2, ismine3

      character(4) suffix

      integer mb,nb,myid,nproc,wantnproc,myrow,mycol
      integer nprow,npcol,icontxt, wantnprocs
      integer lrindx,lcindx,rsrc,csrc,ipos,ii,jj


      integer rsrc1, csrc1, irnc1, icnc1
      integer rsrc2, csrc2, irnc2, icnc2
      integer rsrc3, csrc3, irnc3, icnc3


*     ---------------------------------
*     set default values of input data:
*     ---------------------------------
      nnodex = nmodesx  !assignment to nnodex      
      nnoderho = nnodex/2

      n_prof_flux = 1

      isolve = 1
      ftrap = 1
      eslowev = 3.5e+06
      amu_slo = 4.0
      z_slo = 2.0
      eta_slo = 0.0

      nnode_local = 0
      nnode_overlap = 0

      iprofile = 1
      nboundary = 1
      nprow = 8
      npcol = 8

      nwdot = 0

      lmax = 5
      ibessel = 1
      xnuomg = 0.0

      xnuead = 0.0000E+00
      xnu1ad = 0.0000E+00
      xnu2ad = 0.0000E+00
      rant = 1.6

      dthetant0 = 60.
      dpsiant0 = .05

      psilim = .99
      psiant = 0.95
      psimol = 0.90
      psipne = 0.50
      psipte = .30
      psipti1 = .30
      psipti2 = .30
      psipti3 = .30

      yant = .1

      te0  = 4.2900E+03
      ti0  = 7.0700E+03
      ti02 = 7.0700E+03
      ti03 = 7.0700E+03

      inu = 0
      iprint = 50
      iexact = 1
      delta0 = 0.0000E+00
      xwall = 0.0000E+00
      iroot = 2
      iequat = 1
      igeom = 5
      epszet = 1.0000E-07
      iqx = 4
      iqprof = 1
      iez = 0
      icurve = 2
      nphi = 33
      amu1 = 2.0000E+00
      amu2 = 1.0000E+00
      z1 = 1.0000E+00
      z2 = 1.0000E+00
      eta = 4.5000E-01
      q0 = 1.0
      ekappa = 2.0

      ymax = 0.0
      xnurf = 3.2000E+07
      aplasm = 7.0000E-01
      alim = 100.0
      grad = 0.0

      xn0 = 3.1100E+19
      flat = 0.0000E+00
      b1rat = 7.0000E-01
      b2rat = 1.3000E+00
      curdnx = 0.0000E+00
      curdny = 1.0
      curdnz = 0.0000E+00
      prfin = 0.0
      nstep  = 16
      nabs = 2
      xnuabs = 0.0000E+00
      xbnch = 0.0000E+00
      xleft = -7.000E-01
      xright = 7.0000E-01
      isigma = 1
      nzfun = 1
      qavg0 = 1.0
      iabsorb = 2
      itemp = 0

      xnlim  = 0.0000E+00
      telim  = 0.0000E+00
      tilim  = 0.0000E+00
      ti2lim = 0.0000E+00
      ti3lim = 0.0000E+00


      nfreqm = 1
      dfreq = 0.0000E+00
      nkzm = 1
      dkz = 0.0000E+00
      idens = 0
      r0 = 2.3000E+00
      xnudip = 2.5000E+00
      adip = 0.0000E+00
      efold = 0.0000E+00
      amu3 = 1.2000E+01
      z3 = 6.0000E+00
      eta3 = 4.6600E-02
      xnu3ad = 0.0000E+00


      xdelta = 5.5000E-01
      wdelta = 0.0000E+00
      xdelt2 = -7.000E-02
      wdelt2 = 0.0000E+00
      zeffcd = 2.5000E+00
      rzoom1 = 0.0
      rzoom2 = 0.0
      yzoom1 = 0.0
      yzoom2 = 0.0

      ibackground = 1


      alphan = 1.0
      alphate = 1.0
      alphati = 1.0
      idiag = 5
      jdiag = 4






      idiag = nmodesx / 2
      jdiag = nmodesy / 2




      nrhs = 1
      nmaxe=1


      nkx2 = nmodesx / 2
      nkx1 = - nmodesx / 2 + 1


      nky2 = nmodesy / 2
      nky1 = - nmodesy / 2 + 1
      nnodey = nmodesy

      jequat = nnodey / 2
      icenter = nnodex / 2

      if (qavg0 .ne. 0.0) xiota0 = 1./qavg0

      rholim = sqrt(psilim)

      q = 1.6e-19
      if(te0 .eq. 0.0)te0 = ti0

      t0e = te0
      t0i = ti0
      t0i2 = ti02
      t0i3 = ti03


      teedge = 400.0


  101 format (10i10)


      qhat = qavg0

      eslow = eslowev * q
      t0e = t0e   * q
      t0i = t0i   * q
      t0i2 = t0i2 * q
      t0i3 = t0i3 * q

      telim   = telim  * q
      tilim   = tilim  * q
      ti2lim  = ti2lim * q
      ti3lim  = ti3lim * q


      teedge = teedge * q

      xme = 9.11e-31
      xmh = 1.67e-27
      xmi1 = amu1 * xmh
      xmi2 = amu2 * xmh
      xmi3 = amu3 * xmh
      xmi_slo = amu_slo * xmh

      qi1 = z1 * q
      qi2 = z2 * q
      qi3 = z3 * q
      qi_slo = z_slo * q



      qe = -q
      zi = cmplx(0.0,1.0)
      eps0 = 8.85e-12
      pi = 3.141592654
      twopi = 2. * pi
      xlnlam = 20.0
      xmu0 = 1.26e-06
      clight = 1.0 / sqrt(eps0 * xmu0)
      xmax = rwright - rwleft
      ymax = ytop - ybottom

      xkphi0 = nphi / rt



      omgrf = 2.0 * pi * xnurf
      omgci10 = qi1 * b0 / xmi1
      vthi10 = sqrt(2.0 * t0i / xmi1)
      v0i = vthi10
      rhoi10 = vthi10 / omgci10
      rhoi1overl = rhoi10 / xmax

      vphase = omgrf / xkphi0
      vthe = sqrt(t0e / xme)
      vsound = sqrt(teedge / xmi1)
      wphase = 0.0
      if(vthe .ne. 0.0)wphase = vphase / vthe


      xkthrho = 0.2
      omgestar = xkthrho * vsound / ymax
      xkthdx = 1.0
      domgk = omgestar
      shearedge = domgk / xkthdx


      xk0 = omgrf / clight
      rnz = xkphi0 / xk0
      xn1 = xn0 / z1 * (1.0 - z2 * eta - z3 * eta3 - z_slo * eta_slo)
      eta1 = xn1 / xn0
      xn2 = xn0 * eta
      eta2 = xn2 / xn0
      xn3 = xn0 * eta3

      xn_slo = xn0 * eta_slo


      xwleft = rwleft - rt
      xwright = rwright - rt

*--------------------------------------------
*     Define x mesh: x(i), xprime(i), capr(i)
*     xprime: 0 to xmax
*     x(i) : -xmax / 2.0   to   xmax / 2.0
*--------------------------------------------
      dx = xmax / nnodex
      do i = 1, nnodex
         xprime(i) = (i-1) * dx + dx / 2.0
c--   Note: the code gives slightly smoother results with dx/2.0 added

         x(i) = xprime(i) + xwleft
         capr(i) = rt + x(i)

         xkphi(i) = nphi / capr(i)

      end do


      if(rzoom1 .eq. 0.0)rzoom1 = capr(1)
      if(rzoom2 .eq. 0.0)rzoom2 = capr(nnodex)

      izoom1 = int((rzoom1 - rwleft - dx / 2.0) / dx) + 1
      izoom2 = int((rzoom2 - rwleft - dx / 2.0) / dx) + 1

c      if (myid.eq.0) then
c         write(6, 1312) izoom1
c         write(6, 1312) izoom2
c      end if

*--------------------------------------------
*     Define y mesh: y(j), yprime(j)
*     yprime: 0 to ymax
*     y(j) : -ymax / 2.0   to   ymax / 2.0
*--------------------------------------------
      dy = ymax / nnodey
      do j = 1, nnodey
         yprime(j) = (j-1) * dy + dy / 2.0

c--      Note: the code gives slightly smoother results with dy/2.0 added

         y(j) = yprime(j) + ybottom
      end do



      if(yzoom1 .eq. 0.0)yzoom1 = ybottom
      if(yzoom2 .eq. 0.0)yzoom2 = ytop

      jzoom1 = int((yzoom1 - ybottom - dy / 2.0) / dy) + 1
      jzoom2 = int((yzoom2 - ybottom - dy / 2.0) / dy) + 1


*-----------------------
*     define theta0(i,j)
*-----------------------
      do i = 1, nnodex
         do j = 1, nnodey
            if(y(j) .ne. 0.0 .or. x(i) .ne. 0.0) then
               theta0(i,j) = atan2(y(j), x(i))
               if(theta0(i,j) .ge. 0.0) theta(i,j) = theta0(i,j)
               if(theta0(i,j) .lt. 0.0) theta(i,j) =
     &                                       theta0(i,j) + twopi
            end if
         end do
      end do



*-----------------------------
*     Define rho mesh: rhon(n)
*     rhon: 0 to rhomax
*-----------------------------

c      rhomax = 1.0
      drho = rhomax / (nnoderho - 1)
      do n = 1, nnoderho
         rhon(n) = (n-1) * drho
      end do


      allocate(x_fine(nnodex_fine) )
      allocate(y_fine(nnodey_fine) )
      allocate(capr_fine(nnodex_fine) )
      allocate(rho_fine(nnodex_fine, nnodey_fine) )

*     -----------------------------
*     Define fine x mesh: x_fine(i)
*     x_fine(i) : 0 to xmax
*     -----------------------------
      dx_fine = xmax / nnodex_fine

      do i = 1, nnodex_fine
         x_fine(i) = (i - 1) * dx_fine
     &      + dx_fine / 2.0
         capr_fine(i) = rt + x_fine(i) + xwleft
      end do



*     -----------------------------
*     Define fine y mesh: y_fine(j)
*     y_fine(j) : 0 to ymax
*     -----------------------------
      dy_fine = ymax / nnodey_fine

      do j = 1, nnodey_fine
         y_fine(j) = (j - 1) * dy_fine
     &      + dy_fine / 2.0
      end do

*     ----------------------------------
*     Interpolate rho onto fine 2-D grid
*     ----------------------------------

      do i = 1, nnodex_fine
         do j = 1, nnodey_fine

            call intplt2(x_fine(i), y_fine(j), rho_fine(i,j),
     &         nnodex, nnodey, rho, nxmx, nymx, dx, dy)

         end do
      end do



*     ----------------------------------------------
*     Calculate the differential volume on half mesh:
*     ---------------------------------------------- 
     
      dvol  = 0.0
      darea = 0.0           
     
      do i = 1, nnodex_fine
         do j = 1, nnodey_fine
            n = int(rho_fine(i,j) / drho) + 1
            if(n .le. nnoderho .and. n .ge. 1)then
               dvol(n) = dvol(n)+ dx_fine * dy_fine *twopi *capr_fine(i)
               darea(n) = darea(n)+ dx_fine * dy_fine !JCW * capr_fine(i)/ rt
            end if
         end do
      end do
      


*     --------------------------------------------
*     Calculate the integrated volume on even mesh:
*     --------------------------------------------     
      
      volume = 0.0
            
      volume(1) = 0.0      
      do n = 1, nnoderho - 1
         volume(n+1) = volume(n) + dvol(n)
      end do
           


*     -------------------------
*     Do flux surface averages:
*     -------------------------

      call fluxavg2(redotje, redotjeavg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(redotj1, redotj1avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(redotj2, redotj2avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(redotj3, redotj3avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(redotj4, redotj4avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(redotj5, redotj5avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(redotj6, redotj6avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(redotji, redotjiavg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(xjprl, xjprlavg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(wdote, wdoteavg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(wdoti1, wdoti1avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(wdoti2, wdoti2avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(wdoti3, wdoti3avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(wdoti4, wdoti4avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(wdoti5, wdoti5avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      call fluxavg2(wdoti6, wdoti6avg, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, rt,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)
     
     

     
*     --------------------------------------------
*     Integrate flux averaged quantities over rhon:
*     --------------------------------------------
      call rhograte(rhon, xjprlavg, 1, nnoderho, xjprl_int, 
     &                                               nrhomax, darea)
     
      call rhograte(rhon, redotjeavg, 1, nnoderho, redotje_int, 
     &                                                nrhomax, dvol)
      call rhograte(rhon, redotj1avg, 1, nnoderho, redotji1_int, 
     &                                                nrhomax, dvol)
      call rhograte(rhon, redotj2avg, 1, nnoderho, redotji2_int, 
     &                                                nrhomax, dvol)
      call rhograte(rhon, redotj3avg, 1, nnoderho, redotji3_int, 
     &                                                nrhomax, dvol)
      call rhograte(rhon, redotj4avg, 1, nnoderho, redotji4_int, 
     &                                                nrhomax, dvol)
      call rhograte(rhon, redotj5avg, 1, nnoderho, redotji5_int, 
     &                                                nrhomax, dvol)
      call rhograte(rhon, redotj6avg, 1, nnoderho, redotji6_int, 
     &                                                nrhomax, dvol)
     
     
      call rhograte(rhon, wdoteavg,  1, nnoderho, wdote_int, 
     &                                                nrhomax, dvol)
      call rhograte(rhon, wdoti1avg, 1, nnoderho, wdoti1_int, 
     &                                                nrhomax, dvol) 
      call rhograte(rhon, wdoti2avg, 1, nnoderho, wdoti2_int, 
     &                                                nrhomax, dvol) 
      call rhograte(rhon, wdoti3avg, 1, nnoderho, wdoti3_int, 
     &                                                nrhomax, dvol) 
      call rhograte(rhon, wdoti4avg, 1, nnoderho, wdoti4_int, 
     &                                                nrhomax, dvol) 
      call rhograte(rhon, wdoti5avg, 1, nnoderho, wdoti5_int, 
     &                                                nrhomax, dvol) 
      call rhograte(rhon, wdoti6avg, 1, nnoderho, wdoti6_int, 
     &                                                nrhomax, dvol)
     
     
      wdoti1_dvol = 0.0 
      wdoti2_dvol = 0.0
      wdoti3_dvol = 0.0
      wdoti4_dvol = 0.0
      wdoti5_dvol = 0.0
      wdoti6_dvol = 0.0
      wdote_dvol  = 0.0
     
      redotj1_dvol = 0.0
      redotj2_dvol = 0.0
      redotj3_dvol = 0.0
      redotj4_dvol = 0.0
      redotj5_dvol = 0.0
      redotj6_dvol = 0.0
      redotje_dvol = 0.0
     
      do n = 1, nnoderho - 1

         wdoti1_dvol(n+1) = wdoti1avg(n) * dvol(n) 
         wdoti2_dvol(n+1) = wdoti2avg(n) * dvol(n)
         wdoti3_dvol(n+1) = wdoti3avg(n) * dvol(n)
         wdoti4_dvol(n+1) = wdoti4avg(n) * dvol(n)
         wdoti5_dvol(n+1) = wdoti5avg(n) * dvol(n)
         wdoti6_dvol(n+1) = wdoti6avg(n) * dvol(n)
         wdote_dvol (n+1) = wdoteavg(n)  * dvol(n)
     
         redotj1_dvol(n+1) = redotj1avg(n) * dvol(n)
         redotj2_dvol(n+1) = redotj2avg(n) * dvol(n)
         redotj3_dvol(n+1) = redotj3avg(n) * dvol(n)
         redotj4_dvol(n+1) = redotj4avg(n) * dvol(n)
         redotj5_dvol(n+1) = redotj5avg(n) * dvol(n)
         redotj6_dvol(n+1) = redotj6avg(n) * dvol(n)
         redotje_dvol(n+1) = redotjeavg(n) * dvol(n)
      end do      
                    
*     --------------------------------
*     Open and write SWIM output file:
*     --------------------------------
      if (myid.eq.0) then        
                  
         open(unit=99, file='out_swim', status='unknown',
     &                                              form='formatted')
         
            write(99, 309) nnoderho
            write(99, 310) (rhon(n), n = 1, nnoderho)
         
            write(99, 310) (dvol(n), n = 1, nnoderho - 1)
         
            write(99, 310) (redotje_int(n), n = 1, nnoderho)
            write(99, 310) (redotji1_int(n), n = 1, nnoderho)
            write(99, 310) (redotji2_int(n), n = 1, nnoderho)
            write(99, 310) (redotji3_int(n), n = 1, nnoderho)
            write(99, 310) (redotji4_int(n), n = 1, nnoderho)
            write(99, 310) (redotji5_int(n), n = 1, nnoderho)
            write(99, 310) (redotji6_int(n), n = 1, nnoderho)


            write(99, 310) (wdote_int(n), n = 1, nnoderho)
            write(99, 310) (wdoti1_int(n), n = 1, nnoderho)
            write(99, 310) (wdoti2_int(n), n = 1, nnoderho)
            write(99, 310) (wdoti3_int(n), n = 1, nnoderho)
            write(99, 310) (wdoti4_int(n), n = 1, nnoderho)
            write(99, 310) (wdoti5_int(n), n = 1, nnoderho)
            write(99, 310) (wdoti6_int(n), n = 1, nnoderho)
            
            write(99, 310) (redotje_dvol(n), n = 1, nnoderho)
            write(99, 310) (redotj1_dvol(n), n = 1, nnoderho)
            write(99, 310) (redotj2_dvol(n), n = 1, nnoderho)
            write(99, 310) (redotj3_dvol(n), n = 1, nnoderho)
            write(99, 310) (redotj4_dvol(n), n = 1, nnoderho)
            write(99, 310) (redotj5_dvol(n), n = 1, nnoderho)
            write(99, 310) (redotj6_dvol(n), n = 1, nnoderho)


            write(99, 310) (wdote_dvol(n), n = 1, nnoderho)
            write(99, 310) (wdoti1_dvol(n), n = 1, nnoderho)
            write(99, 310) (wdoti2_dvol(n), n = 1, nnoderho)
            write(99, 310) (wdoti3_dvol(n), n = 1, nnoderho)
            write(99, 310) (wdoti4_dvol(n), n = 1, nnoderho)
            write(99, 310) (wdoti5_dvol(n), n = 1, nnoderho)
            write(99, 310) (wdoti6_dvol(n), n = 1, nnoderho)
                         
         close (99)      
                 
                  
         write(38, 309) nnodex, nnodey
         write(38, 310) psilim
         write(38, 310) (x(i), i = 1, nnodex)
         write(38, 310) (y(j), j = 1, nnodey)
         write(38, 310) (capr(i), i = 1, nnodex)

         write(38, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((theta(i, j), i = 1, nnodex), j = 1, nnodey)
        

         write(38, 309) nnoderho
         write(38, 310) (rhon(n), n = 1, nnoderho)
         
         write(38, 310) (dvol(n), n = 1, nnoderho)
         write(38, 310) (volume(n), n = 1, nnoderho)
         
         write(38, 310) (redotjeavg(n), n = 1, nnoderho)
         write(38, 310) (redotj1avg(n), n = 1, nnoderho)
         write(38, 310) (redotj2avg(n), n = 1, nnoderho)
         write(38, 310) (redotj3avg(n), n = 1, nnoderho)
         write(38, 310) (redotj4avg(n), n = 1, nnoderho)
         write(38, 310) (redotj5avg(n), n = 1, nnoderho)
         write(38, 310) (redotj6avg(n), n = 1, nnoderho)
         write(38, 310) (redotjiavg(n), n = 1, nnoderho)

         write(38, 310) (xjprlavg(n), n = 1, nnoderho)
         write(38, 310) (xjprl_int(n), n = 1, nnoderho)
         
         write(38, 310) (wdoteavg(n), n = 1, nnoderho)
         write(38, 310) (wdoti1avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti2avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti3avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti4avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti5avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti6avg(n), n = 1, nnoderho)
         
                

         write(38, 310) (redotje_int(n),  n = 1, nnoderho)
         write(38, 310) (redotji1_int(n), n = 1, nnoderho)
         write(38, 310) (redotji2_int(n), n = 1, nnoderho)
         write(38, 310) (redotji3_int(n), n = 1, nnoderho)
         write(38, 310) (redotji4_int(n), n = 1, nnoderho)
         write(38, 310) (redotji5_int(n), n = 1, nnoderho)
         write(38, 310) (redotji6_int(n), n = 1, nnoderho)

         write(38, 310) (wdote_int(n),  n = 1, nnoderho)
         write(38, 310) (wdoti1_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti2_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti3_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti4_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti5_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti6_int(n), n = 1, nnoderho)


         write(38, 310) (redotje_dvol(n), n = 1, nnoderho)
         write(38, 310) (redotj1_dvol(n), n = 1, nnoderho)
         write(38, 310) (redotj2_dvol(n), n = 1, nnoderho)
         write(38, 310) (redotj3_dvol(n), n = 1, nnoderho)
         write(38, 310) (redotj4_dvol(n), n = 1, nnoderho)
         write(38, 310) (redotj5_dvol(n), n = 1, nnoderho)
         write(38, 310) (redotj6_dvol(n), n = 1, nnoderho)


         write(38, 310) (wdote_dvol(n),  n = 1, nnoderho)
         write(38, 310) (wdoti1_dvol(n), n = 1, nnoderho)
         write(38, 310) (wdoti2_dvol(n), n = 1, nnoderho)
         write(38, 310) (wdoti3_dvol(n), n = 1, nnoderho)
         write(38, 310) (wdoti4_dvol(n), n = 1, nnoderho)
         write(38, 310) (wdoti5_dvol(n), n = 1, nnoderho)
         write(38, 310) (wdoti6_dvol(n), n = 1, nnoderho)
         
      end if
     
      do n = 1, nnoderho
         xkti3avg(n) = xkti3avg(n) / q
      end do


  310 format(1p,6e12.4)
 8310 format(1p,6e14.6)
  309 format(10i10)
 9310 format(1p,7e12.4)
  311 format(1p,10e12.4)
 3117 format(1p,11e12.4)

 1314 format (2i10, 1p,11e12.4)
13149 format (4i10, 1p,8e12.4)
 1313 format(10i10)
 1311 format(1p,9e12.4)
   10 format(i10,1p,4e10.3,i10,1p,e10.3)
 1010 format(1f4.0,4f8.3,3f7.3,1f8.3,1f9.3,2f8.3)
 1009 format(3x,"        frequency  = ",1p,e12.4," hertz ")
 1012 format(3x,"             omgrf = ",1p,e12.4," hertz ")
 2012 format(3x,"2*pi*rt/Ly * real part of impedance (resistance) = ",
     1   1p,e12.4," ohms     ")
 2013 format(3x,
     1   "2*pi*rt/Ly * imaginary part of impedance (reactance) = ",
     1   1p,e12.4," ohms     ")
 1014 format(3x,"               xkz = ",1p,e12.4," m-1   ")
 1321 format(3x,"         vph / vth = ",1p,e12.4,"       ")
 1391 format(3x," critical shear(0) = ",1p,e12.4," s-1   ")
 1392 format(3x," critical shear(a) = ",1p,e12.4," s-1   ")
 1393 format(3x,"         mu neo(a) = ",1p,e12.4," s-1   ")
 1322 format(3x,"               vph = ",1p,e12.4,"       ")
 1323 format(3x,"               vth = ",1p,e12.4,"       ")
 1714 format(3x,"               xk0 = ",1p,e12.4," m-1   ")
 1021 format(3x,"        n parallel = ",1p,e12.4,"       ")
 1812 format(3x,"                rt = ",1p,e12.4," m     ")
 1822 format(3x,"            aplasm = ",1p,e12.4," m     ")
 1823 format(3x,"              rant = ",1p,e12.4," m     ")
 1809 format(3x,"                b0 = ",1p,e12.4," T     ")
 1813 format(3x,"              xn10 = ",1p,e12.4," m-3   ")
 6813 format(3x,"              xne0 = ",1p,e12.4," m-3   ")
 1814 format(3x,"              xn20 = ",1p,e12.4," m-3   ")
 1834 format(3x,"              xn30 = ",1p,e12.4," m-3   ")
 6834 format(3x,"              eta1 = ",1p,e12.4,"       ")
 6835 format(3x,"              eta2 = ",1p,e12.4,"       ")
 6836 format(3x,"              eta3 = ",1p,e12.4,"       ")
 1815 format(3x,"               te0 = ",1p,e12.4," eV    ")
 1821 format(3x,"               ti0 = ",1p,e12.4," eV    ")
 1016 format(3x," xnue/omgrf ad hoc = ",1p,e12.4,"       ")
 1017 format(3x," xnu1/omgrf ad hoc = ",1p,e12.4,"       ")
 1018 format(3x," xnu2/omgrf ad hoc = ",1p,e12.4,"       ")
 1013 format(3x,"              nphi = ",i12,"       ")
 7013 format(3x,"           nmodesx = ",i12,"       "/
     1       3x,"           nmodesy = ",i12,"       ")

 7113 format(3x,"             nwdot = ",i12,"       ")

 7014 format(3x,"              lmax = ",i12,"       ")
 7015 format(3x,"           ibessel = ",i12,"       ")
 7115 format(3x,"             nzfun = ",i12,"       ")
 7016 format(3x,"         rhoi1 / L = ",1p,e12.4,"       ")
 7017 format(3x,"             rhoi1 = ",1p,e12.4," m     ")
 7217 format(3x,"             qavg0 = ",1p,e12.4,"       ")
 1020 format(3x,"            nnodex = ",i12, "       "/
     1       3x,"            nnodey = ",i12, "       ")

 3013 format(3x,"                i0 = ",i12,"       ")
30131 format(3x,"             ileft = ",i12,"       ")
30132 format(3x,"            iright = ",i12,"       ")
 3014 format(3x," xnuii(0)/omgti(0) = ",1p,e12.4," s-1   ")
 3015 format(3x,"           vthi(0) = ",1p,e12.4," m/s   ")
 3016 format(3x,"          omgti(0) = ",1p,e12.4," s-1   ")
 3017 format(3x,"           xnup(0) = ",1p,e12.4," s-1   ")
 3018 format(3x,"          eps**1.5 = ",1p,e12.4,"       ")
 7116 format(3x,"         xnu / omg = ",1p,e12.4,"       ")

 1309 format(3x,
     1   "             total power absorbed = ",1p,e12.4," watts/m" )
11091 format(3x,
     1   "       total power absorbed (old) = ",1p,e12.4," watts/m" )

 1109 format(
     1   3x, "      power absorbed by electrons = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       "/
     1   3x, "  power absorbed by majority ions = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       "/
     1   3x, "  power absorbed by minority ions = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       "/
     1   3x, "power absorbed by 3rd ion species = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by slowing species = ",1e12.5,
     &   " watts/m ", " = ", f10.4, " %       "/
     1   3x, "             total power absorbed = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       ")

71109 format(
     1   3x, "      power absorbed by electrons = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       "/
     1   3x, "  power absorbed by majority ions = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       "/
     1   3x, "  power absorbed by minority ions = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       "/
     1   3x, "power absorbed by 3rd ion species = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       "/
     1   3x, "             total power absorbed = ",1e12.5,
     1   " watts/m ", " = ", f10.4, " %       ")

81109 format(
     &   3x, "  power absorbed by majority ions = ",1e12.5,
     &   " watts/m "/
     1   3x, "  power absorbed by minority ions = ",1e12.5,
     1   " watts/m "/
     1   3x, "      power absorbed by electrons = ",1e12.5,
     1   " watts/m ")

 1112 format(
     1   3x,"      power absorbed by electrons = ",1f12.4," %       "/
     1   3x,"  power absorbed by majority ions = ",1f12.4," %       "/
     1   3x,"  power absorbed by minority ions = ",1f12.4," %       "/
     1   3x,"power absorbed by 3rd ion species = ",1f12.4," %       "/
     1   3x,"             total power absorbed = ",1f12.4," %       ")
 1189 format(
     1   3x,"                             ref1 = ",f12.3," %"/
     1   3x,"                            trans = ",f12.3," %"/
     1   3x,"                           trans2 = ",f12.3," %"/
     1   3x,"                             conv = ",f12.3," %"/
     1   3x,"                            conv2 = ",f12.3," %"/
     1   3x,"  power absorbed by minority ions = ",f12.3," %"/
     1   3x,"  power absorbed by majority ions = ",f12.3," %"/
     1   3x,"power absorbed by 3rd ion species = ",f12.3," %"/
     1   3x,"      power absorbed by electrons = ",f12.3," %"/
     1   3x,"             total power absorbed = ",f12.3," %"/
     1   3x,"                        total sum = ",f12.3," %")
 1289 format(
     1   3x, '                         x(iedge) =',1p,e12.4,' m'/
     1   3x, '                               bz =',1p,e12.4,' T'/
     1   3x, '               real(eps parallel) =',1p,e12.4,'  '/
     1   3x, '     estimate: real(eps parallel) =',1p,e12.4,'  ' /
     1   3x, '                real (eps left)   =',1p,e12.4,'  ' /
     1   3x, '      estimate: real (eps left)   =',1p,e12.4,'  ' /
     1   3x, '                real (eps right)  =',1p,e12.4,'  ' /
     1   3x, '     estimate:  real (eps right)  =',1p,e12.4,'  '/
     1   3x, '                               xn =',1p,e12.4,' m-3'/
     1   3x, '                               LN =',1p,e12.4,' m'/
     1   3x, '                        mod2 (ez) =',1p,e12.4,' (V/m)**2'/
     1   3x, '                            L rf  =',1p,e12.4,' m'/
     1   3x, '                x=omgrf/omgci(i)  =',1p,e12.4,'    '/
     1   3x, '                             xkti =',1p,e12.4,' deg'/
     1   3x, '                             xkte =',1p,e12.4,' deg')

 1290 format(
     1   3x,"             alpha rf parallel    = ",1p,e12.4,"  "/
     1   3x," estimate of alpha rf parallel    = ",1p,e12.4,"  "/
     1   3x,"                 alpha rf left    = ",1p,e12.4,"  "/
     1   3x,"     estimate of alpha rf left    = ",1p,e12.4,"  "/
     1   3x,"                 alpha rf right   = ",1p,e12.4,"  "/
     1   3x,"     estimate of alpha rf right   = ",1p,e12.4,"  "/
     1   3x,"                 alpha rf total   = ",1p,e12.4,"  "/
     1   3x," electron ponderomotive potential = ",1p,e12.4," eV"/
     1   3x,"      ion ponderomotive potential = ",1p,e12.4," eV")
      
 1110 format(
     1   3x,'            real power emitted = ',1p,e12.4,' watts/m',/
     1   3x,'       imaginary power emitted = ',1p,e12.4,' watts/m' )
 1113 format(3x,'driven current per meter = ',1p,e12.4,' Amps/m')
 1114 format(3x,'current driven per watt = ',1p,e12.4,' A/W')
 1215 format(3x,'current drive efficiency = ',1p,e12.4,' A/W/m**2')
 1216 format(3x,'total RF power = ',1p,e12.4,' Watts')
 1217 format(3x,'total driven current = ',1p,e12.4,' Amps')

 1220 format(3x,'total x force on antenna = ',1p,e12.4,' Nt/m'/
     &       3x,'total y force on antenna = ',1p,e12.4,' Nt/m'/
     &       3x,'total z force on antenna = ',1p,e12.4,' Nt/m')


 1120 format(3x,'total x force on antenna = ',1p,e12.4,' Nt/m'/
     &       3x,'total x force on  plasma = ',1p,e12.4,' Nt/m'/)

 1121 format(3x,'total y force on antenna = ',1p,e12.4,' Nt/m'/
     &       3x,'total y force on  plasma = ',1p,e12.4,' Nt/m'/)

 1122 format(3x,'total z force on antenna = ',1p,e12.4,' Nt/m'/
     &       3x,'total z force on  plasma = ',1p,e12.4,' Nt/m'/)

 1123 format(3x,'total poloidal force on plasma = ',1p,e12.4,' Nt/m'/)


  163 format("0")


 1218 format(3x,'pscale = ',1p,e12.4,' ')
 1219 format(3x,'gamma = ',1p,e12.4,' ')
 1111 format(3x," At x = ",1p,e12.4," m",
     1       3x," ifail = ",i5)
  162 format("1")
  169 format(" ")

 2162 format(1p,8e12.4)
 2163 format(i5, 1p,11e12.3)
 2165 format(3i10,5e10.3)
 1002 format(2i10,7e10.3)
11313 format(2i10, 1p,8e12.4)
 1000 format(1i10,7e10.3)

 1001 format(8e10.3)


  164 format(" ",1x," i  ",  3x," R(x) ",
     1           6x,"  x   ",6x," btau ",6x,"bzeta ",6x,"re om1",
     1           6x,"re om2",6x,"re om3",6x," bmod ",6x," xne  ",
     &           6x,"  Te  ")

 9164 format(" ",1x," j  ",  3x," R(x) ",
     1           6x,"  y   ",6x," btau ",6x,"bzeta ",6x,"re om1",
     1           6x,"re om2",6x,"re om3",6x," bmod ",6x," xne  ",
     &           6x,"  Te  ")



 1164 format(" ",3x," R(x) ",

     1           6x,"  x   ", 6x,"vymean",6x,"mu*vyx",
     1           5x,"mu*vzxx",6x,"mu*vzx",
     1           6x,"      ",6x,"      ",6x,"      ",
     1           6x,"      ")

 4000 format(" ",3x,"  x    ",5x,"re k**2",5x,"re k**2",
     1           5x,"re k**2",5x,"re k**2",5x,"re k**2",
     1           5x,"re k**2",5x,"re k**2",5x,"re k**2")

 7000 format(" ",3x,"  x   ",6x," amp1 ",6x," amp2 ",
     1           6x," amp3 ",6x," amp4 ",6x," amp5 ",
     1           6x," amp6 ",6x," amp7 ",6x," amp8 ")

 7100 format(" ",3x,"  x   ",6x,"  sx1 ",6x,"  sx2 ",
     1           6x,"  sx3 ",6x,"  sx4 ",6x,"  sx5 ",
     1           6x,"  sx6 ",6x,"  sx7 ",6x,"  sx8 ",
     1           6x,"sx sum")
 4002 format(" ",3x,"  x   ",6x,"re dt1",6x,"re dt2",
     1           6x,"re dt3",6x,"re dt4",6x,"re dt5",
     1           6x,"re dt6",6x,"re dt7",6x,"re dt8")

 4001 format(" ",3x,"  x   ", 5x,"im k**2",5x,"im k**2",
     1           5x,"im k**2",5x,"im k**2",5x,"im k**2",
     1           5x,"im k**2",5x,"im k**2",5x,"im k**2")

 4003 format(" ",3x,"  x   ",6x,"im dt1",6x,"im dt2",
     1           6x,"im dt3",6x,"im dt4",6x,"im dt5",
     1           6x,"im dt6",6x,"im dt7",6x,"im dt8")


 1650 format(3x,"dispersion relation             ")
 7650 format(3x,"amplitude of modes              ")
 7750 format(3x,"flux of power carried by modes  ")
 1651 format(3x,"dispersion relation check       ")
 1115 format(3x,"        kh = ",1p,e12.4," m-1   ")
 1116 format(3x,"        lp = ",1p,e12.4," m     ")
 1117 format(3x,"      eps2 = ",1p,e12.4," tesla ")
 1118 format(3x,"avg iotabr = ",1p,e12.4,"       ")
 1119 format(2x,"ellipticity = ",1p,e12.4,"       ")

 5555 continue

  899 format('total cpu time used =',f9.3," min")

 9930 format(3x, 'electron power and flow')
 1213 format(9x, 'i', 6x, 'R', 9x, 'wdote', 7x, 'fye')
  933 format(3x, 'species #3 ion power and flow')
 1315 format(9x, 'i', 6x, 'R', 9x, 'wdoti3', 7x, 'fyi3')
  932 format(3x, 'minority ion power and flow')
 1214 format(9x, 'i', 6x, 'R', 9x, 'wdoti2', 7x, 'fyi2')
  931 format(3x, 'majority ion power and flow')
 1319 format(9x, 'i', 6x, 'R', 9x, 'wdoti1', 7x, 'fyi1')

 9321 format(3x, 'Electric field:')
 1394 format(9x, 'i', 6x, 'R', 9x, 'Re Ex', 7x, 'Im Ex',
     1                         7x, 'Re Ey', 7x, 'Im Ey',
     1                         7x, 'Re Ez', 7x, 'Im Ez')


 5002 continue


 1312 format(1i10, 1p,8e12.4)
 2313 format(1p,8e12.4)
 2314 format(1p,e11.4, 1p,e11.4)

      deallocate(x_fine)
      deallocate(y_fine)
      deallocate(capr_fine)
      deallocate(rho_fine)

      return
      end

c
c***************************************************************************
c


      subroutine intplt(xgiv, ygiv, fout, nx, ny, f, 
     &   nxmax, nymax, dx, dy)

      implicit none

      integer nx, ny, nxmax, nymax, n, m
      real x, y, xgiv, ygiv, dx, dy, zeta, eta
      complex f(nxmax, nymax), fout, a, b, c, d


      x = abs(xgiv)
      y = abs(ygiv)
      fout = 0.0

      n = int(x / dx) + 1
      m = int(y / dy) + 1
               
      if (n .ge. nx) return
      if (m .ge. ny) return

      zeta = (x - (n - 1) * dx) / dx
      eta  = (y - (m - 1) * dy) / dy

      a = f(n, m)
      b = f(n+1 ,m) - f(n, m)
      c = f(n, m+1) - f(n, m)
      d = f(n+1, m+1) + f(n, m) - f(n+1, m) - f(n, m+1)

      fout = a + b * zeta + c * eta + d * zeta * eta

      return

  900 format(1p,6e12.4)
  901 format(10i10)      
  910 format(1p,8e12.4)
  920 format(2i10)
      end

c
c***************************************************************************
c
      subroutine intplt_re(xgiv, ygiv, fout, nx, ny, f, 
     &   nxmax, nymax, dx, dy)

      implicit none

      integer nx, ny, nxmax, nymax, n, m
      real x, y, xgiv, ygiv, dx, dy, zeta, eta
      real f(nxmax, nymax), fout, a, b, c, d


      x = abs(xgiv)
      y = abs(ygiv)
      fout = 0.0

      n = int(x / dx) + 1
      m = int(y / dy) + 1
               
      if (n .ge. nx) return
      if (m .ge. ny) return

      zeta = (x - (n - 1) * dx) / dx
      eta  = (y - (m - 1) * dy) / dy

      a = f(n, m)
      b = f(n+1 ,m) - f(n, m)
      c = f(n, m+1) - f(n, m)
      d = f(n+1, m+1) + f(n, m) - f(n+1, m) - f(n, m+1)

      fout = a + b * zeta + c * eta + d * zeta * eta

      return

  900 format(1p,6e12.4)
  901 format(10i10)      
  910 format(1p,8e12.4)
  920 format(2i10)
      end
      
c
c***************************************************************************
c
      subroutine intplt_to_cql3d(xgiv, ygiv, fout, nx, ny, nz, f, 
     &   nxmax, nymax, nzmax, xa, ya, k)

      implicit none

      integer nx, ny, nz, nxmax, nymax, n, m, nzmax, k
      real x, y, xgiv, ygiv, zeta, eta, dx, dy, xmax, ymax, xmin, ymin
      real f(nxmax, nymax, nzmax), fout, a, b, c, d
      real xa(nxmax), ya(nymax)
      
      x = xgiv
      y = ygiv
      fout = 0.0      

      xmin = xa(1)
      xmax = xa(nx)      
      ymin = ya(1)
      ymax = ya(ny)
      
      if(xgiv .lt. xmax .and. xgiv .gt. xmin .and. 
     &   ygiv .lt. ymax .and. ygiv .gt. ymin) then
      
         dx = (xmax - xmin) / (nx - 1)
         dy = (ymax - ymin) / (ny - 1)

         n = int((x - xmin) / dx) + 1
         m = int((y - ymin) / dy) + 1
               
         zeta = (x - (xmin + (n - 1) * dx)) / dx
         eta  = (y - (ymin + (m - 1) * dy)) / dy

         a = f(n, m, k)
         b = f(n+1 ,m, k) - f(n, m, k)
         c = f(n, m+1, k) - f(n, m, k)
         d = f(n+1, m+1, k) + f(n, m, k) - f(n+1, m, k) - f(n, m+1, k)

         fout = a + b * zeta + c * eta + d * zeta * eta
         
      end if

      return

  900 format(1p,6e12.4)
  901 format(10i10)      
  910 format(1p,8e12.4)
  920 format(2i10)
      end


c
c***************************************************************************
c

      subroutine intplt2(xgiv, ygiv, fout, nx, ny, f, nxmax, nymax, dx,
     &   dy)

      implicit none

      integer nx, ny, nxmax, nymax, i, j
      real x, y, xgiv, ygiv, dx, dy, p, q, a, b, c, d, fout
      real a0, a1, a2, a3, a4, a5
      real f(nxmax, nymax)


      x = abs(xgiv)
      y = abs(ygiv)
      fout = 100000.0

      i = int(x / dx) + 1
      j = int(y / dy) + 1

      if (i .ge. nx) return
      if (j .ge. ny) return
      if (i .le. 1) return
      if (j .le. 1) return

      p = (x - (i - 1) * dx) / dx
      q = (y - (j - 1) * dy) / dy

      a0 = .5 * q * (q - 1.)
      a1 = .5 * p * (p - 1.)
      a2 = 1.0 - p**2 - q**2 + p * q
      a3 = .5 * p * (p + 1. - 2. * q)
      a4 = .5 * q * (q + 1. - 2. * p)
      a5 = p * q


      fout = f(i, j-1) * a0 + f(i-1, j) * a1 + f(i,j) * a2
     &     + f(i+1, j) * a3 + f(i, j+1) * a4 + f(i+1, j+1) * a5

      return

  900 format(1p,6e12.4)
  910 format(1p,8e12.4)
  920 format(2i10)
      end
      
c
c***************************************************************************
c

      subroutine fluxavg2(f, favg, nxdim, nydim, nrhodim,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, r0,
     &   x_fine, y_fine, rho_fine,
     &   nnodex_fine, nnodey_fine, capr_fine, dx_fine, dy_fine)

      implicit none

      integer nxdim, nydim, nrhodim, nnodex, nnodey, nnoderho
      integer n, i, j, nnodex_fine, nnodey_fine

      real f(nxdim, nydim), favg(nrhodim), r0, drho, dx, dy, pi, twopi

      real x_fine(nnodex_fine), y_fine(nnodey_fine),
     &   rho_fine(nnodex_fine, nnodey_fine),
     &   capr_fine(nnodex_fine), dx_fine, dy_fine

      real, dimension(:,:),   allocatable :: f_fine
      real, dimension(:),     allocatable :: vol
      real, dimension(:),     allocatable :: fvol

      allocate(f_fine(nnodex_fine, nnodey_fine) )
      allocate(fvol(nrhodim) )
      allocate( vol(nrhodim) )

*     --------------------------------
*     Interpolate f onto fine 2-D grid
*     --------------------------------

      do i = 1, nnodex_fine
         do j = 1, nnodey_fine

            call intplt2(x_fine(i), y_fine(j), f_fine(i,j),
     &         nnodex, nnodey, f, nxdim, nydim, dx, dy)

         end do
      end do
      
      pi = 3.141592654
      twopi = 2.0 * pi

      fvol = 0.0
      vol = 0.0
      favg = 0.0
             
      do i = 1, nnodex_fine
         do j = 1, nnodey_fine
            n = int(rho_fine(i,j) / drho) + 1
            if(n .le. nnoderho)then
               fvol(n) = fvol(n) + dx_fine * dy_fine 
     &                              * twopi * capr_fine(i) * f_fine(i,j)
               vol(n) =  vol(n) + dx_fine * dy_fine *twopi *capr_fine(i)
            end if
         end do
      end do


      do n = 1, nnoderho 
         favg(n) = 0.0
         if(vol(n) .ne. 0.0)favg(n) = fvol(n) / vol(n)
      end do
      
      do n = 2, nnoderho - 1
         if(favg(n) .eq. 0.0)then
            favg(n) = (favg(n-1) + favg(n+1)) / 2.0
         end if
      end do

      if (favg(1) .eq. 0.0) favg(1) = favg(2)
      if (favg(nnoderho) .eq. 0.0) favg(nnoderho) = favg(nnoderho - 1)

      deallocate(f_fine)
      deallocate(fvol)
      deallocate(vol)

  100 format (1i10, 1p,8e12.4)
  102 format (2i10)
  
      return
      end

c
c***************************************************************************
c


