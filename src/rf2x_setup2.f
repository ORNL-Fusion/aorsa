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

      integer::  nxdim, nydim, myid
      integer n

      integer nxmx, nymx, ieq

      integer mkdim1, mkdim2

      integer nkdim1, nkdim2, nkx1, nkx2, nldim, nldim3,
     1   nky1, nky2, iant

      integer nrhomax
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



      real rhomax

      real gaussian_ne, gaussian_te, gaussian_ti1, gaussian_ti2,
     &   gaussian_ti3


      integer nmodesx, nmodesy,
     &    inu, iprint, iexact,
     &    nstep, nabs,
     &    isigma, itemp,
     &    nfreqm,  nkzm,
     &    idens,  iabsorb,
     &    nzfun, nnodex, nnodey, i, j,
     &    iflag, liw, lw, nrhs, icenter

      integer nnoderho

      integer m, nphi, iflag_gammab

      real xant, te0,
     &    xwall, xnwall,
     &    epszet, amu1, amu2, z1, z2, eta,
     &    b0, rt, ytop, ybottom, aplasm, xnlim, 
     &    xn0, flat, b1rat, b2rat, curdnx, curdny, curdnz,
     &    xbnch, 
     &    rhoplasm,
     &    dfreq,  dkz, reomg1, reomg2, reomg3,
     &    r0, adip, efold,
     &    amu3, z3, eta3,
     &    q0, prfin,
     &    alim, grad, qavg0, ymax,
     &    ekappa, xiota0, rholim, psilim, psimol, yant,
     &    rwleft, rwright, xwleft, psi_lim, psi1

      real xktau, xkrho, rant, xnphi
      real rlim

      real domgk, 
     &   vthi10, vthi, rnz, xn2, eta2,
     &   eta1, xmh, xme, 
     &   clight, xmu0, xlnlam,
     &   omgci10, xmax, qe,i3, pi, eps0, xn3, qi3,
     &   costh, sinth, radius, rnx,  rny, rnphi, twopi

      real capr(nxmx), xprime(nxmx), x(nxmx), dx, dxc


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

      real fyavg(nrhomax), qhat, omgte, omgti

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
      real rhomavg(nrhomax)

      real drhodx, drhodxx, drhody, drhodyy, gradrho

      real dxpsi, dxxpsi, dxypsi, dypsi, dyypsi
        

      real rho(nxdim, nydim), theta(nxmx, nymx),
     &     rhohatx(nxmx, nymx), rhohaty(nxmx, nymx),
     &     theta0(nxmx, nymx),
     &     bx(nxmx, nymx), by(nxmx, nymx), bz(nxmx, nymx),
     &     btau(nxmx, nymx), bzeta(nxmx, nymx),
     &     dxdth(nxmx, nymx), dzdth(nxmx, nymx), xntau(nxmx, nymx),
     &     xn(nxmx, nymx), xkte(nxmx, nymx),
     &     omgce(nxmx, nymx),
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

      real betan, betan2, betan3, betate, betati, betati2,
     &    betati3, taue

      real yprimec(nymx), ycourse(nymx),
     &     yprime(nymx), y(nymx), dy, dyc


!------------------end of declarations----------------------------------

*     ---------------------------------
*     set default values of input data:
*     ---------------------------------

      if (allocated(capr_fine)) then
         deallocate(capr_fine)
      end if
      nnoderho = nmodesx/2 !50  !JCW bad magic number

      rant = 1.6

      psilim = .99
      psimol = 0.90

      yant = .1

      inu = 0
      iprint = 50
      iexact = 1
      xwall = 0.0000E+00

      epszet = 1.0000E-07
      nphi = 33
      amu1 = 2.0000E+00
      amu2 = 1.0000E+00
      z1 = 1.0000E+00
      z2 = 1.0000E+00
      eta = 4.5000E-01
      q0 = 1.0
      ekappa = 2.0

      ymax = 0.0
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
      xbnch = 0.0000E+00

      isigma = 1
      nzfun = 1
      qavg0 = 1.0
      iabsorb = 2
      itemp = 0

      xnlim  = 0.0000E+00

      nfreqm = 1
      dfreq = 0.0000E+00
      nkzm = 1
      dkz = 0.0000E+00
      idens = 0
      r0 = 2.3000E+00
      adip = 0.0000E+00
      efold = 0.0000E+00
      amu3 = 1.2000E+01
      z3 = 6.0000E+00
      eta3 = 4.6600E-02


      nrhs = 1

      nkx2 = nmodesx / 2
      nkx1 = - nmodesx / 2 + 1
      nnodex = nmodesx


      nnodey = nmodesy

      icenter = nnodex / 2

      if (qavg0 .ne. 0.0) xiota0 = 1./qavg0

      rholim = sqrt(psilim)


      qhat = qavg0


      eps0 = 8.85e-12
      pi = 3.141592654
      twopi = 2. * pi
      xlnlam = 20.0
      xmu0 = 1.26e-06
      clight = 1.0 / sqrt(eps0 * xmu0)
      xmax = rwright - rwleft
      ymax = ytop - ybottom

      xn2 = xn0 * eta
      eta2 = xn2 / xn0
      xn3 = xn0 * eta3

      xwleft = rwleft - rt

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

      end do


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
#ifdef DEBUG
      write(*,*) 'rhomax',rhomax,drho
#endif
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
     

  310 format(1p,6e12.4)
  309 format(10i10)
  311 format(1p,10e12.4)



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
            if (n .le. 0) then
               write(*,*) 'fvol',i,j,drho ,rho_fine(i,j),n
            end if
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


