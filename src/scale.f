
      subroutine disp_ono(yplot)
    
      use size_mod

      implicit none

      integer nmax, npoints, nline, i, nsum, n, kmax, kn, jmid
      integer pgopen, pgbeg, ier, nlevmax, iflag, numb
      integer nxmx, nymx, nnodex, nnodey, j, jplot

      parameter (nxmx = nmodesmax)
      parameter (nymx = mmodesmax)
      parameter (nmax = nmodesmax)

      parameter (nlevmax = 101)

      real xgrid(nmax), ygrid(nmax), f(nmax), dx, xmin, xmax, capf(nmax)
      real capr(nmax), ff(101), rhoplasm, dy, yplot, ybottom, yprime
      complex xkperp2_slow(nxmx, nymx), xkperp2_fast(nxmx, nymx)
      complex xkperp_slow(nxmx, nymx), xkperp_fast(nxmx, nymx)
      complex P(nxmx, nymx)
      real xkprl(nxmx, nymx)

      complex xkperp2_slow_plot(nxmx), xkperp2_fast_plot(nxmx)
      complex xkperp_slow_plot(nxmx), xkperp_fast_plot(nxmx)

      real rho(nxmx, nymx)

      real f3(nmax)
      real y(nmax), k, capk, c, xdeg, yh, L, kr, nr
      real factrl, xn, an, diff1, diff2

      real ans_simpson, ans_trapezoidal, ans_analytic, ymax, ymin
      real dummy, time, t1, t2, second1, pi, x, a, e, dydx, sum

      character(32):: title
      character(32):: titlx
      character(32):: tityl
      character(32):: tityr

      character(32):: titx
      character(32):: tity
      character(32):: titz

      integer nblack, nred, nyellow, ngreen, nblue, ncyan, nmagenta,
     &   nwhite, norange

      common/boundcom/rhoplasm

      open(unit=130, file='Ono_disp', status='unknown',
     &      form='formatted')

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      numb = 100


      pi = 3.14159
      e =  2.71828
      a = 3.0
      k = .0015
      capk = 6000.
      f3 = 0.

      read(130, 309) nnodex, nnodey, jmid
      read(130, 310) rhoplasm
      read(130, 310) (capr(i), i = 1, nnodex)
      read(130, 310) (y(j), j = 1, nnodey)
      read(130, 310) ((xkperp2_slow(i, j), i = 1, nnodex),
     &      j = 1, nnodey)
      read(130, 310) ((xkperp2_fast(i, j), i = 1, nnodex),
     &      j = 1, nnodey)
      read(130, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      read(130, 310) ((xkprl(i,j), i = 1, nnodex), j = 1, nnodey)
      read(130, 310) ((P(i, j), i = 1, nnodex), j = 1, nnodey)

      xgrid = capr
      ygrid = y
      npoints = nnodex

      xmin = xgrid(1)
      xmax = xgrid(nnodex)

      jplot = jmid

      dy = y(2) - y(1)
      ybottom = y(1)
c      yplot = 0.05
      yprime = yplot - ybottom
      jplot = yprime / dy + 1
      write(6, *) "jplot = ", jplot
      write(6, *) "yplot = ", yplot


!     -------------------------------
!     Open the pgplot graphics device
!     -------------------------------

#ifdef GIZA
      IER = PGBEG(0, 'disper.pdf', 1, 1)
#else
      IER = PGBEG(0, 'disper.ps/vcps', 1, 1)
#endif

      IF (IER.NE.1) STOP

      call PGSCH (1.5)


!     -------------------------------------------
!     Plot contours of k**2 - slow root
!     -------------------------------------------
      titx = 'R (m)'
      tity = 'Z (m)'
      title = 'Re(k**2) - slow root'
      call ezconc(capr, y, real(xkperp2_slow), ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &   nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity)


      title = 'k_{parallel}'
      call ezconc(capr, y, xkprl, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &   nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity)


      title = 'Re(P)'
      call ezconc(capr, y, real(P), ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &   nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity)

      title = 'Im(P)'
      call ezconc(capr, y, aimag(P), ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &   nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity)

      title = 'Im(k**2) - slow root'
      call ezconc(capr, y, aimag(xkperp2_slow), ff, nnodex, nnodey,numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &   nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity)

!     ------------------------
!     Plot cut in y at yplot
!     ------------------------
      do i = 1, nnodex
         xkperp2_slow_plot(i) = xkperp2_slow(i, jplot)
         xkperp2_fast_plot(i) = xkperp2_fast(i, jplot)
      end do

      xkperp_slow_plot = csqrt(xkperp2_slow_plot)
      xkperp_fast_plot = csqrt(xkperp2_fast_plot)

      titlx  = 'R (m)'
      tityl = 'kperp2 (m-2)'
      tityr = ' '
      title = 'fast and slow roots'

      call ezplot2_s(title, titlx, tityl, tityr, xgrid,
     &   real(xkperp2_slow_plot), real(xkperp2_fast_plot), f3,
     &   npoints, nmax, ymax, ymin, xmin, xmax)

      tityl = 'kperp (m-1)'

      call ezplot2_b(title, titlx, tityl, tityr, xgrid,
     &   real(xkperp_slow_plot), real(xkperp_fast_plot), f3,
     &   npoints, nmax, ymax, ymin, xmin, xmax)

c      title = 'slow root'
c     call ezplot1_s(title, titlx, tityl, tityr, xgrid,
c     &   real(xkperp2_slow_plot), f3,
c     &   npoints, nmax, ymax, ymin, xmin, xmax, nred)

      tityl = 'kperp2 (m-2)'
      title = 'fast root'
      call ezplot1_s(title, titlx, tityl, tityr, xgrid,
     &   real(xkperp2_fast_plot), f3,
     &   npoints, nmax, ymax, ymin, xmin, xmax, nblue)

!     -------------------------------------------
!     Plot contours of k**2 - fast root
!     -------------------------------------------

      title = 'Re(k**2) - fast root'
      call ezconc(capr, y, real(xkperp2_fast), ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &   nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity)

      title = 'Im(k**2) - fast root'
      call ezconc(capr, y, aimag(xkperp2_fast), ff, nnodex, nnodey,numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &   nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity)


      xkperp_slow = csqrt(xkperp2_slow)
      xkperp_fast = csqrt(xkperp2_fast)

      do i = 1, nnodex
         xkperp_slow_plot(i) = xkperp_slow(i, jplot)
         xkperp_fast_plot(i) = xkperp_fast(i, jplot)
      end do

!     -------------------------------------------------
!     write the function xkperp_fast_plot(i) to screen
!     -------------------------------------------------
      do i = 1, npoints
         write(6, '(i10, 1p,5e15.5)') i, xgrid(i), xkperp_fast_plot(i)
      end do

      title = 'Re(k) - fast and slow roots'
      titlx  = 'R (m)'
      tityl = 'Re(k) (m-1)'
      tityr = ' '

      call ezplot2_sym(title, titlx, tityl, tityr, xgrid,
     &   real(xkperp_slow_plot), real(xkperp_fast_plot), f3,
     &   npoints, nmax, 0.0, 0.0, xmin, xmax)



      title = 'Re(k) - fast and slow (blowup)'
      titlx  = 'R (m)'
      tityl = 'Re(k) (m-1)'
      tityr = ' '

      call ezplot2_sym(title, titlx, tityl, tityr, xgrid,
     &   real(xkperp_slow_plot), real(xkperp_fast_plot), f3,
     &   npoints, nmax, 2000., -2000., xmin, xmax)



      title = 'Re(k) - slow root'
      call ezconc(capr, y, real(xkperp_slow), ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &   nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity)



!     -------------------------
!     Close the graphics device
!     -------------------------
      call pgclos




      close (130)

  310 format(1p,6e12.4)
 3310 format(1p,6e18.10)
 8310 format(1p,6e14.6)
  309 format(10i10)

      return
      end
!
!*******************************************************************************
!


      subroutine ezplot1_s(title, titx, titll, titlr, x1, y1, y3,
     &   nr, nrmax,
     &   ymax, ymin, xmin, xmax, ncolor)

      implicit none

      integer nr, nrmax
      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y3(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin, ymax, xmin, xmax

      character(32):: title
      character(32):: titx
      character(32):: titll
      character(32):: titlr


      integer nplot1,ncollab, ncolion,ncolbox, ncyan, ncolelec,
     &   ncolln2, ncollin, ncolbrd

      integer nblack, nred, nyellow, ngreen, naqua, npink, nwheat,
     &   ngrey, nbrown, nblue, nblueviolet, ncyan1, ncolor
      integer nturquoise, nmagenta, nsalmon, nwhite, norange, ncolln3


      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      ncolln3 = ngreen

      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = xmax
      rhomin = xmin


! Advance plotter to a new page, define coordinate range of graph and draw axes


      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX   ('BCNST', 0.0, 0,  'BCNST', 0.0, 0)

      CALL PGSCI(nblack)


! Label the axes (note use of \u and \d for raising exponent).

      call pglab(titx, titll, title)

! Plot the line graph.

      CALL PGSCI(ncolor)
      call pgline(nr, x1, y1)

  300 format (1p,9e11.3)

      CALL PGSCI(nblack)
      call pgline(nr, x1, y3)

      return
      end



!
!********************************************************************
!
      subroutine ezplot2_s(title, titlb, titll, titlr, x1, y1, y2, y3,
     &   nr, nrmax, ymax, ymin, xmin, xmax)

      implicit none

      integer nr, nrmax, n

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax), y2(nrmax), y3(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax, xmin, xmax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, ncolelec,
     &   ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, nwheat,
     &   ngrey,nbrown,nblue,nblueviolet,ncyan1
      integer nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)

      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = xmax
      rhomin = xmin


      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


      CALL PGSCI(nblack)


! Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)

      CALL PGSCI(nred)
      call pgline(nr, x1, y1)


      CALL PGSCI(nblue)
      call pgline(nr, x1, y2)


      CALL PGSCI(nblack)
      call pgline(nr, x1, y3)


  300 format (1p,9e11.3)
      return
      end


!
!********************************************************************
!
      subroutine ezplot2_b(title, titlb, titll, titlr, x1, y1, y2, y3,
     &   nr, nrmax, ymax, ymin, xmin, xmax)

      implicit none

      integer nr, nrmax, n

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax), y2(nrmax), y3(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax, xmin, xmax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, ncolelec,
     &   ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, nwheat,
     &   ngrey,nbrown,nblue,nblueviolet,ncyan1
      integer nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y2, nrmax, nr, y2min, y2max)

      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y2max * 2.0
      ymin = y2min

      rhomax = xmax
      rhomin = xmin


      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


      CALL PGSCI(nblack)


! Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)

      CALL PGSCI(nred)
      call pgline(nr, x1, y1)


      CALL PGSCI(nblue)
      call pgline(nr, x1, y2)


      CALL PGSCI(nblack)
      call pgline(nr, x1, y3)


  300 format (1p,9e11.3)
      return
      end

!
!*****************************************************************
!
      subroutine ezplot2_sym(title, titlb, titll, titlr, x1, y1, y2, y3,
     &   nr, nrmax, ymax_giv, ymin_giv, xmin, xmax)

      implicit none

      integer nr, nrmax, n
      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax), y2(nrmax), y3(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax, xmin, xmax, ymax_giv, ymin_giv

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, ncolelec,
     &   ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, nwheat,
     &   ngrey,nbrown,nblue,nblueviolet,ncyan1
      integer nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)

!      write(6, 300) y2min, y2max
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max * 1.1
      ymin = y1min

      if(ymax_giv .ne. 0.0)ymax = ymax_giv
      ymin = -ymax

      rhomax = xmax
      rhomin = xmin


      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


      CALL PGSCI(nblack)


! Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)

      CALL PGSCI(nred)
      call pgline(nr, x1, y1)
      call pgline(nr, x1, -y1)


      CALL PGSCI(nblue)
      call pgline(nr, x1, y2)
      call pgline(nr, x1, -y2)


      CALL PGSCI(nblack)
      call pgline(nr, x1, y3)


  300 format (1p,9e11.3)
      return
      end

!
!********************************************************************
!
