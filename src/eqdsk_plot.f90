!
      subroutine trace_plot
      use size_mod

      implicit none


      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titx
      character(32):: tity
      character(32):: titz
      character(32):: titlb

      real:: logmax
      integer:: ibackground, nkx1, nkx2, nky1, nky2, n, m, &
         nkxplt, nkyplt

      integer::  ncolln10, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, ncollin, ncollab, ncolion, &
         ncolelec, norange

      integer:: nnoderho, iflag, nnoderho_half, nnoderho2_half

      integer:: pgopen, pgbeg, ier

      integer:: iostatval

      integer, parameter::  nlevmax = 101, nxmx = nmodesmax, &
           nymx = mmodesmax,      nrhomax = nxmx, &
           nkdim1 = - nxmx / 2,   nkdim2 =   nxmx / 2, &
           mkdim1 = - nxmx / 2,   mkdim2 =   nxmx / 2, &
           nkpltdim = 2 * nkdim2, mkpltdim = 2 * mkdim2, &
           nxeqdmax = 257,        nyeqdmax = 257

      real:: capr_x(6000), capz_x(6000), dx, dy
      integer:: n_phi, n_phi_max, number_points


      real rhoeqdsk(nxeqdmax), psigrid(nxeqdmax), agrid(nxeqdmax), &
         fpsi(nxeqdmax), dfpsida(nxeqdmax), &
         qpsi(nxeqdmax), dqpsida(nxeqdmax)

      real xkxsav(nkpltdim), xkysav(mkpltdim), pscale
      real rhon(nrhomax), vyi1avg(nrhomax),  vyi2avg(nrhomax), &
           rhon_half(nrhomax)
      real capr_bpol_midavg(nrhomax), bmod_midavg(nrhomax), &
           dldbavg(nrhomax)



      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
         ncolion,ncolelec,ncollin,ncollab,ncolln2
      common/boundcom/rhoplasm
      common/zoom/ xmaxz, xminz, ymaxz

      real xmaxz, xminz, ymaxz
      real fmin, fmax, fminre, fmaxre, fminim, fmaxim, fmin1, &
         fmax1, fmax2, fmax3, fmaxt, fmin2, fmin3, fmint
      real x(nxmx), capr(nxmx), y(nymx)

      real, dimension(:,:), allocatable :: xkti, xkte, xn, &
         rho, xjy, bmod, bmod_mid, capr_bpol_mid2, &
         xjx, psi, xiota, qsafety, btau, bzeta, &
         dbxdx, dbydx, dbzdx, dbxdy, dbydy, dbzdy, dbdx, dbdy, &
         gradprlb, dxxbxn, dxxbyn, dxxbzn, &
         dxybxn, dxybyn, dxybzn, dyybxn, dyybyn, dyybzn, &
         dxxmodb, dxymodb, dyymodb, bx, by, bz, bz_dum, dldb_tot12, &
         freal, fimag, fmod, rho_tor2d


      real xnmid(nxmx), xktimid(nxmx),  xktemid(nxmx), qmid(nxmx)
      real bpmid(nxmx), xiotamid(nxmx)
      real fmidre(nxmx), fmidim(nxmx), fmid1(nxmx), &
           fmid2(nxmx), fmid3(nxmx), fmidt(nxmx)
      real ff(101)
      real q, omgrf, xk0, n0, clight, xmu0, eps0, rhoplasm
      real temax, temin, timin, tmin, tmax, timax, caprmaxp, &
         caprmin, caprminp, caprmax, xnmax, xnmin, qmin, qmax
      real bpmin, bpmax
      integer:: nnodex, j, i, nnodey, numb, jmid, nxeqd

      namelist/plotin/ibackground, xminz, xmaxz, ymaxz, logmax

      allocate ( xkti(nxmx, nymx), xkte(nxmx, nymx), xn(nxmx, nymx), &
         rho(nxmx, nymx), &
         xjy(nxmx, nymx), bmod(nxmx, nymx), &
         bmod_mid(nxmx, nymx), capr_bpol_mid2(nxmx, nymx), &
         xjx(nxmx, nymx), psi(nxmx, nymx), &
         xiota(nxmx, nymx), qsafety(nxmx, nymx), &
         btau(nxmx, nymx), bzeta(nxmx, nymx), &
         dbxdx(nxmx, nymx), dbydx(nxmx, nymx), dbzdx(nxmx, nymx), &
         dbxdy(nxmx, nymx), dbydy(nxmx, nymx), dbzdy(nxmx, nymx), &
         dbdx(nxmx, nymx),  dbdy(nxmx, nymx), &
         gradprlb(nxmx, nymx), &
         dxxbxn(nxmx, nymx), dxxbyn(nxmx, nymx), dxxbzn(nxmx, nymx), &
         dxybxn(nxmx, nymx), dxybyn(nxmx, nymx), dxybzn(nxmx, nymx), &
         dyybxn(nxmx, nymx), dyybyn(nxmx, nymx), dyybzn(nxmx, nymx), &
         dxxmodb(nxmx, nymx), dxymodb(nxmx, nymx), dyymodb(nxmx,nymx), &
         bx(nxmx, nymx), by(nxmx, nymx), bz(nxmx, nymx), &
         bz_dum(nxmx, nymx), dldb_tot12(nxmx, nymx), &
         freal(nxmx, nymx), fimag(nxmx, nymx), fmod(nxmx, nymx), &
         rho_tor2d(nxmx, nymx) )


!--set default values of input data:
      ibackground = 1
      xminz = 0.51
      xmaxz = 0.62
      ymaxz = 0.24
      logmax = 2.0

! <<  input data
!      read (63, plotin)

      open(unit=116,file='trace',status='old',form='formatted')


!      eps0 = 8.85e-12
!      xmu0 = 1.26e-06
!      clight = 1./sqrt(eps0 * xmu0)
!      xk0 = omgrf / clight
      q = 1.6e-19


!--For white background set ibackground = 0
!      ibackground = 0

!--For black background set ibackground = 1
!      ibackground = 1

!      read(63, 309)ibackground


      read(116, 309) nnodex, nnodey
      read(116, 310) rhoplasm
      read(116, 310) (x(i), i = 1, nnodex)
      read(116, 310) (y(j), j = 1, nnodey)
      read(116, 310) (capr(i), i = 1, nnodex)

      read(116, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)

      read(116, 309) n_phi_max
      read(116, 310) (capr_x(n_phi), n_phi = 1, n_phi_max)
      read(116, 310) (capz_x(n_phi), n_phi = 1, n_phi_max)


      close (116)


      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      ncollin=ncyan
      ncolln2=nblueviolet
      ncolln3=ngreen
      ncolln4=naqua
      ncolln5=nyellow
      ncolbrd=nwheat
      ncollab=ncyan

      ncolbox=nred
      ncolbrd=nwheat
!     ncolbox=nblack
!     ncolbrd=nblack

!     ncolion=nblue
!     ncolelec=nred
      ncolion=naqua
      ncolelec=nyellow

      if (ibackground.eq.1)then
         ncolbox=nred
         ncolbrd=nwheat
      end if

      if (ibackground.eq.0)then
         ncolbox=nblack
         ncolbrd=nblack
      end if

! Open graphics device

#ifdef GIZA
      IER = PGBEG(0, 'trace.pdf', 1, 1)
#else
      IER = PGBEG(0, 'trace.ps/vcps', 1, 1)
#endif

      IF (IER.NE.1) STOP
      call PGSCF(2)
      call PGSCH (1.5)


!--Contour plots using 16 contour levels:

      if (ibackground.eq.1)then
!        black background
         ncollin = ncyan
         ncolln2 = nblueviolet
      end if

      if (ibackground.eq.0)then
!        white background
         ncollin = nblue
         ncolln2 = nyellow
      end if

      titx = 'R (m)'
      tity = 'Y (m)'

      numb = 20
      title = 'poloidal rho surfaces'
      call ezconc(capr, y, rho, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)


      CALL PGSCI(nred)
      call pgline(n_phi_max, capr_x, capz_x)


! Close the graphics device.

      call pgclos

      deallocate ( xkti, xkte, xn, &
         rho, xjy, bmod, bmod_mid, capr_bpol_mid2, &
         xjx, psi, xiota, qsafety, btau, bzeta, &
         dbxdx, dbydx, dbzdx, dbxdy, dbydy, dbzdy, dbdx, dbdy, &
         gradprlb, dxxbxn, dxxbyn, dxxbzn, &
         dxybxn, dxybyn, dxybzn, dyybxn, dyybyn, dyybzn, &
         dxxmodb, dxymodb, dyymodb, bx, by, bz, bz_dum, dldb_tot12, &
         freal, fimag, fmod, rho_tor2d )





  310 format(1p6e12.4)
  309 format(10i10)
  311 format(1p10e12.4)
 1312 format(i10,1p8e12.4)
 1313 format(2i10,1p8e12.4)
 9310 format(1p7e12.4)
   10 format(i10,1p4e10.3,i10,1pe10.3)

      return
      end
!
!********************************************************************
!

!
      subroutine eqdsk_plot
      use size_mod

      implicit none


      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titx
      character(32):: tity
      character(32):: titz
      character(32):: titlb

      real:: logmax
      integer:: ibackground, nkx1, nkx2, nky1, nky2, n, m, &
         nkxplt, nkyplt

      integer:: iostatval

      integer::  ncolln10, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, ncollin, ncollab, ncolion, &
         ncolelec, norange

      integer:: nnoderho, iflag, nnoderho_half, nnoderho2_half

      integer:: pgopen, pgbeg, ier

      integer, parameter::  nlevmax = 101, nxmx = nmodesmax, &
           nymx = mmodesmax,      nrhomax = nxmx, &
           nkdim1 = - nxmx / 2,   nkdim2 =   nxmx / 2, &
           mkdim1 = - nxmx / 2,   mkdim2 =   nxmx / 2, &
           nkpltdim = 2 * nkdim2, mkpltdim = 2 * mkdim2, &
           nxeqdmax = 257,        nyeqdmax = 257

      real:: capr_x(6000), capz_x(6000), dx, dy
      integer:: n_phi, n_phi_max, number_points


      real rhoeqdsk(nxeqdmax), psigrid(nxeqdmax), agrid(nxeqdmax), &
         fpsi(nxeqdmax), dfpsida(nxeqdmax), &
         qpsi(nxeqdmax), dqpsida(nxeqdmax)

      real xkxsav(nkpltdim), xkysav(mkpltdim), pscale
      real rhon(nrhomax), vyi1avg(nrhomax),  vyi2avg(nrhomax), &
           rhon_half(nrhomax)
      real capr_bpol_midavg(nrhomax), bmod_midavg(nrhomax), &
           dldbavg(nrhomax)



      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
         ncolion,ncolelec,ncollin,ncollab,ncolln2
      common/boundcom/rhoplasm
      common/zoom/ xmaxz, xminz, ymaxz

      real xmaxz, xminz, ymaxz
      real fmin, fmax, fminre, fmaxre, fminim, fmaxim, fmin1, &
         fmax1, fmax2, fmax3, fmaxt, fmin2, fmin3, fmint
      real x(nxmx), capr(nxmx), y(nymx)



!      real xkti(nxmx, nymx), xkte(nxmx, nymx), xn(nxmx, nymx)
!      real rho(nxmx, nymx)
!      real xjy(nxmx, nymx), bmod(nxmx, nymx)
!      real bmod_mid(nxmx, nymx), capr_bpol_mid2(nxmx, nymx)
!      real xjx(nxmx, nymx), psi(nxmx, nymx)
!      real xiota(nxmx, nymx), qsafety(nxmx, nymx)
!      real btau(nxmx, nymx), bzeta(nxmx, nymx)
!      real dbxdx(nxmx, nymx), dbydx(nxmx, nymx), dbzdx(nxmx, nymx),
!     .     dbxdy(nxmx, nymx), dbydy(nxmx, nymx), dbzdy(nxmx, nymx),
!     .     dbdx(nxmx, nymx),  dbdy(nxmx, nymx)
!      real gradprlb(nxmx, nymx)
!      real dxxbxn(nxmx, nymx), dxxbyn(nxmx, nymx), dxxbzn(nxmx, nymx),
!     .     dxybxn(nxmx, nymx), dxybyn(nxmx, nymx), dxybzn(nxmx, nymx),
!     .     dyybxn(nxmx, nymx), dyybyn(nxmx, nymx), dyybzn(nxmx, nymx),
!     .     dxxmodb(nxmx, nymx), dxymodb(nxmx, nymx), dyymodb(nxmx,nymx),
!     .     bx(nxmx, nymx), by(nxmx, nymx), bz(nxmx, nymx),
!     .     dldb_tot12(nxmx, nymx)
!      real freal(nxmx, nymx), fimag(nxmx, nymx), fmod(nxmx, nymx)
!      real rho_tor2d(nxmx, nymx)

      real, dimension(:,:), allocatable :: xkti, xkte, xn, &
         rho, xjy, bmod, bmod_mid, capr_bpol_mid2, &
         xjx, psi, xiota, qsafety, btau, bzeta, &
         dbxdx, dbydx, dbzdx, dbxdy, dbydy, dbzdy, dbdx, dbdy, &
         gradprlb, dxxbxn, dxxbyn, dxxbzn, &
         dxybxn, dxybyn, dxybzn, dyybxn, dyybyn, dyybzn, &
         dxxmodb, dxymodb, dyymodb, bx, by, bz, bz_dum, dldb_tot12, &
         freal, fimag, fmod, rho_tor2d




      real xnmid(nxmx), xktimid(nxmx),  xktemid(nxmx), qmid(nxmx)
      real bpmid(nxmx), xiotamid(nxmx)
      real fmidre(nxmx), fmidim(nxmx), fmid1(nxmx), &
           fmid2(nxmx), fmid3(nxmx), fmidt(nxmx)
      real ff(101)
      real q, omgrf, xk0, n0, clight, xmu0, eps0, rhoplasm
      real temax, temin, timin, tmin, tmax, timax, caprmaxp, &
         caprmin, caprminp, caprmax, xnmax, xnmin, qmin, qmax
      real bpmin, bpmax
      integer:: nnodex, j, i, nnodey, numb, jmid, nxeqd
!      complex zi

      namelist/plotin/ibackground, xminz, xmaxz, ymaxz, logmax

      allocate ( xkti(nxmx, nymx), xkte(nxmx, nymx), xn(nxmx, nymx), &
         rho(nxmx, nymx), &
         xjy(nxmx, nymx), bmod(nxmx, nymx), &
         bmod_mid(nxmx, nymx), capr_bpol_mid2(nxmx, nymx), &
         xjx(nxmx, nymx), psi(nxmx, nymx), &
         xiota(nxmx, nymx), qsafety(nxmx, nymx), &
         btau(nxmx, nymx), bzeta(nxmx, nymx), &
         dbxdx(nxmx, nymx), dbydx(nxmx, nymx), dbzdx(nxmx, nymx), &
         dbxdy(nxmx, nymx), dbydy(nxmx, nymx), dbzdy(nxmx, nymx), &
         dbdx(nxmx, nymx),  dbdy(nxmx, nymx), &
         gradprlb(nxmx, nymx), &
         dxxbxn(nxmx, nymx), dxxbyn(nxmx, nymx), dxxbzn(nxmx, nymx), &
         dxybxn(nxmx, nymx), dxybyn(nxmx, nymx), dxybzn(nxmx, nymx), &
         dyybxn(nxmx, nymx), dyybyn(nxmx, nymx), dyybzn(nxmx, nymx), &
         dxxmodb(nxmx, nymx), dxymodb(nxmx, nymx), dyymodb(nxmx,nymx), &
         bx(nxmx, nymx), by(nxmx, nymx), bz(nxmx, nymx), &
         bz_dum(nxmx, nymx), dldb_tot12(nxmx, nymx), &
         freal(nxmx, nymx), fimag(nxmx, nymx), fmod(nxmx, nymx), &
         rho_tor2d(nxmx, nymx) )

!--set default values of input data:
      ibackground = 1
      xminz = 0.51
      xmaxz = 0.62
      ymaxz = 0.24
      logmax = 2.0

! <<  input data
!      read (63, plotin)

      open(unit=138,file='out138',status='old',form='formatted')
      open(unit=638,file='fpm2',status ='unknown', form='formatted')

      open(unit=145,file='Bfield_2D.vtk',status='unknown', &
                                             form='formatted')

!      open(unit=53,file='out53',status='unknown',form='formatted')
!      open(unit=58,file='out58',status='unknown',form='formatted')
!      open(unit=52,file='out52',status='unknown',form='formatted')
!      open(unit=59,file='out59',status='unknown',form='formatted')

!      open(unit=51,file='rho',status='unknown',form='formatted')
!      open(unit=54,file='movie_spx',status='unknown',form='formatted')
!      open(unit=55,file='movie_spy',status='unknown',form='formatted')
!      open(unit=56,file='movie_spz',status='unknown',form='formatted')


!      open(unit=57,file='emodes',status='unknown',form='formatted')

!      open(unit=60,file='test',status='unknown',form='formatted')



!      zi = cmplx(0.0, 1.0)
!      eps0 = 8.85e-12
!      xmu0 = 1.26e-06
!      clight = 1./sqrt(eps0 * xmu0)
!      xk0 = omgrf / clight
      q = 1.6e-19


!--For white background set ibackground = 0
!      ibackground = 0

!--For black background set ibackground = 1
!      ibackground = 1

!      read(63, 309)ibackground

      read(138, 309) nnodex, nnodey, nxeqd
      read(138, 310) rhoplasm

      ! assert that nxeqdmax >= nxeqd

      if(nxeqd.gt.nxeqdmax)then
          write(*,*) 'ERROR'
          write(*,*) 'Value of nxeqdmax is too small in eqdsk_plot.f90'
          write(*,*) 'nxeqdmax = ', nxeqdmax
          write(*,*) 'nxeq = ', nxeqd
          stop 1
      end if

      read(138, 310) (x(i), i = 1, nnodex)
      read(138, 310) (y(j), j = 1, nnodey)
      read(138, 310) (capr(i), i = 1, nnodex)


      read(138, 310) ((bmod(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((psi(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((qsafety(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((btau(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((bzeta(i, j), i = 1, nnodex), j = 1, nnodey)

      read(138, 310) ((dbxdx(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dbydx(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dbzdx(i, j), i = 1, nnodex), j = 1, nnodey)

      read(138, 310) ((dbxdy(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dbydy(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dbzdy(i, j), i = 1, nnodex), j = 1, nnodey)

      read(138, 310) ((dbdx(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dbdy(i, j), i = 1, nnodex), j = 1, nnodey)

      write(*,*) nxeqd, nxeqdmax
      read(138, 310) (psigrid(i), i = 1, nxeqd)
      read(138, 310, iostat=iostatval) (rhoeqdsk(i), i = 1, nxeqd)
      if(iostatval>0)then
        write(*,*) 'iostat = ', iostatval
      end if

      read(138, 310, iostat=iostatval) (qpsi(i), i = 1, nxeqd)
      if(iostatval>0)then
        write(*,*) 'iostat = ', iostatval
      end if

      read(138, 310) ((dxxbxn(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dxxbyn(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dxxbzn(i, j), i = 1, nnodex), j = 1, nnodey)

      read(138, 310) ((dxybxn(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dxybyn(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dxybzn(i, j), i = 1, nnodex), j = 1, nnodey)

      read(138, 310) ((dyybxn(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dyybyn(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dyybzn(i, j), i = 1, nnodex), j = 1, nnodey)

      read(138, 310) ((dxxmodb(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dxymodb(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((dyymodb(i, j), i = 1, nnodex), j = 1, nnodey)

      read(138, 310) ((gradprlb(i, j), i = 1, nnodex), j = 1, nnodey)

      read(138, 310) ((bmod_mid(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((capr_bpol_mid2(i, j), i = 1,nnodex), j= 1,nnodey)

      read(138, 309) nnoderho
      read(138, 310) (rhon(n), n = 1, nnoderho)
      read(138, 310) (capr_bpol_midavg(n), n = 1, nnoderho)
      read(138, 310) (bmod_midavg(n), n = 1, nnoderho)

      read(138, 310) ((rho_tor2d(i, j), i = 1, nnodex), j = 1, nnodey)

      read(138, 309) n_phi_max
      read(138, 310) (capr_x(n_phi), n_phi = 1, n_phi_max)
      read(138, 310) (capz_x(n_phi), n_phi = 1, n_phi_max)

      read(138, 310) ((dldb_tot12(i, j), i = 1, nnodex), &
                                        j = 1, nnodey)

      read(138, 310) (dldbavg(n), n = 1, nnoderho)

      read(138, 310) ((bx(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((by(i, j), i = 1, nnodex), j = 1, nnodey)
      read(138, 310) ((bz(i, j), i = 1, nnodex), j = 1, nnodey)

      close (138)



      write(638, 309) nnodex, nnodey
      write(638, 310) ((bx(i, j), i = 1, nnodex), j = 1, nnodey)
      write(638, 310) ((by(i, j), i = 1, nnodex), j = 1, nnodey)
      write(638, 310) ((bz(i, j), i = 1, nnodex), j = 1, nnodey)
      close (638)


!     --------------------------------------
!     Write 2D.vtk file "structured points"
!     --------------------------------------
      number_points = nnodex * nnodey
      dx = capr(2) - capr(1)
      dy = y(2) - y(1)

      write(145, 2840)
 2840 format('# vtk DataFile Version 2.0')

      write(145, 2845)
 2845 format('2D vector B field')

      write(145, 2846)
 2846 format('ASCII')

      write(145, 2847)
 2847 format('DATASET STRUCTURED_POINTS')

      write(145, 2841) nnodex,  nnodey, 1
 2841 format('DIMENSIONS', 3i8)

      write(145, 2842) capr(1), y(1), 0
 2842 format('ORIGIN', 2f12.3, 1i8)

      write(145, 2843) dx, dy, 1
 2843 format('SPACING', 2f9.5, 1i8)

      write(145, 2844) number_points
 2844 format('POINT_DATA', i10)

      write(145, 2848)
 2848 format('SCALARS Bx float 1')

      write(145, 2849)
 2849 format('LOOKUP_TABLE default')

      write (145, 3411) ((bx(i,j), i = 1, nnodex),  j = 1, nnodey)

 3410 format(1p4e10.2)
 3411 format(6e17.4)



      write(145, 3853)
 3853 format('SCALARS By float 1')
      write(145, 2849)
      write (145, 3411) ((by(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(145, 6853)
 6853 format('SCALARS Bz float 1')
      write(145, 2849)
      write (145, 3411) ((bz(i,j), i = 1, nnodex),  j = 1, nnodey)



      write(145, 2850)
 2850 format('VECTORS B_vector float')

      bz_dum = 0.0

      write(145, 3412) ((bx(i,j), by(i,j), bz_dum(i,j), i = 1, nnodex), &
                                                        j = 1, nnodey)

 3412 format(3f17.4)


      close (145)


      nnoderho_half = nnoderho - 1

      do n = 1, nnoderho_half
       rhon_half(n) = (rhon(n) + rhon(n+1)) / 2.
      end do


!      write(6, 310) (capr(i), i = 1, nnodex)
!      write(6, 310) ((bmod(i, j), i = 1, nnodex), j = 1, nnodey)

!      write(52,309) nnodex, nnodey
!      write(52,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)



!      write(58,309) nnodex, nnodey
!      write(58,310) ((bmod(i, j), i = 1, nnodex), j = 1, nnodey)


      jmid = nnodey / 2
      do i = 1, nnodex
         xnmid(i) = xn(i, jmid)
         xktemid(i) = xkte(i, jmid) / q
         xktimid(i) = xkti(i, jmid) / q
         qmid(i) = qsafety(i, jmid)
         xiotamid(i) = 1.0 / qmid(i)
         bpmid(i) = btau(i, jmid)
      end do



      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      ncollin=ncyan
      ncolln2=nblueviolet
      ncolln3=ngreen
      ncolln4=naqua
      ncolln5=nyellow
      ncolbrd=nwheat
      ncollab=ncyan

      ncolbox=nred
      ncolbrd=nwheat
!     ncolbox=nblack
!     ncolbrd=nblack

!     ncolion=nblue
!     ncolelec=nred
      ncolion=naqua
      ncolelec=nyellow

      if (ibackground.eq.1)then
         ncolbox=nred
         ncolbrd=nwheat
      end if

      if (ibackground.eq.0)then
         ncolbox=nblack
         ncolbrd=nblack
      end if

! Open graphics device

#ifdef GIZA
      IER = PGBEG(0, 'eqdsk_setup.pdf', 1, 1)
#else
      IER = PGBEG(0, 'eqdsk_setup.ps/vcps', 1, 1)
#endif

      IF (IER.NE.1) STOP
      call PGSCF(2)
      call PGSCH (1.5)


!--Contour plots using 16 contour levels:

      if (ibackground.eq.1)then
!        black background
         ncollin = ncyan
         ncolln2 = nblueviolet
      end if

      if (ibackground.eq.0)then
!        white background
         ncollin = nblue
         ncolln2 = nyellow
      end if

      titx = 'R (m)'
      tity = 'Y (m)'

      numb = 20
      title = 'poloidal rho surfaces'
      call ezconc(capr, y, rho, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)


      CALL PGSCI(nred)
      call pgline(n_phi_max, capr_x, capz_x)


      title = 'toroidal rho surfaces'
      call ezconc(capr, y, rho_tor2d, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)



      title = 'psi surfaces'
!      write(15,*) 'psi test ',psi
!      write(6,*) 'psi test ',psi
      call ezconc(capr, y, psi, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

!      numb = 99

      title = 'dl/B surfaces'
      call ezconc(capr, y, dldb_tot12, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)


      numb = 20

      title = 'q surfaces'
      call ezconc(capr, y, qsafety, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)


      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = 1.0 / bmod(i,j)
         end do
      end do

!      numb = 10

      title = 'Mod 1/B surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)


      title = 'bmod_mid surfaces'
      call ezconc(capr, y, bmod_mid, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)


      title = 'capr bpol_{mid2} surfaces'
      call ezconc(capr, y, capr_bpol_mid2, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)



      do i = 1, nnodex
         fmid1(i) = bx(i, jmid)
      end do

      title= 'Midplane B_x'
      titll= 'B (T)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmid1, &
         nnodex, nxmx)


      do i = 1, nnodex
         fmid1(i) = by(i, jmid)
      end do

      title= 'Midplane B_y'
      titll= 'B (T)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmid1, &
         nnodex, nxmx)


      do i = 1, nnodex
         fmid1(i) = bz(i, jmid)
      end do

      title= 'Midplane B_z'
      titll= 'B (T)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmid1, &
         nnodex, nxmx)


      title= 'Flux average bmod_mid'
      titll= 'bmod_mid (T)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, bmod_midavg, &
          nnoderho_half, nrhomax)


      title= 'Flux average capr*bpol'
      titll= 'capr*bpol (mT)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, capr_bpol_midavg, &
           nnoderho_half, nrhomax)

!      do n = 1, nnoderho
!         write(6,*)n, dldbavg(n)
!      end do


      title= 'Integral dl/B'
      titll= 'Integral dl/B (m)'
      titlr='       '

      call ezplot1_0(title, titll, titlr, rhon_half, dldbavg, &
          nnoderho_half, nrhomax)







      numb = 15



      title = 'B poloidal'
      call ezconc(capr, y, btau, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      title='3-D plot of B poloidal'
      titz='Bpol'
!      call ezcon3d(capr, y, btau, ff, nnodex, nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity, titz)


      title = 'B toroidal'
      call ezconc(capr, y, bzeta, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      title='3-D plot of B toroidal'
      titz='Btor'
!      call ezcon3d(capr, y, bzeta, ff, nnodex, nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity, titz)



      title = 'grad parallel B'
      call ezconc(capr, y, gradprlb, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      title='3-D plot of gradprlb'
      titz='gradprlb'
!      call ezcon3d(capr, y, gradprlb, ff, nnodex, nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity, titz)



!
!--plot q profile in midplane
!

      call a1mnmx(capr, nxmx, nnodex, caprmin, caprmax)
      caprminp = caprmin * 1.01
      caprmaxp = caprmax * .99

      call a2mnmx(qmid, nxmx, nnodex, &
         capr, caprminp, caprmaxp, qmin, qmax)
      qmin = 0.0
!     qmax = 5.0

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (caprmin, caprmax, qmin, qmax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BNST', 0.0, 0)
      call pgmtxt('b', 3.2, 0.5, 0.5, 'R (m)')
      call pgmtxt('t', 2.0, 0.5, 0.5, 'q profile in midplane')

      CALL PGSCI(nred)
      call pgmtxt('l', 2.0, 0.5, 0.5, 'q midplane')
      call pgline(nnodex, capr, qmid)

      CALL PGSCI(nblue)
      call pgsls(2)
      call pgline(nnodex, capr, xiotamid)

      call a1mnmx(bpmid, nxmx, nnodex, bpmin, bpmax)
      bpmin = 0.0

      CALL PGSWIN (caprmin, caprmax, bpmin, bpmax)
      CALL PGSCI(nblack)
      CALL PGBOX  (' ', 0.0, 0, 'CMST', 0.0, 0)

      CALL PGSCI(ngreen)
      call pgsls(3)
      call pgline(nnodex, capr, bpmid)
      call PGMTXT ('r', 2.0, 0.5, 0.5, 'btau')

      CALL PGSCI(nblack)
      call pgsls(1)



      numb = 50

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dbxdx(i,j)
         end do
      end do

      title = 'dbxdx surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dbydx(i,j)
         end do
      end do

      title = 'dbydx surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dbzdx(i,j)
         end do
      end do

      title = 'dbzdx surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dbxdy(i,j)
         end do
      end do

      title = 'dbxdy surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dbydy(i,j)
         end do
      end do

      title = 'dbydy surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dbzdy(i,j)
         end do
      end do

      title = 'dbzdy surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dbdx(i,j)
         end do
      end do

      title = 'dbdx surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dbdy(i,j)
         end do
      end do

      title = 'dbdy surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dxxbxn(i,j)
         end do
      end do

      title = 'dxxbxn surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dxxbyn(i,j)
         end do
      end do

      title = 'dxxbyn surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dxxbzn(i,j)
         end do
      end do

      title = 'dxxbzn surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)




      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dxybxn(i,j)
         end do
      end do

      title = 'dxybxn surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dxybyn(i,j)
         end do
      end do

      title = 'dxybyn surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dxybzn(i,j)
         end do
      end do

      title = 'dxybzn surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)




      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dyybxn(i,j)
         end do
      end do

      title = 'dyybxn surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dyybyn(i,j)
         end do
      end do

      title = 'dyybyn surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dyybzn(i,j)
         end do
      end do

      title = 'dyybzn surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)




      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dxxmodb(i,j)
         end do
      end do

      title = 'dxxmodb surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dxymodb(i,j)
         end do
      end do

      title = 'dxymodb surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = dyymodb(i,j)
         end do
      end do

      title = 'dyymodb surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag,1.)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity,1.)


! Close the graphics device.

      call pgclos

      deallocate ( xkti, xkte, xn, &
         rho, xjy, bmod, bmod_mid, capr_bpol_mid2, &
         xjx, psi, xiota, qsafety, btau, bzeta, &
         dbxdx, dbydx, dbzdx, dbxdy, dbydy, dbzdy, dbdx, dbdy, &
         gradprlb, dxxbxn, dxxbyn, dxxbzn, &
         dxybxn, dxybyn, dxybzn, dyybxn, dyybyn, dyybzn, &
         dxxmodb, dxymodb, dyymodb, bx, by, bz, bz_dum, dldb_tot12, &
         freal, fimag, fmod, rho_tor2d )





  310 format(1p6e12.4)
  309 format(10i10)
  311 format(1p10e12.4)
 1312 format(i10,1p8e12.4)
 1313 format(2i10,1p8e12.4)
 9310 format(1p7e12.4)
   10 format(i10,1p4e10.3,i10,1pe10.3)

      return
      end
!
!********************************************************************
!

!----------------------------------------------------------------------------!
! Subroutine a2dmnmx_eq
!----------------------------------------------------------------------------!
! Minimum and the maximum elements of a 2-d array.
!----------------------------------------------------------------------------!

      subroutine a2dmnmx_eq(f, nxmax, nymax, nx, ny, fmin, fmax)

      implicit none

      integer:: nx, ny, i, j, nxmax, nymax
      real f(nxmax, nymax), fmin, fmax

      fmax=f(2, 1)
      fmin=fmax
      do 23000 j=1, ny
!         do 23002 i=1, nx
         do 23002 i=2, nx-1
            fmax=amax1(fmax, f(i, j))
            fmin=amin1(fmin, f(i, j))
23002    continue
23000 continue

      return
 2201 format(2i5,1p8e12.4)
      end


!
!***************************************************************************
!


      subroutine ezplot1_0(title, titll, titlr, x1, y1, nr, nrmax)

      implicit none

      integer:: nr,nrmax
      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
       nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
       nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
       ncolion,ncolelec,ncollin,ncolln2,ncollab

      ncolln3=ngreen

      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = x1(nr)
      rhomin = x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1

      ymin = 0.0

! Advance plotter to a new page, define coordinate range of graph and draw axes

!      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


! Label the axes (note use of \u and \d for raising exponent).

      call pglab('rho', titll, title)

! Plot the line graph.

      call pgline(nr, x1, y1)

  300 format (1p9e11.3)

      return
      end
!
!***************************************************************************
!
