      subroutine plotws()
      !> Not called
      use size_mod

      implicit none

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb
      character(32):: titx
      character(32):: tity
      character(32):: titz

      integer:: ncolln10, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, ncollin, ncollab, ncolion,
     &   ncolelec, norange

      integer:: nlmx, n, ibackground

      parameter (nlmx = 256)

      integer:: npoints, nnodelb, lnwidth
      integer:: pgopen, pgbeg, ier, np, nl
      real:: xklsavp(nlmx), Elnp(nlmx)
      real:: lb(nlmx), lbprime(nlmx)
      complex El(nlmx)


      open(unit=938,file='out938',status='old',form='formatted')

      read(938, 309) nnodelb
      read(938, 310) (lb(nl), nl = 1, nnodelb)
      read(938, 310) (lbprime(nl), nl = 1, nnodelb)
      read(938, 310) (El(nl), nl = 1, nnodelb)

      read(938, 309) npoints
      read(938, 310) (xklsavp(np), np = 1, npoints)
      read(938, 310) (Elnp(np), np = 1, npoints)

      lnwidth = 2


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

c Open graphics device


#ifdef GIZA      
      IER = PGBEG(0, 'scale.pdf', 1, 1)
#else
      IER = PGBEG(0, 'scale.ps/vcps', 1, 1)
#endif

      IF (IER.NE.1) STOP

      call PGSCF(2)
      call PGSCH (1.5)
      CALL PGSCI(1)
      CALL PGSLW(lnwidth)

c
c--plot El along field line
c
      title= 'E(l) along B'
      titll= 'Mod El (V/m)'
      titlr='       '
      titlb='l(m)'

      call ezplot1q(title, titll, titlr, titlb, lb, real(El),
     &    nnodelb, nlmx)

c
c--plot El spectrum
c
      title= 'Spectrum of E(k) along B'
      titll= 'Mod E(k) (V/m)'
      titlr='       '
      titlb='kl (m-1)'

      call ezplot1q(title, titll, titlr, titlb, xklsavp, Elnp,
     &    npoints, nlmx)

      call pgclos


  310 format(1p,6e12.4)
  309 format(10i10)

      close (938)

      return
      end
c
c********************************************************************
c

      subroutine fieldws(dfquotient, rmin_zoom, rmax_zoom)

      use size_mod
      use aorsa2din_mod, only: eqtype, nzeta_wdoti
      
      implicit none

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb
      character(32):: titx
      character(32):: tity
      character(32):: titz

      real:: dfquotient

      real:: logmax, ycut, dy, xmax, ymax, xmi, E_eV, vperp_mks,
     &   vpara_mks, vperp_cgs, uperp_1kev, duperp, dz, dx
      real:: exkmin, exkmax, prfin, rmin_zoom, rmax_zoom

      integer:: pgopen, pgbeg, ier, nmid, mmid, it, mmax, nmax
      integer:: n_theta_max, n_u_max, n_psi_max, idiag, jdiag

      integer:: ndisti1, ndisti2, ndisti3, number_points, k, nnodez
      integer:: nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, nlevmax
      integer:: ibackground, nkx1, nkx2, nky1, nky2, n, m,
     &   nkpltdim, mkpltdim, nkxplt, nkyplt, ipage, n1, n2, n3, n4
      integer:: i_psi, i_psi1, i_psi2, i_psi3, i_psi4, i_psi5, i_psi6
      integer:: i_psi_array(6), i_psi_index
      integer:: nkpr, ierr

      integer:: nuper, nupar, i_uperp, i_upara, nz
      real:: uminpara, umaxpara, vc_cgs, vpara_cgs

      integer::  ncolln10, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, ncollin, ncollab, ncolion,
     &   ncolelec, norange

      integer:: nrhomax, nnoderho, nnoderho2, iflag, nnoderho_half,
     &   nnoderho2_half, mnodetheta, nthetamax, ntable, mtable


      parameter (nlevmax = 101)
      parameter (ntable = 1000)
      parameter (mtable = 1000)

      parameter (nxmx = nmodesmax)
      parameter (nymx = mmodesmax)

      parameter (nrhomax = nxmx)
      parameter (nthetamax = nymx)

      parameter (nkdim1 = - nxmx / 2)
      parameter (nkdim2 =   nxmx / 2)

      parameter (mkdim1 = - nymx / 2)
      parameter (mkdim2 =   nymx / 2)

      parameter (nkpltdim = 2 * nkdim2)
      parameter (mkpltdim = 2 * mkdim2)

      parameter (n_theta_max = 150)
      parameter (n_u_max = 150)
      parameter (n_psi_max = 150)

      real:: u(n_u_max), theta_u(n_theta_max)
      real:: f_cql(n_theta_max, n_u_max, n_psi_max)
      real:: f_cql_2d(n_u_max, n_theta_max)
      real:: f_cql_1d_1(n_u_max), f_cql_1d_2(n_u_max)
      real:: f_cql_1d_3(n_u_max), f_cql_1d_4(n_u_max)
      real:: f_cql_1d_5(n_u_max), f_cql_1d_6(n_u_max)
      real:: f_cql_1d_7(n_u_max)


      integer:: n_theta_(n_psi_max)
      real:: theta_(n_theta_max, n_psi_max)

      integer:: n_theta, n_u, n_psi
      integer:: i_theta, i_u
      integer:: i_theta1, i_theta2, i_theta3, i_theta4,
     &     i_theta5, i_theta6, i_theta7
      integer, parameter :: r15 = selected_real_kind(15)
      real:: vc, r0

      real, dimension(:),   allocatable :: UPERP, UPARA

      real, dimension(:,:), allocatable :: bqlavg_i1_2d
      real, dimension(:,:), allocatable :: cqlavg_i1_2d

      real, dimension(:,:), allocatable :: E_kick_2d
      real, dimension(:,:), allocatable :: f_cql_cart_2d
      real, dimension(:,:,:), allocatable :: f_cql_cart

      real, dimension(:), allocatable :: wperp1_cql
      real, dimension(:), allocatable :: wpar1_cql

      real, dimension(:), allocatable :: wperp2_cql
      real, dimension(:), allocatable :: wpar2_cql

      real, dimension(:,:,:), allocatable :: bqlavg_i1
      real, dimension(:,:,:), allocatable :: cqlavg_i1
      real, dimension(:,:,:), allocatable :: eqlavg_i1
      real, dimension(:,:,:), allocatable :: fqlavg_i1


      real:: capr_bpol_mid(nrhomax), capr_bpol_midavg(nrhomax)
      real:: capr_bpol_mid2(nxmx, nymx), bmod_midavg(nrhomax)
      real:: bmod_mid(nxmx, nymx), bratio(nxmx, nymx)

      real:: xkxsav(nkpltdim), xkysav(mkpltdim)
      real:: wdoti1avg(nrhomax), wdoti2avg(nrhomax), wdoti3avg(nrhomax)
      real:: wdoti4avg(nrhomax), wdoti5avg(nrhomax), wdoti6avg(nrhomax)
      real:: zdummy(3)

      real:: xjprl_int(nrhomax)

      real:: redotje_int(nrhomax),
     &     redotj1_int(nrhomax), redotj2_int(nrhomax),
     &     redotj3_int(nrhomax), redotj4_int(nrhomax),
     &     redotj5_int(nrhomax), redotj6_int(nrhomax)

      real:: wdote_int(nrhomax),
     &     wdot1_int(nrhomax), wdot2_int(nrhomax),
     &     wdot3_int(nrhomax), wdot4_int(nrhomax),
     &     wdot5_int(nrhomax), wdot6_int(nrhomax)

      real:: wdoti1_dvol(nrhomax), wdoti2_dvol(nrhomax),
     &     wdoti3_dvol(nrhomax), wdoti4_dvol(nrhomax),
     &     wdoti5_dvol(nrhomax), wdoti6_dvol(nrhomax),
     &     wdote_dvol(nrhomax)

      real:: redotj1_dvol(nrhomax), redotj2_dvol(nrhomax),
     &     redotj3_dvol(nrhomax), redotj4_dvol(nrhomax),
     &     redotj5_dvol(nrhomax), redotj6_dvol(nrhomax),
     &     redotje_dvol(nrhomax)


      real:: wdoteavg_int(nrhomax), wdoti1avg_int(nrhomax),
     &                            wdoti2avg_int(nrhomax),
     &                            wdoti3avg_int(nrhomax),
     &                            wdoti4avg_int(nrhomax),
     &                            wdoti5avg_int(nrhomax),
     &                            wdoti6avg_int(nrhomax)

      real:: wdote_ql_int(nrhomax), wdoti1_ql_int(nrhomax),
     &                            wdoti2_ql_int(nrhomax),
     &                            wdoti3_ql_int(nrhomax),
     &                            wdoti4_ql_int(nrhomax),
     &                            wdoti5_ql_int(nrhomax),
     &                            wdoti6_ql_int(nrhomax)


      real:: redotjeavg_int(nrhomax), redotj1avg_int(nrhomax),
     &                              redotj2avg_int(nrhomax),
     &                              redotj3avg_int(nrhomax),
     &                              redotj4avg_int(nrhomax),
     &                              redotj5avg_int(nrhomax),
     &                              redotj6avg_int(nrhomax)

      real:: rhon(nrhomax), thetam(nthetamax)

      real:: xnavg(nrhomax), xn1avg(nrhomax), xkteavg(nrhomax),
     &     xktiavg(nrhomax),
     &     xkti1avg(nrhomax), xkti2avg(nrhomax), xkti3avg(nrhomax),
     &     xn2avg(nrhomax), xn3avg(nrhomax), xk3avg(nrhomax),
     &     xna_sloavg(nrhomax),
     &     vyi1avg(nrhomax),  vyi2avg(nrhomax),
     &     dvol(nrhomax), rhon_half(nrhomax), volume(nrhomax)
      real:: rhon_save(nrhomax), rhon_half_save(nrhomax)
      real:: redotj2avg_save(nrhomax), wdoti2avg_save(nrhomax)

      real:: xkti4avg(nrhomax),  xkti5avg(nrhomax), xkti6avg(nrhomax)
      real:: wdote_ql(nrhomax), wdoti1_ql(nrhomax), wdoti2_ql(nrhomax),
     &     wdoti3_ql(nrhomax), wdoti4_ql(nrhomax), wdoti5_ql(nrhomax),
     &     wdoti6_ql(nrhomax), dldbavg(nrhomax), gradprlb2_avg(nrhomax)
      real:: xn4avg(nrhomax), xn5avg(nrhomax), xn6avg(nrhomax)
      real:: vyavg(nrhomax), dvydrho(nrhomax), wdoteavg(nrhomax)
      real:: fz0i1avg(nrhomax), fz0i2avg(nrhomax), fz0i3avg(nrhomax),
     &     fz0eavg(nrhomax), fz0avg(nrhomax)

      real:: fz0i4avg(nrhomax), fz0i5avg(nrhomax), fz0i6avg(nrhomax)

      real:: gpsi_avg(nrhomax), kpsi_avg(nrhomax), xjhat(nrhomax),
     &     muhat_avg(nrhomax), nu_star_avg(nrhomax), xkhat(nrhomax)

      real:: redotjeavg(nrhomax), redotj1avg(nrhomax),
     &     redotj2avg(nrhomax), redotj3avg(nrhomax),
     &     redotj4avg(nrhomax),
     &     redotj5avg(nrhomax), redotj6avg(nrhomax),
     &     redotjsavg(nrhomax),
     &     redotjtavg(nrhomax), xjprlavg(nrhomax),
     &     redotjiavg(nrhomax), ipsi_avg(nrhomax)

      real:: fz0_int(nrhomax)

      real:: xnmid(nxmx), xktimid(nxmx),  xktemid(nxmx), qmid(nxmx)
      real:: bpmid(nxmx), xiotamid(nxmx), xjymid(nxmx)
      real:: xjxmid(nxmx), xjzmid(nxmx)

      real:: xjxmin, xjxmax
      real:: xjymin, xjymax
      real:: xjzmin, xjzmax

      real:: fmidre(nxmx), fmidim(nxmx), fmid1(nxmx),
     &     fmid2(nxmx), fmid3(nxmx), fmid4(nxmx), fmid5(nxmx),
     &     fmids(nxmx), fmidt(nxmx), fmid6(nxmx)
      real:: x(nxmx), capr(nxmx), y(nymx), xkphi(nxmx)
      real:: mod_Eplus_mid(nxmx),  mod_Eminus_mid(nxmx),
     &   mod_Eb_mid(nxmx), mod_e_mid(nxmx)


      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     &   ncolion,ncolelec,ncollin,ncollab,ncolln2

      common/boundcom/rhoplasm
      common/zoom/ xmaxz, xminz, ymaxz, yminz

      real:: xmaxz, xminz, ymaxz, yminz

      real, dimension(:,:), allocatable :: exkmod,
     &     eykmod, ezkmod, exklog, eyklog, ezklog

      real:: fmodm(mkpltdim)
      complex:: z0_table1(ntable, mtable)
      complex:: z1_table1(ntable, mtable)
      complex:: z2_table1(ntable, mtable)
      real:: zetai_table(ntable)
      real:: dKdL_table(mtable)

      complex, dimension(:,:), allocatable :: exk, eyk, ezk

      real:: fmin, fmax, fminre, fmaxre, fminim, fmaxim, fmin1, isq2,
     &   fmax1, fmax2, fmax3, fmaxs, fmaxt, fmin2, fmin3, fmins,fmint

      real:: exklogmin, exklogmax,
     &     eyklogmin, eyklogmax,
     &     ezklogmin, ezklogmax

      real, dimension(:,:), allocatable ::  xn, xkti, xkte, rho, theta,
     &   psi, xjy, bmod, xjx, xiota, qsafety, ipsi, btau, bzeta,
     &   freal, fimag, fmod, mod_Eplus, mod_Eminus, mod_Eb,
     &   mod_Ealpha, mod_Ebeta, mod_E,
     &   dldb_tot12, reex_dx, reey_dx, reez_dx, ximex_dx,
     &   ximey_dx, ximez_dx, spx, spy, spz, reomg1a, reomg2a,
     &   reomg3a, capr_bpol, pressi, redotj1, redotj2, redotj3,
     &   redotje, redotjt, redotjs, wdoti1, wdoti2, wdoti3,
     &   wdoti4, wdoti5, wdoti6, wdote, wdott, fz0e, fz0i1,
     &   fz0i2, fz0i3, fz0, fz0i4, fz0i5, fz0i6, gpsi,  kpsi,
     &   omgexb, uzeta, utheta, fpsi0, ftheta0, muhat, nu_star,
     &   redotj4, redotj5, redotj6, divq, capd, capd_plot, xkb,
     &   xkb_plot, reomglha, xjz


      real, dimension(:,:), allocatable :: dxxuyy, dyyuzz

      real, dimension(:,:), allocatable :: gradprlb

      complex, dimension(:,:), allocatable :: xkperp_cold, acold,
     &   bcold, ccold,
     &   ex, ey, ez, bxwave, bywave, bzwave, xkperp_cold2,
     &   ealpha, ebeta, eb, eplus, eminus, eplus_flux_plot,
     &   eminus_flux_plot, xkperp_flux_plot, ntilda_e

      complex, dimension(:,:), allocatable :: xkperp2_slow,
     &        xkperp2_fast, P_a
      real, dimension(:,:), allocatable :: xkprl_a

      complex, dimension(:,:), allocatable :: xjpxe, xjpye, xjpze


      complex, dimension(:,:), allocatable :: eplus_flux,
     &   eminus_flux,  xkperp_flux

      real:: bmod_flux(nrhomax, nthetamax)
      real:: capr_flux(nrhomax,nthetamax), capz_flux(nrhomax, nthetamax)


      real:: bmod_plot(nrhomax)


      real:: ff(101)
      real:: q, omgrf, xk0, n0, clight, xmu0, eps0, rhoplasm
      real:: period, pi, time, dt
      real:: temax, temin, timin, tmin, tmax, timax, caprmaxp,
     &   caprmin, caprminp, caprmax, xnmax, xnmin, qmin, qmax
      real:: bpmin, bpmax, scalex
      integer:: nnodex, j, i, nnodey, numb, jmid, jcut, jmid_sav
      complex:: zi


!--------------End declarations-------------------
      scalex = 1.0
      if (eqtype=='mirror') scalex=10.0 !JCW mirror maybe make namelist input


      
      allocate (xkperp_cold(nxmx, nymx), acold(nxmx, nymx),
     &   bcold(nxmx, nymx), ccold(nxmx, nymx),
     &   ex(nxmx, nymx), ey(nxmx, nymx), ez(nxmx, nymx),
     &   bxwave(nxmx, nymx), bywave(nxmx, nymx), bzwave(nxmx, nymx),
     &   ealpha(nxmx, nymx), ebeta(nxmx, nymx), eb(nxmx, nymx),
     &   eplus(nxmx, nymx), eminus(nxmx, nymx),
     &   eplus_flux_plot(nxmx, nymx), xkperp_cold2(nxmx, nymx),
     &   eminus_flux_plot(nxmx, nymx),
     &   xkperp_flux_plot(nxmx, nymx), ntilda_e(nxmx, nymx)  )

      allocate (xjpxe(nxmx, nymx), xjpye(nxmx, nymx), xjpze(nxmx, nymx))

      allocate (xkperp2_slow(nxmx, nymx), xkperp2_fast(nxmx, nymx) )
      allocate (xkprl_a(nxmx, nymx) )
      allocate (P_a(nxmx, nymx) )

      allocate (dxxuyy(nxmx, nymx), dyyuzz(nxmx, nymx))

      allocate (gradprlb(nxmx, nymx))

      allocate ( exkmod(nkpltdim, mkpltdim),
     &   eykmod(nkpltdim, mkpltdim),
     &   ezkmod(nkpltdim, mkpltdim),
     &   exklog(nkpltdim, mkpltdim),
     &   eyklog(nkpltdim, mkpltdim),
     &   ezklog(nkpltdim, mkpltdim)  )

      allocate ( eplus_flux(nrhomax, nthetamax),
     &   eminus_flux(nrhomax, nthetamax),
     &   xkperp_flux(nrhomax, nthetamax)  )

      allocate ( exk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &        eyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &        ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)  )


      allocate ( xn(nxmx, nymx),
     &   xkti(nxmx, nymx), xkte(nxmx, nymx),
     &   rho(nxmx, nymx), theta(nxmx, nymx), psi(nxmx, nymx),
     &   xjy(nxmx, nymx), bmod(nxmx, nymx),
     &   xjx(nxmx, nymx), xjz(nxmx, nymx),
     &   xiota(nxmx, nymx), qsafety(nxmx, nymx), ipsi(nxmx, nymx),
     &   btau(nxmx, nymx), bzeta(nxmx, nymx),
     &   freal(nxmx, nymx), fimag(nxmx, nymx), fmod(nxmx, nymx),
     &   mod_Eplus(nxmx, nymx), mod_Eminus(nxmx, nymx),
     &   mod_Eb(nxmx, nymx), mod_Ealpha(nxmx, nymx),
     &   mod_Ebeta(nxmx, nymx), mod_E(nxmx, nymx),
     &   dldb_tot12(nxmx, nymx), reex_dx(nxmx, nymx),
     &     reey_dx(nxmx, nymx), reez_dx(nxmx, nymx) ,stat=ierr)

      allocate(
     &   ximex_dx(nxmx, nymx), ximey_dx(nxmx, nymx),
     &   ximez_dx(nxmx, nymx),
     &   spx(nxmx, nymx), spy(nxmx, nymx), spz(nxmx, nymx),
     &   reomg1a(nxmx, nymx), reomg2a(nxmx, nymx), reomg3a(nxmx, nymx),
     &   capr_bpol(nxmx, nymx), pressi(nxmx, nymx),
     &   redotj1(nxmx, nymx), redotj2(nxmx, nymx), redotj3(nxmx, nymx),
     &   redotje(nxmx, nymx), redotjt(nxmx, nymx), redotjs(nxmx, nymx),
     &   wdoti1(nxmx, nymx), wdoti2(nxmx, nymx), wdoti3(nxmx, nymx),
     &   wdoti4(nxmx, nymx), wdoti5(nxmx, nymx), wdoti6(nxmx, nymx),
     &   wdote(nxmx, nymx), wdott(nxmx, nymx),
     &   fz0e(nxmx, nymx), fz0i1(nxmx, nymx),
     &     fz0i2(nxmx, nymx), fz0i3(nxmx, nymx), fz0(nxmx, nymx))
      allocate(
     &   fz0i4(nxmx, nymx), fz0i5(nxmx, nymx), fz0i6(nxmx, nymx),
     &   gpsi(nxmx, nymx),  kpsi(nxmx, nymx),
     &   omgexb(nxmx, nymx), uzeta(nxmx, nymx), utheta(nxmx, nymx),
     &   fpsi0(nxmx, nymx), ftheta0(nxmx, nymx),
     &   muhat(nxmx, nymx), nu_star(nxmx, nymx),
     &   redotj4(nxmx, nymx), redotj5(nxmx, nymx), redotj6(nxmx, nymx),
     &   divq(nxmx, nymx), reomglha(nxmx, nymx) )

      allocate(capd(nxmx, nkpltdim))
      allocate(capd_plot(nxmx, nkpltdim))

      allocate(xkb(nxmx, nkpltdim))
      allocate(xkb_plot(nxmx, nkpltdim))

!JCW init to zero
c--set default values of input data:

      ibackground = 0

      xminz = .58
      xmaxz = .72
      ymaxz = .30
      yminz = -99.99

      logmax = 2.0
      ipage = 2
      numb = 20
      ycut = -0.0

      i_psi1 = 2
      i_psi2 = 10
      i_psi3 = 20
      i_psi4 = 30
      i_psi5 = 40
      i_psi6 = 50

      i_psi_array(1) = 2
      i_psi_array(2) = 10
      i_psi_array(3) = 20
      i_psi_array(4) = 30
      i_psi_array(5) = 40
      i_psi_array(6) = 50


      if(yminz .eq. -99.99)yminz = -ymaxz


      open(unit=38,file='out38',status='old',form='formatted')

      open(unit=140,file='Efield_2D.vtk',status='unknown',
     &                                       form='formatted')
      open(unit=245,file='Poynting_2D.vtk',status='unknown',
     &                                       form='formatted')
      open(unit=141,file='E_kicks_2D.vtk',status='unknown',
     &                                       form='formatted')
      open(unit=142,file='Bql_avg_2D.vtk',status='unknown',
     &                                       form='formatted')
      open(unit=242,file='Cql_avg_2D.vtk',status='unknown',
     &                                       form='formatted')

      open(unit=143,file='capd.vtk',status='unknown',
     &                                       form='formatted')

      open(unit=144,file='Eb_spectrum.vtk',status='unknown',
     &                                       form='formatted')

      open(unit=51,file='rho',status='unknown',form='formatted')

      open(unit=66,file='bharvey',status='unknown',form='formatted')
      open(unit=57,file='swain',status='unknown',form='formatted')
      open(unit=62,file='murakami',status='unknown', form='formatted')
      open(unit=72,file='bertelli',status='unknown', form='formatted')
      open(unit=65,file='mchoi',status='unknown', form='formatted')
      open(unit=67,file='mchoi2',status='unknown', form='formatted')
      open(unit=69,file='mchoi3',status='unknown', form='formatted')
      open(unit=68,file='zowens',status='unknown', form='formatted')

      open(unit=64,file='E_lab_frame',status='unknown',form='formatted')


      zi = cmplx(0.0, 1.0)
      eps0 = 8.85e-12
      xmu0 = 1.26e-06
      clight = 1./sqrt(eps0 * xmu0)
      xk0 = omgrf / clight
      q = 1.6e-19

      fmid1 = 0.0
      fmid2 = 0.0
      fmid3 = 0.0
      fmid4 = 0.0
      fmid5 = 0.0
      fmid6 = 0.0

      read(38, 309) nmax, mmax
      read(38, 310) (zetai_table(n), n = 1, nmax)
      read(38, 310) (dKdL_table(m), m = 1, mmax)
      read(38, 310) ((z0_table1(n, m), n = 1, nmax),  m = 1, mmax)
      read(38, 310) ((z1_table1(n, m), n = 1, nmax),  m = 1, mmax)
      read(38, 310) ((z2_table1(n, m), n = 1, nmax),  m = 1, mmax)

      read(38, 309) nnodex, nnodey, jmid, idiag, jdiag
      read(38, 310) rhoplasm, prfin, omgrf
      read(38, 310) (x(i), i = 1, nnodex)
      read(38, 310) (y(j), j = 1, nnodey)
      read(38, 310) (capr(i), i = 1, nnodex)
      read(38, 310) (xkphi(i), i = 1, nnodex)

      xmax = x(nnodex) - x(1)
      ymax = y(nnodex) - y(1)


      read(38, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((theta(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((xn(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((xkte(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((xkti(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((xkperp_cold(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((xkperp_cold2(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((acold(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((bcold(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((ccold(i, j), i = 1, nnodex), j = 1, nnodey)

      read (38, 310) ((reomg1a(i, j), i = 1, nnodex), j = 1, nnodey)
      read (38, 310) ((reomg2a(i, j), i = 1, nnodex), j = 1, nnodey)
      read (38, 310) ((reomg3a(i, j), i = 1, nnodex), j = 1, nnodey)

      read (38, 310) ((reomglha(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((xjx(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((xjy(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((xjz(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)

      read (38, 310) ((xjpxe(i,j), i = 1, nnodex), j = 1, nnodey)
      read (38, 310) ((xjpye(i,j), i = 1, nnodex), j = 1, nnodey)
      read (38, 310) ((xjpze(i,j), i = 1, nnodex), j = 1, nnodey)



      write(51,309) nnodex, nnodey
      write(51,310) rhoplasm
      write(51,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)



      read(38, 309) nkx1, nkx2
      read(38, 309) nky1, nky2

      nkxplt = 2 * nkx2
      nkyplt = 2 * nky2

      read (38, 310) (xkxsav(n), n = 1, nkxplt)
      read (38, 310) (xkysav(m), m = 1, nkyplt)

      read(38, 310) ((exkmod(n, m), n = 1, nkxplt), m = 1, nkyplt)
      read(38, 310) ((eykmod(n, m), n = 1, nkxplt), m = 1, nkyplt)
      read(38, 310) ((ezkmod(n, m), n = 1, nkxplt), m = 1, nkyplt)



c
c--calculate log exkmod, eykmod, ezkmod
c
      do n = 1, nkxplt
         do m = 1, nkyplt


            exklog(n, m) = 0.0
            eyklog(n, m) = 0.0
            ezklog(n, m) = 0.0

            if(abs(exkmod(n,m)) .ne. 0.0)
     &         exklog(n,m)=alog10(abs(exkmod(n,m)))
            if(abs(eykmod(n,m)) .ne. 0.0)
     &         eyklog(n,m)=alog10(abs(eykmod(n,m)))
            if(abs(ezkmod(n,m)) .ne. 0.0)
     &         ezklog(n,m)=alog10(abs(ezkmod(n,m)))
         end do
      end do




      call a2dmnmx_r4(exklog, nkpltdim, mkpltdim, nkxplt, nkyplt,
     &    exklogmin, exklogmax)

      call a2dmnmx_r4(eyklog, nkpltdim, mkpltdim, nkxplt, nkyplt,
     &    eyklogmin, eyklogmax)

      call a2dmnmx_r4(ezklog, nkpltdim, mkpltdim, nkxplt, nkyplt,
     &    ezklogmin, ezklogmax)


c--   Normalize logs to a maximum value of logmax
      do n = 1, nkxplt
         do m = 1, nkyplt
            if(exklog(n,m) .ne. 0.0)
     &         exklog(n,m) = exklog(n,m) + (logmax - exklogmax)
            if(eyklog(n,m) .ne. 0.0)
     &         eyklog(n,m) = eyklog(n,m) + (logmax - eyklogmax)
            if(ezklog(n,m) .ne. 0.0)
     &         ezklog(n,m) = ezklog(n,m) + (logmax - ezklogmax)
         end do
      end do


      read(38, 310) ((redotje(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj1(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj2(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj3(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotjt(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((bmod(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((qsafety(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((ebeta(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((btau(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((bzeta(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((spx(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((spy(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((spz(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((exk(n, m), n = nkx1, nkx2), m = nky1, nky2)
      read(38, 310) ((eyk(n, m), n = nkx1, nkx2), m = nky1, nky2)
      read(38, 310) ((ezk(n, m), n = nkx1, nkx2), m = nky1, nky2)

c      write(6, *)"ealphak(32,32) = ", exk(32,32)
c      write(6, *)"ebetak(32,32)  = ", eyk(32,32)
c      write(6, *)"ebk(32,32)     = ", ezk(32,32)

      read(38, 310) ((wdote(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti1(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti2(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdott(i, j), i = 1, nnodex), j = 1, nnodey)

      read (38, 310) ((divq(i, j), i = 1, nnodex), j = 1, nnodey)
      read (38, 310) ((capd(i, n), i = 1, nnodex), n = 1, nkxplt)
      read (38, 310) ((xkb(i, n), i = 1, nnodex), n = 1, nkxplt)


      do n = 1, nkxplt
         write(6, *) "n = ",n,"capd(nnodex/2, n) = ",capd(nnodex/2, n)
      end do

      read(38, 309) nnoderho
      read(38, 310) (rhon(n), n = 1, nnoderho)
      read(38, 310) (gradprlb2_avg(n), n = 1, nnoderho)

      read(38, 309) mnodetheta
      read(38, 310) (thetam(m), m = 1, mnodetheta)

      read(38, 310) (dvol(n), n = 1, nnoderho)
      read(38, 310) (volume(n), n = 1, nnoderho)

      read(38, 310) (xnavg(n), n = 1, nnoderho)
      read(38, 310) (wdoteavg(n), n = 1, nnoderho)
      read(38, 310) (wdoti1avg(n), n = 1, nnoderho)
      read(38, 310) (wdoti2avg(n), n = 1, nnoderho)
      read(38, 310) (vyavg(n), n = 1, nnoderho)
      read(38, 310) (dvydrho(n), n = 1, nnoderho)

      read(38, 310) (redotjeavg(n), n = 1, nnoderho)
      read(38, 310) (redotj1avg(n), n = 1, nnoderho)
      read(38, 310) (redotj2avg(n), n = 1, nnoderho)
      read(38, 310) (redotj3avg(n), n = 1, nnoderho)
      read(38, 310) (redotjtavg(n), n = 1, nnoderho)

      read(38, 310) ((redotjs(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti3(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) (wdoti3avg(n), n = 1, nnoderho)

      read(38, 310) ((fz0e(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i1(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i2(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i3(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (fz0eavg(n), n = 1, nnoderho)
      read(38, 310) (fz0i1avg(n), n = 1, nnoderho)
      read(38, 310) (fz0i2avg(n), n = 1, nnoderho)
      read(38, 310) (fz0i3avg(n), n = 1, nnoderho)
      read(38, 310) (fz0avg(n), n = 1, nnoderho)

      read(38, 310) (xjprlavg(n), n = 1, nnoderho)

      read(38, 310) (xjprl_int(n), n = 1, nnoderho)
      read(38, 310) (fz0_int(n), n = 1, nnoderho)

      read(38, 310) (xjhat(n), n = 1, nnoderho)
      read(38, 310) (xkhat(n), n = 1, nnoderho)

      read(38, 310) (gpsi_avg(n), n = 1, nnoderho)
      read(38, 310) ((gpsi(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (kpsi_avg(n), n = 1, nnoderho)
      read(38, 310) ((kpsi(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (muhat_avg(n), n = 1, nnoderho)
      read(38, 310) ((muhat(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (nu_star_avg(n), n = 1, nnoderho)
      read(38, 310) ((nu_star(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (ipsi_avg(n), n = 1, nnoderho)
      read(38, 310) ((ipsi(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((omgexb(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((uzeta(i,j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((utheta(i,j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (xn1avg(n), n = 1, nnoderho)
      read(38, 310) (xn2avg(n), n = 1, nnoderho)
      read(38, 310) (xn3avg(n), n = 1, nnoderho)
      read(38, 310) (xna_sloavg(n), n = 1, nnoderho)
      read(38, 310) (xkteavg(n), n = 1, nnoderho)
      read(38, 310) (xktiavg(n), n = 1, nnoderho)
      read(38, 310) (xkti2avg(n), n = 1, nnoderho)
      read(38, 310) (xkti3avg(n), n = 1, nnoderho)

      read(38, 310) (redotjsavg(n), n = 1, nnoderho)

      read(38, 310) ((bmod_mid(i,j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) (bmod_midavg(n), n = 1, nnoderho)

      read(38, 310) ((capr_bpol_mid2(i, j), i = 1,nnodex),
     &                                                    j = 1,nnodey)
      read(38, 310) (capr_bpol_midavg(n), n = 1, nnoderho)

      read(38, 310) ((redotj4(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj5(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj6(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (redotj4avg(n), n = 1, nnoderho)
      read(38, 310) (redotj5avg(n), n = 1, nnoderho)
      read(38, 310) (redotj6avg(n), n = 1, nnoderho)

      read(38, 310) ((wdoti4(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti5(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti6(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (wdoti4avg(n), n = 1, nnoderho)
      read(38, 310) (wdoti5avg(n), n = 1, nnoderho)
      read(38, 310) (wdoti6avg(n), n = 1, nnoderho)

      read(38, 310) ((fz0i4(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i5(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i6(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (fz0i4avg(n), n = 1, nnoderho)
      read(38, 310) (fz0i5avg(n), n = 1, nnoderho)
      read(38, 310) (fz0i6avg(n), n = 1, nnoderho)

      read(38, 310) (xn4avg(n), n = 1, nnoderho)
      read(38, 310) (xn5avg(n), n = 1, nnoderho)
      read(38, 310) (xn6avg(n), n = 1, nnoderho)

      read(38, 310) (xkti4avg(n), n = 1, nnoderho)
      read(38, 310) (xkti5avg(n), n = 1, nnoderho)
      read(38, 310) (xkti6avg(n), n = 1, nnoderho)

      read(38, 310) (wdoteavg_int(n),  n = 1, nnoderho)
      read(38, 310) (wdoti1avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti2avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti3avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti4avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti5avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti6avg_int(n), n = 1, nnoderho)

      read(38, 310) (wdote_ql_int(n),  n = 1, nnoderho)
      read(38, 310) (wdoti1_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti2_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti3_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti4_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti5_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti6_ql_int(n), n = 1, nnoderho)

      read(38, 310) (redotjeavg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj1avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj2avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj3avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj4avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj5avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj6avg_int(n), n = 1, nnoderho)




*     -----------------------
*     Calculate bratio(i,j)
*     -----------------------
      do i = 1, nnodex
         do j = 1, nnodey
            bratio(i,j) = bmod_mid(i,j) / bmod(i,j)
            if (bratio(i,j) .gt. 1.0) bratio(i,j) = 1.0
         end do
      end do



      nnoderho_half = nnoderho - 1

      do n = 1, nnoderho_half
         rhon_half(n) = (rhon(n) + rhon(n+1)) / 2.
      end do

c      do n = 1, nnoderho
c        rhon_save(n) = rhon(n)
c        rhon_half_save(n) = rhon_half(n)
c      end do

      rhon_save = rhon
      rhon_half_save = rhon_half
      redotj2avg_save = redotj2avg
      wdoti2avg_save = wdoti2avg


*     ---------------------
*     write swain (57) file
*     ---------------------
      write(57, 309) nnoderho
      do n = 1, nnoderho
         write(57,1312)n, rhon(n), redotjeavg(n), redotj1avg(n),
     &              redotj2avg(n), redotj3avg(n), redotj4avg(n),
     &              redotj5avg(n), xjprlavg(n)
      end do


*     ---------------------
*     write murakami (62) file
*     ---------------------
      write(62, 309) nnoderho
      do n = 1, nnoderho
         write(62,1312)n, rhon(n), wdoteavg(n), wdoti1avg(n),
     &              wdoti2avg(n), wdoti3avg(n), wdoti4avg(n),
     &              wdoti5avg(n), wdoti6avg(n), xjprlavg(n)
      end do

*     ---------------------
*     write mchoi2 (67) file
*     ---------------------
c      write(67, 309) nnoderho_half
c      do n = 1, nnoderho_half
c         write (67, 1312) n, rhon_half(n), redotj2avg(n),
c     &      wdoti2avg(n)
c      end do

c      jmid_sav = jmid
c      jmid = 135
      do i = 1, nnodex
         xnmid(i) = xn(i, jmid)
         xktemid(i) = xkte(i, jmid) / q
         xktimid(i) = xkti(i, jmid) / q
         qmid(i) = qsafety(i, jmid)
         xiotamid(i) = 1.0 / qmid(i)
         bpmid(i) = btau(i, jmid)
         xjxmid(i) = xjx(i, jmid)
         xjymid(i) = xjy(i, jmid)
         xjzmid(i) = xjz(i, jmid)

      end do
c      jmid = jmid_sav



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
c     ncolbox=nblack
c     ncolbrd=nblack

c     ncolion=nblue
c     ncolelec=nred
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

c Open graphics device

#ifdef GIZA
      IER = PGBEG(0, 'aorsa2d.pdf', 1, 1)
#else
      IER = PGBEG(0, 'aorsa2d.ps/vcps', 1, 1)
#endif

      IF (IER.NE.1) STOP


c
c--plot plasma profiles
c
      call PGSCH (1.5)

      CALL PGSCI(1)
      CALL PGSLW(3)

*     -------------------------------
*     plot 2D contours of Z2_table table
*     -------------------------------
      titx = '1 / zeta'
      tity = 'dK/dL'

      title = 'Real z0 table1'
      call ezconc(zetai_table, dKdL_table, real(z0_table1),
     &     ff, nmax, mmax, numb, ntable, mtable, nlevmax,
     &     title, titx, tity, iflag,0.25)

      title = 'Imag z0 table1'
      call ezconc(zetai_table, dKdL_table, aimag(z0_table1),
     &    ff, nmax, mmax,
     &   numb, ntable, mtable, nlevmax, title, titx, tity, iflag,0.25)

      title = 'Real z1 table1'
      call ezconc(zetai_table, dKdL_table, real(z1_table1),
     &   ff, nmax, mmax,
     &   numb, ntable, mtable, nlevmax, title, titx, tity, iflag,0.25)

      title = 'Imag z1 table1'
      call ezconc(zetai_table, dKdL_table, aimag(z1_table1),
     &    ff, nmax, mmax,
     &   numb, ntable, mtable, nlevmax, title, titx, tity, iflag,0.25)

      title = 'Real z2 table1'
      call ezconc(zetai_table, dKdL_table, real(z2_table1),
     &   ff, nmax, mmax,
     &   numb, ntable, mtable, nlevmax, title, titx, tity, iflag,0.25)

      title = 'Imag z2 table1'
      call ezconc(zetai_table, dKdL_table, aimag(z2_table1),
     &    ff, nmax, mmax,
     &   numb, ntable, mtable, nlevmax, title, titx, tity, iflag,0.25)

      titx = 'R (m)'
      tity = 'Z (m)'

      title = '|B|_{mid} surfaces'
      call ezconc(capr, y, bmod_mid, ff, nnodex, nnodey, 50,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      title = 'Bratio surfaces'
      call ezconc(capr, y, bratio, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = 'I(psi) m-T'
      call ezconc(capr, y, ipsi, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title= 'I(psi) = R B_{phi}'
      titll= 'I(psi) m-T'
      titlr='       '
      call ezplot0(title, titll, titlr, rhon, ipsi_avg,
     &    nnoderho, nrhomax)


      title = 'R*B_{pol}-mid2 surfaces'
      call ezconc(capr, y, capr_bpol_mid2, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      title= '<|B|_{mid}>'
      titll= 'bmod_{mid} (T)'
      titlr='       '

      call ezplot0(title, titll, titlr, rhon_half, bmod_midavg,
     &    nnoderho_half, nrhomax)

      title= 'parallel gradient of B'
      titll= 'gradprlb avg (T/m)'
      titlr='       '
!JCW note this avg~0 for mirrors, maybe vs Z instead?
      call ezplot0(title, titll, titlr, rhon_half, gradprlb2_avg,
     &    nnoderho_half, nrhomax)


      title= '<R*B_{pol}>'
      titll= 'capr_bpol (mT)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, capr_bpol_midavg,
     &    nnoderho_half, nrhomax)


      !-------------- plot Plasma Profiles -----------
      call a1mnmx(capr, nxmx, nnodex, caprmin, caprmax)
      caprminp = caprmin * 1.01
      caprmaxp = caprmax * .99
      call a2mnmx(xnmid, nxmx, nnodex,
     &   capr, caprminp, caprmaxp, xnmin, xnmax)

      xnmax = 2.0 * xnmax

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (caprmin, caprmax, xnmin, xnmax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BNST', 0.0, 0)
      call pgmtxt('b', 3.2, 0.5, 0.5, 'R (m)')
      call pgmtxt('t', 2.0, 0.5, 0.5, 'Plasma profiles')

      CALL PGSCI(nblue)
      call pgmtxt('l', 2.0, 0.5, 0.5, ' n [m^{-3}]')
      call pgline(nnodex, capr, xnmid)

      call a1mnmx(xktemid, nxmx, nnodex, temin, temax)
      call a1mnmx(xktimid, nxmx, nnodex, timin, timax)
      tmax = max(timax, temax)
      tmin = 0.0

      CALL PGSWIN (caprmin, caprmax, tmin, tmax)
      CALL PGSCI(nblack)
      CALL PGBOX  (' ', 0.0, 0, 'CMST', 0.0, 0)

      CALL PGSCI(nred)
      call pgsls(2)
      call pgline(nnodex, capr, xktemid)
      call PGMTXT ('r', 2.0, 0.5, 0.5, 'T (eV)')

      CALL PGSCI(ngreen)
      call pgsls(3)
      call pgline(nnodex, capr, xktimid)

      call a1mnmx(xjxmid, nxmx, nnodex, xjxmin, xjxmax)
      call a1mnmx(xjymid, nxmx, nnodex, xjymin, xjymax)
      call a1mnmx(xjzmid, nxmx, nnodex, xjzmin, xjzmax)

      CALL PGSWIN (caprmin, caprmax, xjymin, xjymax)
      CALL PGSCI(nblack)
      call pgsls(1)
      call pgline(nnodex, capr, xjxmid)
      call pgline(nnodex, capr, xjymid)
      call pgline(nnodex, capr, xjzmid)

      CALL PGSCI(nblack)
      call pgsls(1)


c
c--plot midplane density vs capr on zoomed grid
c

      call a2mnmx(xnmid, nxmx, nnodex,
     &   capr, rmin_zoom, rmax_zoom, xnmin, xnmax)

      title= 'Edge density in midplane'
      titll= 'density (m^{-3})'
      titlr= ''
      titlb = 'R (m)'

c      call ezzoom1(title, titll, titlr, titlb, capr, xnmid,
c     &    fmid1, fmid2, fmid3, fmid4, fmid5,
c     &    nnodex, nxmx, xnmin, xnmax, rmin_zoom, rmax_zoom)

      call ezzoom6(title, titll, titlr, titlb, capr, xnmid,
     &    xktemid, xktimid, xjymid, fmid4, fmid5,
     &    nnodex, nxmx, rmin_zoom, rmax_zoom) !jcw other J components?

*     -----------------------
*     plot kphi = nphi / R(i) xxxx
*     -----------------------
      title = 'k_{\phi} = n_{phi} / R'
      titll = 'kphi [m^{-1}]'
      call ezplot1q(title, titll, titlr, titlb, capr, xkphi,
     &   nnodex, nxmx)


      title= 'Flux surface average density'
      titll= 'density [m^{-3}]'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xnavg,
     &   nnoderho_half, nrhomax)

      
      title= 'Differential volume element'
      titll= 'dVol (m3)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, dvol,
     &                                          nnoderho_half, nrhomax)


      title= 'Integrated volume'
      titll= 'volume (m3)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon, volume, nnoderho, nrhomax)



      do n = 1, nnoderho
         xkteavg(n) = xkteavg(n) / q
         xktiavg(n) = xktiavg(n) / q
         xkti2avg(n) = xkti2avg(n) / q
         xkti3avg(n) = xkti3avg(n) / q
      end do

      title= '<electron temperature>'
      titll= 'kTe (eV)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xkteavg,
     &    nnoderho_half, nrhomax)

      title= 'Flux average ion temperature'
      titll= 'kTi (eV)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xktiavg,
     &    nnoderho_half, nrhomax)


      title= 'Flux average beam density'
      titll= 'density (m^{-3})'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xn3avg,
     &   nnoderho_half, nrhomax)


      title= 'Flux average beam temperature'
      titll= 'kTi_{fast} (eV)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xkti3avg,
     &    nnoderho_half, nrhomax)


      title= 'Flux average minority density'
      titll= 'density (m^{-3})'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xn2avg,
     &   nnoderho_half, nrhomax)

      title= 'Flux average majority density'
      titll= 'density (m^{-3})'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xn1avg,
     &   nnoderho_half, nrhomax)

      title= 'Flux average minority temperature'
      titll= 'kT2 (eV)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xkti2avg,
     &    nnoderho_half, nrhomax)

      title = 'Contour plot of density'
      call ezconc(capr, y, xn, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title= 'Flux surface average redotj'
      titll= 'redotj (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half,
     &   redotj1avg, redotj2avg, redotj3avg, redotj4avg, redotj5avg,
     &   redotj6avg, redotjeavg,  nnoderho_half, nrhomax)

      !redo with > 0 only
      title= 'Flux surface average redotj'
      titll= 'redotj (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot70(title, titll, titlr, titlb, rhon_half,
     &   redotj1avg, redotj2avg, redotj3avg, redotj4avg, redotj5avg,
     &   redotj6avg, redotjeavg,  nnoderho_half, nrhomax)

      title= 'Integrated redotj'
      titll= 'P (Watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon,
     &   redotj1avg_int, redotj2avg_int, redotj3avg_int, redotj4avg_int,
     &   redotj5avg_int, redotj6avg_int, redotjeavg_int,
     &   nnoderho, nrhomax)

c      write(15, *)
c      write(15, *) 'Flux surface driven current'
c      write(15, *)
c      write(15, *) '        n      rho       J (A/m2)     I (A)'
c      write(15, *)

c      do n = 1, nnoderho_half
c        write (15, 1312) n, rhon_half(n), xjprlavg(n), xjprl_int(n)
c        write (6, 1312) n, rhon_half(n), xjprlavg(n), xjprl_int(n)
c      end do

      title= 'Flux surface driven current'
      titll= 'xjprl (Amps/m2)'
      titlr= 'I (Amps)'
      titlb= 'rho'

      call ezplot2(title, titll, titlr, titlb, rhon_half, xjprlavg,
     &    xjprl_int, nnoderho_half, nrhomax)

      if(xjprl_int(nnoderho) .lt. 0.0) then


         titll= '-xjprl (Amps/m2)'
         titlr= '-I (Amps)'

         xjprlavg = - xjprlavg
         xjprl_int = - xjprl_int

         call ezplot2(title, titll, titlr, titlb, rhon_half, xjprlavg,
     &      xjprl_int, nnoderho_half, nrhomax)

         titll= 'xjprl (A/cm2/MW)'

         xjprlavg = 0.0
         if (prfin .ne. 0.0)
     &      xjprlavg = xjprlavg  / prfin * 1.0e+06 * 1.0e-04

         call ezplot1(title, titll, titlr, rhon_half, xjprlavg,
     &      nnoderho_half, nrhomax)

         titll= 'I (kA)'
         xjprl_int  = xjprl_int / 1.0e+03
         call ezplot1(title, titll, titlr, rhon_half, xjprl_int,
     &      nnoderho_half, nrhomax)

      end if

      if(xjprl_int(nnoderho) .gt. 0.0) then

         titll= 'xjprl (Amps/m2)'
         titlr= 'I (Amps)'

         xjprlavg = xjprlavg
         xjprl_int = xjprl_int

         titll= 'xjprl (A/cm2/MW)'

         xjprlavg = 0.0
         if (prfin .ne. 0.0)
     &      xjprlavg = xjprlavg  / prfin * 1.0e+06 * 1.0e-04

         call ezplot1(title, titll, titlr, rhon_half, xjprlavg,
     &      nnoderho_half, nrhomax)

         titll= 'I (kA)'
         xjprl_int  = xjprl_int / 1.0e+03
         call ezplot1(title, titll, titlr, rhon_half, xjprl_int,
     &      nnoderho_half, nrhomax)

      end if

      title= 'Flux surface average wdot'
      titll= 'wdot (Watts/m3)'
      titlr= '       '
      titlb= 'rho'

      ! ------- plot with ztable usage -------------!
      if (nzeta_wdoti/=0) then
      call ezplot7(title, titll, titlr, titlb, rhon_half,
     &   wdoti1avg, wdoti2avg, wdoti3avg, wdoti4avg, wdoti5avg,
     &   wdoti6avg, wdoteavg, nnoderho_half, nrhomax)

      title= 'Flux surface average wdote'
      titll= 'wdote (W/cm3/MW)'
      titlr= '       '
      titlb= 'rho'

!      wdoteavg = 0.0  !JCW why?

!      if (prfin .ne. 0.0)
!     &   wdoteavg = wdoteavg / prfin * 1.0e+06 * 1.0e-06

!      call ezplot1_red(title, titll, titlr, rhon_half, wdoteavg,
!     &    nnoderho_half, nrhomax)

!      title= 'Integrated wdot'
!      titll= 'P (Watts)'
!      titlr= '       '
!      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon,
     &   wdoti1avg_int, wdoti2avg_int, wdoti3avg_int, wdoti4avg_int,
     &   wdoti5avg_int, wdoti6avg_int, wdoteavg_int, nnoderho, nrhomax)

      title= 'Flux surface average wdoti1'
      titll= 'wdoti1 (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot1(title, titll, titlr, rhon_half, wdoti1avg,
     &    nnoderho_half, nrhomax)


      title= 'Flux surface average wdoti2'
      titll= 'wdoti2 (watts/m3)'
      titlr= '       '
      titlb= 'rho'
      
      call ezplot1(title, titll, titlr, rhon_half, wdoti2avg,
     &    nnoderho_half, nrhomax)

      end if
      !-------------------end!
      
      title= 'Flux average Toroidal force'
      titll= 'force (Nt/m^{3})'
      titlr='       '

      call ezplot7(title, titll, titlr, titlb, rhon_half,
     &    fz0i1avg, fz0i2avg, fz0i3avg, fz0i4avg, fz0i5avg,
     &    fz0i6avg, fz0eavg, nnoderho_half, nrhomax)

      !JCW gets to here
      
      title= 'Flux average toroidal force'
      titll= 'force (Nt/m3)'
      titlr= 'integrated force (Nt/m)'
      titlb= 'rho'

      call ezplot2(title, titll, titlr, titlb, rhon_half, fz0avg,
     &    fz0_int, nnoderho_half, nrhomax)


      title= 'Flux surface average Vy'
      titll= 'Vy (m/s)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, vyavg,
     &    nnoderho_half, nrhomax)

*     -------------------
*     plot G(psi) in 2-D
*     -------------------

      title = 'Toroidal angular speed G(psi)'
      call ezconc(capr, y, gpsi, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

*     -----------------------
*     plot G(psi) in midplane
*     -----------------------

      do i = 1, nnodex
         fmidre(i) = gpsi(i, jmid)
      end do

      title = 'Toroidal angular speed G(psi)'
      titll = 'Gpsi (sec-1)'
      titlr = '      '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre,
     &   nnodex, nxmx)


*     -------------------
*     plot K(psi) in 2-D
*     -------------------

      title = 'Poloidal speed/B_{pol} = K(psi)'
      call ezconc(capr, y, kpsi, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

*     -----------------------
*     plot K(psi) in midplane
*     -----------------------

      do i = 1, nnodex
         fmidre(i) = kpsi(i, jmid)
      end do

      title = 'Poloidal speed / B_pol = K(psi)'
      titll = 'K(psi) (m/s/T)'
      titlr = '      '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre,
     &   nnodex, nxmx)


*     -----------------
*     plot uzeta in 2-D
*     -----------------

      title = 'Toroidal velocity (uzeta)'
      call ezconc(capr, y, uzeta, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

*     -----------------------
*     plot uzeta in midplane
*     -----------------------

      do i = 1, nnodex
         fmidre(i) = uzeta(i, jmid)
      end do

      title = 'Toroidal velocity (uzeta)'
      titll = 'uzeta (m/sec)'
      titlr = '  '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre,
     &   nnodex, nxmx)




*     -----------------
*     plot utheta in 2-D
*     -----------------

      title = 'Poloidal velocity (utheta)'
      call ezconc(capr, y, utheta, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

*     -----------------------
*     plot utheta in midplane
*     -----------------------
      do i = 1, nnodex
         fmidre(i) = utheta(i, jmid)
      end do

      title = 'Poloidal velocity (utheta)'
      titll = 'utheta (m/sec)'
      titlr = '  '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre,
     &   nnodex, nxmx)









*     -------------------
*     plot omgexb in 2-D
*     -------------------

      title = 'ExB shearing rate (omgexb)'
      call ezconc(capr, y, omgexb, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

*     -----------------------
*     plot omgexb in midplane
*     -----------------------

      do i = 1, nnodex
         fmidre(i) = omgexb(i, jmid)
      end do

      title = 'ExB shearing rate (omgexb)'
      titll = 'omgexb (sec-1)'
      titlr = '      '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre,
     &   nnodex, nxmx)



      title= 'xjhat(rho)'
      titll= 'xjhat (kg m-2 s-1)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xjhat,
     &    nnoderho_half, nrhomax)

      if (eqtype=='tokamak') then !these are NaN for mirror, check
      title= 'xkhat(rho)'
      titll= 'xkhat (kg m-2 s-1)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xkhat,
     &    nnoderho_half, nrhomax)

      title= 'poloidal speed <K(rho)>'
      titll= 'K(rho) (m s-1 T-1)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, kpsi_avg,
     &    nnoderho_half, nrhomax)


      title= 'toroidal angular speed <G(\rho)>'
      titll= 'G(\rho) (s^{-1})'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, gpsi_avg,
     &    nnoderho_half, nrhomax)
      end if
      
      title= 'normalized viscosity <mu hat>'
      titll= 'mu hat (kg/m**3/s-1)'
      titlr='       '
      call ezplot0(title, titll, titlr, rhon_half, muhat_avg,
     &    nnoderho_half, nrhomax)


      title= 'collisionality (nu^{*})'
      titll= 'nu*'
      titlr='       '
      call ezplot0(title, titll, titlr, rhon_half, nu_star_avg,
     &    nnoderho_half, nrhomax)




*     -------------------
*     plot nu_star in 2-D
*     -------------------

      title = 'collisionality  (nu^{*})'
      call ezconc(capr, y, nu_star, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



c
c--plot edotj in midplane
c
      do i = 1, nnodex
         fmidre(i) = redotje(i,jmid)
         fmidim(i) = redotj2(i,jmid)
         fmid1(i)  = redotj1(i,jmid)
         fmid2(i)  = redotj2(i,jmid)
         fmid3(i)  = redotj3(i,jmid)
         fmid4(i)  = redotj4(i,jmid)
         fmid5(i)  = redotj5(i,jmid)
         fmidt(i)  = redotjt(i,jmid)
      end do

      title= 'J dot E in midplane'
      titll= 'Re edotj (W/m3)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot7(title, titll, titlr, titlb, capr,
     &    fmid1, fmid2, fmid3, fmid4, fmid5, fmid6, fmidre,
     &    nnodex, nxmx)

c
c--plot wdot at jcut
c

      dy = y(2) - y(1)
      jcut = int((ycut - y(1) - dy / 2.0) / dy) + 1

      if (ycut .eq. 0.0)jcut = jmid

      do i = 1, nnodex
         fmidre(i) = wdote(i, jcut)
         fmidim(i) = wdoti2(i, jcut)
         fmid1(i)  = wdoti1(i, jcut)
         fmid2(i)  = wdoti2(i, jcut)
         fmid3(i)  = wdoti3(i, jcut)
         fmid4(i)  = wdoti4(i, jcut)
         fmid5(i)  = wdoti5(i, jcut)
         fmidt(i)  = wdott(i, jcut)
      end do

      title= 'Midplane Wdot'
      titll= 'Wdot (W/m3)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot70(title, titll, titlr, titlb, capr,
     &    fmid1, fmid2, fmid3, fmid4, fmid5, fmid6, fmidre,
     &    nnodex, nxmx)


      do i = 1, nnodex
         fmid1(i) = bmod(i, jcut)
      end do

      title= 'Midplane Mod B'
      titll= 'Mod B (T)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmid1,
     &   nnodex, nxmx)



c
c--plot F_toroidal at y_cut
c
      dy = y(2) - y(1)
      jcut = int((ycut - y(1) - dy / 2.0) / dy) + 1

      if (ycut .eq. 0.0)jcut = jmid

      do i = 1, nnodex
         fmidre(i) = fz0e(i, jcut)
         fmidim(i) = fz0i2(i, jcut)
         fmid1(i)  = fz0i1(i, jcut)
         fmid2(i)  = fz0i2(i, jcut)
         fmid3(i)  = fz0i3(i, jcut)
         fmid4(i)  = fz0i4(i, jcut)
         fmid5(i)  = fz0i5(i, jcut)
         fmidt(i)  = fz0(i, jcut)
      end do

      title= 'Toroidal force'
      titll= 'F_toroidal (Nt/m3)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot7(title, titll, titlr, titlb, capr,
     &    fmid1, fmid2, fmid3, fmid4, fmid5, fmid6, fmidre,
     &    nnodex, nxmx)



c--Contour plots using 16 contour levels:

      if (ibackground.eq.1)then
c        black background
         ncollin = ncyan
         ncolln2 = nblueviolet
      end if

      if (ibackground.eq.0)then
c        white background
         ncollin = nblue
         ncolln2 = nyellow
      end if

      titx = 'R (m)'
      tity = 'Z (m)'

c
c--plot q profile in midplane
c
      call a2mnmx(qmid, nxmx, nnodex,
     &   capr, caprminp, caprmaxp, qmin, qmax)
      qmin = 0.0
      qmax = 5.0

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


      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = rho(i,j)**2
            psi(i,j) = rho(i,j)**2
         end do
      end do

      title = 'Contour plot of psi'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

*     --------------------
*     plot psi in midplane
*     --------------------

      do i = 1, nnodex
         fmidre(i) = fmod(i, jmid)
      end do

      title = 'psi = rho**2 in midplane'
      titll = 'psi'
      titlr = '      '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre,
     &   nnodex, nxmx)

*     ---------------------
*     plot Mod 1/B surfaces
*     ---------------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = 1.0 / bmod(i,j)
         end do
      end do

      title = 'Mod 1/B surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



      if (eqtype/='mirror') then
         title = 'Contour plot of theta'
         call ezconc(capr, y, theta, ff, nnodex, nnodey, 28,
     &        nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
         if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &        nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &        scalex)

      end if
      
c--   Plot xkperp_cold:
      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(xkperp_cold(i,j))
            fimag(i,j) = aimag(xkperp_cold(i,j))
         end do
      end do


      title = 'Re xkperp cold'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, 1, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = 'Im xkperp cold'
      call ezconc(capr, y, fimag, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, 1, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


*     --------------------------------------
*     plot xkperp_cold(i,j) in midplane
*     --------------------------------------
      do i = 1, nnodex
         fmidre(i) = real(xkperp_cold(i,jmid))
         fmidim(i) = aimag(xkperp_cold(i,jmid))
      end do

      title = 'Midplane xkperp cold'
      titll = 'Re xkperp (1/m)'
      titlr = 'Im xkperp (1/m)'
      titlb = 'R (m)'

      call ezplot2(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)

*     --------------------------------------
*     plot xkperp_cold2(i,j) in midplane
*     --------------------------------------
      do i = 1, nnodex
         fmidre(i) = real(xkperp_cold2(i,jmid))
         fmidim(i) = aimag(xkperp_cold2(i,jmid))
      end do

      title = 'Midplane xkperp cold2'
      titll = 'Re xkperp^{2} (1/m)'
      titlr = 'Im xkperp^{2} (1/m)'
      titlb = 'R (m)'

      call ezplot2(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)


c--   Plot xkperp_cold2 in 2D:
      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(xkperp_cold2(i,j))
            fimag(i,j) = aimag(xkperp_cold2(i,j))
         end do
      end do


      title = 'Re xkperp cold2'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, 1, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


c--   Plot Acold in 2D:
      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(acold(i,j))
            fimag(i,j) = aimag(acold(i,j))
         end do
      end do

      title = 'Re Acold (resonance)'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, 1,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, 1, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



c--   Plot Bcold in 2D:(slow wave cutoff)
      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(bcold(i,j))
            fimag(i,j) = aimag(bcold(i,j))
         end do
      end do

      title = 'Re Bcold (slow wave cutoff)'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, 1,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, 1, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



c--   Plot Ccold in 2D:(fast wave cutoff)
      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ccold(i,j))
            fimag(i,j) = aimag(ccold(i,j))
         end do
      end do

      title = 'Re Ccold (fast wave cutoff)'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, 1,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, 1, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

c      call ezconz(capr, y, freal, ff, nnodex, nnodey, 1,
c     &   nxmx, nymx, nlevmax, title, titx, tity)


*     -----------------------
*     plot resonance surfaces
*     -----------------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = reomg1a(i,j)
         end do
      end do

      title = 'Majority resonant surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, 15,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb,nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = reomg2a(i,j)
         end do
      end do

      title = 'Minority resonant surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, 15,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = reomg3a(i,j)
         end do
      end do

      title = 'Species 3 resonant surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, 15,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = reomglha(i,j)
         end do
      end do

      title = 'Lower hybrid resonant surface'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, 1,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      title = 'Contour plot of Janty'
      call ezconc(capr, y, xjy, ff, nnodex, nnodey, 2,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = 'Contour plot of Jantx'
      call ezconc(capr, y, xjx, ff, nnodex, nnodey, 2,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = 'Contour plot of Jantz'
      call ezconc(capr, y, xjz, ff, nnodex, nnodey, 2,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ealpha(i,j))
            fimag(i,j) = aimag(ealpha(i,j))
         end do
      end do


      title = 'Real E alpha'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

c      call ezconz(capr, y, freal, ff, nnodex, nnodey, numb,
c     &   nxmx, nymx, nlevmax, title, titx, tity)


!     --------------------------------------
!     Write 2D.vtk file "structured points"
!     --------------------------------------
      number_points = nnodex * nnodey
      dx = capr(2) - capr(1)
      dy = y(2) - y(1)

      write(140, 2840)
 2840 format('# vtk DataFile Version 2.0')

      write(140, 2845)
 2845 format('Real E_alpha')

      write(140, 2846)
 2846 format('ASCII')

      write(140, 2847)
 2847 format('DATASET STRUCTURED_POINTS')

      write(140, 2841) nnodex,  nnodey, 1
 2841 format('DIMENSIONS', 3i8)

      write(140, 2842) capr(1), y(1), 0
 2842 format('ORIGIN', 2f12.3, 1i8)

      write(140, 2843) dx, dy, 1
 2843 format('SPACING', 2f9.5, 1i8)

      write(140, 2844) number_points
 2844 format('POINT_DATA', i10)

      write(140, 2848)
 2848 format('SCALARS Real_E_alpha float 1')

      write(140, 2849)
 2849 format('LOOKUP_TABLE default')

      write (140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

 3410 format(1p,4e10.2)
 3411 format(6e17.4)


      write(140, 3853)
 3853 format('SCALARS Acold float 1')
      write(140, 2849)
      write (140, 3411)((real(acold(i,j)),i = 1, nnodex), j = 1, nnodey)

      write(140, 6853)
 6853 format('SCALARS Bcold float 1')
      write(140, 2849)
      write (140, 3411)((real(bcold(i,j)),i = 1, nnodex), j = 1, nnodey)

      write(140, 6854)
 6854 format('SCALARS Ccold float 1')
      write(140, 2849)
      write (140, 3411)((real(ccold(i,j)),i = 1, nnodex), j = 1, nnodey)


      title = 'Imag E alpha'
      call ezconc(capr, y, fimag, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ebeta(i,j))
            fimag(i,j) = aimag(ebeta(i,j))
         end do
      end do
      title = 'Real Ebeta'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      title = 'Imag Ebeta'
      call ezconc(capr, y, fimag, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(eb(i,j))
            fimag(i,j) = aimag(eb(i,j))
         end do
      end do
      title = 'Real Eb'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      call pseudo(capr*scalex, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      write(140, 3966)
 3966 format('SCALARS Real_eb float 1')
      write(140, 2849)
      write(140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3967)
 3967 format('SCALARS Imag_eb float 1')
      write(140, 2849)
      write(140, 3411) ((fimag(i,j), i = 1, nnodex),  j = 1, nnodey)


      title = 'Imag Eb'
      call ezconc(capr, y, fimag, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

c
c--plot plasma current (electrons)
c

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(xjpxe(i,j))
            fimag(i,j) = aimag(xjpxe(i,j))
         end do
      end do
      title = 'Real Je_{alpha}'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(xjpye(i,j))
            fimag(i,j) = aimag(xjpye(i,j))
         end do
      end do
      title = 'Real Je_{beta}'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(xjpze(i,j))
            fimag(i,j) = aimag(xjpze(i,j))
         end do
      end do
      title = 'Real Je_{b}'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

c
c--plot E alpha in midplane
c

c      jmid = jdiag

      do i = 1, nnodex
         fmidre(i) = real(ealpha(i,jmid))
         fmidim(i) = aimag(ealpha(i,jmid))
      end do

      title = 'Ealpha'
      titll = 'Re Ealpha (V/m)'
      titlr = 'Im Ealpha (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)

c
c--plot Ebeta in midplane
c
      do i = 1, nnodex
         fmidre(i) = real(ebeta(i,jmid))
         fmidim(i) = aimag(ebeta(i,jmid))
      end do

      title = 'Ebeta'
      titll = 'Re Ebeta (V/m)'
      titlr = 'Im Ebeta (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)


*     ----------------------
*     plot Eplus in midplane
*     ----------------------
      ISQ2 = SQRT(0.5)
      do i = 1, nnodex
         do j = 1, nnodey
            eplus(i,j)  = isq2 * (ealpha(i,j) + zi * ebeta(i,j))

            fmod(i,j)   = conjg(eplus(i,j)) * eplus(i,j)
            mod_Eplus(i,j) = sqrt(fmod(i,j))
         end do
      end do


      do i = 1, nnodex
         fmidre(i) = real(eplus(i,jmid))
         fmidim(i) = aimag(eplus(i,jmid))
         mod_Eplus_mid(i) = mod_Eplus(i,jmid)
      end do

      title = 'Eplus'
      titll = 'Re Eplus (V/m)'
      titlr = 'Im Eplus (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)



*     -----------------------
*     plot Eminus in midplane
*     -----------------------

      do i = 1, nnodex
         do j = 1, nnodey
            eminus(i,j) = isq2 * (ealpha(i,j) - zi * ebeta(i,j))

            fmod(i,j)   = conjg(eminus(i,j)) * eminus(i,j)
            mod_Eminus(i,j) = sqrt(fmod(i,j))
         end do
      end do

      do i = 1, nnodex
         fmidre(i) = real(eminus(i,jmid))
         fmidim(i) = aimag(eminus(i,jmid))
         mod_Eminus_mid(i) = mod_Eminus(i,jmid)
      end do

      title = 'E_{-}'
      titll = 'Re Eminus (V/m)'
      titlr = 'Im Eminus (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)


*     ----------------------
*     plot Eb in midplane
*     ----------------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = conjg(eb(i,j)) * eb(i,j)
            mod_Eb(i,j) = sqrt(fmod(i,j))
c            if(i .eq. 70 .and. j .eq. 64)then
c               write(6, *) 'R(70) = ', capr(i)
c               write(6, *) 'Z(64) = ', y(j)
c               write(6, *) 'mod_Eb(70, 64) = ', mod_Eb(i,j)
c            end if
         end do
      end do

      do i = 1, nnodex
         fmidre(i) = real(eb(i,jmid))
         fmidim(i) = aimag(eb(i,jmid))
         mod_Eb_mid(i) = mod_Eb(i,jmid)
      end do

      title = 'E parallel'
      titll = 'Re E parallel (V/m)'
      titlr = 'Im E parallel (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)


*     -----------------------
*     write bharvey (66) file
*     ------------------------
      write(66,309) nnodex, nnodey
      write(66,310) (capr(i), i = 1, nnodex)
      write(66,310) (y(j), j = 1, nnodey)

      write(66,310) ((eplus(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((eminus(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)

      write(66,310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)

*     -----------------------
*     write zowens (68) file
*     ------------------------
      write(68,309) nnodex, nnodey
      write(68,310) (capr(i), i = 1, nnodex)
      write(68,310) (y(j), j = 1, nnodey)

      write(68,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      write(68,310) ((bratio(i, j), i = 1, nnodex), j = 1, nnodey)



      ncollin = nblue
      ncolln2 = nblueviolet

      write(140, 8969)
 8969 format('SCALARS redotj float 1')
      write(140, 2849)
      write(140, 3411) ((redotjt(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 8949)
 8949 format('SCALARS redotje float 1')
      write(140, 2849)
      write(140, 3411) ((redotje(i,j), i = 1, nnodex),  j = 1, nnodey)


      write(140, 5949)
 5949 format('SCALARS redotj1 float 1')
      write(140, 2849)
      write(140, 3411) ((redotj1(i,j), i = 1, nnodex),  j = 1, nnodey)


      write(140, 3949)
 3949 format('SCALARS redotj2 float 1')
      write(140, 2849)
      write(140, 3411) ((redotj2(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 13949)
13949 format('SCALARS redotj3 float 1')
      write(140, 2849)
      write(140, 3411) ((redotj3(i,j), i = 1, nnodex),  j = 1, nnodey)
c
c***********************************************************
c
      ncollin = nred
      ncolln2 = nmagenta

      title = '1/2 Re JedotE'
      call ezconc(capr, y, redotje, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



c
c***********************************************************
c
      ncollin = ncyan
      ncolln2 = nblueviolet

      title = '1/2 Re J1dotE'
      call ezconc(capr, y, redotj1, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



      title = '1/2 Re J2dotE'

      call ezconc(capr, y, redotj2, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)

      if (iflag .eq. 0) call boundary  (capr, y, rho, ff,
     &     nnodex, nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

c
c**********************************************************
c
      ncollin = npink
      ncolln2 = ngrey

      title = '1/2 Re J3dotE'
      call ezconc(capr, y, redotj3, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


c
c***********************************************************
c
      ncollin = ngreen
      ncolln2 = nyellow

      title = '1/2 Re J4dotE'
      call ezconc(capr, y, redotj4, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = '1/2 Re J5dotE'
      call ezconc(capr, y, redotj5, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = '1/2 Re J6dotE'
      call ezconc(capr, y, redotj6, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

c
c***********************************************************
c

      ncollin = ngreen
      ncolln2 = nyellow

      title = '1/2 Re JsdotE'
      call ezconc(capr, y, redotjs, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

c
c***********************************************************
c
      if (ibackground.eq.1)then !  black background
         ncollin = nblue
         ncolln2 = nblueviolet
      end if

      if (ibackground.eq.0)then ! white background
         ncollin = nblue
         ncolln2 = nyellow
      end if

      title = '1/2 Re JdotE'
      call ezconc(capr, y, redotjt, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = 'div Q'
      call ezconc(capr, y, divq, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      if (dfquotient .eq. 0.0) then
         numb = 20
      else
         numb = 1
      end if

      titx='R (m)'
      tity='kx (1/m)'

      title='Contour plot of D(x, kx)'

      call ezconc2(capr, xkxsav, capd, ff, nnodex, nkxplt, numb,
     &   nxmx, nkpltdim, nlevmax, title, titx, tity, iflag, dfquotient)


c      CALL PGENV(0., 1.0, 0., 1.0, 1, -2)
c      CALL FREDDY(capd, nnodex, nkxplt, 1.0, 25.0)

      title='Contour plot of xkb(x, kx)'

      call ezconc2(capr, xkxsav, xkb, ff, nnodex, nkxplt, numb,
     &   nxmx, nkpltdim, nlevmax, title, titx, tity, iflag, dfquotient)


!        ---------------------------------------------
!        Write capd.vtk file "structured points"
!        ---------------------------------------------
         capd_plot = capd * 1.0e-15
         xkb_plot = xkb
         number_points = nnodex * nkxplt
         dx = capr(2) - capr(1)
         dy = xkxsav(2) - xkxsav(1)
         write(143, 2840)
         write(143, 4861)
 4861    format('Contour plot of D(x, kx)')
         write(143, 2846)
         write(143, 2847)
         write(143, 2841) nnodex,  nkxplt, 1
         write(143, 2842) capr(1), xkxsav(1), 0
         write(143, 2843) dx, dy, 1
         write(143, 2844) number_points

         write(143, 4862)
 4862    format('SCALARS capd float 1')
         write(143, 2849)
         write(143, 3411) ((capd_plot(i,j), i = 1, nnodex),
     &                                         j = 1, nkxplt)


         write(143, 4873)
 4873    format('SCALARS xkb float 1')
         write(143, 2849)
         write(143, 3411) ((xkb_plot(i,j), i = 1, nnodex),
     &                                         j = 1, nkxplt)

      numb = 20


      ncollin=nred
      ncolln2=nmagenta

      title = 'wdote'
      call ezconc(capr, y,wdote, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



      write(140, 3852)
 3852 format('SCALARS wdote float 1')
      write(140, 2849)
      write (140, 3411) ((wdote(i,j), i = 1, nnodex),  j = 1, nnodey)



      if (ibackground.eq.1)then
c        black background
         ncollin = ncyan
         ncolln2 = nblueviolet
      end if

      if (ibackground.eq.0)then
c        white background
         ncollin = nblue
         ncolln2 = nyellow
      end if

      title = 'wdoti1'
      call ezconc(capr, y,wdoti1, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      title = 'wdoti2'
      call ezconc(capr, y,wdoti2, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      write(140, 3849)
 3849 format('SCALARS wdoti2 float 1')
      write(140, 2849)
      write (140, 3411) ((wdoti2(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3850)
 3850 format('SCALARS wdoti1 float 1')
      write(140, 2849)
      write (140, 3411) ((wdoti1(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 13859)
13859 format('SCALARS wdoti3 float 1')
      write(140, 2849)
      write (140, 3411) ((wdoti3(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 13860)
13860 format('SCALARS wdoti4 float 1')
      write(140, 2849)
      write (140, 3411) ((wdoti4(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3869)
 3869 format('SCALARS wdot_tot float 1')
      write(140, 2849)
      write (140, 3411) ((wdott(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3851)
 3851 format('SCALARS rho float 1')
      write(140, 2849)
      write (140, 3411) ((rho(i,j), i = 1, nnodex),  j = 1, nnodey)





      title = 'wdoti3'
      call ezconc(capr, y,wdoti3, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = 'wdoti4'
      call ezconc(capr, y,wdoti4, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = 'wdoti5'
      call ezconc(capr, y,wdoti5, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      title = 'wdot_tot'
      call ezconc(capr, y,wdott, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


#ifdef skip
      title = 'Toroidal force'
      call ezconc(capr, y, fz0, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)
#endif


      write(140, 3857)
 3857 format('SCALARS fz0 float 1')
      write(140, 2849)
      write (140, 3411) ((fz0(i,j), i = 1, nnodex),  j = 1, nnodey)



      titx = 'kx (m-1)'
      tity = 'ky (m-1)'
      title = 'Mod E alpha(kx, ky)'
      call ezconc(xkxsav, xkysav, exkmod, ff, nkxplt, nkyplt, numb,
     &     nkpltdim, mkpltdim, nlevmax, title, titx, tity,
     &     iflag, 1.0)

      title='3-D plot of Mod E alpha(kx,ky)'
      titz='Mod E alpha(kx,ky)'
      call ezcon3d(xkxsav, xkysav, exkmod, ff, nkxplt, nkyplt, numb,
     &   nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)


      nmid = nkxplt / 2
      do m = 1, nkyplt
         fmodm(m) = exklog(nmid, m)
      end do
c
c--plot exklog vs ky
c
      title= 'Mod E alpha(ky) Spectrum'
      titll= 'Mod E alpha (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, tity, xkysav, fmodm,
     &    nkyplt, mkpltdim, exkmin, exkmax)







      mmid = nkyplt / 2
      do n = 1, nkxplt
         fmodm(n) = exklog(n, mmid)
      end do

c
c--plot exklog vs kx
c
      title= 'Mod E alpha(kx) Spectrum'
      titll= 'Mod E alpha (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, titx, xkxsav, fmodm,
     &    nkxplt, nkpltdim, exkmin, exkmax)




c--   plot log of E alpha

      titx = 'kx (m-1)'
      tity = 'ky (m-1)'
      title = 'log E alpha(kx, ky)'
      call ezconc(xkxsav, xkysav, exklog, ff, nkxplt, nkyplt, numb,
     &     nkpltdim, mkpltdim, nlevmax, title, titx, tity,
     &     iflag, scalex)

      title='3-D plot of log E alpha(kx,ky)'
      titz='log E alpha(kx,ky)'
      call ezcon3d(xkxsav, xkysav, exklog, ff, nkxplt, nkyplt, numb,
     &   nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)



      title = 'Mod E beta(kx, ky)'
      call ezconc(xkxsav, xkysav, eykmod, ff, nkxplt, nkyplt, numb,
     &     nkpltdim, mkpltdim, nlevmax, title, titx, tity,
     &     iflag, scalex)

      title='3-D plot of Mod E beta(kx,ky)'
      titz='Mod E beta(kx,ky)'
      call ezcon3d(xkxsav, xkysav, eykmod, ff, nkxplt, nkyplt, numb,
     &   nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)


      do m = 1, nkyplt
         fmodm(m) = eyklog(nmid, m)
      end do
c
c--plot eyklog vs ky
c
      title= 'Mod E beta(ky) Spectrum'
      titll= 'Mod E beta (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, tity, xkysav, fmodm,
     &    nkyplt, mkpltdim, exkmin, exkmax)



      mmid = nkyplt / 2
      do n = 1, nkxplt
         fmodm(n) = eyklog(n, mmid)
      end do

c
c--plot eyklog vs kx
c
      title= 'Mod E beta(kx) Spectrum'
      titll= 'Mod E beta (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, titx, xkxsav, fmodm,
     &    nkxplt, nkpltdim, exkmin, exkmax)



c--   plot log of E beta

      title = 'log E beta(kx, ky)'
      call ezconc(xkxsav, xkysav, eyklog, ff, nkxplt, nkyplt, numb,
     &     nkpltdim, mkpltdim, nlevmax, title, titx, tity,
     &     iflag, scalex)

      title='3-D plot of log E beta(kx,ky)'
      titz='log E beta(kx,ky)'
      call ezcon3d(xkxsav, xkysav, eyklog, ff, nkxplt, nkyplt, numb,
     &   nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)



      title = 'Mod E b(kx, ky)'
      call ezconc(xkxsav, xkysav, ezkmod, ff, nkxplt, nkyplt, numb,
     &     nkpltdim, mkpltdim, nlevmax, title, titx, tity,
     &     iflag, scalex)


!        ---------------------------------------------
!        Write Eb_spectrum.vtk file "structured points"
!        ---------------------------------------------
         number_points = nkxplt * nkyplt
         dx = xkxsav(2) - xkxsav(1)
         dy = xkysav(2) - xkysav(1)
         write(144, 2840)
         write(144, 4863)
 4863    format('Contour plot of Eb(kx, ky)')
         write(144, 2846)
         write(144, 2847)
         write(144, 2841) nkxplt,  nkyplt, 1
         write(144, 2842) xkxsav(1), xkysav(1), 0
         write(144, 2843) dx, dy, 1
         write(144, 2844) number_points

         write(144, 4864)
 4864    format('SCALARS Eb_spectrum float 1')
         write(144, 2849)
         write(144, 3411) ((ezkmod(i,j), i = 1, nkxplt),
     &                                         j = 1, nkyplt)

         write(144, 4865)
 4865    format('SCALARS Ex_spectrum float 1')
         write(144, 2849)
         write(144, 3411) ((exkmod(i,j), i = 1, nkxplt),
     &                                         j = 1, nkyplt)

         write(144, 4866)
 4866    format('SCALARS Ey_spectrum float 1')
         write(144, 2849)
         write(144, 3411) ((eykmod(i,j), i = 1, nkxplt),
     &                                         j = 1, nkyplt)

      title='3-D plot of Mod E b(kx,ky)'
      titz='Mod E b(kx,ky)'
      call ezcon3d(xkxsav, xkysav, ezkmod, ff, nkxplt, nkyplt, numb,
     &   nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)


      do m = 1, nkyplt
         fmodm(m) = ezklog(nmid, m)
      end do
c
c--plot ezklog vs ky
c
      title= 'Mod Eb(ky) Spectrum'
      titll= 'Mod Eb (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, tity, xkysav, fmodm,
     &    nkyplt, mkpltdim, exkmin, exkmax)



      mmid = nkyplt / 2
      do n = 1, nkxplt
         fmodm(n) = ezklog(n, mmid)
      end do

c
c--plot ezklog vs kx
c
      title= 'Mod Eb(kx) Spectrum'
      titll= 'Mod Eb (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, titx, xkxsav, fmodm,
     &    nkxplt, nkpltdim, exkmin, exkmax)


c--   plot log of Eb

      title = 'log E b(kx, ky)'
      call ezconc(xkxsav, xkysav, ezklog, ff, nkxplt, nkyplt, numb,
     &     nkpltdim, mkpltdim, nlevmax, title, titx, tity,
     &     iflag, scalex)

      title='3-D plot of log E b(kx,ky)'
      titz='Mod E b(kx,ky)'
      call ezcon3d(xkxsav, xkysav, ezklog, ff, nkxplt, nkyplt, numb,
     &   nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)





      titx = 'R (m)'
      tity = 'Z (m)'
      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = conjg(ealpha(i,j)) * ealpha(i,j)
            mod_Ealpha(i,j) = sqrt(fmod(i,j))
         end do
      end do

      title = 'Mod E alpha'
      call ezconc(capr, y, mod_Ealpha, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      title='3-D plot of Mod E alpha'
      titz='Mod E alpha'
      call ezcon3d(capr, y, mod_Ealpha, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, titz)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = conjg(ebeta(i,j)) * ebeta(i,j)
            mod_Ebeta(i,j) = sqrt(fmod(i,j))
         end do
      end do

      title = 'Mod E beta'
      call ezconc(capr, y, mod_Ebeta, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title='3-D plot of Mod E beta'
      titz='Mod Ebeta'
      call ezcon3d(capr, y, mod_Ebeta, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, titz)

*     ------------
*     plot mod(Eb)
*     ------------


      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = conjg(eb(i,j)) * eb(i,j)
            mod_Eb(i,j) = sqrt(fmod(i,j))
         end do
      end do

      title = ' Mod Eb'
      call ezconc(capr, y, mod_Eb, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title='3-D plot of Mod Eb'
      titz='Mod Eb'
      call ezcon3d(capr, y, mod_Eb, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, titz)

      write(140, 3956)
 3956 format('SCALARS mod_eb float 1')
      write(140, 2849)
      write(140, 3411) ((mod_Eb(i,j), i = 1, nnodex),  j = 1, nnodey)

*     ----------------
*     calculate mod(E)
*     ----------------
      do i = 1, nnodex
         do j = 1, nnodey
            mod_E(i,j)   = sqrt(conjg(ealpha(i,j)) * ealpha(i,j)
     &                        + conjg(ebeta(i,j)) * ebeta(i,j)
     &                        + conjg(eb(i,j)) * eb(i,j) )
         end do
      end do

      do i = 1, nnodex
         mod_E_mid(i) = mod_E(i,jmid)
      end do

      write(140, 3996)
 3996 format('SCALARS mod_e float 1')
      write(140, 2849)
      write(140, 3411) ((mod_E(i,j), i = 1, nnodex),  j = 1, nnodey)


*     ----------------
*     plot mod(Eplus)
*     ----------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = conjg(eplus(i,j)) * eplus(i,j)
            mod_Eplus(i,j) = sqrt(fmod(i,j))
c            if(i .eq. 70 .and. j .eq. 64)then
c               write(6, *) 'R(70) = ', capr(i)
c               write(6, *) 'Z(64) = ', y(j)
c               write(6, *) 'mod_Eplus(70, 64) = ', mod_Eplus(i,j)
c            end if
         end do
      end do

      title = ' Mod(E_{+})'
      call ezconc(capr, y, mod_Eplus, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      title='3-D plot of Mod Eplus'
      titz='Mod Eplus'
      call ezcon3d(capr, y, mod_Eplus, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, titz)


      write(140, 3954)
 3954 format('SCALARS mod_eplus float 1')
      write(140, 2849)
      write(140, 3411) ((mod_Eplus(i,j), i = 1, nnodex),  j = 1, nnodey)

*     ----------------
*     plot mod(Eminus)
*     ----------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = conjg(eminus(i,j)) * eminus(i,j)
            mod_Eminus(i,j) = sqrt(fmod(i,j))
c            if(i .eq. 70 .and. j .eq. 64)then
c               write(6, *) 'R(70) = ', capr(i)
c               write(6, *) 'Z(64) = ', y(j)
c               write(6, *) 'mod_Eminus(70, 64) = ', mod_Eminus(i,j)
c            end if
         end do
      end do

      title = ' Mod(Eminus)'
      call ezconc(capr, y, mod_Eminus, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      title='3-D plot of Mod Eminus'
      titz='Mod Eminus'
      call ezcon3d(capr, y, mod_Eminus, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, titz)


      write(140, 3955)
 3955 format('SCALARS mod_eminus float 1')
      write(140, 2849)
      write(140, 3411) ((mod_Eminus(i,j), i = 1, nnodex),  j = 1,nnodey)


      write(64,309) nnodex, nnodey
      write(64,310) (capr(i), i = 1, nnodex)
      write(64,310) (y(j), j = 1, nnodey)
      write(64,310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
      write(64,310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
      write(64,310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)


*     ---------------
*     plot real Eplus
*     ---------------

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j)   = real(eplus(i,j))
            fimag(i,j)   = aimag(eplus(i,j))
         end do
      end do

      title = 'real(E_{+})'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      write(140, 3950)
 3950 format('SCALARS re_eplus float 1')
      write(140, 2849)
      write(140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3951)
 3951 format('SCALARS im_eplus float 1')
      write(140, 2849)
      write(140, 3411) ((fimag(i,j), i = 1, nnodex),  j = 1, nnodey)



*     ----------------
*     plot real Eminus
*     ----------------

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j)   = real(eminus(i,j))
            fimag(i,j)   = aimag(eminus(i,j))
         end do
      end do

      title = 'Real(E_{minus})'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      write(140, 3952)
 3952 format('SCALARS re_eminus float 1')
      write(140, 2849)
      write(140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3953)
 3953 format('SCALARS im_eminus float 1')
      write(140, 2849)
      write(140, 3411) ((fimag(i,j), i = 1, nnodex),  j = 1, nnodey)


      write(140, 3854)
 3854 format('SCALARS omgexb float 1')
      write(140, 2849)
      write (140, 3411) ((omgexb(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3855)
 3855 format('SCALARS uzeta float 1')
      write(140, 2849)
      write (140, 3411) ((uzeta(i,j), i = 1, nnodex),  j = 1, nnodey)


      title = 'B toroidal'
      call ezconc(capr, y, bzeta, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

c      title='3-D plot of B toroidal'
c      titz='Btor'
c      call ezcon3d(capr, y, bzeta, ff, nnodex, nnodey, numb,
c     &   nxmx, nymx, nlevmax, title, titx, tity, titz)

c      do i = 1, nnodex
c         write(6, 2163)i, capr(i), x(i), btau(i, 16),
c     &      bzeta(i,16)
c      end do
 2163 format(i5, 1p,11e12.3)



      title = 'B poloidal'
      call ezconc(capr, y, btau, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title='3-D plot of B poloidal'
      titz='Bpol'
      call ezcon3d(capr, y, btau, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, titz)



c      title='3-D plot of safety factor'
c      titz='q'
c      call ezcon3d(capr, y, qsafety, ff, nnodex, nnodey, numb,
c     &   nxmx, nymx, nlevmax, title, titx, tity, titz)





c      title = 'Spx'
c      call ezconc(capr, y, spx, ff, nnodex, nnodey, numb,
c     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
c      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
c     &   nnodey, numb,  nxmx, nymx, nlevmax, title, titx, tity,
c     &     scalex)

c      title = 'Spy'
c      call ezconc(capr, y, spy, ff, nnodex, nnodey, numb,
c     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
c      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
c     &   nnodey, numb,  nxmx, nymx, nlevmax, title, titx, tity,
c     &     scalex)

c      title = 'Spz'
c      call ezconc(capr, y, spz, ff, nnodex, nnodey, numb,
c     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
c      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
c     &   nnodey, numb,  nxmx, nymx, nlevmax, title, titx, tity,
c     &     scalex)


*     --------------
*     read from rf2x
*     --------------

      read(38, 309) nnodex, nnodey
      read(38, 310) rhoplasm
      read(38, 310) (x(i), i = 1, nnodex)
      read(38, 310) (y(j), j = 1, nnodey)
      read(38, 310) (capr(i), i = 1, nnodex)

      read(38, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((theta(i, j), i = 1, nnodex), j = 1, nnodey)


      read(38, 309) nnoderho2
      read(38, 310) (rhon(n), n = 1, nnoderho2)
c      read(38, 310) (gradprlb2_avg(n), n = 1, nnoderho2)

      read(38, 310) (dvol(n), n = 1, nnoderho2)
      read(38, 310) (volume(n), n = 1, nnoderho2)

      read(38, 310) (redotjeavg(n), n = 1, nnoderho2)
      read(38, 310) (redotj1avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj2avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj3avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj4avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj5avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj6avg(n), n = 1, nnoderho2)
      read(38, 310) (redotjiavg(n), n = 1, nnoderho2)

      read(38, 310) (xjprlavg(n), n = 1, nnoderho2)
      read(38, 310) (xjprl_int(n), n = 1, nnoderho2)

      read(38, 310) (wdoteavg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti1avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti2avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti3avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti4avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti5avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti6avg(n), n = 1, nnoderho2)



      read(38, 310) (redotje_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj1_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj2_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj3_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj4_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj5_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj6_int(n), n = 1, nnoderho2)

      read(38, 310) (wdote_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot1_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot2_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot3_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot4_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot5_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot6_int(n), n = 1, nnoderho2)

      read(38, 310) (redotje_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj1_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj2_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj3_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj4_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj5_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj6_dvol(n), n = 1, nnoderho2)


      read(38, 310) (wdote_dvol(n),  n = 1, nnoderho2)
      read(38, 310) (wdoti1_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti2_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti3_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti4_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti5_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti6_dvol(n), n = 1, nnoderho2)



*     --------------
*     read from main
*     --------------

      read(38, 310) (wdote_ql(n),  n = 1, nnoderho)
      read(38, 310) (wdoti1_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti2_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti3_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti4_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti5_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti6_ql(n), n = 1, nnoderho)



      read(38, 310) (dldbavg(n), n = 1, nnoderho)

      read(38, 310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((ntilda_e(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((fpsi0(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((ftheta0(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((pressi(i,j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((dldb_tot12(i, j), i = 1, nnodex),
     &                                  j = 1, nnodey)

      read (38, 309) ndisti1, ndisti2, ndisti3

      read (38, 310) ((eplus_flux_plot(i, j), i = 1, nnodex),
     &                                         j = 1, nnodey)

      read (38, 310) ((eminus_flux_plot(i, j), i = 1, nnodex),
     &                                         j = 1, nnodey)

      read (38, 310) ((xkperp_flux_plot(i, j), i = 1, nnodex),
     &                                         j = 1, nnodey)

      read (38, 310) ((eplus_flux(n, m), n = 1, nnoderho),
     &                                   m = 1, mnodetheta)
      read (38, 310) ((eminus_flux(n, m), n = 1, nnoderho),
     &                                    m = 1, mnodetheta)

      read (38, 310) ((xkperp_flux(n, m), n = 1, nnoderho),
     &                                    m = 1, mnodetheta)


      read (38, 310) ((capr_flux(n,m),   n = 1, nnoderho),
     &                                   m = 1, mnodetheta)
      read (38, 310) ((capz_flux(n,m),   n = 1, nnoderho),
     &                                   m = 1, mnodetheta)

      read (38, 310) ((bmod_flux(n,m),   n = 1, nnoderho),
     &                                   m = 1, mnodetheta)

      read (38, 310) ((dxxuyy(i,j), i = 1, nnodex), j = 1, nnodey)
      read (38, 310) ((dyyuzz(i,j), i = 1, nnodex), j = 1, nnodey)
      read (38, 310) ((gradprlb(i,j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((xkperp2_slow(i, j), i = 1, nnodex),
     &      j = 1, nnodey)
      read(38, 310) ((xkperp2_fast(i, j), i = 1, nnodex),
     &      j = 1, nnodey)
      read(38, 310) ((xkprl_a(i,j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((P_a(i, j), i = 1, nnodex), j = 1, nnodey)


      write(66,310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

      write(66,310) rhoplasm
      write(66,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)




      title= 'Integral dl/B'
      titll= 'Integral dl/B (m)'
      titlr='       '

      call ezplot1_0(title, titll, titlr, rhon_half_save, dldbavg,
     &    nnoderho_half, nrhomax)



      title= 'Quasilinear wdot'
      titll= 'wdot (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot70(title, titll, titlr, titlb, rhon_half_save,
     &    wdoti1_ql, wdoti2_ql, wdoti3_ql, wdoti4_ql, wdoti5_ql,
     &    wdoti6_ql, wdote_ql, nnoderho_half, nrhomax)

*     ----------------------
*     write mchoi2 (67) file
*     ----------------------
      write(67, 309) nnoderho_half
      do n = 1, nnoderho_half
         write (67, 1312) n, rhon_half_save(n), redotj2avg_save(n),
     &      wdoti2avg_save(n), wdoti2_ql(n)
      end do


      title= 'Integrated quasilinear wdot'
      titll= 'P (Watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half_save,
     &   wdoti1_ql_int, wdoti2_ql_int, wdoti3_ql_int, wdoti4_ql_int,
     &   wdoti5_ql_int, wdoti6_ql_int, wdote_ql_int, nnoderho_half,
     &   nrhomax)

      title= 'Quasilinear wdoti1'
      titll= 'wdoti1 (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot10(title, titll, titlr, rhon_half_save, wdoti1_ql,
     &    nnoderho_half, nrhomax)

      title= 'Quasilinear wdoti2'
      titll= 'wdoti2 (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot10(title, titll, titlr, rhon_half_save, wdoti2_ql,
     &    nnoderho_half, nrhomax)


      title= 'Quasilinear wdoti3'
      titll= 'wdoti2 (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot10(title, titll, titlr, rhon_half_save, wdoti3_ql,
     &    nnoderho_half, nrhomax)



*     ----------------------------------
*     Plot on fine mesh from rf2x_setup2
*     ----------------------------------

      nnoderho2_half = nnoderho2 - 1

      do n = 1, nnoderho2_half
         rhon_half(n) = (rhon(n) + rhon(n+1)) / 2.0
c         write(6, 1312)n, rhon(n), rhon_half(n)
      end do


      title= 'Flux surface average redotj'
      titll= 'redotj (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half,
     &    redotj1avg, redotj2avg, redotj3avg, redotj4avg, redotj5avg,
     &    redotj6avg, redotjeavg, nnoderho2_half, nrhomax)

      nnoderho2_half = nnoderho2 - 1

      do n = 1, nnoderho2_half
         rhon_half(n) = (rhon(n) + rhon(n+1)) / 2.0
c         write(6, 1312)n, rhon(n), rhon_half(n)
      end do

      title= 'Flux surface average redotj'
      titll= 'redotj (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot70(title, titll, titlr, titlb, rhon_half,
     &   redotj1avg, redotj2avg, redotj3avg, redotj4avg, redotj5avg,
     &   redotj6avg, redotjeavg,  nnoderho2_half, nrhomax)

      title= 'Integrated Real(EdotJ)'
      titll= 'P(watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon,
     &    redotj1_int, redotj2_int, redotj3_int, redotj4_int,
     &    redotj5_int, redotj6_int, redotje_int, nnoderho2, nrhomax)


      title= 'Real(EdotJ) * dvol'
      titll= 'P(watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon,
     &    redotj1_dvol, redotj2_dvol, redotj3_dvol, redotj4_dvol,
     &    redotj5_dvol, redotj6_dvol, redotje_dvol, nnoderho2, nrhomax)

      title= 'Flux surface average wdot'
      titll= 'wdot (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half,
     &   wdoti1avg, wdoti2avg, wdoti3avg, wdoti4avg, wdoti5avg,
     &   wdoti6avg, wdoteavg, nnoderho2_half, nrhomax)

      write(72, 309) nnoderho
      do n = 1, nnoderho
         write(72,1312)n, rhon_half(n), wdoteavg(n), wdoti1avg(n),
     &              wdoti2avg(n), wdoti3avg(n), wdoti4avg(n),
     &              wdoti5avg(n), wdoti6avg(n)
      end do


      title= 'Integrated Wdot'
      titll= 'P (watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon,
     &    wdot1_int, wdot2_int, wdot3_int, wdot4_int, wdot5_int,
     &    wdot6_int, wdote_int, nnoderho2, nrhomax)

      title= 'Wdot * dvol'
      titll= 'P (watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon,
     &    wdoti1_dvol, wdoti2_dvol, wdoti3_dvol, wdoti4_dvol,
     &    wdoti5_dvol, wdoti6_dvol, wdote_dvol, nnoderho2, nrhomax)

      title= 'Flux surface average redotj'
      titll= 'redotj (watts/m3)'
      titlr='       '
      titlb= 'rho'

      call ezplot2p(title, titll, titlr, titlb, rhon_half, redotjeavg,
     &   redotjiavg, nnoderho2_half, nrhomax)

      write(15, *)
      write(15, *) 'Flux surface driven current'
      write(15, *)
      write(15, *) '        n      rho       J (A/m2)     I (A)'
      write(15, *)

      write(6, *)
      write(6, *) 'Flux surface driven current'
      write(6, *)
      write(6, *) '        n      rho       J (A/m2)     I (A)'
      write(6, *)

      do n = 1, nnoderho2_half
         write (15, 1312) n, rhon_half(n), xjprlavg(n), xjprl_int(n)
         write (6, 1312) n, rhon_half(n), xjprlavg(n), xjprl_int(n)
      end do

c      write(15, *)
c      write(15, *) 'Mod E in the midplane'
c      write(15, *)
c      write(15, *) '        i      R(i)       mod_E_mid      '
c      write(15, *)

c      write(6, *)
c      write(6, *) 'Mod E in the midplane'
c      write(6, *)
c      write(6, *) '        i      R(i)       mod_E_mid      '
c      write(6, *)

c      do i = 1, nnodex
c        write (15, 1312) i, capr(i), mod_E_mid(i)
c        write (6, 1312) i, capr(i), mod_E_mid(i)
c      end do

      title= 'Flux surface driven current'
      titll= 'xjprl (Amps/m2)'
      titlr= 'I (Amps)'
      titlb= 'rho'

      call ezplot2(title, titll, titlr, titlb, rhon_half, xjprlavg,
     &    xjprl_int, nnoderho2_half, nrhomax)


      if(xjprl_int(nnoderho2) .lt. 0.0) then

         titll= '-xjprl (Amps/m2)'
         titlr= '-I (Amps)'

         xjprlavg = - xjprlavg
         xjprl_int = - xjprl_int

         call ezplot2(title, titll, titlr, titlb, rhon_half, xjprlavg,
     &      xjprl_int, nnoderho2_half, nrhomax)

         titll= 'xjprl (MA/m2/MW)'

         xjprlavg = 0.0
         if(prfin .ne. 0.0) xjprlavg = xjprlavg  / prfin

         call ezplot1(title, titll, titlr, rhon_half, xjprlavg,
     &      nnoderho2_half, nrhomax)

         titll= 'I (kA)'
         xjprl_int  = xjprl_int / 1.0e+03
         call ezplot1(title, titll, titlr, rhon_half, xjprl_int,
     &      nnoderho2_half, nrhomax)

      end if

      if(xjprl_int(nnoderho2) .gt. 0.0) then

         titll= 'xjprl (Amps/m2)'
         titlr= 'I (Amps)'

         xjprlavg = xjprlavg
         xjprl_int = xjprl_int

         titll= 'xjprl (MAmps/m2/MW)'

         xjprlavg = 0.0
         if (prfin .ne. 0.0) xjprlavg = xjprlavg  / prfin

         call ezplot1(title, titll, titlr, rhon_half, xjprlavg,
     &      nnoderho2_half, nrhomax)

         titll= 'I (kA)'
         xjprl_int  = xjprl_int / 1.0e+03
         call ezplot1(title, titll, titlr, rhon_half, xjprl_int,
     &      nnoderho2_half, nrhomax)

      end if



      title= 'Volume element dvol'
      titll= 'dvol (m3)'
      titlr= '       '
      call ezplot1(title, titll, titlr, rhon_half, dvol,
     &     nnoderho2_half, nrhomax)

      title= 'Integrated volume'
      titll= 'volume (m3)'
      titlr= '       '
      call ezplot1(title, titll, titlr, rhon, volume,
     &     nnoderho2, nrhomax)

*     --------------------------------------
*     Reset rhon back to its original value
*     --------------------------------------
      do n = 1, nnoderho
         rhon(n) = rhon_save(n)
         rhon_half(n) = rhon_half_save(n)
      end do

*     ----------------------------------------------------------------
*     Read quasilinear diffusion coefficients from file out_cql3d.coef
*     ----------------------------------------------------------------


      if (ndisti2 .eq. 1 .or. ndisti1 .eq. 1) then

         if (ndisti1 .eq. 1)
     &      open(unit=42,file='out_cql3d.coef1', status='unknown',
     &                                             form='formatted')

         if (ndisti2 .eq. 1)
     &      open(unit=42,file='out_cql3d.coef2', status='unknown',
     &                                             form='formatted')

         read (42, 309) nuper
         read (42, 309) nupar
         read (42, 309) nnoderho

         allocate( UPERP(nuper) )
         allocate( UPARA(nupar) )
         allocate( f_cql_cart(nuper, nupar, nnoderho) )

         allocate( wperp1_cql(nnoderho) )
         allocate( wpar1_cql(nnoderho) )

         allocate( wperp2_cql(nnoderho) )
         allocate( wpar2_cql(nnoderho) )

         allocate( bqlavg_i1(nuper, nupar, nnoderho) )
         allocate( cqlavg_i1(nuper, nupar, nnoderho) )
         allocate( eqlavg_i1(nuper, nupar, nnoderho) )
         allocate( fqlavg_i1(nuper, nupar, nnoderho) )

         allocate( bqlavg_i1_2d(nupar, nuper) )
         allocate( cqlavg_i1_2d(nupar, nuper) )

         allocate( E_kick_2d (nupar, nnoderho) )
         allocate( f_cql_cart_2d (nupar, nuper) )


         read (42, 3310) vc_cgs
         read (42, 3310) UminPara, UmaxPara

         read (42, 3310) (rhon(n), n = 1, nnoderho)
         read (42, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
         read (42, 3310) (upara(i_upara), i_upara = 1, nupar)

         read (42, 3310) (((bqlavg_i1(i_uperp, i_upara, n),
     &     i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) (((cqlavg_i1(i_uperp, i_upara, n),
     &     i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) (((eqlavg_i1(i_uperp, i_upara, n),
     &     i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) (((fqlavg_i1(i_uperp, i_upara, n),
     &     i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) xmi

         close (42)

!        --------------------------------------------
!        read data for plotting f(u_perp, u_parallel)
!        --------------------------------------------

         if (ndisti1 .eq. 1) then
         open(unit=237,file='out237',status='old',form='formatted')

         read  (237, 309) n_u
         read  (237, 309) n_psi
         read  (237, 310) vc

         read  (237, 310) (u(i_u), i_u = 1, n_u)
         read  (237, 309) (n_theta_(i_psi), i_psi = 1, n_psi)
         read  (237, 310) ((theta(i_theta, i_psi),
     &          i_theta = 1, n_theta_(i_psi)), i_psi = 1, n_psi)

         read  (237, 310) (((f_cql(i_theta, i_u, i_psi),
     &     i_theta = 1, n_theta_(i_psi)), i_u = 1, n_u),
     &     i_psi = 1, n_psi)


         read (237, 309) nuper
         read (237, 309) nupar
         read (237, 309) nnoderho

         read (237, 310) (uperp(i_uperp), i_uperp = 1, nuper)
         read (237, 310)  (upara(i_upara), i_upara = 1, nupar)

         read (237, 310) (((f_cql_cart(i_uperp, i_upara, i_psi),
     &     i_uperp = 1, nuper), i_upara = 1, nupar),
     &     i_psi = 1, nnoderho)

         read (237, 310) (wperp1_cql(i_psi), i_psi = 1, nnoderho)
         read (237, 310) (wpar1_cql(i_psi),  i_psi = 1, nnoderho)

         close (237)
         end if



         if (ndisti2 .eq. 1) then
         open(unit=238,file='out238',status='old',form='formatted')

         read  (238, 309) n_u
         read  (238, 309) n_psi
         read  (238, 310) vc

         read  (238, 310) (u(i_u), i_u = 1, n_u)
         read  (238, 309) (n_theta_(i_psi), i_psi = 1, n_psi)
         read  (238, 310) ((theta(i_theta, i_psi),
     &          i_theta = 1, n_theta_(i_psi)), i_psi = 1, n_psi)

         read  (238, 310) (((f_cql(i_theta, i_u, i_psi),
     &     i_theta = 1, n_theta_(i_psi)), i_u = 1, n_u),
     &     i_psi = 1, n_psi)


         read (238, 309) nuper
         read (238, 309) nupar
         read (238, 309) nnoderho

         read (238, 310) (uperp(i_uperp), i_uperp = 1, nuper)
         read (238, 310)  (upara(i_upara), i_upara = 1, nupar)

         read (238, 310) (((f_cql_cart(i_uperp, i_upara, i_psi),
     &     i_uperp = 1, nuper), i_upara = 1, nupar),
     &     i_psi = 1, nnoderho)

         read (238, 310) (wperp2_cql(i_psi), i_psi = 1, nnoderho)
         read (238, 310) (wpar2_cql(i_psi),  i_psi = 1, nnoderho)

         close (238)
         end if

 4319    format(i10, 1p,1e16.8)

         title= 'Wperp1 and Wpar1 (keV)'
         titll= 'W (keV)'
         titlr='       '
         titlb = 'rho'

         call ezplot2q(title, titll, titlr, titlb, rhon_half,
     &      wperp1_cql, wpar1_cql, nnoderho_half, nrhomax)


         title= 'Wperp2 and Wpar2 (keV)'
         titll= 'W (keV)'
         titlr='       '
         titlb = 'rho'

         call ezplot2q(title, titll, titlr, titlb, rhon_half,
     &      wperp2_cql, wpar2_cql, nnoderho_half, nrhomax)


!       ---------------------------------
!       2D plots of f(u, theta = const)
!       ---------------------------------
        titll = 'log f(u)'
        titlr = '    '
        titx = 'u'

        do i_psi_index = 1, 6

           if(i_psi_index .eq. 1)title = 'log f_psi1(u, theta = const)'
           if(i_psi_index .eq. 2)title = 'log f_psi2(u, theta = const)'
           if(i_psi_index .eq. 3)title = 'log f_psi3(u, theta = const)'
           if(i_psi_index .eq. 4)title = 'log f_psi4(u, theta = const)'
           if(i_psi_index .eq. 5)title = 'log f_psi5(u, theta = const)'
           if(i_psi_index .eq. 6)title = 'log f_psi6(u, theta = const)'

           i_psi = i_psi_array(i_psi_index)

           i_theta1 = 1
           i_theta2 = int(1./6. * n_theta_(i_psi))
           i_theta3 = int(2./6. * n_theta_(i_psi))
           i_theta4 = int(3./6. * n_theta_(i_psi))
           i_theta5 = int(4./6. * n_theta_(i_psi))
           i_theta6 = int(5./6. * n_theta_(i_psi))
           i_theta7 = n_theta_(i_psi)

           do i_u = 1, n_u
              f_cql_1d_1(i_u) = -40.0
              f_cql_1d_2(i_u) = -40.0
              f_cql_1d_3(i_u) = -40.0
              f_cql_1d_4(i_u) = -40.0
              f_cql_1d_5(i_u) = -40.0
              f_cql_1d_6(i_u) = -40.0
              f_cql_1d_7(i_u) = -40.0

              if (f_cql(i_theta1, i_u, i_psi) .gt. 0.0)
     &        f_cql_1d_1(i_u) = alog10(f_cql(i_theta1, i_u, i_psi))
              if (f_cql(i_theta2, i_u, i_psi) .gt. 0.0)
     &         f_cql_1d_2(i_u) = alog10(f_cql(i_theta2, i_u, i_psi))
              if (f_cql(i_theta3, i_u, i_psi) .gt. 0.0)
     &         f_cql_1d_3(i_u) = alog10(f_cql(i_theta3, i_u, i_psi))
              if (f_cql(i_theta4, i_u, i_psi) .gt. 0.0)
     &         f_cql_1d_4(i_u) = alog10(f_cql(i_theta4, i_u, i_psi))
              if (f_cql(i_theta5, i_u, i_psi) .gt. 0.0)
     &         f_cql_1d_5(i_u) = alog10(f_cql(i_theta5, i_u, i_psi))
              if (f_cql(i_theta6, i_u, i_psi) .gt. 0.0)
     &         f_cql_1d_6(i_u) = alog10(f_cql(i_theta6, i_u, i_psi))
              if (f_cql(i_theta7, i_u, i_psi) .gt. 0.0)
     &         f_cql_1d_7(i_u) = alog10(f_cql(i_theta7, i_u, i_psi))
           end do


           call ezlog1_f(title, titll, titlr, titx, u, f_cql_1d_1,
     &        f_cql_1d_2, f_cql_1d_3, f_cql_1d_4, f_cql_1d_5,
     &        f_cql_1d_6, f_cql_1d_7,
     &        n_u, n_u_max)

        end do




!        ------------------------------------------
!        2D plots of f_cql_cart(u_perp, u_parallel)
!        ------------------------------------------
         titx = 'u_parallel'
         tity = 'u_perp'

         write(6, *)"i_psi1 = ", i_psi1, "   rho1 = ", rhon(i_psi1)
         write(6, *)"i_psi2 = ", i_psi2, "   rho2 = ", rhon(i_psi2)
         write(6, *)"i_psi3 = ", i_psi3, "   rho3 = ", rhon(i_psi3)
         write(6, *)"i_psi4 = ", i_psi4, "   rho4 = ", rhon(i_psi4)
         write(6, *)"i_psi5 = ", i_psi5, "   rho5 = ", rhon(i_psi5)
         write(6, *)"i_psi6 = ", i_psi6, "   rho6 = ", rhon(i_psi6)


         write(15, *)"i_psi1 = ", i_psi1, "   rho1 = ", rhon(i_psi1)
         write(15, *)"i_psi2 = ", i_psi2, "   rho2 = ", rhon(i_psi2)
         write(15, *)"i_psi3 = ", i_psi3, "   rho3 = ", rhon(i_psi3)
         write(15, *)"i_psi4 = ", i_psi4, "   rho4 = ", rhon(i_psi4)
         write(15, *)"i_psi5 = ", i_psi5, "   rho5 = ", rhon(i_psi5)
         write(15, *)"i_psi6 = ", i_psi6, "   rho6 = ", rhon(i_psi6)

         numb = 20

         i_psi = i_psi1
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) =
     &                     f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi1(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         i_psi = i_psi2
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) =
     &                     f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi2(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         i_psi = i_psi3
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) =
     &                     f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi3(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         i_psi = i_psi4
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) =
     &                     f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi4(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         i_psi = i_psi5
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) =
     &                     f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi5(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)


         i_psi = i_psi6
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) =
     &                     f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi6(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)



!        --------------------------------------
!        2D plots of bqlavg(u_perp, u_parallel)
!        --------------------------------------

         i_psi = i_psi1
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) =
     &                     bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi1(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)


!        ---------------------------------------------
!        Write Bql_avg_2D.vtk file "structured points"
!        ---------------------------------------------
         number_points = nupar * nuper
         dx = upara(2) - upara(1)
         dy = uperp(2) - uperp(1)
         write(142, 2840)
         write(142, 4850)
 4850    format('Quasilinear diffusion coefficients')
         write(142, 2846)
         write(142, 2847)
         write(142, 2841) nupar,  nuper, 1
         write(142, 2842) upara(1), uperp(1), 0
         write(142, 2843) dx, dy, 1
         write(142, 2844) number_points

         write(142, 4851)
 4851    format('SCALARS bqlavg_psi1 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)



         i_psi = i_psi2
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) =
     &                     bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi2(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(142, 4852)
 4852    format('SCALARS bqlavg_psi2 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)






         i_psi = i_psi3
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) =
     &                     bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi3(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(142, 4853)
 4853    format('SCALARS bqlavg_psi3 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)





         i_psi = i_psi4
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) =
     &                     bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi4(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(142, 4854)
 4854    format('SCALARS bqlavg_psi4 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)






         i_psi = i_psi5
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) =
     &                     bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi5(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(142, 4855)
 4855    format('SCALARS bqlavg_psi5 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)




         i_psi = i_psi6
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) =
     &                     bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi6(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(142, 4856)
 4856    format('SCALARS bqlavg_psi6 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)





























!        --------------------------------------
!        2D plots of cqlavg(u_perp, u_parallel)
!        --------------------------------------

         i_psi = i_psi1
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlavg_i1_2d(i_upara, i_uperp) =
     &                     cqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do

         title = 'cqlavg_psi1(u_perp, u_parallel)'
         call ezconc(upara, uperp, cqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)


!        ---------------------------------------------
!        Write Cql_avg_2D.vtk file "structured points"
!        ---------------------------------------------
         number_points = nupar * nuper
         dx = upara(2) - upara(1)
         dy = uperp(2) - uperp(1)
         write(242, 2840)
         write(242, 4850)

         write(242, 2846)
         write(242, 2847)
         write(242, 2841) nupar,  nuper, 1
         write(242, 2842) upara(1), uperp(1), 0
         write(242, 2843) dx, dy, 1
         write(242, 2844) number_points

         write(242, 4951)
 4951    format('SCALARS cqlavg_psi1 float 1')
         write(242, 2849)
         write(242, 3411) ((cqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)



         i_psi = i_psi2
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlavg_i1_2d(i_upara, i_uperp) =
     &                     cqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do

         title = 'cqlavg_psi2(u_perp, u_parallel)'
         call ezconc(upara, uperp, cqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(242, 4952)
 4952    format('SCALARS cqlavg_psi2 float 1')
         write(242, 2849)
         write(242, 3411) ((cqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)






         i_psi = i_psi3
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlavg_i1_2d(i_upara, i_uperp) =
     &                     cqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do

         title = 'cqlavg_psi3(u_perp, u_parallel)'
         call ezconc(upara, uperp, cqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(242, 4953)
 4953    format('SCALARS cqlavg_psi3 float 1')
         write(242, 2849)
         write(242, 3411) ((cqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)





         i_psi = i_psi4
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlavg_i1_2d(i_upara, i_uperp) =
     &                     cqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do

         title = 'cqlavg_psi4(u_perp, u_parallel)'
         call ezconc(upara, uperp, cqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(242, 4954)
 4954    format('SCALARS cqlavg_psi4 float 1')
         write(242, 2849)
         write(242, 3411) ((cqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)






         i_psi = i_psi5
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlavg_i1_2d(i_upara, i_uperp) =
     &                     cqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do

         title = 'cqlavg_psi5(u_perp, u_parallel)'
         call ezconc(upara, uperp, cqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(242, 4955)
 4955    format('SCALARS cqlavg_psi5 float 1')
         write(242, 2849)
         write(242, 3411) ((cqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)




         i_psi = i_psi6
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlavg_i1_2d(i_upara, i_uperp) =
     &                     cqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do

         title = 'cqlavg_psi6(u_perp, u_parallel)'
         call ezconc(upara, uperp, cqlavg_i1_2d, ff, nupar, nuper, numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag,1.0)

         write(242, 4956)
 4956    format('SCALARS cqlavg_psi6 float 1')
         write(242, 2849)
         write(242, 3411) ((cqlavg_i1_2d(i,j), i = 1, nupar),
     &                                         j = 1, nuper)




























!        --------------------------------------
!        2D plots energy kicks (eV)
!        --------------------------------------

         write(6,  *) 'energy kicks (eV)'
         write(15, *) 'energy kicks (eV)'

         title = 'Energy kick (eV)'
         titz='Energy kick (eV)'
         numb = 15

         i_psi = i_psi1
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
            vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
               vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
               E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) *
     &            sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs)
     &            * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)



         i_psi = i_psi2
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
            vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
               vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
               E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) *
     &            sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs)
     &            * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)




         i_psi = i_psi3
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
            vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
               vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
               E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) *
     &            sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs)
     &            * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)




         i_psi = i_psi4
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
            vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
               vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
               E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) *
     &            sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs)
     &            * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,numb,
     &      NUPAR, NUPER, nlevmax, title, titx, tity, iflag)



         i_psi = i_psi5
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
            vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
               vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
               E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) *
     &            sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs)
     &            * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,
     &      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)




         i_psi = i_psi6
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
            vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
               vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
               E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) *
     &            sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs)
     &            * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,
     &      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)



!        --------------------------------------
!        2D plots energy kicks for 1keV ion (eV)
!        --------------------------------------
         write(6, *)
         write(6, *) 'energy kicks for 1 keV ion (eV)'

         write(15, *)
         write(15, *) 'energy kicks for 1 keV ion (eV)'

         title = 'Energy kick for 1 keV ion (eV)'
         titx = 'u_parallel'
         tity = 'rho'
         titz='Energy kick (eV)'
         numb = 15


         E_eV = 1000.
         vperp_mks = sqrt(2.0 * E_eV * q / xmi)
         vperp_cgs = vperp_mks * 100.
         uperp_1kev = vperp_cgs / vc_cgs
         duperp = uperp(3) - uperp(2)
         i_uperp = int(uperp_1kev / duperp) + 2

         uperp_1kev = (i_uperp - 1) * duperp
         vperp_cgs = uperp_1kev * vc_cgs
         vperp_mks = vperp_cgs / 100.

         E_eV = 0.5 * xmi * (vperp_mks**2) / q

         write (6, *) 'i_uperp =', i_uperp
         write (6, *) 'E = ', E_eV, ' eV'

         write (15, *) 'i_uperp =', i_uperp
         write (15, *) 'E = ', E_eV, ' eV'

         do i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
            vpara_mks = vpara_cgs * 1.0e-02
            do i_psi = 1, nnoderho
               vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
               E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               E_kick_2d(i_upara, i_psi) = (xmi / q) *
     &            sqrt(2. * bqlavg_i1(i_uperp, i_upara, i_psi)
     &            / vpara_cgs) * 1.0e-04
            end do
         end do

         call ezconcx(upara, rhon, E_kick_2d, ff, nupar, nnoderho, numb,
     &      NUPAR, nnoderho, nlevmax, title, titx, tity, iflag)
         call ezcon3dv(upara, rhon, E_kick_2d, ff, nupar, nnoderho,
     &      numb, NUPAR, nnoderho, nlevmax, title, titx, tity, titz)


!        ---------------------------------------------
!        Write E_kicks_2D.vtk file "structured points"
!        ---------------------------------------------
         number_points = nupar * nnoderho
         dx = upara(2) - upara(1)
         dy = rhon(2) - rhon(1)

         write(141, 2840)
         write(141, 2852)
 2852    format('Energy kicks for 1 keV')
         write(141, 2846)
         write(141, 2847)
         write(141, 2841) nupar,  nnoderho, 1
         write(141, 2842) upara(1), rhon(1), 0
         write(141, 2843) dx, dy, 1
         write(141, 2844) number_points

         write(141, 3848)
 3848    format('SCALARS E_kick_2d float 1')
         write(141, 2849)
         write(141, 3411) ((E_kick_2d(i,j), i = 1, nupar),
     &                                      j = 1, nnoderho)


      end if

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(bxwave(i,j))
            fimag(i,j) = aimag(bxwave(i,j))
         end do
      end do

      titx = 'R (m)'
      tity = 'Z (m)'

      title = 'Real Bx wave'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag, scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(bzwave(i,j))
            fimag(i,j) = aimag(bzwave(i,j))
         end do
      end do


      title = 'Real Bz wave'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = dxxuyy(i,j)
            if(rho(i,j) .gt. 1.0) freal(i,j) = 0.0
         end do
      end do


      title = 'dxxuyy'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)





      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = dyyuzz(i,j)
            if(rho(i,j) .gt. 1.0) freal(i,j) = 0.0
         end do
      end do


      title = 'dyyuzz'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = gradprlb(i,j)
            if(rho(i,j) .gt. 1.0) freal(i,j) = 0.0
         end do
      end do


      title = 'gradprlb'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)




      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ntilda_e(i,j))
            fimag(i,j) = aimag(ntilda_e(i,j))
         end do
      end do


      title = 'Real ntilda_e'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


c
c--plot ntilda_e in midplane
c
      do i = 1, nnodex
         fmidre(i) = real(ntilda_e(i, nnodey / 3 - 1))
         fmidim(i) = aimag(ntilda_e(i, nnodey / 3 - 1))
      end do

      title = 'ntilda_e'
      titll = 'Re ntilda_e (m^{-3})'
      titlr = 'Im ntilda_e (m^{-3})'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)



c
c--plot Bx_wave in midplane
c
      do i = 1, nnodex
         fmidre(i) = real(bxwave(i,jmid))
         fmidim(i) = aimag(bxwave(i,jmid))
      end do

      title = 'Bx wave'
      titll = 'Re Bx wave (T)'
      titlr = 'Im Bx wave (T)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)


c
c--plot Bz_wave in midplane
c
      do i = 1, nnodex
         fmidre(i) = real(bzwave(i,jmid))
         fmidim(i) = aimag(bzwave(i,jmid))
      end do

      title = 'Bz wave'
      titll = 'Re Bz wave (T)'
      titlr = 'Im Bz wave (T)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)


*     -------------------
*     plot fpsi in 2-D
*     -------------------

      title = 'Radial force (fpsi0)'
      call ezconc(capr, y, fpsi0, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = 'Poloidal force (ftheta0)'
      call ezconc(capr, y, ftheta0, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      title = 'pressure surfaces'
      call ezconc(capr, y, pressi, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      numb = 20

      title = 'dl/B surfaces'
      call ezconc(capr, y, dldb_tot12, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)


      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(eplus_flux_plot(i,j))
            fimag(i,j) = aimag(eplus_flux_plot(i,j))
         end do
      end do


      title = 'E_{+} flux plot'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)




      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(eminus_flux_plot(i,j))
            fimag(i,j) = aimag(eminus_flux_plot(i,j))
         end do
      end do


      title = 'E_{-} flux plot'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)





      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(xkperp_flux_plot(i,j))
            fimag(i,j) = aimag(xkperp_flux_plot(i,j))
         end do
      end do


      title = 'xkperp flux plot'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &   nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &     nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)



*     --------------------------------
*     plot eplus_flux_plot(i,j) in midplane
*     --------------------------------
      do i = 1, nnodex
         fmidre(i) = real(eplus_flux_plot(i,jmid))
         fmidim(i) = aimag(eplus_flux_plot(i,jmid))
      end do

      title = 'E_{+} flux plot'
      titll = 'Re Eplus (V/m)'
      titlr = 'Im Eplus (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)


*     --------------------------------
*     plot eminus_flux_plot(i,j) in midplane
*     --------------------------------
      do i = 1, nnodex
         fmidre(i) = real(eminus_flux_plot(i,jmid))
         fmidim(i) = aimag(eminus_flux_plot(i,jmid))
      end do

      title = 'E_{--} flux plot'
      titll = 'Re Eminus (V/m)'
      titlr = 'Im Eminus (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)



*     --------------------------------------
*     plot xkperp_flux_plot(i,j) in midplane
*     --------------------------------------
      do i = 1, nnodex
         fmidre(i) = real(xkperp_flux_plot(i,jmid))
         fmidim(i) = aimag(xkperp_flux_plot(i,jmid))
      end do

      title = 'xkperp flux plot'
      titll = 'Re xkperp s(1/m)'
      titlr = 'Im xkperp (1/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim,
     &   nnodex, nxmx)



      do n = 1, nnoderho
c         bmod_plot(n) = bmod_flux(n, 1) / bmod_flux(1, 1)
          bmod_plot(n) = bmod_flux(n, 1)
      end do

      title= 'Mod B (theta = 0)'
      titll= 'bmod flux (T)'
      titlr='       '
      call ezplot0(title, titll, titlr, rhon, bmod_plot,
     &    nnoderho, nrhomax)





*     ---------------------
*     write mchoi (65) file
*     ---------------------
      write(65, 309) nnoderho, mnodetheta
      do n = 1, nnoderho
         do m = 1, mnodetheta
            write (65, 1313) n, m, rhon(n), thetam(m),
     &         capr_flux(n, m), capz_flux(n, m), eplus_flux(n,m),
     &                                          eminus_flux(n,m),
     &                                          xkperp_flux(n,m)

            freal(n,m) = real(eplus_flux(n,m))
            fimag(n,m) = aimag(eplus_flux(n,m))
         end do
      end do

*     ---------------------
*     write mchoi2 (67) file
*     ---------------------
c      write(67, 309) nnoderho_half
c      do n = 1, nnoderho_half
c         write (67, 1312) n, rhon_half_save(n), redotj2avg(n),
c     &      wdoti2avg(n), wdoti2_ql(n)
c      end do


      title = 'E_{+} flux'
      titx   = 'rho'
      tity = 'theta'
      call ezconc(rhon, thetam, freal, ff, nnoderho, mnodetheta, numb,
     &     nrhomax, nthetamax, nlevmax, title, titx, tity,
     &     iflag, scalex)


*     ---------------------
*     write mchoi3 (69) file
*     ---------------------

      write(69, 309) nnodex, nnodey
      do i = 1, nnodex
         do j = 1, nnodey
            write (69, 1313) i, j, capr(i), y(j), wdoti2(i,j)
         end do
      end do



      do n = 1, nnoderho
         do m = 1, mnodetheta
            freal(n,m) = real(eminus_flux(n,m))
            fimag(n,m) = aimag(eminus_flux(n,m))
         end do
      end do

      title = 'E_{--} flux'
      titx   = 'rho'
      tity = 'theta'
      call ezconc(rhon, thetam, freal, ff, nnoderho, mnodetheta, numb,
     &     nrhomax, nthetamax, nlevmax, title, titx, tity,
     &     iflag, scalex)



      write(140, 3856)
 3856 format('SCALARS fpsi0 float 1')
      write(140, 2849)
      write (140, 3411) ((fpsi0(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3858)
 3858 format('SCALARS ftheta0 float 1')
      write(140, 2849)
      write (140, 3411) ((ftheta0(i,j), i = 1, nnodex),  j = 1, nnodey)


      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ntilda_e(i,j)) / 1.0e+19
            fimag(i,j) = aimag(ntilda_e(i,j)) / 1.0e+19
         end do
      end do

      write(140, 3859)
 3859 format('SCALARS ntilda_e_real float 1')
      write(140, 2849)
      write (140, 3411) ((freal(i,j), i = 1, nnodex), j = 1, nnodey)


      write(140, 3960)
 3960 format('SCALARS divq float 1')
      write(140, 2849)
      write (140, 3411) ((divq(i,j), i = 1, nnodex), j = 1, nnodey)

      write(140, 3961)
 3961 format('SCALARS reomg1 float 1')
      write(140, 2849)
      write (140, 3411) ((reomg1a(i,j), i = 1, nnodex), j = 1, nnodey)

      write(140, 3962)
 3962 format('SCALARS reomg2 float 1')
      write(140, 2849)
      write (140, 3411) ((reomg2a(i,j), i = 1, nnodex), j = 1, nnodey)

      write(140, 3963)
 3963 format('SCALARS reomglh float 1')
      write(140, 2849)
      write (140, 3411) ((reomglha(i,j), i = 1, nnodex), j = 1, nnodey)

      write(140, 3964)
 3964 format('SCALARS Re_kperp2_fast float 1')
      write(140, 2849)
      write (140, 3411) ((real(xkperp2_fast(i,j)), i = 1, nnodex),
     &   j = 1, nnodey)

      write(140, 3984)
 3984 format('SCALARS Im_kperp2_fast float 1')
      write(140, 2849)
      write (140, 3411) ((aimag(xkperp2_fast(i,j)), i = 1, nnodex),
     &   j = 1, nnodey)

      write(140, 3965)
 3965 format('SCALARS Re_kperp2_slow float 1')
      write(140, 2849)
      write (140, 3411) ((real(xkperp2_slow(i,j)), i = 1, nnodex),
     &   j = 1, nnodey)

      write(140, 3985)
 3985 format('SCALARS Im_kperp2_slow float 1')
      write(140, 2849)
      write (140, 3411) ((aimag(xkperp2_slow(i,j)), i = 1, nnodex),
     &   j = 1, nnodey)

      write(140, 3968)
 3968 format('SCALARS xkprl_a float 1')
      write(140, 2849)
      write (140, 3411) ((xkprl_a(i,j), i = 1, nnodex), j = 1, nnodey)

      write(140, 3969)
 3969 format('SCALARS Re_P_a float 1')
      write(140, 2849)
      write (140, 3411) ((real(P_a(i,j)), i = 1, nnodex), j = 1, nnodey)

      write(140, 3989)
 3989 format('SCALARS Re_P_a float 1')
      write(140, 2849)
      write (140, 3411) ((real(P_a(i,j)), i = 1, nnodex), j = 1, nnodey)



      freal = xn / 1.0e+19
      write(140, 2168)
 2168 format('SCALARS xn float 1')
      write(140, 2849)
      write (140, 3411) ((freal(i,j), i = 1, nnodex), j = 1, nnodey)


      freal = xkte / q
      write(140, 2169)
 2169 format('SCALARS xkte float 1')
      write(140, 2849)
      write (140, 3411) ((freal(i,j), i = 1, nnodex), j = 1, nnodey)

!      write(140, 2850)
! 2850 format('VECTORS Sp_vector float 1')
!      write(140, 2849)
!      write(140, 3412) ((spx(i,j), spy(i,j), spz(i,j),
!     &                                    i = 1, nnodex), j = 1, nnodey)

 3412 format(3f17.4)

!     --------------------------------------
!     Write 2D.vtk file "structured points"
!     --------------------------------------
      write(245, 2840)
      write(245, 2855)
 2855 format('2D Poynting vector field')

      write(245, 2846)
      write(245, 2847)
      write(245, 2841) nnodex,  nnodey, 1
      write(245, 2842) capr(1), y(1), 0
      dx = capr(2) - capr(1)
      dy = y(2) - y(1)
      write(245, 2843) dx, dy, 1
      write(245, 2844) number_points


      write(245, 2858)
 2858 format('SCALARS Spx float 1')
      write(245, 2849)
      write (245, 3411) ((Spx(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(245, 3873)
 3873 format('SCALARS Spy float 1')
      write(245, 2849)
      write (245, 3411) ((Spy(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(245, 7853)
 7853 format('SCALARS Spz float 1')
      write(245, 2849)
      write (245, 3411) ((Spz(i,j), i = 1, nnodex),  j = 1, nnodey)

       write(245, 2870)
 2870 format('VECTORS Sp_vector float 1')
      write(245, 2849)
      write(245, 3412) ((spx(i,j), spy(i,j), spz(i,j),
     &                                    i = 1, nnodex), j = 1, nnodey)

      close (245)


 5000 continue

c Close the graphics device.

      call pgclos

c Open new graphics device


#ifdef GIZA      
      IER = PGBEG(0, 'movie.pdf', 1, 1)
#else
      IER = PGBEG(0, 'movie.ps/vcps', 1, 1)
#endif
      
      IF (IER.NE.1) STOP

      call PGSCH (1.5)


      titx = 'R (m)'
      tity = 'Z (m)'


      pi = 3.14159
      period = 2.0 * pi / omgrf
      dt = period / 10.
      time = 0.0

      do it = 1, 50
         time = time + dt

         do i = 1, nnodex !?do concurrent or slicing indices?
            do j = 1, nnodey
               freal(i,j) = real(eb(i,j) * exp(-zi * omgrf * time))
            end do
         end do

         title = 'Real Eb(t)'

c        black and white movie:
         call ezconc(capr, y, freal, ff, nnodex, nnodey, numb,
     &      nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)


c        color movie:
c         call pseudo(capr, y, freal, ff, nnodex, nnodey, numb,
c     &      nxmx, nymx, nlevmax, title, titx, tity, iflag,scalex)

         if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
     &        nnodey, numb,    nxmx, nymx, nlevmax, title, titx, tity,
     &     scalex)

      end do

      call pgclos

      close (140)


      close (38)

      close (141)
      close (142)
      close (242)
      close (143)
      close (144)

      close (51)
c      close (61)
      close (66)
      close (68)
      close (57)
      close (62)
      close (65)
      close (67)
      close (64)

      if (ndisti2 .eq. 1 .or. ndisti1 .eq. 1) then
         deallocate( UPERP )
         deallocate( UPARA )

         deallocate( bqlavg_i1 )
         deallocate( cqlavg_i1 )
         deallocate( eqlavg_i1 )
         deallocate( fqlavg_i1 )

         deallocate( bqlavg_i1_2d )
         deallocate( cqlavg_i1_2d )

         deallocate( E_kick_2d  )
         deallocate( f_cql_cart)
         deallocate( f_cql_cart_2d)

         deallocate (wperp1_cql)
         deallocate (Wpar1_cql)

         deallocate (wperp2_cql)
         deallocate (Wpar2_cql)
      end if

      deallocate (xkperp_cold, acold, xkperp_cold2, bcold, ccold,
     &   ex, ey, ez, bxwave, bywave, bzwave,
     &   ealpha, ebeta, eb, eplus, eminus, eplus_flux_plot,
     &   eminus_flux_plot, xkperp_flux_plot, ntilda_e )

      deallocate (xkperp2_slow, xkperp2_fast)
      deallocate (xkprl_a, P_a)

      deallocate (xjpxe, xjpye, xjpze)

      deallocate ( exkmod,
     &     eykmod, ezkmod, exklog, eyklog, ezklog )

      deallocate ( eplus_flux, eminus_flux, xkperp_flux)

      deallocate ( exk, eyk, ezk )

      deallocate (xn, xkti, xkte, rho, theta, xjz,
     &   psi, xjy, bmod, xjx, xiota, qsafety, ipsi, btau, bzeta,
     &   freal, fimag, fmod, mod_Eplus, mod_Eminus, mod_Eb,
     &   mod_Ealpha, mod_Ebeta, Mod_E,
     &   dldb_tot12, reex_dx, reey_dx, reez_dx, ximex_dx,
     &   ximey_dx, ximez_dx, spx, spy, spz, reomg1a, reomg2a,
     &   reomg3a, capr_bpol, pressi, redotj1, redotj2, redotj3,
     &   redotje, redotjt, redotjs, wdoti1, wdoti2, wdoti3,
     &   wdoti4, wdoti5, wdoti6, wdote, wdott, fz0e, fz0i1,
     &   fz0i2, fz0i3, fz0, fz0i4, fz0i5, fz0i6, gpsi,  kpsi,
     &   omgexb, uzeta, utheta, fpsi0, ftheta0, muhat, nu_star,
     &   redotj4, redotj5, redotj6, divq, reomglha )

c      deallocate (dxuxx, dxuxy, dxuxz,
c     &   dxuyx, dxuyy, dxuyz, dxuzx, dxuzy, dxuzz,
c     &   dyuxx, dyuxy, dyuxz, dyuyx, dyuyy, dyuyz,
c     &   dyuzx, dyuzy, dyuzz, dyyuxx, dyyuxy, dyyuxz,
c     &   dyyuyx, dyyuyy, dyyuyz, dyyuzx, dyyuzy, dyyuzz,
c     &   dxyuxx, dxyuxy, dxyuxz, dxyuyx, dxyuyy, dxyuyz,
c     &   dxyuzx, dxyuzy, dxyuzz, dxxuxx, dxxuxy, dxxuxz,
c     &   dxxuyx, dxxuyy, dxxuyz, dxxuzx, dxxuzy, dxxuzz)

      deallocate (dxxuyy, dyyuzz)

      deallocate (gradprlb)
      deallocate (capd)
      deallocate (capd_plot)

      deallocate (xkb)
      deallocate (xkb_plot)

  310 format(1p,6e12.4)
  309 format(10i10)
  311 format(1p,10e12.4)
 3310 format(1p,6e18.10)
 1312 format(i10,1p,9e12.4)
 1313 format(2i10,1p,10e12.4)
 9310 format(1p,7e12.4)
   10 format(i10,1p,4e10.3,i10,1p,e10.3)

      return
      end
c
c********************************************************************
c

c Subroutine a1mnmx
c----------------------------------------------------------------------------!
c Minimum and the maximum elements of a 1-d array.
c----------------------------------------------------------------------------!

      subroutine a1mnmx(f, nxmax, nx, fmin, fmax)

      implicit none

      integer:: nx, i, nxmax
      real:: f(nxmax), fmin, fmax

      fmax=f(1)
      fmin=fmax
         do i=1, nx
            fmax=max(fmax, f(i))
            fmin=min(fmin, f(i))
         end do

      return
 2201 format(2i5,1p,8e12.4)
      end

c
c********************************************************************
c
c Subroutine a2mnmx
c----------------------------------------------------------------------------!
c Minimum and the maximum elements of a 1-d array between x1 and x2
c----------------------------------------------------------------------------!

      subroutine a2mnmx(f, nxmax, nx, x, x1, x2, fmin, fmax)

      implicit none

      integer:: nx, i, nxmax
      real:: f(nxmax), fmin, fmax, x1, x2
      real:: x(nxmax)

      fmax=0.0
      fmin=fmax
      do i = 1, nx
         if(x(i).gt.x1.and.x(i).lt.x2)then
            fmax=max(fmax, f(i))
            fmin=min(fmin, f(i))
         endif
      end do

      return

      end

c
c********************************************************************
c
c----------------------------------------------------------------------------!
c Subroutine a2dmnmx_r4
c----------------------------------------------------------------------------!
c Minimum and the maximum elements of a 2-d array.
c----------------------------------------------------------------------------!

      subroutine a2dmnmx_r4(f, nxmax, nymax, nx, ny, fmin, fmax)

      implicit none

      integer:: nx, ny, i, j, nxmax, nymax
      real:: f(nxmax, nymax), fmin, fmax

      fmax=f(2, 1)
      fmin=fmax

      do j=1, ny
c         do 23002 i=1, nx
         do i=2, nx-1
            fmax=max(fmax, f(i, j))
            fmin=min(fmin, f(i, j))
         end do
      end do

      return
 2201 format(2i5,1p,8e12.4)
      end

c
c********************************************************************
c
      subroutine ezconc(r, theta, f, flevel, nr, nth, nlevel,
     &   nrmax, nthmax, nlevmax, title, titx, tity, iflag, scale)

      implicit none

!     real, intent(in), optional:: scale
      real:: scale
      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &     nndm1, norange

      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &    ncolelec, iflag,rcol

      real:: fmin,ymax,df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, theta, r, flevel, f, xmin, ymin, xmax

      real:: tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character(8):: xopt,yopt
      character(32):: title
      character(32):: titx
      character(32):: tity

      real:: scalex

      scalex=scale
!      if (present(scale)) scalex=scale
      
      
      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      dx = (r(2) - r(1))*scalex
      dy = theta(2) - theta(1)

      tr(1) = r(1)*scalex - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      call a1mnmx(r*scalex, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
      !This doesn't work because work size is greater than used size
!      xmin = MINVAL(r)
!      xmax = MAXVAL(r)
!      ymin = MINVAL(theta)
!      ymax = MAXVAL(theta)

c--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)
!      fmin = MINVAL(f) !(2:nrmax-1,:))
!      fmax = MAXVAL(f) !(2:nrmax-1,:))

#ifdef DEBUG
      write(6,*) 'ezconc title: ', title,
     & ymin,ymax,xmin,xmax,fmin,fmax
#endif
      iflag = 0
      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
         write(6, *)"title = ", title
         write(6, *)"fmax = ", fmax
         write(6, *)"fmin = ", fmin
         iflag = 1
         return
      end if

      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
c        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel==1) then
         flevel(1)=1.0e-03
      else if (nlevel==4 .or. nlevel==15) then
         flevel(1:nlevel) = (/(i*1.0, i=1,nlevel)/)
      end if

c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) exit
         nlevlt = nlevlt + 1
      end do

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps

      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

c--Advance graphics frame and get ready to plot

c     ! call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 1, 0)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then
         call pgsci(nred)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &        nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then
!         CALL PGSFS(1)
!         DO I=1, Nlevel-1
!            RCOL = 0.5+0.5*REAL(I-1)/REAL(NLEVEL-1)
!            CALL PGSCR(I+10, RCOL, RCOL, RCOL)
!            CALL PGSCI(I+10)
!            CALL PGCONF(f, nrmax, nthmax, 1, nr, 1, nth,
!     &        flevel(I),flevel(I+1),TR)
!         END DO
         call pgsci(nblue)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &       nlevgt, tr)
      endif

      call pgsci(nblack)
      call pglab(titx, tity, title)

  310 format(1p,6e12.4)
  312 format(i10, 1p,6e12.4)

      return
      end
c
c*********************************************************************
c
!JCW only used by LH plot
      subroutine ezconc1(r, theta, f, flevel, nr, nth, nlevel,
     &   nrmax, nthmax, nlevmax, title, titx, tity, iflag,scale)

      implicit none

      real:: scale      
      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &   nndm1, norange

      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &    ncolelec, iflag

      real:: fmin,ymax,df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, theta, r, flevel, f, xmin, ymin, xmax

      real:: tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character(8):: xopt,yopt
      character(32):: title
      character(32):: titx
      character(32):: tity

      real:: scalex
      scalex=scale
            
      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      dx = (r(2) - r(1))*scalex
      dy = theta(2) - theta(1)

      tr(1) = r(1)*scalex - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      call a1mnmx(r*scalex, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
!      xmin = MINVAL(r)
!      xmax = MAXVAL(r)
!      ymin = MINVAL(theta)
!      ymax = MAXVAL(theta)

c--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

c      write(6, *)"fmax = ", fmax, "   fmin = ", fmin
c      write(15,*)"fmax = ", fmax, "   fmin = ", fmin

      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
c         write(6, *)"fmax = ", fmax
c         write(6, *)"fmin = ", fmin
         iflag = 1
         return
      end if

      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
c        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if (nlevel .eq. 1) flevel(1) = 1.0

      if(nlevel.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif


      if(nlevel.eq.15)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
         flevel(5)=5.0
         flevel(6)=6.0
         flevel(7)=7.0
         flevel(8)=8.0
         flevel(9)=9.0
         flevel(10)=10.0
         flevel(11)=11.0
         flevel(12)=12.0
         flevel(13)=13.0
         flevel(14)=14.0
         flevel(15)=15.0
      endif

c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) exit
         nlevlt = nlevlt + 1
      end do

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps

      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

c--Advance graphics frame and get ready to plot

      call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 1, 0)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

         call pgsci(nyellow)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &       nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then

         call pgsci(nblue)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &       nlevgt, tr)
      endif

      call pgsci(nblack)
      call pglab(titx, tity, title)

  310 format(1p,6e12.4)
  312 format(i10, 1p,6e12.4)


      return
      end
c
c*********************************************************************
c



      subroutine ezconc2(r, theta, f, flevel, nr, nth, nlevel,
     &   nrmax, nthmax, nlevmax, title, titx, tity, iflag, dfquotient)

      implicit none

      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &   nndm1, norange

      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &    ncolelec, iflag

      real:: dfquotient

      real:: fmin,ymax,df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, theta, r, flevel, f, xmin, ymin, xmax

      real:: tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character(8):: xopt,yopt
      character(32):: title
      character(32):: titx
      character(32):: tity

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      dx = r(2) - r(1)
      dy = theta(2) - theta(1)

      tr(1) = r(1) - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
!      xmin = MINVAL(r)
!      xmax = MAXVAL(r)
!      ymin = MINVAL(theta)
!      ymax = MAXVAL(theta)

c--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

#ifdef DEBUG
      write(*,*) 'title: ', title
      write(6, *)"fmax = ", fmax, "   fmin = ", fmin
#endif

      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
         iflag = 1
         return
      end if

      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
c        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do


      if(nlevel .eq. 1)then
c         flevel(1) = fmin + df / 100000.
c         flevel(1) = fmin + df /  35000.
          flevel(1) = fmin + df /  dfquotient
      endif


      if(nlevel.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif


      if(nlevel.eq.10)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
         flevel(5)=5.0
         flevel(6)=6.0
         flevel(7)=7.0
         flevel(8)=8.0
         flevel(9)=9.0
         flevel(10)=10.0
      endif

c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) exit
         nlevlt = nlevlt + 1
      end do

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps

      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

c--Advance graphics frame and get ready to plot

      call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 0, 0)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

         call pgsci(nyellow)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &       nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then

         call pgsci(nblue)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &       nlevgt, tr)
      endif

      call pgsci(nblack)
      call pglab(titx, tity, title)

  310 format(1p,6e12.4)
  312 format(i10, 1p,6e12.4)


      return
      end
c
c*********************************************************************
c
      subroutine ezconz(r, theta, f, flevel, nr, nth, nlevel,
     &   nrmax, nthmax, nlevmax, title, titx, tity)

      implicit none

      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &   nndm1

      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &     ncolelec, j

      real:: fmin,ymax,df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, theta, r, flevel, f, xmin, ymin, xmax

      real:: xmaxz, xminz, ymaxz, yminz

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character(8):: xopt,yopt
      character(32):: title
      character(32):: titx
      character(32):: tity

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     &    nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     &    nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     &    ncolion,ncolelec,ncollin,ncollab,ncolln2
      common/zoom/ xmaxz, xminz, ymaxz, yminz


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
!      xmin = MINVAL(r)
!      xmax = MAXVAL(r)
!      ymin = MINVAL(theta)
!      ymax = MAXVAL(theta)

c--set up contour levels

c      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)


      fmax=f(2, 1)
      fmin=fmax
      do j = 1, nth
         do i = 2, nr - 1
            if (r(i) .gt. xminz .and. r(i) .le. xmaxz .and.
     &         theta(j) .gt. yminz .and. theta(j) .le. ymaxz) then
               fmax = max(fmax, f(i, j))
               fmin = min(fmin, f(i, j))
            end if
         end do
      end do

      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)return


      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
c        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel .eq. 1)flevel(1) = 1.0e-03

      if(nlevel.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif

c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) exit
         nlevlt = nlevlt + 1
      end do

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


c--Advance graphics frame and get ready to plot
c     call pladv(0)
c     call plcol(ncolbox)

      xmin = xminz
      xmax = xmaxz
      ymax = ymaxz
      ymin = yminz

c     call plenv(xmin, xmax, ymin, ymax, 1, 0)
c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps
c     call pllsty(1)
c     call plvpas(0.15,0.85,0.15,0.85,0.0)
c     call plwind(xpmin,xpmax,ypmin,ypmax)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0
c     call plbox(xopt,xtick,nxsub,yopt,ytick,nysub)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

c        call plwid(0.1)
c     call pllsty(4)
c     call plcol(ncolln2)
c        call plcol(nyellow)

c     call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
c     &      nlevlt, r, theta)
      endif

      if(nlevgt .gt. 0) then

c        call plwid(10)
c     call pllsty(1)
c     call plcol(ncollin)
c        call plcol(nblue)

c     call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
c     &      nlevgt, r, theta)
      endif

c     call plcol(ncollab)
c     call pllab(titx,tity,title)
      return
      end
c
c*********************************************************************
c
      subroutine ezcon3d(r,theta,f,flevel,nr,nth,nlevel,
     &   nrmax,nthmax,nlevmax,title,titx,tity,titz)

      implicit none

      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &   nndm1

      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &    ncolelec

      real:: fmin,ymax,df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, xmin, ymin, xmax

      real:: r(nrmax), theta(nthmax)
      real:: f(nrmax, nthmax)
      real:: flevel(nlevmax)
      character(8):: xopt,yopt
      character(32):: title
      character(32):: titx
      character(32):: tity
      character(32):: titz

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     &    nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     &    nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     &    ncolion,ncolelec,ncollin,ncollab,ncolln2

      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
!      xmin = MINVAL(r)
!      xmax = MAXVAL(r)
!      ymin = MINVAL(theta)
!      ymax = MAXVAL(theta)

c--set up contour levels
      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)return
      df = abs(fmax - fmin) / float(nlevel)

c      write(6, 100)fmin, fmax
  100 format(1p,8e12.4)


c--Advance graphics frame and get ready to plot
c     call plcol(ncolbox)
c     call plenv(-2.5, 2.5, -2.5, 4.0, 0, -2)
c     call plenv(-2.3, 1.8, -1.5, 3.5, 0, -2)

c     call plw3d(2.,4.,3.,xmin,xmax,ymin,ymax,fmin,fmax,30.,30.)
c     call plschr(0.,.50)
c     call plbox3('binstu',titx,0.0,0,'binstu',
c    *     tity,0.0,0,'binstu',titz,0.0,0)
c     call plcol(ngreen)
c     call plmesh(r,theta,f,nr,nth,3,nrmax)
c     call plschr(0.,1.0)
c     call plcol(ncollab)
c     call plmtex('t',1.0,0.5,0.5,title)
      return
      end

c
c********************************************************************
c
      subroutine boundary(r, theta, f, flevel, nr, nth, nlevel,
     &   nrmax, nthmax, nlevmax, title, titx, tity, scale)

      implicit none

!     real, intent(in), optional:: scale
      real:: scale
      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &   nndm1, norange

      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &    ncolelec, nlevelb

      real:: fmin,ymax,df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, theta, r, flevel, f, xmin, ymin, xmax, rhoplasm

      real:: tr(6), dx, dy, rin

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character(8):: xopt,yopt
      character(32):: title
      character(32):: titx
      character(32):: tity

      real:: scalex
      common/boundcom/rhoplasm

      scalex=scale
!      if (present(scale)) scalex=scale
      

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      dx = (r(2) - r(1))*scalex
      dy = theta(2) - theta(1)

      tr(1) = r(1)*scalex - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      nlevelb = 1

      call a1mnmx(r*scalex, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
!      xmin = MINVAL(r)
!      xmax = MAXVAL(r)
!      ymin = MINVAL(theta)
!      ymax = MAXVAL(theta)

c--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

      if(nlevelb.eq.1)then
         flevel(1)= rhoplasm
      endif

c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevelb
         if(flevel(i) .gt. 0.) exit
         nlevlt = nlevlt + 1
      end do

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevelb - nlevlt

c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps


      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

c--Advance graphics frame and get ready to plot

      call pgsci(nblack)
      call pgsls(1)

      CALL PGVSTD
      CALL PGWNAD (xmin, xmax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &       nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &       nlevgt, tr)
      endif

      return
      end
c
c*********************************************************************
c
      subroutine ezcon3dv(r,theta,f,flevel,nr,nth,nlevel,
     &   nrmax,nthmax,nlevmax,title,titx,tity,titz)

      implicit none

      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &   nndm1

      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &     ncolelec

      real:: fmin,ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, xmin, ymin, xmax


      real:: r(nrmax), theta(nthmax)
      real:: f(nrmax, nthmax)
      real:: flevel(nlevmax)
      character(8):: xopt,yopt
      character(32):: title
      character(32):: titx
      character(32):: tity
      character(32):: titz

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     &    nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     &    nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     &    ncolion,ncolelec,ncollin,ncollab,ncolln2

      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
!      xmin = MINVAL(r)
!      xmax = MAXVAL(r)
!      ymin = MINVAL(theta)
!      ymax = MAXVAL(theta)

c--set up contour levels
      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)return
      df = abs(fmax - fmin) / float(nlevel)

c      write(6, 100)fmin, fmax
  100 format(1p,8e12.4)


c--Advance graphics frame and get ready to plot
c     call plcol(ncolbox)
c     call plenv(-2.5, 2.5, -2.5, 4.0, 0, -2)
c        call plenv(-2.3, 1.8, -1.5, 3.5, 0, -2)
c     call plw3d(3.,3.,3.,xmin,xmax,ymin,ymax,fmin,fmax,30.,30.)

c     call plschr(0.,.50)
c     call plbox3('binstu',titx,0.0,0,'binstu',
c     *     tity,0.0,0,'binstu',titz,0.0,0)
c     call plcol(ngreen)
c     call plmesh(r,theta,f,nr,nth,3,nrmax)
c     call plschr(0.,1.0)
c     call plcol(ncollab)
c     call plmtex('t',1.0,0.5,0.5,title)
      return
      end

c
c********************************************************************
c

      subroutine ezplot3(title, titll, titlr, x1, y1, y2, y3, y4,
     &    nr, nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax)
      real:: y1max, y2max, y3max, y4max
      real:: y1min, y2min, y3min, y4min
      real:: ymin, ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)

      ymax = max(y1max, y2max, y3max, y4max)
c      ymax = 1.0e+06

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=(ymax+0.1)*1.1       !JCW July 2023. If nzeta_wdot[ei]=0, all yn are zero so add 0.1 in case
      !ymax=ymax*1.1
c      ymax = .1
      ymin=0.0



c     call pladv(0)
c     call pllsty(1)
c     call plcol(ncolbox)
c     call plvpor(0.15,0.85,0.15,0.85)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

c     call plcol(nblue)
c     call plline(nr,x1,y1)
c     call plptex(0.65*rhomax,.95*ymax,1.0,0.0,0.0,'minority ions')

c     call plcol(nred)
c     call plline(nr,x1,y2)
c     call plptex(0.65*rhomax,.89*ymax,1.0,0.0,0.0,'electrons')

c     call plcol(ncyan)
c     call plline(nr,x1,y3)
c     call plptex(0.65*rhomax,.83*ymax,1.0,0.0,0.0,'majority ions')

c     call plcol(ngreen)
c     call plline(nr, x1, y4)
c     call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'ion3')


c     call plcol(ncollin)
c     call plmtex('l',5.0,0.5,0.5,titll)
c     call plmtex('b',3.2,0.5,0.5,'rho')
c     call plmtex('t',2.0,0.5,0.5,title)

c     call pllsty(1)
c     call plcol(ncyan)
c     call plpoin(n12,x2,y3,4)
c     call pllsty(1)
c     call plcol(ncolbox)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox(' ',0.0,0,'cstv',0.0,0)
c     call plcol(3)
c     call pllsty(1)
c     call plline(nr,x1,y2)
c     call plpoin(n12,x2,y4,4)

c     call pllsty(1)
c     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p,9e11.3)
      return
      end

c
c***************************************************************************
c

      subroutine ezplot5p(title, titll, titlr, x1, y1,
     &   y2, y3, y4, y5,
     &   nr, nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax),
     &     y5(nrmax)
      real:: y1max, y2max, y3max, y4max, y5max
      real:: y1min, y2min, y3min, y4min, y5min
      real:: ymin, ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)

      ymax = max(y1max, y2max, y3max, y4max, y5max)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=(ymax+0.1)*1.1       !JCW July 2023. If nzeta_wdot[ei]=0, all yn are zero so add 0.1 in case
      !ymax=ymax*1.1
      ymin=0.0



c     call pladv(0)
c     call pllsty(1)
c     call plcol(ncolbox)
c     call plvpor(0.15,0.85,0.15,0.85)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

c     call plcol(nblue)
c     call plline(nr,x1,y1)
c     call plptex(0.65*rhomax,.95*ymax,1.0,0.0,0.0,'minority ions')

c     call plcol(nred)
c     call plline(nr,x1,y2)
c     call plptex(0.65*rhomax,.89*ymax,1.0,0.0,0.0,'electrons')

c     call plcol(ncyan)
c     call plline(nr,x1,y3)
c     call plptex(0.65*rhomax,.83*ymax,1.0,0.0,0.0,'majority ions')

c     call plcol(ngreen)
c     call plline(nr, x1, y4)
c     call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'ion3')

c     call plcol(ngreen)
c     call plline(nr, x1, y5)
c     call plptex(0.65*rhomax,.71*ymax,1.0,0.0,0.0,'fast ions')


c     call plcol(ncollin)
c     call plmtex('l',5.0,0.5,0.5,titll)
c     call plmtex('b',3.2,0.5,0.5,'rho')
c     call plmtex('t',2.0,0.5,0.5,title)

c     call pllsty(1)
c     call plcol(ncyan)
c     call plpoin(n12,x2,y3,4)
c     call pllsty(1)
c     call plcol(ncolbox)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox(' ',0.0,0,'cstv',0.0,0)
c     call plcol(3)
c     call pllsty(1)
c     call plline(nr,x1,y2)
c     call plpoin(n12,x2,y4,4)

c     call pllsty(1)
c     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p,9e11.3)
      return
      end

c
c***************************************************************************
c


      subroutine ezplot7(title, titll, titlr, titlb, x1,
     &   y1, y2, y3, y4, y5, y6, ye,
     &   nr, nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax),
     &       y5(nrmax), y6(nrmax), ye(nrmax)
      real:: y1max, y2max, y3max, y4max, y5max, y6max, yemax
      real:: y1min, y2min, y3min, y4min, y5min, y6min, yemin
      real:: ymin, ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &   ncolelec, ncolln2, ncollin, ncolbrd

      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &     nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

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
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)
      call a1mnmx(y6, nrmax, nr, y6min, y6max)
      call a1mnmx(ye, nrmax, nr, yemin, yemax)

      ymax = max(y1max, y2max, y3max, y4max, y5max, y6max, yemax)
      ymin = min(y1min, y2min, y3min, y4min, y5min, y6min, yemin)
c      ymin=0.0
#ifdef DEBUG      
      write(*,*) 'DEBUG ezplot7: ', title, ymin, ymax
#endif
      if (ymin==0. .and. ymax==0.) then
         write(*,*) 'ezplot7 no yrange: ',title,ymin,ymax
         return
      end if

      rhomax=x1(nr)
      rhomin=x1(1)
!      ymax=(ymax+0.1)*1.1 !JCW July 2023. If nzeta_wdot[ei]=0, all yn are zero so add 0.1 in case

c      ymax = 2.8e+06

c Advance plotter to a new page, define coordinate range of graph and draw axes
c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

c Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)

c Plot the line graph.

      call pgsci(nred)
      call pgline(nr, x1, ye)

      call pgsci(ncyan)
      call pgline(nr, x1, y1)

      call pgsci(nblue)
      call pgline(nr, x1, y2)

      call pgsci(ngreen)
      call pgline(nr, x1, y3)

      call pgsci(nmagenta)
      call pgline(nr, x1, y4)

      call pgsci(norange)
      call pgline(nr, x1, y5)

      call pgsci(nyellow)
      call pgline(nr, x1, y6)

      call pgsci(nblack)

  300 format (1p,9e11.3)

      return
      end

c
c***************************************************************************
c

      subroutine ezplot70(title, titll, titlr, titlb, x1,
     &   y1, y2, y3, y4, y5, y6, ye,
     &   nr, nrmax)
      ! plot seven lines set min to 0
      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax),
     &     y5(nrmax), y6(nrmax), ye(nrmax)
      real:: y1max, y2max, y3max, y4max, y5max, y6max, yemax
      real:: y1min, y2min, y3min, y4min, y5min, y6min, yemin
      real:: ymin, ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb


      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd

      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &     nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

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
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)
      call a1mnmx(y6, nrmax, nr, y6min, y6max)
      call a1mnmx(ye, nrmax, nr, yemin, yemax)

      ymax = max(y1max, y2max, y3max, y4max, y5max, y6max, yemax)
      ymin = min(y1min, y2min, y3min, y4min, y5min, y6min, yemin)
#ifdef DEBUG
      write(6,*) 'DEBUG ezplot70 ', title,ymin,ymax
      if (ymax /= ymax) then
          write(6,*) 'Nan found in ymax'
          return
      end if
#endif
      ymin=0.0

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=(ymax+0.1)*1.1 !JCW July 2023. If nzeta_wdot[ei]=0, all yn are zero so add 0.1 in case

c      ymax = 2.8e+06

c Advance plotter to a new page, define coordinate range of graph and draw axes
c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

c Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)

c Plot the line graph.

      call pgsci(nred)
      call pgline(nr, x1, ye)

      call pgsci(ncyan)
      call pgline(nr, x1, y1)

      call pgsci(nblue)
      call pgline(nr, x1, y2)

      call pgsci(ngreen)
      call pgline(nr, x1, y3)

      call pgsci(nyellow)
      call pgline(nr, x1, y4)

      call pgsci(nmagenta)
      call pgline(nr, x1, y5)

      call pgsci(norange)
      call pgline(nr, x1, y6)

      call pgsci(nblack)

  300 format (1p,9e11.3)

      return
      end

c
c***************************************************************************
c
      subroutine ezplot2(title, titll, titlr, titlb,
     &     x1, y1, y2, nr, nrmax)
      !plot two 1D curves, black and  green

      implicit none

      integer:: nr, nrmax, n

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax),y1(nrmax), y2(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb


      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

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
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      y2min = y2min * 1.1
      y2max = y2max * 1.1      
      ymax = max(y1max, y2max)
      ymin = min(y1min, y2min)

      ymax = y1max
      ymin = y1min

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1
      if (abs(ymin-ymax) .lt. 1e-3) ymax=ymin+1.e-3
#ifdef DEBUG
      write(*,*) 'ezplot2:  title: ', title, ymin,ymax
#endif
      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BNST', 0.0, 0)
      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)

      CALL PGSCI(nblue)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)
      call pgline(nr, x1, y1)

      if(y2max .ne. 0.0 .and. y2min .ne. 0.0) then
         CALL PGSWIN (rhomin, rhomax, y2min, y2max)
         CALL PGSCI(nblack)
         CALL PGBOX  (' ', 0.0, 0, 'CMST', 0.0, 0)

         CALL PGSCI(ngreen)
         call pgline(nr, x1, y2)
         call PGMTXT ('r', 2.0, 0.5, 0.5, titlr)
         CALL PGSCI(nblack)
      else
         CALL PGSCI(nblack)
         CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
      end if


  300 format (1p,9e11.3)
      return
      end
c
c***************************************************************************
c

      subroutine ezplot2p(title, titll, titlr, titlb,
     &                                       x1, y1, y2, nr, nrmax)

      implicit none

      integer:: nr, nrmax, n

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax),y1(nrmax), y2(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb


      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

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
      call a1mnmx(y2, nrmax, nr, y2min, y2max)

      ymax = max(y1max, y2max)
      ymin = min(y1min, y2min)

      if(ymax .eq. 0.0 .and. ymin .eq. 0.0)return

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)

      CALL PGSCI(nred)
      call pgline(nr, x1, y1)

      CALL PGSCI(nblue)
      call pgline(nr, x1, y2)

      CALL PGSCI(nblack)

  300 format (1p,9e11.3)
      return
      end
c
c***************************************************************************
c


      subroutine ezplot2q(title, titll, titlr, titlb,
     &                                       x1, y1, y2, nr, nrmax)

      implicit none

      integer:: nr, nrmax, n

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax),y1(nrmax), y2(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

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
      call a1mnmx(y2, nrmax, nr, y2min, y2max)

      ymax = max(y1max, y2max)
      ymin = min(y1min, y2min)
      
#ifdef DEBUG
      write(*,*) 'DEBUG ezplot2q:  title: ', title, ymin,ymax
#endif
      if(ymax .eq. 0.0 .and. ymin .eq. 0.0) return

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin
#ifdef DEBUG
      write(*,*) 'DEBUG ranges ezplot2q: ', rhomin,rhomax,ymin,ymax
      if (ymax < 1e-20) then
         write(*,*) 'ymax too small',y1,y2
         return
      end if
#endif
      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)
      call pgmtxt('r', 2.0, 0.5, 0.5, titlr)

      CALL PGSCI(nblue)
      call pgline(nr, x1, y1)

      CALL PGSCI(nmagenta)
      call pgline(nr, x1, y2)

      CALL PGSCI(nblack)

  300 format (1p,9e11.3)
      return
      end
c
c***************************************************************************
c

      subroutine ezplot4(title, titll, titlr, x1, y1, y2, y3, y4,
     &   nr, nrmax)

      implicit none
      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax),y1(nrmax),y2(nrmax),y3(nrmax),y4(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: y4max, y4min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      call a1mnmx(y2,nrmax,nr,y2min,y2max)
      call a1mnmx(y3,nrmax,nr,y3min,y3max)
      call a1mnmx(y4,nrmax,nr,y4min,y4max)

      ymax = max(y1max,y2max,y3max,y4max)
      ymin = min(y1min,y2min,y3min,y4min)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin * 1.1



c     call pladv(0)
c     call pllsty(1)
c     call plcol(ncolbox)
c     call plvpor(0.15,0.85,0.15,0.85)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

c     call plcol(nblue)
c     call plline(nr,x1,y1)
c     call plptex(0.65*rhomax,.95*ymax,1.0,0.0,0.0,'minority ions')

c     call plcol(nred)
c     call plline(nr,x1,y2)
c     call plptex(0.65*rhomax,.89*ymax,1.0,0.0,0.0,'electrons')

c     call plcol(ncyan)
c     call plline(nr,x1,y3)
c     call plptex(0.65*rhomax,.83*ymax,1.0,0.0,0.0,'majority ions')

c     call plcol(ngreen)
c     call plline(nr,x1,y4)
c     call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'fast ions')

c     call plcol(ncollin)
c     call plmtex('l',5.0,0.5,0.5,titll)
c     call plmtex('b',3.2,0.5,0.5,'rho')
c     call plmtex('t',2.0,0.5,0.5,title)

c     call pllsty(1)
c     call plcol(ncyan)
c     call plpoin(n12,x2,y3,4)
c     call pllsty(1)
c     call plcol(ncolbox)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox(' ',0.0,0,'cstv',0.0,0)
c     call plcol(3)
c     call pllsty(1)
c     call plline(nr,x1,y2)
c     call plpoin(n12,x2,y4,4)

c     call pllsty(1)
c     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p,9e11.3)
      return
      end
c
c***************************************************************************
c
      subroutine ezplot5(title,titll,titlr,x1,
     &   y1,y2,y3,y4,y5,y6,nr,nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax),y1(nrmax),y2(nrmax),y3(nrmax)
      real:: y4(nrmax),y5(nrmax),y6(nrmax)
      real:: y1max,y2max,y3max,y4max,y5max,y6max
      real:: y1min,y2min,y3min,y4min,y5min,y6min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: ncolln4,ncolln5,ncolln6
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen
      ncolln4=nred
      ncolln5=ngrey
      ncolln6=nbrown


      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      call a1mnmx(y2,nrmax,nr,y2min,y2max)
      call a1mnmx(y3,nrmax,nr,y3min,y3max)
      call a1mnmx(y4,nrmax,nr,y4min,y4max)
      call a1mnmx(y5,nrmax,nr,y5min,y5max)
      call a1mnmx(y6,nrmax,nr,y6min,y6max)
      ymax = max(y1max,y2max,y3max,y4max,y5max,y6max)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=(ymax+0.1)*1.1 !JCW July 2023. If nzeta_wdot[ei]=0, all yn are zero so add 0.1 in case
      ymin=0.0



c     call pladv(0)
c     call pllsty(1)
c     call plcol(ncolbox)
c     call plvpor(0.15,0.85,0.15,0.85)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

c     call plcol(ncollin)
c     call pllsty(3)
c     call plline(nr,x1,y1)
c     call plptex(0.6*rhomax,.95*ymax,1.0,0.0,0.0,'Mode -2')

c     call plcol(ncolln2)
c     call pllsty(2)
c     call plline(nr,x1,y2)
c     call plptex(0.6*rhomax,.90*ymax,1.0,0.0,0.0,'Mode -1')

c     call plcol(ncolln3)
c     call pllsty(1)
c     call plline(nr,x1,y3)
c     call plptex(0.6*rhomax,.85*ymax,1.0,0.0,0.0,'Mode 0')

c     call plcol(ncolln4)
c     call pllsty(2)
c     call plline(nr,x1,y4)
c     call plptex(0.6*rhomax,.80*ymax,1.0,0.0,0.0,'Mode 1')

c     call plcol(ncolln5)
c     call pllsty(3)
c     call plline(nr,x1,y5)
c     call plptex(0.6*rhomax,.75*ymax,1.0,0.0,0.0,'Mode 2')



c     call plcol(ncollin)
c     call plmtex('l',5.0,0.5,0.5,titll)
c     call plmtex('b',3.2,0.5,0.5,'rho')
c     call plmtex('t',2.0,0.5,0.5,title)

c     call pllsty(1)
c     call plcol(ncyan)
c     call plpoin(n12,x2,y3,4)
c     call pllsty(1)
c     call plcol(ncolbox)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox(' ',0.0,0,'cstv',0.0,0)
c     call plcol(3)
c     call pllsty(1)
c     call plline(nr,x1,y2)
c     call plpoin(n12,x2,y4,4)

c     call pllsty(1)
c     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p,9e11.3)
      return
      end

c
c***************************************************************************
c
      subroutine ezplot1(title, titll, titlr, x1, y1, nr, nrmax)

      implicit none

      integer:: nr, nrmax, n

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab

      ncolln3=ngreen

      call a1mnmx(y1, nrmax, nr, y1min, y1max)

#ifdef DEBUG
      write(*,*) 'ezplot1:  title: ', title, y1min,y1max
#endif
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0) then
         return
      end if
      
      ymax = y1max
      ymin = y1min

      rhomax = x1(nr)
      rhomin = x1(1)

c      write(6, *) "rhomin = ", rhomin
c      write(6, *) "rhomax = ", rhomax

      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1
      if (abs(ymin-ymax) .lt. 1e-3) ymax=ymin+1.e-3
c Advance plotter to a new page, define coordinate range of graph and draw axes

c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


c Label the axes (note use of \u and \d for raising exponent).

      call pglab('rho', titll, title)

c Plot the line graph.

      call pgline(nr, x1, y1)

  300 format (1p,9e11.3)

      return
      end

c
c***************************************************************************
c
      subroutine ezplot1_red(title, titll, titlr, x1, y1, nr, nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab

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
      if (abs(ymin-ymax) .lt. 1e-3) ymax=ymin+1.e-3
c Advance plotter to a new page, define coordinate range of graph and draw axes

c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


c Label the axes (note use of \u and \d for raising exponent).

      call pglab('rho', titll, title)

c Plot the line graph.

      call pgsci(nred)
      call pgline(nr, x1, y1)
      call pgsci(nblack)

  300 format (1p,9e11.3)

      return
      end
c
c***************************************************************************
c

      subroutine ezplot10(title, titll, titlr, x1, y1, nr, nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab

      ncolln3=ngreen

      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = x1(nr)
      rhomin = x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      ymin = 0.0

      if (ymin .le. 0.0) ymin = ymin * 1.1
      if (abs(ymin-ymax) .lt. 1e-3) ymax=ymin+1.e-3
c Advance plotter to a new page, define coordinate range of graph and draw axes

c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


c Label the axes (note use of \u and \d for raising exponent).

      call pglab('rho', titll, title)

c Plot the line graph.

      call pgline(nr, x1, y1)

  300 format (1p,9e11.3)

      return
      end
c
c***************************************************************************
c

      subroutine ezplot0(title, titll, titlr, x1, y1, nr, nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab

      ncolln3=ngreen

      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = x1(nr)
      rhomin = x1(1)
c      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1
      if (abs(ymin-ymax) .lt. 1e-3) ymax=ymin+1.e-3
c      ymin = 0.0

c Advance plotter to a new page, define coordinate range of graph and draw axes

c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


c Label the axes (note use of \u and \d for raising exponent).

      call pglab('rho', titll, title)

c Plot the line graph.

      call pgline(nr, x1, y1)

  300 format (1p,9e11.3)

      return
      end
c
c***************************************************************************
c


      subroutine ezplot1q(title, titll, titlr, titlb, x1, y1, nr, nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer:: iflag
      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab

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
      if (abs(ymin-ymax) .lt. 1e-3) ymax=ymin+1.e-3

#ifdef DEBUG
      write(*,*) 'ezplot1q title: ', title
      write(6, *)"ymax = ", ymax, "   ymin = ", ymin
#endif

      iflag = 0
      if(ymax == 0.0 .and. ymin == 0.0) then
         ymax = 1.0
         iflag = 1
      end if
      
c Advance plotter to a new page, define coordinate range of graph and draw axes

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


c Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)

c Plot the line graph.
      call pgsci(nblue)
      call pgline(nr, x1, y1)
      call pgsci(nblack)

  300 format (1p,9e11.3)

      return
      end
c
c***************************************************************************
c


      subroutine ezplot1p(title, titll, titlr, x1, y1, nr, nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax),y1(nrmax)
      real:: y1max,y2max,y3max,y1min,y2min,y3min
      real:: ymin,ymax

      character(32):: title
      character(32):: titll
      character(32):: titlr

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd
      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     & nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     & nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     & ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      if(y1max .eq. 0.0 .and. y2min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1

c     call pladv(0)
c     call pllsty(1)
c     call plcol(ncolbox)
c     call plvpor(0.15,0.85,0.15,0.85)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

c     call plcol(ncollin)
c     call plline(nr,x1,y1)


c     call plcol(ncollin)
c     call plmtex('l',5.0,0.5,0.5,titll)
c     call plmtex('b',3.2,0.5,0.5,'rho')
c     call plmtex('t',2.0,0.5,0.5,title)

c     call pllsty(1)
c     call plcol(ncyan)

c     call pllsty(1)
c     call plcol(ncolbox)
c     call plwind(rhomin,rhomax,ymin,ymax)
c     call plbox(' ',0.0,0,'cstv',0.0,0)


  300 format (1p,9e11.3)
      return
      end
c
c***************************************************************************
c
      subroutine ezconcx(r, theta, f, flevel, nr, nth, nlevel0,
     &   nrmax, nthmax, nlevmax, title, titx, tity, iflag)

      implicit none

      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &   nndm1, norange

      real:: fmin,ymax,df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, theta, r, flevel, f, xmin, ymin, xmax, fact
      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &    ncolelec, iflag, nlevel0

      real:: tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character(8):: xopt,yopt
      character(32):: title
      character(32):: titx
      character(32):: tity

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      dx = r(2) - r(1)
      dy = theta(2) - theta(1)

      tr(1) = r(1) - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
!      xmin = MINVAL(r)
!      xmax = MAXVAL(r)
!      ymin = MINVAL(theta)
!      ymax = MAXVAL(theta)

c--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

#ifdef DEBUG
      write(*,*) 'title: ', title
#endif
      write(6, *)"fmax = ", fmax, "   fmin = ", fmin
      write(15,*)"fmax = ", fmax, "   fmin = ", fmin

      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
c         write(6, *)"fmax = ", fmax
c         write(6, *)"fmin = ", fmin
         iflag = 1
         return
      end if

      df = abs(fmax - fmin) / float(nlevel0)
      do i = 1, nlevel0
c        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel0.eq.1)flevel(1)=1.0e-03

      if(nlevel0.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif





      if (nlevel0 .eq. 21) then

c         nlevel = nlevel0 * 2.0

c         fact = 2.0
c         i = 1
c         flevel(i) = fmin

c         do i = 2, nlevel / 2
c            fact = fact * 1.2
c            flevel(i) = flevel(i-1) / fact
c         end do

c         fact = 2.0
c         i = nlevel / 2 + 1
c         flevel(i) = fmax

c         do i = nlevel / 2 + 2, nlevel
c            fact = fact * 1.2
c            flevel(i) = flevel(i-1) / fact
c         end do

         nlevel = nlevel0
         fact = 2.0
         i = 1
         flevel(i) = fmax

         do i = 2, nlevel
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do

      else

         nlevel = nlevel0

      end if








c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) exit
         nlevlt = nlevlt + 1
      end do

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps

      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

c--Advance graphics frame and get ready to plot

      call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 1, 0)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

         call pgsci(nyellow)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &       nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then

         call pgsci(nblue)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &       nlevgt, tr)
      endif

      call pgsci(nblack)
      call pglab(titx, tity, title)

  310 format(1p,6e12.4)
  312 format(i10, 1p,6e12.4)


      return
      end
c
c*********************************************************************
c

      subroutine ezlog1(title, titll, titlr, titlb, x1, y1,
     &   nr, nrmax, ymin, ymax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax)
      real:: y1max, y2max, y3max, y4max, y5max, y6max
      real:: y1min, y2min, y3min, y4min, y5min, y6min
      real:: ymin, ymax, tmin, tmax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd

      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      rhomax=x1(nr)
      rhomin=x1(1)

      call a1mnmx(y1, nrmax, nr, ymin, ymax)
      if (abs(ymin-ymax) .lt. 1e-3) ymax=ymin+1.e-3
c Advance plotter to a new page, define coordinate range of graph and draw axes
c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)

      CALL PGBOX  ('BCNST', 0.0, 0, 'LBCNST', 0.0, 0)

      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)

      call pgsci(nblue)
      call pgline(nr, x1, y1)



      CALL PGSCI(nblack)
      call pgsls(1)

  300 format (1p,9e11.3)
      return
      end

c
c***************************************************************************
c

      subroutine ezlog1_f(title, titll, titlr, titlb, x1, y1,
     &   y2, y3, y4, y5, y6, y7,
     &   nr, nrmax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real:: x1(nrmax), y1(nrmax)
      real:: y2(nrmax), y3(nrmax), y4(nrmax)
      real:: y5(nrmax), y6(nrmax), y7(nrmax)
      real:: y1max, y2max, y3max, y4max, y5max, y6max
      real:: y1min, y2min, y3min, y4min, y5min, y6min
      real:: ymin, ymax, tmin, tmax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd

      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      rhomax=x1(nr)
      rhomin=x1(1)

      call a1mnmx(y1, nrmax, nr, ymin, ymax)

      ymin = -10.0

c Advance plotter to a new page, define coordinate range of graph and draw axes
c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)

      CALL PGBOX  ('BCNST', 0.0, 0, 'LBCNST', 0.0, 0)

      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)

      call pgsci(nblue)
      call pgline(nr, x1, y1)
      call pgline(nr, x1, y2)
      call pgline(nr, x1, y3)
      call pgline(nr, x1, y4)
      call pgline(nr, x1, y5)
      call pgline(nr, x1, y6)
      call pgline(nr, x1, y7)


      CALL PGSCI(nblack)
      call pgsls(1)

  300 format (1p,9e11.3)
      return
      end

c
c***************************************************************************
c



      subroutine ezconpx(r, theta, f, flevel, nr, nth, nlevel0,
     &   nrmax, nthmax, nlevmax, xg, yg, title, r0)

      implicit none

      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &   nndm1

      real:: fmin,ymax,df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, xmin, ymin, xmax
      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &   ncolelec, j, n, m, nth1, nlevel0
      real:: xmaxz, xminz, ymaxz, r0, fact

      real:: r(nrmax),theta(nthmax)
      real:: f(nrmax,nthmax)
      real:: xg(nrmax,nthmax)
      real:: yg(nrmax,nthmax)
      real:: flevel(nlevmax)

      character(8):: xopt,yopt
      character(32):: title


      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     &   ncolion,ncolelec,ncollin,ncollab,ncolln2



c--set up transformation
      do n = 1, nr
         do m = 1, nth
            xg(n,m) = r(n) * cos(theta(m)) + r0
            yg(n,m) = r(n) * sin(theta(m))
         end do
         xg(n, nth+1) = xg(n,1)
         yg(n, nth+1) = yg(n,1)
      end do
      nth1 = nth + 1
      call a2dmnmx(xg, nrmax, nthmax, nr, nth, xmin, xmax)
      call a2dmnmx(yg, nrmax, nthmax, nr, nth, ymin, ymax)



c--set up contour levels
      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)
      df=abs(fmax-fmin)/float(nlevel0)

c        write (6, *) "df = ", df

      if(df .eq. 0.0)return

      do i = 1, nlevel0
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel0.eq.1)flevel(1)=1.0e-03
      if(nlevel0.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif

      if (nlevel0 .eq. 20) then

         nlevel = nlevel0 * 2.0

         fact = 2.0
         i = 1
         flevel(i) = fmin

         do i = 2, nlevel / 2
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do


         fact = 2.0
         i = nlevel / 2 + 1
         flevel(i) = fmax

         do i = nlevel / 2 + 2, nlevel
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do

      else
         nlevel = nlevel0
      end if


c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0

      do i = 1, nlevel
         if(flevel(i) .gt. 0.) exit
         nlevlt = nlevlt + 1
      end do

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt



c--Advance graphics frame and get ready to plot

c      call plcol(ncolbox)
c      call plenv(xmin, xmax, ymin, ymax, 0, 0)

      call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 1, 0)

c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps


c      call pllsty(1)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then
c         call plwid(1)
c         call pllsty(2)
c         call plcol(ncolln2)  #plplot routines
c         call pgsci(nyellow)
c         call plcon2(f, nrmax, nthmax, 1, nr, 1, nth1, flevel(ilevlt),
c     &       nlevlt, xg, yg)
      endif

      if(nlevgt .gt. 0) then
c         call plwid(1)
c         call pllsty(1)
c         call plcol(ncollin)
c         call plcol(nblue)
c         call plcon2(f,nrmax,nthmax,1,nr,1,nth1,flevel(ilevgt),
c     &   nlevgt, xg, yg)
      endif
c      call plwid(1)


c      call plcol(ncollab)
c      call plcol(nblue)
c      call pllab('u_parallel','u_perp', title)


  310 format(1p,6e12.4)
  311 format(10i10)

      return
      end
c
c*********************************************************************
c

      subroutine pseudo(r, theta, f, flevel, nr, nth, nlevel,
     &   nrmax, nthmax, nlevmax, title, titx, tity, iflag)

      implicit none

      integer:: nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     &   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     &   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     &   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     &   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     &   nndm1, norange

      real:: fmin,ymax,df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     &   ypmin, theta, r, flevel, f, xmin, ymin, xmax
      integer:: nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     &   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     &    ncolelec, iflag

      real:: tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character(8):: xopt,yopt
      character(32):: title
      character(32):: titx
      character(32):: tity

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      dx = r(2) - r(1)
      dy = theta(2) - theta(1)

      tr(1) = r(1) - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
!      xmin = MINVAL(r)
!      xmax = MAXVAL(r)
!      ymin = MINVAL(theta)
!      ymax = MAXVAL(theta)

c--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

c      write(6, *)"fmax = ", fmax, "   fmin = ", fmin
c      write(15,*)"fmax = ", fmax, "   fmin = ", fmin

      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
c         write(6, *)"fmax = ", fmax
c         write(6, *)"fmin = ", fmin
         iflag = 1
         return
      end if


c--Advance graphics frame and get ready to plot

      call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 1, 0)

      call pgscir(1, nlevel)
      call palett(2, 1., 0.5)
      call pgimag(f, nrmax, nthmax, 1, nr, 1, nth, fmin, fmax, tr)

!      CALL PGSVP (0.15,0.85,0.85,0.95)'
      call palett(2, 1., 0.5)
      CALL PGWEDG('RI', 1.0, 2.5, FMIN, FMAX, 'Values')

      call pgsci(nblack)
      call pglab(titx, tity, title)

  310 format(1p,6e12.4)
  312 format(i10, 1p,6e12.4)


      return
      end
c
c*********************************************************************
c


      subroutine ezzoom1(title, titll, titlr, titlb, x1, y1,
     &   y2, y3, y4, y5, y6,
     &   nr, nrmax, ymin, ymax, rhomin, rhomax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax, xzmin, xnmin, xnmax, rhomin, rhomax
      real:: x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax),
     &     y5(nrmax), y6(nrmax)
      real:: y1max, y2max, y3max, y4max, y5max, y6max
      real:: y1min, y2min, y3min, y4min, y5min, y6min
      real:: ymin, ymax, tmin, tmax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd

      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      if (abs(ymin-ymax) .lt. 1e-3) ymax=ymin+1.e-3
c Advance plotter to a new page, define coordinate range of graph and draw axes
c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)


c Plot the line graph.

      call pgsci(nblue)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)
      call pgline(nr, x1, y1)

c      tmax = max(y2max, y3max)
c      tmin = 0.0

c      CALL PGSWIN (rhomin, rhomax, tmin, tmax)
c      CALL PGSCI(nblack)
c      CALL PGBOX  (' ', 0.0, 0, 'CMST', 0.0, 0)

c      CALL PGSCI(nred)
c      call pgsls(2)
c      call pgline(nr, x1, y2)
c      call PGMTXT ('r', 2.0, 0.5, 0.5, titlr)


      call pgsci(nblack)
      call pgsls(1)

  300 format (1p,9e11.3)

      return
      end

c
c*********************************************************************
c


      subroutine ezzoom6(title, titll, titlr, titlb, x1, y1,
     &   y2, y3, y4, y5, y6,
     &   nr, nrmax, rhomin, rhomax)

      implicit none

      integer:: nr,nrmax

      real:: xzmax, xzmin, xnmin, xnmax, rhomin, rhomax
      real:: x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax),
     &     y5(nrmax), y6(nrmax)
      real:: y1max, y2max, y3max, y4max, y5max, y6max
      real:: y1min, y2min, y3min, y4min, y5min, y6min
      real:: ymin, ymax, tmin, tmax

      character(32):: title
      character(32):: titll
      character(32):: titlr
      character(32):: titlb

      integer:: nplot1,ncollab, ncolion,ncolbox, ncyan,
     &    ncolelec, ncolln2, ncollin, ncolbrd

      integer:: nblack,nred,nyellow, ngreen,naqua,npink,
     &   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     &   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

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
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)
      call a1mnmx(y6, nrmax, nr, y6min, y6max)

      ymax = max(y1max, y2max, y3max, y4max, y5max, y6max)
      ymin = min(y1min, y2min, y3min, y4min, y5min, y6min)
c      ymin=0.0

      ymax = ymax * 0.5
      if (abs(ymin-ymax) .lt. 1e-3) ymax=ymin+1.e-3
c Advance plotter to a new page, define coordinate range of graph and draw axes
c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BNST', 0.0, 0)
      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)


c Plot the line graph.

      call pgsci(nblue)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)
      call pgline(nr, x1, y1)

      tmax = max(y2max, y3max) * 0.5
      tmin = 0.0

      CALL PGSWIN (rhomin, rhomax, tmin, tmax)
      CALL PGSCI(nblack)
      CALL PGBOX  (' ', 0.0, 0, 'CMST', 0.0, 0)

      CALL PGSCI(nred)
      call pgsls(1)
      call pgline(nr, x1, y2)
      call PGMTXT ('r', 2.0, 0.5, 0.5, titlr)

      CALL PGSCI(ngreen)
      call pgsls(3)
      call pgline(nr, x1, y3)

      CALL PGSWIN (rhomin, rhomax, y4min, y4max)
      CALL PGSCI(nblack)


      call pgsci(nblack)
      call pgsls(1)
      call pgline(nr, x1, y4)



      call pgsci(nblack)
      call pgsls(1)

  300 format (1p,9e11.3)

      return
      end

c
c***************************************************************************
c
C-----------------------------------------------------------------------
      SUBROUTINE FREDDY(ARRAY,KX,NY,SIZE,ANGLE)
      INTEGER:: KX, NY
      REAL:: ARRAY(KX,NY), SIZE, ANGLE
C
C Draws isometric plot of array
C
      REAL:: FMAX,FMIN,DELTAX,DELTAY,DELTAV,SINE,PEAK,X,DX,HEIGHT
      INTEGER:: I,J,KI,KJ,NX,MX,MY,STEP,LEFT,RIGHT,IT,MN,INCX
      LOGICAL:: VISBLE
      COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
      MN = KX*NY
      NX = KX
C     Check array size:
      IF(NX.LT.2 .OR. NY.LT.2) RETURN
      FMAX = ARRAY(1,1)
      FMIN = FMAX
      DO J=1,NY
          DO I=1,NX
              FMIN = MIN(ARRAY(I,J),FMIN)
              FMAX = MAX(ARRAY(I,J),FMAX)
           END DO
      END DO
      DELTAX = SIZE/(NX+NY)
      SINE = SIN(ANGLE/58.)
      DELTAY = DELTAX*SINE
      HEIGHT = SIZE*(1.-ABS(SINE))
      DELTAV = HEIGHT
      FMAX = FMAX-FMIN
      IF(FMAX.LT.0.0001) FMAX = DELTAV
      DELTAV = DELTAV/FMAX
      MX = NX+1
      MY = NY+1
      STEP = MX
C
C Start PGPLOT buffering.
C
      CALL PGBBUF
C
C Work our way down the Y axis, then up the X axis,
C calculating the Y plotter coordinates for each
C column of the plot, doing the hidden-line suppression
C at the same time.
C
      DO 50 J=1,NY
          KJ = MY-J
          KI = 1
C               ( KI,KJ are coordinates of bottom of column)
          ARRAY(KI,KJ) = DELTAY*(KI+KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
   30     PEAK = ARRAY(KI,KJ)
   40     KI = KI+1
          KJ = KJ+1
          IF(KI.GT.NX .OR. KJ.GT.NY) GOTO 50
          ARRAY(KI,KJ) = DELTAY*(KI+KJ) + DELTAV*(ARRAY(KI,KJ)-FMIN)
          IF(ARRAY(KI,KJ).GT.PEAK) GOTO 30
          IF(ARRAY(KI,KJ).LE.PEAK) ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
          GOTO 40
   50    CONTINUE
C
C Now to work our way up the X axis
C
      DO 80 I=2,NX
          KI = I
          KJ = 1
          ARRAY(KI,KJ) = DELTAY*(KI+KJ)+DELTAV*(ARRAY(KI,KJ)-FMIN)
   60     PEAK = ARRAY(KI,KJ)
   70     KI = KI+1
          KJ = KJ+1
          IF(KI.GT.NX .OR. KJ.GT.NY) GOTO 80
          ARRAY(KI,KJ) = DELTAY*(KI+KJ)+DELTAV*(ARRAY(KI,KJ)-FMIN)
          IF(ARRAY(KI,KJ).GT.PEAK) GOTO 60
          IF(ARRAY(KI,KJ).LE.PEAK) ARRAY(KI,KJ) = -ABS(ARRAY(KI,KJ))
          GOTO 70
   80 CONTINUE
C
C Draw a line along the bottom of the vertical faces
C
      CALL PGMOVE(DELTAX*(NX+NY-2), DELTAY*(MX))
      CALL PGDRAW(DELTAX*(NY-1),    DELTAY*2)
      CALL PGDRAW(0.0,              DELTAY*MY)
C
C Array is now ready for plotting.  If a point is
C positive, then it is to be plotted at that Y
C coordinate; if it is negative, then it is
C invisible, but at minus that Y coordinate (the point
C where the line heading towards it disappears has to
C be determined by finding the intersection of it and
C the cresting line).
C
C Plot rows:
C
      DO J=1,NY,2
          KJ = MY-J
          DX = DELTAX*(J-2)
          X = DX+DELTAX
          CALL PGMOVE(X,DELTAY*(KJ+1))
          CALL PGDRAW(X,ARRAY(1,KJ))
          VISBLE = .TRUE.
          DO I=2,NX
              RIGHT = I+NX*(KJ-1)
              LEFT = RIGHT-1
              IT = RIGHT
              X = DX+DELTAX*I
              CALL FREDGO(ARRAY,MN)
          END DO
C
C Now at far end of row so come back
C
          KJ = KJ-1
          IF(KJ.LE.0) GOTO 170
          VISBLE = ARRAY(NX,KJ).GE.0.0
          DX = DELTAX*(NX+J)
          IF(VISBLE) CALL PGMOVE(DX-DELTAX,ARRAY(NX,KJ))
          DELTAX = -DELTAX
          DO I=2,NX
              KI = MX-I
              LEFT = KI+NX*(KJ-1)
              RIGHT = LEFT+1
              IT = LEFT
              X = DX+DELTAX*I
              CALL FREDGO(ARRAY,MN)
          END DO
C
          X = DX+DELTAX*NX
          IF(.NOT.VISBLE) CALL PGMOVE(X,ARRAY(1,KJ))
          CALL PGDRAW(X,DELTAY*(KJ+1))
C               (set DELTAX positive for return trip)
          DELTAX = -DELTAX
      END DO
C
C Now do the columns:
C as we fell out of the last DO-loop we do the
C columns in ascending-X order
C
      INCX = 1
      KI = 1
C               (set DELTAX -ve since scanning R to L)
  120 DX = DELTAX*(KI+NY-1)
      DELTAX = -DELTAX
      X = DX+DELTAX
      CALL PGMOVE(X,ARRAY(1,1))
  130 VISBLE = .TRUE.
      DO J=2,NY
          LEFT = KI+NX*(J-1)
          RIGHT = LEFT-NX
          IT = LEFT
          X = DX+DELTAX*J
          CALL FREDGO(ARRAY,MN)
      END DO
C
C At far end, increment X and check still inside array
C
      KI = KI+INCX
      IF(KI.LE.0 .OR. KI.GT.NX) GOTO 180
      VISBLE = ARRAY(KI,NY).GE.0.0
      DELTAX = -DELTAX
      DX = DELTAX*(KI-2)
      X = DX+DELTAX
      IF(VISBLE) CALL PGMOVE(X,ARRAY(KI,NY))
      DO J=2,NY
          KJ = MY-J
          RIGHT = KI+NX*(KJ-1)
          LEFT = RIGHT+NX
          IT = RIGHT
          X = DX+DELTAX*J
          CALL FREDGO(ARRAY,MN)
      END DO

      X = DX+DELTAX*NY
      IF(.NOT.VISBLE) CALL PGMOVE(X,ARRAY(KI,1))
      IF(KI.EQ.1) GOTO 180
      CALL PGDRAW(X,DELTAY*(KI+1))
      KI = KI+INCX
      IF(KI.GT.NX) GOTO 180
      IF(KI.EQ.1) GOTO 120
  160 DELTAX = -DELTAX
      DX = DELTAX*(1-KI-NY)
      X = DX+DELTAX
      CALL PGMOVE(X,DELTAY*(KI+1))
      CALL PGDRAW(X,ARRAY(KI,1))
      GOTO 130
C
C Do columns backwards because ended rows at far end of X
C
  170 KI = NX
      INCX = -1
      DX = DELTAX*(KI+NY)
      GOTO 160
C
C
  180 CALL PGEBUF
      END
C-----------------------------------------------------------------------
      SUBROUTINE FREDGO(ARRAY,MN)
      INTEGER:: MN
      REAL:: ARRAY(MN)
C
      INTEGER:: STEP,LEFT,RIGHT,IT,NX
      LOGICAL:: VISBLE
      REAL:: AL,AR,BL,EM,XX,X,Y,DELTAX
      COMMON /FREDCM/ DELTAX,X,STEP,LEFT,RIGHT,IT,NX,VISBLE
C
C Test visibility
C
      IF(ARRAY(IT).LT.0.0) GOTO 80
C
C This point is visible - was last?
C
      IF(VISBLE) GOTO 50
C
C No: calculate point where this line vanishes
C
   10 IF(LEFT.LE.NX .OR. MOD(LEFT-1,NX).EQ.0 .OR.
     &     RIGHT.LE.NX .OR. MOD(RIGHT-1,NX).EQ.0) GOTO 100
      AL = ABS(ARRAY(LEFT))
      AR = ABS(ARRAY(RIGHT))
      IF(ARRAY(LEFT).LT.0.0) GOTO 70
C               Right-hand point is crested
   20 RIGHT = RIGHT-STEP
      IF(ARRAY(RIGHT).LT.0.0) GOTO 20
C               Left-hand end of cresting line is either
C               RIGHT+NX or RIGHT-1
      LEFT = RIGHT+NX
      IF(ARRAY(LEFT).LT.0.0) LEFT = RIGHT-1
C
C               RIGHT and LEFT index into the endpoints of the
C               cresting line
   30 BL = ABS(ARRAY(LEFT))
      EM = ABS(ARRAY(RIGHT))-BL
      XX = EM-AR+AL
      IF(ABS(XX).LT.0.0001) GOTO 60
      XX = (AL-BL)/XX
   40 Y = EM*XX+BL
      IF(DELTAX.GT.0.0) XX = 1.0-XX
      XX = X-XX*DELTAX
      IF(VISBLE) GOTO 90
C               Drawing a line from an invisible point
C               to a visible one
      CALL PGMOVE(XX,Y)
      VISBLE = .TRUE.
   50 CALL PGDRAW(X,ARRAY(IT))
      RETURN
C
   60 XX = 0.5
      GOTO 40
C
C Left-hand point crested
C
   70 LEFT = LEFT-STEP
      IF(ARRAY(LEFT).LT.0.0) GOTO 70
C
C Right-hand end of cresting line is either LEFT+1 or LEFT-NX
C
      RIGHT = LEFT+1
      IF(ARRAY(RIGHT).LT.0.0) RIGHT = LEFT-NX
      GOTO 30
C
C This point is invisible; if last one was too, then forget it;
C else draw a line towards it
C
   80 IF(.NOT.VISBLE) RETURN
      GOTO 10
C
   90 CALL PGDRAW(XX,Y)
  100 VISBLE = .FALSE.
      RETURN
      END


!--------------------------------------------------------------------------------
      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
C-----------------------------------------------------------------------
C Set a "palette" of colors in the range of color indices used by
C PGIMAG.
C-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      END
