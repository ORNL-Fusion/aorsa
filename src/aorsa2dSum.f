*------------------
*     AORSA2D_SUM:
*------------------

!      OPTIONS / check = no underflow

      subroutine aorsa2dSum(myid)

      use size_mod
      use aorsa2din_mod

      implicit none

      integer:: kp1, km1
      integer:: nchmax, jhalf, khalf, kshift
      integer:: iezrant, nr, k, nt1, nt2, jmid,
     &   nant, nphi3d, n, nphi, nstrpol, nphi_fpm, nt_max
      integer:: nphi3d_1quarter, nphi3d_half, nphi3d_3quarter

      integer:: i, j
      integer:: nxdim, nydim, myid
      integer:: nnodex, nnodey, nphi1, nphi2, nphi1_dum, nphi2_dum
      integer:: nxplot, nyplot, ir, nt, imax, jmax, kmax

      integer:: pgopen, pgbeg, ier, numb, iflag
      integer:: number_points
      integer:: nnoderho

      integer, parameter:: nlevmax = 101, nxmx = nmodesmax,
     & nymx = mmodesmax, nrhomax = nxmx * 2, nxplot_dim   = 1000,
     &     nyplot_dim   = 1000, nphmax = 2048,
     &     ntmax = nphmax / 2, ntmin = -ntmax

      character(32):: title, titll, titlr, titlb
      character(32):: titx, tity, titz, tityl, tityr

      real, dimension(:),   allocatable :: UPERP, UPARA
      real, dimension(:,:,:), allocatable :: bqlavg_i1
      real, dimension(:,:,:), allocatable :: cqlavg_i1
      real, dimension(:,:,:), allocatable :: eqlavg_i1
      real, dimension(:,:,:), allocatable :: fqlavg_i1

      real, dimension(:,:,:), allocatable :: bqlsum_i1
      real, dimension(:,:,:), allocatable :: cqlsum_i1
      real, dimension(:,:,:), allocatable :: eqlsum_i1
      real, dimension(:,:,:), allocatable :: fqlsum_i1

      real, dimension(:,:,:), allocatable :: bqlavg_i2
      real, dimension(:,:,:), allocatable :: cqlavg_i2
      real, dimension(:,:,:), allocatable :: eqlavg_i2
      real, dimension(:,:,:), allocatable :: fqlavg_i2

      real, dimension(:,:,:), allocatable :: bqlsum_i2
      real, dimension(:,:,:), allocatable :: cqlsum_i2
      real, dimension(:,:,:), allocatable :: eqlsum_i2
      real, dimension(:,:,:), allocatable :: fqlsum_i2

      integer:: i_uperp, i_upara, i_psi, i_psi_eq
      real:: vce_mks,vc1_mks,vc2_mks, vc3_mks, vc4_mks, vc5_mks, vc6_mks
      real:: vce_cgs,vc1_cgs,vc2_cgs, vc3_cgs, vc4_cgs, vc5_cgs, vc6_cgs
      real:: UminPara, UmaxPara, xmi1, xmi2
      real:: term1_coef, term2_coef, diff_coef


      real:: xnurf, pi, omgrf, q, width, eps0, xmu0, clight, qe, absqe,
     &     xk0, phimin, phimax, dphi, dphi3d, r0, xltcm
      real:: xi0, separ0, rhoplasm, rmax, rmin
      real:: dphase0, phase0, phi(nphmax), capr_min, capr_max
      real:: phi3d(nphimx), phi3d_shift(nphimx), zk3d(nphimx), dxplot
      real:: zk(nphmax), dyplot, dx, dy, emax
      real:: workr(nphmax), worki(nphmax), capr(nxmx), capz(nymx)
      real:: rho(nxmx, nymx), drho
      real:: xnea(nxmx, nymx)
      real:: ff(101), dummy, tmin, t1, second1, ptot, pcito2, pcrto2
      real:: capr_giv, phi3d_giv, fout, dcapr, x_giv, ptotal, xjtot
      real:: ptotal_edge, ptotal_core, fedge, fcore

      real:: pcte, pcti1, pcti2, pcti3, pcti4, pcti5, pcti6, pctt, pt

      real::  pcedotje, pcedotj1, pcedotj2, pcedotj3, pcedotj4,
     &      pcedotj5, pcedotj6, pcedotjs, pcedotjt, pedotjt


      complex:: work(nphmax), zi, csum, curent(nphmax),cn(ntmin : ntmax)
      complex:: cur3d(nphimx), jfun, cursum(nphmax), cexpkz
      complex:: div_x, div_y, div_z, div_j


      complex, dimension(:,:), allocatable :: ealpha
      complex, dimension(:,:), allocatable :: ebeta
      complex, dimension(:,:), allocatable :: eb
      complex, dimension(:,:), allocatable :: xjpx
      complex, dimension(:,:), allocatable :: xjpy
      complex, dimension(:,:), allocatable :: xjpz

      real, dimension(:,:), allocatable :: xjx
      real, dimension(:,:), allocatable :: xjy
      real, dimension(:,:), allocatable :: xjz

      real, dimension(:,:), allocatable :: bx
      real, dimension(:,:), allocatable :: by
      real, dimension(:,:), allocatable :: bz

      complex, dimension(:,:), allocatable :: xjpxe_lab
      complex, dimension(:,:), allocatable :: xjpye_lab
      complex, dimension(:,:), allocatable :: xjpze_lab

      real, dimension(:,:,:), allocatable :: sR_3d
      real, dimension(:,:,:), allocatable :: sZ_3d
      real, dimension(:,:,:), allocatable :: sphi_3d
      real:: Sx, Sy, Sz

      complex, dimension(:,:), allocatable :: ntilda_e

      complex, dimension(:,:), allocatable :: eplus
      complex, dimension(:,:), allocatable :: eminus

      complex, dimension(:,:), allocatable :: ex, ey, ez,
     &   bxwave, bywave, bzwave


      real:: rhon(nrhomax)
      real:: ntilda_max

      real:: wdoteavg(nrhomax),
     &     wdot1avg(nrhomax), wdot2avg(nrhomax),
     &     wdot3avg(nrhomax), wdot4avg(nrhomax),
     &     wdot5avg(nrhomax), wdot6avg(nrhomax)

      real::  wdote_ql(nrhomax), wdoti1_ql(nrhomax), wdoti2_ql(nrhomax),
     &     wdoti3_ql(nrhomax), wdoti4_ql(nrhomax), wdoti5_ql(nrhomax)

      real:: redotjeavg(nrhomax),
     &     redotj1avg(nrhomax), redotj2avg(nrhomax),
     &     redotj3avg(nrhomax), redotj4avg(nrhomax),
     &     redotj5avg(nrhomax), redotj6avg(nrhomax)

      real:: xjprlavg(nrhomax), xjprl_sum(nrhomax)

      real:: wdotesum(nrhomax),
     &     wdot1sum(nrhomax), wdot2sum(nrhomax),
     &     wdot3sum(nrhomax), wdot4sum(nrhomax),
     &     wdot5sum(nrhomax), wdot6sum(nrhomax)

      real:: wdotesum_ql(nrhomax),
     &     wdot1sum_ql(nrhomax), wdot2sum_ql(nrhomax),
     &     wdot3sum_ql(nrhomax), wdot4sum_ql(nrhomax),
     &     wdot5sum_ql(nrhomax)

      real:: redotjesum(nrhomax), redotjisum(nrhomax),
     &     redotj1sum(nrhomax), redotj2sum(nrhomax),
     &     redotj3sum(nrhomax), redotj4sum(nrhomax),
     &     redotj5sum(nrhomax), redotj6sum(nrhomax)

      real:: redotje_int(nrhomax),
     &     redotji1_int(nrhomax), redotji2_int(nrhomax),
     &     redotji3_int(nrhomax), redotji4_int(nrhomax),
     &     redotji5_int(nrhomax), redotji6_int(nrhomax)


      real:: wdoti1_dvol(nrhomax), wdoti2_dvol(nrhomax),
     &     wdoti3_dvol(nrhomax), wdoti4_dvol(nrhomax),
     &     wdoti5_dvol(nrhomax), wdoti6_dvol(nrhomax),
     &     wdote_dvol(nrhomax)

      real:: redotj1_dvol(nrhomax), redotj2_dvol(nrhomax),
     &     redotj3_dvol(nrhomax), redotj4_dvol(nrhomax),
     &     redotj5_dvol(nrhomax), redotj6_dvol(nrhomax),
     &     redotje_dvol(nrhomax)

      real:: wdote_int(nrhomax),
     &     wdoti1_int(nrhomax), wdoti2_int(nrhomax),
     &     wdoti3_int(nrhomax), wdoti4_int(nrhomax),
     &     wdoti5_int(nrhomax), wdoti6_int(nrhomax)

      real:: wdote_int_ql(nrhomax),
     &     wdoti1_int_ql(nrhomax), wdoti2_int_ql(nrhomax),
     &     wdoti3_int_ql(nrhomax), wdoti4_int_ql(nrhomax),
     &     wdoti5_int_ql(nrhomax), wdoti6_int_ql(nrhomax)

      real:: xjprl_int(nrhomax)

      real::  dvol(nrhomax), darea(nrhomax)


      real, dimension (:, :), allocatable :: freal
      real, dimension (:, :), allocatable :: fimag
      real, dimension (:, :), allocatable :: power
      real, dimension (:, :), allocatable :: fluct
      real, dimension (:, :), allocatable :: capr_plot

      real, dimension (:, :), allocatable :: freal_phi
      real, dimension (:, :), allocatable :: fimag_phi
      real, dimension (:, :), allocatable :: power_phi
      real, dimension (:, :), allocatable :: capr_plot_phi

      real:: fmidre(nxmx), fmidim(nxmx)

      real:: cnmod2(nphimx),xnphi(nphimx), pabs(nphimx), jdriven(nphimx)
      real:: spa(nphimx), spa_cold
      real:: pabs_weight, pabs_sum, pscale, j_driven_weight
      real:: jdriven_weight, jdriven_sum

      real, dimension (:), allocatable :: xplot
      real, dimension (:), allocatable :: yplota
      real, dimension (:), allocatable :: xplotm
      real, dimension (:), allocatable :: yplotm

      real, dimension (:, :), allocatable :: frealp
      real, dimension (:, :), allocatable :: fimagp
      real, dimension (:, :), allocatable :: powerp
      real, dimension (:, :), allocatable :: capr_plotp

      real, dimension (:, :, :), allocatable :: ealpha_sum_real
      real, dimension (:, :, :), allocatable :: ealpha_sum_imag

      real, dimension (:, :, :), allocatable :: eb_sum_real
      real, dimension (:, :, :), allocatable :: eb_sum_imag


      real, dimension (:, :, :), allocatable :: eb_sum_mod
      real, dimension (:, :, :), allocatable :: ealpha_sum_mod
      real, dimension (:, :, :), allocatable :: ebeta_sum_mod
      real, dimension (:, :, :), allocatable :: mod_e_sum

      complex, dimension(:, :, :), allocatable :: ealpha_sum
      complex, dimension(:, :, :), allocatable :: ebeta_sum
      complex, dimension(:, :, :), allocatable :: eb_sum

      real, dimension (:, :, :), allocatable :: bR_3d
      real, dimension (:, :, :), allocatable :: bphi_3d
      real, dimension (:, :, :), allocatable :: bZ_3d

      complex, dimension(:, :, :), allocatable :: ntilda
      real, dimension(:, :, :), allocatable :: ntilda_real

      complex, dimension(:, :, :), allocatable :: xjpx_sum
      complex, dimension(:, :, :), allocatable :: xjpy_sum
      complex, dimension(:, :, :), allocatable :: xjpz_sum

      complex, dimension(:, :, :), allocatable :: xjx_sum
      complex, dimension(:, :, :), allocatable :: xjy_sum
      complex, dimension(:, :, :), allocatable :: xjz_sum
      real, dimension(:, :, :), allocatable :: xjmod_sum

      complex, dimension(:, :, :), allocatable :: xjpxe_lab_sum
      complex, dimension(:, :, :), allocatable :: xjpye_lab_sum
      complex, dimension(:, :, :), allocatable :: xjpze_lab_sum

      complex, dimension(:, :, :), allocatable :: bxwave_sum
      complex, dimension(:, :, :), allocatable :: bywave_sum
      complex, dimension(:, :, :), allocatable :: bzwave_sum

      complex, dimension(:, :, :), allocatable :: ex_sum
      complex, dimension(:, :, :), allocatable :: ey_sum
      complex, dimension(:, :, :), allocatable :: ez_sum

      complex, dimension(:, :, :), allocatable :: ntilda_e_sum

      real, dimension (:, :, :), allocatable :: redotj
      real, dimension (:, :, :), allocatable :: rho_sum

      integer::  nsum

      external jfun

      t1 = second1(dummy)

      open(unit=63, file='aorsa2d.in', status='old', form='formatted')
      rewind(63)


      write(6, *) "myid_sum = ", myid

      nphi1 = -50
      nphi2 =  50

      write(6, *) "nphi1 = ", nphi1
      write(6, *) "nphi2 = ", nphi2

*     ------------------------
*     Read namelist input data
*     ------------------------
      read (63, aorsa2din)
      close (63)

*     ------------------------
*     Write namelist input data
*     ------------------------
      write (6, aorsa2din)

*     -------------------------------------------
*     Set the remainder of the nphi_array to zero
*     -------------------------------------------
      do i = nphi_number + 1, nphimx
         nphi_array(i) = 0.0
      end do

      write(6, *) "AORSA2D: version_number = ", version_number
      write(6, *) "nphi_number = ", nphi_number
      write(6, *) "nphi_array = ", nphi_array
      write(6, *) "phase_array = ", phase_array


      open(unit=16,file='out16',status='unknown',form='formatted')
      open(unit=17,file='bharvey2_3d',status='unknown',form='formatted')


      write(16, *) "AORSA2D: version_number = ", version_number
      write(16, *) "nphi_number = ", nphi_number
      write(16, *) "nphi_array = ", nphi_array
      write(16, *) "phase_array = ", phase_array


      pi = 3.141592654
      zi = cmplx(0., 1.)
      q = 1.6e-19
      width = pi / 4.

      omgrf = 2.0 * pi * freqcy

      write(6, *) "omgrf = ", omgrf

      eps0 = 8.85e-12
      xmu0 = 1.26e-06
      clight = 1.0 / sqrt(eps0 * xmu0)
      qe = -1.6e-19
      absqe = abs(qe)
      xk0 = omgrf / clight

      write(6,* ) "qe = ", qe

      nstrpol = 1
      separ0 = 0.0000E+00
      dphase0 = 0.0000E+00
      phase0 = 0.0000E+00

      nphi  = 2048
      nphi3d = 200

      nphi3d_1quarter = nphi3d / 4 * 1
      nphi3d_half     = nphi3d / 4 * 2
      nphi3d_3quarter = nphi3d / 4 * 3


*     ----------
*     input data
*     ----------
*     ------------------------
*     Open and read file "fpm"
*     ------------------------
      open(unit=34, file='fpm', status='old', form= 'formatted')
      open(unit=166, file='bharvey_3d', status='old', form= 'formatted')

      read(34, 309) nphi1_dum, nphi2_dum

      read(34,309) nmodesx, nmodesy
      read(34,310) rwleft, rwright, ytop, ybottom

      read(34, 309) nnodex, nnodey
      read(34, 310) (capr(i), i = 1, nnodex)
      read(34, 310) (capz(j), j = 1, nnodey)

      read(34, 310) rhoplasm, rt, b0

      rhoplasm = .99
      read(34, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      read(34, 310) ((xnea(i, j), i = 1, nnodex), j = 1, nnodey)

      read (34, 309) nuper
      read (34, 309) nupar
      read (34, 309) nnoderho

      allocate( UPERP(nuper) )
      allocate( UPARA(nupar) )

      allocate( xplot(nxplot_dim) )
      allocate( yplota(nyplot_dim) )
      allocate( xplotm(nxplot_dim) )
      allocate( yplotm(nyplot_dim) )

      allocate( ealpha_sum_real(nxmx, nymx, nphimx) )
      allocate( ealpha_sum_imag(nxmx, nymx, nphimx) )

      allocate( eb_sum_real(nxmx, nymx, nphimx) )
      allocate( eb_sum_imag(nxmx, nymx, nphimx) )

      allocate( eb_sum_mod(nxmx, nymx, nphimx) )
      allocate( ealpha_sum_mod(nxmx, nymx, nphimx) )

      allocate( ebeta_sum_mod(nxmx, nymx, nphimx)  )
      allocate( mod_e_sum(nxmx, nymx, nphimx)  )

      allocate( bR_3d(nxmx, nymx, nphimx) )
      allocate( bphi_3d(nxmx, nymx, nphimx) )
      allocate( bZ_3d(nxmx, nymx, nphimx) )

      allocate( ealpha_sum(nxmx, nymx, nphimx) )
      allocate( ebeta_sum(nxmx, nymx, nphimx) )
      allocate( eb_sum(nxmx, nymx, nphimx) )

      allocate(ntilda(nxmx, nymx, nphimx) )
      allocate(ntilda_real(nxmx, nymx, nphimx) )

      allocate( xjpx_sum(nxmx, nymx, nphimx) )
      allocate( xjpy_sum(nxmx, nymx, nphimx) )
      allocate( xjpz_sum(nxmx, nymx, nphimx) )

      allocate( xjx_sum(nxmx, nymx, nphimx) )
      allocate( xjy_sum(nxmx, nymx, nphimx) )
      allocate( xjz_sum(nxmx, nymx, nphimx) )
      allocate( xjmod_sum(nxmx, nymx, nphimx) )

      allocate( xjpxe_lab_sum(nxmx, nymx, nphimx) )
      allocate( xjpye_lab_sum(nxmx, nymx, nphimx) )
      allocate( xjpze_lab_sum(nxmx, nymx, nphimx) )

      allocate( bxwave_sum(nxmx, nymx, nphimx) )
      allocate( bywave_sum(nxmx, nymx, nphimx) )
      allocate( bzwave_sum(nxmx, nymx, nphimx) )

      allocate( ex_sum(nxmx, nymx, nphimx) )
      allocate( ey_sum(nxmx, nymx, nphimx) )
      allocate( ez_sum(nxmx, nymx, nphimx) )

      allocate( ntilda_e_sum(nxmx, nymx, nphimx) )

      allocate( redotj(nxmx, nymx, nphimx)  )
      allocate( rho_sum(nxmx, nymx, nphimx) )

      allocate( frealp(nxplot_dim, nyplot_dim) )
      allocate( fimagp(nxplot_dim, nyplot_dim) )
      allocate( powerp(nxplot_dim, nyplot_dim) )
      allocate( capr_plotp(nxplot_dim, nyplot_dim) )

      allocate( freal(nxmx, nymx) )
      allocate( fimag(nxmx, nymx) )
      allocate( capr_plot(nxmx, nymx) )
      allocate( power(nxmx, nymx) )
      allocate( fluct(nxmx, nymx) )

      allocate( freal_phi(nxmx, nphimx) )
      allocate( fimag_phi(nxmx, nphimx) )
      allocate( capr_plot_phi(nxmx, nphimx) )
      allocate( power_phi(nxmx, nphimx) )

      allocate( ealpha(nxmx, nymx) )
      allocate( ebeta(nxmx, nymx) )
      allocate( eb(nxmx, nymx) )

      allocate( xjpx(nxmx, nymx) )
      allocate( xjpy(nxmx, nymx) )
      allocate( xjpz(nxmx, nymx) )

      allocate( xjx(nxmx, nymx) )
      allocate( xjy(nxmx, nymx) )
      allocate( xjz(nxmx, nymx) )

      allocate( xjpxe_lab(nxmx, nymx) )
      allocate( xjpye_lab(nxmx, nymx) )
      allocate( xjpze_lab(nxmx, nymx) )

      allocate( sR_3d(nxmx, nymx, nphimx) )
      allocate( sZ_3d(nxmx, nymx, nphimx) )
      allocate( sphi_3d(nxmx, nymx, nphimx) )

      allocate( ntilda_e(nxmx, nymx) )

      allocate( eplus(nxmx, nymx) )
      allocate( eminus(nxmx, nymx) )

      allocate (ex(nxmx, nymx) )
      allocate (ey(nxmx, nymx) )
      allocate (ez(nxmx, nymx) )
      allocate (bxwave(nxmx, nymx) )
      allocate (bywave(nxmx, nymx) )
      allocate (bzwave(nxmx, nymx) )

      allocate (bx(nxmx, nymx) )
      allocate (by(nxmx, nymx) )
      allocate (bz(nxmx, nymx) )



      if(ndisti1 .eq. 1)read (34, 3310) vc1_cgs
      if(ndisti2 .eq. 1)read (34, 3310) vc2_cgs

      read (34, 3310) UminPara, UmaxPara

      read (34, 3310) (rhon(n), n = 1, nnoderho)
      read (34, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
      read (34, 3310) (upara(i_upara), i_upara = 1, nupar)

      if(ndisti1 .eq. 1) then
         allocate( bqlavg_i1(nuper, nupar, nnoderho) )
         allocate( cqlavg_i1(nuper, nupar, nnoderho) )
         allocate( eqlavg_i1(nuper, nupar, nnoderho) )
         allocate( fqlavg_i1(nuper, nupar, nnoderho) )

         allocate( bqlsum_i1(nuper, nupar, nnoderho) )
         allocate( cqlsum_i1(nuper, nupar, nnoderho) )
         allocate( eqlsum_i1(nuper, nupar, nnoderho) )
         allocate( fqlsum_i1(nuper, nupar, nnoderho) )

         bqlavg_i1 = 0.0
         cqlavg_i1 = 0.0
         eqlavg_i1 = 0.0
         fqlavg_i1 = 0.0

         bqlsum_i1 = 0.0
         cqlsum_i1 = 0.0
         eqlsum_i1 = 0.0
         fqlsum_i1 = 0.0
      end if


      if(ndisti2 .eq. 1) then
         allocate( bqlavg_i2(nuper, nupar, nnoderho) )
         allocate( cqlavg_i2(nuper, nupar, nnoderho) )
         allocate( eqlavg_i2(nuper, nupar, nnoderho) )
         allocate( fqlavg_i2(nuper, nupar, nnoderho) )

         allocate( bqlsum_i2(nuper, nupar, nnoderho) )
         allocate( cqlsum_i2(nuper, nupar, nnoderho) )
         allocate( eqlsum_i2(nuper, nupar, nnoderho) )
         allocate( fqlsum_i2(nuper, nupar, nnoderho) )

         bqlavg_i2 = 0.0
         cqlavg_i2 = 0.0
         eqlavg_i2 = 0.0
         fqlavg_i2 = 0.0

         bqlsum_i2 = 0.0
         cqlsum_i2 = 0.0
         eqlsum_i2 = 0.0
         fqlsum_i2 = 0.0
      end if



*     --------------------------------------
*     Set up arrays for toroidal polar plots
*     --------------------------------------
      nxplot = 1000
      nyplot = 1000

      call a1mnmx_dp(capr, nxmx, nnodex, capr_min, capr_max)

      dcapr = (capr_max - capr_min) / (nnodex - 1)

      rmax = capr_max * 1.1
      rmin = -rmax



      dxplot = (rmax - rmin) / (nxplot - 1)
      dyplot = (rmax - rmin) / (nyplot - 1)


      do i = 1, nxplot
         xplot(i) = rmin + (i-1) * dxplot
         xplotm(i) = -xplot(i)
      end do

c      write(6, *)"xplotm(i) = "
c      write (6, 310) (xplotm(i),  i = 1, nxplot)

      do j = 1, nyplot
         yplota(j) = rmin + (j-1) * dyplot
         yplotm(j) = -yplota(j)
      end do



      r0 = rt

      phase_array = phase_array * pi / 180.
      phi0 = pi
      xi0 = 57.385758

      do  n = 1, nstrap
         amplt(n) = 1.0
         write(6, 1053)amplt(n)
      end do

      phimin = 0.0
      phimax = 2.0 * pi
      dphi   = (phimax - phimin) / nphi
      dphi3d = (phimax - phimin) / nphi3d

      zmin = rant * phimin
      zmax = rant * phimax
      xltcm = xlt * 100.0


*     --------------------------------------------------
*     Antenna current on fine mesh of nphi = 2048 points
*     --------------------------------------------------
      do k = 1, nphi
         phi(k) = (k-1) * dphi
         zk(k) = rant * phi(k)

         curent(k) = jfun(zk(k))
      end do


*     ------------------------------------------------------
*     Antenna current on course mesh of nphi3d (~200) points
*     ------------------------------------------------------
      do k = 1, nphi3d
         phi3d(k) =         0. + (k - 1) * dphi3d
         phi3d_shift(k) = - pi + (k - 1) * dphi3d

         zk3d(k)  = rant * phi3d(k)
         cur3d(k) = jfun(zk3d(k))
      end do



*     ------------------------------------------
*     do fft on fine mesh giving nphi=2048 modes
*     ------------------------------------------
      call fftn2(curent, work, nphi, nphmax, ntmin, ntmax, cn)


*     ------------------------------------------------------
*     reconstruct current on fine mesh of nphi = 2048 points
*     ------------------------------------------------------

      call fftin2(cursum, work, nphi, nphmax, ntmin, ntmax, cn)

      nt2 = nphi / 2
      nt1 = nt2 - nphi + 1

      write(16, 162)
      write(16, 1163)
      write(16, 163)

      do k = 1, nphi

         csum = 0.0
         do nphi = nt1, nt2
            csum = csum + cn(nphi) * cexp(zi * nphi * phi(k))
         end do

         write(16,2001)k, phi(k), curent(k), cursum(k), csum
c        write( 6,2001)k, phi(k), curent(k), cursum(k), csum
      end do


*     ---------------------------------------------------------
*     discard all modes with nphi .lt. nphi1 or nphi .gt. nphi2
*     ---------------------------------------------------------
      do nphi = ntmin, ntmax
         if (nphi .lt. nphi1 .or. nphi .gt. nphi2) cn(nphi) = 0.0
      end do


*     ---------------------------------------------------
*     reconstruct current on mesh of nphi3d (~200) points
*     ---------------------------------------------------
      cursum = 0.0

c      call fftin2(cursum, work, nphi3d, nphmax, ntmin, ntmax, cn)

      nt2 = nphi3d / 2
      nt1 = nt2 - nphi3d + 1

      write(16, 162)
      write(16, 1163)
      write(16, 163)

      do k = 1, nphi3d

         csum = 0.0
         do nphi = nt1, nt2
            csum = csum + cn(nphi) * cexp(zi * nphi * phi3d(k))
         end do

c         write(16,2001)k, phi3d(k), cur3d(k), cursum(k), csum
         write(16,2001)k, phi3d(k), cur3d(k), csum
         write(6, 2001)k, phi3d(k), cur3d(k), csum
      end do



      ealpha_sum = 0.0
      ebeta_sum = 0.0
      eb_sum = 0.0

      ntilda = 0.0
      ntilda_real = 0.0

      xjpx_sum = 0.0
      xjpy_sum = 0.0
      xjpz_sum = 0.0

      xjx_sum = 0.0
      xjy_sum = 0.0
      xjz_sum = 0.0

      bxwave_sum = 0.0
      bywave_sum = 0.0
      bzwave_sum = 0.0

      ex_sum = 0.0
      ey_sum = 0.0
      ez_sum = 0.0

      xjpxe_lab_sum = 0.0
      xjpye_lab_sum = 0.0
      xjpze_lab_sum = 0.0

      ntilda_e_sum = 0.0

      pabs_sum = 0.0
      jdriven_sum = 0.0

      wdotesum = 0.0
      wdot1sum = 0.0
      wdot2sum = 0.0
      wdot3sum = 0.0
      wdot4sum = 0.0
      wdot5sum = 0.0
      wdot6sum = 0.0

      wdotesum_ql = 0.0
      wdot1sum_ql = 0.0
      wdot2sum_ql = 0.0
      wdot3sum_ql = 0.0
      wdot4sum_ql = 0.0
      wdot5sum_ql = 0.0

      redotjesum = 0.0
      redotj1sum = 0.0
      redotj2sum = 0.0
      redotj3sum = 0.0
      redotj4sum = 0.0
      redotj5sum = 0.0
      redotj6sum = 0.0

      xjprl_sum = 0.0


      write(6, *)
      write(16,*)
      write(6, *) "       nphi       nt     xnphi      cnmod2      pabs
     &    jdriven         spa"
      write(16,*) "       nphi       nt     xnphi      cnmod2      pabs
     &    jdriven         spa"
      write(6, *)
      write(16,*)

      write(17, 310) xi0, xlt

      nt = 0
      do nphi = nphi1, nphi2
         nt = nt +1
         cnmod2(nt) = conjg(cn(nphi)) * cn(nphi) * (xi0 / xlt)**2
c        if(nphi .eq. 16) cnmod2(nt) = 0.0
         xnphi(nt) = real(nphi)
         write(17, 123) nt, nphi, cn(nphi), cnmod2(nt)
      end do

      close (17)

c      read(166, 309) nnodex, nnodey
c      read(166, 310) (capr(i), i = 1, nnodex)
c      read(166, 310) (capz(j), j = 1, nnodey)
c      read(166, 310) psilim
c      read(166, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)

c      read(166, 309) nphi_number
c      read(166, 309) (nphi_array(n), n = 1, nphi_number)

*     ----------------
*     Sum over nphi's
*     ----------------
*     -------------------------------------------------------------
*     Leave all modes out of the sum except the pre-determined ones
*     -------------------------------------------------------------
      do 9000 nsum = 1, nphi_number

         nphi = nphi_array(nsum)
         nt =  nphi - nphi1 + 1

c         read(166,310) ((eplus(i, j), i = 1, nnodex), j = 1, nnodey)
c         read(166,310) ((eminus(i, j),  i = 1, nnodex), j = 1, nnodey)
c         read(166,310) ((eb(i, j),     i = 1, nnodex), j = 1, nnodey)

c        read(166,310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
c         read(166,310) ((ey(i, j),  i = 1, nnodex), j = 1, nnodey)
c         read(166,310) ((ez(i, j),     i = 1, nnodex), j = 1, nnodey)

c         read(166,310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
c         read(166,310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
c         read(166,310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)



         read(34, 309) nphi_fpm

         read(34, 310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34, 310) ((ebeta(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34, 310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)

         read(34,310) ((xjpx(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34,310) ((xjpy(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34,310) ((xjpz(i, j), i = 1, nnodex), j = 1, nnodey)

         read(34,310) ((xjx(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34,310) ((xjy(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34,310) ((xjz(i, j), i = 1, nnodex), j = 1, nnodey)

         read(34,310) ((xjpxe_lab(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34,310) ((xjpye_lab(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34,310) ((xjpze_lab(i, j), i = 1, nnodex), j = 1, nnodey)

         read(34, 310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34, 310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34, 310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

         read(34, 310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34, 310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
         read(34, 310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)

         read(34,310) ((ntilda_e(i, j), i = 1, nnodex), j = 1, nnodey)

         read(34,310) ptot, pcrto2, pcito2, xjtot, spa_cold

         read(34, 309) nnoderho
         read(34, 310) (rhon(n), n = 1, nnoderho)

         read(34, 310) (redotjeavg(n), n = 1, nnoderho)
         read(34, 310) (redotj1avg(n), n = 1, nnoderho)
         read(34, 310) (redotj2avg(n), n = 1, nnoderho)
         read(34, 310) (redotj3avg(n), n = 1, nnoderho)
         read(34, 310) (redotj4avg(n), n = 1, nnoderho)
         read(34, 310) (redotj5avg(n), n = 1, nnoderho)
         read(34, 310) (redotj6avg(n), n = 1, nnoderho)

         read(34, 310) (wdoteavg(n), n = 1, nnoderho)
         read(34, 310) (wdot1avg(n), n = 1, nnoderho)
         read(34, 310) (wdot2avg(n), n = 1, nnoderho)
         read(34, 310) (wdot3avg(n), n = 1, nnoderho)
         read(34, 310) (wdot4avg(n), n = 1, nnoderho)
         read(34, 310) (wdot5avg(n), n = 1, nnoderho)
         read(34, 310) (wdot6avg(n), n = 1, nnoderho)

         read(34, 310) (wdote_ql(n),  n = 1, nnoderho)
         read(34, 310) (wdoti1_ql(n), n = 1, nnoderho)
         read(34, 310) (wdoti2_ql(n), n = 1, nnoderho)
         read(34, 310) (wdoti3_ql(n), n = 1, nnoderho)
         read(34, 310) (wdoti4_ql(n), n = 1, nnoderho)
         read(34, 310) (wdoti5_ql(n), n = 1, nnoderho)

         read(34, 310) (xjprlavg(n),  n = 1, nnoderho)



         if(ndisti1   .eq. 1) then
            read (34, 3310) (((bqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            read (34, 3310) (((cqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            read (34, 3310) (((eqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            read (34, 3310) (((fqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            read (34, 3310) xmi1
         end if

         if(ndisti2   .eq. 1) then
            read (34, 3310) (((bqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            read (34, 3310) (((cqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            read (34, 3310) (((eqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            read (34, 3310) (((fqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            read (34, 3310) xmi1
         end if



               pabs(nt)    = ptot
               jdriven(nt) = xjtot
               spa(nt) = spa_cold

               if (cnmod2(nt) .eq. 0.0) then
                  pabs(nt)    = 0.0
                  jdriven(nt) = 0.0
                  spa(nt) = 0.0
               end if


               pabs_weight    = pabs(nt)    * cnmod2(nt) * 2. * pi * r0
               jdriven_weight = jdriven(nt) * cnmod2(nt)

               pabs_sum    = pabs_sum    + pabs_weight
               jdriven_sum = jdriven_sum + jdriven_weight

               write(6,123)nphi, nt, xnphi(nt), cnmod2(nt), pabs(nt),
     &                                       jdriven(nt), spa(nt)
               write(16,123)nphi, nt, xnphi(nt), cnmod2(nt), pabs(nt),
     &                                       jdriven(nt), spa(nt)

               do n = 1, nnoderho
                  wdotesum(n) = wdotesum(n) + wdoteavg(n) * cnmod2(nt)
                  wdot1sum(n) = wdot1sum(n) + wdot1avg(n) * cnmod2(nt)
                  wdot2sum(n) = wdot2sum(n) + wdot2avg(n) * cnmod2(nt)
                  wdot3sum(n) = wdot3sum(n) + wdot3avg(n) * cnmod2(nt)
                  wdot4sum(n) = wdot4sum(n) + wdot4avg(n) * cnmod2(nt)
                  wdot5sum(n) = wdot5sum(n) + wdot5avg(n) * cnmod2(nt)
                  wdot6sum(n) = wdot6sum(n) + wdot6avg(n) * cnmod2(nt)

                  wdotesum_ql(n) = wdotesum_ql(n)+wdote_ql(n)*cnmod2(nt)
                  wdot1sum_ql(n)= wdot1sum_ql(n)+wdoti1_ql(n)*cnmod2(nt)
                  wdot2sum_ql(n)= wdot2sum_ql(n)+wdoti2_ql(n)*cnmod2(nt)
                  wdot3sum_ql(n)= wdot3sum_ql(n)+wdoti3_ql(n)*cnmod2(nt)
                  wdot4sum_ql(n)= wdot4sum_ql(n)+wdoti4_ql(n)*cnmod2(nt)
                  wdot5sum_ql(n)= wdot5sum_ql(n)+wdoti5_ql(n)*cnmod2(nt)

                  redotjesum(n) = redotjesum(n)+redotjeavg(n)*cnmod2(nt)
                  redotj1sum(n) = redotj1sum(n)+redotj1avg(n)*cnmod2(nt)
                  redotj2sum(n) = redotj2sum(n)+redotj2avg(n)*cnmod2(nt)
                  redotj3sum(n) = redotj3sum(n)+redotj3avg(n)*cnmod2(nt)
                  redotj4sum(n) = redotj4sum(n)+redotj4avg(n)*cnmod2(nt)
                  redotj5sum(n) = redotj5sum(n)+redotj5avg(n)*cnmod2(nt)
                  redotj6sum(n) = redotj6sum(n)+redotj6avg(n)*cnmod2(nt)

                  xjprl_sum(n) = xjprl_sum(n) + xjprlavg(n) * cnmod2(nt)


                  if(ndisti1  .eq. 1) then
                     do i_uperp =  1, nuper
                        do i_upara = 1, nupar
                           bqlsum_i1(i_uperp, i_upara, n) =
     &                        bqlsum_i1(i_uperp, i_upara, n)
     &                      + bqlavg_i1(i_uperp, i_upara, n)* cnmod2(nt)

                           cqlsum_i1(i_uperp, i_upara, n) =
     &                        cqlsum_i1(i_uperp, i_upara, n)
     &                      + cqlavg_i1(i_uperp, i_upara, n)* cnmod2(nt)

                           eqlsum_i1(i_uperp, i_upara, n) =
     &                        eqlsum_i1(i_uperp, i_upara ,n)
     &                      + eqlavg_i1(i_uperp, i_upara, n)* cnmod2(nt)

                           fqlsum_i1(i_uperp, i_upara, n) =
     &                        fqlsum_i1(i_uperp, i_upara, n)
     &                      + fqlavg_i1(i_uperp, i_upara, n)* cnmod2(nt)
                        end do
                     end do
                  end if

                  if(ndisti2  .eq. 1) then
                     do i_uperp =  1, nuper
                        do i_upara = 1, nupar
                           bqlsum_i2(i_uperp, i_upara, n) =
     &                     bqlsum_i2(i_uperp, i_upara, n)
     &                   + bqlavg_i2(i_uperp, i_upara, n) * cnmod2(nt)

                           cqlsum_i2(i_uperp, i_upara, n) =
     &                     cqlsum_i2(i_uperp, i_upara, n)
     &                   + cqlavg_i2(i_uperp, i_upara, n) * cnmod2(nt)

                           eqlsum_i2(i_uperp, i_upara, n) =
     &                     eqlsum_i2(i_uperp, i_upara ,n)
     &                   + eqlavg_i2(i_uperp, i_upara, n) * cnmod2(nt)

                           fqlsum_i2(i_uperp, i_upara, n) =
     &                     fqlsum_i2(i_uperp, i_upara, n)
     &                   + fqlavg_i2(i_uperp, i_upara, n) * cnmod2(nt)
                        end do
                     end do
                  end if

               end do


               do k = 1, nphi3d

                  cexpkz = xi0 / xlt * cn(nphi) * cexp(zi*nphi*phi3d(k))

                  do i = 1, nnodex
                     do j = 1, nnodey

                        ealpha_sum(i,j,k)=ealpha_sum(i,j,k)
     &                                              +ealpha(i,j)*cexpkz
                        ebeta_sum(i,j,k) =ebeta_sum(i,j,k)
     &                                              +ebeta(i,j) *cexpkz
                        eb_sum(i,j,k)    =eb_sum(i,j,k)
     &                                               +  eb(i,j) *cexpkz

                        xjpx_sum(i,j,k)=xjpx_sum(i,j,k)+xjpx(i,j)*cexpkz
                        xjpy_sum(i,j,k)=xjpy_sum(i,j,k)+xjpy(i,j)*cexpkz
                        xjpz_sum(i,j,k)=xjpz_sum(i,j,k)+xjpz(i,j)*cexpkz

                        xjx_sum(i,j,k)=xjx_sum(i,j,k)+xjx(i,j)*cexpkz
                        xjy_sum(i,j,k)=xjy_sum(i,j,k)+xjy(i,j)*cexpkz
                        xjz_sum(i,j,k)=xjz_sum(i,j,k)+xjz(i,j)*cexpkz

                        bxwave_sum(i,j,k) = bxwave_sum(i,j,k)
     &                                            + bxwave(i,j) * cexpkz
                        bywave_sum(i,j,k) = bywave_sum(i,j,k)
     &                                            + bywave(i,j) * cexpkz
                        bzwave_sum(i,j,k) = bzwave_sum(i,j,k)
     &                                            + bzwave(i,j) * cexpkz

                        ex_sum(i,j,k) = ex_sum(i,j,k) + ex(i,j) * cexpkz
                        ey_sum(i,j,k) = ey_sum(i,j,k) + ey(i,j) * cexpkz
                        ez_sum(i,j,k) = ez_sum(i,j,k) + ez(i,j) * cexpkz

                        xjpxe_lab_sum(i,j,k)=xjpxe_lab_sum(i,j,k)
     &                                  + xjpxe_lab(i,j) * cexpkz
                        xjpye_lab_sum(i,j,k)=xjpye_lab_sum(i,j,k)
     &                                   +xjpye_lab(i,j) * cexpkz
                        xjpze_lab_sum(i,j,k)=xjpze_lab_sum(i,j,k)
     &                                   +xjpze_lab(i,j) * cexpkz

                        ntilda_e_sum(i,j,k) = ntilda_e_sum(i,j,k)
     &                                     + ntilda_e(i,j) * cexpkz

                        rho_sum(i,j,k)  = rho(i,j)


                     end do
                  end do

               end do

 9000 continue
*     ----------------------
*     end of sum over nphi's
*     ----------------------

      open(unit = 638, file ='fpm2', status ='old', form ='formatted')
      read(638, 309) nnodex, nnodey
      read(638, 310) ((bx(i, j), i = 1, nnodex), j = 1, nnodey)
      read(638, 310) ((by(i, j), i = 1, nnodex), j = 1, nnodey)
      read(638, 310) ((bz(i, j), i = 1, nnodex), j = 1, nnodey)


      close (638)
      close (34)
      close (166)


*     -------------------------------------------
*     scale fields, currents, and powers to prfin
*     -------------------------------------------
      pscale = 1.0
      if(pabs_sum .ne. 0.0)pscale = prfin / pabs_sum

      write(6, *) "pscale = ", pscale
      write(16, *)"pscale = ", pscale



      do k = 1, nphi3d
         do i = 1, nnodex
            do j = 1, nnodey

               bR_3d(i,j,k) = bx(i,j)
               bZ_3d(i,j,k) = by(i,j)
               bphi_3d(i,j,k) = signbz * bz(i,j)

*     -----------------------------------------------------
*     Poynthing flux in cylindrical coordinates (R, phi, Z)
*     -----------------------------------------------------
               sR_3d(i,j,k) = 1. / (2. * xmu0) * real(
     &                         conjg(ez_sum(i,j,k)) * bywave_sum(i,j,k)
     &                       - conjg(ey_sum(i,j,k)) * bzwave_sum(i,j,k))
               sphi_3d(i,j,k) = signbz * 1. / (2. * xmu0) * real(
     &                         conjg(ey_sum(i,j,k)) * bxwave_sum(i,j,k)
     &                       - conjg(ex_sum(i,j,k)) * bywave_sum(i,j,k))
               sZ_3d(i,j,k) = 1. / (2. * xmu0) * real(
     &                         conjg(ex_sum(i,j,k)) * bzwave_sum(i,j,k)
     &                       - conjg(ez_sum(i,j,k)) * bxwave_sum(i,j,k))



*     -----------------------------------------------------
*     Poynting flux in Cartesian (x, y, z) with phi > - phi
*     -----------------------------------------------------
c               Sx = 1. / (2. * xmu0) * real(
c     .                         conjg(ey_sum(i,j,k)) * bzwave_sum(i,j,k)
c     .                       - conjg(ez_sum(i,j,k)) * bywave_sum(i,j,k))
c               Sy = 1. / (2. * xmu0) * real(
c     .                         conjg(ez_sum(i,j,k)) * bxwave_sum(i,j,k)
c     .                       - conjg(ex_sum(i,j,k)) * bzwave_sum(i,j,k))
c               Sz = 1. / (2. * xmu0) * real(
c     .                         conjg(ex_sum(i,j,k)) * bywave_sum(i,j,k)
c     .                       - conjg(ey_sum(i,j,k)) * bxwave_sum(i,j,k))
c               sR_3d(i,j,k)   = - Sx
c               sZ_3d(i,j,k)   = - Sy
c               sphi_3d(i,j,k) = - Sz


               ealpha_sum(i,j,k) =  ealpha_sum(i,j,k) * sqrt(pscale)
               ebeta_sum(i,j,k)  =  ebeta_sum(i,j,k)  * sqrt(pscale)
               eb_sum(i,j,k)     =  eb_sum(i,j,k)     * sqrt(pscale)

               ex_sum(i,j,k) =  ex_sum(i,j,k) * sqrt(pscale)
               ey_sum(i,j,k) =  ey_sum(i,j,k) * sqrt(pscale)
               ez_sum(i,j,k) =  ez_sum(i,j,k) * sqrt(pscale)

               xjpx_sum(i,j,k) = xjpx_sum(i,j,k) * sqrt(pscale)
               xjpy_sum(i,j,k) = xjpy_sum(i,j,k) * sqrt(pscale)
               xjpz_sum(i,j,k) = xjpz_sum(i,j,k) * sqrt(pscale)

               xjx_sum(i,j,k) = xjx_sum(i,j,k) * sqrt(pscale)
               xjy_sum(i,j,k) = xjy_sum(i,j,k) * sqrt(pscale)
               xjz_sum(i,j,k) = xjz_sum(i,j,k) * sqrt(pscale)

               xjmod_sum(i,j,k) = sqrt(
     &               conjg(xjx_sum(i,j,k)) * xjx_sum(i,j,k)
     &             + conjg(xjy_sum(i,j,k)) * xjy_sum(i,j,k)
     &             + conjg(xjz_sum(i,j,k)) * xjz_sum(i,j,k) )

               bxwave_sum(i,j,k) = bxwave_sum(i,j,k) * sqrt(pscale)
               bywave_sum(i,j,k) = bywave_sum(i,j,k) * sqrt(pscale)
               bzwave_sum(i,j,k) = bzwave_sum(i,j,k) * sqrt(pscale)

               xjpxe_lab_sum(i,j,k) = xjpxe_lab_sum(i,j,k) *sqrt(pscale)
               xjpye_lab_sum(i,j,k) = xjpye_lab_sum(i,j,k) *sqrt(pscale)
               xjpze_lab_sum(i,j,k) = xjpze_lab_sum(i,j,k) *sqrt(pscale)

               ntilda_e_sum(i,j,k) = ntilda_e_sum(i,j,k) *sqrt(pscale)

               ntilda(i,j,k) = ntilda_e_sum(i,j,k)
               ntilda_real(i,j,k) = real(ntilda(i,j,k)) / 1.0e+20



            end do
         end do
      end do


      pabs_sum    = pabs_sum    * pscale
      jdriven_sum = jdriven_sum * pscale

      do n = 1, nnoderho
         wdotesum(n) = wdotesum(n) * pscale
         wdot1sum(n) = wdot1sum(n) * pscale
         wdot2sum(n) = wdot2sum(n) * pscale
         wdot3sum(n) = wdot3sum(n) * pscale
         wdot4sum(n) = wdot4sum(n) * pscale
         wdot5sum(n) = wdot5sum(n) * pscale
         wdot6sum(n) = wdot6sum(n) * pscale

         wdotesum_ql(n) = wdotesum_ql(n) * pscale
         wdot1sum_ql(n) = wdot1sum_ql(n) * pscale
         wdot2sum_ql(n) = wdot2sum_ql(n) * pscale
         wdot3sum_ql(n) = wdot3sum_ql(n) * pscale
         wdot4sum_ql(n) = wdot4sum_ql(n) * pscale
         wdot5sum_ql(n) = wdot5sum_ql(n) * pscale

         redotjesum(n) = redotjesum(n) * pscale
         redotj1sum(n) = redotj1sum(n) * pscale
         redotj2sum(n) = redotj2sum(n) * pscale
         redotj3sum(n) = redotj3sum(n) * pscale
         redotj4sum(n) = redotj4sum(n) * pscale
         redotj5sum(n) = redotj5sum(n) * pscale
         redotj6sum(n) = redotj6sum(n) * pscale

         xjprl_sum(n)  = xjprl_sum(n)  * pscale

         if(ndisti1  .eq. 1) then
            do i_uperp =  1, nuper
               do i_upara = 1, nupar
                 bqlsum_i1(i_uperp, i_upara, n) =
     &                     bqlsum_i1(i_uperp, i_upara, n) * pscale
                 cqlsum_i1(i_uperp, i_upara, n) =
     &                     cqlsum_i1(i_uperp, i_upara, n) * pscale
                 eqlsum_i1(i_uperp, i_upara, n) =
     &                     eqlsum_i1(i_uperp, i_upara ,n) * pscale
                 fqlsum_i1(i_uperp, i_upara, n) =
     &                     fqlsum_i1(i_uperp, i_upara, n) * pscale
               end do
            end do
         end if

         if(ndisti2  .eq. 1) then
            do i_uperp =  1, nuper
               do i_upara = 1, nupar
                 bqlsum_i2(i_uperp, i_upara, n) =
     &                     bqlsum_i2(i_uperp, i_upara, n) * pscale
                 cqlsum_i2(i_uperp, i_upara, n) =
     &                     cqlsum_i2(i_uperp, i_upara, n) * pscale
                 eqlsum_i2(i_uperp, i_upara, n) =
     &                     eqlsum_i2(i_uperp, i_upara ,n) * pscale
                 fqlsum_i2(i_uperp, i_upara, n) =
     &                     fqlsum_i2(i_uperp, i_upara, n) * pscale
               end do
            end do
         end if

      end do


*     ---------------------------------------------------
*     Calculate the differential volume elements, dvol(n):
*     ---------------------------------------------------
      dx = capr(2) - capr(1)
      dy = capz(2) - capz(1)
      drho = rhon(2) - rhon(1)

      dvol  = 0.0
      darea = 0.0

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1
            if(n .le. nnoderho .and. n .ge. 1)then
               dvol(n)  =  dvol(n) + dx * dy * 2.0 * pi * capr(i)
               darea(n) =  darea(n)+ dx * dy * capr(i)/ rt
            end if
         end do
      end do

      dvol(1) = 0.0



*     --------------------------------------------
*     Integrate flux averaged quantities over rhon:
*     --------------------------------------------

      call rhograte(rhon, redotjesum, 1, nnoderho, redotje_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj1sum, 1, nnoderho, redotji1_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj2sum, 1, nnoderho, redotji2_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj3sum, 1, nnoderho, redotji3_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj4sum, 1, nnoderho, redotji4_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj5sum, 1, nnoderho, redotji5_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj6sum, 1, nnoderho, redotji6_int,
     &   nrhomax, dvol)

      pedotjt = redotje_int(nnoderho)  + redotji1_int(nnoderho)
     &        + redotji2_int(nnoderho) + redotji3_int(nnoderho)
     &        + redotji4_int(nnoderho) + redotji5_int(nnoderho)
     &        + redotji6_int(nnoderho)


*     ---------------------------------------------
*     Calculate species fractions from Estar dot J:
*     ---------------------------------------------
      pcedotje = redotje_int(nnoderho)  / pedotjt * 100.
      pcedotj1 = redotji1_int(nnoderho) / pedotjt * 100.
      pcedotj2 = redotji2_int(nnoderho) / pedotjt * 100.
      pcedotj3 = redotji3_int(nnoderho) / pedotjt * 100.
      pcedotj4 = redotji4_int(nnoderho) / pedotjt * 100.
      pcedotj5 = redotji5_int(nnoderho) / pedotjt * 100.
      pcedotj6 = redotji6_int(nnoderho) / pedotjt * 100.

      pcedotjt = pcedotje + pcedotj1 + pcedotj2 + pcedotj3
     &         + pcedotj4 + pcedotj5 + pcedotj6

 1109 format(
     &   3x, "      power absorbed by electrons = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "  power absorbed by majority ions = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "  power absorbed by minority ions = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 3rd ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 4th ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 5th ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 6th ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "             total power absorbed = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       ")



71667 format(' Species absorption from Estar dot J:')

         write(16, 169)
         write(16, 71667)
         write(16, 1109) redotje_int(nnoderho),  pcedotje,
     &                   redotji1_int(nnoderho), pcedotj1,
     &                   redotji2_int(nnoderho), pcedotj2,
     &                   redotji3_int(nnoderho), pcedotj3,
     &                   redotji4_int(nnoderho), pcedotj4,
     &                   redotji5_int(nnoderho), pcedotj5,
     &                   redotji6_int(nnoderho), pcedotj6,
     &                   pedotjt, pcedotjt

         write(6, 169)
         write(6, 71667)
         write(6, 1109) redotje_int(nnoderho),  pcedotje,
     &                   redotji1_int(nnoderho), pcedotj1,
     &                   redotji2_int(nnoderho), pcedotj2,
     &                   redotji3_int(nnoderho), pcedotj3,
     &                   redotji4_int(nnoderho), pcedotj4,
     &                   redotji5_int(nnoderho), pcedotj5,
     &                   redotji6_int(nnoderho), pcedotj6,
     &                   pedotjt, pcedotjt



      call rhograte(rhon, wdotesum, 1, nnoderho, wdote_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot1sum, 1, nnoderho, wdoti1_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot2sum, 1, nnoderho, wdoti2_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot3sum, 1, nnoderho, wdoti3_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot4sum, 1, nnoderho, wdoti4_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot5sum, 1, nnoderho, wdoti5_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot6sum, 1, nnoderho, wdoti6_int,
     &   nrhomax, dvol)


      call rhograte(rhon, wdotesum_ql, 1, nnoderho, wdote_int_ql,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot1sum_ql, 1, nnoderho, wdoti1_int_ql,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot2sum_ql, 1, nnoderho, wdoti2_int_ql,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot3sum_ql, 1, nnoderho, wdoti3_int_ql,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot4sum_ql, 1, nnoderho, wdoti4_int_ql,
     &   nrhomax, dvol)
      call rhograte(rhon, wdot5sum_ql, 1, nnoderho, wdoti5_int_ql,
     &   nrhomax, dvol)



      call rhograte(rhon, xjprl_sum, 1, nnoderho, xjprl_int,
     &   nrhomax, darea)

      pt = wdote_int(nnoderho)  + wdoti1_int(nnoderho)
     &   + wdoti2_int(nnoderho) + wdoti3_int(nnoderho)
     &   + wdoti4_int(nnoderho) + wdoti5_int(nnoderho)
     &   + wdoti6_int(nnoderho)


c--   Calculate species fractions from Wdot:
      if (pt .ne. 0.0) then
         pcte =   wdote_int(nnoderho)  / pt * 100.
         pcti1 =  wdoti1_int(nnoderho) / pt * 100.
         pcti2 =  wdoti2_int(nnoderho) / pt * 100.
         pcti3 =  wdoti3_int(nnoderho) / pt * 100.
         pcti4 =  wdoti4_int(nnoderho) / pt * 100.
         pcti5 =  wdoti5_int(nnoderho) / pt * 100.
         pcti6 =  wdoti6_int(nnoderho) / pt * 100.
      end if

      pctt = pcte + pcti1 + pcti2 + pcti3
     &            + pcti4 + pcti5 + pcti6

      write(16, 169)
      write(16, 81667)
      write(16, 71109) wdote_int(nnoderho),  pcte,
     &                 wdoti1_int(nnoderho), pcti1,
     &                 wdoti2_int(nnoderho), pcti2,
     &                 wdoti3_int(nnoderho), pcti3,
     &                 wdoti4_int(nnoderho), pcti4,
     &                 wdoti5_int(nnoderho), pcti5,
     &                 wdoti6_int(nnoderho), pcti6,
     &                 pt,  pctt

      write(6, 169)
      write(6, 81667)
      write(6, 71109) wdote_int(nnoderho),  pcte,
     &                 wdoti1_int(nnoderho), pcti1,
     &                 wdoti2_int(nnoderho), pcti2,
     &                 wdoti3_int(nnoderho), pcti3,
     &                 wdoti4_int(nnoderho), pcti4,
     &                 wdoti5_int(nnoderho), pcti5,
     &                 wdoti6_int(nnoderho), pcti6,
     &                 pt,  pctt




      pt = wdote_int_ql(nnoderho)  + wdoti1_int_ql(nnoderho)
     &   + wdoti2_int_ql(nnoderho) + wdoti3_int_ql(nnoderho)
     &   + wdoti4_int_ql(nnoderho) + wdoti5_int_ql(nnoderho)


c--   Calculate species fractions from Wdot:
      if (pt .ne. 0.0) then
         pcte =   wdote_int_ql(nnoderho)  / pt * 100.
         pcti1 =  wdoti1_int_ql(nnoderho) / pt * 100.
         pcti2 =  wdoti2_int_ql(nnoderho) / pt * 100.
         pcti3 =  wdoti3_int_ql(nnoderho) / pt * 100.
         pcti4 =  wdoti4_int_ql(nnoderho) / pt * 100.
         pcti5 =  wdoti5_int_ql(nnoderho) / pt * 100.

      end if

      pctt = pcte + pcti1 + pcti2 + pcti3
     &            + pcti4 + pcti5


      write(16, 169)
      write(16, 81668)
      write(16, 71108) wdote_int_ql(nnoderho),  pcte,
     &                 wdoti1_int_ql(nnoderho), pcti1,
     &                 wdoti2_int_ql(nnoderho), pcti2,
     &                 wdoti3_int_ql(nnoderho), pcti3,
     &                 wdoti4_int_ql(nnoderho), pcti4,
     &                 wdoti5_int_ql(nnoderho), pcti5,
     &                 pt,  pctt

      write(6, 169)
      write(6, 81668)
      write(6, 71108) wdote_int_ql(nnoderho),  pcte,
     &                 wdoti1_int_ql(nnoderho), pcti1,
     &                 wdoti2_int_ql(nnoderho), pcti2,
     &                 wdoti3_int_ql(nnoderho), pcti3,
     &                 wdoti4_int_ql(nnoderho), pcti4,
     &                 wdoti5_int_ql(nnoderho), pcti5,
     &                 pt,  pctt



  169 format(" ")
81667 format(' Species absorption from Wdot:')
81668 format(' Species absorption from quasilinear Wdot:')
71109 format(
     &   3x, "      power absorbed by electrons = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "  power absorbed by majority ions = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "  power absorbed by minority ions = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 3rd ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 4th ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 5th ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 6th ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "             total power absorbed = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       ")

71108 format(
     &   3x, "      power absorbed by electrons = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "  power absorbed by majority ions = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "  power absorbed by minority ions = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 3rd ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 4th ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "power absorbed by 5th ion species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "             total power absorbed = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       ")


      write(6,  1217) xjprl_int(nnoderho)
      write(16, 1217) xjprl_int(nnoderho)
 1217 format(3x,'total driven current = ',1p,e12.4,' Amps')


*     -------------------------------------------
*     Calculate density fluctuation from 3D div J
*     -------------------------------------------
      do k = 1, nphi3d
         do i = 2, nnodex - 1
            do j = 2, nnodey - 1

c               xnea(i, j) = 1.277E+20

               if(k .ne. nphi3d) kp1 = k + 1
               if(k .eq. nphi3d) kp1 = 1
               if(k .ne. 1) km1 = k - 1
               if(k .eq. 1) km1 = nphi3d

               div_x = (capr(i+1) * xjpxe_lab_sum(i+1,j,k)
     &                - capr(i-1) * xjpxe_lab_sum(i-1,j,k))
     &                / dx / capr(i)
               div_y = (xjpye_lab_sum(i,j+1,k) - xjpye_lab_sum(i,j-1,k))
     &                / dy
               div_z = (xjpze_lab_sum(i,j,kp1) - xjpze_lab_sum(i,j,km1))
     &                / dphi /capr(i)

               div_j = div_x + div_y + div_z

c              ntilda(i,j,k) = - zi / (omgrf * qe) * div_j
c              ntilda_real(i,j,k) = real(ntilda(i,j,k))

            end do
         end do
      end do

      ntilda_max = maxval(ntilda_real)
      write(6, *) "ntilda_max = ", ntilda_max
      write(16, *) "ntilda_max = ", ntilda_max

*     ---------------------------
*     Calculate power from edotj:
*     ---------------------------

      do k = 1, nphi3d
         do i = 1, nnodex
            do j = 1, nnodey
               redotj(i,j,k)= 0.5 * real(
     &                      conjg(ealpha_sum(i,j,k)) * xjpx_sum(i,j,k)
     &                    + conjg(ebeta_sum(i,j,k))  * xjpy_sum(i,j,k)
     &                    + conjg(eb_sum(i,j,k))     * xjpz_sum(i,j,k) )
            end do
         end do
      end do

*     --------------------------------
*     Integrate edotj over 3-D volume:
*     --------------------------------
      call sgrate_3d(capr, capz, phi3d, redotj,
     &               1, nnodex, 1, nnodey, 1, nphi3d,
     &               ptotal, capr, nxmx, nymx, nphimx)

      call sgrate_3d_edge(capr, capz, phi3d, redotj,
     &               1, nnodex, 1, nnodey, 1, nphi3d,
     &               ptotal_edge, capr, nxmx, nymx, nphimx, rho)

      call sgrate_3d_core(capr, capz, phi3d, redotj,
     &               1, nnodex, 1, nnodey, 1, nphi3d,
     &               ptotal_core, capr, nxmx, nymx, nphimx, rho)

      fedge = ptotal_edge / ptotal * 100.0
      fcore = ptotal_core / ptotal * 100.0


      write(6, *) "ptotal = ", ptotal, "Watts"

      write(6, *) "jdriven = ", jdriven_sum, "Amps"

      write(16, *) "ptotal = ", ptotal, "Watts"
      write(16, *) "jdriven = ", jdriven_sum, "Amps"


*     ----------------------------------------
*     Open and write file "out36" for plotting
*     ----------------------------------------
      open(unit=36, file='out36', status='unknown', form= 'formatted')

*     -----------------------------------------------
*     Open file 46, "naoto" for synthetic diagnostic
*     -----------------------------------------------
      open(unit=46, file='naoto', status='unknown', form= 'formatted')



      open(unit=140,file='Efield2D_sum.vtk',status='unknown',
     &                                               form='formatted')
      open(unit=141,file='E_Sum_rect.vtk',status='unknown',
     &                                               form='formatted')
      open(unit=142,file='E_Sum_toroidal.vtk',status='unknown',
     &                                               form='formatted')
      open(unit=143,file='E_sum_3D.vtk', status='unknown',
     &                                               form='formatted')




!     ----------------------------------
!     take real part and mod of E_sum_3d
!     ----------------------------------
      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nphi3d

               ealpha_sum_real(i,j,k) = real(ealpha_sum(i, j, k))

               ealpha_sum_imag(i,j,k) = aimag(ealpha_sum(i, j, k))

               eb_sum_real(i,j,k) = real(eb_sum(i, j, k))

               eb_sum_imag(i,j,k) = aimag(eb_sum(i, j, k))

               eb_sum_mod(i,j,k) = sqrt(conjg(eb_sum(i, j, k))
     &                                          * eb_sum(i, j, k))

               ealpha_sum_mod(i,j,k) = sqrt(conjg(ealpha_sum(i, j, k))
     &                                          * ealpha_sum(i, j, k))

               ebeta_sum_mod(i,j,k) = sqrt(conjg(ebeta_sum(i, j, k))
     &                                         * ebeta_sum(i, j, k))

               mod_E_sum(i,j,k)   =
     &                 sqrt(conjg(ealpha_sum(i,j,k)) * ealpha_sum(i,j,k)
     &                    + conjg(ebeta_sum(i,j,k)) * ebeta_sum(i,j,k)
     &                    + conjg(eb_sum(i,j,k)) * eb_sum(i,j,k) )

               redotj(i,j,k)= 0.5 * real(
     &                     conjg(ealpha_sum(i,j,k)) * xjpx_sum(i,j,k)
     &                   + conjg(ebeta_sum(i,j,k))  * xjpy_sum(i,j,k)
     &                   + conjg(eb_sum(i,j,k))     * xjpz_sum(i,j,k))

            end do
         end do
      end do





      nt_max = nphi2 - nphi1 + 1

      write (36, 309) nphi1, nphi2, nt_max, nnoderho
      write (36, 310) prfin

      write (36, 310) (xnphi(nt),  nt = 1, nt_max)
      write (36, 310) (cnmod2(nt), nt = 1, nt_max)
      write (36, 310) (pabs(nt),   nt = 1, nt_max)
      write (36, 310) (jdriven(nt),   nt = 1, nt_max)
c      write (36, 310) (spa(nt),   nt = 1, nt_max)


      write (36, 310) (rhon(n),    n = 1, nnoderho)
      write (36, 310) (redotj2sum(n), n = 1, nnoderho)
      write (36, 310) (redotjesum(n), n = 1, nnoderho)
      write (36, 310) (redotj1sum(n), n = 1, nnoderho)
      write (36, 310) (redotj3sum(n), n = 1, nnoderho)
      write (36, 310) (redotj4sum(n), n = 1, nnoderho)
      write (36, 310) (redotj5sum(n), n = 1, nnoderho)

      write (36, 310) (redotji2_int(n), n = 1, nnoderho)
      write (36, 310) (redotje_int(n),  n = 1, nnoderho)
      write (36, 310) (redotji1_int(n), n = 1, nnoderho)
      write (36, 310) (redotji3_int(n), n = 1, nnoderho)
      write (36, 310) (redotji4_int(n), n = 1, nnoderho)
      write (36, 310) (redotji5_int(n), n = 1, nnoderho)

      write (36, 310) (xjprl_sum(n), n = 1, nnoderho)
      write (36, 310) (xjprl_int(n), n = 1, nnoderho)

      write (36, 310) (wdot2sum(n), n = 1, nnoderho)
      write (36, 310) (wdotesum(n),  n = 1, nnoderho)
      write (36, 310) (wdot1sum(n), n = 1, nnoderho)
      write (36, 310) (wdot3sum(n), n = 1, nnoderho)
      write (36, 310) (wdot4sum(n), n = 1, nnoderho)
      write (36, 310) (wdot5sum(n), n = 1, nnoderho)


      nxdim = nxmx
      nydim = nymx

      redotjisum = redotj1sum + redotj2sum + redotj3sum
     &           + redotj4sum + redotj5sum + redotj6sum

c      call run_rf2x(nmodesx, nmodesy, rwleft, rwright,
c     .   ytop, ybottom, myid, nxdim, nydim,
c     .   rt, b0, rho,
c     .   redotjesum, redotjisum,
c     .   redotj1sum, redotj2sum,
c     .   redotj3sum, redotj4sum,
c     .   redotj5sum, redotj6sum,
c     .   xjprl_sum, wdotesum,
c     .   wdot1sum,  wdot2sum,
c     .   wdot3sum,  wdot4sum,
c     .   wdot5sum,  wdot6sum)


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

         wdoti1_dvol(n+1) = wdot1sum(n) * dvol(n)
         wdoti2_dvol(n+1) = wdot2sum(n) * dvol(n)
         wdoti3_dvol(n+1) = wdot3sum(n) * dvol(n)
         wdoti4_dvol(n+1) = wdot4sum(n) * dvol(n)
         wdoti5_dvol(n+1) = wdot5sum(n) * dvol(n)
         wdoti6_dvol(n+1) = wdot6sum(n) * dvol(n)
         wdote_dvol (n+1) = wdotesum(n)  * dvol(n)

         redotj1_dvol(n+1) = redotj1sum(n) * dvol(n)
         redotj2_dvol(n+1) = redotj2sum(n) * dvol(n)
         redotj3_dvol(n+1) = redotj3sum(n) * dvol(n)
         redotj4_dvol(n+1) = redotj4sum(n) * dvol(n)
         redotj5_dvol(n+1) = redotj5sum(n) * dvol(n)
         redotj6_dvol(n+1) = redotj6sum(n) * dvol(n)
         redotje_dvol(n+1) = redotjesum(n) * dvol(n)
      end do

*     --------------------------------
*     Open and write SWIM output file:
*     --------------------------------

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





      write (36, 310) (wdot2sum_ql(n), n = 1, nnoderho)
      write (36, 310) (wdotesum_ql(n), n = 1, nnoderho)
      write (36, 310) (wdot1sum_ql(n), n = 1, nnoderho)
      write (36, 310) (wdot3sum_ql(n), n = 1, nnoderho)
      write (36, 310) (wdot4sum_ql(n), n = 1, nnoderho)
      write (36, 310) (wdot5sum_ql(n), n = 1, nnoderho)

      write (36, 310) (wdoti2_int(n), n = 1, nnoderho)
      write (36, 310) (wdote_int(n),  n = 1, nnoderho)
      write (36, 310) (wdoti1_int(n), n = 1, nnoderho)
      write (36, 310) (wdoti3_int(n), n = 1, nnoderho)
      write (36, 310) (wdoti4_int(n), n = 1, nnoderho)
      write (36, 310) (wdoti5_int(n), n = 1, nnoderho)

      write (36, 310) (dvol(n), n = 1, nnoderho)

      write (36, 309) nnodex, nnodey, nphi3d
      write (36, 310) (capr(i),  i = 1, nnodex)
      write (36, 310) (capz(j),  j = 1, nnodey)
      write (36, 310) (phi3d(k), k = 1, nphi3d)
      write (36, 310) (phi3d_shift(k), k = 1, nphi3d)
      write (36, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      write (36, 310) rhoplasm, rt

      write (36, 310) (((ealpha_sum(i, j, k), i = 1, nnodex),
     &                                       j = 1, nnodey),
     &                                       k = 1, nphi3d)

      write (36, 310) (((ealpha_sum_mod(i, j, k), i = 1, nnodex),
     &                                           j = 1, nnodey),
     &                                           k = 1, nphi3d)

      write (36, 310) (((ebeta_sum(i, j, k), i = 1, nnodex),
     &                                       j = 1, nnodey),
     &                                       k = 1, nphi3d)

      write (36, 310) (((ebeta_sum_mod(i, j, k), i = 1, nnodex),
     &                                           j = 1, nnodey),
     &                                           k = 1, nphi3d)

      write (36, 310) (((eb_sum(i, j, k), i = 1, nnodex),
     &                                       j = 1, nnodey),
     &                                       k = 1, nphi3d)


      write (36, 310) (((redotj(i, j, k), i = 1, nnodex),
     &                                    j = 1, nnodey),
     &                                    k = 1, nphi3d)






!     -------------------------------------------
!     Find maximum value of ebeta_sum_mod(i,j,1)
!     ------------------------------------------

      emax = 0.0

      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nphi3d
               if(ebeta_sum_mod(i,j,k) .gt. emax) then
                  emax = ebeta_sum_mod(i,j,k)
                  imax = i
                  jmax = j
                  kmax = k
               end if
            end do
         end do
      end do

      write(16, *) 'emax = ', emax
      write(16, *) 'imax = ', imax
      write(16, *) 'jmax = ', jmax
      write(16, *) 'kmax = ', kmax

      write(6, *) 'emax = ', emax
      write(6, *) 'imax = ', imax
      write(6, *) 'jmax = ', jmax
      write(6, *) 'kmax = ', kmax




!     ------------------------------------------
!     Write E_sum_3D.vtk file "rectilinear_grid"
!     ------------------------------------------
      number_points = nnodex * nnodey * nphi3d

      write(143, 2840)
 2840 format('# vtk DataFile Version 2.0')

      write(143, 2945)
 2945 format('3D electric field')

      write(143, 2846)
 2846 format('ASCII')

      write(143, 2847)
 2847 format('DATASET RECTILINEAR_GRID')

      write(143, 2841)nnodex, nnodey, nphi3d
 2841 format('DIMENSIONS', 3i8)

      write(143, 2942) nnodex
 2942 format('X_COORDINATES', i8, ' float')
      write (143, 3411) (capr(i), i = 1, nnodex)

      write(143, 2943) nnodey
 2943 format('Y_COORDINATES', i8, ' float')
      write (143, 3411) (capz(j), j = 1, nnodey)

      phi3d = signbz * phi3d

      write(143, 2951) nphi3d
 2951 format('Z_COORDINATES', i8, ' float')
      write (143, 3411) (phi3d(k), k = 1, nphi3d)

      write(143, 2844) number_points
 2844 format('POINT_DATA', i10)

      write(143, 2969)
 2969 format('SCALARS redotj float 1')
      write(143, 2849)
      write (143, 3411) (((redotj(i,j,k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)

      write(143, 2948)
 2948 format('SCALARS mod_E_sum float 1')
      write(143, 2849)
      write (143, 3411) (((mod_E_sum(i, j, k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)


      write(143, 2888)
 2888 format('SCALARS ealpha_sum_mod float 1')
      write(143, 2849)
      write (143, 3411) (((ealpha_sum_mod(i,j,k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)

      write(143, 2949)
 2949 format('SCALARS ealpha_sum_real float 1')
      write(143, 2849)
      write (143, 3411) (((ealpha_sum_real(i, j, k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)

      write(143, 2950)
 2950 format('SCALARS ealpha_sum_imag float 1')
      write(143, 2849)
      write (143, 3411) (((ealpha_sum_imag(i, j, k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)

      write(143, 2958)
 2958 format('SCALARS rho_sum float 1')
      write(143, 2849)
      write (143, 3411) (((rho_sum(i, j, k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)

      write(143, 2988)
 2988 format('SCALARS ntilda_real float 1')
      write(143, 2849)
      write (143, 3411) (((ntilda_real(i, j, k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)

      write(143, 2989)
 2989 format('SCALARS xjmod_sum float 1')
      write(143, 2849)
      write (143, 3411) (((xjmod_sum(i, j, k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)

      write(143, 2850)
 2850 format('VECTORS B_3d_vector float')

      write(143, 3412) (((bR_3d(i,j,k), bphi_3d(i,j,k), bZ_3d(i,j,k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)

 3412 format(3f17.4)

      write(143, 2851)
 2851 format('VECTORS Sp_3d_vector float')


      write(143, 3412) (((sR_3d(i,j,k), sphi_3d(i,j,k), sZ_3d(i,j,k),
     &   i = 1, nnodex), j = 1, nnodey), k = 1, nphi3d)

      close(143)




 3410 format(1p4e10.2)
 3411 format(6e16.4)


      if(ndisti1  .eq. 1) then

*        ----------------------------------------------------------------
*        Write quasilinear diffusion coefficients to file out_cql3d.coef1
*        ----------------------------------------------------------------
         open(unit=41, file='out_cql3d.coef1',
     &                               status='unknown', form='formatted')
         write (41, 309) nuper
         write (41, 309) nupar
         write (41, 309) nnoderho

         write (41, 3310) vc1_cgs
         write (41, 3310) UminPara, UmaxPara

         write (41, 3310) (rhon(n), n = 1, nnoderho)
         write (41, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
         write (41, 3310) (upara(i_upara), i_upara = 1, nupar)


         write (41, 3310) (((bqlsum_i1(i_uperp, i_upara, n),
     &      i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         write (41, 3310) (((cqlsum_i1(i_uperp, i_upara, n),
     &      i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         write (41, 3310) (((eqlsum_i1(i_uperp, i_upara, n),
     &      i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         write (41, 3310) (((fqlsum_i1(i_uperp, i_upara, n),
     &      i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         write (41, 3310) xmi1

         close (41)
         call dql_write_nc('out_cql3d.coef1','out_cql3d.coef1.nc')

      end if

      if(ndisti2  .eq. 1) then
*        ----------------------------------------------------------------
*        Write quasilinear diffusion coefficients to file out_cql3d.coef2
*        ----------------------------------------------------------------
         open(unit=42, file='out_cql3d.coef2',
     &                               status='unknown', form='formatted')
         write (42, 309) nuper
         write (42, 309) nupar
         write (42, 309) nnoderho

         write (42, 3310) vc2_cgs
         write (42, 3310) UminPara, UmaxPara

         write (42, 3310) (rhon(n), n = 1, nnoderho)
         write (42, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
         write (42, 3310) (upara(i_upara), i_upara = 1, nupar)

         write (42, 3310) (((bqlsum_i2(i_uperp, i_upara, n),
     &      i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         write (42, 3310) (((cqlsum_i2(i_uperp, i_upara, n),
     &      i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         write (42, 3310) (((eqlsum_i2(i_uperp, i_upara, n),
     &      i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         write (42, 3310) (((fqlsum_i2(i_uperp, i_upara, n),
     &      i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         write (42, 3310) xmi2

         close (42)
         call dql_write_nc('out_cql3d.coef2','out_cql3d.coef2.nc')

      end if



      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ealpha_sum(i, j, 1))
            fimag(i,j) = aimag(ealpha_sum(i, j, 1))
            power(i,j) = real(redotj(i, j, 1))
            fluct(i,j) = real(ntilda(i, j, 1)) / 1.0e+20
         end do
      end do

!     -----------------------------------------------
!     Write Efield2D_sum.vtk file "structured points"
!     -----------------------------------------------
      number_points = nnodex * nnodey

      write(140, 2840)

      write(140, 2845)
 2845 format('Real E_alpha')

      write(140, 2846)

c      write(140, 2847)
       write(140, 3847)
 3847 format('DATASET STRUCTURED_POINTS')

      write(140, 2841) nnodex,  nnodey, 1


      write(140, 2842) capr(1), capz(1), 0
 2842 format('ORIGIN', 2f9.3, 1i8)

      write(140, 2843) dx, dy, 1
 2843 format('SPACING', 2f9.5, 1i8)

      write(140, 2844) number_points

      write(140, 2848)
 2848 format('SCALARS Real_E_alpha float 1')

      write(140, 2849)
 2849 format('LOOKUP_TABLE default')

      write (140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 8854)
 8854 format('SCALARS Imag_E_alpha float 1')
      write(140, 2849)
      write (140, 3411) ((fimag(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3851)
 3851 format('SCALARS rho float 1')
      write(140, 2849)
      write (140, 3411) ((rho(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3852)
 3852 format('SCALARS redotj float 1')
      write(140, 2849)
      write (140, 3411) ((power(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3882)
 3882 format('SCALARS fluctuation float 1')
      write(140, 2849)
      write (140, 3411) ((fluct(i,j), i = 1, nnodex),  j = 1, nnodey)

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(eb_sum(i, j, 1))
            fimag(i,j) = aimag(eb_sum(i, j, 1))
         end do
      end do

      write(140, 88854)
88854 format('SCALARS Real_Eb float 1')
      write(140, 2849)
      write (140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 98854)
98854 format('SCALARS Imag_Eb float 1')
      write(140, 2849)
      write (140, 3411) ((fimag(i,j), i = 1, nnodex),  j = 1, nnodey)

*     -----------------------------------------------------------
*     plot electric field in poloidal cross-section (one quarter)
*     -----------------------------------------------------------

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ealpha_sum(i, j, nphi3d_1quarter))
            fimag(i,j) = aimag(ealpha_sum(i, j, nphi3d_1quarter))
            power(i,j) = real(redotj(i, j, nphi3d_1quarter))
         end do
      end do


      write(140, 3861)
 3861 format('SCALARS ealpha_re_1_4th float 1')
      write(140, 2849)
      write (140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)


*     --------------------------------------------------------
*     plot electric field in poloidal cross-section (one half)
*     --------------------------------------------------------

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ealpha_sum(i, j, nphi3d_half))
            fimag(i,j) = aimag(ealpha_sum(i, j, nphi3d_half))
            power(i,j) = real(redotj(i, j, nphi3d_half))
         end do
      end do



      write(140, 3862)
 3862 format('SCALARS ealpha_re_2_4th float 1')
      write(140, 2849)
      write (140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)



*     -----------------------------------------------------------
*     plot electric field in poloidal cross-section (three quarter)
*     -----------------------------------------------------------

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ealpha_sum(i, j, nphi3d_3quarter))
            fimag(i,j) = aimag(ealpha_sum(i, j, nphi3d_3quarter))
            power(i,j) = real(redotj(i, j, nphi3d_3quarter))
         end do
      end do



      write(140, 3863)
 3863 format('SCALARS ealpha_re_3_4th float 1')
      write(140, 2849)
      write (140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)


      close (140)




      jhalf = nnodey / 2
      khalf = nphi3d / 2

      write(6, *) "jhalf = ", jhalf
      write(6, *) "nnodex = ", nnodex
      write(6, *) "nphi3d = ", nphi3d
      write(6, *) "signbz = ", signbz

      write(16, *) "jhalf = ", jhalf
      write(16, *) "nnodex = ", nnodex
      write(16, *) "nphi3d = ", nphi3d
      write(16, *) "signbz = ", signbz


      do i = 1, nnodex
         do k = 1, nphi3d

            if (k .lt. khalf) kshift = k + khalf
            if (k .ge. khalf) kshift = k - khalf + 1

            freal_phi(i, k) = real(ealpha_sum(i, jhalf, kshift))
            fimag_phi(i, k) = aimag(ealpha_sum(i, jhalf, kshift))
            power_phi(i, k) = real(redotj(i, jhalf, kshift))
            capr_plot_phi(i, k) = capr(i)

         end do
      end do



!     ---------------------------------------------
!     Write E_rect_2D.vtk file "structured points"
!     --------------------------------------------
      number_points = nnodex * nphi3d
      dx = capr(2) - capr(1)
      dy = phi3d_shift(2) - phi3d_shift(1)

      write(141, 2840)
      write(141, 2852)
 2852 format('Real E alpha')

      write(141, 2846)
c      write(141, 2847)

      write(141, 6848)
 6848 format('DATASET STRUCTURED_POINTS')

      write(141, 2841) nnodex, nphi3d, 1
      write(141, 2842) capr(1), phi3d_shift(1), 0
      write(141, 2843) dx, dy, 1
      write(141, 2844) number_points

      write(141, 3848)
 3848 format('SCALARS Re_E_alpha float 1')
      write(141, 2849)
      write (141, 3411) ((freal_phi(i, k), i = 1, nnodex), k = 1,nphi3d)

      write(141, 3849)
 3849 format('SCALARS Im_E_alphafloat 1')
      write(141, 2849)
      write (141, 3411) ((fimag_phi(i, k), i = 1, nnodex), k = 1,nphi3d)

      write(141, 3859)
 3859 format('SCALARS Re_EdotJ float 1')
      write(141, 2849)
      write (141, 3411) ((power_phi(i, k), i = 1, nnodex), k = 1,nphi3d)


      close (141)



*     --------------------------------------
*     plot in toroidal cross-section (polar)
*     --------------------------------------

*     ---------------------------------------
*     Interpolate E_alpha onto Cartesian grid
*     ---------------------------------------
      do i = 1, nxplot
         do j = 1, nyplot

            capr_giv = sqrt(xplot(i)**2 + yplota(j)**2)

            if(capr_giv .ge. capr_min .and. capr_giv .le. capr_max) then

               x_giv = capr_giv - capr_min
               phi3d_giv = atan2(yplota(j), xplot(i))
               if(phi3d_giv .lt. 0.0) phi3d_giv = phi3d_giv + 2.0 * pi

               call intplt_phi(x_giv, phi3d_giv, frealp(i, j), nnodex,
     &                 nphi3d, freal_phi, nxmx, nphimx, dcapr, dphi3d)

               call intplt_phi(x_giv, phi3d_giv, fimagp(i, j), nnodex,
     &                 nphi3d, fimag_phi, nxmx, nphimx, dcapr, dphi3d)

               call intplt_phi(x_giv, phi3d_giv, capr_plotp(i,j),nnodex,
     &               nphi3d, capr_plot_phi, nxmx, nphimx, dcapr, dphi3d)

               call intplt_phi(x_giv, phi3d_giv, powerp(i, j), nnodex,
     &                  nphi3d, power_phi, nxmx, nphimx, dcapr, dphi3d)

            end if

         end do
      end do

      write(6, *) "capr_max = ", capr_max
      write(6, *) "capr_min = ", capr_min

      write(16, *) "capr_max = ", capr_max
      write(16, *) "capr_min = ", capr_min


      write (36, 309) nxplot, nyplot
      write (36, 310) capr_max, capr_min
      write (36, 310) (xplotm(i),  i = 1, nxplot)
      write (36, 310) (yplotm(j),  j = 1, nyplot)
      write (36, 310) ((frealp(i, j), i = 1, nxplot), j = 1, nyplot)
      write (36, 310) ((fimagp(i, j), i = 1, nxplot), j = 1, nyplot)
      write (36, 310) ((capr_plotp(i, j), i = 1, nxplot), j = 1, nyplot)
      write (36, 310) ((powerp(i, j), i = 1, nxplot), j = 1, nyplot)

      write (36, 310) (((ntilda(i, j, k), i = 1, nnodex),
     &                                    j = 1, nnodey),
     &                                    k = 1, nphi3d)

      close (36)


*     -----------------------------------------------
*     Write file 46, "naoto" for synthetic diagnostic
*     -----------------------------------------------

      write (46, 309) nnodex, nnodey, nphi3d
      write (46, 310) (capr(i),  i = 1, nnodex)
      write (46, 310) (capz(j),  j = 1, nnodey)
      write (46, 310) (phi3d(k), k = 1, nphi3d)

      write (46, 310) (((ealpha_sum(i, j, k), i = 1, nnodex),
     &                                       j = 1, nnodey),
     &                                       k = 1, nphi3d)

      write (46, 310) (((ebeta_sum(i, j, k), i = 1, nnodex),
     &                                       j = 1, nnodey),
     &                                       k = 1, nphi3d)

      write (46, 310) (((eb_sum(i, j, k), i = 1, nnodex),
     &                                    j = 1, nnodey),
     &                                    k = 1, nphi3d)

      write (46, 310) (((ntilda(i, j, k), i = 1, nnodex),
     &                                    j = 1, nnodey),
     &                                    k = 1, nphi3d)
      close (46)




!     ------------------------------------------------
!     Write E_toroidal_2D.vtk file "structured points"
!     ------------------------------------------------
      number_points = nxplot * nyplot
      dx = xplotm(2) - xplotm(1)
      dy = yplotm(2) - yplotm(1)



      write(142, 2840)

      write(142, 2853)
 2853 format('Real E_toroidal_2D')

      write(142, 2846)
c      write(142, 2847)

      write(142, 6849)
 6849 format('DATASET STRUCTURED_POINTS')

      write(142, 2841) nxplot, nyplot, 1
      write(142, 2842) xplotm(1), yplotm(1), 0
      write(142, 2843) dx, dy, 1
      write(142, 2844) number_points

      write(142, 4848)
 4848 format('SCALARS E_alpha_re float 1')
      write(142, 2849)
      write (142, 3411) ((frealp(i, j), i = 1, nxplot), j = 1, nyplot)

      write(142, 3854)
 3854 format('SCALARS E_alpha_im float 1')
      write(142, 2849)
      write (142, 3411) ((fimagp(i, j), i = 1, nxplot), j = 1, nyplot)

      write(142, 8858)
 8858 format('SCALARS Real_EdotJ float 1')
      write(142, 2849)
      write(142, 3411) ((powerp(i, j), i = 1, nxplot), j = 1, nyplot)

      write(142, 3853)
 3853 format('SCALARS capr_plot float 1')
      write(142, 2849)
      write(142, 3411)((capr_plotp(i,j), i = 1, nxplot), j = 1, nyplot)

      close (142)






      write(6, 163)
      write(6, 1009)ptotal

      write(16, 163)
      write(16, 1009)ptotal

      write(6, *) "fraction absorbed in edge = ", fedge, " %"
      write(6, *) "fraction absorbed in core = ", fcore, " %"
      write(16, *) "fraction absorbed in edge = ", fedge, " %"
      write(16, *) "fraction absorbed in core = ", fcore, " %"


      write(6, 163)
      write(6, 2011)jdriven_sum

      write(16, 163)
      write(16, 2011)jdriven_sum




      deallocate( xplot )
      deallocate( yplota)
      deallocate( xplotm )
      deallocate( yplotm )

      deallocate( ealpha_sum_real )
      deallocate( ealpha_sum_imag )

      deallocate( eb_sum_real )
      deallocate( eb_sum_imag )

      deallocate( ealpha_sum_mod )

      deallocate( ebeta_sum_mod )
      deallocate( mod_E_sum)

      deallocate( bR_3d )
      deallocate( bphi_3d )
      deallocate( bZ_3d )

      deallocate (ealpha_sum)
      deallocate( ebeta_sum )
      deallocate( eb_sum )

      deallocate( ntilda )
      deallocate( ntilda_real )

      deallocate( xjpx_sum )
      deallocate( xjpy_sum )
      deallocate( xjpz_sum )

      deallocate( xjx_sum )
      deallocate( xjy_sum )
      deallocate( xjz_sum )
      deallocate( xjmod_sum )

      deallocate( bxwave_sum )
      deallocate( bywave_sum )
      deallocate( bzwave_sum )

      deallocate( ex_sum )
      deallocate( ey_sum )
      deallocate( ez_sum )

      deallocate( xjpxe_lab_sum )
      deallocate( xjpye_lab_sum )
      deallocate( xjpze_lab_sum )

      deallocate( ntilda_e_sum )

      deallocate( redotj )
      deallocate( rho_sum )

      deallocate( frealp )
      deallocate( fimagp )
      deallocate( powerp )
      deallocate( capr_plotp )


      deallocate( freal )
      deallocate( fimag )
      deallocate( capr_plot )
      deallocate( power )
      deallocate( fluct )

      deallocate( freal_phi )
      deallocate( fimag_phi )
      deallocate( capr_plot_phi )
      deallocate( power_phi )

      deallocate( ealpha )
      deallocate( ebeta )
      deallocate( eb )

      deallocate( xjpx )
      deallocate( xjpy )
      deallocate( xjpz )

      deallocate( xjx )
      deallocate( xjy )
      deallocate( xjz )

      deallocate( bx )
      deallocate( by )
      deallocate( bz )

      deallocate( xjpxe_lab )
      deallocate( xjpye_lab )
      deallocate( xjpze_lab )

      deallocate( sR_3d )
      deallocate( sZ_3d )
      deallocate( sphi_3d )

      deallocate( ntilda_e )

      deallocate( eplus )
      deallocate( eminus )

      deallocate (ex)
      deallocate (ey)
      deallocate (ez)
      deallocate (bxwave)
      deallocate (bywave)
      deallocate (bzwave)

      tmin = (second1(dummy) - t1) / 60.
      write(6 , 8847) tmin
      write(16, 8847) tmin
 8847 format('time to sum modes =', f9.3, ' min')

      call plot(ndisti2)

      close (16)

  309 format(10i10)
 3310 format(1p6e18.10)
11001 format(1e16.8,2e15.8)
99998 format(" ",f12.4, 1p9e12.4)
99997 format(i10)
88887 format(2i10)
  899 format(" ",5x,"      cpu = ",1p,e12.4," min"/
     1       " ",5x,"       io = ",1p,e12.4," min"/
     1       " ",5x,"      sys = ",1p,e12.4," min"/
     1       " ",5x,"      mem = ",1p,e12.4," min"/
     1       " ",5x,"    total = ",1p,e12.4," min")
 2000 format(1i10,1e10.3,1i10,1e10.3)
 2001 format(1i10,1p9e12.4)
 2517 format(1i10,1e10.3,1i10,1e10.3)
 2211 format(8i10)
 2220 format(3i10,3e10.3,2i10)
 2230 format(2i10,7e10.3)
 2240 format(1i10,7e10.3)
 2260 format(1i10,1e10.3,1i10,1e10.3)
 2280 format(2i10,5e10.3,1i10)
 1009 format(3x,"   power absorbed  = ",1p,e12.4," Watts   ")
  300 format(i10/i10,4e10.3,i10)
  310 format(1p6e12.4)
  311 format(1p10e12.4)
  122 format(1i10,1p9e12.4)
  123 format(2i10,1p9e12.4)
 1051 format(6i10)
 1052 format(1i10,6e12.4)
 1053 format(6e12.4)
 1013 format(3x,"         lt  = ",1p,e12.4," m       ")
 1014 format(3x,"         wd  = ",1p,e12.4," m       ")
 1015 format(3x,"       dphi  = ",1p,e12.4," radians ")
 3015 format(3x,"       phi0  = ",1p,e12.4," radians ")
 3016 format(3x,"        xi0  = ",1p,e12.4," Amps    ")
c     skip a line
  163 format("0")
 1163 format(" ",5x,"  k   ",2x," phi  ",6x,"re cur",4x," im cur ",
     1       4x," re jsum",4x," im jsum",
     1       4x," re csum",4x," im csum",
     1       4x,"        ",4x,"        ")
 1164 format(" ",5x,"nphi  ",2x,"  k"// ,6x," pabs ",5x," weight ",
     1       5x," re pc  ",4x," im pc  ",
     1       4x," pabs wt",4x," gammacd",
     1       4x,"  damp1 ",4x,"        ")
 1165 format(" ",5x,"  k   ",2x," phi  ",6x,"re cur",4x," im cur ",
     1       4x,"re jthsm",4x,"im jthsm",
     1       4x," re epsi",4x," im epsi",
     1       4x,"        ",4x,"        ")
  162 format("1")
 2011 format(3x,"total fast wave current driven = ",1p,e12.4," Amps   ")
 2016 format(3x,"total ohmic current = ",1p,e12.4," Amps    ")
 2017 format(3x,"total boot strap current = ",1p,e12.4," Amps    ")
 7014 format(3x,"current driven per watt = ",1p,e12.4," Amps/watt")
 7015 format(3x,"current drive efficiency = ",1p,e12.4,
     1   " Amps/watt/m**2")
 7016 format(3x," average single pass damping = ",1p,e12.4,
     1   "               ")
 2009 format(3x,"power absorbed by electrons = ",1p,e12.4," watts   ")
 2012 format(3x,"power absorbed by alphas = ",1p,e12.4," watts   ")
 2013 format(3x,"power absorbed by majority = ",1p,e12.4," watts   ")
 2014 format(3x,"power absorbed by minority = ",1p,e12.4," watts   ")
 2015 format(3x,"power absorbed by impurity = ",1p,e12.4," watts   ")
 1008 format(3x,"      phi antenna  = ",1p,e12.4," pi      ")
 1010 format(3x,"   real power transferred from antenna  = ",
     1   1p,e12.4," watts   "/
     1      ,3x,"   imaginary power transferred from antenna  = ",
     1   1p,e12.4," watts   ")
 1011 format(3x,"real part of antenna impedance (resistance) = ",
     1   1p,e12.4," ohms    "/
     1      ,3x,"imaginary part of antenna impedance (reactance) = ",
     1   1p,e12.4," ohms    ")
 1016 format(3x,"total power coupled at antenna = ",
     1   1p,e12.4," watts   ")
 1002 format(2i10,7e10.3)
 2165 format(3i10,5e10.3)
 1101 format(2i10,5e10.3,1i10)
 1001 format(8e10.3)
 1000 format(1i10,7e10.3)
11009 format(3x,"        frequency  = ",1p,e12.4," hertz ")
81822 format(3x,"            xwall  = ",1p,e12.4," m     ")
81823 format(3x,"           xnwall  = ",1p,e12.4," m-3   ")
81824 format(3x,"             amu1  = ",1p,e12.4," amu   ")
81825 format(3x,"             amu2  = ",1p,e12.4," amu   ")
81826 format(3x,"             amu3  = ",1p,e12.4," amu   ")
98824 format(3x,"             amu4  = ",1p,e12.4," amu   ")
98825 format(3x,"             amu5  = ",1p,e12.4," amu   ")
98826 format(3x,"             amu6  = ",1p,e12.4," amu   ")
98827 format(3x,"             amu7  = ",1p,e12.4," amu   ")
81827 format(3x,"               z1  = ",1p,e12.4," q proton")
81828 format(3x,"               z2  = ",1p,e12.4," q proton")
81829 format(3x,"               z3  = ",1p,e12.4," q proton")
98834 format(3x,"               z4  = ",1p,e12.4," q proton")
98835 format(3x,"               z5  = ",1p,e12.4," q proton")
98836 format(3x,"               z6  = ",1p,e12.4," q proton")
98837 format(3x,"               z7  = ",1p,e12.4," q proton")
81830 format(3x,"              eta  = ",1p,e12.4,"       ")
81831 format(3x,"             eta3  = ",1p,e12.4,"       ")
91834 format(3x,"             eta4  = ",1p,e12.4,"       ")
91835 format(3x,"             eta5  = ",1p,e12.4,"       ")
91836 format(3x,"             eta6  = ",1p,e12.4,"       ")
91837 format(3x,"             eta7  = ",1p,e12.4,"       ")
81832 format(3x,"            xnlim  = ",1p,e12.4," m-3   ")
81833 format(3x,"            qavg0  = ",1p,e12.4,"       ")
81834 format(3x,"           scrape  = ",1p,e12.4," m     ")
81835 format(3x,"           alphac  = ",1p,e12.4,"       ")
81836 format(3x,"           alphan  = ",1p,e12.4,"       ")
81837 format(3x,"           alphat  = ",1p,e12.4,"       ")
81838 format(3x,"           zeffcd  = ",1p,e12.4," q proton")
81839 format(3x,"            dpsi0  = ",1p,e12.4,"       ")
81840 format(3x,"            telim  = ",1p,e12.4," eV    ")
81841 format(3x,"            tilim  = ",1p,e12.4," eV    ")
91841 format(3x,"           ti2lim  = ",1p,e12.4," eV    ")
91842 format(3x,"           ti3lim  = ",1p,e12.4," eV    ")
91844 format(3x,"           ti4lim  = ",1p,e12.4," eV    ")
91845 format(3x,"           ti5lim  = ",1p,e12.4," eV    ")
91846 format(3x,"           ti6lim  = ",1p,e12.4," eV    ")
91847 format(3x,"           ti7lim  = ",1p,e12.4," eV    ")
81842 format(3x,"           deltap  = ",1p,e12.4,"       ")
81843 format(3x,"                p  = ",1p,e12.4,"       ")
81845 format(3x,"            q07qa  = ",1p,e12.4,"       ")
 1012 format(3x,"             omgrf = ",1p,e12.4," hertz ")
 1714 format(3x,"               xk0 = ",1p,e12.4," m-1   ")
 1021 format(3x,"                nz = ",1p,e12.4,"       ")
 1921 format(3x,"           w phase = ",1p,e12.4,"       ")
 1812 format(3x,"                rt = ",1p,e12.4," m     ")
 1835 format(3x,"            ekappa = ",1p,e12.4,"       ")
 1822 format(3x,"            aplasm = ",1p,e12.4," m     ")
 1823 format(3x,"              xant = ",1p,e12.4," m     ")
 1809 format(3x,"                b0 = ",1p,e12.4," T     ")
81844 format(3x,"              xne0 = ",1p,e12.4," m-3   ")
 1813 format(3x,"              xn10 = ",1p,e12.4," m-3   ")
 1814 format(3x,"              xn20 = ",1p,e12.4," m-3   ")
 1834 format(3x,"              xn30 = ",1p,e12.4," m-3   ")
 9834 format(3x,"              xn40 = ",1p,e12.4," m-3   ")
 9835 format(3x,"              xn50 = ",1p,e12.4," m-3   ")
 9836 format(3x,"              xn60 = ",1p,e12.4," m-3   ")
 9837 format(3x,"              xn70 = ",1p,e12.4," m-3   ")
 8834 format(3x,"              xn40 = ",1p,e12.4," m-3   ")
 8835 format(3x,"              xn50 = ",1p,e12.5," m-3   ")
 1815 format(3x,"               te0 = ",1p,e12.4," eV    ")
 1821 format(3x,"               ti0 = ",1p,e12.4," eV    ")
 1922 format(3x,"              ti02 = ",1p,e12.4," eV    ")
 1923 format(3x,"              ti03 = ",1p,e12.4," eV    ")
 1924 format(3x,"              ti04 = ",1p,e12.4," eV    ")
 1925 format(3x,"              ti05 = ",1p,e12.4," eV    ")
 1926 format(3x,"              ti06 = ",1p,e12.4," eV    ")
 1927 format(3x,"              ti07 = ",1p,e12.4," eV    ")
 1017 format(3x,"        xnu1/omgrf = ",1p,e12.4,"       ")
 1018 format(3x,"        xnu2/omgrf = ",1p,e12.4,"       ")
 1019 format(3x,"        xnu3/omgrf = ",1p,e12.4,"       ")
 1054 format(3x,"           flandau = ",1p,e12.4,"       ")
 1022 format(3x,"                nl = ",i12,"       ")
71022 format(3x,"                nr = ",i12,"       ")
71023 format(3x,"            ntheta = ",i12,"       ")
71024 format(3x,"             mdiag = ",i12,"       ")
71034 format(3x,"            mpdiag = ",i12,"       ")
73031 format(3x,"             jharm = ",i12,"       ")
73032 format(3x,"              lmax = ",i12,"       ")
73033 format(3x,"               iez = ",i12,"       ")
73034 format(3x,"             isort = ",i12,"       ")
71025 format(3x,"             ndiag = ",i12,"       ")
71026 format(3x,"             ikprl = ",i12,"       ")
71027 format(3x,"             igeom = ",i12,"       ")
71028 format(3x,"             nphi1 = ",i12,"       ")
71029 format(3x,"             nphi2 = ",i12,"       ")
71030 format(3x,"             ncoll = ",i12,"       ")
91031 format(3x,"            ncoll7 = ",i12,"       ")
71031 format(3x,"             imode = ",i12,"       ")
31012 format(3x,"    power absorbed by electrons = ",
     1   1p,e12.4,"  %      ")
31013 format(3x,"    power absorbed by majority  = ",
     1   1p,e12.4,"  %      ")
31014 format(3x,"    power absorbed by minority  = ",
     1   1p,e12.4,"  %      ")
31015 format(3x,"    power absorbed by species 3 = ",
     1   1p,e12.4,"  %      ")
31016 format(3x,"    power absorbed by species 4 = ",
     1   1p,e12.4,"  %      ")
31017 format(3x,"    power absorbed by species 5 = ",
     1   1p,e12.4,"  %      ")
31018 format(3x,"    power absorbed by species 6 = ",
     1   1p,e12.4,"  %      ")
31019 format(3x,"    power absorbed by species 7 = ",
     1   1p,e12.4,"  %      ")

 5000 continue




      return
      end

c
c********************************************************************
c



      complex function jfun(z)

      use aorsa2din_mod

      implicit none

      integer:: n

      real:: zn, pi, z0, zdiff2, phin, zdiff, zdiff1,
     &   z, dphi

      complex:: zi, fun1

      pi = 3.141592654
      zi = cmplx(0., 1.)

      z0 = -wd * (real(nstrap) - 1.) / 2.


      jfun = cmplx(0., 0.)

      do n = 1, nstrap
         zn = z0 + (n - 1) * wd

         if(zn .lt. zmin)zn = zn + zmax
         if(zn .gt. zmax)zn = zn - zmax

         zdiff = z - zn

         if(zn .eq. 0) then
            zdiff1 = z - zn
            zdiff2 = zmax - z
            zdiff  = amin1(zdiff1, zdiff2)
         endif

         phin = phase_array(n)
c        write (6, *) "phin = ", phin

         jfun = jfun + cexp(zi * phin) * fun1(zdiff) * amplt(n)

      end do

c      call exit

      return
      end
c
c***************************************************************************
c
      complex function fun1(z)

      use aorsa2din_mod

      implicit none

      real:: zabs, z, dphi

      zabs=abs(z)

      fun1=cmplx(1.,0.)

      if(zabs .gt. xlt/2.)fun1=cmplx(0.0, 0.0)

      return
      end


c
c***************************************************************************
c
      subroutine fftn2(ftheta, work, ntheta, nthmax, mpmin, mpmax, cn)

      implicit none

      integer:: nhalf, mpol, m, nthmax, ntheta, mpmin, mpmax

      complex:: ftheta(nthmax), work(nthmax),  cn(mpmin : mpmax)

      nhalf = ntheta / 2

      do m = 1, ntheta
         work(m) = ftheta(m)
      end do

      call four1(work, ntheta, -1)

      do m = 1, ntheta
         mpol = m - 1
         if(mpol .gt. nhalf)mpol = mpol - ntheta
         cn(mpol) = work(m) / ntheta
      end do

      return
      end

c
c***************************************************************************
c

      subroutine dgrate(x, f, nx1, nx2, ans, nxmax)

      implicit none

      integer:: n, nx1, nx2, nxmax
      real:: f(nxmax), x(nxmax), ans

      ans=0.0

      do n=nx1,nx2-1
         ans=ans+(f(n)+f(n+1))/2.*(x(n+1)-x(n))
      end do

      return
      end


c
c***************************************************************************
c
      subroutine fftin2(ftheta, work, ntheta, nthmax, mpmin, mpmax, cn)

      implicit none

      integer:: nhalf, mpol, m, nthmax, ntheta, mpmin, mpmax

      complex:: ftheta(nthmax), work(nthmax), cn(mpmin:mpmax)

      nhalf = ntheta / 2

      do m = 1, ntheta
         mpol = m - 1
         if(mpol .gt. nhalf) mpol = mpol - ntheta
         work(m) = cn(mpol)
      end do

      call four1(work, ntheta, +1)

      do m = 1, ntheta
         ftheta(m) = work(m)
      end do

      return
      end

c
c***************************************************************************
c

      subroutine sgrate_3d(x, y, phi, f, nx1, nx2, ny1, ny2,
     &   nphi1, nphi2, ans, capr, nxmax, nymax, nphimx)

      implicit none

      integer:: i, j, k, nx1, nx2, ny1, ny2, nphi1, nphi2
      integer:: nxmax, nymax, nphimx

      real:: f(nxmax, nxmax, nphimx), capr(nxmax),
     &      x(nxmax), y(nymax), phi(nphimx)

      real:: ans, dx, dy, dphi

      ans = 0.0

      do k = nphi1, nphi2 - 1
         dphi = phi(k+1) - phi(k)

         do i = nx1, nx2 - 1
            dx = x(i+1) - x(i)

            do j = ny1, ny2 - 1
               dy = y(j+1) - y(j)

               ans = ans + dx * dy * dphi *
     &                           (capr(i)   * f(i, j, k)
     &                          + capr(i+1) * f(i+1, j, k)
     &                          + capr(i)   * f(i, j+1, k)
     &                          + capr(i+1) * f(i+1, j+1, k)

     &                          + capr(i)   * f(i, j, k+1)
     &                          + capr(i+1) * f(i+1, j, k+1)
     &                          + capr(i)   * f(i, j+1, k+1)
     &                          + capr(i+1) * f(i+1, j+1, k+1)  )/ 8.0
            end do
         end do
      end do


      return
      end

c
c***************************************************************************
c

      subroutine sgrate_3d_edge(x, y, phi, f, nx1, nx2, ny1, ny2,
     &   nphi1, nphi2, ans, capr, nxmax, nymax, nphimx, rho)

      implicit none

      integer:: i, j, k, nx1, nx2, ny1, ny2, nphi1, nphi2
      integer:: nxmax, nymax, nphimx

      real:: f(nxmax, nxmax, nphimx), capr(nxmax),
     &      x(nxmax), y(nymax), phi(nphimx), rho(nxmax, nymax)

      real:: ans, dx, dy, dphi

      ans = 0.0

      do k = nphi1, nphi2 - 1
         dphi = phi(k+1) - phi(k)

         do i = nx1, nx2 - 1
            dx = x(i+1) - x(i)

            do j = ny1, ny2 - 1
               dy = y(j+1) - y(j)

               if (rho(i,j) .gt. 1.0) then

                  ans = ans + dx * dy * dphi *
     &                           (capr(i)   * f(i, j, k)
     &                          + capr(i+1) * f(i+1, j, k)
     &                          + capr(i)   * f(i, j+1, k)
     &                          + capr(i+1) * f(i+1, j+1, k)

     &                          + capr(i)   * f(i, j, k+1)
     &                          + capr(i+1) * f(i+1, j, k+1)
     &                          + capr(i)   * f(i, j+1, k+1)
     &                          + capr(i+1) * f(i+1, j+1, k+1)  )/ 8.0
                end if

            end do
         end do
      end do


      return
      end

c
c***************************************************************************
c

      subroutine sgrate_3d_core(x, y, phi, f, nx1, nx2, ny1, ny2,
     &   nphi1, nphi2, ans, capr, nxmax, nymax, nphimx, rho)

      implicit none

      integer:: i, j, k, nx1, nx2, ny1, ny2, nphi1, nphi2
      integer:: nxmax, nymax, nphimx

      real:: f(nxmax, nxmax, nphimx), capr(nxmax),
     &      x(nxmax), y(nymax), phi(nphimx), rho(nxmax, nymax)

      real:: ans, dx, dy, dphi

      ans = 0.0

      do k = nphi1, nphi2 - 1
         dphi = phi(k+1) - phi(k)

         do i = nx1, nx2 - 1
            dx = x(i+1) - x(i)

            do j = ny1, ny2 - 1
               dy = y(j+1) - y(j)

               if (rho(i,j) .le. 1.0) then

                  ans = ans + dx * dy * dphi *
     &                           (capr(i)   * f(i, j, k)
     &                          + capr(i+1) * f(i+1, j, k)
     &                          + capr(i)   * f(i, j+1, k)
     &                          + capr(i+1) * f(i+1, j+1, k)

     &                          + capr(i)   * f(i, j, k+1)
     &                          + capr(i+1) * f(i+1, j, k+1)
     &                          + capr(i)   * f(i, j+1, k+1)
     &                          + capr(i+1) * f(i+1, j+1, k+1)  )/ 8.0
                end if

            end do
         end do
      end do


      return
      end


c
c***************************************************************************
c
      subroutine psigr3(f, npsi, ntheta, nphi, ans, psi, theta, phi,
     1   rootg, r0, nrmax, nthmax, nphimx)

      data pi/3.141592654/

      dimension f(nrmax,nthmax,nphimx), psi(nrmax),
     1   theta(nthmax), phi(nphimx),
     1   rootg(nrmax,nthmax)

      ans=0.0

      do k=1,nphi

         if(k.eq.2)ansz1=ansz
         if(k.gt.1)anszm1=ansz
         if(k.gt.1)dphi=phi(k)-phi(k-1)

         ansz=0.0

         do n=1,npsi-1
            dpsi=psi(n+1)-psi(n)

            do m=1,ntheta-1
               dtheta=theta(m+1)-theta(m)
               ansz=ansz+dpsi*dtheta
     &            *(rootg(n,m)*f(n,m,k)
     &            +rootg(n+1,m)*f(n+1,m,k)
     &            +rootg(n,m+1)*f(n,m+1,k)
     &            +rootg(n+1,m+1)*f(n+1,m+1,k))/4.
            end do
         end do

         m=ntheta
         dtheta=2.0*pi-theta(m)

         do n=1,npsi-1
            dpsi=psi(n+1)-psi(n)
            ansz = ansz + dpsi * dtheta
     &        *(rootg(n,m) * f(n,m,k)
     &         +rootg(n+1,m) * f(n+1,m,k)
     &         +rootg(n,1) * f(n,1,k)
     &         +rootg(n+1,1) * f(n+1,1,k))/4.0
         end do

         if(k.gt.1)ans=ans+(anszm1+ansz)/2.0*dphi

      end do

      k=nphi+1
      dphi=2.*pi-phi(k-1)
      ans=ans+(ansz+ansz1)/2.0*dphi

      return
      end
c
c***************************************************************************
c
      subroutine psig3c(f,npsi,ntheta,nphi,ans,psi,theta,phi,
     &   rootg,r0,nrmax,nthmax,nphimx)
      real, parameter:: pi=3.141592654
      dimension f(nrmax,nthmax,nphimx),psi(nrmax),
     &   theta(nthmax),phi(nphimx),
     &   rootg(nrmax,nthmax)
      complex:: f,ansz,ans,anszm1,ansz1
      ans=0.0
      do k=1,nphi
      if(k.eq.2)ansz1=ansz
      if(k.gt.1)anszm1=ansz
      if(k.gt.1)dphi=phi(k)-phi(k-1)
      ansz=0.0
      do n=1,npsi-1
        dpsi=psi(n+1)-psi(n)
        do m=1,ntheta-1
          dtheta=theta(m+1)-theta(m)
          ansz=ansz+dpsi*dtheta
     &      *(rootg(n,m)*f(n,m,k)
     &      +rootg(n+1,m)*f(n+1,m,k)
     &      +rootg(n,m+1)*f(n,m+1,k)
     &      +rootg(n+1,m+1)*f(n+1,m+1,k))/4.
        end do
      end do
      m=ntheta
      dtheta=2.0*pi-theta(m)
      do n=1,npsi-1
        dpsi=psi(n+1)-psi(n)
        ansz=ansz+dpsi*dtheta
     &     *(rootg(n,m)*f(n,m,k)
     &     +rootg(n+1,m)*f(n+1,m,k)
     &     +rootg(n,1)*f(n,1,k)
     &     +rootg(n+1,1)*f(n+1,1,k))/4.0
      end do
      if(k.gt.1)ans=ans+(anszm1+ansz)/2.0*dphi
      end do
      k=nphi+1
      dphi=2.*pi-phi(k-1)
      ans=ans+(ansz+ansz1)/2.0*dphi
      return
      end
c
c***************************************************************************
c
      subroutine intplc(rgiv,thegiv,fl,nr,nth,f,nrmax,nthmax,
     &   dr,dtheta,rwall)
      dimension f(nrmax,nthmax)
      complex:: fl,f,a,b,c,d
      if(rgiv.ge.rwall)return
      r=abs(rgiv)
      theta=abs(thegiv)
      fl=0.0
      n=int(r/dr)+1
      m=int(theta/dtheta)+1
      if(n.ge.nr)return
      if(m.eq.nth)mp1=1
      if(m.lt.nth)mp1=m+1
      zeta=(r-(n-1)*dr)/dr
      eta=(theta-(m-1)*dtheta)/dtheta
      a=f(n,m)
      b=f(n+1,m)-f(n,m)
      c=f(n,mp1)-f(n,m)
      d=f(n+1,mp1)+f(n,m)-f(n+1,m)-f(n,mp1)
      fl=a+b*zeta+c*eta+d*zeta*eta
  300 format(1p6e12.4)
 1000 format(1p8e12.4)
  301 format(2i10)
c     call exit
      return
      end
c
c***************************************************************************
c
      subroutine intp3c(rgiv,thegiv,fl,nr,nth,f,nrmax,nthmax,
     &   nphimx,k,dr,dtheta,rwall)
      dimension f(nrmax,nthmax,nphimx)
      complex:: fl,f,a,b,c,d
      if(rgiv.ge.rwall)return
      r=abs(rgiv)
      theta=abs(thegiv)
      fl=0.0
      n=int(r/dr)+1
      m=int(theta/dtheta)+1
      if(n.ge.nr)return
      if(m.eq.nth)mp1=1
      if(m.lt.nth)mp1=m+1
      zeta=(r-(n-1)*dr)/dr
      eta=(theta-(m-1)*dtheta)/dtheta
      a=f(n,m,k)
      b=f(n+1,m,k)-f(n,m,k)
      c=f(n,mp1,k)-f(n,m,k)
      d=f(n+1,mp1,k)+f(n,m,k)-f(n+1,m,k)-f(n,mp1,k)
      fl=a+b*zeta+c*eta+d*zeta*eta
  300 format(1p6e12.4)
 1000 format(1p8e12.4)
  301 format(2i10)
c     call exit
      return
      end
c
c***************************************************************************
c
      subroutine intp3d(rgiv,thegiv,fl,nr,nth,f,nrmax,nthmax,
     &   nphimx,k,dr,dtheta,rwall)
      dimension f(nrmax,nthmax,nphimx)
      if(rgiv.ge.rwall)return
      r=abs(rgiv)
      theta=abs(thegiv)
      fl=0.0
      n=int(r/dr)+1
      m=int(theta/dtheta)+1
      if(n.ge.nr)return
      if(m.eq.nth)mp1=1
      if(m.lt.nth)mp1=m+1
      zeta=(r-(n-1)*dr)/dr
      eta=(theta-(m-1)*dtheta)/dtheta
      a=f(n,m,k)
      b=f(n+1,m,k)-f(n,m,k)
      c=f(n,mp1,k)-f(n,m,k)
      d=f(n+1,mp1,k)+f(n,m,k)-f(n+1,m,k)-f(n,mp1,k)
      fl=a+b*zeta+c*eta+d*zeta*eta
  300 format(1p6e12.4)
 1000 format(1p8e12.4)
  301 format(2i10)
c     call exit
      return
      end
c
c***************************************************************************
c
      subroutine intp2d(rgiv,thegiv,fl,nr,nth,f,nrmax,nthmax,
     &   nphimx,k,dr,dtheta,rwall)
      dimension f(nrmax,nthmax)
      if(rgiv.ge.rwall)return
      r=abs(rgiv)
      theta=abs(thegiv)
      fl=0.0
      n=int(r/dr)+1
      m=int(theta/dtheta)+1
      if(n.ge.nr)return
      if(m.eq.nth)mp1=1
      if(m.lt.nth)mp1=m+1
      zeta=(r-(n-1)*dr)/dr
      eta=(theta-(m-1)*dtheta)/dtheta
      a=f(n,m)
      b=f(n+1,m)-f(n,m)
      c=f(n,mp1)-f(n,m)
      d=f(n+1,mp1)+f(n,m)-f(n+1,m)-f(n,mp1)
      fl=a+b*zeta+c*eta+d*zeta*eta
  300 format(1p6e12.4)
 1000 format(1p8e12.4)
  301 format(2i10)
c     call exit
      return
      end
c
c***************************************************************************
c

      subroutine intplt_real(xgiv, ygiv, fout, nx, ny, f, nxmax, nymax,
     &   dx, dy)

      implicit none

      integer:: nx, ny, nxmax, nymax, n, m
      real:: x, y, xgiv, ygiv, dx, dy, zeta, eta, a, b, c, d, fout
      real:: f(nxmax, nymax)


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

  900 format(1p6e12.4)
  910 format(1p8e12.4)
  920 format(2i10)
      end



c
c***************************************************************************
c

      subroutine intplt_phi(xgiv, ygiv, fout, nx, ny, f, nxmax, nymax,
     ,   dx, dy)

      implicit none

      integer:: nx, ny, nxmax, nymax, n, m, mp1
      real:: x, y, xgiv, ygiv, dx, dy, zeta, eta, a, b, c, d, fout
      real:: f(nxmax, nymax)

      x = xgiv
      y = ygiv

      fout = 0.0

      n = int(x / dx) + 1
      m = int(y / dy) + 1

      if (n .ge. nx) return

      if (m .eq. ny) mp1 = 1
      if (m .lt. ny) mp1 = m + 1

      zeta = (x - (n - 1) * dx) / dx
      eta  = (y - (m - 1) * dy) / dy

      a = f(n, m)
      b = f(n+1 ,m) - f(n, m)
      c = f(n, mp1) - f(n, m)
      d = f(n+1, mp1) + f(n, m) - f(n+1, m) - f(n, mp1)

      fout = a + b * zeta + c * eta + d * zeta * eta

      return

  900 format(1p6e12.4)
  910 format(1p8e12.4)
  920 format(2i10)
      end

c
c***************************************************************************
c

      subroutine intp1d(rgiv,fl,nr,f,nrmax,dr,rwall)
      dimension f(nrmax)
      if(rgiv.gt.rwall)return
      r=abs(rgiv)
      fl=0.0
      n=int(r/dr)+1
      if(n.ge.nr)return
      zeta=(r-(n-1)*dr)/dr
      a=f(n)
      b=f(n+1)-f(n)
      fl=a+b*zeta
  300 format(1p6e12.4)
 1000 format(1p8e12.4)
  301 format(2i10)
      return
      end


c
c***************************************************************************
c


c Subroutine a1mnmx_dp
c----------------------------------------------------------------------------!
c Minimum and the maximum elements of a 1-d array.
c----------------------------------------------------------------------------!

      subroutine a1mnmx_dp(f, nxmax, nx, fmin, fmax)

      implicit none

      integer:: nx, i, nxmax
      real:: f(nxmax), fmin, fmax

      fmax=f(1)
      fmin=fmax
      do i=1, nx
            fmax=amax1(fmax, f(i))
            fmin=amin1(fmin, f(i))
c           write(6, *)fmax, fmin
      end do

      return
 2201 format(2i5,1p,8e12.4)
      end


c***************************************************************************
c


