
*----------------------------------------------------------------------
*     AORSA2D: version 20 (04-28-2010)
*     This version includes odd-order derivative in the conductivity
*     Special eqdsk_setup for Alcator-C lower hybrid
*        (eqdsk_special .ne. 0)
*     Does not change the sign of psi for abnormal eqdsk's
*     Uses simple but slower load and solve
*        - good out to the outer box wall
*----------------------------------------------------------------------


      program aorsa2dMain

      use size_mod

      use ql_myra_mod

      use wdot_mod
      use wdot_sum_mod


      use profile_mod

      use aorsa2din_mod
c      use zfun_hilbert

      implicit none

      integer, parameter :: real32=selected_real_kind(6,37)
      complex :: l2norm_res, l2norm_rhs
      integer :: ieer ! fft error flag
      integer :: n_theta_max, n_theta, n_psi_max, nt, k, nbessj
      integer :: n_upper, n_lower, m_upper, m_lower, n1
      integer :: jdisp, mdisp, idum, jdum, norbit
      integer :: pgopen, pgbeg, ier, np
      integer :: npoints, nlen_exitp, nlen_exitm, nlen_exit, n_len

      integer :: nphi1 = -50
      integer :: nphi2 =  50

      complex xjpxh, xjpyh, xjpzh

      CHARACTER(LEN=60) :: rowfmt

      real :: zetai_max, zetai_min, zetai_step, zetai, phil,
     &    dKdL_max, dKdL_min, dKdL_step, modEln, modEln2, cosarg, sinarg

      real :: zetai_max_table, zetai_min_table,
     &    dKdL_max_table, dKdL_min_table, psiorb, sin_antang, cos_antang
      real :: antang_rad,xn00,xn10,xn20,xn30, xn40, xn50, xn60, xn_slo0

c      CHARACTER(128) :: root

      character(16):: title
      character(16):: titll
      character(16):: titlr
      character(16):: titlb
      character(16):: titlx
      character(16):: titly
      character(16):: titlz
      character(16):: tityl
      character(16):: tityr


      real :: i_ant, z_load, ranval, nmin1, nmin2, nmin3, nmin4, nmin5
      real :: caprh,fn, fm, nperp, nperp_max, fnmax, fmmax, drho1, nmin6
      real :: Tmini, Tmini2, Tmini3, Tmini4, Tmini5, Tmini6, rho0, rhod1
      real :: Tmin1, Tmin2, Tmin3, Tmin4, Tmin5, Tmin6, Temin0
      real :: nminslo, ans, spa_cold, spa_slow, spa_fast, c_ehst, bmaxa
      real :: xmut2,eta0, ceta, damp, r1, xm1, c1, afit, yt, cfit, xnexp
      real :: xjtild, signkz
      real :: ydiff, ydiff_min, xkprl, xkprl_max, xkprl_min, xkprl_orig
      real :: dKdL_giv, abs_dKdL, sgn_dKdL, sgn_vprl, dKdL_giv_1
      real :: d2KdL2_giv
      real :: dupara, duperp, dvperp, dvpara, vpara_in
      real :: vperp_in, xmi_slo, ntilda_e_max
      real :: vce_mks_cql3d, vc1_mks_cql3d, vc2_mks_cql3d
      real :: vce_cgs_cql3d, vc1_cgs_cql3d, vc2_cgs_cql3d, vc3_cgs_cql3d
      real :: vc4_cgs_cql3d
      real :: bqlavg_i2_out, cqlavg_i2_out, eqlavg_i2_out, fqlavg_i2_out
      real :: diff_coef
      real :: xkalp, xkbet
      real :: xkprl_neg_max, xkprl_pos_min
      real :: sgn_kprl, akprl, gammab, y0, yb, fgam, xkprl_eff, descrim
      real :: xk_cutoff, term2, rmaxis, zmaxis, z0, ptot_wdot, xj
      real :: term1_coef, term2_coef, akprl_min

      integer :: jmid, nxdim, nydim, l, lmaxdim, n_ono, m_ono
      integer :: nnodex_loc, nx_overlap, nnodex_eff
      integer :: nnodey_loc, ny_overlap, nnodey_eff
      integer :: istart, ifinish, jstart, jfinish
     &  i_seg, j_seg, iseg1, iseg2, jseg1, jseg2, i0, j0

      integer :: nxmx,nymx,ieq, nlmx, nkldim1, nkldim2, nkl, nkl1, nkl2
      integer nmax, mmax
      real :: zeta_max, zeta_min, dzeta, dKdL, dkprldl, d2kprldl2
      real :: a0, a1, b1, psil, c0, xkprl_phi, xkprl_perp, psil0
      real :: rt_init, b0_init

      complex z0n, z1n, z2n
      integer :: incX, global_row, global_col
      integer :: ndfmax, lwork, lrwork,
     &    ninteg, nd, izoom1, izoom2, jzoom1, jzoom2, ndf, nmaxe, irnc
      integer :: nrow, ncol, norder

      integer :: mkdim1, mkdim2

      integer :: nkdim1, nkdim2, nkx1, nkx2, nldim, nldim3,
     &   nky1, nky2, iant

      real :: kperp_max, valfven, kalfven, xky_ant
      real ::  kperp_max_actual
      real :: kperprhoi1, kperprhoi2, kperprhoi3, kperprhoi4, kperprhoi5

      real :: eslow, qi_slo, xmi_sl, yant_max,
     &    xn_slo, a, epsa, epsa_potato, rmaxa, rmina

      real :: flimiter, parabola, dfdx, dfdy, caprmax, caprmin, rmax
      real :: caprl, xprimel, yprimel, zprimel

      real :: tmem, tsys, tio, ttotal, time0, time, cpu, dummy, second1

      real :: tmin, gflops, gflopsp, ops, teev, t1ev, t2ev, t3ev
      real :: t4ev, t5ev, t6ev

      real :: sqx, signb, gausspsi, dpsiant, theta_antr

      real :: dxsqx, dxxsqx, dxysqx, dysqx, dyysqx, UminPara, UmaxPara

      real :: vce_mks,vc1_mks,vc2_mks,vc3_mks, vc4_mks, vc5_mks, vc6_mks
      real :: vce_cgs,vc1_cgs,vc2_cgs,vc3_cgs, vc4_cgs, vc5_cgs, vc6_cgs

      integer ::  nrhomax, nthetamax, ndim
      integer, parameter ::
     &     ntable    = 1000,
     &     mtable    = 1000,
     &     norb_dim  = 4000,
     &     norb_dim2 = norb_dim*2

      real :: capr_xp(norb_dim), capz_xp(norb_dim), phin_xp(norb_dim)
      real :: capr_xm(norb_dim), capz_xm(norb_dim), phin_xm(norb_dim)
      real :: len_xm(norb_dim),  len_xp(norb_dim)
      real :: b_capr(norb_dim2),  c_capr(norb_dim2),  d_capr(norb_dim2)
      real :: b_capz(norb_dim2),  c_capz(norb_dim2),  d_capz(norb_dim2)
      real :: b_phin(norb_dim2),  c_phin(norb_dim2),  d_phin(norb_dim2)

      real :: len_x(norb_dim2),  capr_x(norb_dim2),
     &     capz_x(norb_dim2), phin_x(norb_dim2)

      real :: ispline


      complex z0_table1(ntable, mtable)
      complex z1_table1(ntable, mtable)
      complex z2_table1(ntable, mtable)
      real :: zetai_table(ntable), dKdL_table(mtable)

      parameter (n_theta_max = 200)
      parameter (n_psi_max = 200)

      parameter (nxmx = nmodesmax)
      parameter (nymx = mmodesmax)

      parameter (nlmx = 512)
      parameter (nkldim1 = -nlmx / 2)
      parameter (nkldim2 =  nlmx / 2)

      parameter (nrhomax = nmodesmax * 2)
      parameter (nthetamax = nmodesmax * 2)

      parameter (nkdim1 = - nmodesmax / 2)
      parameter (nkdim2 =   nmodesmax / 2)

      parameter (mkdim1 = - mmodesmax / 2)
      parameter (mkdim2 =   mmodesmax / 2)

      parameter (nldim  = nxmx * nymx)
      parameter (nldim3 = 3 * nldim)
      parameter (ndfmax = nldim3)

      parameter (lmaxdim = 99)

      parameter (ndim = 256)

      integer :: nsum

      real :: rhod(ndim), xneavg(ndim), tekev(ndim), tikev(ndim),
     &     omega(ndim), zeffavg(ndim), xn_beam(ndim), e_beam(ndim),
     &     drho_data, rhot, ti2kev(ndim), xn_maj(ndim), xn_min(ndim),
     &     xn_min_cgs(ndim), xn_maj_cgs(ndim), xneavg_cgs(ndim),
     &     xn_beam_cgs(ndim)

      real :: rhod2(ndim), xneavg2(ndim), tekev2(ndim), tikev2(ndim),
     &    xn_beam2(ndim), e_beam2(ndim), zeffavg2(ndim)

      integer, parameter :: n_psi_dim = 64

      integer :: n_psi

      real :: dens(n_psi_dim), rho_a(n_psi_dim)

      real, dimension(:),   allocatable :: UPERP, UPARA
      real, dimension(:),   allocatable :: VPERP, VPARA
      real, dimension(:),   allocatable :: UPERP_work, UPARA_work

      real, dimension(:,:), allocatable :: DFDUPERe, DFDUPARe

      real, dimension(:,:), allocatable :: DFDUPER1, DFDUPAR1
      real, dimension(:,:), allocatable :: DFDUPER2, DFDUPAR2
      real, dimension(:,:), allocatable :: DFDUPER3, DFDUPAR3
      real, dimension(:,:), allocatable :: DFDUPER4, DFDUPAR4
      real, dimension(:,:), allocatable :: DFDUPER5, DFDUPAR5
      real, dimension(:,:), allocatable :: DFDUPER6, DFDUPAR6

      real, dimension(:),   allocatable :: FPERPe
      real, dimension(:),   allocatable :: FPERP1
      real, dimension(:),   allocatable :: FPERP2
      real, dimension(:),   allocatable :: FPERP3
      real, dimension(:),   allocatable :: FPERP4
      real, dimension(:),   allocatable :: FPERP5
      real, dimension(:),   allocatable :: FPERP6


      real, dimension(:,:), allocatable :: f

      real, dimension(:,:,:), allocatable :: fe_cql_cart
      real, dimension(:,:,:), allocatable :: dfe_cql_uprp
      real, dimension(:,:,:), allocatable :: dfe_cql_uprl

      real, dimension(:,:,:), allocatable :: f1_cql_cart
      real, dimension(:,:,:), allocatable :: df1_cql_uprp
      real, dimension(:,:,:), allocatable :: df1_cql_uprl

      real, dimension(:,:,:), allocatable :: f2_cql_cart
      real, dimension(:,:,:), allocatable :: df2_cql_uprp
      real, dimension(:,:,:), allocatable :: df2_cql_uprl

      real, dimension(:,:,:), allocatable :: f3_cql_cart
      real, dimension(:,:,:), allocatable :: df3_cql_uprp
      real, dimension(:,:,:), allocatable :: df3_cql_uprl

      real, dimension(:,:,:), allocatable :: f4_cql_cart
      real, dimension(:,:,:), allocatable :: df4_cql_uprp
      real, dimension(:,:,:), allocatable :: df4_cql_uprl

      real, dimension(:,:,:), allocatable :: f5_cql_cart
      real, dimension(:,:,:), allocatable :: df5_cql_uprp
      real, dimension(:,:,:), allocatable :: df5_cql_uprl

      real, dimension(:,:,:), allocatable :: f6_cql_cart
      real, dimension(:,:,:), allocatable :: df6_cql_uprp
      real, dimension(:,:,:), allocatable :: df6_cql_uprl

      real, dimension(:,:,:), allocatable :: bqlavg_e
      real, dimension(:,:,:), allocatable :: cqlavg_e
      real, dimension(:,:,:), allocatable :: eqlavg_e
      real, dimension(:,:,:), allocatable :: fqlavg_e

      real, dimension(:,:,:), allocatable :: bqlavg_i1
      real, dimension(:,:,:), allocatable :: cqlavg_i1
      real, dimension(:,:,:), allocatable :: eqlavg_i1
      real, dimension(:,:,:), allocatable :: fqlavg_i1

      real, dimension(:,:,:), allocatable :: bqlavg_i2
      real, dimension(:,:,:), allocatable :: cqlavg_i2
      real, dimension(:,:,:), allocatable :: eqlavg_i2
      real, dimension(:,:,:), allocatable :: fqlavg_i2

      real, dimension(:,:,:), allocatable :: bqlavg_i3
      real, dimension(:,:,:), allocatable :: cqlavg_i3
      real, dimension(:,:,:), allocatable :: eqlavg_i3
      real, dimension(:,:,:), allocatable :: fqlavg_i3

      real, dimension(:,:,:), allocatable :: bqlavg_i4
      real, dimension(:,:,:), allocatable :: cqlavg_i4
      real, dimension(:,:,:), allocatable :: eqlavg_i4
      real, dimension(:,:,:), allocatable :: fqlavg_i4

      real, dimension(:,:,:), allocatable :: bqlavg_i5
      real, dimension(:,:,:), allocatable :: cqlavg_i5
      real, dimension(:,:,:), allocatable :: eqlavg_i5
      real, dimension(:,:,:), allocatable :: fqlavg_i5

      real, dimension(:,:,:), allocatable :: bqlavg_i6
      real, dimension(:,:,:), allocatable :: cqlavg_i6
      real, dimension(:,:,:), allocatable :: eqlavg_i6
      real, dimension(:,:,:), allocatable :: fqlavg_i6

      real, dimension(:,:,:), allocatable :: bqlavg_work
      real, dimension(:,:,:), allocatable :: cqlavg_work
      real, dimension(:,:,:), allocatable :: eqlavg_work
      real, dimension(:,:,:), allocatable :: fqlavg_work

      integer :: i_uperp, i_upara, i_psi, i_psi_eq

      real :: factl_sav(0:lmaxdim), factrl

      character(1):: trans

      integer :: info
      complex b, zi, cexpkxky
      complex alpha

      real :: norm2, norm2_max, alphae, alpha1

      complex fdk, fek, ffk, fgk, fak, fpk, frk, fqk, fsk, detd

      real :: xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2),
     &     xkperp, xketa, xkperp2, cs, xkperp2_eicw, diff,
     &     ptot, ptot_core, ptot_edge, pcito2, pcrto2, powtot, pscale,
     &     cosalp, sinalp, t1, gaussian, psimax, fedge, fcore,
     &     xnprl, rhoant, gaussantx, shapey,
     &     dthetant, gaussantth, gaussiant, gaussthetap,
     &     constantth, abstheta0, xnustar, thetamax, dtheta

      real :: wdot_edge, wdot_core, f_wdot_edge, f_wdot_core

      real :: t1e, t1i1, t1i2, t1i3, t1i4, t1i5, t1i6

      real :: gaussian_ne,  gaussian_te,  gaussian_ti1, gaussian_ti2,
     &     gaussian_ti3, gaussian_ti4, gaussian_ti5, gaussian_ti6


      real :: shape, shapen,  shapen2,  shapen3, shapen_slo,
     &   shapen4, shapen5, shapen6,
     &   shapete, shapeti, shapeti2, shapeti3,
     &   shapeti4, shapeti5, shapeti6

      real :: frho, dxfrho, dxxfrho, dxyfrho, dyfrho, dyyfrho

c     complex xx(nkdim1 : nkdim2, 1 : nxmx),
c     &        yy(mkdim1 : mkdim2, 1 : nymx)
c      complex xx_inv(nkdim1 : nkdim2, 1 : nxmx),
c     &        yy_inv(mkdim1 : mkdim2, 1 : nymx)

      complex, dimension(:,:), allocatable :: xx, yy, xx_inv, yy_inv

      real :: dldb_tot12(nxmx, nymx)

      real :: theta_(n_theta_max, n_psi_max)
      integer :: n_theta_(n_psi_max)

      integer :: ilayer,nlayer
      integer :: mask(nxmx, nymx), mask2(nxmx,nymx)

*     ---------------------
*     2D allocatable arrays
*     ---------------------

      complex, dimension(:,:), allocatable :: ealphak
      complex, dimension(:,:), allocatable :: ebetak
      complex, dimension(:,:), allocatable :: ebk

      real, dimension(:,:), allocatable :: ealphakmod
      real, dimension(:,:), allocatable :: ebetakmod
      real, dimension(:,:), allocatable :: ebkmod

      complex, dimension(:,:), allocatable :: fdksav
      complex, dimension(:,:), allocatable :: feksav
      complex, dimension(:,:), allocatable :: ffksav
      complex, dimension(:,:), allocatable :: fgksav
      complex, dimension(:,:), allocatable :: faksav
      complex, dimension(:,:), allocatable :: fpksav
      complex, dimension(:,:), allocatable :: frksav
      complex, dimension(:,:), allocatable :: fqksav
      complex, dimension(:,:), allocatable :: fsksav

      complex, dimension(:,:), allocatable :: fdksav2d
      complex, dimension(:,:), allocatable :: feksav2d
      complex, dimension(:,:), allocatable :: ffksav2d
      complex, dimension(:,:), allocatable :: fgksav2d
      complex, dimension(:,:), allocatable :: faksav2d
      complex, dimension(:,:), allocatable :: fpksav2d
      complex, dimension(:,:), allocatable :: frksav2d
      complex, dimension(:,:), allocatable :: fqksav2d
      complex, dimension(:,:), allocatable :: fsksav2d

      real, dimension(:,:), allocatable :: capd
      real, dimension(:,:), allocatable :: capd2d

      real, dimension(:,:), allocatable :: xkb
      real, dimension(:,:), allocatable :: xkb2d

      integer :: nnodex, nnodey, i, j,
     &    jequat, iflag, liw, lw, nrhs, icenter
      integer :: nnodelb, nnodelb0, nl
      real :: lbmax, dlb, xl, yl, zl

      integer :: nnodexm1, nnodeym1
      integer :: nnodex_p, nnodey_p
      integer :: nnoderho, mnodetheta
      integer :: n, m, nphi, ndata, ndata2
      real :: xant, delta, xnurf, psio, psimag, psi_tor_max
      real :: telimj, tilimj, ti2limj, ti3limj,
     &           ti4limj, ti5limj, ti6limj

      real :: rhoplasm, reomg1, reomg2, reomg3,
     &    reomg4, reomg5, reomg6,
     &    r0, xiota0, rholim,
     &    deltay, abs_deltay, delta_theta,
     &    xwleft, xwright, psi_lim, psi1
      real :: xkphi(nxmx), xktau, xkrho, xkphi0, xkphia, xnphi, rlim
      real :: xkthrho, wphase, vsound, domgk, xkthdx, omgestar, rhoi10,
     &   v0i, vthi10, vthe, vthi, vphase, rhoi1overl, rnz, eta2,
     &   xk0, shearedge, eta1, xn1, xmi1, xmh, xme, qi1, xmi2,
     &   xmi3, t0i, t0i2, t0i3, t0e, q, teedge, clight, xmu0, xlnlam,
     &   omgci10, omgrf, xmax, qe, qi2, pi, twopi, eps0, qi3,
     &   costh, sinth, radius, rnx,  rny, rnphi
      real :: ymaxp, yminp, xminp, xmaxp
      real :: rhoi20, vthi20, omgci20
      real :: rhoi30, vthi30, omgci30
      real :: rhoi40, vthi40, omgci40
      real :: rhoi50, vthi50, omgci50
      real :: xmi4, t0i4, qi4
      real :: xmi5, t0i5, qi5
      real :: xmi6, t0i6, qi6
      real :: xjantx, xjanty, xjantz, xjant
      complex :: djxdx, djydy

      real :: xprimec(nxmx), caprc(nxmx), xcourse(nxmx), capr(nxmx),
     &   xprime(nxmx), x(nxmx), dx, dxc, drdx(nxmx), d2rdx2(nxmx)
      real :: dbxdx, dbxdy, dbydx, dbydy, dbphidx, dbphidy,
     &   d2bxdx2, d2bydx2, d2bzdx2,
     &   d2bxdxy, d2bydxy, d2bzdxy,
     &   d2bxdy2, d2bydy2, d2bzdy2,
     &   d2bphirdx2, d2bphirdxy, dbphirdx,
     &   bhatgradbx, bhatgradby, bhatgradbz,
     &   bhatgrad2bx, bhatgrad2by, bhatgrad2bz

      complex :: dxdkwxx, dxdkwxy, dxdkwxz,
     &        dxdkwyx, dxdkwyy, dxdkwyz,
     &        dxdkwzx, dxdkwzy, dxdkwzz

      complex :: dydkwxx, dydkwxy, dydkwxz,
     &        dydkwyx, dydkwyy, dydkwyz,
     &        dydkwzx, dydkwzy, dydkwzz

      complex :: dxdkwxxp, dxdkwxyp, dxdkwxzp,
     &        dxdkwyxp, dxdkwyyp, dxdkwyzp,
     &        dxdkwzxp, dxdkwzyp, dxdkwzzp

      complex :: dydkwxxp, dydkwxyp, dydkwxzp,
     &        dydkwyxp, dydkwyyp, dydkwyzp,
     &        dydkwzxp, dydkwzyp, dydkwzzp

      real :: capr_flux(nrhomax, nthetamax),
     &   capz_flux(nrhomax, nthetamax),xgiv, ygiv
      complex :: eplus_flux(nrhomax, nthetamax), fout
      complex :: eminus_flux(nrhomax, nthetamax)
      real :: bmod_flux(nrhomax, nthetamax)

      real :: rhon(nrhomax), thetam(nthetamax), rho_in, theta_in,
     &   capr_out, capz_out,
     &   wdoti1avg(nrhomax), wdoti2avg(nrhomax),
     &   wdoti3avg(nrhomax), wdoti4avg(nrhomax),
     &   dvol(nrhomax), darea(nrhomax), rn(nrhomax), rh(nrhomax),
     &   dn(nrhomax), dh(nrhomax), dr, an(nrhomax), bn(nrhomax),
     &   cn(nrhomax), yn(nrhomax), epsn(nrhomax),
     &   wdoteavg(nrhomax), wdotavg(nrhomax), drho, dpsi,
     &   capr_bpol_mid(nrhomax), capr_bpol_midavg(nrhomax),
     &   bmod_midavg(nrhomax), fvol(nrhomax), dvolp(nrhomax),
     &   volume(nrhomax)
      real :: zeffavg1(nrhomax)

      real :: ank(nrhomax), bnk(nrhomax), cnk(nrhomax), ynk(nrhomax)


      real ::  wdote_ql(nrhomax), wdoti1_ql(nrhomax), wdoti2_ql(nrhomax)
      real ::  wdoti3_ql(nrhomax),wdoti4_ql(nrhomax), wdoti5_ql(nrhomax)
      real ::  wdoti6_ql(nrhomax)

      real :: wdoti5avg(nrhomax), wdoti6avg(nrhomax)

      real :: xnavg(nrhomax), xn1avg(nrhomax), xn2avg(nrhomax),
     &     xn3avg(nrhomax), xna_sloavg(nrhomax)
      real :: xn4avg(nrhomax), xn5avg(nrhomax), xn6avg(nrhomax)


      real :: xkteavg(nrhomax), xktiavg(nrhomax),  xkti2avg(nrhomax),
     &     xkti3avg(nrhomax)
      real :: xkti4avg(nrhomax), xkti5avg(nrhomax), xkti6avg(nrhomax)
      real :: dldbavg(nrhomax)


      real :: fyavg(nrhomax), qhat, omgte, omgti,
     &     vthe0, vthi0, xnuee, xnuii, xnu7omg

      real :: redotj1avg(nrhomax), redotj2avg(nrhomax),
     &     redotjeavg(nrhomax), redotj3avg(nrhomax),
     &     redotj4avg(nrhomax),
     &     redotj5avg(nrhomax), redotj6avg(nrhomax),
     &     redotjsavg(nrhomax),
     &     redotjtavg(nrhomax), xjprlavg(nrhomax)

      real :: xjprl_int(nrhomax), fyp_int(nrhomax), fz0_int(nrhomax)

      real :: wdote_ql_int(nrhomax), wdoti1_ql_int(nrhomax),
     &                            wdoti2_ql_int(nrhomax),
     &                            wdoti3_ql_int(nrhomax),
     &                            wdoti4_ql_int(nrhomax),
     &                            wdoti5_ql_int(nrhomax),
     &                            wdoti6_ql_int(nrhomax)

      real :: wdoteavg_int(nrhomax), wdoti1avg_int(nrhomax),
     &                            wdoti2avg_int(nrhomax),
     &                            wdoti3avg_int(nrhomax),
     &                            wdoti4avg_int(nrhomax),
     &                            wdoti5avg_int(nrhomax),
     &                            wdoti6avg_int(nrhomax)


      real :: redotjeavg_int(nrhomax), redotj1avg_int(nrhomax),
     &                              redotj2avg_int(nrhomax),
     &                              redotj3avg_int(nrhomax),
     &                              redotj4avg_int(nrhomax),
     &                              redotj5avg_int(nrhomax),
     &                              redotj6avg_int(nrhomax)

      real :: fypavg(nrhomax), fypi1avg(nrhomax), fypi2avg(nrhomax),
     &     fypeavg(nrhomax), fypi3avg(nrhomax)

      real :: fypi4avg(nrhomax), fypi5avg(nrhomax), fypi6avg(nrhomax)


      real :: fz0avg(nrhomax), fz0i1avg(nrhomax), fz0i2avg(nrhomax),
     &     fz0eavg(nrhomax), fz0i3avg(nrhomax)

      real :: fz0i4avg(nrhomax), fz0i5avg(nrhomax), fz0i6avg(nrhomax)


      real :: capr_fzeta_avg(nrhomax), capr2_avg(nrhomax)
      real :: gpsi_avg(nrhomax), xjhat(nrhomax), dgdpsi_avg(nrhomax)
      real :: kpsi_avg(nrhomax), xkhat(nrhomax)

      real :: dpdpsi(nrhomax)
      real :: epsi_avg(nrhomax), dedpsi_avg(nrhomax)
      real :: epsig_avg(nrhomax), epsip_avg(nrhomax), epsir_avg(nrhomax)

      real :: dgdpsi(nxmx, nymx), dedpsi(nxmx, nymx)
      real :: kimag_slow(nxmx), kimag_fast(nxmx), kimag_cold(nxmx)
      real :: vyavg(nrhomax),  vyi1avg(nrhomax),  vyi2avg(nrhomax)
      real :: muhat_avg(nrhomax), nu_star_avg(nrhomax),
     &     pressiavg(nrhomax),
     &     rhom1avg(nrhomax), qsafety_avg(nrhomax), ipsi_avg(nrhomax),
     &     gamma_avg(nrhomax), beta_avg(nrhomax), pitch_avg(nrhomax)

      real :: taup_b, taup_pl, taup_ps, taup, tauii, taup_inv, angle

      real :: dvydrho(nrhomax), psi_dim_avg(nrhomax),
     &     fpsi1_avg(nrhomax), ftheta1_avg(nrhomax),
     &     gradprlb_avg(nrhomax),
     &     bdotf_avg(nrhomax), bmod2_avg(nrhomax)

      real :: drhodx, drhodxx, drhody, drhodyy, gradrho

      real :: dxpsi, dxxpsi, dxypsi, dypsi, dyypsi

      real, dimension(:,:), allocatable :: psi, rho, theta,
     &     rhohatx, rhohaty, theta0, psi_dim, bx, by, bz,
     &     btau, bzeta, dxdth, dzdth, xntau, xkte, xkti,
     &     xkti2, xkti3, xkti4, xkti5, xkti6, xn1a, xnea, xn2a,
     &     xn3a, bmod, omgce, xn4a, xn5a, xn6a,
     &     omgci1, omgci2, omgci3, omgci4, omgci5, omgci6,
     &     omgpe2, bpol, capr_bpol, omgp12, omgp22, omgp32,
     &     omgp42, omgp52, omgp62, xiota, qsafety, bmod_mid,
     &     capr_bpol_mid2, rho_tor2d, psi_tor2d, omglh, pitch,
     &     xnuomea, xnuomia

      real :: reomgpi, reomglh

      real :: dxbx, dxxbx, dxybx, dybx, dyybx
      real :: dxby, dxxby, dxyby, dyby, dyyby
      real :: dxbz, dxxbz, dxybz, dybz, dyybz
      real :: dxbmod, dxxbmod, dxybmod, dybmod, dyybmod
      real :: rhome, rhomi1, rhomi2, rhomi3, rhomslo, prod
      real :: rhomi4, rhomi5, rhomi6, omglh2i
      real :: presse, pressi1, pressi2, pressi3,
     &                                 pressi4, pressi5, pressi6
      real :: dbdx, dbdy, d2bdx2, d2bdy2
      real :: dpsidx, dpsidy, d2psidx2, d2psidy2, gradpsi, Ix, Iy, Iz
      real :: dxbxn, dybxn, dxxbxn, dyybxn, dxybxn
      real :: dxbyn, dybyn, dxxbyn, dyybyn, dxybyn
      real :: dxbzn, dybzn, dxxbzn, dyybzn, dxybzn
      complex :: dezdx, dezdy, deydx, dexdy
      real :: re_acold, rho_giv, the_giv, isq2

      real, dimension(:,:), allocatable :: dxuxx, dxuxy, dxuxz,
     &   dxuyx, dxuyy, dxuyz, dxuzx, dxuzy, dxuzz,
     &   dyuxx, dyuxy, dyuxz, dyuyx, dyuyy, dyuyz,
     &   dyuzx, dyuzy, dyuzz, dyyuxx, dyyuxy, dyyuxz,
     &   dyyuyx, dyyuyy, dyyuyz, dyyuzx, dyyuzy, dyyuzz,
     &   dxyuxx, dxyuxy, dxyuxz, dxyuyx, dxyuyy, dxyuyz,
     &   dxyuzx, dxyuzy, dxyuzz, dxxuxx, dxxuxy, dxxuxz,
     &   dxxuyx, dxxuyy, dxxuyz, dxxuzx, dxxuzy, dxxuzz,
     &   spx, spy, spz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz,
     &   bmod2, rhom1, bxn, byn, bzn, pressi, muhat, nu_star,
     &   ipsi, xna_slo, gradprlb, omgp2_slo, omgci_slo, zeff,
     &   rho_pol2d, psi_pol2d, drhodr, drhodz,
     &   dthedr, dthedz, reomg1a, reomg2a, reomg3a, xjx, xjy, xjz,
     &   ntilda_e_real, reomglha

      real :: cur_mod, cur_mod_max, npar
      integer :: i_cur_max, j_cur_max

      real, dimension(:,:), allocatable :: kappa_hatx, kappa_haty,
     &   kappa_hatz, rc, bhatx, bhaty, bhatz

      real :: xjy_max, xjx_max, xjz_max, mod_kappa
      real :: xjy_min, xjx_min, xjz_min

      complex, dimension(:,:), allocatable :: acold, bcold, ccold,
     &   xkperp_cold, xkperp_cold2, ex, ey, ez, bxwave, bywave, bzwave,
     &   eplus_flux_plot, eminus_flux_plot, xkperp_flux_plot,
     &   ealpha, ebeta, eb, eplus, eminus, ealphakx, ebetakx,
     &   ealphaky, ebetaky, ebkx, ebky, rho_ant, rho_pla, ntilda_e,
     &   xb, xc, xd

      complex, dimension(:,:), allocatable :: xkperp2_slow,
     &        xkperp2_fast, P_a
      real, dimension(:,:), allocatable :: xkprl_a

      complex :: xkperp_flux(nrhomax, nthetamax)

c      complex :: uxxk(nkdim1 : nkdim2, mkdim1 : mkdim2),
c     &        uxyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
c     &        uxzk(nkdim1 : nkdim2, mkdim1 : mkdim2),
c     &        uyxk(nkdim1 : nkdim2, mkdim1 : mkdim2),
c     &        uyyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
c     &        uyzk(nkdim1 : nkdim2, mkdim1 : mkdim2),
c     &        uzxk(nkdim1 : nkdim2, mkdim1 : mkdim2),
c     &        uzyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
c     &        uzzk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      complex, dimension(:,:), allocatable :: uxxk, uxyk,
     &  uxzk, uyxk, uyyk, uyzk, uzxk, uzyk, uzzk


      complex, dimension(:,:), allocatable :: xjpxe, xjpye, xjpze,
     &   xjpx1, xjpy1, xjpz1,xjpx2, xjpy2, xjpz2, xjpx3,
     &   xjpy3, xjpz3, xjpx4, xjpy4, xjpz4, xjpx5, xjpy5,
     &   xjpz5, xjpx6, xjpy6, xjpz6


      complex, dimension(:,:), allocatable :: xjpx_ehst, xjpy_ehst,
     &   xjpz_ehst, xj_slox, xj_sloy, xj_sloz,
     &   xjpx, xjpy, xjpz, xjpx_lab, xjpy_lab, xjpz_lab,  pc

      complex, dimension(:,:), allocatable ::  xjpxe_lab, xjpye_lab,
     &      xjpze_lab

      real :: fye, fyi1, fyi2, fyi3, fy, fyi4, fyi5, fyi6
      real :: denom, xjtot, fyptot, fz0tot
      real :: fxtot_ant, fytot_ant, fztot_ant
      real :: fxtot_plasma, fytot_plasma, fztot_plasma


      real, dimension(:,:), allocatable :: wdote, wdoti1,
     &   wdoti2, wdoti3, wdot,
     &   wdoti4, wdoti5, wdoti6, fype, fypi1,
     &   fypi2, fypi3, fyp, fypi4, fypi5, fypi6, fpol0e, fpol0i1,
     &   fpol0i2, fpol0i3, fpol0i4, fpol0i5, fpol0i6, fx0e, fx0i1,
     &   fx0i2, fx0i3,
     &   fx0i4, fx0i5, fx0i6, fy0e , fy0i1, fy0i2, fy0i3,
     &   fy0i4, fy0i5, fy0i6, fz0e , fz0i1,
     &   fz0i2, fz0i3,
     &   fz0i4, fz0i5, fz0i6, bdotf, capr_fzeta, capr2,
     &   uzeta, utheta, gpsi, omgexb, jhat, kpsi, epsi, fprl0e,
     &   fprl0i1,
     &   fprl0i2, fprl0i3, fprl0i4, fprl0i5, fprl0i6, fx0, fy0, fz0,
     &   fpsi0, fpsi1, ftheta0, ftheta1, fpol1e , fpol1i1,
     &   fpol1i2, fpol1i3, fpol1i4, fpol1i5, fpol1i6, dfpol1e,
     &   dfpol1i1, dfpol1i2, dfpol1i3, dfpol1i4, dfpol1i5, dfpol1i6,
     &   redotj, divq,
     &   pcre, pcim, fx_ant, fy_ant, fz_ant, fx_pla, fy_pla, fz_pla,
     &   rekxdotj, rekydotj,
     &   rekzdotj, xjprl, redotj1, redotj2, redotj3, redotj4, redotj5,
     &   redotj6, redotje, redotjt, redotji, redotjs, redotj_ehst



      real :: pcedotj1, pcedotje, pcedotj2, pcedotjt, pcedotj3, pcedotjs
      real :: pcedotj4, pcedotj5, pcedotj6

      real :: pedotj1, pedotje, pedotj2, pedotjt, pedotj3, pedotjs
      real :: pedotj4, pedotj5, pedotj6

      real :: p, pi1, pi2, pit, pi3, pe, pt, pi4, pi5, pi6

      real :: p_ql, pi1_ql, pi2_ql, pit_ql, pi3_ql, pe_ql,
     &            pt_ql, pi4_ql, pi5_ql, pi6_ql

      real :: pcti1, pcte, pcti2, pctt, pcti3, pcti4, pcti5, pcti6

      real :: pcti1_ql, pcte_ql, pcti2_ql, pctt_ql, pcti3_ql, pcti4_ql,
     &     pcti5_ql, pcti6_ql

      complex :: ekl, xintc

      complex :: sigexx, sigexy, sigexz,
     &        sigeyx, sigeyy, sigeyz,
     &        sigezx, sigezy, sigezz

      complex :: sige0, sige1, sige2,
     &        sige3, sige4, sige5

      complex :: sige0_new, sige1_new, sige2_new,
     &        sige3_new, sige4_new, sige5_new

      complex :: sigexx_new, sigexy_new, sigexz_new,
     &        sigeyx_new, sigeyy_new, sigeyz_new,
     &        sigezx_new, sigezy_new, sigezz_new

      complex :: sig1xx, sig1xy, sig1xz,
     &        sig1yx, sig1yy, sig1yz,
     &        sig1zx, sig1zy, sig1zz

      complex :: sig2xx, sig2xy, sig2xz,
     &        sig2yx, sig2yy, sig2yz,
     &        sig2zx, sig2zy, sig2zz

      complex :: sig3xx, sig3xy, sig3xz,
     &        sig3yx, sig3yy, sig3yz,
     &        sig3zx, sig3zy, sig3zz

      complex :: sig4xx, sig4xy, sig4xz,
     &        sig4yx, sig4yy, sig4yz,
     &        sig4zx, sig4zy, sig4zz

      complex :: sig5xx, sig5xy, sig5xz,
     &        sig5yx, sig5yy, sig5yz,
     &        sig5zx, sig5zy, sig5zz

      complex :: sig6xx, sig6xy, sig6xz,
     &        sig6yx, sig6yy, sig6yz,
     &        sig6zx, sig6zy, sig6zz

      complex :: sigsloxx, sigsloxy, sigsloxz,
     &        sigsloyx, sigsloyy, sigsloyz,
     &        sigslozx, sigslozy, sigslozz

      complex
     &     sigxx, sigxy, sigxz,
     &     sigyx, sigyy, sigyz,
     &     sigzx, sigzy, sigzz

      complex
     &     xkxx, xkxy, xkxz,
     &     xkyx, xkyy, xkyz,
     &     xkzx, xkzy, xkzz
      complex :: scap, capk_perp, capk_x, left, right

      complex :: xk1xx, xk1xy, xk1xz,
     &        xk1yx, xk1yy, xk1yz,
     &        xk1zx, xk1zy, xk1zz

      complex :: xk2xx, xk2xy, xk2xz,
     &        xk2yx, xk2yy, xk2yz,
     &        xk2zx, xk2zy, xk2zz


      complex
     &     dxx, dxy, dxz,
     &     dyx, dyy, dyz,
     &     dzx, dzy, dzz

      complex
     &     bxx, bxy, bxz,
     &     byx, byy, byz,
     &     bzx, bzy, bzz


      real :: yprimec(nymx), ycourse(nymx),
     &     yprime(nymx), y(nymx), dy, dyc

      real :: lb(nlmx), lbprime(nlmx), capr_lb(nlmx), capz_lb(nlmx),
     &                                              phin_lb(nlmx)
      complex :: El(nlmx),Elsum(nlmx),workl(nlmx), Eln(nkldim1:nkldim2)
      complex :: Eln_chk(nkldim1 : nkldim2)
      complex :: Esum

      complex :: Jxxsum, Jxysum, Jxzsum,
     &        Jyxsum, Jyysum, Jyzsum,
     &        Jzxsum, Jzysum, Jzzsum

      complex :: J0sum, J1sum, J2sum,
     &        J3sum, J4sum, J5sum

      real :: xklsav(nkldim1 : nkldim2)
      real :: xklsavp(nlmx), Elnp(nlmx)
      real :: arg



* ------------------------------
* storage for parallel scalapack
* ------------------------------
      INTEGER ::            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER ::            CTXT_, M_, N_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6 )
      INTEGER ::            RSRC_, CSRC_, LLD_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

      integer ::            p_amat_dim,p_brhs_dim,p_ipiv_dim
      integer ::            p_amat_size


      logical use_cur_mod
      parameter(use_cur_mod=.true.)

      complex, allocatable :: p_amat(:), p_amat_copy(:)
      integer, allocatable :: p_ipiv(:)
      complex, allocatable :: p_brhs(:), p_brhs_copy(:)

      real :: err_max,err,b_nrm2,err_nrm2
      complex, dimension(ndfmax) ::brhs2,brhsk,brhs_tmp
      integer :: icnc
      logical :: ismine, ismine1

      integer :: numroc
      external numroc

      integer :: lld,nrow_local,ncol_local

      character(4):: suffix

      integer :: mb,nb,myid,nproc,wantnproc,myrow,mycol
      integer :: icontxt, wantnprocs
      integer :: lrindx,lcindx,rsrc,csrc,ipos,ii,jj
      integer :: desc_amat(dlen_), desc_brhs(dlen_)

      integer :: rsrc1, csrc1, irnc1, icnc1

!efd-begin
      integer, parameter :: idebug=0
      integer :: ia,ja, iastart,iaend, jastart,jaend
      integer :: ia_first,ja_first,   iproc
c      integer, allocatable, dimension(:) ::
c     &         itable,jtable,ntable,mtable
      integer :: lrindx_base,lcindx_base
      integer :: lrindx1,lcindx1
      logical :: isok, found
      integer, parameter :: undefined=-987654321

*       ---------------------------------------------------
*       local variables for transformation to real :: space
*       ---------------------------------------------------
        integer :: ib,jb
        integer :: Locp, Locq, irow, icol, nia,nja
        integer :: niastart, niaend
        integer :: num_org, num_new, niasize
        integer :: org_nrow, org_ncol, mmb, nnb
        integer :: new_nrow, new_ncol
        integer :: num_inside, mm, ni, nn, iia, jja, niu, miu

        logical, dimension(:,:), allocatable :: is_inside
c       integer, dimension(:), allocatable :: new_to_org
c       integer, dimension(:,:), allocatable :: niabegin_all
c       integer, dimension(:,:), allocatable :: isize_all
c       integer, dimension(:,:,:), allocatable :: descBtmp_all
c       integer, dimension(:,:,:), allocatable :: descbrhs_all
        integer, dimension(DLEN_) :: descBtmp, descbrhs

        real :: rho_thresh



        complex :: aij, beta
c       complex, dimension(:,:), allocatable :: brhs
        real,    dimension(:,:), allocatable :: work
        complex, dimension(:,:), allocatable :: workc

c       complex, dimension(:,:), allocatable :: Btmp

        complex, dimension(ndfmax) :: row,rowk

        complex, dimension(ndfmax) :: workn
        real, dimension(ndfmax) :: rwork
        real :: pzlange, anorm, rcond

        logical, parameter :: use_fft = .true.


*       ---------------------------------------------------------
*       further permutation to avoid load imbalance due to
*       assigning cells in resonance region to the same processor
*       ---------------------------------------------------------
        logical, parameter :: need_shuffle = .true.
        integer :: ip, nrow3,ilen
        real, allocatable, dimension(:) :: dtemp2
        integer, allocatable, dimension(:) :: itemp2, iperm

      common/upcom/akprl_min


      integer :: indxg2p,indxg2l,indxl2g
      external indxg2p,indxg2l,indxl2g

! DLG - check to make sure the problem size you are running is <= the
! compile time parameters, i.e.,

      if (nmodesx .gt. nmodesmax .or. nmodesy .gt. mmodesmax) then
        stop 'ERROR - nmodesx and nmodesy must be <= nmodesmax and '//
     &        'mmodesmax respectively (see src/size_mode.f90)'
      end if

*     ---------------
*     allocate arrays
*     ---------------

      allocate ( xx(nkdim1 : nkdim2, 1 : nxmx),
     &           yy(mkdim1 : mkdim2, 1 : nymx) )
      allocate ( xx_inv(nkdim1 : nkdim2, 1 : nxmx),
     &           yy_inv(mkdim1 : mkdim2, 1 : nymx) )

      allocate (work(nxmx, nymx))
      allocate (workc(nxmx, nymx))

      allocate ( psi(nxmx, nymx),  rho(nxmx, nymx), theta(nxmx, nymx),
     &   rhohatx(nxmx, nymx), rhohaty(nxmx, nymx), theta0(nxmx, nymx),
     &   psi_dim(nxmx, nymx), bx(nxmx, nymx), by(nxmx, nymx),
     &   bz(nxmx, nymx), btau(nxmx, nymx), bzeta(nxmx, nymx),
     &   dxdth(nxmx, nymx),  dzdth(nxmx, nymx), xntau(nxmx, nymx),
     &   xkte(nxmx, nymx), xkti(nxmx, nymx), xkti2(nxmx, nymx),
     &   xkti3(nxmx, nymx), xkti4(nxmx, nymx), xkti5(nxmx, nymx),
     &   xkti6(nxmx, nymx), xn1a(nxmx, nymx), xnea(nxmx, nymx),
     &   xn2a(nxmx, nymx), xn3a(nxmx, nymx), bmod(nxmx, nymx),
     &   omgce(nxmx, nymx), xn4a(nxmx, nymx), xn5a(nxmx, nymx),
     &   xn6a(nxmx, nymx), omgci1(nxmx, nymx), omgci2(nxmx, nymx),
     &   omgci3(nxmx, nymx), omgci4(nxmx, nymx), omgci5(nxmx, nymx),
     &   omgci6(nxmx, nymx), omgpe2(nxmx, nymx), bpol(nxmx, nymx),
     &   capr_bpol(nxmx, nymx), omgp12(nxmx, nymx), omgp22(nxmx, nymx),
     &   omgp32(nxmx, nymx), omgp42(nxmx, nymx), omgp52(nxmx, nymx),
     &   omgp62(nxmx, nymx), xiota(nxmx, nymx), qsafety(nxmx, nymx),
     &   bmod_mid(nxmx, nymx), capr_bpol_mid2(nxmx, nymx),
     &   rho_tor2d(nxmx, nymx), psi_tor2d(nxmx, nymx), omglh(nxmx,nymx),
     &   pitch(nxmx, nymx), xnuomea(nxmx, nymx), xnuomia(nxmx, nymx))

      allocate (acold(nxmx, nymx), bcold(nxmx, nymx), ccold(nxmx, nymx),
     &   xkperp_cold(nxmx, nymx),  xkperp_cold2(nxmx, nymx),
     &   ex(nxmx, nymx), ey(nxmx, nymx), ez(nxmx, nymx),
     &   bxwave(nxmx, nymx), bywave(nxmx, nymx), bzwave(nxmx, nymx),
     &   eplus_flux_plot(nxmx, nymx),
     &   eminus_flux_plot(nxmx, nymx),
     &   xkperp_flux_plot(nxmx, nymx),
     &   ealpha(nxmx, nymx), ebeta(nxmx, nymx), eb(nxmx, nymx),
     &   eplus(nxmx, nymx), eminus(nxmx, nymx),
     &   ealphakx(nxmx, nymx), ebetakx(nxmx, nymx),
     &   ealphaky(nxmx, nymx), ebetaky(nxmx, nymx),
     &   ebkx(nxmx, nymx), ebky(nxmx, nymx),
     &   rho_ant(nxmx, nymx), rho_pla(nxmx, nymx), ntilda_e(nxmx,nymx),
     &   xb(nxmx, nymx), xc(nxmx, nymx), xd(nxmx, nymx),
     &   ntilda_e_real(nxmx, nymx)  )

      allocate (xkperp2_slow(nxmx, nymx), xkperp2_fast(nxmx, nymx) )
      allocate (xkprl_a(nxmx, nymx) )
      allocate (P_a(nxmx, nymx) )

      allocate (kappa_hatx(nxmx, nymx), kappa_haty(nxmx, nymx),
     &   kappa_hatz(nxmx, nymx), rc(nxmx, nymx),
     &   bhatx(nxmx, nymx), bhaty(nxmx, nymx), bhatz(nxmx, nymx))

      allocate (dxuxx(nxmx, nymx), dxuxy(nxmx, nymx), dxuxz(nxmx, nymx),
     &   dxuyx(nxmx, nymx), dxuyy(nxmx, nymx), dxuyz(nxmx, nymx),
     &   dxuzx(nxmx, nymx), dxuzy(nxmx, nymx), dxuzz(nxmx, nymx),
     &   dyuxx(nxmx, nymx), dyuxy(nxmx, nymx), dyuxz(nxmx, nymx),
     &   dyuyx(nxmx, nymx), dyuyy(nxmx, nymx), dyuyz(nxmx, nymx),
     &   dyuzx(nxmx, nymx), dyuzy(nxmx, nymx), dyuzz(nxmx, nymx),
     &   dyyuxx(nxmx, nymx), dyyuxy(nxmx, nymx), dyyuxz(nxmx, nymx),
     &   dyyuyx(nxmx, nymx), dyyuyy(nxmx, nymx), dyyuyz(nxmx, nymx),
     &   dyyuzx(nxmx, nymx), dyyuzy(nxmx, nymx), dyyuzz(nxmx, nymx),
     &   dxyuxx(nxmx, nymx), dxyuxy(nxmx, nymx), dxyuxz(nxmx, nymx),
     &   dxyuyx(nxmx, nymx), dxyuyy(nxmx, nymx), dxyuyz(nxmx, nymx),
     &   dxyuzx(nxmx, nymx), dxyuzy(nxmx, nymx), dxyuzz(nxmx, nymx),
     &   dxxuxx(nxmx, nymx), dxxuxy(nxmx, nymx), dxxuxz(nxmx, nymx),
     &   dxxuyx(nxmx, nymx), dxxuyy(nxmx, nymx), dxxuyz(nxmx, nymx),
     &   dxxuzx(nxmx, nymx), dxxuzy(nxmx, nymx), dxxuzz(nxmx, nymx),
     &   spx(nxmx, nymx), spy(nxmx, nymx), spz(nxmx, nymx),
     &   uxx(nxmx, nymx), uxy(nxmx, nymx), uxz(nxmx, nymx),
     &   uyx(nxmx, nymx), uyy(nxmx, nymx), uyz(nxmx, nymx),
     &   uzx(nxmx, nymx), uzy(nxmx, nymx), uzz(nxmx, nymx),
     &   bmod2(nxmx, nymx), rhom1(nxmx, nymx),
     &   bxn(nxmx, nymx), byn(nxmx, nymx), bzn(nxmx, nymx),
     &   pressi(nxmx, nymx), muhat(nxmx, nymx), nu_star(nxmx, nymx),
     &   ipsi(nxmx, nymx),
     &   xna_slo(nxmx, nymx), gradprlb(nxmx, nymx),
     &   omgp2_slo(nxmx, nymx), omgci_slo(nxmx, nymx),
     &   zeff(nxmx, nymx),
     &   rho_pol2d(nxmx, nymx), psi_pol2d(nxmx, nymx),
     &   drhodr(nxmx, nymx), drhodz(nxmx, nymx),
     &   dthedr(nxmx, nymx),  dthedz(nxmx, nymx),
     &   reomg1a(nxmx, nymx), reomg2a(nxmx, nymx), reomg3a(nxmx, nymx),
     &   xjx(nxmx, nymx), xjy(nxmx, nymx), xjz(nxmx, nymx),
     &   reomglha(nxmx, nymx))


      allocate (wdote(nxmx, nymx), wdoti1(nxmx, nymx),
     &   wdoti2(nxmx, nymx), wdoti3(nxmx, nymx), wdot(nxmx, nymx),
     &   wdoti4(nxmx, nymx), wdoti5(nxmx, nymx), wdoti6(nxmx, nymx),
     &   fype(nxmx, nymx), fypi1(nxmx, nymx),
     &   fypi2(nxmx, nymx), fypi3(nxmx, nymx), fyp(nxmx, nymx),
     &   fypi4(nxmx, nymx), fypi5(nxmx, nymx), fypi6(nxmx, nymx),
     &   fpol0e (nxmx, nymx), fpol0i1(nxmx, nymx),
     &   fpol0i2(nxmx, nymx), fpol0i3(nxmx, nymx),
     &   fpol0i4(nxmx, nymx), fpol0i5(nxmx, nymx), fpol0i6(nxmx,nymx),
     &   fx0e (nxmx, nymx), fx0i1(nxmx, nymx),
     &   fx0i2(nxmx, nymx), fx0i3(nxmx, nymx),
     &   fx0i4(nxmx, nymx), fx0i5(nxmx, nymx), fx0i6(nxmx, nymx),
     &   fy0e (nxmx, nymx), fy0i1(nxmx, nymx),
     &   fy0i2(nxmx, nymx), fy0i3(nxmx, nymx),
     &   fy0i4(nxmx, nymx), fy0i5(nxmx, nymx), fy0i6(nxmx, nymx),
     &   fz0e (nxmx, nymx), fz0i1(nxmx, nymx),
     &   fz0i2(nxmx, nymx), fz0i3(nxmx, nymx),
     &   fz0i4(nxmx, nymx), fz0i5(nxmx, nymx), fz0i6(nxmx, nymx),
     &   bdotf(nxmx, nymx),
     &   capr_fzeta(nxmx, nymx), capr2(nxmx, nymx),
     &   uzeta(nxmx, nymx), utheta(nxmx, nymx),
     &   gpsi(nxmx, nymx), omgexb(nxmx, nymx),
     &   jhat(nxmx, nymx),
     &   kpsi(nxmx, nymx), epsi(nxmx, nymx),
     &   fprl0e (nxmx, nymx), fprl0i1(nxmx, nymx),
     &   fprl0i2(nxmx, nymx), fprl0i3(nxmx, nymx),
     &   fprl0i4(nxmx, nymx), fprl0i5(nxmx, nymx), fprl0i6(nxmx, nymx),
     &   fx0(nxmx, nymx), fy0(nxmx, nymx), fz0(nxmx, nymx),
     &   fpsi0(nxmx, nymx), fpsi1(nxmx, nymx),
     &   ftheta0(nxmx, nymx), ftheta1(nxmx, nymx),
     &   fpol1e (nxmx, nymx), fpol1i1(nxmx, nymx),
     &   fpol1i2(nxmx, nymx), fpol1i3(nxmx, nymx),
     &   fpol1i4(nxmx, nymx), fpol1i5(nxmx, nymx), fpol1i6(nxmx, nymx),
     &   dfpol1e (nxmx, nymx), dfpol1i1(nxmx, nymx),
     &   dfpol1i2(nxmx, nymx), dfpol1i3(nxmx, nymx),
     &   dfpol1i4(nxmx,nymx), dfpol1i5(nxmx,nymx), dfpol1i6(nxmx,nymx),
     &   redotj(nxmx, nymx), pcre(nxmx, nymx), pcim(nxmx, nymx),
     &   divq(nxmx, nymx),
     &   fx_ant(nxmx, nymx), fy_ant(nxmx, nymx), fz_ant(nxmx, nymx),
     &   fx_pla(nxmx, nymx), fy_pla(nxmx, nymx), fz_pla(nxmx, nymx),
     &   rekxdotj(nxmx, nymx), rekydotj(nxmx, nymx),
     &   rekzdotj(nxmx, nymx),
     &   xjprl(nxmx, nymx),
     &   redotj1(nxmx, nymx), redotj2(nxmx, nymx), redotj3(nxmx, nymx),
     &   redotj4(nxmx, nymx), redotj5(nxmx, nymx), redotj6(nxmx, nymx),
     &   redotje(nxmx, nymx), redotjt(nxmx, nymx), redotji(nxmx, nymx),
     &   redotjs(nxmx, nymx),
     &   redotj_ehst(nxmx, nymx)  )

      allocate ( uxxk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &         uxyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &         uxzk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &         uyxk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &         uyyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &         uyzk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &         uzxk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &         uzyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     &         uzzk(nkdim1 : nkdim2, mkdim1 : mkdim2) )


      allocate (xjpxe(nxmx, nymx), xjpye(nxmx, nymx), xjpze(nxmx, nymx),
     &   xjpx1(nxmx, nymx), xjpy1(nxmx, nymx), xjpz1(nxmx, nymx),
     &   xjpx2(nxmx, nymx), xjpy2(nxmx, nymx), xjpz2(nxmx, nymx),
     &   xjpx3(nxmx, nymx), xjpy3(nxmx, nymx), xjpz3(nxmx, nymx),
     &   xjpx4(nxmx, nymx), xjpy4(nxmx, nymx), xjpz4(nxmx, nymx),
     &   xjpx5(nxmx, nymx), xjpy5(nxmx, nymx), xjpz5(nxmx, nymx),
     &   xjpx6(nxmx, nymx), xjpy6(nxmx, nymx), xjpz6(nxmx, nymx)  )


      allocate (xjpx_ehst(nxmx, nymx), xjpy_ehst(nxmx, nymx),
     &   xjpz_ehst(nxmx, nymx), xj_slox(nxmx, nymx),
     &   xj_sloy(nxmx, nymx), xj_sloz(nxmx, nymx),
     &   xjpx(nxmx, nymx), xjpy(nxmx, nymx), xjpz(nxmx, nymx),
     &   xjpx_lab(nxmx, nymx), xjpy_lab(nxmx, nymx),
     &   xjpz_lab(nxmx, nymx), pc(nxmx, nymx) )

      allocate (xjpxe_lab(nxmx, nymx), xjpye_lab(nxmx, nymx),
     &   xjpze_lab(nxmx, nymx) )

      allocate ( fdksav2d(nxmx, nymx) )
      allocate ( feksav2d(nxmx, nymx) )
      allocate ( ffksav2d(nxmx, nymx) )
      allocate ( fgksav2d(nxmx, nymx) )
      allocate ( faksav2d(nxmx, nymx) )
      allocate ( fpksav2d(nxmx, nymx) )
      allocate ( frksav2d(nxmx, nymx) )
      allocate ( fqksav2d(nxmx, nymx) )
      allocate ( fsksav2d(nxmx, nymx) )

      allocate ( fdksav(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( feksav(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( ffksav(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( fgksav(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( faksav(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( fpksav(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( frksav(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( fqksav(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( fsksav(nkdim1 : nkdim2, mkdim1 : mkdim2) )

      allocate(capd(nxmx, nkdim1 : nkdim2))
      allocate(capd2d(nxmx, nmodesmax))

      allocate(xkb(nxmx, nkdim1 : nkdim2))
      allocate(xkb2d(nxmx, nmodesmax))

      allocate ( ealphakmod(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( ebetakmod(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( ebkmod(nkdim1 : nkdim2, mkdim1 : mkdim2) )

      allocate ( ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2) )
      allocate ( ebk(nkdim1 : nkdim2, mkdim1 : mkdim2) )



*     --------------------------
*     setup parallel environment
*     --------------------------

*     -----------------------------------------
*     some environment may not require mpi_init
*     -----------------------------------------
      call mpi_init( info )

      open(unit=63, file='aorsa2d.in',
     &          status='old', form='formatted')
      rewind(63)


*     ------------------------
*     Read namelist input data
*     ------------------------
      read (63, aorsa2din)

      rt_init = rt
      b0_init = b0

      if (z2_electron .eq. 0) then
         use_new_z2 = .false.
         use_fourier_z2 = .false.
      end if

      if (z2_electron .eq. 1) then
         use_new_z2 = .true.
         use_fourier_z2 = .false.
      end if

      if (z2_electron .eq. 2) then
         use_new_z2 = .false.
         use_fourier_z2 = .true.
      end if


*     -------------------------------------------
*     Set the remainder of the nphi_array to zero
*     -------------------------------------------
      do i = nphi_number + 1, nphimx
         nphi_array(i) = 0.0
      end do


*     ---------------------------------------------------------
*     if (nphi_number .gt. 1), set prfin=0.0 so that pscale=1.0
*     ---------------------------------------------------------
      if(nphi_number .gt. 1) prfin=0.0

c      idiag = nmodesx / 2
c      jdiag = nmodesy / 2

      idiag = 104
      jdiag = 90

      nbessj = lmax + 2

      call blacs_pinfo( myid, nproc )
      if (nproc .lt. 1) then
         write(6,*) '** blacs_pinfo returns:myid,nproc',myid,nproc
         write(15,*) '** blacs_pinfo returns:myid,nproc',myid,nproc
         stop '** blacs not setup '
      endif



      if (myid .eq. 0) then
         open(unit=53, file='out_fields', status='unknown',
     &      form='formatted')
         rewind (53)

         open(unit=963, file='namelist_input',
     &                               status='unknown', form='formatted')

         open(unit=15, file='out15', status='unknown', form='formatted')
         open(unit=166,file='bharvey_3d',status='unknown',
     &                                                 form='formatted')
         if(nphi_sum_only .neqv. .true.)
     &      open(unit=34, file='fpm',status='unknown', form='formatted')

*        ------------------------
*        Write namelist input data
*        ------------------------
         write (963, aorsa2din)
      end if

*     -----------------------
*     write fpm (34) file
*     ------------------------
      if (myid.eq.0 .and. nphi_sum_only .neqv. .true.) then
         write(34,309) nphi1, nphi2
      end if


      call blacs_get(-1,0, icontxt)

*     --------------------
*     setup processor grid
*     --------------------
*     -------------------------------------------------------
*     1 x n processor grid would simplify pivoting
*     n x 1 processor grid would simplify eqnx,eqny,eqnz
*     sqrt(n) x sqrt(n) normally would give best performance
*     -------------------------------------------------------
      if (myid .eq. 0) then
         write(6,*) "setting up processor grid"
         write(6, *) "nphi_sum_only = ", nphi_sum_only
      end if

      call blacs_gridinit( icontxt, 'column-major', nprow,npcol)
      call blacs_gridinfo( icontxt, nprow, npcol, myrow, mycol)





      if(nphi_sum_only .eqv. .true.)go to 9001

*     ----------------
*     loop over nphi's
*     ----------------
      nt = 0

*     -------------------------------------------------------------
*     Leave all modes out of the sum except the pre-determined ones
*     -------------------------------------------------------------
      do 9000 nsum = 1, nphi_number

         fdksav = 0.0
         feksav = 0.0
         ffksav = 0.0

         fgksav = 0.0
         faksav = 0.0
         fpksav = 0.0

         frksav = 0.0
         fqksav = 0.0
         fsksav = 0.0

             fdksav2d = 0.0
         feksav2d = 0.0
         ffksav2d = 0.0

         fgksav2d = 0.0
         faksav2d = 0.0
         fpksav2d = 0.0

         frksav2d = 0.0
         fqksav2d = 0.0
         fsksav2d = 0.0

      capd = 0.0
      capd2d = 0.0

      xkb = 0.0
      xkb2d = 0.0

      nphi = nphi_array(nsum)
      nt = nt + 1

      signkz =  1.0
      if(nphi .lt. 0) signkz = -1.0

      ceta = 8.27
      xnexp = 2.48
      cfit = .0987
      afit = 12.3
      damp = 23.83/(0.678 + zeffcd)




      time0=second1(dummy)

*     ----------------------------------------------------------
*     simplifying assumption that all processors can open
*     same files and perform read
*     no need to re-broadcast input data
*     ----------------------------------------------------------
      if (myid .eq. 0) then
          write(6,*) 'blacs started: nprow,npcol,nproc ',
     &            nprow,npcol,nproc
c          write(15,*) 'blacs started: nprow,npcol,nproc ',
c     &            nprow,npcol,nproc
      endif

c      write(6, *) 'myid = ', myid


*     --------------------------------------------------
*     may need to open different files like 'out38.001'
*     for processor '001'
*     --------------------------------------------------
      suffix = '.000'
      suffix(4:4) = char(ichar('0') + mod( myid,    10))
      suffix(3:3) = char(ichar('0') + mod( myid/10, 10))
      suffix(2:2) = char(ichar('0') + mod( myid/100,10))

      open(unit=54, file='fields_fourier', status='unknown',
     &      form='formatted')
      rewind (54)


      if (myid.eq.0) then
         open(unit=38,file='out38',status='unknown',form='formatted')
         open(unit=130, file='Ono_disp', status='unknown',
     &      form='formatted')

c         open(unit=29,file='out29',status='unknown',form='formatted')
c         open(unit=30,file='out30',status='unknown',form='formatted')

         open(unit=59,file='movie_eb',status='unknown',form='formatted')
         open(unit=60,file='movie_ealpha',status='unknown',
     &      form='formatted')
         open(unit=51,file='rho',status='unknown',form='formatted')
         open(unit=69,file='movie_wdot',status='unknown',
     &        form='formatted')

         open(unit=65,file='umat',status='unknown', form='formatted')


c         open(unit=32,file='diagnostic',status='unknown',
c     &      form='formatted')


         open(unit=28, file='fields_local',   status='unknown',
     &      form='formatted')

         if (ndiste .eq. 1)  open(unit=40, file='out_cql3d.coefe',
     &                            status='unknown', form='formatted')
         if (ndisti1 .eq. 1) open(unit=41, file='out_cql3d.coef1',
     &                            status='unknown', form='formatted')
         if (ndisti2 .eq. 1) open(unit=42, file='out_cql3d.coef2',
     &                            status='unknown', form='formatted')
         if (ndisti3 .eq. 1) open(unit=43, file='out_cql3d.coef3',
     &                            status='unknown', form='formatted')
         if (ndisti4 .eq. 1) open(unit=44, file='out_cql3d.coef4',
     &                            status='unknown', form='formatted')


      endif

      nrhs = 1
      nmaxe=1


      nkx2 = nmodesx / 2
      nkx1 = - nmodesx / 2 + 1
      nnodex = nmodesx
      nnodey = nmodesy

      nnodelb = nmodeslb
      nnodelb0 = nnodelb / 2 + 1
      nkl2 = nnodelb / 2
      nkl1= nkl2 - nnodelb + 1



      nnodexm1 = nnodex - 1
      nnodeym1 = nnodey - 1

      nnoderho = nnodex / 2
      mnodetheta = nnodey / 2

      nky2 = nmodesy / 2
      nky1 = - nmodesy / 2 + 1


c      jequat  = nnodey / 2
      icenter = nnodex / 2


      if (qavg0 .ne. 0.0) xiota0 = 1./qavg0

      rholim = sqrt(psilim)

      q = 1.6e-19
      if(te0 .eq. 0.0)te0 = ti0

      t0e = te0
      t0i = ti0
      t0i2 = ti02
      t0i3 = ti03
      t0i4 = ti04
      t0i5 = ti05
      t0i6 = ti06

      teedge = 400.0

      qhat = qavg0

      eslow = eslowev * q
      t0e = t0e   * q
      t0i = t0i   * q
      t0i2 = t0i2 * q
      t0i3 = t0i3 * q
      t0i4 = t0i4 * q
      t0i5 = t0i5 * q
      t0i6 = t0i6 * q

      teedge = teedge * q

      xme = 9.11e-31
      xmh = 1.67e-27
      xmi1 = amu1 * xmh
      xmi2 = amu2 * xmh
      xmi3 = amu3 * xmh
      xmi4 = amu4 * xmh
      xmi5 = amu5 * xmh
      xmi6 = amu6 * xmh
      xmi_slo = amu_slo * xmh

      qi1 = z1 * q
      qi2 = z2 * q
      qi3 = z3 * q
      qi4 = z4 * q
      qi5 = z5 * q
      qi6 = z6 * q
      qi_slo = z_slo * q

      qe = -q
      zi = cmplx(0.0,1.0)
      eps0 = 8.85e-12

      pi = 3.141592654
      twopi = 2.0 * pi

      xlnlam = 20.0
      xmu0 = 1.26e-06
      clight = 1.0 / sqrt(eps0 * xmu0)


      xkphi0 = nphi / rt_init
      xkphia = nphi / rant
      if(xkphi0 .eq. 0.0)xkphi0 = 1.0e-05


      xnurf = freqcy

      omgrf = 2.0 * pi * xnurf


      if (plot_only .eqv. .true.)go to 9002

      vphase = omgrf / xkphi0
      vthe = sqrt(t0e / xme)
      vsound = sqrt(teedge / xmi1)
      wphase = 0.0
      if(vthe .ne. 0.0)wphase = vphase / vthe


      xkthrho = 0.2

      xkthdx = 1.0



      xk0 = omgrf / clight
      rnz = xkphi0 / xk0


      if (myid .eq. 0) then

         write(15, *) "AORSA2D: version_number = ", version_number
         write(15, *) "z2_electron = ", z2_electron
         write(15, *) "use_new_Z2 = ", use_new_Z2
         write(15, *) "use_fourier_Z2 = ", use_fourier_Z2

         write(15, *) "nphi_number = ", nphi_number
         write(15, *) "nphi_array = ", nphi_array


         write(6, *) "AORSA2D: version_number = ", version_number
         write(6, *) "z2_electron = ", z2_electron
         write(6, *) "use_new_Z2 = ", use_new_Z2
         write(6, *) "use_fourier_Z2 = ", use_fourier_Z2

         write(6, *) "nphi_number = ", nphi_number
         write(6, *) "nphi_array = ", nphi_array


         write (6, *)
         write (6, *) "xkperp_cutoff = ", xkperp_cutoff
         write (6, *) "damping       = ", damping
         write (6, *) "upshift = ", upshift
         write (6, *)
         write (6, *) "n_prof_flux = ", n_prof_flux
         write (6, *) "eqdsk       = ", eqdsk
         write (6, *) "netcdf_file1 = ", netcdf_file1
         write (6, *) "netcdf_file2 = ", netcdf_file2
         write (6, *)
         write (6, *) "nnodex  = ", nnodex
         write (6, *) "nnodey  = ", nnodey
         write (6, *) "nnodelb  = ", nnodelb
         write (6, *)


         write (6, *)
         write (6, *) "nprow   = ", nprow
         write (6, *) "npcol   = ", npcol
         write (6, *)
         write (6, *) "i_antenna = ", i_antenna
         write (6, *) "prfin   = ", prfin
         write (6, *) "freqcy  = ", freqcy
         write (6, *) "nphi    = ", nphi
         write (6, *) "antlen  = ", antlen
         write (6, *) "antlc   = ", antlc
         write (6, *) "rant    = ", rant
         write (6, *) "antang    = ", antang

         write (6, *) "yant    = ", yant
         write (6, *)
         write (6, *) "delta0  = ", delta0
         write (6, *) "xnuome  = ", xnuome
         write (6, *) "xnuomi  = ", xnuomi
         write (6, *)
         write (6, *) "lmax    = ", lmax
         write (6, *) "lmaxe    = ", lmaxe
         write (6, *) "zeffcd  = ", zeffcd
         write (6, *)
         write (6, *) "xn0     = ", xn0
         write (6, *) "xnlim   = ", xnlim
         write (6, *) "alphan  = ", alphan
         write (6, *) "betan   = ", betan
         write (6, *) "te0     = ", te0
         write (6, *) "telim   = ", telim
         write (6, *) "alphate = ", alphate
         write (6, *) "betate  = ", betate
         write (6, *)
         write (6, *) "ndisti1 = ", ndisti1
         write (6, *) "amu1    = ", amu1
         write (6, *) "z1      = ", z1
         write (6, *) "ti0     = ", ti0
         write (6, *) "tilim   = ", tilim
         write (6, *) "alphati = ", alphati
         write (6, *) "betati  = ", betati
         write (6, *)
         write (6, *) "ndisti2 = ", ndisti2
         write (6, *) "amu2    = ", amu2
         write (6, *) "z2      = ", z2
         write (6, *) "xn2     = ", xn2
         write (6, *) "xn2lim  = ", xn2lim
         write (6, *) "alphan2 = ", alphan2
         write (6, *) "betan2  = ", betan2
         write (6, *) "ti02    = ", ti02
         write (6, *) "ti2lim  = ", ti2lim
         write (6, *) "alphati2= ", alphati2
         write (6, *) "betati2 = ", betati2
         write (6, *)
         write (6, *) "ndisti3 = ", ndisti3
         write (6, *) "amu3    = ", amu3
         write (6, *) "z3      = ", z3
         write (6, *) "xn3     = ", xn3
         write (6, *) "xn3lim  = ", xn3lim
         write (6, *) "alphan3 = ", alphan3
         write (6, *) "betan3  = ", betan3
         write (6, *) "ti03    = ", ti03
         write (6, *) "ti3lim  = ", ti3lim
         write (6, *) "alphati3= ", alphati3
         write (6, *) "betati3 = ", betati3
         write (6, *)
         write (6, *) "ndisti4 = ", ndisti4
         write (6, *) "amu4    = ", amu4
         write (6, *) "z4      = ", z4
         write (6, *) "xn4     = ", xn4
         write (6, *) "xn4lim  = ", xn4lim
         write (6, *) "alphan4 = ", alphan4
         write (6, *) "betan4  = ", betan4
         write (6, *) "ti04    = ", ti04
         write (6, *) "ti4lim  = ", ti4lim
         write (6, *) "alphati4= ", alphati4
         write (6, *) "betati4 = ", betati4
         write (6, *)
         write (6, *) "ndisti5 = ", ndisti5
         write (6, *) "amu5    = ", amu5
         write (6, *) "z5      = ", z5
         write (6, *) "xn5     = ", xn5
         write (6, *) "xn5lim  = ", xn5lim
         write (6, *) "alphan5 = ", alphan5
         write (6, *) "betan5  = ", betan5
         write (6, *) "ti05    = ", ti05
         write (6, *) "ti5lim  = ", ti5lim
         write (6, *) "alphati5= ", alphati5
         write (6, *) "betati5 = ", betati5
         write (6, *)


         write (15, *)
         write (15, *) "xkperp_cutoff = ", xkperp_cutoff
         write (15, *) "damping       = ", damping
         write (15, *) "upshift = ", upshift
         write (15, *)
         write (15, *) "n_prof_flux = ", n_prof_flux
         write (15, *) "eqdsk       = ", eqdsk
         write (15, *) "netcdf_file1 = ", netcdf_file1
         write (15, *) "netcdf_file2 = ", netcdf_file2
         write (15, *)
         write (15, *) "nnodex  = ", nnodex
         write (15, *) "nnodey  = ", nnodey
         write (15, *) "nnodelb  = ", nnodelb

         write (15, *)
         write (15, *) "nprow   = ", nprow
         write (15, *) "npcol   = ", npcol
         write (15, *)
         write (15, *) "i_antenna = ", i_antenna
         write (15, *) "prfin     = ", prfin
         write (15, *) "freqcy    = ", freqcy
         write (15, *) "nphi      = ", nphi
         write (15, *) "antlen    = ", antlen
         write (15, *) "rant      = ", rant
         write (15, *) "yant      = ", yant
         write (15, *) "antang    = ", antang
         write (15, *)
         write (15, *) "delta0  = ", delta0
         write (15, *) "xnuome  = ", xnuome
         write (15, *) "xnuomi  = ", xnuomi
         write (15, *)
         write (15, *) "lmax    = ", lmax
         write (15, *) "lmaxe    = ", lmaxe

         write (15, *) "zeffcd  = ", zeffcd
         write (15, *)
         write (15, *) "xn0     = ", xn0
         write (15, *) "xnlim   = ", xnlim
         write (15, *) "alphan  = ", alphan
         write (15, *) "betan   = ", betan
         write (15, *) "te0     = ", te0
         write (15, *) "telim   = ", telim
         write (15, *) "alphate = ", alphate
         write (15, *) "betate  = ", betate
         write (15, *)
         write (15, *) "ndisti1 = ", ndisti1
         write (15, *) "amu1    = ", amu1
         write (15, *) "z1      = ", z1
         write (15, *) "ti0     = ", ti0
         write (15, *) "tilim   = ", tilim
         write (15, *) "alphati = ", alphati
         write (15, *) "betati  = ", betati
         write (15, *)
         write (15, *) "ndisti2 = ", ndisti2
         write (15, *) "amu2    = ", amu2
         write (15, *) "z2      = ", z2
         write (15, *) "xn2     = ", xn2
         write (15, *) "xn2lim  = ", xn2lim
         write (15, *) "alphan2 = ", alphan2
         write (15, *) "betan2  = ", betan2
         write (15, *) "ti02    = ", ti02
         write (15, *) "ti2lim  = ", ti2lim
         write (15, *) "alphati2= ", alphati2
         write (15, *) "betati2 = ", betati2
         write (15, *)
         write (15, *) "ndisti3 = ", ndisti3
         write (15, *) "amu3    = ", amu3
         write (15, *) "z3      = ", z3
         write (15, *) "xn3     = ", xn3
         write (15, *) "xn3lim  = ", xn3lim
         write (15, *) "alphan3 = ", alphan3
         write (15, *) "betan3  = ", betan3
         write (15, *) "ti03    = ", ti03
         write (15, *) "ti3lim  = ", ti3lim
         write (15, *) "alphati3= ", alphati3
         write (15, *) "betati3 = ", betati3
         write (15, *)
         write (15, *) "ndisti4 = ", ndisti4
         write (15, *) "amu4    = ", amu4
         write (15, *) "z4      = ", z4
         write (15, *) "xn4     = ", xn4
         write (15, *) "xn4lim  = ", xn4lim
         write (15, *) "alphan4 = ", alphan4
         write (15, *) "betan4  = ", betan4
         write (15, *) "ti04    = ", ti04
         write (15, *) "ti4lim  = ", ti4lim
         write (15, *) "alphati4= ", alphati4
         write (15, *) "betati4 = ", betati4
         write (15, *)
         write (15, *) "ndisti5 = ", ndisti5
         write (15, *) "amu5    = ", amu5
         write (15, *) "z5      = ", z5
         write (15, *) "xn5     = ", xn5
         write (15, *) "xn5lim  = ", xn5lim
         write (15, *) "alphan5 = ", alphan5
         write (15, *) "betan5  = ", betan5
         write (15, *) "ti05    = ", ti05
         write (15, *) "ti5lim  = ", ti5lim
         write (15, *) "alphati5= ", alphati5
         write (15, *) "betati5 = ", betati5
         write (15, *)



c         write(6, 7013)nmodesx, nmodesy
c         write(6, 7113)nwdot


c        write(6, *)"nprow = ", nprow
c        write(6, *)"npcol = ", npcol
c        write(6, *)"n_bin = ", n_bin

c         write(6, *) "i_write = ", i_write
c        write(15, *)"i_write = ", i_write


         write(6, 1815)te0
         write(6, 1821)ti0

         write(6, 7217)qavg0
         write(6, 7014)lmax
         write(6,17014)lmaxe

         write(6, 7015)ibessel
         write(6, 7115)nzfun
         write(6, 71160)xnuome
         write(6, 71161)xnuomi
c         write(6, 7016)rhoi1overl
c         write(6, 7017)rhoi10
         write(6, 7217)qavg0
         write(6, 1013)nphi
         write(6, 1321)wphase
         write(6, 1323)vthe
         write(6, 1021)rnz
         write(6, 1882)rt_init
         write(6, 1822)aplasm
         write(6, 1823)rant
         write(6, 1889)b0_init
         write(6, 6813)xn0
         write(6, 1813)xn1
         write(6, 1814)xn2
         write(6, 1834)xn3

         write(6, 1012)omgrf
         write(6, 1009)xnurf
         write(6, 1013)nphi
         write(6, 1321)wphase
         write(6, 1322)vphase
         write(6, 1323)vthe
         write(6, 1714)xk0
         write(6, 1016)xnuead
         write(6, 1017)xnu1ad
         write(6, 1018)xnu2ad
         write(6, 1020)nnodex, nnodey

         write(6, *) "iprofile = ", iprofile
         write(6, *) "psilim = ", psilim


         write(15, 162)
c         write(15, 7013)nmodesx, nmodesy
         write(15, 7113)nwdot


c         write(15,*)"nprow = ", nprow
c         write(15,*)"npcol = ", npcol

         write(15, 1815)te0
         write(15, 1821)ti0

         write(15, 7217)qavg0
         write(15, 7014)lmax
         write(15,17014)lmaxe

         write(15, 7015)ibessel
         write(15, 7115)nzfun
         write(15, 71160)xnuome
         write(15, 71161)xnuomi
c         write(15, 7016)rhoi1overl
c         write(15, 7017)rhoi10
         write(15, 7217)qavg0
         write(15, 1013)nphi
         write(15, 1321)wphase
         write(15, 1323)vthe
         write(15, 1021)rnz
         write(15, 1882)rt_init
         write(15, 1822)aplasm
         write(15, 1823)rant
         write(15, 1889)b0_init
         write(15, 6813)xn0
         write(15, 1813)xn1
         write(15, 1814)xn2
         write(15, 1834)xn3

         write(15, 1012)omgrf
         write(15, 1009)xnurf
         write(15, 1013)nphi
         write(15, 1321)wphase
         write(15, 1322)vphase
         write(15, 1323)vthe
         write(15, 1714)xk0
         write(15, 1016)xnuead
         write(15, 1017)xnu1ad
         write(15, 1018)xnu2ad
         write(15, 1020)nnodex, nnodey

         write(15, *) "iprofile = ", iprofile
         write(15, *) "psilim = ", psilim
      endif

      telimj   = telim  * q
      tilimj   = tilim  * q
      ti2limj  = ti2lim * q
      ti3limj  = ti3lim * q
      ti4limj  = ti4lim * q
      ti5limj  = ti5lim * q
      ti6limj  = ti6lim * q


*------------------------------------
*     Store factorials from 0 to lmax
*------------------------------------
      do l = 0, lmax
         factl_sav(l) = factrl(l)
      end do


      t1 = second1(dummy)

*------------------------------------------
*     Flux surfaces from EQDSK (ga version)
*------------------------------------------
      call blacs_barrier(icontxt, 'All')


      if(igeom .eq. 5) then

         if(eqdsk_special .eq. 0)
     &      call eqdsk_setup(myid, eqdsk, nmodesx, nmodesy,
     &      rwleft, rwright, ytop, ybottom,
     &      rmaxis, zmaxis, b0, psio, psimag, psi_tor_max,
     &      bxn, byn, bzn, bmod, psi_pol2d, rho_pol2d, qsafety,
     &      bmod_mid,  capr_bpol_mid2,  capr_bpol_mid, rho_tor2d,
     &      i_psi_eq, dldb_tot12, dldbavg, n_prof_flux, rhomax)

         if(eqdsk_special .ne. 0)
     &      call eqdsk_setup2(myid, eqdsk, nmodesx, nmodesy,
     &      rwleft, rwright, ytop, ybottom,
     &      rmaxis, zmaxis, b0, psio, psimag, psi_tor_max,
     &      bxn, byn, bzn, bmod, psi_pol2d, rho_pol2d, qsafety,
     &      bmod_mid,  capr_bpol_mid2,  capr_bpol_mid, rho_tor2d,
     &      i_psi_eq, dldb_tot12, dldbavg, n_prof_flux, rhomax)

         call blacs_barrier(icontxt, 'All')


 9318    format(a128)
         rt = rmaxis

         r0 = rmaxis
         z0 = zmaxis

         if (myid .eq. 0)then
            write (6, *) "rwleft  = ", rwleft
            write (6, *) "rwright = ", rwright
            write (6, *) "ybottom = ", ybottom
            write (6, *) "ytop    = ", ytop

            write (15, *)
            write (15, *) "rwleft  = ", rwleft
            write (15, *) "rwright = ", rwright
            write (15, *) "ybottom = ", ybottom
            write (15, *) "ytop    = ", ytop

            write (6, *) "eqdsk       = ", eqdsk
            write (15, *) "eqdsk       = ", eqdsk
         end if






*        ----------------------------------------------------
*        Default:  if n_prof_flux equals 0, use poloidal flux
*        ----------------------------------------------------

            do i = 1, nnodex
               do j = 1, nnodey
                  rho(i,j) = rho_pol2d(i,j)
                  psi(i,j) = psi_pol2d(i,j)
                  psi_dim(i,j) = psi(i,j) * psio
               end do
            end do


*        --------------------------------------------------------------
*        Alternate: if n_prof_flux does not equal 0,  use toroidal flux
*        Note:  Here we use toroidal flux  inside rho_pol = 1.0
*                       and poloidal flux outside rho_pol = 1.0!!
*        --------------------------------------------------------------
         if(n_prof_flux .ne. 0)then
            do i = 1, nnodex
               do j = 1, nnodey
                  psi_tor2d(i,j) = rho_tor2d(i,j)**2

                  if(rho(i,j) .lt. 1.0)then
                     rho(i,j) = rho_tor2d(i,j)
                     psi(i,j) = psi_tor2d(i,j)
                     psi_dim(i,j) = psi(i,j) * psi_tor_max
                  end if

               end do
            end do
         end if

*        -------------------------------------------------------
*        Recalculate rnz and vphase with major radius from eqdsk
*        --------------------------------------------- ---------
         xkphi0 = nphi / rt
         if(xkphi0 .eq. 0.0)xkphi0 = 1.0e-05

         vphase = omgrf / xkphi0

         wphase = 0.0
         if(vthe .ne. 0.0)wphase = vphase / vthe
         rnz = xkphia / xk0


         if(myid .eq. 0)then
            if(n_prof_flux .eq. 0)write(6, *)'use sqrt(poloidal flux)'
            if(n_prof_flux .ne. 0)write(6, *)'use sqrt(toroidal flux)'

            write(15, 1321)wphase
            write(15, 1323)vthe
c            write(15, 1022)rnz

            write(6, 1321)wphase
            write(6, 1323)vthe
c            write(6, 1022)rnz

            write(6, 1812)rt
            write(6, 1809)b0
            write(6, *) "rmaxis = ", rmaxis
            write(6, *) "zmaxis = ", zmaxis

            if(n_prof_flux .eq. 0)write(15,*)'use sqrt(poloidal flux)'
            if(n_prof_flux .ne. 0)write(15,*)'use sqrt(toroidal flux)'
            write(15, 1812)rt
            write(15, 1809)b0
         end if


         omgci10 = qi1 * b0 / xmi1
         vthi10 = sqrt(2.0 * t0i / xmi1)
         rhoi10 = vthi10 / omgci10
         if(myid .eq. 0)write(15, *)"rhoi10 = ", rhoi10

         omgci20 = qi2 * b0 / xmi2
         vthi20 = sqrt(2.0 * t0i2 / xmi2)
         rhoi20 = vthi20 / omgci20
         if(myid .eq. 0)write(15, *)"rhoi20 = ", rhoi20

         omgci30 = qi3 * b0 / xmi3
         vthi30 = sqrt(2.0 * t0i3 / xmi3)
         rhoi30 = vthi30 / omgci30
         if(myid .eq. 0)write(15, *)"rhoi30 = ", rhoi30

         omgci40 = qi4 * b0 / xmi4
         vthi40 = sqrt(2.0 * t0i4 / xmi4)
         rhoi40 = vthi40 / omgci40
         if(myid .eq. 0)write(15, *)"rhoi40 = ", rhoi40

         omgci50 = qi5 * b0 / xmi5
         vthi50 = sqrt(2.0 * t0i5 / xmi5)
         rhoi50 = vthi50 / omgci50
         if(myid .eq. 0)write(15, *)"rhoi50 = ", rhoi50


         do i = 1, nnodex
            do j = 1, nnodey
               btau(i,j) = sqrt(bxn(i,j)**2 + byn(i,j)**2)
               bpol(i,j) = btau(i,j) * bmod(i,j)
               bzeta(i,j) = bzn(i,j)
               pitch(i,j) = btau(i,j) / bzeta(i,j)


*              ---------------------------
*              Calculate rotation matrix U
*              ---------------------------

               sqx = sqrt(1.0 - bxn(i,j)**2)

               uxx(i, j) =   sqx
               uxy(i, j) = - bxn(i, j) * byn(i, j) / sqx
               uxz(i, j) = - bxn(i, j) * bzn(i, j) / sqx
               uyx(i, j) =   0.0
               uyy(i, j) =   bzn(i, j) / sqx
               uyz(i, j) = - byn(i, j) / sqx
               uzx(i, j) =   bxn(i, j)
               uzy(i, j) =   byn(i, j)
               uzz(i, j) =   bzn(i, j)

            end do
         end do

      end if
*     -------------------
*     End of igeom .eq. 5
*     -------------------

c----------------------------------------------------
c     Solovev flux surfaces with analytic derivatives:
c----------------------------------------------------
      if(igeom .eq. 2)then

         psi_lim = xiota0 * b0 / 2.0 *
     &       ( (rlim**2 - rt**2)**2 / 4. / rt**2 )


         do i = 1, nnodex
            do j = 1, nnodey
               psi1 = b0 * xiota0 / 2. *
     1             ((capr(i) * y(j) / rt / ekappa)**2
     1             + (capr(i)**2 - rt**2)**2 / 4. / rt**2 )

               psi(i,j) = psi1 / psi_lim
               if(psi(i,j) .eq. 0.0)psi(i,j) = 1.0e-10


               rho(i,j) = sqrt(psi(i,j))


               dxpsi = b0 * xiota0 / (2. * psi_lim) *
     &             (2. * capr(i) * (y(j) / rt / ekappa)**2
     &             + (capr(i)**2 - rt**2) * capr(i) / rt**2 )

               dxxpsi = b0 * xiota0 / (2. * psi_lim) *
     &             (2. * (y(j) / rt / ekappa)**2
     &             + 3.0 * capr(i)**2 / rt**2 - 1.0 )

               dxypsi = b0 * xiota0 / psi_lim *
     &             2.0 * capr(i) * y(j) / (rt * ekappa)**2

               dypsi = b0 * xiota0 / psi_lim *
     &             capr(i)**2 * y(j) / (rt * ekappa)**2

               dyypsi = b0 * xiota0 / psi_lim *
     &                       (capr(i) / rt / ekappa)**2




               denom = ekappa / 2. / capr(i) *(capr(i)**2 - rt**2)
               if(y(j) .ne. 0.0 .or. denom .ne. 0.0) then
                  theta0(i,j) = atan2(y(j), denom)
                  if(theta0(i,j) .ge. 0.0) theta(i,j) = theta0(i,j)
                  if(theta0(i,j) .lt. 0.0) theta(i,j) =
     1                                       theta0(i,j) + 2.0 * pi
               end if




               bx(i,j) = -b0 * xiota0 * capr(i) * y(j)
     1             / rt**2 / ekappa**2
               by(i,j) = xiota0 * b0 / 2.0 *
     1                ( 2.0 * y(j)**2 / rt**2 / ekappa**2
     1                 + capr(i)**2 / rt**2 - 1.0)



               gaussian =  exp(-psi(i,j) / psipne)

c--            use for regular runs (default):
               if (iqprof .eq. 1)then
                  frho = gaussian
                  dxfrho = - frho / psipne * dxpsi
                  dyfrho = - frho / psipne * dypsi
                  dxxfrho = dxfrho**2 / frho - frho / psipne * dxxpsi
                  dyyfrho = dyfrho**2 / frho - frho / psipne * dyypsi
                  dxyfrho = dxfrho * dyfrho / frho
     &                                        - frho / psipne * dxypsi
               end if

c--            use for TAE modes:
               if (iqprof .eq. 2)then
                  frho = gaussian**(0.5)
               end if

               bx(i,j) = bx(i,j) * frho
               by(i,j) = by(i,j) * frho
               bz(i,j) = b0 * rt / capr(i)

               dxbx = bx(i,j) * (1. / capr(i) + dxfrho / frho)

               dybx = - b0 * xiota0 * capr(i)
     &             / rt**2 / ekappa**2 * frho + bx(i,j) / frho * dyfrho

               dxxbx = dxbx**2 / bx(i,j) + bx(i,j) * (-1.0 / capr(i)**2
     &            - (dxfrho / frho)**2 + dxxfrho / frho)

               dyybx = - b0 * xiota0 * capr(i)
     &             / rt**2 / ekappa**2 * dyfrho + dybx * dyfrho / frho
     &             + bx(i,j) * (- dyfrho**2 / frho**2 + dyyfrho / frho)

               dxybx = dxbx * dybx / bx(i,j) + bx(i,j) * (
     &            - dxfrho * dyfrho / frho**2 + dxyfrho / frho)



               dxby = xiota0 * b0 * capr(i) / rt**2 * frho
     &                                       + by(i,j) / frho * dxfrho
               dyby = xiota0 * b0 * 2.0 * y(j) / rt**2
     &                    / ekappa**2 * frho + by(i,j) / frho * dyfrho

               dxxby =  xiota0 * b0 / rt**2 * (frho + capr(i) * dxfrho)
     &            + dxby * dxfrho / frho  + by(i,j) *
     &                        ( - (dxfrho / frho)**2  + dxxfrho / frho)
               dyyby =  xiota0 * b0 * 2.0 / rt**2 / ekappa**2
     &                * (frho + y(j) * dyfrho) + dyby * dyfrho / frho
     &            + by(i,j) * (- (dyfrho / frho)**2  + dyyfrho / frho)
               dxyby =  xiota0 * b0 * capr(i) / rt**2 * dyfrho +
     &             dyby * dxfrho / frho + by(i,j) *
     &                (- dxfrho * dyfrho / frho**2  + dxyfrho / frho)


               dxbz = - bz(i,j) / capr(i)
               dxxbz = - dxbz / capr(i) + bz(i,j) / capr(i)**2
               dybz = 0.0
               dyybz = 0.0
               dxybz = 0.0



               bmod(i, j) = sqrt(bx(i, j)**2 + by(i, j)**2
     1                      + bz(i, j)**2)


               dxdth(i, j) = - y(j) / ekappa
               dzdth(i, j) = ekappa * (capr(i)**2 - rt**2) /
     1            (2. *capr(i)) + y(j)**2 / (capr(i) * ekappa)
               xntau(i, j) = sqrt(dxdth(i, j)**2 + dzdth(i, j)**2)
               if(xntau(i,j) .eq. 0.0)xntau(i,j) = 1.0e-08


               btau(i,j) = 1.0 / xntau(i, j) *
     1               (dxdth(i, j) * bx(i, j) + dzdth(i, j) * by(i, j))
     1               / bmod(i,j)
               bpol(i,j) = btau(i,j) * bmod(i,j)
               bzeta(i,j) = bz(i,j) / bmod(i,j)

               bxn(i,j) = bx(i,j) / bmod(i,j)
               byn(i,j) = by(i,j) / bmod(i,j)
               bzn(i,j) = bz(i,j) / bmod(i,j)

               dxbmod = bxn(i,j) * dxbx
     &                + byn(i,j) * dxby
     &                + bzn(i,j) * dxbz
               dybmod = bxn(i,j) * dybx
     &                + byn(i,j) * dyby
     &                + bzn(i,j) * dybz

               dxbxn = bxn(i,j) * (dxbx / bx(i,j) - dxbmod / bmod(i,j))
               dybxn = bxn(i,j) * (dybx / bx(i,j) - dybmod / bmod(i,j))

               dxbyn = byn(i,j) * (dxby / by(i,j) - dxbmod / bmod(i,j))
               dybyn = byn(i,j) * (dyby / by(i,j) - dybmod / bmod(i,j))

               dxbzn = bzn(i,j) * (dxbz / bz(i,j) - dxbmod / bmod(i,j))
               dybzn = bzn(i,j) * (dybz / bz(i,j) - dybmod / bmod(i,j))



               dxxbmod = dxbxn * dxbx
     &                 + dxbyn * dxby
     &                 + dxbzn * dxbz
     &                 + bxn(i,j) * dxxbx
     &                 + byn(i,j) * dxxby
     &                 + bzn(i,j) * dxxbz
               dyybmod = dybxn * dybx
     &                 + dybyn * dyby
     &                 + dybzn * dybz
     &                 + bxn(i,j) * dyybx
     &                 + byn(i,j) * dyyby
     &                 + bzn(i,j) * dyybz
               dxybmod = dybxn * dxbx
     &                 + dybyn * dxby
     &                 + dybzn * dxbz
     &                 + bxn(i,j) * dxybx
     &                 + byn(i,j) * dxyby
     &                 + bzn(i,j) * dxybz




               dxxbxn = dxbxn**2 / bxn(i,j)
     &             + bxn(i,j) * (-dxbx**2 / bx(i,j)**2 + dxxbx / bx(i,j)
     &             + (dxbmod / bmod(i,j))**2  - dxxbmod / bmod(i,j))
               dyybxn = dybxn**2 / bxn(i,j)
     &             + bxn(i,j) * (-dybx**2 / bx(i,j)**2 + dyybx / bx(i,j)
     &             + (dybmod / bmod(i,j))**2  - dyybmod / bmod(i,j))
               dxybxn = dxbxn * dybxn / bxn(i,j)
     &             + bxn(i,j) * (- dxbx * dybx / bx(i,j)**2
     &             + dxybx / bx(i,j)
     &             + dxbmod * dybmod / bmod(i,j)**2
     &             - dxybmod / bmod(i,j))


               dxxbyn = dxbyn**2 / byn(i,j)
     &             + byn(i,j) * (-dxby**2 / by(i,j)**2 + dxxby / by(i,j)
     &             + (dxbmod / bmod(i,j))**2  - dxxbmod / bmod(i,j))
               dyybyn = dybyn**2 / byn(i,j)
     &             + byn(i,j) * (-dyby**2 / by(i,j)**2 + dyyby / by(i,j)
     &             + (dybmod / bmod(i,j))**2  - dyybmod / bmod(i,j))
               dxybyn = dxbyn * dybyn / byn(i,j)
     &             + byn(i,j) * (- dxby * dyby / by(i,j)**2
     &             + dxyby / by(i,j)
     &             + dxbmod * dybmod / bmod(i,j)**2
     &             - dxybmod / bmod(i,j))


               dxxbzn = dxbzn**2 / bzn(i,j)
     &             + bzn(i,j) * (-dxbz**2 / bz(i,j)**2 + dxxbz / bz(i,j)
     &             + (dxbmod / bmod(i,j))**2  - dxxbmod / bmod(i,j))
               dyybzn = dybzn**2 / bzn(i,j)
     &             + bzn(i,j) * (-dybz**2 / bz(i,j)**2 + dyybz / bz(i,j)
     &             + (dybmod / bmod(i,j))**2  - dyybmod / bmod(i,j))
               dxybzn = dxbzn * dybzn / bzn(i,j)
     &             + bzn(i,j) * (- dxbz * dybz / bz(i,j)**2
     &             + dxybz / bz(i,j)
     &             + dxbmod * dybmod / bmod(i,j)**2
     &             - dxybmod / bmod(i,j))


               xiota(i,j) = btau(i,j)/ bzeta(i,j) * capr(i) / xntau(i,j)
               if(xiota(i,j) .eq. 0.0) xiota(i,j) = 1.0e-06
               qsafety(i,j) = 1.0 / xiota(i,j)

c               if(j .eq. jequat)write(6, 1312) i, capr(i), btau(i,j),
c     1             bzeta(i,j)

               sqx = sqrt(1.0 - bxn(i,j)**2)

               dxsqx = - bxn(i,j) / sqx * dxbxn
               dysqx = - bxn(i,j) / sqx * dybxn
               dxxsqx =  dxsqx * (dxbxn / bxn(i,j) - dxsqx / sqx)
     &                  - bxn(i,j) / sqx * dxxbxn
               dyysqx =  dysqx * (dybxn / bxn(i,j) - dysqx / sqx)
     &                  - bxn(i,j) / sqx * dyybxn
               dxysqx =  dxsqx * (dybxn / bxn(i,j) - dysqx / sqx)
     &                  - bxn(i,j) / sqx * dxybxn


               uxx(i, j) =   sqx
               uxy(i, j) = - bxn(i, j) * byn(i, j) / sqx
               uxz(i, j) = - bxn(i, j) * bzn(i, j) / sqx
               uyx(i, j) =   0.0
               uyy(i, j) =   bzn(i, j) / sqx
               uyz(i, j) = - byn(i, j) / sqx
               uzx(i, j) =   bxn(i, j)
               uzy(i, j) =   byn(i, j)
               uzz(i, j) =   bzn(i, j)

               dxuxx(i,j)  = dxsqx
               dxxuxx(i,j) = dxxsqx
               dyuxx(i,j)  = dysqx
               dyyuxx(i,j) = dyysqx
               dxyuxx(i,j) = dxysqx


               dxuxy(i,j)  = uxy(i,j) * (dxbxn / bxn(i,j)
     &                                 + dxbyn / byn(i,j) - dxsqx / sqx)
               dyuxy(i,j)  = uxy(i,j) * (dybxn / bxn(i,j)
     &                                 + dybyn / byn(i,j) - dysqx / sqx)
               dxxuxy(i,j) = dxuxy(i,j)**2 / uxy(i,j)
     &               + uxy(i,j) * (dxxbxn / bxn(i,j)
     &                           + dxxbyn / byn(i,j)
     &                           - dxxsqx / sqx  )
     &               - uxy(i,j) * (dxbxn**2 / bxn(i,j)**2
     &                           + dxbyn**2 / byn(i,j)**2
     &                           - dxsqx**2 / sqx**2 )
               dyyuxy(i,j) = dyuxy(i,j)**2 / uxy(i,j)
     &               + uxy(i,j) * (dyybxn / bxn(i,j)
     &                           + dyybyn / byn(i,j)
     &                           - dyysqx / sqx  )
     &               - uxy(i,j) * (dybxn**2 / bxn(i,j)**2
     &                           + dybyn**2 / byn(i,j)**2
     &                           - dysqx**2 / sqx**2 )
               dxyuxy(i,j) = dxuxy(i,j) * dyuxy(i,j) / uxy(i,j)
     &               + uxy(i,j) * (dxybxn / bxn(i,j)
     &                           + dxybyn / byn(i,j)
     &                           - dxysqx / sqx  )
     &               - uxy(i,j) * (dxbxn * dybxn / bxn(i,j)**2
     &                          +  dxbyn * dybyn / byn(i,j)**2
     &                          -  dxsqx * dysqx / sqx**2   )


               dxuxz(i,j)  = uxz(i,j) * (dxbxn / bxn(i,j)
     &                                 + dxbzn / bzn(i,j) - dxsqx / sqx)
               dyuxz(i,j)  = uxz(i,j) * (dybxn / bxn(i,j)
     &                                 + dybzn / bzn(i,j) - dysqx / sqx)
               dxxuxz(i,j) = dxuxz(i,j)**2 / uxz(i,j)
     &               + uxz(i,j) * (dxxbxn / bxn(i,j)
     &                           + dxxbzn / bzn(i,j)
     &                           - dxxsqx / sqx  )
     &               - uxz(i,j) * (dxbxn**2 / bxn(i,j)**2
     &                           + dxbzn**2 / bzn(i,j)**2
     &                           - dxsqx**2 / sqx**2 )
               dyyuxz(i,j) = dyuxz(i,j)**2 / uxz(i,j)
     &               + uxz(i,j) * (dyybxn / bxn(i,j)
     &                           + dyybzn / bzn(i,j)
     &                           - dyysqx / sqx  )
     &               - uxz(i,j) * (dybxn**2 / bxn(i,j)**2
     &                           + dybzn**2 / bzn(i,j)**2
     &                           - dysqx**2 / sqx**2 )
               dxyuxz(i,j) = dxuxz(i,j) * dyuxz(i,j) / uxz(i,j)
     &               + uxz(i,j) * (dxybxn / bxn(i,j)
     &                           + dxybzn / bzn(i,j)
     &                           - dxysqx / sqx  )
     &               - uxz(i,j) * (dxbxn * dybxn / bxn(i,j)**2
     &                          +  dxbzn * dybzn / bzn(i,j)**2
     &                          -  dxsqx * dysqx / sqx**2   )


               dxuyx(i,j)  = 0.0
               dxxuyx(i,j) = 0.0
               dyuyx(i,j)  = 0.0
               dyyuyx(i,j) = 0.0
               dxyuyx(i,j) = 0.0

               dxuyy(i,j) = uyy(i,j) * (dxbzn / bzn(i,j) - dxsqx / sqx)
               dyuyy(i,j) = uyy(i,j) * (dybzn / bzn(i,j) - dysqx / sqx)
               dxxuyy(i,j) = dxuyy(i,j)**2 / uyy(i,j)
     &            + uyy(i,j) * (dxxbzn / bzn(i,j) - dxxsqx / sqx)
     &            - uyy(i,j) * (dxbzn**2 / bzn(i,j)**2
     &                        - dxsqx**2 / sqx**2        )
               dyyuyy(i,j) = dyuyy(i,j)**2 / uyy(i,j)
     &            + uyy(i,j) * (dyybzn / bzn(i,j) - dyysqx / sqx)
     &            - uyy(i,j) * (dybzn**2 / bzn(i,j)**2
     &                        - dysqx**2 / sqx**2        )
               dxyuyy(i,j) = dyuyy(i,j) * dxuyy(i,j) / uyy(i,j)
     &            + uyy(i,j) * (dxybzn / bzn(i,j) - dxysqx / sqx)
     &            - uyy(i,j) * (dxbzn * dybzn / bzn(i,j)**2
     &                        - dxsqx * dysqx / sqx**2        )


               dxuyz(i,j) = uyz(i,j) * (dxbyn / byn(i,j) - dxsqx / sqx)
               dyuyz(i,j) = uyz(i,j) * (dybyn / byn(i,j) - dysqx / sqx)
               dxxuyz(i,j) = dxuyz(i,j)**2 / uyz(i,j)
     &            + uyz(i,j) * (dxxbyn / byn(i,j) - dxxsqx / sqx)
     &            - uyz(i,j) * (dxbyn**2 / byn(i,j)**2
     &                        - dxsqx**2 / sqx**2        )
               dyyuyz(i,j) = dyuyz(i,j)**2 / uyz(i,j)
     &            + uyz(i,j) * (dyybyn / byn(i,j) - dyysqx / sqx)
     &            - uyz(i,j) * (dybyn**2 / byn(i,j)**2
     &                        - dysqx**2 / sqx**2        )
               dxyuyz(i,j) = dxuyz(i,j) * dyuyz(i,j) / uyz(i,j)
     &            + uyz(i,j) * (dxybyn / byn(i,j) - dxysqx / sqx)
     &            - uyz(i,j) * (dxbyn * dybyn / byn(i,j)**2
     &                        - dxsqx * dysqx / sqx**2        )


               dxuzx(i,j)  = dxbxn
               dxxuzx(i,j) = dxxbxn
               dyuzx(i,j)  = dybxn
               dyyuzx(i,j) = dyybxn
               dxyuzx(i,j) = dxybxn


               dxuzy(i,j)  = dxbyn
               dxxuzy(i,j) = dxxbyn
               dyuzy(i,j)  = dybyn
               dyyuzy(i,j) = dyybyn
               dxyuzy(i,j) = dxybyn


               dxuzz(i,j)  = dxbzn
               dxxuzz(i,j) = dxxbzn
               dyuzz(i,j)  = dybzn
               dyyuzz(i,j) = dyybzn
               dxyuzz(i,j) = dxybzn

*              -----
*              Exact:
*              -----
c               gradprlb(i,j) = uzx(i,j) * dxbmod + uzy(i,j) * dybmod

*              -----------------------
*              Brambilla approximation:
*              -----------------------
               sinth = y(j) / sqrt(x(i)**2 + y(j)**2)
               gradprlb(i,j) = bmod(i,j) / capr(i)
     &                                       * abs(btau(i,j) * sinth)

*              -------------------------------
*              Constant gradient scale length:
*              -------------------------------
c               gradprlb(i,j) = bmod(i,j) / capr(i)



               if (nzfun .eq. 0)gradprlb(i,j) = 1.0e-10


               drhodx = 0.5 / sqrt(psi(i,j)) * dxpsi
               drhody = 0.5 / sqrt(psi(i,j)) * dypsi

               gradrho = sqrt(drhodx**2 + drhody**2)
               if (gradrho .eq. 0.0) gradrho = 1.0e-08

               rhohatx(i,j) = drhodx / gradrho
               rhohaty(i,j) = drhody / gradrho

            end do
         end do


      end if
*     --------------------------------------------
*     End of igeom .eq. 2 and analytic derivatives
*     --------------------------------------------


      xwleft = rwleft - rt
      xwright = rwright - rt

      akprl_min = 0.25 * abs(nphi / rt)
c      akprl_min = abs(nphi / rt)





*     --------------------------------
*     Distribution function from CQL3D
*     --------------------------------

!     -----------------
!     allocate arrays
!     -----------------

      if(ndisti1   .eq. 0 .and.
     &   ndisti2   .eq. 0 .and.
     &   ndisti3   .eq. 0 .and.
     &   ndisti4   .eq. 0 .and.
     &   ndisti5   .eq. 0 .and.
     &   ndisti6   .eq. 0 .and.
     &   ndiste    .eq. 0)   then

         allocate( UPERP(nuper) )
         allocate( UPARA(nupar) )
         allocate( VPERP(nuper) )
         allocate( VPARA(nupar) )

         allocate( UPERP_work(nuper) )
         allocate( UPARA_work(nupar) )

         allocate( f(nuper, nupar) )

         allocate( dfdupere(nuper, nupar) )
         allocate( dfdupare(nuper, nupar) )
         allocate( fperpe(nuper) )

         allocate( dfduper1(nuper, nupar) )
         allocate( dfdupar1(nuper, nupar) )
         allocate( fperp1(nuper) )

         allocate( dfduper2(nuper, nupar) )
         allocate( dfdupar2(nuper, nupar) )
         allocate( fperp2(nuper) )

         allocate( dfduper3(nuper, nupar) )
         allocate( dfdupar3(nuper, nupar) )
         allocate( fperp3(nuper) )

         allocate( dfduper4(nuper, nupar) )
         allocate( dfdupar4(nuper, nupar) )
         allocate( fperp4(nuper) )

         allocate( dfduper5(nuper, nupar) )
         allocate( dfdupar5(nuper, nupar) )
         allocate( fperp5(nuper) )

         allocate( dfduper6(nuper, nupar) )
         allocate( dfdupar6(nuper, nupar) )
         allocate( fperp6(nuper) )




         UPERP = 0.0
         UPARA = 0.0
         UPERP_work = 0.0
         UPARA_work = 0.0

         f = 0.0

         dfdupere = 0.0
         dfdupare = 0.0

         dfduper1 = 0.0
         dfdupar1 = 0.0

         dfduper2 = 0.0
         dfdupar2 = 0.0

         dfduper3 = 0.0
         dfdupar3 = 0.0

         dfduper4 = 0.0
         dfdupar4 = 0.0

         dfduper5 = 0.0
         dfdupar5 = 0.0

         dfduper6 = 0.0
         dfdupar6 = 0.0

      end if

      if(ndisti1   .ne. 0 .or.
     &   ndisti2   .ne. 0 .or.
     &   ndisti3   .ne. 0 .or.
     &   ndisti4   .ne. 0 .or.
     &   ndisti5   .ne. 0 .or.
     &   ndisti6   .ne. 0 .or.
     &   ndiste    .ne. 0)   then


!        -----------------
!        allocate arrays
!        -----------------
         allocate( UPERP(nuper) )
         allocate( UPARA(nupar) )
         allocate( VPERP(nuper) )
         allocate( VPARA(nupar) )
         allocate( UPERP_work(nuper) )
         allocate( UPARA_work(nupar) )


         allocate( f(nuper, nupar) )

         allocate( dfdupere(nuper, nupar) )
         allocate( dfdupare(nuper, nupar) )
         allocate( fperpe(nuper) )

         allocate( dfduper1(nuper, nupar) )
         allocate( dfdupar1(nuper, nupar) )
         allocate( fperp1(nuper) )

         allocate( dfduper2(nuper, nupar) )
         allocate( dfdupar2(nuper, nupar) )
         allocate( fperp2(nuper) )

         allocate( dfduper3(nuper, nupar) )
         allocate( dfdupar3(nuper, nupar) )
         allocate( fperp3(nuper) )

         allocate( dfduper4(nuper, nupar) )
         allocate( dfdupar4(nuper, nupar) )
         allocate( fperp4(nuper) )

         allocate( dfduper5(nuper, nupar) )
         allocate( dfdupar5(nuper, nupar) )
         allocate( fperp5(nuper) )

         allocate( dfduper6(nuper, nupar) )
         allocate( dfdupar6(nuper, nupar) )
         allocate( fperp6(nuper) )


         UPERP = 0.0
         UPARA = 0.0
         UPERP_work = 0.0
         UPARA_work = 0.0

         f = 0.0

         dfdupere = 0.0
         dfdupare = 0.0

         dfduper1 = 0.0
         dfdupar1 = 0.0

         dfduper2 = 0.0
         dfdupar2 = 0.0

         dfduper3 = 0.0
         dfdupar3 = 0.0

         dfduper4 = 0.0
         dfdupar4 = 0.0

         dfduper5 = 0.0
         dfdupar5 = 0.0

         dfduper6 = 0.0
         dfdupar6 = 0.0

         allocate( fe_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( dfe_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( dfe_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f1_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df1_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df1_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f2_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df2_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df2_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f3_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df3_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df3_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f4_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df4_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df4_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f5_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df5_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df5_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f6_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df6_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df6_cql_uprl(nuper, nupar, n_psi_dim) )



         fe_cql_cart = 0.0
         dfe_cql_uprp = 0.0
         dfe_cql_uprl = 0.0

         f1_cql_cart = 0.0
         df1_cql_uprp = 0.0
         df1_cql_uprl = 0.0

         f2_cql_cart = 0.0
         df2_cql_uprp = 0.0
         df2_cql_uprl = 0.0

         f3_cql_cart = 0.0
         df3_cql_uprp = 0.0
         df3_cql_uprl = 0.0

         f4_cql_cart = 0.0
         df4_cql_uprp = 0.0
         df4_cql_uprl = 0.0

         f5_cql_cart = 0.0
         df5_cql_uprp = 0.0
         df5_cql_uprl = 0.0

         f6_cql_cart = 0.0
         df6_cql_uprp = 0.0
         df6_cql_uprl = 0.0

         call blacs_barrier(icontxt, 'All')



         if (ndisti1 .eq. 1) then

            if(myid .eq. 0) call cql3d_setup1(netcdf_file1,nuper,nupar,
     &                                xmi1, enorm_factor, vc1_mks_cql3d)
            vc1_cgs_cql3d = vc1_mks_cql3d * 100.

            call blacs_barrier(icontxt, 'All')

            open(unit=50, file='cql3d.out', status='unknown',
     &                                                 form='formatted')

!           ---------------------------------------------
!           Read data from CQL3D in u_perp and u_parallel
!           ---------------------------------------------

            rewind (50)

            read (50, 309) nuper
            read (50, 309) nupar
            read (50, 309) n_psi

            read (50, 3310) vc1_mks
            read (50, 3310) UminPara, UmaxPara

            read (50, 3310) (rho_a(i_psi), i_psi = 1, n_psi)
            read (50, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            read (50, 3310) (upara(i_upara), i_upara = 1, nupar)

            read (50, 3310) (((f1_cql_cart(i_uperp, i_upara, i_psi),
     &       i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

            read (50, 3310) (((df1_cql_uprp(i_uperp, i_upara, i_psi),
     &       i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

            read (50, 3310) (((df1_cql_uprl(i_uperp, i_upara, i_psi),
     &       i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

             close(50)

         end if

         call blacs_barrier(icontxt, 'All')

         if (ndisti2 .eq. 1) then

            if(myid.eq.0) call cql3d_setup2(netcdf_file2, nuper, nupar,
     &                                xmi2, enorm_factor, vc2_mks_cql3d)

            vc2_cgs_cql3d = vc2_mks_cql3d * 100.

            call blacs_barrier(icontxt, 'All')

            open(unit=50,file='cql3d.out', status='unknown',
     &                                                 form='formatted')

!           ---------------------------------------------
!           Read data from CQL3D in u_perp and u_parallel
!           ---------------------------------------------

            rewind (50)

            read (50, 309) nuper
            read (50, 309) nupar
            read (50, 309) n_psi

            read (50, 3310) vc2_mks
            read (50, 3310) UminPara, UmaxPara

            read (50, 3310) (rho_a(i_psi), i_psi = 1, n_psi)
            read (50, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            read (50, 3310) (upara(i_upara), i_upara = 1, nupar)

            read (50, 3310) (((f2_cql_cart(i_uperp, i_upara, i_psi),
     &       i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

            read (50, 3310) (((df2_cql_uprp(i_uperp, i_upara, i_psi),
     &       i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

            read (50, 3310) (((df2_cql_uprl(i_uperp, i_upara, i_psi),
     &       i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)
            close(50)

         end if



         if(myid .eq. 0)then
            WRITE (6,*)
            WRITE (6,*) "nuper = ",  nuper
            WRITE (6,*) "nupar = ",  nupar
            WRITE (6,*) "n_psi = ",  n_psi

            WRITE (6,*)
            WRITE (6,*) "vce_mks = ", vce_mks
            WRITE (6,*) "vc1_mks = ", vc1_mks
            WRITE (6,*) "vc2_mks = ", vc2_mks
            WRITE (6,*) "vc3_mks = ", vc3_mks
            WRITE (6,*) "vc4_mks = ", vc4_mks
            WRITE (6,*) "vc5_mks = ", vc5_mks
            WRITE (6,*) "vc6_mks = ", vc6_mks

            write(6, *)
            write(6, *) "rho/a(i_psi)"
            write(6, 310) (rho_a(i_psi), i_psi = 1, n_psi)

            write(6, *)
            write(6, *) "uperp ="
            write(6, 310) (uperp(i_uperp), i_uperp = 1, nuper)

            write(6, *)
            write(6, *) "upara ="
            write(6, 310) (upara(i_upara), i_upara = 1, nupar)



            WRITE (15,*)
            WRITE (15,*) "nuper = ",  nuper
            WRITE (15,*) "nupar = ",  nupar
            WRITE (15,*) "n_psi = ",  n_psi

            WRITE (15,*)
            WRITE (15,*) "vce_mks = ", vce_mks
            WRITE (15,*) "vc1_mks = ", vc1_mks
            WRITE (15,*) "vc2_mks = ", vc2_mks
            WRITE (15,*) "vc3_mks = ", vc3_mks
            WRITE (15,*) "vc4_mks = ", vc4_mks
            WRITE (15,*) "vc5_mks = ", vc5_mks
            WRITE (15,*) "vc6_mks = ", vc6_mks

            write(15, *)
            write(15, *) "rho/a(i_psi)"
            write(15, 310) (rho_a(i_psi), i_psi = 1, n_psi)

            write(15, *)
            write(15, *) "uperp ="
            write(15, 310) (uperp(i_uperp), i_uperp = 1, nuper)

            write(15, *)
            write(15, *) "upara ="
            write(15, 310) (upara(i_upara), i_upara = 1, nupar)
         end if
      else

c        Allocate as single element arrays so compiler warnings can be
c        enabled. This should really be fixed properly, but lets just
c        get the compiler warnings turned on first

         allocate( fe_cql_cart (1,1,1) )
         allocate( dfe_cql_uprp(1,1,1) )
         allocate( dfe_cql_uprl(1,1,1) )

         allocate( f1_cql_cart (1,1,1) )
         allocate( df1_cql_uprp(1,1,1) )
         allocate( df1_cql_uprl(1,1,1) )

         allocate( f2_cql_cart (1,1,1) )
         allocate( df2_cql_uprp(1,1,1) )
         allocate( df2_cql_uprl(1,1,1) )

         allocate( f3_cql_cart (1,1,1) )
         allocate( df3_cql_uprp(1,1,1) )
         allocate( df3_cql_uprl(1,1,1) )

         allocate( f4_cql_cart (1,1,1) )
         allocate( df4_cql_uprp(1,1,1) )
         allocate( df4_cql_uprl(1,1,1) )

         allocate( f5_cql_cart (1,1,1) )
         allocate( df5_cql_uprp(1,1,1) )
         allocate( df5_cql_uprl(1,1,1) )

         allocate( f6_cql_cart (1,1,1) )
         allocate( df6_cql_uprp(1,1,1) )
         allocate( df6_cql_uprl(1,1,1) )

      end if

      xmax = rwright - rwleft
      ymax = ytop - ybottom

      if(ndisti1   .ne. 0 .or.
     &   ndisti2   .ne. 0 .or.
     &   ndisti3   .ne. 0 .or.
     &   ndisti4   .ne. 0 .or.
     &   ndisti5   .ne. 0 .or.
     &   ndisti6   .ne. 0 .or.
     &   ndiste    .ne. 0)   then

         deallocate( fe_cql_cart )
         deallocate( f1_cql_cart )
         deallocate( f2_cql_cart )
         deallocate( f3_cql_cart )
         deallocate( f4_cql_cart )
         deallocate( f5_cql_cart )
         deallocate( f6_cql_cart )
      end if

*--------------------------------------------
*--   Define x mesh: x(i), xprime(i), capr(i)
*--------------------------------------------
c--   xprime: 0 to xmax
c--   x(i) : -xmax / 2.0   to   xmax / 2.0
      dx = xmax / nnodex



      do i = 1, nnodex
         xprime(i) = (i-1) * dx
     &     + dx / 2.0
c--   Note: the code gives slightly smoother results with dx/2.0 added
         x(i) = xprime(i) + xwleft
         capr(i) = rt + x(i)
c         write(6, 1312)i, x(i), xprime(i), capr(i)

         xkphi(i) = nphi / capr(i)

      end do

      do i = 2, nnodex - 1
         drdx(i) = (capr(i+1) - capr(i-1) ) / (2.0 * dx)
         d2rdx2(i) = (capr(i+1) - 2.0 * capr(i) + capr(i-1)) / dx**2
      end do
      drdx(1) = (capr(2) - capr(1)) / dx
      d2rdx2(1) = 0.0
      drdx(nnodex) = (capr(nnodex) - capr(nnodex -1)) / dx
      d2rdx2(nnodex) = 0.0


      if(rzoom1 .eq. 0.0)rzoom1 = capr(1)
      if(rzoom2 .eq. 0.0)rzoom2 = capr(nnodex)

      izoom1 = int((rzoom1 - rwleft - dx / 2.0) / dx) + 1
      izoom2 = int((rzoom2 - rwleft - dx / 2.0) / dx) + 1
c      write(6, 1312) izoom1
c      write(6, 1312) izoom2

*     -----------------------------------
*     Define y mesh: y(j), yprime(j)
*     -----------------------------------
c--   yprime: 0 to ymax
c--   y(j) : -ymax / 2.0   to   ymax / 2.0
      dy = ymax / nnodey


      do j = 1, nnodey
         yprime(j) = (j-1) * dy
     &      + dy / 2.0
c--      Note: the code gives slightly smoother results with dy/2.0 added
         y(j) = yprime(j) + ybottom

c         write(6, 1312)j, y(j), yprime(j)
      end do





      if(yzoom1 .eq. 0.0)yzoom1 = ybottom
      if(yzoom2 .eq. 0.0)yzoom2 = ytop

      jzoom1 = int((yzoom1 - ybottom - dy / 2.0) / dy) + 1
      jzoom2 = int((yzoom2 - ybottom - dy / 2.0) / dy) + 1

      jdisp = int((ydisp - ybottom - dy / 2.0) / dy) + 1
      if(myid .eq. 0)write(6, *)  "ydisp = ", ydisp, "jdisp = ", jdisp

c      write(6, 1312) jzoom1
c      write(6, 1312) jzoom2


*--------------------------------------
*     find jmid = j in equatorial plane
*--------------------------------------
      ydiff_min = 100.0
      do j = 1, nnodey
         ydiff = abs(y(j) - zmaxis)
         if (ydiff .lt. ydiff_min) then
            jmid = j
            ydiff_min = ydiff
         end if
      end do

      jequat = jmid


*-------------------------------------------
*     define theta(i,j)
*        theta0(i,j) : -pi / 2.0  to  pi / 2.0
*        theta(i,j)  :       0.0  to  2.0 * pi
*-------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey
            if(y(j) .ne. 0.0 .or. x(i) .ne. 0.0) then
               theta0(i,j) = atan2(y(j), x(i))
               if(theta0(i,j) .ge. 0.0) theta(i,j) = theta0(i,j)
               if(theta0(i,j) .lt. 0.0) theta(i,j) =
     1                                       theta0(i,j) + 2.0 * pi
            end if
         end do
      end do

*----------------------------------------------------
*     Calculate derivatives of rho(i,j) and theta(i,j)
*-----------------------------------------------------
      drhodr = 0.
      drhodz = 0.
      dthedr = 0.
      dthedz = 0.

      do i = 2, nnodex - 1
         do j = 2, nnodey - 1
            drhodr(i,j) = (rho(i+1, j) - rho(i-1, j)) / (2.0 * dx)
            drhodz(i,j) = (rho(i, j+1) - rho(i, j-1)) / (2.0 * dy)
            dthedr(i,j) = (theta(i+1, j) - theta(i-1, j)) / (2.0 * dx)
            dthedz(i,j) = (theta(i, j+1) - theta(i, j-1)) / (2.0 * dy)
         end do
      end do


*-----------------------------
*     Define rho mesh: rhon(n)
*------------------------------
c--   rhon: 0 to rhomax
c      rhomax = 1.0
      drho = rhomax / (nnoderho - 1)
      do n = 1, nnoderho
         rhon(n) = (n-1) * drho
      end do

*---------------------------------
*     Define theta mesh: thetam(m)
*---------------------------------
c--   thetam: 0 to 2pi
      thetamax = 2.0 * pi
      dtheta = thetamax / (mnodetheta - 1)
      do m = 1, mnodetheta
         thetam(m) = (m-1) * dtheta
         if(myid .eq. 0)write(6, *)m, thetam(m)
      end do

*---------------------------------------------------
*     Invert mesh to get capr(nnoderho, nnodetheta)
*                    and capz(nnoderho, nnodetheta)
*---------------------------------------------------
      do n = 1, nnoderho
         do m = 1, mnodetheta
            rho_in = rhon(n)
            theta_in = thetam(m)

            call invert(rho_in, theta_in, capr_out, capz_out,
     &         nnodex, nnodey, nxmx, nymx, rho, theta, capr, y, myid,
     &         drhodr, drhodz, dthedr, dthedz, xj)


            capr_flux(n,m) = capr_out
            capz_flux(n,m) = capz_out
c           if(myid .eq. 0)write(6,*)n, m, capr_flux(n,m), xj
         end do
      end do


      delta = sqrt(dx**2 + dy**2)

      xant = rant - rt

      antang_rad = antang / 180. * pi
      sin_antang = sin(antang_rad)
      cos_antang = cos(antang_rad)

      iant=int((rant - rwleft) / dx) + 1
      if(rant .ne. 0.0)psiant = psi(iant, jequat)


c--note curden is in Amps per meter of toroidal length (2.*pi*rt).
      xjantx = curdnx / dx
      xjanty = curdny / dx
      xjantz = curdnz / dx

      xjant=sqrt(xjantx**2 + xjanty**2 + xjantz**2)

      rlim   = rt + alim






*     ----------------------------------------
*     Set up r mesh for circular flux surfaces
*     ----------------------------------------
      j = jequat
      do i = 2, nnodex -1
c         if (myid .eq. 0)write(6, 1312)i, rho(i,j), x(i), xprime(i),
c     &       capr(i)
         if(rho(i,j) .lt. 1. .and. rho(i-1,j) .ge. 1.)caprmin = capr(i)
         if(rho(i,j) .gt. 1. .and. rho(i-1,j) .le. 1.)caprmax =capr(i-1)
      end do

      rmax = (caprmax - caprmin) / 2.0

      if(myid .eq. 0)write(6, 311)caprmax, caprmin, rmax



      do n = 1, nnoderho
         rn(n) = rhon(n) * rmax
         epsn(n) = rn(n) / rt
         dn(n) = rmax**2 / taue
c         if (myid .eq. 0)write(6, 1312)n, rn(n), rhon(n), dn(n)
c         if (myid .eq. 0)write(15, 1312)n, rn(n), rhon(n), dn(n)
      end do

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1
            if (n .eq. 1) then
               i0 = i
               j0 = j
            end if
         end do
      end do

      if(i0 .eq. 0 .and. j0 .eq. 0) then
         do i = 1, nnodex
            do j = 1, nnodey
               n = int(rho(i,j) / drho) + 1
               if (n .eq. 2) then
                  i0 = i
                  j0 = j
               end if
            end do
         end do

      end if


      dr = drho * rmax

      do n = 1, nnoderho - 1
         rh(n) = (rn(n) + rn(n+1)) / 2.0
         dh(n) = (dn(n) + dn(n+1)) / 2.0
c         if (myid .eq. 0)write(6, 1312)n, rh(n), dh(n)
      end do
      rh(nnoderho) = 0.0
      dh(nnoderho) = 0.0



*     -------------------------------------
*     Calculate capr * bpol in the midplane
*     -------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey
            capr_bpol(i,j) = capr(i) * bpol(i,j)
            ipsi(i,j)= capr(i) * bzeta(i,j) * bmod(i,j)
         end do
      end do

      do i = 1, nnodex
         do j = 1, nnodey
            call midplane(i, j, capr_bpol, capr_bpol_mid2(i,j),
     &          rho, nxmx, nymx, nnodex, nnodey, capr, rt, 0.0, jmid)
           end do
      end do

      call polavg(capr_bpol_mid2, capr_bpol_mid, rho, nxmx, nymx,
     &   nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvolp, fvol)


      if (myid .eq. 0)then

         write(6, *) "jmid = ", jmid
         write(15, *) "jmid = ", jmid

         write(6, *)  "psio = ", psio
         write(15, *) "psio = ", psio
         write(6, *)  "psi_tor_max = ", psi_tor_max
         write(15, *) "psi_tor_max = ", psi_tor_max

         write(15, *)"rt = ", rt
         write(15, *)"b0 = ", b0
         write(15, *)"dx = ", dx
         write(15, *)"dy = ", dy
         write(15, *)"drho = ", drho

         write(15, *)
         write(15, *) "     i    capr  rhoij  capr_bpol  capr_bpol_mid2"
         write(15, *)

         write(6, *)
         write(6, *)  "     i    capr  rhoij  capr_bpol  capr_bpol_mid2"
         write(6, *)

         do i = 1, nnodex
            write(6,  1312)i, capr(i), rho(i, jequat),
     &                   capr_bpol(i, jequat), capr_bpol_mid2(i, jequat)
            write(15, 1312)i, capr(i), rho(i, jequat),
     &                   capr_bpol(i, jequat), capr_bpol_mid2(i, jequat)
         end do

         write(15, *)
         write(15, *) "     n   capr_bpol_mid"

         write(6, *)
         write(6, *)  "     n   capr_bpol_mid"

         do n = 1, nnoderho
            write(6,  1312)n, capr_bpol_mid(n)
            write(15, 1312)n, capr_bpol_mid(n)
         end do

      end if




*     -----------------------
*     Calculate bmod_mid(i,j)
*     -----------------------
      do i = 1, nnodex
         do j = 1, nnodey
            call midplane(i, j, bmod, bmod_mid(i,j), rho,
     &                nxmx, nymx, nnodex, nnodey, capr, rt, b0, jmid)
           end do
      end do

      call polavg(bmod_mid, bmod_midavg, rho, nxmx, nymx,
     &   nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvolp, fvol)


c      if (myid .eq. 0)then

c         write(15, *)
c         write(15, *) "     n        bmod_midavg"

c         write(6, *)
c         write(6, *)  "     n        bmod_midavg"


c         do n = 1, nnoderho
c            write(6,  1312)n, bmod_midavg(n)
c            write(15, 1312)n, bmod_midavg(n)
c         end do

c      end if



c      if (myid .eq. 0)then

c         write(6, *)
c         write(15, *)

c         do i = 1, nnodex
c            do j = 1, nnodey
c              if(j .eq. jequat) then
c                 write(6,  1312)i, x(i), bmod(i, j), bmod_mid(i, j)
c                 write(15, 1312)i, x(i), bmod(i, j), bmod_mid(i, j)
c              end if
c            end do
c         end do

c         do i = 1, nnodex
c            do j = 1, nnodey
c               if (i .eq. 5 .and. j .eq. 5)then
c                  write(6, *) "bmod_mid(i, j) = ", bmod_mid(i, j)
c                  write(6, *) "bxn(i, j) = ", bxn(i, j)
c                  write(6, *) "byn(i, j) = ", byn(i, j)
c                  write(6, *) "bzn(i, j) = ", bzn(i, j)
c               end if
c            end do
c         end do

c      end if



*     ---------------
*     Antenna current
*     ---------------

      theta_antr = theta_ant / 180. * pi

      dpsiant = dpsiant0
      dthetant = dthetant0 / 360. * 2.0 * pi
      yant_max = antlen / 2.0


      do i = 1, nnodex
         do j = 1, nnodey

            call deriv_x(psi, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dpsidx, d2psidx2)
            call deriv_y(psi, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dpsidy, d2psidy2)

            gradpsi = sqrt(dpsidx**2 + dpsidy**2)
            if (gradpsi .le. 1.0e-06)gradpsi = 1.0e-06


            delta_theta = theta0(i,j) - theta_antr
            gaussantth = exp(-delta_theta**2 / dthetant**2)

*           ------------------------
*           Curved antenna (default)
*           ------------------------
            if(curved_ant .eqv. .true.) then
               gausspsi = exp(-( psi(i,j) - psiant)**2 / dpsiant**2)
            end if

*           ----------------
*           Straight antenna
*           -----------------
            if(curved_ant .neqv. .true.) then
               gausspsi   = exp(-(capr(i) -   rant)**2 / dpsiant**2)
            end if

            shapey = 0.0
            deltay = y(j) - yant

            if(i_antenna .eq. 0 .or. i_antenna .eq. 1) then

               if (capr(i) .gt. rt) then

*                 ------------------------------------------------
*                 if(i_antenna .eq. 0) antenna current is Gaussian
*                 ------------------------------------------------

                  if (i_antenna .eq. 0) then
                     shapey = exp(- deltay**2 / yant_max**2)
                  end if


*                 ---------------------------------------------------
*                 if(i_antenna .eq. 1) antenna current is cos(ky * y)  (DEFAULT)
*                 ---------------------------------------------------
c                  antlc = pi / (xk0 * antlen)
                  if (i_antenna.eq.1 .and. abs(deltay).lt.yant_max)then
                     shapey = cos(xk0 * antlc * deltay)
                  end if


               end if

               Ix = - dpsidy / gradpsi * cos_antang
               Iy =   dpsidx / gradpsi * cos_antang
               Iz =   sin_antang

               xjx(i,j) = xjantx * shapey  * gausspsi * Ix
               xjy(i,j) = xjanty * shapey  * gausspsi * Iy
               xjz(i,j) = xjantz * shapey  * gausspsi * Iz

            end if

            if(i_antenna .eq. 2) then
               call random_number(ranval)
               xjx(i,j) = 0
               xjy(i,j) = ranval * 1.0e-06
               xjz(i,j) = 0
            end if

            if(i_antenna .eq. 3) then
               call random_number(ranval)
               xjx(i,j) = ranval * 1.0e-06
               xjy(i,j) = 0
               xjz(i,j) = 0
            end if

            if(i_antenna .eq. 4) then
               call random_number(ranval)
               xjx(i,j) = 0
               xjy(i,j) = 0
               xjz(i,j) = ranval * 1.0e-06
            end if


         end do
      end do

*     ----------------------------
*     find maximum antenna current
*     ----------------------------
      cur_mod_max = 0.0
      do i = 1, nnodex
         do j = 1, nnodey
            cur_mod = sqrt(xjx(i,j)**2 + xjy(i,j)**2 + xjz(i,j)**2)
            if (cur_mod .gt. cur_mod_max) then
               cur_mod_max = cur_mod
               i_cur_max = i
               j_cur_max = j
            end if
         end do
      end do

      xkphia = nphi / capr(i_cur_max)
      npar = xkphia / xk0

      if(myid .eq. 0)then
         write (15, 1022) npar
         write ( 6, 1022) npar
      end if


*    ----------------
*    plasma profiles:
*    ----------------
      do i = 1, nnodex
         do j = 1, nnodey


            if(psimol .eq. 1.0)flimiter = 1.0
            if(psimol .ne. 1.0)
     &         flimiter = 1.0 / (1.0 + (psi(i,j) / psimol)**16)

            if (iprofile .eq. 3) then

               shapen = 0.0
               shapen2 = 0.0
               shapen3 = 0.0
               shapen4 = 0.0
               shapen5 = 0.0
               shapen6 = 0.0
               shapen_slo = 0.0

               shapete = 0.0
               shapeti = 0.0
               shapeti2 = 0.0
               shapeti3 = 0.0
               shapeti4 = 0.0
               shapeti5 = 0.0
               shapeti6 = 0.0


               if(rho(i,j) .le. rholim)then

                  shapen  = 1.0 - rho(i,j)**betan
                  shapen2 = 1.0 - rho(i,j)**betan2
                  shapen3 = 1.0 - rho(i,j)**betan3
                  shapen4 = 1.0 - rho(i,j)**betan4
                  shapen5 = 1.0 - rho(i,j)**betan5
                  shapen6 = 1.0 - rho(i,j)**betan6
                  shapen_slo =  1.0 - rho(i,j)**betan_slo

                  shapete  = 1.0 - rho(i,j)**betate
                  shapeti  = 1.0 - rho(i,j)**betati
                  shapeti2 = 1.0 - rho(i,j)**betati2
                  shapeti3 = 1.0 - rho(i,j)**betati3
                  shapeti4 = 1.0 - rho(i,j)**betati4
                  shapeti5 = 1.0 - rho(i,j)**betati5
                  shapeti6 = 1.0 - rho(i,j)**betati6


               end if

               xnea(i,j) =
     &            xnlim  + (xn0 - xnlim)  * shapen**alphan   * flimiter
               xn2a(i,j) =
     &            xn2lim + (xn2 - xn2lim) * shapen2**alphan2 * flimiter
               xn3a(i,j) =
     &            xn3lim + (xn3 - xn3lim) * shapen3**alphan3 * flimiter
               xn4a(i,j) =
     &            xn4lim + (xn4 - xn4lim) * shapen4**alphan4 * flimiter
               xn5a(i,j) =
     &            xn5lim + (xn5 - xn5lim) * shapen5**alphan5 * flimiter
               xn6a(i,j) =
     &            xn6lim + (xn6 - xn6lim) * shapen6**alphan6 * flimiter
               xna_slo(i,j) =
     &            xnslolim + (xnslo - xnslolim) * shapen_slo**alphan_slo
     &                                                        * flimiter

               nmin2 = xn2lim / xnlim * nmin
               nmin3 = xn3lim / xnlim * nmin
               nmin4 = xn4lim / xnlim * nmin
               nmin5 = xn5lim / xnlim * nmin
               nmin6 = xn6lim / xnlim * nmin
               nminslo = xnslolim / xnlim * nmin


               if(rho(i,j) .gt. rholim) then
                  drho1 = rho(i,j) - rholim
                  xnea(i,j) = nmin +
     &                  (xnlim - nmin) * exp(- drho1 / scrape)
                  xn2a(i,j) = nmin2 +
     &                  (xn2lim - nmin2) * exp(- drho1 / scrape)
                  xn3a(i,j) = nmin3 +
     &                  (xn3lim - nmin3) * exp(- drho1 / scrape)
                  xn4a(i,j) = nmin4 +
     &                  (xn4lim - nmin4) * exp(- drho1 / scrape)
                  xn5a(i,j) = nmin5 +
     &                  (xn5lim - nmin6) * exp(- drho1 / scrape)
                  xn6a(i,j) = nmin6 +
     &                  (xn6lim - nmin6) * exp(- drho1 / scrape)
                  xna_slo(i,j) = nminslo +
     &                  (xnslolim - nminslo) * exp(- drho1 / scrape)
              end if


               xn1a(i, j) = (xnea(i, j) - z2 * xn2a(i, j)
     &                                  - z3 * xn3a(i, j)
     &                                  - z4 * xn4a(i, j)
     &                                  - z5 * xn5a(i, j)
     &                                  - z6 * xn6a(i, j)
     &                               - z_slo * xna_slo(i, j)   ) / z1

               if (xn1a(i, j) .le. 0.0) then
                  xn1a(i, j) = 1.0e-10
                  if(myid .eq.0) write(6, *) xnea(i,j), xn1a(i,j),
     &               xn2a(i,j), xn3a(i,j), xn4a(i,j), xn5a(i,j),
     &               xn6a(i,j), xna_slo(i,j)
                  stop 3
               end if



               eta2 = xn2 / xn0
               eta3 = xn3 / xn0
               eta4 = xn4 / xn0
               eta5 = xn5 / xn0
               eta6 = xn6 / xn0
               eta_slo = xnslo / xn0

               eta1 = 1.0 / z1 * (1.0 - z2 * eta2 - z3 * eta3
     &                    - z4 * eta4 - z5 * eta5 - z6 * eta6
     &                                            - z_slo * eta_slo)
               xn1 = eta1 * xn0

               xkte(i,j) =
     &            telimj + (t0e - telimj) * shapete**alphate * flimiter
               xkti(i,j) =
     &            tilimj + (t0i - tilimj) * shapeti**alphati * flimiter
               xkti2(i,j) =
     &            ti2limj + (t0i2 - ti2limj) * shapeti2**alphati2
     &                                                      * flimiter
               xkti3(i,j) =
     &            ti3limj + (t0i3 - ti3limj) * shapeti3**alphati3
     &                                                      * flimiter
               xkti4(i,j) =
     &            ti4limj + (t0i4 - ti4limj) * shapeti4**alphati4
     &                                                      * flimiter
               xkti5(i,j) =
     &            ti5limj + (t0i5 - ti5limj) * shapeti5**alphati5
     &                                                      * flimiter
               xkti6(i,j) =
     &            ti6limj + (t0i6 - ti6limj) * shapeti6**alphati6
     &                                                      * flimiter

               Temin0 = Temin * q * 1.0e+3
               Tmini =  tilimj  / telimj * Temin0
               Tmini2 = ti2limj / telimj * Temin0
               Tmini3 = ti3limj / telimj * Temin0
               Tmini4 = ti4limj / telimj * Temin0
               Tmini5 = ti5limj / telimj * Temin0
               Tmini6 = ti6limj / telimj * Temin0

               if(rho(i,j) .gt. rholim) then
                  xkte(i,j)  = telimj
                  xkti(i,j)  = tilimj
                  xkti2(i,j) = ti2limj
                  xkti3(i,j) = ti3limj
                  xkti4(i,j) = ti4limj
                  xkti5(i,j) = ti5limj
                  xkti6(i,j) = ti6limj
                  drho1 = rho(i,j) - rholim
c                  xkte(i,j) = Temin0 +
c     &                  (telimj - Temin0) * exp(- drho1 / scrape)
c                 xkti(i,j) = Tmini +
c     &                  (tilimj - Tmini) * exp(- drho1 / scrape)
c                 xkti2(i,j) = Tmini2 +
c     &                  (ti2limj - Tmini2) * exp(- drho1 / scrape)
c                 xkti3(i,j) = Tmini3 +
c     &                  (ti3limj - Tmini3) * exp(- drho1 / scrape)
c                 xkti4(i,j) = Tmini4 +
c     &                  (ti4limj - Tmini4) * exp(- drho1 / scrape)
c                 xkti5(i,j) = Tmini5 +
c     &                  (ti5limj - Tmini5) * exp(- drho1 / scrape)
c                 xkti6(i,j) = Tmini6 +
c     &                  (ti6limj - Tmini6) * exp(- drho1 / scrape)
               end if

            end if
*          --------------------
*           end iprofile .eq. 3:
*          --------------------


            if (iprofile .eq. 1 .or. iprofile .eq. 2) then


               xn1 = xn0 / z1 * (1.0 - z2 * eta  - z3 * eta3
     &                   - z4 * eta4 - z5 * eta5 - z6 * eta6
     &                                             - z_slo * eta_slo)
               eta1 = xn1 / xn0
               xn20 = xn0 * eta
               eta2 = xn20 / xn0
               xn3 = xn0 * eta3
               xn4 = xn0 * eta4
               xn5 = xn0 * eta5
               xn6 = xn0 * eta6

               xn_slo = xn0 * eta_slo

               gaussian_ne  = exp(-psi(i,j) / psipne)
               gaussian_te  = exp(-psi(i,j) / psipte)
               gaussian_ti1 = exp(-psi(i,j) / psipti1)
               gaussian_ti2 = exp(-psi(i,j) / psipti2)
               gaussian_ti3 = exp(-psi(i,j) / psipti3)
               gaussian_ti4 = exp(-psi(i,j) / psipti4)
               gaussian_ti5 = exp(-psi(i,j) / psipti5)
               gaussian_ti6 = exp(-psi(i,j) / psipti6)

               parabola = (1. - psi(i,j) / psilim)
               if (parabola .le. 0.0) parabola = 0.0


               if(iprofile .eq. 1) shape = gaussian_ne
               if(iprofile .eq. 2) shape = parabola * flimiter


               xnea(i,j)  = xnlim +  (xn0 - xnlim) * shape**alphan
               xn2a(i, j) = xnea(i, j) * eta2
               xn3a(i, j) = xnea(i, j) * eta3
               xn4a(i, j) = xnea(i, j) * eta4
               xn5a(i, j) = xnea(i, j) * eta5
               xn6a(i, j) = xnea(i, j) * eta6
               xna_slo(i, j) = xnea(i,j) * eta_slo

               xn1a(i, j) = (xnea(i, j) - z2 * xn2a(i, j)
     &                                  - z3 * xn3a(i, j)
     &                                  - z4 * xn4a(i, j)
     &                                  - z5 * xn5a(i, j)
     &                                  - z6 * xn6a(i, j)
     &                               - z_slo * xna_slo(i, j)   ) / z1

               if (xn1a(i, j) .le. 0.0) xn1a(i, j) = 1.0e-10


               xkte(i,j)  = telimj +  (t0e - telimj) * shape**alphate
               xkti(i,j)  = tilimj +  (t0i - tilimj) * shape**alphati
               xkti2(i,j) = ti2limj + (t0i2 - ti2limj) * shape**alphati
               xkti3(i,j) = ti3limj + (t0i3 - ti3limj) * shape**alphati
               xkti4(i,j) = ti4limj + (t0i4 - ti4limj) * shape**alphati
               xkti5(i,j) = ti5limj + (t0i5 - ti5limj) * shape**alphati
               xkti6(i,j) = ti6limj + (t0i6 - ti6limj) * shape**alphati

            end if

         end do
      end do



*     ---------------------------------------
*     read profile data from the plasma state:
*     ---------------------------------------
      if (iprofile .eq. 5 .or. iprofile .eq. 6) then
*        ----------------------------------------------------------
*        iprofile = 5: numerical profiles for core but analytic SOL
*        iprofile = 6: numerical profiles for both core and SOL
*        ----------------------------------------------------------
*        ------------------
*        read profile data:
*        ------------------
         rewind (63)
         read(63, state)

         ndata = s_nrho_n


*        ------------------------------------------------------------
*        Optional: change from stored energy to equivalent kT/e (keV)
*        ------------------------------------------------------------
c        do n = 1, ndata
c            s_t_s(n,2) = (s_t_s(n,2) + s_t_s(n,3)) * 1.0 e+06
c     &                         / q / s_n_s(n,2) * 2./ 3./ 1.e+03
c            s_t_s(n,3) = s_t_s(n,3) / q / s_n_s(n,3) * 2./ 3./ 1.e+03
c            s_t_s(n,4) = s_t_s(n,4) / q / s_n_s(n,4) * 2./ 3./ 1.e+03
c            s_t_s(n,0) = s_t_s(n,0) * 2.5
c        end do



*        -------------------------
*        Optional: modify profiles
*        -------------------------

         go to 5009
         rho0 = 0.90
         n1 = 14
         do n = 1, ndata
            rhod(n) = s_rho_n_grid(n)

c           s_n_s(n,0) = s_n_s(n,0) / (1. + (rhod(n) / rho0)**n1)
c           s_n_s(n,1) = s_n_s(n,1) / (1. + (rhod(n) / rho0)**n1)
c           s_n_s(n,2) = s_n_s(n,2) / (1. + (rhod(n) / rho0)**n1)
c           s_n_s(n,3) = s_n_s(n,3) / (1. + (rhod(n) / rho0)**n1)
c           s_n_s(n,4) = s_n_s(n,4) / (1. + (rhod(n) / rho0)**n1)

            s_n_s(n,0) = s_n_s(n,0) * 0.54112554
            s_n_s(n,1) = s_n_s(n,1) * 0.54112554


c           s_t_s(n,0) = s_t_s(n,0) / (1. + (rhod(n) / rho0)**n1)
c           s_t_s(n,1) = s_t_s(n,1) / (1. + (rhod(n) / rho0)**n1)
c           s_t_s(n,2) = s_t_s(n,2) / (1. + (rhod(n) / rho0)**n1)
c           s_t_s(n,3) = s_t_s(n,3) / (1. + (rhod(n) / rho0)**n1)
c           s_t_s(n,4) = s_t_s(n,4) / (1. + (rhod(n) / rho0)**n1)

            s_t_s(n,0) = s_t_s(n,0) * 0.45801527
            s_t_s(n,1) = s_t_s(n,1) * 0.45801527
         end do
 5009    continue


         if(myid .eq. 0)then
            write(6, *)
            write(6, *) "ndata = ", s_nrho_n
            write(6, 66164)
66164       format(" " ,5x,"   n  ",3x," rhod ",
     &                 6x,"  ne  ",6x,"  te  ",6x,"  ti1 ",6x," ti2  ",
     &                 6x,"  ti3 ",6x," zeff ",6x,"  n1  ",6x,"  n2  ",
     &                 6x,"  n3  ")
            write(6, *)


            write(15, *)
            write(15, *) "ndata = ", s_nrho_n
            write(15, 66164)
            write(15, *)

            do n = 1, s_nrho_n
               write(6, 1312)n, s_rho_n_grid(n),
     &             s_n_s(n,0), s_t_s(n,0), s_t_s(n,1), s_t_s(n,2),
     &             s_t_s(n,3), s_zeff(n),  s_n_s(n,1), s_n_s(n,2),
     &             s_n_s(n,3)
               write(15, 1312)n, s_rho_n_grid(n),
     &             s_n_s(n,0), s_t_s(n,0), s_t_s(n,1), s_t_s(n,2),
     &             s_t_s(n,3), s_zeff(n),  s_n_s(n,1), s_n_s(n,2),
     &             s_n_s(n,3)
            end do

            ndata = s_nrho_n

            do n = 1, ndata

               rhod(n) = s_rho_n_grid(n)

               tekev(n)  = s_t_s(n,0)
               tikev(n)  = s_t_s(n,1)
               ti2kev(n) = s_t_s(n,2)
               e_beam(n) = s_t_s(n,3)

               xneavg(n) = s_n_s(n,0)
               xn_min(n) = s_n_s(n,2)
               xn_maj(n) = s_n_s(n,1)
               xn_beam(n)= s_n_s(n,3)





               if(eta .gt. 0.0)then
                  xn_min(n) = eta * xneavg(n)
               end if

               if(eta3 .gt. 0.0)then
                  xn_beam(n) = eta3 * xneavg(n)
               end if

               xn_maj(n) = (xneavg(n) - z2 * xn_min(n)
     &                                - z3 * s_n_s(n,3)
     &                                - z4 * s_n_s(n,4)
     &                                - z5 * s_n_s(n,5)
     &                                - z6 * s_n_s(n,6) ) / z1

               xneavg_cgs(n) = xneavg(n) * 1.0e-06
               xn_min_cgs(n) = xn_min(n) * 1.0e-06
               xn_maj_cgs(n) = xn_maj(n) * 1.0e-06
               xn_beam_cgs(n) = xn_beam(n) * 1.0e-06

               zeffavg(n) = s_zeff(n)

c              rhod(n) = sqrt(rhod(n))
            end do

c           tekev  = tekev  / 3.
c           tikev  = tikev  / 3.
c           ti2kev = ti2kev / 3.

            write(6, *)
            write(6, *) "rhod = "
            write(6, 100) (rhod(n), n = 1, ndata)
            write(6, *)

            write(6, *) "xneavg = "
            write(6, 100) (xneavg(n), n = 1, ndata)
            write(6, *)

            write(6, *) "xn_majority = "
            write(6, 100) (xn_maj(n), n = 1, ndata)
            write(6, *)

            write(6, *) "xn_minority = "
            write(6, 100) (xn_min(n), n = 1, ndata)
            write(6, *)

            write(6, *) "xn_beam = "
            write(6, 100) (xn_beam(n), n = 1, ndata)
            write(6, *)

            write(6, *) "xn_ion4 = "
            write(6, 100) (s_n_s(n,4), n = 1, ndata)
            write(6, *)

            write(6, *) "xn_ion5 = "
            write(6, 100) (s_n_s(n,5), n = 1, ndata)
            write(6, *)

            write(6, *) "xn_ion6 = "
            write(6, 100) (s_n_s(n,6), n = 1, ndata)
            write(6, *)

            write(6, *) "zeffavg = "
            write(6, 100) (zeffavg(n), n = 1, ndata)
            write(6, *)

c           tekev = tekev * 2.1 / 3.54298
c           tikev = tikev * 2.1 / 2.83439

            write(6, *) "tekev = "
            write(6, 100) (tekev(n), n = 1, ndata)
            write(6, *)

            write(6, *) "tikev = "
            write(6, 100) (tikev(n), n = 1, ndata)
            write(6, *)

            write(6, *) "T_ion2 (keV) = "
            write(6, 100) (s_t_s(n,2), n = 1, ndata)
            write(6, *)

            write(6, *) "T_ion3 (keV) = "
            write(6, 100) (s_t_s(n,3), n = 1, ndata)
            write(6, *)

            write(6, *) "T_ion4 (keV) = "
            write(6, 100) (s_t_s(n,4), n = 1, ndata)
            write(6, *)

            write(6, *) "T_ion5 (keV) = "
            write(6, 100) (s_t_s(n,5), n = 1, ndata)
            write(6, *)

            write(6, *) "T_ion6 (keV) = "
            write(6, 100) (s_t_s(n,6), n = 1, ndata)
            write(6, *)

            write(6, *) "xn_majority_cgs = "
            write(6, 100) (xn_maj_cgs(n), n = 1, ndata)
            write(6, *)

            write(6, *) "xn_minority_cgs = "
            write(6, 100) (xn_min_cgs(n), n = 1, ndata)
            write(6, *)

            write(6, *) "xneavg_cgs = "
            write(6, 100) (xneavg_cgs(n), n = 1, ndata)
            write(6, *)


            write(6, *) "ti2kev = "
            write(6, 100) (ti2kev(n), n = 1, ndata)
            write(6, *)

            write(6, *) "xme = ", xme
            write(6, *) "xmi1 = ", xmi1
            write(6, *) "xmi2 = ", xmi2
            write(6, *) "xmi3 = ", xmi3
            write(6, *) "xmi4 = ", xmi4
            write(6, *) "xmi5 = ", xmi5
            write(6, *) "xmi6 = ", xmi6

            write(6, *) "qe = ", qe
            write(6, *) "qi1 = ", qi1
            write(6, *) "qi2 = ", qi2
            write(6, *) "qi3 = ", qi3
            write(6, *) "qi4 = ", qi4
            write(6, *) "qi5 = ", qi5
            write(6, *) "qi6 = ", qi6

            write(15, *)
            write(15, *) "rhod = "
            write(15, 100) (rhod(n), n = 1, ndata)
            write(15, *)

            write(15, *) "xneavg = "
            write(15, 100) (xneavg(n), n = 1, ndata)
            write(15, *)

            write(15, *) "xn_majority = "
            write(15, 100) (xn_maj(n), n = 1, ndata)
            write(15, *)

            write(15, *) "xn_minority = "
            write(15, 100) (xn_min(n), n = 1, ndata)
            write(15, *)

            write(15, *) "xn_beam = "
            write(15, 100) (xn_beam(n), n = 1, ndata)
            write(15, *)

            write(15, *) "zeffavg = "
            write(15, 100) (zeffavg(n), n = 1, ndata)
            write(15, *)

            write(15, *) "tekev = "
            write(15, 100) (tekev(n), n = 1, ndata)
            write(15, *)

            write(15, *) "tikev = "
            write(15, 100) (tikev(n), n = 1, ndata)
            write(15, *)

            write(15, *) "e_beam = "
            write(15, 100) (e_beam(n), n = 1, ndata)
            write(15, *)

            write(15, *) "xn_majority_cgs = "
            write(15, 100) (xn_maj_cgs(n), n = 1, ndata)
            write(15, *)

            write(15, *) "xn_minority_cgs = "
            write(15, 100) (xn_min_cgs(n), n = 1, ndata)
            write(15, *)

            write(15, *) "xneavg_cgs = "
            write(15, 100) (xneavg_cgs(n), n = 1, ndata)
            write(15, *)


            write(15, *) "ti2kev = "
            write(15, 100) (ti2kev(n), n = 1, ndata)
            write(15, *)

            write(15, *) "xme = ", xme
            write(15, *) "xmi1 = ", xmi1
            write(15, *) "xmi2 = ", xmi2
            write(15, *) "xmi3 = ", xmi3
            write(15, *) "xmi4 = ", xmi4
            write(15, *) "xmi5 = ", xmi5
            write(15, *) "xmi6 = ", xmi6

            write(15, *) "qe = ", qe
            write(15, *) "qi1 = ", qi1
            write(15, *) "qi2 = ", qi2
            write(15, *) "qi3 = ", qi3
            write(15, *) "qi4 = ", qi4
            write(15, *) "qi5 = ", qi5
            write(15, *) "qi6 = ", qi6



         end if
  100    format(1p,6e12.3)


         s_t_s = q * s_t_s * 1.0e3    !input is kev TEMPS IN JOULES
         Temin0 = q * Temin * 1.0e3    !input is kev TEMPS IN JOULES

         nmin1 = s_n_s(ndata, 1) / s_n_s(ndata, 0) * nmin
         nmin2 = s_n_s(ndata, 2) / s_n_s(ndata, 0) * nmin
         nmin3 = s_n_s(ndata, 3) / s_n_s(ndata, 0) * nmin
         nmin4 = s_n_s(ndata, 4) / s_n_s(ndata, 0) * nmin
         nmin5 = s_n_s(ndata, 5) / s_n_s(ndata, 0) * nmin
         nmin6 = s_n_s(ndata, 6) / s_n_s(ndata, 0) * nmin

         Tmin1 = s_t_s(ndata, 1) / s_t_s(ndata, 0) * Temin0
         Tmin2 = s_t_s(ndata, 2) / s_t_s(ndata, 0) * Temin0
         Tmin3 = s_t_s(ndata, 3) / s_t_s(ndata, 0) * Temin0
         Tmin4 = s_t_s(ndata, 4) / s_t_s(ndata, 0) * Temin0
         Tmin5 = s_t_s(ndata, 5) / s_t_s(ndata, 0) * Temin0
         Tmin6 = s_t_s(ndata, 6) / s_t_s(ndata, 0) * Temin0



*        -----------------------------
*        deposit data on 2D flux grid:
*        -----------------------------

         xnea = 1.0e-10
c        if(nboundary .eq. 0)xnea = s_n_s(ndata,0)

         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,0),
     &        xnea(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xkte = 1.0e-10 * q
c        if(nboundary .eq. 0)xkte = s_t_s(ndata,0)

         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,0),
     &        xkte(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xn1a = 1.0e-10
c        if(nboundary .eq. 0)xn1a = s_n_s(ndata,1)

         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,1),
     &        xn1a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xkti = 1.0e-10 *q
c        if(nboundary .eq. 0)xkti = s_t_s(ndata,1)


         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,1),
     &        xkti(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xn2a = 1.0e-10
c        if(nboundary .eq. 0)xn2a = s_n_s(ndata,2)

         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,2),
     &        xn2a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xkti2 = 1.0e-10 *q
c        if(nboundary .eq. 0)xkti2 = s_t_s(ndata,2)

         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,2),
     &        xkti2(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xn3a = 1.0e-10
c        if(nboundary .eq. 0)xn3a = s_n_s(ndata,3)

         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,3),
     &        xn3a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xkti3 = 1.0e-10 *q
c        if(nboundary .eq. 0)xkti3 = s_t_s(ndata,3)

         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,3),
     &        xkti3(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xn4a = 1.0e-10
c        if(nboundary .eq. 0)xn4a = s_n_s(ndata,4)

         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,4),
     &        xn4a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xkti4 = 1.0e-10 *q
c        if(nboundary .eq. 0)xkti4 = s_t_s(ndata,4)

         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,4),
     &        xkti4(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xn5a = 1.0e-10
c        if(nboundary .eq. 0)xn5a = s_n_s(ndata,5)

         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,5),
     &        xn5a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xkti5 = 1.0e-10 *q
c        if(nboundary .eq. 0)xkti5 = s_t_s(ndata,5)

         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,5),
     &        xkti5(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xn6a = 1.0e-10
c        if(nboundary .eq. 0)xn6a = s_n_s(ndata,6)

         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,6),
     &        xn6a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))

         xkti6 = 1.0e-10 *q
c        if(nboundary .eq. 0)xkti6 = s_t_s(ndata,6)

         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,6),
     &        xkti6(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))


c        if(nboundary .eq. 0)zeff = s_zeff(ndata)

         call flux_to_rz(nnodex,nnodey,s_zeff(1:s_nrho_n),
     &        zeff(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n),
     &        s_nrho_n, rho(1:nnodex,1:nnodey))


         if (nboundary .eq. 0) then
            if (iprofile .eq. 5)then
*           -------------------------
*           analytic scrape off layer
*           -------------------------

            do i = 1, nnodex
               do j = 1, nnodey
                  drho1 = rho(i,j) - rholim
                  if(rho(i,j) .gt. rholim) then

                     xna_slo(i,j) = 0.0e-10

                     xnea(i,j) = nmin +
     &                  (s_n_s(ndata, 0) - nmin) * exp(- drho1 / scrape)
                     xn1a(i,j) = nmin1 +
     &                  (s_n_s(ndata, 1) - nmin1) * exp(- drho1/ scrape)
                     xn2a(i,j) = nmin2 +
     &                  (s_n_s(ndata, 2) - nmin2) * exp(- drho1/ scrape)
                     xn3a(i,j) = nmin3 +
     &                  (s_n_s(ndata, 3) - nmin3) * exp(- drho1/ scrape)
                     xn4a(i,j) = nmin4 +
     &                  (s_n_s(ndata, 4) - nmin4) * exp(- drho1/ scrape)
                     xn5a(i,j) = nmin5 +
     &                  (s_n_s(ndata, 5) - nmin5) * exp(- drho1/ scrape)
                     xn6a(i,j) = nmin6 +
     &                  (s_n_s(ndata, 6) - nmin6) * exp(- drho1/ scrape)

c                     xnuomea(i,j) = xnuome * exp(+ drho1/ scrape)
c                    xnuomia(i,j) = xnuomi * exp(+ drho1/ scrape)

                     xnuomea(i,j) = xnuome
                     xnuomia(i,j) = xnuomi

c                     xkte(i,j)  = s_t_s(ndata,0)
c                     xkti(i,j)  = s_t_s(ndata,1)
c                     xkti2(i,j) = s_t_s(ndata,2)
c                     xkti3(i,j) = s_t_s(ndata,3)
c                     xkti4(i,j) = s_t_s(ndata,4)
c                     xkti5(i,j) = s_t_s(ndata,5)
c                     xkti6(i,j) = s_t_s(ndata,6)
                     xkte(i,j) = Temin0 +
     &                  (s_t_s(ndata, 0) -Temin0) * exp(- drho1/ scrape)
                     xkti(i,j) = Tmin1 +
     &                  (s_t_s(ndata, 1) - Tmin1) * exp(- drho1/ scrape)
                     xkti2(i,j) = Tmin2 +
     &                  (s_t_s(ndata, 2) - Tmin2) * exp(- drho1/ scrape)
                     xkti3(i,j) = Tmin3 +
     &                  (s_t_s(ndata, 3) - Tmin3) * exp(- drho1/ scrape)
                     xkti4(i,j) = Tmin4 +
     &                  (s_t_s(ndata, 4) - Tmin4) * exp(- drho1/ scrape)
                     xkti5(i,j) = Tmin5 +
     &                  (s_t_s(ndata, 5) - Tmin5) * exp(- drho1/ scrape)
                     xkti6(i,j) = Tmin6 +
     &                  (s_t_s(ndata, 6) - Tmin6) * exp(- drho1/ scrape)
                     zeff(i,j) = s_zeff(ndata)
                  end if
               end do
            end do
            end if
         end if

         if (nboundary .eq. 1) then
            do i = 1, nnodex
               do j = 1, nnodey
                  if(rho(i,j) .gt. rholim) then
                     xnea(i,j) = 1.0e-10
                     zeff(i,j) = 1.0e-10
                     xn3a(i,j) = 1.0e-10
                     xna_slo(i,j) = 0.0e-10
                     xkte(i,j) = 1.0e-10 * q
                     xkti(i,j) = 1.0e-10 * q
                     xkti3(i,j) = 1.0e-10 * q
                  end if
               end do
            end do
         end if

*        ---------------------------------
*        Minority ion (species 2) profiles:
*        ---------------------------------
*        for now leave the species 2 minority option open
*        if eta comes in as zero or small, use namelist values
*        else create a minority

         if(eta .ge. 1.0e-6) then
            do i = 1, nnodex
               do j = 1, nnodey
                  xn2a(i,j) = eta * xnea(i,j)

*                 ----------------
*                 Cold minority H:
*                 ----------------
*                  xkti2(i,j) = xkti(i,j)

*                 ----------------
*                 Hot minority H:
*                 ----------------
c                 xkti2(i,j) = xkti3(i,j)

               end do
            end do

            do i = 1, nnodex
               do j = 1, nnodey
                  xn1a(i, j) = (xnea(i, j) - z2 * xn2a(i, j)
     &                                  - z3 * xn3a(i, j)
     &                                  - z4 * xn4a(i, j)
     &                                  - z5 * xn5a(i, j)
     &                                  - z6 * xn6a(i, j)
     &                               - z_slo * xna_slo(i, j)   ) / z1

                  if (xn1a(i, j) .le. 0.0) xn1a(i, j) = 1.0e-10
               end do
            end do

         end if

*        --------------------------------------------------
*        Adjust ni1a for charge neutrality on all processors:
*        --------------------------------------------------
         do i = 1, nnodex
            do j = 1, nnodey
               xn1a(i, j) = (xnea(i, j) - z2 * xn2a(i, j)
     &                                  - z3 * xn3a(i, j)
     &                                  - z4 * xn4a(i, j)
     &                                  - z5 * xn5a(i, j)
     &                                  - z6 * xn6a(i, j)
     &                               - z_slo * xna_slo(i, j)   ) / z1

               if (xn1a(i, j) .le. 0.0) xn1a(i, j) = 1.0e-10
            end do
         end do



         xn00 = xnea(i0, j0)
         xn10 = xn1a(i0, j0)
         xn20 = xn2a(i0, j0)
         xn30 = xn3a(i0, j0)
         xn_slo0 = xna_slo(i0, j0)
         eta1 = xn10 / xn00
         eta2 = xn20 / xn00
         eta3 = xn30 / xn00
         eta_slo = xn_slo0 / xn00

         if(myid .eq. 0)then
            write(6, 6834)eta1
            write(6, 6835)eta2
            write(6, 6836)eta3

            write(15, 6834)eta1
            write(15, 6835)eta2
            write(15, 6836)eta3
         end if


         call fluxavg3(xnea, xnavg, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(zeff, zeffavg1, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xn1a, xn1avg, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xn2a, xn2avg, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xn3a, xn3avg, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xna_slo, xna_sloavg, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xkte, xkteavg, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xkti, xktiavg, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xkti2, xkti2avg, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xkti3, xkti3avg, rho, nxmx, nymx, nrhomax,
     &      nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)


         xn00 = xnea(i0, j0)
         xn10 = xn1a(i0, j0)
         xn20 = xn2a(i0, j0)
         xn30 = xn3a(i0, j0)
         xn40 = xn4a(i0, j0)
         xn50 = xn5a(i0, j0)
         xn60 = xn6a(i0, j0)
         xn_slo0 = xna_slo(i0, j0)

         eta1 = xn10 / xn00
         eta2 = xn20 / xn00
         eta3 = xn30 / xn00
         eta4 = xn40 / xn00
         eta5 = xn50 / xn00
         eta6 = xn60 / xn00
         eta_slo = xn_slo0 / xn00

      end if

*     --------------------------
*     End of iprofile = 5 option
*     --------------------------



c     write (6, *) "b0 = ", b0
c     write (6, *) "xmu0 = ", xmu0
c     write (6, *) "xn1 = ", xn1
c     write (6, *) "xmi1 = ", xmi1

      valfven = sqrt(b0**2/(xmu0 * xn1 * xmi1))
      kalfven = omgrf / valfven

      kperprhoi1 = 0.0
      kperprhoi2 = 0.0
      kperprhoi3 = 0.0
      kperprhoi4 = 0.0
      kperprhoi5 = 0.0

      if (eta1 .ne. 0.0) kperprhoi1 = kalfven * rhoi10
      if (eta2 .ne. 0.0) kperprhoi2 = kalfven * rhoi20
      if (eta3 .ne. 0.0) kperprhoi3 = kalfven * rhoi30
      if (eta4 .ne. 0.0) kperprhoi4 = kalfven * rhoi40
      if (eta5 .ne. 0.0) kperprhoi5 = kalfven * rhoi50


*     -----------
*     Mask array:
*     -----------
      do i = 1, nnodex
         do j = 1, nnodey

            mask(i,j) = 1
            if (psi(i,j) .gt. psimask) mask(i,j) = 0

         end do
      end do

*     ----------------------------------------------
*     Calculate the differential volume on half mesh:
*     ----------------------------------------------

      dvol  = 0.0
      darea = 1.0e-05

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1
            if(n .le. nnoderho)then
               dvol(n)  =  dvol(n) + dx * dy * 2.0 * pi * capr(i)
               darea(n) = darea(n) + dx * dy * capr(i) / rt
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




*     ------------------------
*     Do flux surface averages:
*     ------------------------

      call fluxavg(xnea, xnavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn1a, xn1avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn2a, xn2avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn3a, xn3avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn4a, xn4avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn5a, xn5avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn6a, xn6avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xna_slo, xna_sloavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(xkte, xkteavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti, xktiavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti2, xkti2avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti3, xkti3avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti4, xkti4avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti5, xkti5avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti6, xkti6avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


c      call polavg(dldb_tot12, dldbavg, rho, nxmx, nymx, nrhomax,
c     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvolp, fvol)


      if(myid .eq. 0) then
         write(6, *)
         write(6, *)"        i    taup        tauii        nustar",
     &                "     epsa         xnea        xkti"
         write(6, *)

         write(15,*)
         write(15,*)"        i    taup        tauii        nustar",
     &                "     epsa         xnea        xkti"
         write(15,*)
      end if

      do i = 1, nnodex
         do j = 1, nnodey

            rhome  = xnea(i, j) * xme
            rhomi1 = xn1a(i, j) * xmi1
            rhomi2 = xn2a(i, j) * xmi2
            rhomi3 = xn3a(i, j) * xmi3
            rhomi4 = xn4a(i, j) * xmi4
            rhomi5 = xn5a(i, j) * xmi5
            rhomi6 = xn6a(i, j) * xmi6

            rhom1(i,j) = xn1a(i, j) * xmi1


            presse  = xnea(i, j) * xkte(i, j)
            pressi1 = xn1a(i, j) * xkti(i, j)
            pressi2 = xn2a(i, j) * xkti2(i, j)
            pressi3 = xn3a(i, j) * xkti3(i, j)
            pressi4 = xn4a(i, j) * xkti4(i, j)
            pressi5 = xn5a(i, j) * xkti5(i, j)
            pressi6 = xn6a(i, j) * xkti6(i, j)

c           pressi(i, j) = presse + pressi1 pressi2 + pressi3
c     &                             + pressi4 pressi5 + pressi6
            pressi(i, j) = pressi1

            xnuee = 2.859e-27 * xnea(i,j) * xlnlam /
     &                      (sqrt(xme) * (xkte(i,j)/q)**1.5)

            xnuii = 2.859e-27 * xn1a(i,j) * xlnlam /
     &                      (sqrt(xmi1) * (xkti(i,j)/q)**1.5)

            if(xnuii .ne. 0.0)then

               vthi0 = sqrt(2.0 * xkti(i,j) / xmi1)
               vthe0 = sqrt(2.0 * xkte(i,j) / xme)

               tauii = 1.0 / xnuii
               omgti =  vthi0 / (rt * qhat)
               omgte =  vthe0 / (rt * qhat)
               xnu7omg = xnuii / omgti

*              -----------------------------
*              Assume circular flux surfaces
*              -----------------------------

               a = sqrt((x(i)-r0)**2 + (y(j)-z0)**2)

               rmaxa = rt + a
               rmina = rt - a

               epsa = (rmaxa - rmina) / (rmaxa + rmina)
               if( epsa .eq. 0.0) epsa = 1.0e-08

               xnustar = xnuii / (omgti * epsa**1.5)
               nu_star(i, j) = xnustar

*              -----------------------------------------------
*              bananna (collisionless) regime:   xnuii/omgti << eps**1.5
*              -----------------------------------------------
*              taup_b = 1.29 * epsa**1.5 * tauii

*              -------------------------------------------------
*              plateau regime:    eps**1.5  << xnuii/omgti << 1
*              -------------------------------------------------
*              taup_pl = 2.0 / sqrt(pi) / omgti

*              --------------------------------------------
*              Pfirsch - Schluter regime:  xnuii/omgti >> 1
*              --------------------------------------------
*              taup_ps = 1.0 / ( omgti * tauii) * (1.0 / omgti)


*              ----------------------------------------
*              Mult-regime model:  arbitray xnuii/omgti
*              ----------------------------------------
               taup_inv = 1.5 * omgti * xnustar / (1.0 + xnustar)
     &                             * 1.0 / (1.0 + xnustar * epsa**1.5)
               taup = 1.0 / taup_inv

               if(myid .eq. 0) then
                  if(j .eq. 16) then
                     write(6, 1312) i, taup, tauii, xnustar, epsa,
     &                  xnea(i,j), xkti(i,j)
                     write(15,1312) i, taup, tauii, xnustar, epsa,
     &                  xnea(i,j), xkti(i,j)
                  end if
               end if

*              ---------------------
*              poloidal viscosities:
*              ---------------------
               muhat(i, j) = epsa**2 / taup

            end if





            signb = b0 / abs(b0)

            omgce(i, j)  = qe  * bmod(i, j) / xme  * signb
            omgci1(i, j) = qi1 * bmod(i, j) / xmi1 * signb
            omgci2(i, j) = qi2 * bmod(i, j) / xmi2 * signb
            omgci3(i, j) = qi3 * bmod(i, j) / xmi3 * signb
            omgci4(i, j) = qi4 * bmod(i, j) / xmi4 * signb
            omgci5(i, j) = qi5 * bmod(i, j) / xmi5 * signb
            omgci6(i, j) = qi6 * bmod(i, j) / xmi6 * signb
            omgci_slo(i, j) = qi_slo * bmod(i, j) / xmi_slo * signb

            omgpe2(i, j) = xnea(i ,j) * qe**2  / (eps0 * xme)
            omgp12(i, j) = xn1a(i, j) * qi1**2 / (eps0 * xmi1)
            omgp22(i, j) = xn2a(i, j) * qi2**2 / (eps0 * xmi2)
            omgp32(i, j) = xn3a(i, j) * qi3**2 / (eps0 * xmi3)
            omgp42(i, j) = xn4a(i, j) * qi4**2 / (eps0 * xmi4)
            omgp52(i, j) = xn5a(i, j) * qi5**2 / (eps0 * xmi5)
            omgp62(i, j) = xn6a(i, j) * qi6**2 / (eps0 * xmi6)
            omgp2_slo(i, j) = xna_slo(i, j) * qi_slo**2
     &                                         / (eps0 * xmi_slo)

            omglh2i = 1.0 / (omgci1(i,j)**2 + omgp12(i,j))
     &           + 1.0 / abs(omgci1(i,j) * omgce(i,j))


            omglh(i,j) =  sqrt(1.0 / omglh2i)


            xb(i,j) = -zi / omgrf / eps0 * xjx(i,j)
            xc(i,j) = -zi / omgrf / eps0 * xjy(i,j)
            xd(i,j) = -zi / omgrf / eps0 * xjz(i,j)

            capk_perp = 1.0 - omgpe2(i,j) / (omgrf**2 - omgce(i,j)**2)
     &                      - omgp12(i,j) / (omgrf**2 - omgci1(i,j)**2)
     &                      - omgp22(i,j) / (omgrf**2 - omgci2(i,j)**2)
     &                      - omgp32(i,j) / (omgrf**2 - omgci3(i,j)**2)
     &                      - omgp42(i,j) / (omgrf**2 - omgci4(i,j)**2)
     &                      - omgp52(i,j) / (omgrf**2 - omgci5(i,j)**2)
     &                      - omgp62(i,j) / (omgrf**2 - omgci6(i,j)**2)

            capk_x =(omgce(i,j) * omgpe2(i,j)/(omgrf**2-omgce(i,j)**2)
     &            + omgci1(i,j) * omgp12(i,j)/(omgrf**2-omgci1(i,j)**2)
     &            + omgci2(i,j) * omgp22(i,j)/(omgrf**2-omgci2(i,j)**2)
     &            + omgci3(i,j) * omgp32(i,j)/(omgrf**2-omgci3(i,j)**2)
     &            + omgci4(i,j) * omgp42(i,j)/(omgrf**2-omgci4(i,j)**2)
     &            + omgci5(i,j) * omgp52(i,j)/(omgrf**2-omgci5(i,j)**2)
     &            + omgci6(i,j) * omgp62(i,j)/(omgrf**2-omgci6(i,j)**2))
     &             / omgrf

           left  = capk_perp - capk_x
           right = capk_perp + capk_x

            xnphi = bzn(i,j) * xkphi(i) / xk0
c            acold(i,j) = xnphi**2 - capk_perp

            acold(i,j) = capk_perp - xnphi**2
            bcold(i,j) = left - xnphi**2
            ccold(i,j) = right - xnphi**2

            xkperp_cold2(i,j) = bcold(i,j) * ccold(i,j) / acold(i,j)
            xkperp_cold(i,j) = csqrt(xkperp_cold2(i,j))
            if (j .eq. jmid) kimag_cold(i) = aimag(xkperp_cold(i,j))

         end do
      end do

      call dgrate(xprime, kimag_cold, 1, nnodex, ans, nxmx)
      spa_cold = 1.0 - exp(-2.0 * ans)

      if (myid .eq. 0) then

         write(15, *) "spa_cold = ", spa_cold
         write(6,  *) "spa_cold = ", spa_cold

         write (15, *) 'valfven = ', valfven, 'm/s'
         write (6,  *) 'valfven = ', valfven, 'm/s'

         write (15, *) 'kalfven = ', kalfven, 'm-1'
         write (6,  *) 'kalfven = ', kalfven, 'm-1'

         write (15, *) 'kperprhoi1 = ', kperprhoi1
         write (15, *) 'kperprhoi2 = ', kperprhoi2
         write (15, *) 'kperprhoi3 = ', kperprhoi3
         write (15, *) 'kperprhoi4 = ', kperprhoi4
         write (15, *) 'kperprhoi5 = ', kperprhoi5

         write (6, *) 'kperprhoi1 = ', kperprhoi1
         write (6, *) 'kperprhoi2 = ', kperprhoi2
         write (6, *) 'kperprhoi3 = ', kperprhoi3
         write (6, *) 'kperprhoi4 = ', kperprhoi4
         write (6, *) 'kperprhoi5 = ', kperprhoi5

         write (15, *) 'eta1 = ', eta1
         write (15, *) 'eta2 = ', eta2
         write (15, *) 'eta3 = ', eta3
         write (15, *) 'eta4 = ', eta4
         write (15, *) 'eta5 = ', eta5
         write (15, *) 'eta6 = ', eta6

         write (6, *) 'eta1 = ', eta1
         write (6, *) 'eta2 = ', eta2
         write (6, *) 'eta3 = ', eta3
         write (6, *) 'eta4 = ', eta4
         write (6, *) 'eta5 = ', eta5
         write (6, *) 'eta6 = ', eta6

         write(15, 920)
         write(15, 164)
         write(15, 920)

         write(6, 920)
         write(6, 164)
         write(6, 920)

      endif


      do i = 1, nnodex
         do j = 1, nnodey
            reomg1a(i,j) = omgrf / omgci1(i, j)
            reomg2a(i,j) = omgrf / omgci2(i, j)
            reomg3a(i,j) = omgrf / omgci3(i, j)
            reomglha(i,j) = omgrf / omglh(i, j)
         end do
      end do

      j = jequat
      do i = 1, nnodex
         reomg1 = omgrf / omgci1(i, j)
         reomg2 = omgrf / omgci2(i, j)
         reomg3 = omgrf / omgci3(i, j)
         reomg4 = omgrf / omgci4(i, j)
         reomg5 = omgrf / omgci5(i, j)
         reomg6 = omgrf / omgci6(i, j)

         reomgpi = sqrt(omgp12(i,j))
         reomglh = reomglha(i,j)

         teev = xkte(i,j) / q
         re_acold = real(acold(i,j))

         if (myid .eq. 0) then
c            write(6, 2163)i, capr(i), x(i), btau(i, j),
c     &         bmod(i,j), reomg1, reomg2, reomg3, reomg4, reomg5,
c     &         xnea(i,j), re_acold
c            write(15, 2163)i, capr(i), x(i), btau(i, j),
c     &         bmod(i,j), reomg1, reomg2, reomg3, reomg4, reomg5,
c     &         xnea(i,j), re_acold

            write(6, 2163)i, capr(i), x(i), btau(i, j),
     &         bmod(i,j), reomg1, reomg2, reomg3, reomgpi, reomglh,
     &         xnea(i,j), re_acold
            write(15, 2163)i, capr(i), x(i), btau(i, j),
     &         bmod(i,j), reomg1, reomg2, reomg3, reomgpi, reomglh,
     &         xnea(i,j), re_acold
         endif
      end do

      if (myid .eq. 0) then
         write(15, 920)
         write(15, 9164)
         write(15, 920)

         write(6, 920)
         write(6, 9164)
         write(6, 920)
      endif

      i = icenter
      do j = 1, nnodey
         reomg1 = omgrf / omgci1(i, j)
         reomg2 = omgrf / omgci2(i, j)
         reomg3 = omgrf / omgci3(i, j)
         reomg4 = omgrf / omgci4(i, j)
         reomg5 = omgrf / omgci5(i, j)
         reomg6 = omgrf / omgci6(i, j)

         reomgpi = sqrt(omgp12(i,j))
         reomglh = reomglha(i,j)

         teev = xkte(i,j) / q

         if (myid .eq. 0) then
c            write(6, 2163)j, capr(i), y(j), btau(i, j),
c     &         bmod(i,j), reomg1, reomg2, reomg3, reomg4, reomg5,
c     &         xnea(i,j), re_acold
c            write(15, 2163)j, capr(i), y(j), btau(i, j),
c     &         bmod(i,j), reomg1, reomg2, reomg3, reomg4, reomg5,
c     &         xnea(i,j), re_acold

            write(6, 2163)j, capr(i), y(j), btau(i, j),
     &         bmod(i,j), reomg1, reomg2, reomg3, reomgpi, reomglh,
     &         xnea(i,j), re_acold
            write(15, 2163)j, capr(i), y(j), btau(i, j),
     &         bmod(i,j), reomg1, reomg2, reomg3, reomgpi, reomglh,
     &         xnea(i,j), re_acold
         endif
      end do



      if (myid .eq. 0) then
         write(15, 920)
         write(15, 8164)
         write(15, 920)

         write(6, 920)
         write(6, 8164)
         write(6, 920)
      endif

 8164 format(" ",1x," i  ",  3x," R(x) ",
     &           6x,"  x   ",6x," bmod ",6x,"   n1 ",6x,"  T1  ",
     &           6x,"  n2  ",6x,"  T2  ",6x,"   n3 ",6x,"  T3  ",
     &           6x,"  ne  ",6x,"  n4  ")

      j = jequat
      do i = 1, nnodex
         t1ev = xkti(i,j)  / q
         t2ev = xkti2(i,j) / q
         t3ev = xkti3(i,j) / q
         t4ev = xkti4(i,j) / q
         t5ev = xkti5(i,j) / q
         t6ev = xkti6(i,j) / q

         if (myid .eq. 0) then
            write(6, 2163) i, capr(i), x(i), bmod(i,j),
     &           xn1a(i,j), t1ev, xn2a(i,j), t2ev, xn3a(i,j), t3ev,
     &           xnea(i,j), xn4a(i,j)
            write(15, 2163)i, capr(i), x(i), bmod(i,j),
     &           xn1a(i,j), t1ev, xn2a(i,j), t2ev, xn3a(i,j), t3ev,
     &           xnea(i,j), xn4a(i,j)
         endif
      end do


      if(myid .eq. 0)write(6, *) "      n       xkxsav"
      do n = nkx1, nkx2
         xkxsav(n) = 2.0 * pi * n / xmax
         if(myid .eq. 0)write(6, *) n,  xkxsav(n)
      end do



      if(myid .eq. 0)write(6, *) "      m       xkysav"
      do m = nky1, nky2
         xkysav(m) = 2.0 * pi * m / ymax
         if(myid .eq. 0)write(6, *) m,  xkysav(m)
      end do

      mdisp = xkydisp * ymax / (2.0 * pi)

      n_ono = xkx_ono * xmax / (2.0 * pi)
      m_ono = xky_ono * ymax / (2.0 * pi)




      kperp_max_actual = sqrt(xkxsav(nkx2)**2 + xkysav(nky2)**2)

      kperp_max = kperp_max_actual * 2.0


      xk_cutoff = kperp_max_actual * xkperp_cutoff

      if(myid .eq. 0) then
         write (6, *) "kperp_max = ", kperp_max
         write(6, *) "xkydisp = ", xkydisp, "mdisp = ",mdisp
         write(6, *) "xkx_ono = ", xkx_ono, "xky_ono = ", xky_ono
         write (6,  *) "xk_cutoff = ", xk_cutoff
         write (15, *) "xk_cutoff = ", xk_cutoff
      end if


c     -----------------------------
c     Initialize Hilbert Z function
c     -----------------------------
c      call zfun_hil_init(ieer)

c     -------------------------------------
c     setup Smithe's Z function table when nzfun = 3
c     -------------------------------------
      if (nzfun .eq. 3) then

         iflag = 0

         call ztable_setup(myid, iflag)

         call blacs_barrier(icontxt, 'All')

         if(iflag .eq. 1)call ztable_setup(myid, nproc)

      end if


c     ---------------------------
c     precompute xx(n,i), yy(m,j)
c     ---------------------------

      do i = 1, nnodex
         do n = nkx1, nkx2
            xx(n, i) = exp(zi * xkxsav(n) * xprime(i))
            xx_inv(n,i) = 1.0/xx(n,i)
         end do
      end do

      do j = 1, nnodey
         do m = nky1, nky2
            yy(m, j) = exp(zi * xkysav(m) * yprime(j))
            yy_inv(m,j) = 1.0/yy(m,j)
         end do
      end do



*     ----------------------------------------
*     Take numerical derivatives for igeom = 5
*     ----------------------------------------
      if (igeom .eq. 5) then
      do i = 1, nnodex
         do j = 1, nnodey

            call deriv_x(bmod, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dbdx, d2bdx2)
            call deriv_y(bmod, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dbdy, d2bdy2)


*           -----
*           Exact:
*           -----
c            gradprlb(i,j) = abs(uzx(i,j) * dbdx + uzy(i,j) * dbdy)

*           -----------------------
*           Brambilla approximation:
*           -----------------------
            sinth = y(j) / sqrt(x(i)**2 + y(j)**2)
            gradprlb(i,j) = bmod(i,j) / capr(i) * abs(btau(i,j) * sinth)

*           -------------------------------
*           Constant gradient scale length:
*           -------------------------------
c            gradprlb(i,j) = bmod(i,j) / capr(i)


            if (nzfun .eq. 0)gradprlb(i,j) = 1.0e-10

c            if (myid .eq. 0 .and. i .eq. 34) then
c                  write(6,*) "j = ", j,
c     &                       "gradprlb(34,j) = ", gradprlb(34,j)
c            end if


            call deriv_x(rho, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         drhodx, drhodxx)
            call deriv_y(rho, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         drhody, drhodyy)

            gradrho = sqrt(drhodx**2 + drhody**2)
            if (gradrho .eq. 0.0) gradrho = 1.0e-08

            rhohatx(i,j) = drhodx / gradrho
            rhohaty(i,j) = drhody / gradrho


            call deriv_x(uxx, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dxuxx(i,j), dxxuxx(i,j))
            call deriv_x(uxy, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dxuxy(i,j), dxxuxy(i,j))
            call deriv_x(uxz, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dxuxz(i,j), dxxuxz(i,j))

            call deriv_x(uyx, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dxuyx(i,j), dxxuyx(i,j))
            call deriv_x(uyy, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dxuyy(i,j), dxxuyy(i,j))
            call deriv_x(uyz, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dxuyz(i,j), dxxuyz(i,j))

            call deriv_x(uzx, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dxuzx(i,j), dxxuzx(i,j))
            call deriv_x(uzy, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dxuzy(i,j), dxxuzy(i,j))
            call deriv_x(uzz, nxmx, nymx, i, j, nnodex, nnodey, dx,
     &         dxuzz(i,j), dxxuzz(i,j))



            call deriv_y(uxx, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dyuxx(i,j), dyyuxx(i,j))
            call deriv_y(uxy, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dyuxy(i,j), dyyuxy(i,j))
            call deriv_y(uxz, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dyuxz(i,j), dyyuxz(i,j))

            call deriv_y(uyx, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dyuyx(i,j), dyyuyx(i,j))
            call deriv_y(uyy, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dyuyy(i,j), dyyuyy(i,j))
            call deriv_y(uyz, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dyuyz(i,j), dyyuyz(i,j))

            call deriv_y(uzx, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dyuzx(i,j), dyyuzx(i,j))
            call deriv_y(uzy, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dyuzy(i,j), dyyuzy(i,j))
            call deriv_y(uzz, nxmx, nymx, i, j, nnodex, nnodey, dy,
     &         dyuzz(i,j), dyyuzz(i,j))



            call deriv_xy(uxx, nxmx, nymx, i, j, nnodex, nnodey, dx, dy,
     &         dxyuxx(i,j))
            call deriv_xy(uxy, nxmx, nymx, i, j, nnodex, nnodey, dx, dy,
     &         dxyuxy(i,j))
            call deriv_xy(uxz, nxmx, nymx, i, j, nnodex, nnodey, dx, dy,
     &         dxyuxz(i,j))

            call deriv_xy(uyx, nxmx, nymx, i, j, nnodex, nnodey, dx, dy,
     &         dxyuyx(i,j))
            call deriv_xy(uyy, nxmx, nymx, i, j, nnodex, nnodey, dx, dy,
     &         dxyuyy(i,j))
            call deriv_xy(uyz, nxmx, nymx, i, j, nnodex, nnodey, dx, dy,
     &         dxyuyz(i,j))

            call deriv_xy(uzx, nxmx, nymx, i, j, nnodex, nnodey, dx, dy,
     &         dxyuzx(i,j))
            call deriv_xy(uzy, nxmx, nymx, i, j, nnodex, nnodey, dx, dy,
     &         dxyuzy(i,j))
            call deriv_xy(uzz, nxmx, nymx, i, j, nnodex, nnodey, dx, dy,
     &         dxyuzz(i,j))

c           if (myid .eq. 0) then
c              dbxdx = dxuzx(i,j)
c              dbxdy = dyuzx(i,j)
c              dbydx = dxuzy(i,j)
c              dbydy = dyuzy(i,j)
c              dbphidx = dxuzz(i,j)
c              dbphidy = dyuzz(i,j)
c               write(65, 1314) i,j, bxn(i,j), byn(i,j), bzn(i,j),
c     &                              dbxdx, dbydx, dbphidx,
c     &                              dbxdy, dbydy, dbphidy
c           end if


         end do
      end do



      call smooth2d_(nsmooth, gradprlb, nxmx, nymx, nnodex, nnodey,work)
      call smooth2d_(nsmooth, rhohatx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, rhohaty, nxmx, nymx, nnodex, nnodey, work)

c      if (myid .eq. 0) then
c         do j = 1, nnodey
c            write(6,*) "j = ", j, "gradprlb(34,j) = ", gradprlb(34,j)
c        end do
c      end if


c      call smooth2d_(nsmooth, dbdx, nxmx, nymx, nnodex, nnodey, work)
c      call smooth2d_(nsmooth, dbdy, nxmx, nymx, nnodex, nnodey, work)
c      call smooth2d_(nsmooth, d2bdx2, nxmx, nymx, nnodex, nnodey, work)
c      call smooth2d_(nsmooth, d2bdy2, nxmx, nymx, nnodex, nnodey, work)


c      call smooth2d_(nsmooth, drhodx, nxmx, nymx, nnodex, nnodey, work)
c      call smooth2d_(nsmooth, drhodxx, nxmx, nymx, nnodex, nnodey, work)
c      call smooth2d_(nsmooth, drhody, nxmx, nymx, nnodex, nnodey, work)
c      call smooth2d_(nsmooth, drhodyy, nxmx, nymx, nnodex, nnodey, work)



      call smooth2d_(nsmooth, dxuxx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxuxy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxuxz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dxuyx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxuyy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxuyz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dxuzx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxuzy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxuzz, nxmx, nymx, nnodex, nnodey, work)



      call smooth2d_(nsmooth, dxxuxx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxxuxy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxxuxz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dxxuyx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxxuyy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxxuyz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dxxuzx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxxuzy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxxuzz, nxmx, nymx, nnodex, nnodey, work)


      call smooth2d_(nsmooth, dyyuxx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyyuxy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyyuxz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dyyuyx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyyuyy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyyuyz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dyyuzx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyyuzy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyyuzz, nxmx, nymx, nnodex, nnodey, work)


      call smooth2d_(nsmooth, dyuxx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyuxy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyuxz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dyuyx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyuyy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyuyz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dyuzx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyuzy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dyuzz, nxmx, nymx, nnodex, nnodey, work)


      call smooth2d_(nsmooth, dxyuxx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxyuxy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxyuxz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dxyuyx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxyuyy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxyuyz, nxmx, nymx, nnodex, nnodey, work)

      call smooth2d_(nsmooth, dxyuzx, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxyuzy, nxmx, nymx, nnodex, nnodey, work)
      call smooth2d_(nsmooth, dxyuzz, nxmx, nymx, nnodex, nnodey, work)

      do i = 1, nnodex
         do j = 1, nnodey

            dbxdx = dxuzx(i,j)
            dbxdy = dyuzx(i,j)
            dbydx = dxuzy(i,j)
            dbydy = dyuzy(i,j)
            dbphidx = dxuzz(i,j)
            dbphidy = dyuzz(i,j)

            kappa_hatx(i,j) = bxn(i,j) * dbxdx + byn(i,j) * dbxdy
     &                            - bzn(i,j)**2 / capr(i)

            kappa_haty(i,j) = bxn(i,j) * dbydx + byn(i,j) * dbydy

            kappa_hatz(i,j) = bxn(i,j) * dbphidx + byn(i,j) * dbphidy
     &                            + bzn(i,j) * bxn(i,j) / capr(i)

            mod_kappa = sqrt(kappa_hatx(i,j)**2
     &              + kappa_haty(i,j)**2 + kappa_hatz(i,j)**2)

            rc(i,j) = 1.0 / mod_kappa

            bhatx(i,j) = bxn(i,j)
            bhaty(i,j) = byn(i,j)
            bhatz(i,j) = bzn(i,j)

         end do
      end do


      end if

*     ----------------------------------------------
*     End of numerical devivatives for igeom .eq. 5
*     ---------------------------------------------

      time = second1(dummy) - t1
      tmin = time / 60.
      if(myid .eq. 0)then
         write(6, 386)tmin
         write(15, 386) tmin
c        call flush
      end if
  386 format('time to call eqdsk_setup = ',f9.3, " min")

*     -----------------------------------------------------------------
*     Read 2D lookup table for Lee Berry's new Z function (xnuomge = 0)
*     -----------------------------------------------------------------
      open(unit=195,file='Ztable',status='unknown',form='formatted')
      read (195, 3409)  nmax, mmax
      read (195, 3410) (zetai_table(n),n = 1,nmax)
      read (195, 3410) (dKdL_table(m), m = 1,mmax)
      read (195, 3410) ((z0_table1(n,m), n = 1,nmax), m = 1,mmax)
      read (195, 3410) ((z1_table1(n,m), n = 1,nmax), m = 1,mmax)
      read (195, 3410) ((z2_table1(n,m), n = 1,nmax), m = 1,mmax)
      close (195)


 3409 format (2i10)
 3410 format(1p,8e22.14)

      zetai_max_table = zetai_table(nmax)
      zetai_min_table = zetai_table(1)

      dKdL_max_table = dKdL_table(mmax)
      dKdL_min_table = dKdL_table(1)

      if (myid .eq. 0) then
         write(15, *) "dKdL_min_table = ", dKdL_min_table
         write(15, *) "dKdL_max_table = ", dKdL_max_table

         write(6,  *) "dKdL_min_table = ", dKdL_min_table
         write(6,  *) "dKdL_max_table = ", dKdL_max_table

         write(15, *) "zetai_min_table = ", zetai_min_table
         write(15, *) "zetai_max_table = ", zetai_max_table

         write(6, *) "zetai_min_table = ", zetai_min_table
         write(6, *) "zetai_max_table = ", zetai_max_table
      end if


      if (read_solution .eqv. .true.) then
         t1 = second1(dummy)
!        ---------------------------------------------
!        read solution in Fourier space to all processors
!        ---------------------------------------------
         read(54, 309) nkx1, nkx2, nky1, nky2
         read(54, 8310)((ealphak(n, m), n = nkx1, nkx2), m = nky1,nky2)
         read(54, 8310)((ebetak(n, m), n = nkx1, nkx2), m = nky1, nky2)
         read(54, 8310)((ebk(n, m), n = nkx1, nkx2), m = nky1, nky2)
         read(54, 8310)(xkxsav(n), n = nkx1, nkx2)
         read(54, 8310)(xkysav(m), m = nky1, nky2)

         time = second1(dummy) - t1
         tmin = time / 60.
         if(myid .eq. 0)then
            write(6, 6386)tmin
            write(15, 6386) tmin
c           call flush
         end if
 6386    format('time to read solution = ',f9.3, " min")

         go to 6000
      end if





      t1=second1(dummy)

*     ------------------------------------------
*     define scalapack array descriptor for amat
*     ------------------------------------------
      mb = 3*33
      nb = mb

      nrow = 3 * nnodex * nnodey
      ncol = nrow
      norder = nrow

      nrow_local = max(1,numroc( nrow,mb,myrow,0,nprow ))
      ncol_local = max(1,numroc( ncol,nb,mycol,0,npcol ))
      lld = max(1,nrow_local)

      p_amat_dim = nrow_local * ncol_local + 1
      p_amat_size = p_amat_dim * 16 / 1.0e+06

      if (isolve .ne. 0) then


        if (myid.eq.0) then
           write (6, 6839) p_amat_dim, p_amat_size
           write (15,6839) p_amat_dim, p_amat_size
        end if
 6839   format('p_amat_dim = ', i10, ' words,   ',
     &       'p_amat_size = ', i5, ' MBytes')

*     ----------------------------------------
*     allocate storage for scalapack matrices
*     ----------------------------------------

      allocate(p_amat(p_amat_dim))
      p_amat(:) = cmplx(0.0,0.0)

      if(write_for_azzam)then
        allocate(p_amat_copy(p_amat_dim))
        p_amat_copy(:) = cmplx(0.0,0.0)
      end if

      p_brhs_dim = nrow
      p_ipiv_dim = nrow
      allocate( p_brhs(p_brhs_dim), p_ipiv(p_ipiv_dim) )
      if(write_for_azzam)then
        allocate( p_brhs_copy(p_brhs_dim) )
      end if
      p_brhs(:) = cmplx(0.0,0.0)
      p_ipiv(:) = 0


      call descinit( desc_amat, nrow,ncol, mb,nb, 0,0, icontxt,
     &         lld, info )
      if (info.ne.0) then
        write(6,*) '** descinit for amat returns info = ', info
        write(6,*) 'lld = ', lld
        write(15,*) '** descinit for amat returns info = ', info
        write(15,*) 'lld = ', lld
        stop '** error ** '
      endif

      do irnc=1,nrow_local*ncol_local
        p_amat(irnc) = cmplx(0., 0.)
      enddo

      if (lld .gt. p_brhs_dim) then
        write(6,*) 'increase p_brhs_dim to ', lld+1
        write(15,*) 'increase p_brhs_dim to ', lld+1
        stop '** error ** '
      endif


      call descinit( desc_brhs, nrow, 1, mb,nb, 0,0, icontxt,
     &         lld, info )
      if (info.ne.0) then
        write(6,*) '** descinit for brhs returns info ', info
        write(15,*) '** descinit for brhs returns info ', info
        stop '** error ** '
      endif

*     --------------------------------------------------
*     Calculate minimum and maximum of zetai and dKdL
*     --------------------------------------------------
      zetai_max = 0.0
      zetai_min = 0.0

      dKdL_max = 0.0
      dKdL_min = 0.0

      xjx_max = 0.0
      xjy_max = 0.0
      xjz_max = 0.0

      xjx_min = 0.0
      xjy_min = 0.0
      xjz_min = 0.0

      do i = 1, nnodex
         do j = 1, nnodey

            if(xjx(i,j) .gt. xjx_max) xjx_max = xjx(i,j)
            if(xjy(i,j) .gt. xjy_max) xjy_max = xjy(i,j)
            if(xjz(i,j) .gt. xjz_max) xjz_max = xjz(i,j)

            if(xjx(i,j) .lt. xjx_min) xjx_min = xjx(i,j)
            if(xjy(i,j) .lt. xjy_min) xjy_min = xjy(i,j)
            if(xjz(i,j) .lt. xjz_min) xjz_min = xjz(i,j)

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     &                               .or. nboundary .eq. 0)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     xkalp = uxx(i, j) * xkxsav(n)
     &                     + uxy(i, j) * xkysav(m)
     &                     + uxz(i, j) * xkphi(i)
                     xkbet = uyx(i, j) * xkxsav(n)
     &                     + uyy(i, j) * xkysav(m)
     &                     + uyz(i, j) * xkphi(i)
                     xkprl = uzx(i, j) * xkxsav(n)
     &                     + uzy(i, j) * xkysav(m)
     &                     + uzz(i, j) * xkphi(i)


                     xkperp = sqrt(xkalp**2 + xkbet**2)
                     xkprl_orig = xkprl

                     alphae = sqrt(2.0 * xkte(i,j) / xme)
                     zetai = xkprl_orig * alphae / omgrf

                     alpha1 = sqrt(2.0 * xkti(i,j) / xmi1)

                     if(zetai .gt. zetai_max) zetai_max = zetai
                     if(zetai .lt. zetai_min) zetai_min = zetai

                     if (xkprl  .eq. 0.0) xkprl  = 1.0e-08


                     if (xkperp .eq. 0.0) xkperp = 1.0e-08

                     dbxdx = dxuzx(i,j)
                     dbxdy = dyuzx(i,j)

                     dbydx = dxuzy(i,j)
                     dbydy = dyuzy(i,j)

                     dbphidx = dxuzz(i,j)
                     dbphidy = dyuzz(i,j)

                     d2bxdx2 = dxxuzx(i,j)
                     d2bydx2 = dxxuzy(i,j)
                     d2bzdx2 = dxxuzz(i,j)

                     d2bxdxy = dxyuzx(i,j)
                     d2bxdy2 = dyyuzx(i,j)
                     d2bzdxy = dxyuzz(i,j)

                     d2bphiRdx2 = -2. / capr(i) * drdx(i) * dbphidx
     &                  + d2bzdx2 + bzn(i,j) / capr(i)**2 * drdx(i)**2
     &                  - bzn(i,j) / capr(i) * d2rdx2(i)

                     d2bphiRdxy = - 1. / capr(i) * drdx(i) * dbphidy
     &                           + d2bzdx2

                     dbphiRdx = dbphidx - bzn(i,j) / capr(i) * drdx(i)

                     bhatgradbx =  bxn(i,j) * dbxdx
     &                           + byn(i,j) * dbxdy
                     bhatgradby =  bxn(i,j) * dbydx
     &                           + byn(i,j) * dbydy
                     bhatgradbz =  bxn(i,j) * (dbphidx
     &                           - bzn(i,j) / capr(i) * drdx(i))
     &                           + byn(i,j) * dbphidy

                     bhatgrad2bx =  bxn(i,j) * dbxdx**2
     &                            + bxn(i,j)**2 * d2bxdx2
     &                            + dbxdy * (bxn(i,j) * dbydx
     &                                     + byn(i,j) * dbxdx
     &                                     + byn(i,j) * dbydy)
     &                            + 2.0 * bxn(i,j) * byn(i,j) * d2bxdxy
     &                            + byn(i,j)**2 * d2bxdy2

                     bhatgrad2by =  byn(i,j) * dbydy**2
     &                            + bxn(i,j)**2 * d2bydx2
     &                            + dbydx * (bxn(i,j) * dbxdx
     &                                     + bxn(i,j) * dbydy
     &                                     + byn(i,j) * dbxdy)
     &                            + 2.0 * bxn(i,j) * byn(i,j) * d2bydxy
     &                            + byn(i,j)**2 * d2bydy2

                     bhatgrad2bz =  bxn(i,j)**2 * d2bzdx2
     &                            + byn(i,j)**2 * d2bzdy2
     &                            + 2. * bxn(i,j) * byn(i,j) * d2bzdxy
     &                            + dbphiRdx * (bxn(i,j) * dbxdx
     &                                        + byn(i,j) * dbxdy)
     &                            + dbphidy * (bxn(i,j) * dbydx
     &                                       + byn(i,j) * dbydy)


                     dkprldl = xkxsav(n) * bhatgradbx
     &                       + xkysav(m) * bhatgradby
     &                       + xkphi(i)  * bhatgradbz

                     d2kprldl2 = xkxsav(n) * bhatgrad2bx
     &                         + xkysav(m) * bhatgrad2by
     &                         + xkphi(i)  * bhatgrad2bz

                     dKdL_giv = (alphae / omgrf)**2 * dkprldl



                     dKdL_giv_1 = (alpha1 / omgrf)**2 * dkprldl

*                    ------------------------------------
*                    Calculate max and min dK/dL
*                    --------------------------------- --
                     if (dKdL_giv .gt. dKdL_max) dKdL_max = dKdL_giv
                     if (dKdL_giv .lt. dKdL_min) dKdL_min = dKdL_giv
                  end do
               end do

            end if

         end do
      end do

      if (myid .eq. 0) then

         write(15, *) "xjx_max = ", xjx_max
         write(15, *) "xjy_max = ", xjy_max
         write(15, *) "xjz_max = ", xjz_max

         write(15, *) "xjx_min = ", xjx_min
         write(15, *) "xjy_min = ", xjy_min
         write(15, *) "xjz_min = ", xjz_min

         write(15, *) "dKdL_min = ", dKdL_min
         write(15, *) "dKdL_max = ", dKdL_max

         write(6,  *) "dKdL_min = ", dKdL_min
         write(6,  *) "dKdL_max = ", dKdL_max

         write(15, *) "zetai_min = ", zetai_min
         write(15, *) "zetai_max = ", zetai_max

         write(6, *) "zetai_min = ", zetai_min
         write(6, *) "zetai_max = ", zetai_max

      end if

      psiorb = psilim * 10.0

*     ----------------------------------------------------------------
*     Load x, y and z equations for spatial point (i,j) and mode number (n,m)
*     ----------------------------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey

            xnuome = xnuomea(i,j)
            xnuomi = xnuomia(i,j)


            irnc1 = (j-1) * 3 + (i-1) * nnodey * 3 + 1



            if (use_fourier_z2 .eqv. .true.) then

*           ---------------------------------------------------------
*           if nboundary = 1, just do orbits in interior plasma region
*           ---------------------------------------------------------
            if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     &                                    .or. nboundary .eq. 0) then
            if(psi(i,j) .lt. psiorb) then

*              -----------------------------
*              Define mesh along field line:
*              -----------------------------
c--            lb(nl):  0 to lbmax
c--            lbprime(nl): -lbmax / 2 to lbmax / 2

                lbmax = 2.0 * pi * rt * qsafety(i, j)


               dlb = lbmax / nnodelb

               do nl = 1, nnodelb
                  lbprime(nl) = (nl - 1) * dlb
                  lb(nl) = lbprime(nl) - lbmax / 2.0
               end do

               do nkl = nkl1, nkl2
                  xklsav(nkl) = 2.0 * pi * nkl / lbmax
                  if (nkl .eq. 0) xklsav(nkl) = 1.0e-05
               end do


c              if (myid .eq. 0) then

               sgn_vprl = +1.0
               call field_line_trace(sgn_vprl, i, j,
     &            nmodesx, nmodesy,
     &            rwleft, rwright, ytop, ybottom, xprime, yprime,
     &            rmaxis, zmaxis, b0, psio, psimag, psi_tor_max,
     &            bxn, byn, bzn, bmod, capr, y, psi_pol2d, rho_pol2d,
     &            qsafety, dx, dy,
     &            bmod_mid,  capr_bpol_mid2, capr_bpol_mid, rho_tor2d,
     &            i_psi_eq, dldb_tot12, dldbavg,
     &            n_prof_flux, rhomax, norb_dim, lbmax,
     &            len_xp, phin_xp, capr_xp, capz_xp, nlen_exitp)

               sgn_vprl = -1.0
               call field_line_trace(sgn_vprl, i, j,
     &            nmodesx, nmodesy,
     &            rwleft, rwright, ytop, ybottom, xprime, yprime,
     &            rmaxis, zmaxis, b0, psio, psimag, psi_tor_max,
     &            bxn, byn, bzn, bmod, capr, y, psi_pol2d, rho_pol2d,
     &            qsafety, dx, dy,
     &            bmod_mid,  capr_bpol_mid2, capr_bpol_mid, rho_tor2d,
     &            i_psi_eq, dldb_tot12, dldbavg,
     &            n_prof_flux, rhomax, norb_dim, lbmax,
     &            len_xm, phin_xm, capr_xm, capz_xm, nlen_exitm)


               do norbit = 1, nlen_exitm
                  len_x(norbit)  = len_xm(nlen_exitm - norbit + 1)
                  capr_x(norbit)  = capr_xm(nlen_exitm - norbit + 1)
                  capz_x(norbit)  = capz_xm(nlen_exitm - norbit + 1)
                  phin_x(norbit)  = phin_xm(nlen_exitm - norbit + 1)
               end do

               do norbit = nlen_exitm + 1, nlen_exitm + nlen_exitp
                  len_x(norbit) = len_xp(norbit - nlen_exitm)
                  capr_x(norbit) = capr_xp(norbit - nlen_exitm)
                  capz_x(norbit) = capz_xp(norbit - nlen_exitm)
                  phin_x(norbit) = phin_xp(norbit - nlen_exitm)
               end do

               nlen_exit = nlen_exitm + nlen_exitp


*              -----------------------------
*              Calculate spline coeficients
*              -----------------------------
               call spline (len_x, capr_x,
     &                       b_capr, c_capr, d_capr, nlen_exit)
               call spline (len_x, capz_x,
     &                       b_capz, c_capz, d_capz, nlen_exit)
               call spline (len_x, phin_x,
     &                       b_phin, c_phin, d_phin, nlen_exit)

               do nl = 1, nnodelb
                  capr_lb(nl) = ispline(lb(nl), len_x, capr_x,
     &                           b_capr, c_capr, d_capr, nlen_exit)
                  capz_lb(nl) = ispline(lb(nl), len_x, capz_x,
     &                           b_capz, c_capz, d_capz, nlen_exit)
                  phin_lb(nl) = ispline(lb(nl), len_x, phin_x,
     &                           b_phin, c_phin, d_phin, nlen_exit)
               end do


               ! dead block for debugging : JCW 2023 (who didnt code this block)
               go to 7070

c              if(myid .eq. 0) then

                  open(unit=116,file='trace',status='unknown',
     &                                                 form='formatted')
                  write(116, 309) nnodex, nnodey
                  write(116, 310) rhoplasm
                  write(116, 310) (x(idum), idum = 1, nnodex)
                  write(116, 310) (y(jdum), jdum = 1, nnodey)
                  write(116, 310) (capr(idum), idum = 1, nnodex)

                  write(116, 310) ((rho(idum, jdum), idum = 1,nnodex),
     &                 jdum = 1,nnodey)

                  write(116, 309) nnodelb
                  write(116, 310) (capr_lb(nl), nl = 1, nnodelb)
                  write(116, 310) (capz_lb(nl), nl = 1, nnodelb)

                  close (116)
                  call trace_plot
c               end if

               call blacs_barrier(icontxt, 'All')
               stop 3
 7070          continue

            end if    !  end if psi(i,j) .lt. psiorb

            end if    ! end if only do orbits in interior if nboundary = 1:
            end if    ! end if use_fourier_z2 .eq. .true.


            do n = nkx1, nkx2
               do m = nky1, nky2


                  icnc1 = (n - nkx1) * 3 * nmodesy + (m - nky1) * 3 + 1

                  call infog2l(irnc1, icnc1, desc_amat, nprow, npcol,
     &                myrow, mycol, lrindx, lcindx, rsrc1, csrc1)

                  ismine1 = (rsrc1 .eq. myrow) .and. (csrc1 .eq. mycol)


                  if (ismine1) then

                  cexpkxky = xx(n, i) * yy(m, j)

                  sigxx = 0.0
                  sigxy = 0.0
                  sigxz = 0.0
                  sigyx = 0.0
                  sigyy = 0.0
                  sigyz = 0.0
                  sigzx = 0.0
                  sigzy = 0.0
                  sigzz = 0.0

*                 -------------------------------------------------
*                 if nboundary = 1, only do sigma in interior plasma
*                 -------------------------------------------------
                  if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     &                                     .or. nboundary .eq. 0) then

                     if (isigma .eq. 1)then

                        xkalp = uxx(i,j) * xkxsav(n)
     &                        + uxy(i,j) * xkysav(m)
     &                        + uxz(i,j) * xkphi(i)
                        xkbet = uyx(i,j) * xkxsav(n)
     &                        + uyy(i,j) * xkysav(m)
     &                        + uyz(i,j) * xkphi(i)
                        xkprl = uzx(i,j) * xkxsav(n)
     &                        + uzy(i,j) * xkysav(m)
     &                        + uzz(i,j) * xkphi(i)

                        xkperp = sqrt(xkalp**2 + xkbet**2)

                        if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
                        if (xkperp .eq. 0.0) xkperp = 1.0e-08

                        dbxdx = dxuzx(i,j)
                        dbxdy = dyuzx(i,j)

                        dbydx = dxuzy(i,j)
                        dbydy = dyuzy(i,j)

                        dbphidx = dxuzz(i,j)
                        dbphidy = dyuzz(i,j)

                        d2bxdx2 = dxxuzx(i,j)
                        d2bydx2 = dxxuzy(i,j)
                        d2bzdx2 = dxxuzz(i,j)

                        d2bxdxy = dxyuzx(i,j)
                        d2bxdy2 = dyyuzx(i,j)
                        d2bzdxy = dxyuzz(i,j)



                        d2bphiRdx2 = -2. / capr(i) * drdx(i) * dbphidx
     &                  + d2bzdx2 + bzn(i,j) / capr(i)**2 * drdx(i)**2
     &                  - bzn(i,j) / capr(i) * d2rdx2(i)

                        d2bphiRdxy = - 1. / capr(i) * drdx(i) * dbphidy
     &                           + d2bzdx2

                        dbphiRdx = dbphidx - bzn(i,j) / capr(i) *drdx(i)


                        bhatgradbx =  bxn(i,j) * dbxdx
     &                              + byn(i,j) * dbxdy
                        bhatgradby =  bxn(i,j) * dbydx
     &                              + byn(i,j) * dbydy
                        bhatgradbz =  bxn(i,j) * (dbphidx
     &                              - bzn(i,j) / capr(i) * drdx(i))
     &                              + byn(i,j) * dbphidy

                        bhatgrad2bx =  bxn(i,j) * dbxdx**2
     &                               + bxn(i,j)**2 * d2bxdx2
     &                               + dbxdy * (bxn(i,j) * dbydx
     &                                     + byn(i,j) * dbxdx
     &                                     + byn(i,j) * dbydy)
     &                               + 2. * bxn(i,j) * byn(i,j) *d2bxdxy
     &                               + byn(i,j)**2 * d2bxdy2

                        bhatgrad2by =  byn(i,j) * dbydy**2
     &                               + bxn(i,j)**2 * d2bydx2
     &                               + dbydx * (bxn(i,j) * dbxdx
     &                                     + bxn(i,j) * dbydy
     &                                     + byn(i,j) * dbxdy)
     &                               + 2. * bxn(i,j) * byn(i,j) *d2bydxy
     &                               + byn(i,j)**2 * d2bydy2

                        bhatgrad2bz =  bxn(i,j)**2 * d2bzdx2
     &                               + byn(i,j)**2 * d2bzdy2
     &                               + 2. * bxn(i,j) * byn(i,j) *d2bzdxy
     &                               + dbphiRdx * (bxn(i,j) * dbxdx
     &                                        + byn(i,j) * dbxdy)
     &                               + dbphidy * (bxn(i,j) * dbydx
     &                                        + byn(i,j) * dbydy)

                        dkprldl = xkxsav(n) * bhatgradbx
     &                          + xkysav(m) * bhatgradby
     &                          + xkphi(i)  * bhatgradbz

                        d2kprldl2 = xkxsav(n) * bhatgrad2bx
     &                            + xkysav(m) * bhatgrad2by
     &                            + xkphi(i)  * bhatgrad2bz

                        alphae = sqrt(2. * xkte(i,j) / xme)
                        dKdL_giv   = (alphae / omgrf)**2 * dkprldl
                        d2KdL2_giv = (alphae / omgrf)**2 * d2kprldl2

                        alpha1 = sqrt(2. * xkti(i,j) / xmi1)
                        dKdL_giv_1   = (alpha1 / omgrf)**2 * dkprldl
*                       ------------------------------------
*                       Optional: leave out upshift in xkprl
*                       --------------------------------- --
                        if (upshift .eq. 0) xkprl = uzz(i, j) * xkphi(i)


                        if (upshift .eq. -1) then
                           if (xkperp .gt. xk_cutoff)
     &                                       xkprl = uzz(i,j) * xkphi(i)
                        end if

                        if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
                        if (xkperp .eq. 0.0) xkperp = 1.0e-08

                        sgn_kprl = sign(1.0, xkprl)
                        akprl = abs(xkprl)
                        sgn_dKdL = sign(1.0, dKdL_giv)
                        abs_dKdL = abs(dKdL_giv)
!                       -------------------------------------------------
!                       Optional: Don't allow xkprl to be 0 (upshift = -2)
!                       -------------------------------------------------
                        if (upshift .eq. -2)then
                           if (akprl .lt. akprl_min) then
                              xkprl = akprl_min * sgn_kprl
                           end if
                        end if

                        if(xkperp .gt. kperp_max)then
                           write (6, *)"xkperp gt kperp_max in sigmad"
                           write (15, *)"xkperp gt kperp_max in sigmad"
                        end if


                         if (use_fourier_z2 .eqv. .true. .and.
     &                                    psi(i,j) .lt. psiorb) then
*                        ----------------------------------------------
*                        Evaluate AORSA basis function along field line
*                        ----------------------------------------------
                           do nl = 1, nnodelb
                              psil = xkxsav(n) * capr_lb(nl)
     &                              + xkysav(m) * capz_lb(nl)
     &                              + nphi      * phin_lb(nl)
                              El(nl) = exp(zi * psil)
                           end do
*                     -------------------------------------------------
*                    fft on lb mesh giving nnodelb modes (nkl1 to nkl2)
*                     -------------------------------------------------
                           call fftn2(El, workl, nnodelb, nlmx,
     &                                            nkldim1, nkldim2, Eln)

*                    --------------------------------------------------
*                    Evaluate current at lb = 0, or lbprime = lbmax / 2.
*                    --------------------------------------------------
                           nl = nnodelb0

                           Esum = 0.0
                           Jzzsum = 0.0

*                          ------------------------------------------------------
*                          Calculate sigma_e (l = 0 harmonic) for each k in the
*                             Fourier expansion along the field line - xklsav(nkl)
*                          ------------------------------------------------------
                           do nkl = nkl1, nkl2

                              sige3_new = 0.0

                              call sigmad_cql3d_e(i, j, n, m, rho(i,j),
     &                        rho_a,
     &                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     &                        xme, qe, xnea(i,j), xnuome,
     &                        xkte(i,j), omgce(i,j), omgpe2(i,j),
     &                        -0, 0, nzfun, ibessel,
     &                        xkxsav(n), xkysav(m), nphi, capr(i),
     &                        bxn(i,j), byn(i,j), bzn(i,j),
     &                        uxx(i,j), uxy(i,j), uxz(i,j),
     &                        uyx(i,j), uyy(i,j), uyz(i,j),
     &                        uzx(i,j), uzy(i,j), uzz(i,j),
     &                        sigexx, sigexy, sigexz,
     &                        sigeyx, sigeyy, sigeyz,
     &                        sigezx, sigezy, sigezz,
     &                        delta0, ndiste, nupar, nuper, n_psi,
     &                        n_psi_dim, dfdupere, dfdupare,
     &                        UminPara, UmaxPara, UPERP, UPARA,
     &                        vce_mks, dfe_cql_uprp, dfe_cql_uprl,
     &                        nbessj,
     &                        nkperp, zi, eps0, v0i, omgrf, xk0,
     &                        kperp_max, i_sav, j_sav, upshift, damping,
     &                        xk_cutoff, rt, nkx2, nky2, xklsav(nkl),
     &                        xkperp,
     &                        xkalp, xkbet,  z0_table1, z1_table1,
     &                        z2_table1, zetai_table, dKdL_table,
     &                        dKdL_giv, nmax, mmax, use_new_z2, ntable,
     &                        mtable, sige3_new)

                              ekl = cexp(zi * xklsav(nkl) * lbprime(nl))

                              Esum   = Esum   + Eln(nkl) * ekl
                              Jzzsum = Jzzsum + Eln(nkl) * sigezz * ekl
                           end do
                           sige3_new  = Jzzsum / Esum

*                         ----------------------------------------------
*                          Recalulate sigma_e using all harmonics, but
*                          replacing the l = 0 term in sigezz
*                          with the Fourier expanded term sige3_new
*                         ----------------------------------------------

                           call sigmad_cql3d_e(i, j, n, m, rho(i,j),
     &                     rho_a,
     &                     gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     &                     xme, qe, xnea(i,j), xnuome,
     &                     xkte(i,j), omgce(i,j), omgpe2(i,j),
     &                     -lmaxe, lmaxe, nzfun, ibessel,
     &                     xkxsav(n), xkysav(m), nphi, capr(i),
     &                     bxn(i,j), byn(i,j), bzn(i,j),
     &                     uxx(i,j), uxy(i,j), uxz(i,j),
     &                     uyx(i,j), uyy(i,j), uyz(i,j),
     &                     uzx(i,j), uzy(i,j), uzz(i,j),
     &                     sigexx, sigexy, sigexz,
     &                     sigeyx, sigeyy, sigeyz,
     &                     sigezx, sigezy, sigezz,
     &                     delta0, ndiste, nupar, nuper, n_psi,
     &                     n_psi_dim, dfdupere, dfdupare,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vce_mks, dfe_cql_uprp, dfe_cql_uprl, nbessj,
     &                     nkperp, zi, eps0, v0i, omgrf, xk0,
     &                     kperp_max, i_sav, j_sav, upshift, damping,
     &                     xk_cutoff, rt, nkx2, nky2, xkprl, xkperp,
     &                     xkalp, xkbet,  z0_table1, z1_table1,
     &                     z2_table1, zetai_table, dKdL_table,
     &                     dKdL_giv, nmax, mmax, use_new_z2, ntable,
     &                     mtable, sige3_new)

                        end if !  end if use_fourier_z2 = true and psi .lt. psiorb


                        if(use_fourier_z2 .eqv. .false.) then

                           sige3_new = 0.0

                           call sigmad_cql3d_e(i, j, n, m, rho(i,j),
     &                        rho_a, gradprlb(i,j), bmod(i,j),
     &                        bmod_mid(i,j), xme, qe, xnea(i,j),
     &                        xnuome,
     &                        xkte(i,j), omgce(i,j), omgpe2(i,j),
     &                        -lmaxe, lmaxe, nzfun, ibessel,
     &                        xkxsav(n), xkysav(m), nphi, capr(i),
     &                        bxn(i,j), byn(i,j), bzn(i,j),
     &                        uxx(i,j), uxy(i,j), uxz(i,j),
     &                        uyx(i,j), uyy(i,j), uyz(i,j),
     &                        uzx(i,j), uzy(i,j), uzz(i,j),
     &                        sigexx, sigexy, sigexz,
     &                        sigeyx, sigeyy, sigeyz,
     &                        sigezx, sigezy, sigezz,
     &                        delta0, ndiste, nupar, nuper, n_psi,
     &                        n_psi_dim, dfdupere, dfdupare,
     &                        UminPara, UmaxPara, UPERP, UPARA,
     &                        vce_mks, dfe_cql_uprp, dfe_cql_uprl,
     &                        nbessj,
     &                        nkperp, zi, eps0, v0i, omgrf, xk0,
     &                        kperp_max, i_sav, j_sav, upshift, damping,
     &                        xk_cutoff, rt, nkx2, nky2, xkprl,
     &                        xkperp,
     &                        xkalp, xkbet,  z0_table1, z1_table1,
     &                        z2_table1, zetai_table, dKdL_table,
     &                        dKdL_giv, nmax, mmax, use_new_z2,
     &                        ntable, mtable, sige3_new)

                        end if    !  if   use_fourier_z2 .eq. .false.


                        if(use_fourier_z2 .eqv. .true. .and.
     &                                   psi(i,j) .ge. psiorb) then

                           sige3_new = 0.0

                           call sigmad_cql3d_e(i, j, n, m, rho(i,j),
     &                        rho_a, gradprlb(i,j), bmod(i,j),
     &                        bmod_mid(i,j), xme, qe, xnea(i,j),
     &                        xnuome,
     &                        xkte(i,j), omgce(i,j), omgpe2(i,j),
     &                        -lmaxe, lmaxe, nzfun, ibessel,
     &                        xkxsav(n), xkysav(m), nphi, capr(i),
     &                        bxn(i,j), byn(i,j), bzn(i,j),
     &                        uxx(i,j), uxy(i,j), uxz(i,j),
     &                        uyx(i,j), uyy(i,j), uyz(i,j),
     &                        uzx(i,j), uzy(i,j), uzz(i,j),
     &                        sigexx, sigexy, sigexz,
     &                        sigeyx, sigeyy, sigeyz,
     &                        sigezx, sigezy, sigezz,
     &                        delta0, ndiste, nupar, nuper, n_psi,
     &                        n_psi_dim, dfdupere, dfdupare,
     &                        UminPara, UmaxPara, UPERP, UPARA,
     &                        vce_mks, dfe_cql_uprp, dfe_cql_uprl,
     &                        nbessj,
     &                        nkperp, zi, eps0, v0i, omgrf, xk0,
     &                        kperp_max, i_sav, j_sav, upshift, damping,
     &                        xk_cutoff, rt, nkx2, nky2, xkprl,
     &                        xkperp,
     &                        xkalp, xkbet,  z0_table1, z1_table1,
     &                        z2_table1, zetai_table, dKdL_table,
     &                        dKdL_giv, nmax, mmax, .true.,
     &                        ntable, mtable, sige3_new)

                        end if ! use_fourier_z2 .eq. .true. .and. psi.ge. psiorb



                        if(odd_order .ne. 0)
     &                     call odd_order_derivs(i, j, n, m, nxmx, nymx,
     &                     nnodex, nnodey, mkdim1, mkdim2, capr,
     &                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2,
     &                     xme, qe, xnea, xkte, omgce, omgpe2, lmaxe,
     &                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi,
     &                     bxn, byn, bzn,
     &                     uxx, uxy, uxz,
     &                     uyx, uyy, uyz,
     &                     uzx, uzy, uzz,
     &                     sigexx, sigexy, sigexz,
     &                     sigeyx, sigeyy, sigeyz,
     &                     sigezx, sigezy, sigezz,
     &                     myrow, mycol, nprow, npcol, icontxt,
     &                     desc_amat, dlen_, nboundary,
     &                     xx, yy, isigma, xnuome, psi, psilim,
     &                     myid, nproc, delta0, gradprlb, bmod,
     &                     bmod_mid, ndiste, nupar, nuper, n_psi,
     &                     n_psi_dim, dfdupere, dfdupare,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vce_mks, dfe_cql_uprp, dfe_cql_uprl, rho,
     &                     rho_a, nbessj, nkperp,
     &                     zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift,
     &                     damping, xk_cutoff, rt, dx, dy,
     &                     xkalp, xkbet, xkprl)



                        call sigmad_cql3d_1(i, j, n, m, rho(i,j), rho_a,
     &                     gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     &                     xmi1, qi1, xn1a(i,j), xnuomi,
     &                     xkti(i,j), omgci1(i,j), omgp12(i,j),
     &                     -lmax, lmax, nzfun, ibessel,
     &                     xkxsav(n), xkysav(m), nphi, capr(i),
     &                     bxn(i,j), byn(i,j), bzn(i,j),
     &                     uxx(i,j), uxy(i,j), uxz(i,j),
     &                     uyx(i,j), uyy(i,j), uyz(i,j),
     &                     uzx(i,j), uzy(i,j), uzz(i,j),
     &                     sig1xx, sig1xy, sig1xz,
     &                     sig1yx, sig1yy, sig1yz,
     &                     sig1zx, sig1zy, sig1zz,
     &                     delta0, ndisti1, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper1, dfdupar1,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc1_mks, df1_cql_uprp, df1_cql_uprl, nbessj,
     &                     nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav1, j_sav1, upshift, damping, xk_cutoff,
     &                     rt, nkx2, nky2, xkprl, xkperp, xkalp, xkbet,
     &                     z0_table1, z1_table1,z2_table1, zetai_table,
     &                     dKdL_table,
     &                     dKdL_giv_1, nmax, mmax, use_new_z2, ntable,
     &                     mtable)

                        if(odd_order .ne. 0)
     &                     call odd_order_derivs(i, j, n, m, nxmx, nymx,
     &                     nnodex, nnodey, mkdim1, mkdim2, capr,
     &                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2,
     &                     xmi1, qi1, xn1a, xkti, omgci1, omgp12, lmax,
     &                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi,
     &                     bxn, byn, bzn,
     &                     uxx, uxy, uxz,
     &                     uyx, uyy, uyz,
     &                     uzx, uzy, uzz,
     &                     sig1xx, sig1xy, sig1xz,
     &                     sig1yx, sig1yy, sig1yz,
     &                     sig1zx, sig1zy, sig1zz,
     &                     myrow, mycol, nprow, npcol, icontxt,
     &                     desc_amat, dlen_, nboundary,
     &                     xx, yy, isigma, xnuomi, psi, psilim,
     &                     myid, nproc, delta0, gradprlb, bmod,
     &                     bmod_mid, ndisti1, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper1, dfdupar1,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc1_mks, df1_cql_uprp, df1_cql_uprl, rho,
     &                     rho_a, nbessj, nkperp,
     &                     zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift,
     &                     damping, xk_cutoff, rt, dx, dy,
     &                     xkalp, xkbet, xkprl)


                        if (eta2 .ne. 0.0) then
                        call sigmad_cql3d_2(i, j, n, m, rho(i,j), rho_a,
     &                     gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     &                     xmi2, qi2, xn2a(i,j), xnuomi,
     &                     xkti2(i,j), omgci2(i,j), omgp22(i,j),
     &                     -lmax, lmax, nzfun, ibessel,
     &                     xkxsav(n), xkysav(m), nphi, capr(i),
     &                     bxn(i,j), byn(i,j), bzn(i,j),
     &                     uxx(i,j), uxy(i,j), uxz(i,j),
     &                     uyx(i,j), uyy(i,j), uyz(i,j),
     &                     uzx(i,j), uzy(i,j), uzz(i,j),
     &                     sig2xx, sig2xy, sig2xz,
     &                     sig2yx, sig2yy, sig2yz,
     &                     sig2zx, sig2zy, sig2zz,
     &                     delta0, ndisti2, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper2, dfdupar2,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc2_mks, df2_cql_uprp, df2_cql_uprl, nbessj,
     &                     nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav2, j_sav2, upshift, damping, xk_cutoff,
     &                     rt, nkx2, nky2, xkprl, xkperp, xkalp, xkbet)

                        if(odd_order .ne. 0)
     &                     call odd_order_derivs(i, j, n, m, nxmx, nymx,
     &                     nnodex, nnodey, mkdim1, mkdim2, capr,
     &                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2,
     &                     xmi2, qi2, xn2a, xkti2, omgci2, omgp22, lmax,
     &                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi,
     &                     bxn, byn, bzn,
     &                     uxx, uxy, uxz,
     &                     uyx, uyy, uyz,
     &                     uzx, uzy, uzz,
     &                     sig2xx, sig2xy, sig2xz,
     &                     sig2yx, sig2yy, sig2yz,
     &                     sig2zx, sig2zy, sig2zz,
     &                     myrow, mycol, nprow, npcol, icontxt,
     &                     desc_amat, dlen_, nboundary,
     &                     xx, yy, isigma, xnuomi, psi, psilim,
     &                     myid, nproc, delta0, gradprlb, bmod,
     &                     bmod_mid, ndisti2, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper2, dfdupar2,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc2_mks, df2_cql_uprp, df2_cql_uprl, rho,
     &                     rho_a, nbessj, nkperp,
     &                     zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav2, j_sav2, upshift,
     &                     damping, xk_cutoff, rt, dx, dy,
     &                     xkalp, xkbet, xkprl)

                        end if    !   eta2 .ne. 0.0


                        if(eta3 .ne. 0.0) then
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a,
     &                     gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     &                     xmi3, qi3, xn3a(i,j), xnuomi,
     &                     xkti3(i,j), omgci3(i,j), omgp32(i,j),
     &                     -lmax, lmax, nzfun, ibessel,
     &                     xkxsav(n), xkysav(m), nphi, capr(i),
     &                     bxn(i,j), byn(i,j), bzn(i,j),
     &                     uxx(i,j), uxy(i,j), uxz(i,j),
     &                     uyx(i,j), uyy(i,j), uyz(i,j),
     &                     uzx(i,j), uzy(i,j), uzz(i,j),
     &                     sig3xx, sig3xy, sig3xz,
     &                     sig3yx, sig3yy, sig3yz,
     &                     sig3zx, sig3zy, sig3zz,
     &                     delta0, ndisti3, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper3, dfdupar3,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc3_mks, df3_cql_uprp, df3_cql_uprl, nbessj,
     &                     nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift, damping, xk_cutoff,
     &                     rt, nkx2, nky2, xkprl, xkperp, xkalp, xkbet)

                        if(odd_order .ne. 0)
     &                     call odd_order_derivs(i, j, n, m, nxmx, nymx,
     &                     nnodex, nnodey, mkdim1, mkdim2, capr,
     &                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2,
     &                     xmi3, qi3, xn3a, xkti3, omgci3, omgp32, lmax,
     &                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi,
     &                     bxn, byn, bzn,
     &                     uxx, uxy, uxz,
     &                     uyx, uyy, uyz,
     &                     uzx, uzy, uzz,
     &                     sig3xx, sig3xy, sig3xz,
     &                     sig3yx, sig3yy, sig3yz,
     &                     sig3zx, sig3zy, sig3zz,
     &                     myrow, mycol, nprow, npcol, icontxt,
     &                     desc_amat, dlen_, nboundary,
     &                     xx, yy, isigma, xnuomi, psi, psilim,
     &                     myid, nproc, delta0, gradprlb, bmod,
     &                     bmod_mid, ndisti3, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper3, dfdupar3,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc3_mks, df3_cql_uprp, df3_cql_uprl, rho,
     &                     rho_a, nbessj, nkperp,
     &                     zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift,
     &                     damping, xk_cutoff, rt, dx, dy,
     &                     xkalp, xkbet, xkprl)
                        end if    !   eta3 .ne. 0.0

                        if(eta4 .ne. 0.0) then
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a,
     &                     gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     &                     xmi4, qi4, xn4a(i,j), xnuomi,
     &                     xkti4(i,j), omgci4(i,j), omgp42(i,j),
     &                     -lmax, lmax, nzfun, ibessel,
     &                     xkxsav(n), xkysav(m), nphi, capr(i),
     &                     bxn(i,j), byn(i,j), bzn(i,j),
     &                     uxx(i,j), uxy(i,j), uxz(i,j),
     &                     uyx(i,j), uyy(i,j), uyz(i,j),
     &                     uzx(i,j), uzy(i,j), uzz(i,j),
     &                     sig4xx, sig4xy, sig4xz,
     &                     sig4yx, sig4yy, sig4yz,
     &                     sig4zx, sig4zy, sig4zz,
     &                     delta0, ndisti4, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper4, dfdupar4,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc4_mks, df4_cql_uprp, df4_cql_uprl, nbessj,
     &                     nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift, damping, xk_cutoff,
     &                     rt, nkx2, nky2, xkprl, xkperp, xkalp, xkbet)

                        if(odd_order .ne. 0)
     &                     call odd_order_derivs(i, j, n, m, nxmx, nymx,
     &                     nnodex, nnodey, mkdim1, mkdim2, capr,
     &                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2,
     &                     xmi4, qi4, xn4a, xkti4, omgci4, omgp42, lmax,
     &                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi,
     &                     bxn, byn, bzn,
     &                     uxx, uxy, uxz,
     &                     uyx, uyy, uyz,
     &                     uzx, uzy, uzz,
     &                     sig4xx, sig4xy, sig4xz,
     &                     sig4yx, sig4yy, sig4yz,
     &                     sig4zx, sig4zy, sig4zz,
     &                     myrow, mycol, nprow, npcol, icontxt,
     &                     desc_amat, dlen_, nboundary,
     &                     xx, yy, isigma, xnuomi, psi, psilim,
     &                     myid, nproc, delta0, gradprlb, bmod,
     &                     bmod_mid, ndisti4, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper4, dfdupar4,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc4_mks, df4_cql_uprp, df4_cql_uprl, rho,
     &                     rho_a, nbessj, nkperp,
     &                     zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift,
     &                     damping, xk_cutoff, rt, dx, dy,
     &                     xkalp, xkbet, xkprl)

                        end if     !   eta4 .ne. 0.0

                        if(eta5 .ne. 0.0) then
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a,
     &                     gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     &                     xmi5, qi5, xn5a(i,j), xnuomi,
     &                     xkti5(i,j), omgci5(i,j), omgp52(i,j),
     &                     -lmax, lmax, nzfun, ibessel,
     &                     xkxsav(n), xkysav(m), nphi, capr(i),
     &                     bxn(i,j), byn(i,j), bzn(i,j),
     &                     uxx(i,j), uxy(i,j), uxz(i,j),
     &                     uyx(i,j), uyy(i,j), uyz(i,j),
     &                     uzx(i,j), uzy(i,j), uzz(i,j),
     &                     sig5xx, sig5xy, sig5xz,
     &                     sig5yx, sig5yy, sig5yz,
     &                     sig5zx, sig5zy, sig5zz,
     &                     delta0, ndisti5, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper5, dfdupar5,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc5_mks, df5_cql_uprp, df5_cql_uprl, nbessj,
     &                     nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift, damping, xk_cutoff,
     &                     rt, nkx2, nky2, xkprl, xkperp, xkalp, xkbet)

                        if(odd_order .ne. 0)
     &                     call odd_order_derivs(i, j, n, m, nxmx, nymx,
     &                     nnodex, nnodey, mkdim1, mkdim2, capr,
     &                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2,
     &                     xmi5, qi5, xn5a, xkti5, omgci5, omgp52, lmax,
     &                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi,
     &                     bxn, byn, bzn,
     &                     uxx, uxy, uxz,
     &                     uyx, uyy, uyz,
     &                     uzx, uzy, uzz,
     &                     sig5xx, sig5xy, sig5xz,
     &                     sig5yx, sig5yy, sig5yz,
     &                     sig5zx, sig5zy, sig5zz,
     &                     myrow, mycol, nprow, npcol, icontxt,
     &                     desc_amat, dlen_, nboundary,
     &                     xx, yy, isigma, xnuomi, psi, psilim,
     &                     myid, nproc, delta0, gradprlb, bmod,
     &                     bmod_mid, ndisti5, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper5, dfdupar5,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc5_mks, df5_cql_uprp, df5_cql_uprl, rho,
     &                     rho_a, nbessj, nkperp,
     &                     zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift,
     &                     damping, xk_cutoff, rt, dx, dy,
     &                     xkalp, xkbet, xkprl)

                        end if  !  eta5 .ne. 0.0

                        if(eta6 .ne. 0.0) then
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a,
     &                     gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     &                     xmi6, qi6, xn6a(i,j), xnuomi,
     &                     xkti6(i,j), omgci6(i,j), omgp62(i,j),
     &                     -lmax, lmax, nzfun, ibessel,
     &                     xkxsav(n), xkysav(m), nphi, capr(i),
     &                     bxn(i,j), byn(i,j), bzn(i,j),
     &                     uxx(i,j), uxy(i,j), uxz(i,j),
     &                     uyx(i,j), uyy(i,j), uyz(i,j),
     &                     uzx(i,j), uzy(i,j), uzz(i,j),
     &                     sig6xx, sig6xy, sig6xz,
     &                     sig6yx, sig6yy, sig6yz,
     &                     sig6zx, sig6zy, sig6zz,
     &                     delta0, ndisti6, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper6, dfdupar6,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc6_mks, df6_cql_uprp, df6_cql_uprl, nbessj,
     &                     nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift, damping, xk_cutoff,
     &                     rt, nkx2, nky2, xkprl, xkperp, xkalp, xkbet)

                        if(odd_order .ne. 0)
     &                     call odd_order_derivs(i, j, n, m, nxmx, nymx,
     &                     nnodex, nnodey, mkdim1, mkdim2, capr,
     &                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2,
     &                     xmi6, qi6, xn6a, xkti6, omgci6, omgp62, lmax,
     &                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi,
     &                     bxn, byn, bzn,
     &                     uxx, uxy, uxz,
     &                     uyx, uyy, uyz,
     &                     uzx, uzy, uzz,
     &                     sig6xx, sig6xy, sig6xz,
     &                     sig6yx, sig6yy, sig6yz,
     &                     sig6zx, sig6zy, sig6zz,
     &                     myrow, mycol, nprow, npcol, icontxt,
     &                     desc_amat, dlen_, nboundary,
     &                     xx, yy, isigma, xnuomi, psi, psilim,
     &                     myid, nproc, delta0, gradprlb, bmod,
     &                     bmod_mid, ndisti6, nupar, nuper, n_psi,
     &                     n_psi_dim, dfduper6, dfdupar6,
     &                     UminPara, UmaxPara, UPERP, UPARA,
     &                     vc6_mks, df6_cql_uprp, df6_cql_uprl, rho,
     &                     rho_a, nbessj, nkperp,
     &                     zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                     i_sav, j_sav, upshift,
     &                     damping, xk_cutoff, rt, dx, dy,
     &                     xkalp, xkbet, xkprl)

                        end if  !  eta6 .ne. 0.0

                        if(eta_slo .ne. 0.0)
     &                  call sigmah_slow(i, j, n, m,
     &                  xmi_slo, qi_slo, xna_slo(i,j), xnuomi,
     &                  eslow, omgci_slo(i,j), omgp2_slo(i,j),
     &                  -lmax, lmax, nzfun, ibessel,
     &                  xkxsav(n), xkysav(m), nphi, capr(i),
     &                  bxn(i,j), byn(i,j), bzn(i,j),
     &                  uxx(i,j), uxy(i,j), uxz(i,j),
     &                  uyx(i,j), uyy(i,j), uyz(i,j),
     &                  uzx(i,j), uzy(i,j), uzz(i,j),
     &                  sigsloxx, sigsloxy, sigsloxz,
     &                  sigsloyx, sigsloyy, sigsloyz,
     &                  sigslozx, sigslozy, sigslozz,
     &                  xkte(i,j), zeffcd, zi, eps0, v0i, omgrf, xk0,
     &                  kperp_max, i_sav, j_sav)


                     end if  ! end if isigma .eq. 1





                     if (isigma .eq. 0)then

                        call sigmac_stix(i, j, n, m,
     &                  xme, qe, xnea(i,j), xnuome,
     &                  xkte(i,j), omgce(i,j), omgpe2(i,j),
     &                  -lmaxe, lmaxe, nzfun, ibessel,
     &                  xkxsav(n), xkysav(m), nphi, capr(i),
     &                  bxn(i,j), byn(i,j), bzn(i,j),
     &                  uxx(i,j), uxy(i,j), uxz(i,j),
     &                  uyx(i,j), uyy(i,j), uyz(i,j),
     &                  uzx(i,j), uzy(i,j), uzz(i,j),
     &                  sigexx, sigexy, sigexz,
     &                  sigeyx, sigeyy, sigeyz,
     &                  sigezx, sigezy, sigezz,
     &                  delta0, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                  i_sav, j_sav)

                        call sigmac_stix(i, j, n, m,
     &                  xmi1, qi1, xn1a(i,j), xnuomi,
     &                  xkti(i,j), omgci1(i,j), omgp12(i,j),
     &                  -lmax, lmax, nzfun, ibessel,
     &                  xkxsav(n), xkysav(m), nphi, capr(i),
     &                  bxn(i,j), byn(i,j), bzn(i,j),
     &                  uxx(i,j), uxy(i,j), uxz(i,j),
     &                  uyx(i,j), uyy(i,j), uyz(i,j),
     &                  uzx(i,j), uzy(i,j), uzz(i,j),
     &                  sig1xx, sig1xy, sig1xz,
     &                  sig1yx, sig1yy, sig1yz,
     &                  sig1zx, sig1zy, sig1zz,
     &                  delta0, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                  i_sav, j_sav)

                        if(eta2 .ne. 0.0)
     &                  call sigmac_stix(i, j, n, m,
     &                  xmi2, qi2, xn2a(i,j), xnuomi,
     &                  xkti2(i,j), omgci2(i,j), omgp22(i,j),
     &                  -lmax, lmax, nzfun, ibessel,
     &                  xkxsav(n), xkysav(m), nphi, capr(i),
     &                  bxn(i,j), byn(i,j), bzn(i,j),
     &                  uxx(i,j), uxy(i,j), uxz(i,j),
     &                  uyx(i,j), uyy(i,j), uyz(i,j),
     &                  uzx(i,j), uzy(i,j), uzz(i,j),
     &                  sig2xx, sig2xy, sig2xz,
     &                  sig2yx, sig2yy, sig2yz,
     &                  sig2zx, sig2zy, sig2zz,
     &                  delta0, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                  i_sav, j_sav)

                        if(eta3 .ne. 0.0)
     &                  call sigmac_stix(i, j, n, m,
     &                  xmi3, qi3, xn3a(i,j), xnuomi,
     &                  xkti3(i,j), omgci3(i,j), omgp32(i,j),
     &                  -lmax, lmax, nzfun, ibessel,
     &                  xkxsav(n), xkysav(m), nphi, capr(i),
     &                  bxn(i,j), byn(i,j), bzn(i,j),
     &                  uxx(i,j), uxy(i,j), uxz(i,j),
     &                  uyx(i,j), uyy(i,j), uyz(i,j),
     &                  uzx(i,j), uzy(i,j), uzz(i,j),
     &                  sig3xx, sig3xy, sig3xz,
     &                  sig3yx, sig3yy, sig3yz,
     &                  sig3zx, sig3zy, sig3zz,
     &                  delta0, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                  i_sav, j_sav)


                        if(eta4 .ne. 0.0)
     &                  call sigmac_stix(i, j, n, m,
     &                  xmi4, qi4, xn4a(i,j), xnuomi,
     &                  xkti4(i,j), omgci4(i,j), omgp42(i,j),
     &                  -lmax, lmax, nzfun, ibessel,
     &                  xkxsav(n), xkysav(m), nphi, capr(i),
     &                  bxn(i,j), byn(i,j), bzn(i,j),
     &                  uxx(i,j), uxy(i,j), uxz(i,j),
     &                  uyx(i,j), uyy(i,j), uyz(i,j),
     &                  uzx(i,j), uzy(i,j), uzz(i,j),
     &                  sig4xx, sig4xy, sig4xz,
     &                  sig4yx, sig4yy, sig4yz,
     &                  sig4zx, sig4zy, sig4zz,
     &                  delta0, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                  i_sav, j_sav)

                        if(eta5 .ne. 0.0)
     &                  call sigmac_stix(i, j, n, m,
     &                  xmi5, qi5, xn5a(i,j), xnuomi,
     &                  xkti5(i,j), omgci5(i,j), omgp52(i,j),
     &                  -lmax, lmax, nzfun, ibessel,
     &                  xkxsav(n), xkysav(m), nphi, capr(i),
     &                  bxn(i,j), byn(i,j), bzn(i,j),
     &                  uxx(i,j), uxy(i,j), uxz(i,j),
     &                  uyx(i,j), uyy(i,j), uyz(i,j),
     &                  uzx(i,j), uzy(i,j), uzz(i,j),
     &                  sig5xx, sig5xy, sig5xz,
     &                  sig5yx, sig5yy, sig5yz,
     &                  sig5zx, sig5zy, sig5zz,
     &                  delta0, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                  i_sav, j_sav)

                        if(eta6 .ne. 0.0)
     &                  call sigmac_stix(i, j, n, m,
     &                  xmi6, qi6, xn6a(i,j), xnuomi,
     &                  xkti6(i,j), omgci6(i,j), omgp62(i,j),
     &                  -lmax, lmax, nzfun, ibessel,
     &                  xkxsav(n), xkysav(m), nphi, capr(i),
     &                  bxn(i,j), byn(i,j), bzn(i,j),
     &                  uxx(i,j), uxy(i,j), uxz(i,j),
     &                  uyx(i,j), uyy(i,j), uyz(i,j),
     &                  uzx(i,j), uzy(i,j), uzz(i,j),
     &                  sig6xx, sig6xy, sig6xz,
     &                  sig6yx, sig6yy, sig6yz,
     &                  sig6zx, sig6zy, sig6zz,
     &                  delta0, zi, eps0, v0i, omgrf, xk0, kperp_max,
     &                  i_sav, j_sav)

                     end if  ! end if isigma .eq. 0


                     sigxx = sigexx + sig1xx + sig2xx + sig3xx
     &                           + sig4xx + sig5xx + sig6xx + sigsloxx
                     sigxy = sigexy + sig1xy + sig2xy + sig3xy
     &                           + sig4xy + sig5xy + sig6xy + sigsloxy
                     sigxz = sigexz + sig1xz + sig2xz + sig3xz
     &                           + sig4xz + sig5xz + sig6xz + sigsloxz

                     sigyx = sigeyx + sig1yx + sig2yx + sig3yx
     &                           + sig4yx + sig5yx + sig6yx + sigsloyx
                     sigyy = sigeyy + sig1yy + sig2yy + sig3yy
     &                           + sig4yy + sig5yy + sig6yy + sigsloyy
                     sigyz = sigeyz + sig1yz + sig2yz + sig3yz
     &                           + sig4yz + sig5yz + sig6yz + sigsloyz

                     sigzx = sigezx + sig1zx + sig2zx + sig3zx
     &                           + sig4zx + sig5zx + sig6zx + sigslozx
                     sigzy = sigezy + sig1zy + sig2zy + sig3zy
     &                           + sig4zy + sig5zy + sig6zy + sigslozy
                     sigzz = sigezz + sig1zz + sig2zz + sig3zz
     &                           + sig4zz + sig5zz + sig6zz + sigslozz


                     xkxx = 1.0 + zi / (eps0 * omgrf) * sigxx
                     xkxy =       zi / (eps0 * omgrf) * sigxy
                     xkxz =       zi / (eps0 * omgrf) * sigxz

                     xkyx =       zi / (eps0 * omgrf) * sigyx
                     xkyy = 1.0 + zi / (eps0 * omgrf) * sigyy
                     xkyz =       zi / (eps0 * omgrf) * sigyz

                     xkzx =       zi / (eps0 * omgrf) * sigzx
                     xkzy =       zi / (eps0 * omgrf) * sigzy
                     xkzz = 1.0 + zi / (eps0 * omgrf) * sigzz


                     rnx = xkxsav(n) / xk0
                     rny = xkysav(m) / xk0
                     rnphi = xkphi(i) / xk0

                     dxx = (xkxx - rny**2 - rnphi**2) * uxx(i,j)
     &                +  xkyx * uyx(i,j)
     &                +  xkzx * uzx(i,j)
     &                + rnx * (rny * uxy(i,j) + rnphi * uxz(i,j))
     &                - zi * rnphi / xk0 *
     &                         (uxz(i,j) / capr(i) + dxuxz(i,j))
     &                - zi * rny / xk0 * (dxuxy(i,j) - 2. * dyuxx(i,j))
     &                - zi * rnx / xk0 * dyuxy(i,j)
     &                + 1. / xk0**2 * (dyyuxx(i,j) - dxyuxy(i,j))

                     dxy =  xkxy * uxx(i,j)
     &                + (xkyy - rny**2 - rnphi**2) * uyx(i,j)
     &                +  xkzy * uzx(i,j)
     &                + rnx * (rny * uyy(i,j) + rnphi * uyz(i,j))
     &                - zi * rnphi / xk0 *
     &                         (uyz(i,j) / capr(i) + dxuyz(i,j))
     &                - zi * rny / xk0 * (dxuyy(i,j) - 2. * dyuyx(i,j))
     &                - zi * rnx / xk0 * dyuyy(i,j)
     &                + 1. / xk0**2 * (dyyuyx(i,j) - dxyuyy(i,j))

                     dxz =  xkxz * uxx(i,j)
     &                +  xkyz * uyx(i,j)
     &                + (xkzz - rny**2 - rnphi**2) * uzx(i,j)
     &                + rnx * (rny * uzy(i,j) + rnphi * uzz(i,j))
     &                - zi * rnphi / xk0 *
     &                         (uzz(i,j) / capr(i) + dxuzz(i,j))
     &                - zi * rny / xk0 * (dxuzy(i,j) - 2. * dyuzx(i,j))
     &                - zi * rnx / xk0 * dyuzy(i,j)
     &                + 1. / xk0**2 * (dyyuzx(i,j) - dxyuzy(i,j))

                     dyx = (xkxx - rnx**2 - rnphi**2) * uxy(i,j)
     &                +  xkyx * uyy(i,j)
     &                +  xkzx * uzy(i,j)
     &                + rny * (rnx * uxx(i,j) + rnphi * uxz(i,j))
     &                - zi * rny / xk0 *
     &                         (dxuxx(i,j) + uxx(i,j) / capr(i))
     &                - zi * rnphi / xk0 * dyuxz(i,j)
     &                - zi * rnx / xk0 *
     &                (dyuxx(i,j) - uxy(i,j) / capr(i) - 2.* dxuxy(i,j))
     &                + 1. / xk0**2 * (dxxuxy(i,j) - dxyuxx(i,j)
     &                - dyuxx(i,j)/ capr(i) + dxuxy(i,j) / capr(i))

                     dyy =  xkxy * uxy(i,j)
     &                + (xkyy - rnx**2 - rnphi**2) * uyy(i,j)
     &                +  xkzy * uzy(i,j)
     &                + rny * (rnx * uyx(i,j) + rnphi * uyz(i,j))
     &                - zi * rny / xk0 *
     &                         (dxuyx(i,j) + uyx(i,j) / capr(i))
     &                - zi * rnphi / xk0 * dyuyz(i,j)
     &                - zi * rnx / xk0 *
     &                (dyuyx(i,j) - uyy(i,j) / capr(i) - 2.* dxuyy(i,j))
     &                + 1. / xk0**2 * (dxxuyy(i,j) - dxyuyx(i,j)
     &                - dyuyx(i,j)/ capr(i) + dxuyy(i,j) / capr(i))

                     dyz =  xkxz * uxy(i,j)
     &                +  xkyz * uyy(i,j)
     &                + (xkzz - rnx**2 - rnphi**2) * uzy(i,j)
     &                + rny * (rnx * uzx(i,j) + rnphi * uzz(i,j))
     &                - zi * rny / xk0 *
     &                         (dxuzx(i,j) + uzx(i,j) / capr(i))
     &                - zi * rnphi / xk0 * dyuzz(i,j)
     &                - zi * rnx / xk0 *
     &                (dyuzx(i,j) - uzy(i,j) / capr(i) - 2.* dxuzy(i,j))
     &                + 1. / xk0**2 * (dxxuzy(i,j) - dxyuzx(i,j)
     &                - dyuzx(i,j)/ capr(i) + dxuzy(i,j) / capr(i))

                     dzx = (xkxx - rnx**2 - rny**2) * uxz(i,j)
     &                +  xkyx * uyz(i,j)
     &                +  xkzx * uzz(i,j)
     &                + rnphi * (rny * uxy(i,j) + rnx * uxx(i,j))
     &                + zi * rny / xk0 * 2. * dyuxz(i,j)
     &                - zi * rnphi / xk0 *
     &                   (dyuxy(i,j) + dxuxx(i,j) - uxx(i,j) / capr(i))
     &                + zi * rnx / xk0 *
     &                   (uxz(i,j) / capr(i) + 2.* dxuxz(i,j))
     &                - 1. / (xk0**2 * capr(i)) *
     &                   (uxz(i,j)/ capr(i) - dxuxz(i,j))
     &                + 1. / xk0**2  * (dxxuxz(i,j) + dyyuxz(i,j))

                     dzy =  xkxy * uxz(i,j)
     &                + (xkyy - rnx**2 - rny**2) * uyz(i,j)
     &                +  xkzy * uzz(i,j)
     &                + rnphi * (rny * uyy(i,j) + rnx * uyx(i,j))
     &                + zi * rny / xk0 * 2. * dyuyz(i,j)
     &                - zi * rnphi / xk0 *
     &                   (dyuyy(i,j) + dxuyx(i,j) - uyx(i,j) / capr(i))
     &                + zi * rnx / xk0 *
     &                   (uyz(i,j) / capr(i) + 2.* dxuyz(i,j))
     &                - 1. / (xk0**2 * capr(i)) *
     &                   (uyz(i,j)/ capr(i) - dxuyz(i,j))
     &                + 1. / xk0**2  * (dxxuyz(i,j) + dyyuyz(i,j))

                     dzz =  xkxz * uxz(i,j)
     &                +  xkyz * uyz(i,j)
     &                + (xkzz - rnx**2 - rny**2) * uzz(i,j)
     &                + rnphi * (rny * uzy(i,j) + rnx * uzx(i,j))
     &                + zi * rny / xk0 * 2. * dyuzz(i,j)
     &                - zi * rnphi / xk0 *
     &                   (dyuzy(i,j) + dxuzx(i,j) - uzx(i,j) / capr(i))
     &                + zi * rnx / xk0 *
     &                   (uzz(i,j) / capr(i) + 2.* dxuzz(i,j))
     &                - 1. / (xk0**2 * capr(i)) *
     &                   (uzz(i,j)/ capr(i) - dxuzz(i,j))
     &                + 1. / xk0**2  * (dxxuzz(i,j) + dyyuzz(i,j))


                     fdk = dxx * cexpkxky
                     fek = dxy * cexpkxky
                     ffk = dxz * cexpkxky

                     fgk = dyx * cexpkxky
                     fak = dyy * cexpkxky
                     fpk = dyz * cexpkxky

                     frk = dzx * cexpkxky
                     fqk = dzy * cexpkxky
                     fsk = dzz * cexpkxky

                     detd = fsk * (fdk * fak - fek * fgk)
     &                 - fqk * (fdk * fpk - ffk * fgk)
     &                 + frk * (fek * fpk - ffk * fak)

*                    ------------------------------------------------------
*                    plot all orders dispersion kernel D(x, kx) = capd(i, n)
*                    ------------------------------------------------------
                     if(j .eq. jdisp .and. m .eq. mdisp)then

                        capd(i, n) = cabs(detd)

                        xkb(i, n) = uzx(i, j) * xkxsav(n)
     &                            + uzy(i, j) * xkysav(m)
     &                            + uzz(i, j) * xkphi(i)
                     end if

                  end if  ! end if only do sigma in interior when nboundary = 1


*                 --------------------------
*                 metal boundary edge region:
*                 --------------------------
                  if(psi(i,j) .gt. psilim .and. nboundary .eq. 1
     &               .or. nboundary .eq. 0 .and. i .eq. 1
     &               .or. nboundary .eq. 0 .and. i .eq. nnodex
     &               .or. nboundary .eq. 0 .and. j .eq. 1
     &               .or. nboundary .eq. 0 .and. j .eq. nnodey) then

c                      if(psi(i,j) .gt. psilim .and. nboundary .eq. 1
c     &               .or. nboundary .eq. 0 .and. i .le. 2
c     &               .or. nboundary .eq. 0 .and. i .ge. nnodexm1
c     &               .or. nboundary .eq. 0 .and. j .le. 2
c     &               .or. nboundary .eq. 0 .and. j .ge. nnodeym1) then

                     fdk = cexpkxky
                     fek = 0.0
                     ffk = 0.0

                     fgk = 0.0
                     fak = cexpkxky
                     fpk = 0.0

                     frk = 0.0
                     fqk = 0.0
                     fsk = cexpkxky

                     xb(i,j) = 0.0
                     xc(i,j) = 0.0
                     xd(i,j) = 0.0
                  end if  !  end if matal boundary edge region


                  if(i .eq. idiag .and. j .eq. jdiag)then
                     fdksav(n, m) = fdk
                     feksav(n, m) = fek
                     ffksav(n, m) = ffk

                     fgksav(n, m) = fgk
                     faksav(n, m) = fak
                     fpksav(n, m) = fpk

                     frksav(n, m) = frk
                     fqksav(n, m) = fqk
                     fsksav(n, m) = fsk
                  end if   !  if i .eq. idiag .and. j .eq. jdiag



*                 ------------------------------------------
*                 Load the collocation solution matrix, amat
*                 ------------------------------------------
                  ipos = lrindx + (lcindx - 1) * desc_amat(LLD_)
                  p_amat(ipos)     = fdk
                  p_amat(ipos + 1) = fgk
                  p_amat(ipos + 2) = frk

                  if (icnc1 .eq. 1) then
                     p_brhs(ipos)     = xb(i,j)
                     p_brhs(ipos + 1) = xc(i,j)
                     p_brhs(ipos + 2) = xd(i,j)
                  end if

                  ipos = lrindx + (lcindx + 1 - 1) * desc_amat(LLD_)
                  p_amat(ipos)     = fek
                  p_amat(ipos + 1) = fak
                  p_amat(ipos + 2) = fqk

                  ipos = lrindx + (lcindx + 2 - 1) * desc_amat(LLD_)
                  p_amat(ipos)     = ffk
                  p_amat(ipos + 1) = fpk
                  p_amat(ipos + 2) = fsk

                  end if  !   if ismine1

               end do
            end do


         end do
      end do



      time = second1(dummy) - t1
      tmin = time / 60.
      if (myid.eq.0) then
         write(6 ,835) tmin
         write(15,835) tmin
c        call flush
      endif

  835 format('time to load matrix =',f9.3," min")

      call blacs_barrier(icontxt, 'All')


      do n = nkx1, nkx2
         do m = nky1, nky2
            fdksav2d(n - nkx1 + 1, m - nky1 + 1) = fdksav(n,m)
            feksav2d(n - nkx1 + 1, m - nky1 + 1) = feksav(n,m)
            ffksav2d(n - nkx1 + 1, m - nky1 + 1) = ffksav(n,m)

            fgksav2d(n - nkx1 + 1, m - nky1 + 1) = fgksav(n,m)
            faksav2d(n - nkx1 + 1, m - nky1 + 1) = faksav(n,m)
            fpksav2d(n - nkx1 + 1, m - nky1 + 1) = fpksav(n,m)

            frksav2d(n - nkx1 + 1, m - nky1 + 1) = frksav(n,m)
            fqksav2d(n - nkx1 + 1, m - nky1 + 1) = fqksav(n,m)
            fsksav2d(n - nkx1 + 1, m - nky1 + 1) = fsksav(n,m)
         end do
      end do


      do i = 1, nnodex
         do n = nkx1, nkx2
            capd2d(i, n - nkx1 + 1) = capd(i, n)
            xkb2d(i, n - nkx1 + 1) = xkb(i, n)
         end do
      end do


      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fdksav2d,
     &   nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, feksav2d,
     &   nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, ffksav2d,
     &   nxmx, -1, -1)

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fgksav2d,
     &   nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, faksav2d,
     &   nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fpksav2d,
     &   nxmx, -1, -1)

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, frksav2d,
     &   nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fqksav2d,
     &   nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fsksav2d,
     &   nxmx, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex, nmodesmax, capd2d,
     &   nxmx, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex, nmodesmax, xkb2d,
     &   nxmx, -1, -1)


      do n = nkx1, nkx2
         do m = nky1, nky2
            fdksav(n,m) = fdksav2d(n - nkx1 + 1, m - nky1 + 1)
            feksav(n,m) = feksav2d(n - nkx1 + 1, m - nky1 + 1)
            ffksav(n,m) = ffksav2d(n - nkx1 + 1, m - nky1 + 1)

            fgksav(n,m) = fgksav2d(n - nkx1 + 1, m - nky1 + 1)
            faksav(n,m) = faksav2d(n - nkx1 + 1, m - nky1 + 1)
            fpksav(n,m) = fpksav2d(n - nkx1 + 1, m - nky1 + 1)

            frksav(n,m) = frksav2d(n - nkx1 + 1, m - nky1 + 1)
            fqksav(n,m) = fqksav2d(n - nkx1 + 1, m - nky1 + 1)
            fsksav(n,m) = fsksav2d(n - nkx1 + 1, m - nky1 + 1)
         end do
      end do



      do i = 1, nnodex
         do n = nkx1, nkx2
            capd(i, n) = capd2d(i, n - nkx1 + 1)
            xkb(i, n) = xkb2d(i, n - nkx1 + 1)
         end do
      end do


*     --------------------------
*     scale each row by its norm
*     --------------------------
      global_col = 1
      incX = desc_amat(M_)

      do global_row = 1, nrow
        call pdznrm2( ncol, norm2, p_amat, global_row, global_col,
     &       desc_amat, incX )

        alpha = 1.0
        if (norm2 .ne. 0.0) alpha = 1.0/norm2

        call pzscal ( ncol, alpha, p_amat, global_row, global_col,
     &       desc_amat, incX )

        call pzscal ( 1, alpha, p_brhs, global_row, 1,
     &       desc_brhs, desc_brhs(M_))


      enddo


      call blacs_barrier(icontxt, 'All')
      t1 = second1(dummy)

*   write p_amat to file for Azzam

      if((write_for_azzam) .and. (nprow*npcol.eq.1))then

        p_amat_copy = p_amat;
        p_brhs_copy = p_brhs;

        write(6,*) 'Starting p_mat file write for Azzam ...'
        write(6,*) 'ncol = ', ncol
        write(6,*) 'nrow = ', nrow
        WRITE(rowfmt,'(A,i10,A)') '(',ncol,
     &     '(1X,"(",e18.10,",",e18.10,")"))'

        OPEN(UNIT=12, FILE="p_mat_for_azzam.txt",
     &       ACTION="write", STATUS="replace", RECL=(40*ncol+10))
        write(12,'(i6)') ncol
        DO i=1,nrow
          WRITE(12,FMT=rowfmt) (p_amat_copy(i+(j-1)*nrow), j=1,ncol)
        END DO
        CLOSE(UNIT=12)
        OPEN(UNIT=12, FILE="p_brhs_for_azzam.txt",
     &       ACTION="write", STATUS="replace", RECL=(40))
        write(12,'(i6)') ncol
        DO i=1,nrow
          WRITE(12,'(1X,"(",e18.10,",",e18.10,")")') p_brhs_copy(i)
        END DO
        CLOSE(UNIT=12)

        write(6,*) 'DONE'
      end if


*       ---------------------
*       Call scalapack solver
*       ---------------------
      call pzgetrf( nrow, ncol, p_amat, 1, 1, desc_amat, p_ipiv, info)
      if (info.ne.0) then
          write(6,*) 'pzgetrf returns info = ', info
          write(15,*) 'pzgetrf returns info = ', info
          stop '** error ** '
      endif

      call blacs_barrier(icontxt, 'All')

*       ------------------------------
*       Optional: get condition number
*       ------------------------------
       lwork = ndfmax
       lrwork = ndfmax
      anorm = pzlange('1', nrow, ncol, p_amat, 1, 1, desc_amat, rwork)
      call pzgecon('1', nrow, p_amat, 1, 1, desc_amat, anorm, rcond,
     &          workn, lwork, rwork, lrwork, info)

      if (info.ne.0) then
         write(6,*) 'pzgecon returns info = ', info
         write(15,*) 'pzgecon returns info = ', info
         stop '** error ** '
      endif

      call pzgetrs( 'notrans', nrow, 1, p_amat, 1, 1, desc_amat,
     &                p_ipiv, p_brhs, 1,1,desc_brhs, info )

      if (info.ne.0) then
         write(6,*) 'pzgetrs returns info = ', info
         write(15,*) 'pzgetrs returns info = ', info
         stop '** error ** '
      endif

      call blacs_barrier(icontxt, 'All')
      time = second1(dummy) - t1
      time = max(1.0, time)

      if((write_for_azzam) .and. (nprow*npcol.eq.1))then
        write(6,*) 'Starting solution file write for Azzam ...'
        OPEN(UNIT=12, FILE="p_solution_for_azzam.txt",
     &       ACTION="write", STATUS="replace", RECL=(40))
        write(12,'(i6)') ncol
        DO i=1,nrow
          WRITE(12,'(1X,"(",e18.10,",",e18.10,")")') p_brhs(i)
        END DO
        CLOSE(UNIT=12)
        write(6,*) 'DONE'
      end if

c--   Operations for complex matrix factor and solve:
      ops = 8. / 3. * (real(nrow))**3 + 7. * (real(nrow))**2

      gflops  = ops / time / 1.0e+09
      gflopsp = gflops / nproc
      tmin = time / 60.
      if (myid.eq.0) then
        write(6 ,839) nrow, nproc
        write(15,839) nrow, nproc
        write(6 ,833) tmin
        write(15,833) tmin
        write(6 ,837) gflops
        write(15,837) gflops
        write(6 ,838) gflopsp
        write(15,838) gflopsp

        write(6, *) "workn(1) = ", workn(1)
        write(6, *) "rwork(1) = ", rwork(1)
        write(6, *) "rcond = ", rcond

        write(15, *) "workn(1) = ", workn(1)
        write(15, *) "rwork(1) = ", rwork(1)
        write(15, *) "rcond = ", rcond
      endif

  833 format('time taken by ScaLAPACK =',f9.3," min")
  834 format('info =', i10)
  839 format('nrow =', i10, 5x, 'nproc =', i10)
  837 format('operations per sec by pzgetrf & pzgetrs ='
     1   ,f11.3," Gflops/sec")
  838 format('operations per sec per processor ='
     1   ,f9.3," Gflops/sec/processor")

      call blacs_barrier(icontxt, 'All')


*     ------------------
*     broadcast solution
*     ------------------

      brhs2(:) = 0.0
      brhsk(:) = 0.0
      brhs_tmp(:) = 0.0


      icnc = 1
      do irnc=1,nrow
         call infog2l( irnc, icnc, desc_brhs, nprow,npcol,
     &                    myrow,mycol, lrindx, lcindx, rsrc,csrc )
         ismine = (rsrc .eq. myrow) .and. (csrc .eq. mycol)
         if (ismine) then
             ipos = lrindx + (lcindx-1)*desc_amat(LLD_)
             brhs_tmp(irnc) =  p_brhs(ipos)
         endif
      enddo
      call zgsum2d(icontxt, 'All', ' ', nrow, 1, brhs_tmp, nrow, -1, -1)

      brhs2 = brhs_tmp

      if((write_for_azzam) .and. (nprow*npcol.eq.1))then
        write(6,*) 'Starting solution file write for Azzam ...'
        OPEN(UNIT=12, FILE="p_solution2_for_azzam.txt",
     &       ACTION="write", STATUS="replace", RECL=(40))
        write(12,'(i6)') ncol
        DO i=1,nrow
          WRITE(12,'(1X,"(",e18.10,",",e18.10,")")') brhs2(i)
        END DO
        CLOSE(UNIT=12)
        write(6,*) 'DONE'
      end if

      if(write_for_azzam)then
*       Examine solve error
        brhs_tmp = 0
        do i=1,nrow
          do j=1,ncol
            brhs_tmp(i) = brhs_tmp(i) +
     &         p_amat_copy(i+(j-1)*nrow) * p_brhs(j)
          end do
        end do

        l2norm_res = 0;
        l2norm_rhs = 0;
        do i=1,nrow
          l2norm_res = l2norm_res + abs(p_brhs_copy(i)-brhs_tmp(i))**2
          l2norm_rhs = l2norm_rhs + abs(p_brhs_copy(i))**2
*          write(*,*) 'error : ', p_brhs_copy(i)-brhs_tmp(i)
        enddo
        l2norm_res = sqrt(l2norm_res)
        write(*,*) 'l2norm_res : ', l2norm_res
        write(*,*) 'l2norm_rhs : ', l2norm_rhs

        deallocate(p_amat_copy)
        deallocate(p_brhs_copy)
      end if

      deallocate(p_amat)
c      deallocate(Btmp)

      do n = nkx1, nkx2
         do m = nky1, nky2

            irnc = (m - nky1) * 3 + 1 + (n - nkx1) * 3 * nmodesy
            ealphak(n, m) = brhs2(irnc)

            irnc = (m - nky1) * 3 + 2 + (n - nkx1) * 3 * nmodesy
            ebetak(n, m) = brhs2(irnc)

            irnc = (m - nky1) * 3 + 3 + (n - nkx1) * 3 * nmodesy
            ebk(n, m) = brhs2(irnc)

         end do
      end do


      end if

 6000 continue

      call blacs_barrier(icontxt, 'All')




!     -----------------------------------
!     write out solution in Fourier modes
!     -----------------------------------

      if (myid .eq.0) then
         write(54, 309) nkx1, nkx2, nky1, nky2
         write(54, 8310)((ealphak(n, m), n = nkx1, nkx2), m = nky1,nky2)
         write(54, 8310)((ebetak(n, m), n = nkx1, nkx2), m = nky1, nky2)
         write(54, 8310)((ebk(n, m), n = nkx1, nkx2), m = nky1, nky2)
         write(54, 8310)(xkxsav(n), n = nkx1, nkx2)
         write(54, 8310)(xkysav(m), m = nky1, nky2)

         close(54)

      end if





      call checkx(myid, idiag, jdiag, nkx1, nkx2, nky1, nky2,
     1   fdksav, feksav, ffksav, ealphak, ebetak, ebk,
     1   xb, xc, xd, nkdim1, nkdim2, mkdim1, mkdim2,
     1   nnodex, nnodey, nxmx, nymx, xprime, yprime,
     1   idiag, jdiag, 1)

      call checky(myid, idiag, jdiag, nkx1, nkx2, nky1, nky2,
     1   fgksav, faksav, fpksav, ealphak, ebetak, ebk,
     1   xb, xc, xd, nkdim1, nkdim2, mkdim1, mkdim2,
     1   nnodex, nnodey, nxmx, nymx, xprime, yprime,
     1   idiag, jdiag, 1)

      call checkz(myid, idiag, jdiag, nkx1, nkx2, nky1, nky2,
     1   frksav, fqksav, fsksav, ealphak, ebetak, ebk,
     1   xb, xc, xd, nkdim1, nkdim2, mkdim1, mkdim2,
     1   nnodex, nnodey, nxmx, nymx, xprime, yprime,
     1   idiag, jdiag, 1)



      if (anti_alias .eqv. .true.) then
*        --------------------------------------
*        Anti-aliasing filter (two-thirds rule):
*        --------------------------------------
         n_upper = 2. / 3. * nkx2
         n_lower = 2. / 3. * nkx1
         m_upper = 2. / 3. * nky2
         m_lower = 2. / 3. * nky1

         fnmax = 2. / 3. * real(nkx2)
         fmmax = 2. / 3. * real(nky2)


*        ----------------------------------------------
*        Anti-aliasing filter on R only (one-half rule):
*        ----------------------------------------------
c         n_upper = 1. / 2. * nkx2
c         n_lower = 1. / 2. * nkx1
c         m_upper = nky2
c         m_lower = nky1

c         fnmax = 1. / 2. * real(nkx2)
c         fmmax = real(nky2)


*        -----------------------------------------------
*        Anti-aliasing filter on R and Z (one-half rule):
*        -----------------------------------------------
c         n_upper = 1. / 2. * nkx2
c         n_lower = 1. / 2. * nkx1
c         m_upper = 1. / 2. * nky2
c         m_lower = 1. / 2. * nky1

c         fnmax = 1. / 2. * real(nkx2)
c         fmmax = 1. / 2. * real(nky2)


         if (myid .eq. 0) then
         write(6, *) "n_lower = ", n_lower, "   n_upper = ", n_upper
         write(6, *) "m_lower = ", m_lower, "   m_upper = ", m_upper
         write(6, *) "fnmax = ", fnmax, "   fmmax = ", fmmax
         end if

         do n = nkx1, nkx2
            do m = nky1, nky2

               nperp = sqrt((real(n) / fnmax)**2 + (real(m) / fmmax)**2)

c              if (nperp .gt. 1.0) then
               if (n .gt. n_upper .or. n .lt. n_lower .or.
     &            m .gt. m_upper .or. m .lt. m_lower)then

                  ealphak(n, m) = 0.0
                  ebetak(n, m) = 0.0
                  ebk(n, m) = 0.0
               end if

            end do
         end do

      end if


      t1 = second1(dummy)

      ealpha = 0.0
      ebeta  = 0.0
      eb     = 0.0
*     ---------------------------------------
*     Invert Fourier transform to real :: space:
*     ---------------------------------------

c      call sftinv2d_old(xprime, yprime, xkxsav, xkysav, ealpha, ealphak,
c     &   nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2,
c     &   nnodex, nnodey, nkx1, nkx2, nky1, nky2, xx, yy,
c     &   myid, nproc, icontxt)

c      call sftinv2d_old(xprime, yprime, xkxsav, xkysav, ebeta, ebetak,
c     &   nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2,
c     &   nnodex, nnodey, nkx1, nkx2, nky1, nky2, xx, yy,
c     &   myid, nproc, icontxt)

c      call sftinv2d_old(xprime, yprime, xkxsav, xkysav, eb, ebk,
c     &   nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2,
c     &   nnodex, nnodey, nkx1, nkx2, nky1, nky2, xx, yy,
c     &   myid, nproc, icontxt)




      call sftinv_filter(root, xprime, yprime, xkxsav, xkysav,
     &   ealpha, ealphakx, ealphaky,
     &   ealphak, nxmx, nymx,
     &   nkdim1, nkdim2, mkdim1, mkdim2,
     &   nnodex, nnodey, nkx1, nkx2, nky1, nky2, xx, yy,
     &   myid, nproc, icontxt,
     &   uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz,
     &   z1, xkte, xkti, xmi1, omgrf, omgci1, omgce, xkphi,
     &   omgpe2, omgp12, xme, upshift, nphi, rt, 0, capr, y, jmid,
     &   n_ono, m_ono, xkperp2_slow, xkperp2_fast, xkprl_a, P_a,
     &   xk_cutoff)

      call sftinv_filter(root, xprime, yprime, xkxsav, xkysav,
     &   ebeta, ebetakx, ebetaky,
     &   ebetak, nxmx, nymx,
     &   nkdim1, nkdim2, mkdim1, mkdim2,
     &   nnodex, nnodey, nkx1, nkx2, nky1, nky2, xx, yy,
     &   myid, nproc, icontxt,
     &   uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz,
     &   z1, xkte, xkti, xmi1, omgrf, omgci1, omgce, xkphi,
     &   omgpe2, omgp12, xme, upshift, nphi, rt, 0, capr, y, jmid,
     &   n_ono, m_ono, xkperp2_slow, xkperp2_fast, xkprl_a, P_a,
     &   xk_cutoff)

      call sftinv_filter(root, xprime, yprime, xkxsav, xkysav,
     &   eb, ebkx, ebky,
     &   ebk, nxmx, nymx,
     &   nkdim1, nkdim2, mkdim1, mkdim2,
     &   nnodex, nnodey, nkx1, nkx2, nky1, nky2, xx, yy,
     &   myid, nproc, icontxt,
     &   uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz,
     &   z1, xkte, xkti, xmi1, omgrf, omgci1, omgce, xkphi,
     &   omgpe2, omgp12, xme, upshift, nphi, rt, 1, capr, y, jmid,
     &   n_ono, m_ono, xkperp2_slow, xkperp2_fast, xkprl_a, P_a,
     &   xk_cutoff)




      do i = 1, nnodex
         kimag_slow(i) = aimag(csqrt(xkperp2_slow(i,jmid)))
         kimag_fast(i) = aimag(csqrt(xkperp2_fast(i,jmid)))
      end do

      call dgrate(xprime, kimag_slow, 1, nnodex, ans, nxmx)
      spa_slow = 1.0 - exp(-2.0 * ans)

      call dgrate(xprime, kimag_fast, 1, nnodex, ans, nxmx)
      spa_fast = 1.0 - exp(-2.0 * ans)


      if (myid .eq. 0) then
         write(15, *) "spa_slow = ", spa_slow
         write(6,  *) "spa_slow = ", spa_slow
         write(15, *) "spa_fast = ", spa_fast
         write(6,  *) "spa_fast = ", spa_fast
      end if


      time = second1(dummy) - t1
      tmin = time / 60.
      if (myid.eq.0) then
         write(6, 1837) tmin
         write(15, 1837) tmin
      endif

 1837 format('time taken to invert Fourier transform =',f9.3," min")


!     --------------------------------
!     write out solution in real :: space
!     --------------------------------

      if (myid .eq.0) then

         write(53, 1314) nnodex, nnodey
         write(53, 310) (capr(i), i = 1, nnodex)
         write(53, 310) (y(j), j = 1, nnodey)
         write(53, 8310)((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)
         write(53, 8310)((ebeta(i, j), i = 1, nnodex), j = 1, nnodey)
         write(53, 8310)((eb(i, j), i = 1, nnodex), j = 1, nnodey)

         close (53)
      end if




      do n = nkx1, nkx2
         do m = nky1, nky2

            ealphakmod(n, m) = sqrt(conjg(ealphak(n, m))* ealphak(n, m))
            ebetakmod(n, m) = sqrt(conjg(ebetak(n, m)) * ebetak(n, m))
            ebkmod(n, m) = sqrt(conjg(ebk(n, m)) * ebk(n, m))

         end do
      end do




      if (myid .eq. 0) then
         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 940)
         write(15, 920)

         write(6, 910)
         write(6, 930)
         write(6, 920)
         write(6, 940)
         write(6, 920)


         i = idiag
         do j = 1, nnodey
            write(6,  900)i, j, capr(i), y(j),
     &                                 ealpha(i,j), ebeta(i,j), eb(i,j)
            write(15, 900)i, j, capr(i), y(j),
     &                                 ealpha(i,j), ebeta(i,j), eb(i,j)
         end do

         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 940)
         write(15, 920)

         write(6, 910)
         write(6, 930)
         write(6, 920)
         write(6, 940)
         write(6, 920)

         j = jdiag
         do i = 1, nnodex
            write(6,  900)i, j, capr(i), y(j),
     &                                 ealpha(i,j), ebeta(i,j), eb(i,j)
            write(15, 900)i, j, capr(i), y(j),
     &                                 ealpha(i,j), ebeta(i,j), eb(i,j)
         end do
      endif


  900 format(i5, i5, 1p,9e12.3)
  940 format(" " , "   i", "    j",
     &                    5x, "  R(i) ", 5x, "  Z(j) ",
     &                    3x, "re ealp", 5x, "im ealp",
     &                    5x, "re ebet", 5x, "im ebet",
     &                    5x, " re eb ", 5x, " im eb ")
  910 format ("1")
  920 format ("0")
  930 format(3x, 'rf electric field in Stix frame')


*     ----------------------------------------------
*     Calculate E in the Lab frame and eplus, eminus
*     ----------------------------------------------
      isq2 = SQRT(0.5)
      do i = 1, nnodex
         do j = 1, nnodey

            ex(i,j)   = uxx(i,j) * ealpha(i,j)
     &                + uyx(i,j) * ebeta(i,j)
     &                + uzx(i,j) * eb(i,j)

            ey(i,j)   = uxy(i,j) * ealpha(i,j)
     &                + uyy(i,j) * ebeta(i,j)
     &                + uzy(i,j) * eb(i,j)

            ez(i,j)   = uxz(i,j) * ealpha(i,j)
     &                + uyz(i,j) * ebeta(i,j)
     &                + uzz(i,j) * eb(i,j)

            eplus(i,j)  = isq2 * (ealpha(i,j) + zi * ebeta(i,j))
            eminus(i,j) = isq2 * (ealpha(i,j) - zi * ebeta(i,j))

         end do
      end do

*-------------------------------------------------------
*     Interpolate Eplus, Eminus and mod B onto (rho, theta) grid
*-------------------------------------------------------
      do n = 1, nnoderho
         do m = 1, mnodetheta

            xgiv = capr_flux(n, m) - rt - xwleft
            ygiv = capz_flux(n, m) - ybottom

            call intplt(xgiv, ygiv, fout, nnodex, nnodey, eplus,
     &         nxmx, nymx, dx, dy)
            eplus_flux(n, m) = fout

            call intplt(xgiv, ygiv, fout, nnodex, nnodey, eminus,
     &         nxmx, nymx, dx, dy)
            eminus_flux(n, m) = fout

            call intplt(xgiv, ygiv, fout, nnodex, nnodey, xkperp_cold,
     &         nxmx, nymx, dx, dy)
            xkperp_flux(n, m) = fout

            call intplt_re(xgiv, ygiv, fout, nnodex, nnodey, bmod,
     &         nxmx, nymx, dx, dy)
            bmod_flux(n, m) = fout



c           if(myid .eq. 0)write(6,*)n, m, eplus_flux(n,m)

         end do
      end do

*---------------------------------------------------
*     Interpolate back onto (R, Z) grid for plotting
*---------------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey

            rho_giv = rho(i, j)
            the_giv = theta(i,j)

            call intplt(rho_giv, the_giv, fout, nnoderho, mnodetheta,
     &         eplus_flux, nrhomax, nthetamax, drho, dtheta)
            eplus_flux_plot(i, j) = fout

            call intplt(rho_giv, the_giv, fout, nnoderho, mnodetheta,
     &         eminus_flux, nrhomax, nthetamax, drho, dtheta)
            eminus_flux_plot(i, j) = fout

            call intplt(rho_giv, the_giv, fout, nnoderho, mnodetheta,
     &         xkperp_flux, nrhomax, nthetamax, drho, dtheta)
            xkperp_flux_plot(i, j) = fout

c           if(myid .eq. 0)write(6,*)i, j, eplus_flux_plot(i, j)

         end do
      end do


      if (myid .eq. 0) then
         write(15, 910)
         write(15, 1930)
         write(15, 920)
         write(15, 1940)
         write(15, 920)

         write(6, 910)
         write(6, 1930)
         write(6, 920)
         write(6, 1940)
         write(6, 920)

         i = idiag
         do j = 1, nnodey
            write(6,  900)i, j, capr(i), y(j), ex(i,j), ey(i,j), ez(i,j)
            write(15, 900)i, j, capr(i), y(j), ex(i,j), ey(i,j), ez(i,j)
         end do

         write(15, 910)
         write(15, 1930)
         write(15, 920)
         write(15, 1940)
         write(15, 920)

         write(6, 910)
         write(6, 1930)
         write(6, 920)
         write(6, 1940)
         write(6, 920)

         j = jdiag
         do i = 1, nnodex
            write(6,  900)i, j, capr(i), y(j), ex(i,j), ey(i,j), ez(i,j)
            write(15, 900)i, j, capr(i), y(j), ex(i,j), ey(i,j), ez(i,j)
         end do
      endif


 1940 format(" " , "   i", "    j",
     &                    5x, "  R(i) ", 5x, "  Z(j) ",
     &                    4x, " re ex", 6x, " im ex",
     &                    6x, " re ey", 6x, " im ey",
     &                    6x, " re ez", 6x, " im ez")

 1930 format(3x, 'rf electric field in Lab frame')


*     -----------------------------------
*     Calculate minimum and maximum xkprl
*     -----------------------------------
      xkprl_min = 10000.0
      xkprl_max = -10000.0
      xkprl_pos_min = 10000.0
      xkprl_neg_max = -10000.0

      zetai_max = 0.0
      zetai_min = 0.0

      xjx_max = 0.0
      xjy_max = 0.0
      xjz_max = 0.0

      do i = 1, nnodex
         do j = 1, nnodey

            if(xjx(i,j) .gt. xjx_max) xjx_max = xjx(i,j)
            if(xjy(i,j) .gt. xjy_max) xjy_max = xjy(i,j)
            if(xjz(i,j) .gt. xjz_max) xjz_max = xjz(i,j)


            if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     &                               .or. nboundary .eq. 0)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     xkalp = uxx(i, j) * xkxsav(n)
     &                     + uxy(i, j) * xkysav(m)
     &                     + uxz(i, j) * xkphi(i)
                     xkbet = uyx(i, j) * xkxsav(n)
     &                     + uyy(i, j) * xkysav(m)
     &                     + uyz(i, j) * xkphi(i)
                     xkprl = uzx(i, j) * xkxsav(n)
     &                     + uzy(i, j) * xkysav(m)
     &                     + uzz(i, j) * xkphi(i)
                     xkperp = sqrt(xkalp**2 + xkbet**2)
                     xkprl_orig = xkprl

*                    ------------------------------------
*                    Optional: leave out upshift in xkprl
*                    --------------------------------- --
                     if (upshift .eq. 0) xkprl = uzz(i,j) * xkphi(i)
c                     if (upshift .eq. 0) xkprl = nphi / rt

                     if (upshift .eq. -1)then
                        if (xkperp  .gt. xk_cutoff)
     &                                   xkprl = uzz(i,j) * xkphi(i)
                     end if

                     if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
                     if (xkperp .eq. 0.0) xkperp = 1.0e-08

                     sgn_kprl = sign(1.0, xkprl)
                     akprl = abs(xkprl)


*                   ----------------------------------------------
*                   Optional: Don't allow xkprl to be 0 (upshift = -2)
*                   ----------------------------------------------
                    if (upshift .eq. -2) then
                       if (akprl .lt. akprl_min) then
                          xkprl = akprl_min * sgn_kprl
                       end if
                    end if


                     gammab = abs(omgci2(i,j)/ (2.0 * vthi20 * xkprl**2)
     &                                     * gradprlb(i,j) / bmod(i,j))

*                    ------------------------------------------------
*                    Calculate Brambilla's xkrpl_eff using l = 1 only
*                    ------------------------------------------------
                     y0 = 1.5
                     yb = y0

                     if(sgn_kprl .ge. 0.0)then
                        fgam = 1.0
                        if(gammab .gt. 1.0e-05)then
                           yb = y0
                           fgam = (sqrt(1. +  4. * gammab * yb) - 1.)
     &                        / (2. * gammab * yb)
                        endif
                        xkprl_eff = xkprl / fgam
                     end if

                     if(sgn_kprl .lt. 0.0)then
                        fgam = 1.0
                        if(gammab .gt. 1.0e-05)then
                           descrim = 1. - 4. * gammab * y0
                           if (descrim .ge. 0.0) yb =   y0
                           if (descrim .lt. 0.0) yb = - y0
                           fgam = (1. - sqrt(1. -  4. * gammab * yb) )
     &                          / (2. * gammab * yb)
                        endif
                        xkprl_eff = xkprl / fgam
                     end if


                     if (upshift .ne. 0) xkprl = xkprl_eff


                     if(xkprl .gt. xkprl_max) xkprl_max = xkprl
                     if(xkprl .lt. xkprl_min) xkprl_min = xkprl

                     if(xkprl .gt. 0.0)then
                       if(xkprl .lt. xkprl_pos_min)xkprl_pos_min = xkprl
                     end if

                     if(xkprl .lt. 0.0)then
                       if(xkprl .gt. xkprl_neg_max)xkprl_neg_max = xkprl
                     end if

                     alphae = sqrt(2.0 * xkte(i,j) / xme)
                     zetai = xkprl_orig * alphae / omgrf

                     if(zetai .gt. zetai_max) zetai_max = zetai
                     if(zetai .lt. zetai_min) zetai_min = zetai

                  end do
               end do

            end if

         end do
      end do

      if (myid .eq. 0) then

         write(15, *) "xjx_max = ", xjx_max
         write(15, *) "xjy_max = ", xjy_max
         write(15, *) "xjz_max = ", xjz_max

         write(15, *) "xkprl_max = ", xkprl_max
         write(15, *) "xkprl_pos_min = ", xkprl_pos_min
         write(15, *) "xkprl_neg_max = ", xkprl_neg_max
         write(15, *) "xkprl_min = ", xkprl_min

         write(6, *) "xjx_max = ", xjx_max
         write(6, *) "xjy_max = ", xjy_max
         write(6, *) "xjz_max = ", xjz_max

         write(6, *) "xkprl_max = ", xkprl_max
         write(6, *) "xkprl_pos_min = ", xkprl_pos_min
         write(6, *) "xkprl_neg_max = ", xkprl_neg_max
         write(6, *) "xkprl_min = ", xkprl_min

      end if

*     ---------------------------------
*     set artificial absorption to zero
*     ---------------------------------
      if (use_no_damp .eqv. .true.)then
         delta0 = 0.0
         damping = 0.0
      end if

*     --------------------------------------
*     calculate individual species currents:
*     --------------------------------------

      t1 = second1(dummy)

      call current_orbit(xjpxe, xjpye, xjpze,
     &   nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xme, qe, xnea, xkte, omgce, omgpe2, lmaxe,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomea, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, 0, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfdupere, dfdupare,
     &   UminPara, UmaxPara, UPERP, UPARA,
     &   vce_mks, dfe_cql_uprp, dfe_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   xjpx_ehst, xjpy_ehst, xjpz_ehst, nkperp, zi, eps0,
     &   v0i, xk0, kperp_max, i_sav, j_sav, upshift, damping,
     &   xk_cutoff, odd_order, dx, dy, xkphi,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx, d2rdx2,
     &   dxxuzx, dxxuzy, dxxuzz, dxyuzx, dyyuzx, dxyuzz,
     &   z0_table1, z1_table1, z2_table1, zetai_table, dKdL_table,
     &   dKdL_giv, nmax, mmax, use_new_z2, ntable, mtable, nnodelb,
     &   nlmx, nkldim1, nkldim2, rc, kappa_hatx, kappa_haty,kappa_hatz,
     &   bhatx, bhaty, bhatz, rwleft, ybottom, workl, Eln, El, lb,
     &   lbprime, Elsum, xklsav, use_fourier_z2, lbmax,
     &   qsafety, rwright, ytop, nmodesx, nmodesy, xprime, yprime,
     &   rmaxis, zmaxis, psio, psimax, psi_tor_max, i_psi_eq,
     &   dldb_tot12, dldbavg, nrhomax, n_prof_flux, rhomax, norb_dim,
     &   capr_x,  capz_x,  phin_x, capr_xp, capz_xp, phin_xp,
     &   capr_xm, capz_xm, phin_xm, len_xm, len_xp,
     &   len_x, b_capr, c_capr,  d_capr, b_capz,
     &   c_capz, d_capz, b_phin, c_phin, d_phin,
     &   capr_lb, capz_lb, phin_lb, rhoplasm, norb_dim2, psiorb)



      call ntilda_(ntilda_e, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xme, qe, xnea, xkte, omgce, omgpe2, lmaxe,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomea, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, 0, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfdupere, dfdupare,
     &   UminPara, UmaxPara, UPERP, UPARA,
     &   vce_mks, dfe_cql_uprp, dfe_cql_uprl, rho, rho_a, nbessj,
     &   nkperp,
     &   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift,
     &   damping, xk_cutoff, rt, xkphi,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,
     &   z0_table1, z1_table1, z2_table1, zetai_table, dKdL_table,
     &   dKdL_giv, nmax, mmax, use_new_z2, ntable, mtable)



      call current_1(xjpx1, xjpy1, xjpz1, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi1, qi1, xn1a, xkti, omgci1, omgp12, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti1, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper1, dfdupar1,
     &   UminPara, UmaxPara, UPERP, UPARA,
     &   vc1_mks, df1_cql_uprp, df1_cql_uprl, rho, rho_a, nbessj,nkperp,
     &   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav1, j_sav1, upshift,
     &   damping, xk_cutoff, rt, odd_order, dx, dy, xkphi,
     &   z0_table1, z1_table1,z2_table1, zetai_table, dKdL_table,
     &   dKdL_giv_1, nmax, mmax, use_new_z2, ntable, mtable)



      if(eta2 .ne. 0.0) then

         call current_2(xjpx2, xjpy2, xjpz2, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi2, qi2, xn2a, xkti2, omgci2, omgp22, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti2, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper2, dfdupar2,
     &   UminPara, UmaxPara, UPERP, UPARA,
     &   vc2_mks, df2_cql_uprp, df2_cql_uprl, rho, rho_a, nbessj,nkperp,
     &   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav2, j_sav2, upshift,
     &   damping, xk_cutoff, rt, odd_order, dx, dy, xkphi)

      endif

      if(eta3 .ne. 0.0) then

         call current(xjpx3, xjpy3, xjpz3, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi3, qi3, xn3a, xkti3, omgci3, omgp32, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti3, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper3, dfdupar3,
     &   UminPara, UmaxPara, UPERP, UPARA,
     &   vc3_mks, df3_cql_uprp, df3_cql_uprl, rho, rho_a, nbessj,nkperp,
     &   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift,
     &   damping, xk_cutoff, rt, odd_order, dx, dy, xkphi)

      endif

      if(eta4 .ne. 0.0) then

         call current(xjpx4, xjpy4, xjpz4, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi4, qi4, xn4a, xkti4, omgci4, omgp42, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti4, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper4, dfdupar4,
     &   UminPara, UmaxPara, UPERP, UPARA,
     &   vc4_mks, df4_cql_uprp, df4_cql_uprl, rho, rho_a, nbessj,nkperp,
     &   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift,
     &   damping, xk_cutoff, rt, odd_order, dx, dy, xkphi)

      endif

      if(eta5 .ne. 0.0) then

         call current(xjpx5, xjpy5, xjpz5, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi5, qi5, xn5a, xkti5, omgci5, omgp52, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti5, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper5, dfdupar5,
     &   UminPara, UmaxPara, UPERP, UPARA,
     &   vc5_mks, df5_cql_uprp, df5_cql_uprl, rho, rho_a, nbessj,nkperp,
     &   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift,
     &   damping, xk_cutoff, rt, odd_order, dx, dy, xkphi)

      endif

      if(eta6 .ne. 0.0) then

         call current(xjpx6, xjpy6, xjpz6, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi6, qi6, xn6a, xkti6, omgci6, omgp62, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti6, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper6, dfdupar6,
     &   UminPara, UmaxPara, UPERP, UPARA,
     &   vc6_mks, df6_cql_uprp, df6_cql_uprl, rho, rho_a, nbessj,nkperp,
     &   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift,
     &   damping, xk_cutoff, rt, odd_order, dx, dy, xkphi)

      endif


      if(eta_slo .ne. 0.0)
     &   call cur_slo(xj_slox, xj_sloy, xj_sloz, nxmx, nymx,
     &   nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi_slo, qi_slo, xna_slo, eslow, omgci_slo, omgp2_slo, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   xkte, zeffcd, myid, nproc)



      time = second1(dummy) - t1
      tmin = time / 60.
      if (myid .eq. 0) then
         write(6 ,836) tmin
         write(15,836) tmin
      endif

  836 format('time to calculate currents =',f9.3," min")


*     --------------------------------------------
*     calculate total plasma current in Stix frame:
*     --------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey
            xjpx(i,j) = xjpxe(i,j) + xjpx1(i,j) + xjpx2(i,j)+ xjpx3(i,j)
     &                + xjpx4(i,j) + xjpx5(i,j) + xjpx6(i,j)
     &                + xj_slox(i, j)
            xjpy(i,j) = xjpye(i,j) + xjpy1(i,j) + xjpy2(i,j)+ xjpy3(i,j)
     &                + xjpy4(i,j) + xjpy5(i,j) + xjpy6(i,j)
     &                + xj_sloy(i, j)
            xjpz(i,j) = xjpze(i,j) + xjpz1(i,j) + xjpz2(i,j)+ xjpz3(i,j)
     &                + xjpz4(i,j) + xjpz5(i,j) + xjpz6(i,j)
     &                + xj_sloz(i, j)

         end do
      end do

*     -------------------------------------------
*     Calculate total plasma current in Lab frame
*     -------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey

            xjpx_lab(i,j)   = uxx(i,j) * xjpx(i,j)
     &                      + uyx(i,j) * xjpy(i,j)
     &                      + uzx(i,j) * xjpz(i,j)

            xjpy_lab(i,j)   = uxy(i,j) * xjpx(i,j)
     &                      + uyy(i,j) * xjpy(i,j)
     &                      + uzy(i,j) * xjpz(i,j)

            xjpz_lab(i,j)   = uxz(i,j) * xjpx(i,j)
     &                      + uyz(i,j) * xjpy(i,j)
     &                      + uzz(i,j) * xjpz(i,j)


            xjpxe_lab(i,j)   = uxx(i,j) * xjpxe(i,j)
     &                       + uyx(i,j) * xjpye(i,j)
     &                       + uzx(i,j) * xjpze(i,j)

            xjpye_lab(i,j)   = uxy(i,j) * xjpxe(i,j)
     &                       + uyy(i,j) * xjpye(i,j)
     &                       + uzy(i,j) * xjpze(i,j)

            xjpze_lab(i,j)   = uxz(i,j) * xjpxe(i,j)
     &                       + uyz(i,j) * xjpye(i,j)
     &                       + uzz(i,j) * xjpze(i,j)


         end do
      end do


*     ----------------------------------------------------------------
*     Calculate fluctuating charge density (rho_plasma) for electrons:
*     ----------------------------------------------------------------
      do i = 2, nnodex - 1
         do j = 2, nnodey - 1

            djxdx = (xjpxe_lab(i+1, j) - xjpxe_lab(i-1, j)) / (2.0 * dx)
            djydy = (xjpye_lab(i, j+1) - xjpye_lab(i, j-1)) / (2.0 * dy)

            rho_pla(i,j) = (djxdx + djydy + xjpxe_lab(i,j) / capr(i)
     &                     + zi * nphi / capr(i) * xjpze_lab(i,j) )
     &                  / (zi * omgrf)

c           ntilda_e(i,j) = rho_pla(i,j) / qe

         end do
      end do

*     -------------------------------------------------
*     Calculate charge density (rho_ant) on the antenna:
*     -------------------------------------------------
      do i = 2, nnodex - 1
         do j = 2, nnodey - 1

            djxdx = (xjx(i+1, j) - xjx(i-1, j)) / (2.0 * dx)
            djydy = (xjy(i, j+1) - xjy(i, j-1)) / (2.0 * dy)

            rho_ant(i,j) = (djxdx + djydy + xjx(i,j)/capr(i)
     &                     + zi * nphi / capr(i) * xjz(i,j) )
     &                  / (zi * omgrf)

         end do
      end do


      t1 = second1(dummy)

*     --------------------------------
*     Distribution function from CQL3D
*     --------------------------------

!     -----------------
!     allocate arrays
!     -----------------

      allocate( bqlavg_e(nuper, nupar, nnoderho) )
      allocate( cqlavg_e(nuper, nupar, nnoderho) )
      allocate( eqlavg_e(nuper, nupar, nnoderho) )
      allocate( fqlavg_e(nuper, nupar, nnoderho) )
      bqlavg_e = 0.0
      cqlavg_e = 0.0
      eqlavg_e = 0.0
      fqlavg_e = 0.0


      allocate( bqlavg_i1(nuper, nupar, nnoderho) )
      allocate( cqlavg_i1(nuper, nupar, nnoderho) )
      allocate( eqlavg_i1(nuper, nupar, nnoderho) )
      allocate( fqlavg_i1(nuper, nupar, nnoderho) )
      bqlavg_i1 = 0.0
      cqlavg_i1 = 0.0
      eqlavg_i1 = 0.0
      fqlavg_i1 = 0.0

      allocate( bqlavg_work(nuper, nupar, nnoderho) )
      allocate( cqlavg_work(nuper, nupar, nnoderho) )
      allocate( eqlavg_work(nuper, nupar, nnoderho) )
      allocate( fqlavg_work(nuper, nupar, nnoderho) )
      bqlavg_work = 0.0
      cqlavg_work = 0.0
      eqlavg_work = 0.0
      fqlavg_work = 0.0


      if(eta2 .ne. 0.0)then
         allocate( bqlavg_i2(nuper, nupar, nnoderho) )
         allocate( cqlavg_i2(nuper, nupar, nnoderho) )
         allocate( eqlavg_i2(nuper, nupar, nnoderho) )
         allocate( fqlavg_i2(nuper, nupar, nnoderho) )
         bqlavg_i2 = 0.0
         cqlavg_i2 = 0.0
         eqlavg_i2 = 0.0
         fqlavg_i2 = 0.0
      end if

      if(eta3 .ne. 0.0)then
         allocate( bqlavg_i3(nuper, nupar, nnoderho) )
         allocate( cqlavg_i3(nuper, nupar, nnoderho) )
         allocate( eqlavg_i3(nuper, nupar, nnoderho) )
         allocate( fqlavg_i3(nuper, nupar, nnoderho) )
         bqlavg_i3 = 0.0
         cqlavg_i3 = 0.0
         eqlavg_i3 = 0.0
         fqlavg_i3 = 0.0
      end if

      if(eta4 .ne. 0.0)then
         allocate( bqlavg_i4(nuper, nupar, nnoderho) )
         allocate( cqlavg_i4(nuper, nupar, nnoderho) )
         allocate( eqlavg_i4(nuper, nupar, nnoderho) )
         allocate( fqlavg_i4(nuper, nupar, nnoderho) )
         bqlavg_i4 = 0.0
         cqlavg_i4 = 0.0
         eqlavg_i4 = 0.0
         fqlavg_i4 = 0.0
      end if

      if(eta5 .ne. 0.0)then
         allocate( bqlavg_i5(nuper, nupar, nnoderho) )
         allocate( cqlavg_i5(nuper, nupar, nnoderho) )
         allocate( eqlavg_i5(nuper, nupar, nnoderho) )
         allocate( fqlavg_i5(nuper, nupar, nnoderho) )
         bqlavg_i5 = 0.0
         cqlavg_i5 = 0.0
         eqlavg_i5 = 0.0
         fqlavg_i5 = 0.0
      end if

      if(eta6 .ne. 0.0)then
         allocate( bqlavg_i6(nuper, nupar, nnoderho) )
         allocate( cqlavg_i6(nuper, nupar, nnoderho) )
         allocate( eqlavg_i6(nuper, nupar, nnoderho) )
         allocate( fqlavg_i6(nuper, nupar, nnoderho) )
         bqlavg_i6 = 0.0
         cqlavg_i6 = 0.0
         eqlavg_i6 = 0.0
         fqlavg_i6 = 0.0
      end if



      if (nzeta_wdote .gt. 0) then

      ! ------------------ !
      ! --  Electrons   -- !
      ! ------------------ !

      t1e = second1(dummy)


      call ql_myra_write(bqlavg_e, cqlavg_e, eqlavg_e, fqlavg_e,
     &   wdote, fx0e, fy0e, fz0e, darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xme, qe, xnea, xkte, omgce, omgpe2, lmaxe,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomea, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, 0, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfdupere, dfdupare, fperpe,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vce_mks, dfe_cql_uprp, dfe_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdote, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff)


      if (use_new_wdot .eqv. .true.) then
         call wdot_new_maxwellian(z2_electron,
     &   wdote, darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xme, qe, xnea, xkte, omgce, omgpe2, lmaxe,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomea, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, 0, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfdupere, dfdupare, fperpe,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vce_mks, dfe_cql_uprp, dfe_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdote, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff,  z0_table1, z1_table1,
     &   z2_table1, ntable, mtable, zetai_table, dKdL_table,
     &   nmax, mmax)
      end if


      tmin = (second1(dummy) - t1e) / 60.
      if (myid.eq.0) then
         write(6 , 2870) tmin
         write(15, 2870) tmin
      endif
 2870 format('time for electron QL operator =', f9.3, ' min')

      endif

      if (nzeta_wdoti .gt. 0) then


      ! ------------------ !
      ! --Majority ions -- !
      ! ------------------ !
      t1i1 = second1(dummy)

      call ql_myra_write(bqlavg_i1, cqlavg_i1, eqlavg_i1, fqlavg_i1,
     &   wdoti1, fx0i1, fy0i1, fz0i1, darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi1, qi1, xn1a, xkti, omgci1, omgp12, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti1, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper1, dfdupar1, fperp1,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc1_mks, df1_cql_uprp, df1_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff)


      if (use_new_wdot .eqv. .true.) then
         call wdot_new_maxwellian(z2_electron,
     &   wdoti1,
     &   darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi1, qi1, xn1a, xkti, omgci1, omgp12, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti1, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper1, dfdupar1, fperp1,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc1_mks, df1_cql_uprp, df1_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff,  z0_table1, z1_table1,
     &   z2_table1, ntable, mtable, zetai_table, dKdL_table,
     &   nmax, mmax)
      end if

      tmin = (second1(dummy) - t1i1) / 60.
      if (myid.eq.0) then
         write(6 , 2871) tmin
         write(15, 2871) tmin
      endif
 2871 format('time for majority ion QL operator =', f9.3, ' min')


      ! ------------------ !
      ! --Minority ions -- !
      ! ------------------ !

      if(eta2 .ne. 0.0)then
      t1i2 = second1(dummy)
      call ql_myra_write(bqlavg_i2, cqlavg_i2, eqlavg_i2, fqlavg_i2,
     &   wdoti2, fx0i2, fy0i2, fz0i2, darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi2, qi2, xn2a, xkti2, omgci2, omgp22, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti2, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper2, dfdupar2, fperp2,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc2_mks, df2_cql_uprp, df2_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff)

      if (use_new_wdot .eqv. .true.) then
         call wdot_new_maxwellian(z2_electron,
     &   wdoti2,
     &   darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi2, qi2, xn2a, xkti2, omgci2, omgp22, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti2, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper2, dfdupar2, fperp2,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc2_mks, df2_cql_uprp, df2_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff,  z0_table1, z1_table1,
     &   z2_table1, ntable, mtable, zetai_table, dKdL_table,
     &   nmax, mmax)
      end if

      tmin = (second1(dummy) - t1i2) / 60.
      if (myid.eq.0) then
         write(6 , 2872) tmin
         write(15, 2872) tmin
      endif
 2872 format('time for minority ion QL operator =', f9.3, ' min')

      end if

      ! --------------------- !
      ! --Third ion species-- !
      ! --------------------- !

      if(eta3 .ne. 0.0)then
      t1i3 = second1(dummy)
      call ql_myra_write(bqlavg_i3, cqlavg_i3, eqlavg_i3, fqlavg_i3,
     &   wdoti3, fx0i3, fy0i3, fz0i3, darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi3, qi3, xn3a, xkti3, omgci3, omgp32, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti3, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper3, dfdupar3, fperp3,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc3_mks, df3_cql_uprp, df3_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff)

      if (use_new_wdot .eqv. .true.) then
         call wdot_new_maxwellian(z2_electron,
     &   wdoti3,
     &   darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi3, qi3, xn3a, xkti3, omgci3, omgp32, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti3, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper3, dfdupar3, fperp3,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc3_mks, df3_cql_uprp, df3_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff,  z0_table1, z1_table1,
     &   z2_table1, ntable, mtable, zetai_table, dKdL_table,
     &   nmax, mmax)
      end if

      tmin = (second1(dummy) - t1i3) / 60.
      if (myid.eq.0) then
         write(6 , 2873) tmin
         write(15, 2873) tmin
      endif
 2873 format('time for 3rd ion species QL operator =', f9.3, ' min')
      end if




      ! ---------------------- !
      ! --Fourth ion species-- !
      ! ---------------------- !

      if(eta4 .ne. 0.0)then
      t1i4 = second1(dummy)
      call ql_myra_write(bqlavg_i4, cqlavg_i4, eqlavg_i4, fqlavg_i4,
     &   wdoti4, fx0i4, fy0i4, fz0i4, darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi4, qi4, xn4a, xkti4, omgci4, omgp42, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti4, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper4, dfdupar4, fperp4,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc4_mks, df4_cql_uprp, df4_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff)

      if (use_new_wdot .eqv. .true.) then
         call wdot_new_maxwellian(z2_electron,
     &   wdoti4,
     &   darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi4, qi4, xn4a, xkti4, omgci4, omgp42, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti4, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper4, dfdupar4, fperp4,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc4_mks, df4_cql_uprp, df4_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff,  z0_table1, z1_table1,
     &   z2_table1, ntable, mtable, zetai_table, dKdL_table,
     &   nmax, mmax)
      end if

      tmin = (second1(dummy) - t1i4) / 60.
      if (myid.eq.0) then
         write(6 , 2874) tmin
         write(15, 2874) tmin
      endif
 2874 format('time for 4th ion species QL operator =', f9.3, ' min')
      end if


      ! --------------------- !
      ! --Fifth ion species-- !
      ! --------------------- !

      if(eta5 .ne. 0.0)then
      t1i5 = second1(dummy)
      call ql_myra_write(bqlavg_i5, cqlavg_i5, eqlavg_i5, fqlavg_i5,
     &   wdoti5, fx0i5, fy0i5, fz0i5, darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi5, qi5, xn5a, xkti5, omgci5, omgp52, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti5, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper5, dfdupar5, fperp5,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc5_mks, df5_cql_uprp, df5_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff)

      if (use_new_wdot .eqv. .true.) then
         call wdot_new_maxwellian(z2_electron,
     &   wdoti5,
     &   darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi5, qi5, xn5a, xkti5, omgci5, omgp52, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti5, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper5, dfdupar5, fperp5,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc5_mks, df5_cql_uprp, df5_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff,  z0_table1, z1_table1,
     &   z2_table1, ntable, mtable, zetai_table, dKdL_table,
     &   nmax, mmax)
      end if

      tmin = (second1(dummy) - t1i5) / 60.
      if (myid.eq.0) then
         write(6 , 2875) tmin
         write(15, 2875) tmin
      endif
 2875 format('time for 5th ion species QL operator =', f9.3, ' min')

      end if


      ! --------------------- !
      ! --Sixth ion species-- !
      ! --------------------- !

      if(eta6 .ne. 0.0)then
      t1i6 = second1(dummy)
      call ql_myra_write(bqlavg_i6, cqlavg_i6, eqlavg_i6, fqlavg_i6,
     &   wdoti6, fx0i6, fy0i6, fz0i6, darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi6, qi6, xn6a, xkti6, omgci6, omgp62, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti6, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper6, dfdupar6, fperp6,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc6_mks, df6_cql_uprp, df6_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff)

      if (use_new_wdot .eqv. .true.) then
         call wdot_new_maxwellian(z2_electron,
     &   wdoti6,
     &   darea, nrhomax, nnoderho, drho,
     &   dx, dy, rt, nxmx, nymx, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xmi6, qi6, xn6a, xkti6, omgci6, omgp62, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   ealphak, ebetak, ebk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomia, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndisti6, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper6, dfdupar6, fperp6,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc6_mks, df6_cql_uprp, df6_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdoti, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,
     &   n_bin, upshift, i_write, xk_cutoff,  z0_table1, z1_table1,
     &   z2_table1, ntable, mtable, zetai_table, dKdL_table,
     &   nmax, mmax)
      end if


      tmin = (second1(dummy) - t1i6) / 60.
      if (myid.eq.0) then
         write(6 , 2876) tmin
         write(15, 2876) tmin
      endif
 2876 format('time for 6th ion species QL operator =', f9.3, ' min')

      end if

      end if


      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2838) tmin
         write(15, 2838) tmin
      endif

 2838 format('time to calculate quasilinear operator =', f9.3, ' min')





      t1 = second1(dummy)

      ! ------------------ !
      ! --  Electrons   -- !
      ! ------------------ !


      call wdot_qlcheck(wdote_ql,
     &   nnoderho, nrhomax,
     &   bqlavg_e, cqlavg_e, xme, omgrf, xkteavg, 0,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfdupere, dfdupare, fperpe,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vce_mks, dfe_cql_uprp, dfe_cql_uprl, rhon, rho_a, myid,
     &   dldbavg)



      ! ------------------ !
      ! --Majority ions -- !
      ! ------------------ !

      call wdot_qlcheck(wdoti1_ql,
     &   nnoderho, nrhomax,
     &   bqlavg_i1, cqlavg_i1, xmi1, omgrf, xktiavg, ndisti1,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper1, dfdupar1, fperp1,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc1_mks, df1_cql_uprp, df1_cql_uprl, rhon, rho_a, myid,
     &   dldbavg)



      ! ------------------ !
      ! --Minority ions -- !
      ! ------------------ !

      if(eta2 .ne. 0.0)then

      call wdot_qlcheck(wdoti2_ql,
     &   nnoderho, nrhomax,
     &   bqlavg_i2, cqlavg_i2, xmi2, omgrf, xkti2avg, ndisti2,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper2, dfdupar2, fperp2,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc2_mks, df2_cql_uprp, df2_cql_uprl, rhon, rho_a, myid,
     &   dldbavg)

      end if



      ! --------------------- !
      ! --Third ion species-- !
      ! --------------------- !

      if(eta3 .ne. 0.0)then

      call wdot_qlcheck(wdoti3_ql,
     &   nnoderho, nrhomax,
     &   bqlavg_i3, cqlavg_i3, xmi3, omgrf, xkti3avg, ndisti3,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper3, dfdupar3, fperp3,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc3_mks, df3_cql_uprp, df3_cql_uprl, rhon, rho_a, myid,
     &   dldbavg)

      end if

      ! --------------------- !
      ! --Fourth ion species-- !
      ! --------------------- !

      if(eta4 .ne. 0.0)then

      call wdot_qlcheck(wdoti4_ql,
     &   nnoderho, nrhomax,
     &   bqlavg_i4, cqlavg_i4, xmi4, omgrf, xkti4avg, ndisti4,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper4, dfdupar4, fperp4,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc4_mks, df4_cql_uprp, df4_cql_uprl, rhon, rho_a, myid,
     &   dldbavg)

      end if

      ! --------------------- !
      ! --Fifth ion species-- !
      ! --------------------- !

      if(eta5 .ne. 0.0)then

      call wdot_qlcheck(wdoti5_ql,
     &   nnoderho, nrhomax,
     &   bqlavg_i5, cqlavg_i5, xmi5, omgrf, xkti5avg, ndisti5,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper5, dfdupar5, fperp5,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc5_mks, df5_cql_uprp, df5_cql_uprl, rhon, rho_a, myid,
     &   dldbavg)

      end if

      ! --------------------- !
      ! --Sixth ion species-- !
      ! --------------------- !

      if(eta6 .ne. 0.0)then

      call wdot_qlcheck(wdoti6_ql,
     &   nnoderho, nrhomax,
     &   bqlavg_i6, cqlavg_i6, xmi6, omgrf, xkti6avg, ndisti6,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper6, dfdupar6, fperp6,
     &   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     &   vc6_mks, df6_cql_uprp, df6_cql_uprl, rhon, rho_a, myid,
     &   dldbavg)

      end if

c      end if


      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2839) tmin
         write(15, 2839) tmin
      endif

 2839 format('time to check the quasilinear operator =', f9.3, ' min')

      t1 = second1(dummy)

*     -------------------------------------------------------
*     Calculate E* dot J, Poynting flux (sp), and zero order
*     poloidal force (fz0) in real :: space:
*     -------------------------------------------------------

      do i = 1, nnodex
         do j = 1, nnodey

            fx0(i,j) = fx0e(i,j) + fx0i1(i,j) + fx0i2(i,j) + fx0i3(i,j)
     &                           + fx0i4(i,j) + fx0i5(i,j) + fx0i6(i,j)
            fy0(i,j) = fy0e(i,j) + fy0i1(i,j) + fy0i2(i,j) + fy0i3(i,j)
     &                           + fy0i4(i,j) + fy0i5(i,j) + fy0i6(i,j)
            fz0(i,j) = fz0e(i,j) + fz0i1(i,j) + fz0i2(i,j) + fz0i3(i,j)
     &                           + fz0i4(i,j) + fz0i5(i,j) + fz0i6(i,j)


            bdotf(i,j) = bmod(i,j) * (bxn(i,j) * fx0(i,j)
     &                             +  byn(i,j) * fy0(i,j)
     &                             +  bzn(i,j) * fz0(i,j))

            bmod2(i,j) = bmod(i,j)**2

            capr_fzeta(i,j) = capr(i) * fz0(i,j)
            capr2(i,j) = capr(i)**2

            wdot(i,j) = wdote(i,j)
     &                + wdoti1(i,j) + wdoti2(i,j) + wdoti3(i,j)
     &                + wdoti4(i,j) + wdoti5(i,j) + wdoti6(i,j)

            redotj(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx(i, j)
     &                              + conjg(ebeta(i, j)) * xjpy(i, j)
     &                              + conjg(eb(i, j)) * xjpz(i, j) )

            divq(i, j) = redotj(i,j) - wdot(i,j)

            rekxdotj(i, j) = 0.5 / omgrf *
     &                         real(conjg(ealphakx(i, j)) * xjpx(i, j)
     &                             + conjg(ebetakx(i, j)) * xjpy(i, j)
     &                                + conjg(ebkx(i, j)) * xjpz(i, j) )

            rekydotj(i, j) = 0.5 / omgrf *
     &                         real(conjg(ealphaky(i, j)) * xjpx(i,j)
     &                             + conjg(ebetaky(i, j)) * xjpy(i, j)
     &                                + conjg(ebky(i, j)) * xjpz(i, j) )

            rekzdotj(i, j) = xkphi(i) / omgrf * redotj(i, j)



            pc(i, j) = -0.5 * (conjg(ex(i, j)) * xjx(i, j)
     &                       + conjg(ey(i, j)) * xjy(i, j)
     &                       + conjg(ez(i, j)) * xjz(i, j) )

            pcre(i, j) = real(pc(i, j))
            pcim(i, j) = aimag(pc(i, j))


            bxwave(i, j) = 0.0
            bzwave(i, j) = 0.0
            bywave(i, j) = 0.0

c           rhohatx(i,j) = 0.0
c            rhohaty(i,j) = 0.0

            if (i .ne. 1 .and. i .ne. nnodex) then
               if (j .ne. 1 .and. j .ne. nnodey) then

c                   drhodx = (rho(i+1, j) - rho(i-1, j)) / (2.0 * dx)
c                   drhody = (rho(i, j+1) - rho(i, j-1)) / (2.0 * dy)

c                   gradrho = sqrt(drhodx**2 + drhody**2)
c                   if (gradrho .eq. 0.0) gradrho = 1.0e-08

c                   rhohatx(i,j) = drhodx / gradrho
c                   rhohaty(i,j) = drhody / gradrho

                   dezdx = (ez(i+1, j) - ez(i-1, j)) / (2.0 * dx)
                   dezdy = (ez(i, j+1) - ez(i, j-1)) / (2.0 * dy)
                   deydx = (ey(i+1, j) - ey(i-1, j)) / (2.0 * dx)
                   dexdy = (ex(i, j+1) - ex(i, j-1)) / (2.0 * dy)

                   bxwave(i, j) = 1. / (zi * omgrf)
     1                * (zi * xkphi(i) * ey(i,j) - dezdy)
                   bzwave(i, j) = 1. / (zi * omgrf)
     1                * (dexdy - deydx)
                   bywave(i, j) = 1. / (zi * omgrf) * (ez(i,j) / capr(i)
     1                + dezdx - zi * xkphi(i) * ex(i,j))


                   spx(i, j) = 1. / (2. * xmu0) * real(
     1                         conjg(ez(i,j)) * bywave(i,j)
     1                       - conjg(ey(i,j)) * bzwave(i,j)  )
                   spz(i, j) = 1. / (2. * xmu0) * real(
     1                         conjg(ey(i,j)) * bxwave(i,j)
     1                       - conjg(ex(i,j)) * bywave(i,j)  )
                   spy(i, j) = 1. / (2. * xmu0) * real(
     1                         conjg(ex(i,j)) * bzwave(i,j)
     1                       - conjg(ez(i,j)) * bxwave(i,j)  )

c                  dspxdx(i) = 1. / (2. * xmu0 * omgrf) * imag(
c    1                conjg(deydx(i)) * (deydx(i) - zi * xky * ex(i))
c    1              + conjg(ey(i)) * (d2eydx(i) - zi * xky * dexdx(i))
c    1              + conjg(dezdx(i)) * (dezdx(i) - zi * xkz * ex(i))
c    1              + conjg(ez(i)) * (d2ezdx(i) - zi * xkz * dexdx(i)))

c                  divsp(i) = dspxdx(i)
c     1               + zi * xky * spy(i) + zi * xkz * spz(i)
               end if
            end if

*           --------------------------------
*           Calculate radial and theta force:
*           --------------------------------

            fpsi0(i, j)  = fx0i2(i, j) * rhohatx(i, j)
     &                   + fy0i2(i, j) * rhohaty(i, j)

            ftheta0(i,j) = fy0i2(i, j) * rhohatx(i, j)
     &                   - fx0i2(i, j) * rhohaty(i, j)

            if (eta2 .gt. 0) then
            fpsi1(i,j) = 1.0 / (qi2 * xn2a(i, j) )
     &                 * fpsi0(i, j) / (capr(i) * bpol(i, j))
            end if

            ftheta1(i,j) = ftheta0(i,j) * bpol(i,j)

*           -------------------------------------------------------------
*           Calculate driven current from wdote(i,j) for xkprl = nphi / R:
*           First calculate C_ehst
*           --------------------------------------------------------------

c           xkprl = uzx(i, j) * xkxsav(n)
c     &                    + uzy(i, j) * xkysav(m)
c     &                    + uzz(i, j) * xkphi(i)

            xkprl = uzz(i, j) * xkphi(i)

            akprl = abs(xkprl)
            if(akprl .lt. 1.0e-05)akprl = 1.0e-05

            vphase = omgrf / akprl
            if (vphase .ge. clight) vphase = clight

            c_ehst = 0.0

            if(xkte(i, j) .ne. 0.0) then

               xlnlam = 24. - alog(sqrt(xnea(i, j)
     &                     /1.e+06) / (xkte(i, j)/ abs(qe)))
               wphase = vphase / sqrt(xkte(i, j) / xme)

               if (wphase .lt. 10.0 .and. wphase .gt. 0.1)then

*                 -----------------------------
*                 Assume circular flux surfaces
*                 -----------------------------
                     a = sqrt(x(i)**2 + y(j)**2)
                     rmaxa = rt + a
                     rmina = rt - a
                     bmaxa = b0 * rt / rmina
*                 ------------------------------------
*                 End circular flux surface assumption
*                 ------------------------------------

                  epsa = (rmaxa - rmina)/(rmaxa + rmina)

                  xmut2 = 1. - bmod(i, j) / bmaxa
                  eta0 = 8. * wphase**2 / (5. + zeffcd) +
     &                     ceta / zeffcd**.707 + damp / wphase
                  r1 = 1. - epsa**.77 * sqrt(3.5**2 + wphase**2)
     &                     /(3.5 * epsa**.77 + wphase)
                  xm1 = 1.0
                  c1 = 1.0

                  if(ftrap .ne. 0)then
                     if (xmut2 .gt. 1.0E-06) then
                         xm1 = 1.0 + afit * (xmut2**1.5 / wphase**3)
                         yt = (1. - xmut2) * wphase**2 / xmut2
                         if(yt.ge.0.0)c1 = 1. - exp(-(cfit*yt)**xnexp)
                     endif
                  endif

                  xjtild = c1 * xm1 * eta0 * r1

                  c_ehst = -19.19E+15 / xlnlam * (xkte(i, j)
     &                     / abs(qe))/ xnea(i, j) * xjtild * signkz
c                 c_ehst = 1.0
               end if
            end if

            if (vphase .ge. clight .or. akprl .lt. 1.e-05) c_ehst = 0.0

            xjprl(i,j) = c_ehst * wdote(i,j)

         end do
      end do

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2840) tmin
         write(15, 2840) tmin
      endif

 2840 format('time to calculate E* dot J and Poynting flux =', f9.3,
     & ' min')


      t1 = second1(dummy)


*     -------------------------
*     Integrate antenna current
*     -------------------------
      call sgrator(x, y, pcre, 1, nnodex, 1, nnodey, pcrto2,
     &   capr, rt,  nxmx, nymx)
      call sgrator(x, y, pcim, 1, nnodex, 1, nnodey, pcito2,
     &   capr, rt, nxmx, nymx)


*     ---------------------------
*     Integrate 1/2 Re (E* dot J):
*     ---------------------------
      call sgrator(x, y, redotj, 1, nnodex, 1, nnodey, ptot,
     &   capr, rt,  nxmx, nymx)

      call sgrator_core(x, y, redotj, 1, nnodex, 1, nnodey, ptot_core,
     &   capr, rt,  nxmx, nymx, rho)

      call sgrator_edge(x, y, redotj, 1, nnodex, 1, nnodey, ptot_edge,
     &   capr, rt,  nxmx, nymx, rho)

*     ---------------
*     Integrate Wdot:
*     ---------------
      call sgrator(x, y, wdot, 1, nnodex, 1, nnodey, ptot_wdot,
     &   capr, rt,  nxmx, nymx)

      call sgrator_core(x, y, wdot, 1, nnodex, 1, nnodey, wdot_core,
     &   capr, rt,  nxmx, nymx, rho)

      call sgrator_edge(x, y, wdot, 1, nnodex, 1, nnodey, wdot_edge,
     &   capr, rt,  nxmx, nymx, rho)


*     -----------------
*     Integrate Wdot_ql:
*     -----------------

      call rhograte(rhon, wdote_ql,  1, nnoderho, wdote_ql_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti1_ql, 1, nnoderho, wdoti1_ql_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti2_ql, 1, nnoderho, wdoti2_ql_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti3_ql, 1, nnoderho, wdoti3_ql_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti4_ql, 1, nnoderho, wdoti4_ql_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti5_ql, 1, nnoderho, wdoti5_ql_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti6_ql, 1, nnoderho, wdoti6_ql_int,
     &   nrhomax, dvol)


      pe_ql  =  wdote_ql_int(nnoderho)
      pi1_ql = wdoti1_ql_int(nnoderho)
      pi2_ql = wdoti2_ql_int(nnoderho)
      pi3_ql = wdoti3_ql_int(nnoderho)
      pi4_ql = wdoti4_ql_int(nnoderho)
      pi5_ql = wdoti5_ql_int(nnoderho)
      pi6_ql = wdoti6_ql_int(nnoderho)
      pt_ql  = pe_ql + pi1_ql + pi2_ql + pi3_ql
     &                   + pi4_ql +pi5_ql + pi6_ql



*     ------------------
*     Normalize to jdote:
*     ------------------
      if(norm .eq. 1)powtot = 2.0 * pi * rt * ptot

*     -----------------
*     Normalize to wdot:
*     -----------------
      if(norm .eq. 2)powtot = 2.0 * pi * rt * ptot_wdot

*     --------------------
*     Normalize to wdot_ql:
*     --------------------
      if(norm .eq. 3)powtot = pt_ql


      if (powtot .ne. 0.0) pscale = prfin / powtot
      if (powtot .eq. 0.0) pscale = 1.0

      if (prfin .eq. 0.0) pscale = 1.0
      if (pscale .le. 0.0) pscale = 1.0

      if (myid .eq. 0) then
         write(6, 1218) pscale
         write(15, 1218) pscale
      endif

*     -----------------------
*     Scale powers to Prf,in:
*     -----------------------

      pe_ql  = pe_ql  * pscale
      pi1_ql = pi1_ql * pscale
      pi2_ql = pi2_ql * pscale
      pi3_ql = pi3_ql * pscale
      pi4_ql = pi4_ql * pscale
      pi5_ql = pi5_ql * pscale
      pi6_ql = pi6_ql * pscale
      pt_ql  = pt_ql  * pscale


      powtot = powtot * pscale
      ptot = ptot * pscale

      ptot_edge = ptot_edge * pscale
      ptot_core = ptot_core * pscale

      pcrto2 = pcrto2 * pscale
      pcito2 = pcito2 * pscale

      fxtot_ant = fxtot_ant * pscale
      fytot_ant = fytot_ant * pscale
      fztot_ant = fztot_ant * pscale

      fxtot_plasma = fxtot_plasma * pscale
      fytot_plasma = fytot_plasma * pscale
      fztot_plasma = fztot_plasma * pscale


*     ------------------------------------
*     Scale fields and current to Prf, in:
*     ------------------------------------

      xjx = xjx * sqrt(pscale)
      xjy = xjy * sqrt(pscale)
      xjz = xjz * sqrt(pscale)

      do n = nkx1, nkx2
         do m = nky1, nky2
            ealphakmod(n, m) = ealphakmod(n, m) * sqrt(pscale)
            ebetakmod(n, m) = ebetakmod(n, m) * sqrt(pscale)
            ebkmod(n, m) = ebkmod(n, m) * sqrt(pscale)
         end do
      end do

      do n = 1, nnoderho
         do m = 1, mnodetheta
            eplus_flux(n, m) = eplus_flux(n, m) * sqrt(pscale)
            eminus_flux(n, m) = eminus_flux(n, m) * sqrt(pscale)
         end do
      end do

      do i = 1, nnodex
         do j = 1, nnodey

            eplus_flux_plot(i,j) = eplus_flux_plot(i,j) * sqrt(pscale)
            eminus_flux_plot(i,j) = eminus_flux_plot(i,j) * sqrt(pscale)

            ex(i,j) = ex(i,j) * sqrt(pscale)
            ey(i,j) = ey(i,j) * sqrt(pscale)
            ez(i,j) = ez(i,j) * sqrt(pscale)

            bxwave(i,j) = bxwave(i,j) * sqrt(pscale)
            bywave(i,j) = bywave(i,j) * sqrt(pscale)
            bzwave(i,j) = bzwave(i,j) * sqrt(pscale)

            ealpha(i,j) = ealpha(i,j) * sqrt(pscale)
            ebeta(i,j) = ebeta(i,j) * sqrt(pscale)
            eb(i,j) = eb(i,j) * sqrt(pscale)

            eplus(i,j) = eplus(i,j) * sqrt(pscale)
            eminus(i,j) = eplus(i,j) * sqrt(pscale)

            ntilda_e(i,j) = ntilda_e(i,j) * sqrt(pscale)
            ntilda_e_real(i,j) = real(ntilda_e(i,j))

            xjpxe(i,j) = xjpxe(i,j) * sqrt(pscale)
            xjpye(i,j) = xjpye(i,j) * sqrt(pscale)
            xjpze(i,j) = xjpze(i,j) * sqrt(pscale)

            xjpxe_lab(i,j) = xjpxe_lab(i,j) * sqrt(pscale)
            xjpye_lab(i,j) = xjpye_lab(i,j) * sqrt(pscale)
            xjpze_lab(i,j) = xjpze_lab(i,j) * sqrt(pscale)

            xjpx_ehst(i,j) = xjpx_ehst(i,j) * sqrt(pscale)
            xjpy_ehst(i,j) = xjpy_ehst(i,j) * sqrt(pscale)
            xjpz_ehst(i,j) = xjpz_ehst(i,j) * sqrt(pscale)

            xjpx1(i,j) = xjpx1(i,j) * sqrt(pscale)
            xjpy1(i,j) = xjpy1(i,j) * sqrt(pscale)
            xjpz1(i,j) = xjpz1(i,j) * sqrt(pscale)

            xjpx2(i,j) = xjpx2(i,j) * sqrt(pscale)
            xjpy2(i,j) = xjpy2(i,j) * sqrt(pscale)
            xjpz2(i,j) = xjpz2(i,j) * sqrt(pscale)

            xjpx3(i,j) = xjpx3(i,j) * sqrt(pscale)
            xjpy3(i,j) = xjpy3(i,j) * sqrt(pscale)
            xjpz3(i,j) = xjpz3(i,j) * sqrt(pscale)

            xjpx4(i,j) = xjpx4(i,j) * sqrt(pscale)
            xjpy4(i,j) = xjpy4(i,j) * sqrt(pscale)
            xjpz4(i,j) = xjpz4(i,j) * sqrt(pscale)

            xjpx5(i,j) = xjpx5(i,j) * sqrt(pscale)
            xjpy5(i,j) = xjpy5(i,j) * sqrt(pscale)
            xjpz5(i,j) = xjpz5(i,j) * sqrt(pscale)

            xjpx6(i,j) = xjpx6(i,j) * sqrt(pscale)
            xjpy6(i,j) = xjpy6(i,j) * sqrt(pscale)
            xjpz6(i,j) = xjpz6(i,j) * sqrt(pscale)

            xj_slox(i,j) = xj_slox(i,j) * sqrt(pscale)
            xj_sloy(i,j) = xj_sloy(i,j) * sqrt(pscale)
            xj_sloz(i,j) = xj_sloz(i,j) * sqrt(pscale)


            redotj(i,j) = redotj(i,j) * pscale
            divq(i,j) = divq(i,j) * pscale

            pc(i,j) = pc(i,j) * pscale
            pcre(i,j) = pcre(i,j) * pscale
            pcim(i,j) = pcim(i,j) * pscale

            wdote(i,j)  = wdote(i,j)  * pscale
            wdoti1(i,j) = wdoti1(i,j) * pscale
            wdoti2(i,j) = wdoti2(i,j) * pscale
            wdoti3(i,j) = wdoti3(i,j) * pscale
            wdoti4(i,j) = wdoti4(i,j) * pscale
            wdoti5(i,j) = wdoti5(i,j) * pscale
            wdoti6(i,j) = wdoti6(i,j) * pscale

            wdot(i,j)   = wdot(i,j)  * pscale

            xjprl(i,j)  = xjprl(i,j)  * pscale

            fype(i,j)  = fype(i,j)  * pscale
            fypi1(i,j) = fypi1(i,j) * pscale
            fypi2(i,j) = fypi2(i,j) * pscale
            fypi3(i,j) = fypi3(i,j) * pscale
            fypi4(i,j) = fypi4(i,j) * pscale
            fypi5(i,j) = fypi5(i,j) * pscale
            fypi6(i,j) = fypi6(i,j) * pscale
            fyp(i,j)   = fyp(i,j)   * pscale


            fz0e(i,j)  = fz0e(i,j)  * pscale
            fz0i1(i,j) = fz0i1(i,j) * pscale
            fz0i2(i,j) = fz0i2(i,j) * pscale
            fz0i3(i,j) = fz0i3(i,j) * pscale
            fz0i4(i,j) = fz0i4(i,j) * pscale
            fz0i5(i,j) = fz0i5(i,j) * pscale
            fz0i6(i,j) = fz0i6(i,j) * pscale
            fz0(i,j)   = fz0(i,j)   * pscale

            capr_fzeta(i,j) = capr_fzeta(i,j) * pscale

            fpsi0(i, j)  = fpsi0(i, j)  * pscale
            ftheta0(i,j) = ftheta0(i,j) * pscale
            bdotf(i,j) = bdotf(i,j) * pscale

            fpsi1(i,j)= fpsi1(i,j) * pscale
            ftheta1(i,j) = ftheta1(i,j) * pscale

         end do
      end do

      do n = 1, nnoderho
         wdote_ql(n)  = wdote_ql(n)  * pscale
         wdoti1_ql(n) = wdoti1_ql(n) * pscale
         wdoti2_ql(n) = wdoti2_ql(n) * pscale
         wdoti3_ql(n) = wdoti3_ql(n) * pscale
         wdoti4_ql(n) = wdoti4_ql(n) * pscale
         wdoti5_ql(n) = wdoti5_ql(n) * pscale
         wdoti6_ql(n) = wdoti6_ql(n) * pscale

         do niu = 1, nuper
            do miu = 1, nupar

               bqlavg_e(niu, miu, n) = bqlavg_e(niu, miu, n) * pscale
               cqlavg_e(niu, miu, n) = cqlavg_e(niu, miu, n) * pscale
               eqlavg_e(niu, miu, n) = eqlavg_e(niu, miu, n) * pscale
               fqlavg_e(niu, miu, n) = fqlavg_e(niu, miu, n) * pscale


               bqlavg_i1(niu, miu, n) = bqlavg_i1(niu, miu, n) * pscale
               cqlavg_i1(niu, miu, n) = cqlavg_i1(niu, miu, n) * pscale
               eqlavg_i1(niu, miu, n) = eqlavg_i1(niu, miu, n) * pscale
               fqlavg_i1(niu, miu, n) = fqlavg_i1(niu, miu, n) * pscale

               if(eta2 .ne. 0.0)then
               bqlavg_i2(niu, miu, n) = bqlavg_i2(niu, miu, n) * pscale
               cqlavg_i2(niu, miu, n) = cqlavg_i2(niu, miu, n) * pscale
               eqlavg_i2(niu, miu, n) = eqlavg_i2(niu, miu, n) * pscale
               fqlavg_i2(niu, miu, n) = fqlavg_i2(niu, miu, n) * pscale
               end if

               if(eta3 .ne. 0.0)then
               bqlavg_i3(niu, miu, n) = bqlavg_i3(niu, miu, n) * pscale
               cqlavg_i3(niu, miu, n) = cqlavg_i3(niu, miu, n) * pscale
               eqlavg_i3(niu, miu, n) = eqlavg_i3(niu, miu, n) * pscale
               fqlavg_i3(niu, miu, n) = fqlavg_i3(niu, miu, n) * pscale
               end if

               if(eta4 .ne. 0.0)then
               bqlavg_i4(niu, miu, n) = bqlavg_i4(niu, miu, n) * pscale
               cqlavg_i4(niu, miu, n) = cqlavg_i4(niu, miu, n) * pscale
               eqlavg_i4(niu, miu, n) = eqlavg_i4(niu, miu, n) * pscale
               fqlavg_i4(niu, miu, n) = fqlavg_i4(niu, miu, n) * pscale
               end if

               if(eta5 .ne. 0.0)then
               bqlavg_i5(niu, miu, n) = bqlavg_i5(niu, miu, n) * pscale
               cqlavg_i5(niu, miu, n) = cqlavg_i5(niu, miu, n) * pscale
               eqlavg_i5(niu, miu, n) = eqlavg_i5(niu, miu, n) * pscale
               fqlavg_i5(niu, miu, n) = fqlavg_i5(niu, miu, n) * pscale
               end if

               if(eta6 .ne. 0.0)then
               bqlavg_i6(niu, miu, n) = bqlavg_i6(niu, miu, n) * pscale
               cqlavg_i6(niu, miu, n) = cqlavg_i6(niu, miu, n) * pscale
               eqlavg_i6(niu, miu, n) = eqlavg_i6(niu, miu, n) * pscale
               fqlavg_i6(niu, miu, n) = fqlavg_i6(niu, miu, n) * pscale
               end if


            end do
         end do

      end do

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2841) tmin
         write(15, 2841) tmin
      endif

 2841 format('time to do power scaling =', f9.3, ' min')



      t1 = second1(dummy)



*     -------------------------------------------------------------------
*     Convert quasilinear diffusion coefficients to CGS units for
*     Bob Harvey and CQL3D:
*
*     1. divide out the density since CQL includes density in f
*     2. multiply by 100 to change xlamda to cm
*     3. multiply by powers of vc_cgs to change to a dimensional velocity
*     -------------------------------------------------------------------

      vce_cgs = vce_mks * 100.
      vc1_cgs = vc1_mks * 100.
      vc2_cgs = vc2_mks * 100.
      vc3_cgs = vc3_mks * 100.
      vc4_cgs = vc4_mks * 100.
      vc5_cgs = vc5_mks * 100.
      vc6_cgs = vc6_mks * 100.

      do n = 1, nnoderho
         do niu = 1, nuper
            do miu = 1, nupar

               if(xnavg(n) .ne. 0.0)then
               bqlavg_e(niu, miu, n) = bqlavg_e(niu, miu, n) / xnavg(n)
     &            * vce_cgs**4 * 100.
               cqlavg_e(niu, miu, n) = cqlavg_e(niu, miu, n) / xnavg(n)
     &            * vce_cgs**3 * 100.
               eqlavg_e(niu, miu, n) = eqlavg_e(niu, miu, n) / xnavg(n)
     &            * vce_cgs**3 * 100.
               fqlavg_e(niu, miu, n) = fqlavg_e(niu, miu, n) / xnavg(n)
     &            * vce_cgs**2 * 100.

               else
                  bqlavg_e(niu, miu, n) = 0.
                  cqlavg_e(niu, miu, n) = 0.
                  eqlavg_e(niu, miu, n) = 0.
                  fqlavg_e(niu, miu, n) = 0.
               end if

               if(xn1avg(n) .ne. 0.0)then
               bqlavg_i1(niu, miu, n)= bqlavg_i1(niu, miu, n)/ xn1avg(n)
     &            * vc1_cgs**4 * 100.
               cqlavg_i1(niu, miu, n)= cqlavg_i1(niu, miu, n)/ xn1avg(n)
     &            * vc1_cgs**3 * 100.
               eqlavg_i1(niu, miu, n)= eqlavg_i1(niu, miu, n)/ xn1avg(n)
     &            * vc1_cgs**3 * 100.
               fqlavg_i1(niu, miu, n)= fqlavg_i1(niu, miu, n)/ xn1avg(n)
     &            * vc1_cgs**2 * 100.

               else
                  bqlavg_i1(niu, miu, n)= 0.
                  cqlavg_i1(niu, miu, n)= 0.
                  eqlavg_i1(niu, miu, n)= 0.
                  fqlavg_i1(niu, miu, n)= 0.
               end if

               if(eta2 .ne. 0.0)then
               if(xn2avg(n) .ne. 0.0)then
               bqlavg_i2(niu, miu, n)= bqlavg_i2(niu, miu, n)/ xn2avg(n)
     &            * vc2_cgs**4 * 100.
               cqlavg_i2(niu, miu, n)= cqlavg_i2(niu, miu, n)/ xn2avg(n)
     &            * vc2_cgs**3 * 100.
               eqlavg_i2(niu, miu, n)= eqlavg_i2(niu, miu, n)/ xn2avg(n)
     &            * vc2_cgs**3 * 100.
               fqlavg_i2(niu, miu, n)= fqlavg_i2(niu, miu, n)/ xn2avg(n)
     &            * vc2_cgs**2 * 100.

               else
                  bqlavg_i2(niu, miu, n)= 0.
                  cqlavg_i2(niu, miu, n)= 0.
                  eqlavg_i2(niu, miu, n)= 0.
                  fqlavg_i2(niu, miu, n)= 0.
               end if
               end if


               if(eta3 .ne. 0.0)then
               if(xn3avg(n) .ne. 0.0)then
               bqlavg_i3(niu, miu, n)= bqlavg_i3(niu, miu, n)/ xn3avg(n)
     &            * vc3_cgs**4 * 100.
               cqlavg_i3(niu, miu, n)= cqlavg_i3(niu, miu, n)/ xn3avg(n)
     &            * vc3_cgs**3 * 100.
               eqlavg_i3(niu, miu, n)= eqlavg_i3(niu, miu, n)/ xn3avg(n)
     &            * vc3_cgs**3 * 100.
               fqlavg_i3(niu, miu, n)= fqlavg_i3(niu, miu, n)/ xn3avg(n)
     &            * vc3_cgs**2 * 100.

               else
                  bqlavg_i3(niu, miu, n)= 0.
                  cqlavg_i3(niu, miu, n)= 0.
                  eqlavg_i3(niu, miu, n)= 0.
                  fqlavg_i3(niu, miu, n)= 0.
               end if
               end if


               if(eta4 .ne. 0.0)then
               if(xn4avg(n) .ne. 0.0)then
               bqlavg_i4(niu, miu, n)= bqlavg_i4(niu, miu, n)/ xn4avg(n)
     &            * vc4_cgs**4 * 100.
               cqlavg_i4(niu, miu, n)= cqlavg_i4(niu, miu, n)/ xn4avg(n)
     &            * vc4_cgs**3 * 100.
               eqlavg_i4(niu, miu, n)= eqlavg_i4(niu, miu, n)/ xn4avg(n)
     &            * vc4_cgs**3 * 100.
               fqlavg_i4(niu, miu, n)= fqlavg_i4(niu, miu, n)/ xn4avg(n)
     &            * vc4_cgs**2 * 100.

               else
                  bqlavg_i4(niu, miu, n)= 0.
                  cqlavg_i4(niu, miu, n)= 0.
                  eqlavg_i4(niu, miu, n)= 0.
                  fqlavg_i4(niu, miu, n)= 0.
               end if
               end if

               if(eta5 .ne. 0.0)then
               if(xn5avg(n) .ne. 0.0)then
               bqlavg_i5(niu, miu, n)= bqlavg_i5(niu, miu, n)/ xn5avg(n)
     &            * vc5_cgs**4 * 100.
               cqlavg_i5(niu, miu, n)= cqlavg_i5(niu, miu, n)/ xn5avg(n)
     &            * vc5_cgs**3 * 100.
               eqlavg_i5(niu, miu, n)= eqlavg_i5(niu, miu, n)/ xn5avg(n)
     &            * vc5_cgs**3 * 100.
               fqlavg_i5(niu, miu, n)= fqlavg_i5(niu, miu, n)/ xn5avg(n)
     &            * vc5_cgs**2 * 100.

               else
                  bqlavg_i5(niu, miu, n)= 0.
                  cqlavg_i5(niu, miu, n)= 0.
                  eqlavg_i5(niu, miu, n)= 0.
                  fqlavg_i5(niu, miu, n)= 0.
               end if
               end if

               if(eta6 .ne. 0.0)then
               if(xn6avg(n) .ne. 0.0)then
               bqlavg_i6(niu, miu, n)= bqlavg_i6(niu, miu, n)/ xn6avg(n)
     &            * vc6_cgs**4 * 100.
               cqlavg_i6(niu, miu, n)= cqlavg_i6(niu, miu, n)/ xn6avg(n)
     &            * vc6_cgs**3 * 100.
               eqlavg_i6(niu, miu, n)= eqlavg_i6(niu, miu, n)/ xn6avg(n)
     &            * vc6_cgs**3 * 100.
               fqlavg_i6(niu, miu, n)= fqlavg_i6(niu, miu, n)/ xn6avg(n)
     &            * vc6_cgs**2 * 100.

               else
                  bqlavg_i6(niu, miu, n)= 0.
                  cqlavg_i6(niu, miu, n)= 0.
                  eqlavg_i6(niu, miu, n)= 0.
                  fqlavg_i6(niu, miu, n)= 0.
               end if
               end if


            end do
         end do
      end do

*     --------------------------------------------------------
*     Interpolate quasilinear coefficients onto the CQL3D grid
*     --------------------------------------------------------

      if(enorm_factor .ne. 0.0) then

          if(eta1 .ne. 0.0 .and. ndisti1 .eq. 1)then

            vperp = uperp * vc1_mks
            vpara = upara * vc1_mks

            do n = 1, nnoderho
               do i_uperp = 1, nuper
                  do i_upara = 1, nupar
                     vperp_in = uperp(i_uperp) * vc1_mks_cql3d
                     vpara_in = upara(i_upara) * vc1_mks_cql3d

                     call intplt_to_cql3d(vperp_in, vpara_in,
     &                        bqlavg_work(i_uperp, i_upara, n),
     &                        nuper, nupar, nnoderho, bqlavg_i1,
     &                        nuper, nupar, nnoderho, vperp, vpara, n)

                     call intplt_to_cql3d(vperp_in, vpara_in,
     &                        cqlavg_work(i_uperp, i_upara, n),
     &                        nuper, nupar, nnoderho, cqlavg_i1,
     &                        nuper, nupar, nnoderho, vperp, vpara, n)

                     call intplt_to_cql3d(vperp_in, vpara_in,
     &                        eqlavg_work(i_uperp, i_upara, n),
     &                        nuper, nupar, nnoderho, eqlavg_i1,
     &                        nuper, nupar, nnoderho, vperp, vpara, n)

                     call intplt_to_cql3d(vperp_in, vpara_in,
     &                        fqlavg_work(i_uperp, i_upara, n),
     &                        nuper, nupar, nnoderho, fqlavg_i1,
     &                        nuper, nupar, nnoderho, vperp, vpara, n)
                  end do
               end do
            end do

            bqlavg_i1 = bqlavg_work
            cqlavg_i1 = cqlavg_work
            eqlavg_i1 = eqlavg_work
            fqlavg_i1 = fqlavg_work

         end if



         if(eta2 .ne. 0.0 .and. ndisti2 .eq. 1) then

            vperp = uperp * vc2_mks
            vpara = upara * vc2_mks

            do n = 1, nnoderho
               do i_uperp = 1, nuper
                  do i_upara = 1, nupar
                     vperp_in = uperp(i_uperp) * vc2_mks_cql3d
                     vpara_in = upara(i_upara) * vc2_mks_cql3d

                     call intplt_to_cql3d(vperp_in, vpara_in,
     &                        bqlavg_work(i_uperp, i_upara, n),
     &                        nuper, nupar, nnoderho, bqlavg_i2,
     &                        nuper, nupar, nnoderho, vperp, vpara, n)

                     call intplt_to_cql3d(vperp_in, vpara_in,
     &                        cqlavg_work(i_uperp, i_upara, n),
     &                        nuper, nupar, nnoderho, cqlavg_i2,
     &                        nuper, nupar, nnoderho, vperp, vpara, n)

                     call intplt_to_cql3d(vperp_in, vpara_in,
     &                        eqlavg_work(i_uperp, i_upara, n),
     &                        nuper, nupar, nnoderho, eqlavg_i2,
     &                        nuper, nupar, nnoderho, vperp, vpara, n)

                     call intplt_to_cql3d(vperp_in, vpara_in,
     &                        fqlavg_work(i_uperp, i_upara, n),
     &                        nuper, nupar, nnoderho, fqlavg_i2,
     &                        nuper, nupar, nnoderho, vperp, vpara, n)
                  end do
               end do
            end do

            bqlavg_i2 = bqlavg_work
            cqlavg_i2 = cqlavg_work
            eqlavg_i2 = eqlavg_work
            fqlavg_i2 = fqlavg_work

         end if


      end if





      if (myid .eq. 0) then

*        -------------------------------------------------------------------------
*        Write electron quasilinear diffusion coefficients to file out_cql3d.coefe
*        -------------------------------------------------------------------------
         if(ndiste .eq. 1)then

            write (40, 309) nuper
            write (40, 309) nupar
            write (40, 309) nnoderho

            write (40, 3310) vce_cgs
            write (40, 3310) UminPara, UmaxPara

            write (40, 3310) (rhon(n), n = 1, nnoderho)
            write (40, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (40, 3310) (upara(i_upara), i_upara = 1, nupar)


            write (40, 3310) (((bqlavg_e(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (40, 3310) (((cqlavg_e(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (40, 3310) (((eqlavg_e(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (40, 3310) (((fqlavg_e(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (40, 3310) xme

            close (40)
            call dql_write_nc('out_cql3d.coefe','out_cql3d.coefe.nc')

         end if

*        ----------------------------------------------------------------
*        Write quasilinear diffusion coefficients to file out_cql3d.coef1
*        ----------------------------------------------------------------
         if(ndisti1 .eq. 1)then

            write (41, 309) nuper
            write (41, 309) nupar
            write (41, 309) nnoderho

            write (41, 3310) vc1_cgs_cql3d
            write (41, 3310) UminPara, UmaxPara

            write (41, 3310) (rhon(n), n = 1, nnoderho)
            write (41, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (41, 3310) (upara(i_upara), i_upara = 1, nupar)


            write (41, 3310) (((bqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (41, 3310) (((cqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (41, 3310) (((eqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (41, 3310) (((fqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (41, 3310) xmi1

            close (41)
            call dql_write_nc('out_cql3d.coef1','out_cql3d.coef1.nc')

         end if


*        ----------------------------------------------------------------
*        Write quasilinear diffusion coefficients to file out_cql3d.coef2
*        ----------------------------------------------------------------
         if(eta2 .ne. 0.0 .and. ndisti2 .eq. 1)then

            write (42, 309) nuper
            write (42, 309) nupar
            write (42, 309) nnoderho

            write (42, 3310) vc2_cgs_cql3d
            write (42, 3310) UminPara, UmaxPara

            write (42, 3310) (rhon(n), n = 1, nnoderho)
            write (42, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (42, 3310) (upara(i_upara), i_upara = 1, nupar)

            write (42, 3310) (((bqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (42, 3310) (((cqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (42, 3310) (((eqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (42, 3310) (((fqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (42, 3310) xmi2

            close (42)
            call dql_write_nc('out_cql3d.coef2','out_cql3d.coef2.nc')

         end if

*        ----------------------------------------------------------------
*        Write quasilinear diffusion coefficients to file out_cql3d.coef3
*        ----------------------------------------------------------------
         if(eta3 .ne. 0.0 .and. ndisti3 .eq. 1)then

            write (43, 309) nuper
            write (43, 309) nupar
            write (43, 309) nnoderho

            write (43, 3310) vc3_cgs_cql3d
            write (43, 3310) UminPara, UmaxPara

            write (43, 3310) (rhon(n), n = 1, nnoderho)
            write (43, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (43, 3310) (upara(i_upara), i_upara = 1, nupar)

            write (43, 3310) (((bqlavg_i3(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (43, 3310) (((cqlavg_i3(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (43, 3310) (((eqlavg_i3(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (43, 3310) (((fqlavg_i3(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (43, 3310) xmi3

            close (43)
            call dql_write_nc('out_cql3d.coef3','out_cql3d.coef3.nc')

         end if


*        ----------------------------------------------------------------
*        Write quasilinear diffusion coefficients to file out_cql3d.coef4
*        ----------------------------------------------------------------
         if(eta4 .ne. 0.0 .and. ndisti4 .eq. 1)then

            write (44, 309) nuper
            write (44, 309) nupar
            write (44, 309) nnoderho

            write (44, 3310) vc4_cgs_cql3d
            write (44, 3310) UminPara, UmaxPara

            write (44, 3310) (rhon(n), n = 1, nnoderho)
            write (44, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (44, 3310) (upara(i_upara), i_upara = 1, nupar)

            write (44, 3310) (((bqlavg_i4(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (44, 3310) (((cqlavg_i4(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (44, 3310) (((eqlavg_i4(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (44, 3310) (((fqlavg_i4(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (44, 3310) xmi4

            close (44)
            call dql_write_nc('out_cql3d.coef4','out_cql3d.coef4.nc')


         end if


      end if


      tmin = (second1(dummy) - t1) / 60.

      if (myid.eq.0) then
         write(6 , 2842) tmin
         write(15, 2842) tmin
      endif

 2842 format('time to write quasilinear operator in cgs units =', f9.3,
     &    ' min')


      t1 = second1(dummy)
*     ------------------------------------------------------------------------
*     calculate individual species edotj's and driven current on solution mesh
*     ------------------------------------------------------------------------

      do i = 1, nnodex
         do j = 1, nnodey
            redotje(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpxe(i, j)
     &                               + conjg(ebeta(i, j)) * xjpye(i, j)
     &                               + conjg(eb(i, j)) * xjpze(i, j) )
            redotj_ehst(i, j) =
     &                 0.5 * real(conjg(ealpha(i, j)) * xjpx_ehst(i, j)
     &                           + conjg(ebeta(i, j)) * xjpy_ehst(i, j)
     &                            + conjg(eb(i, j)) * xjpz_ehst(i, j) )

            redotj1(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx1(i, j)
     &                               + conjg(ebeta(i, j)) * xjpy1(i, j)
     &                               + conjg(eb(i, j)) * xjpz1(i, j) )

            redotj2(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx2(i, j)
     &                               + conjg(ebeta(i, j)) * xjpy2(i, j)
     &                               + conjg(eb(i, j)) * xjpz2(i, j) )

            redotj3(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx3(i, j)
     &                               + conjg(ebeta(i, j)) * xjpy3(i, j)
     &                               + conjg(eb(i, j)) * xjpz3(i, j) )

            redotj4(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx4(i, j)
     &                               + conjg(ebeta(i, j)) * xjpy4(i, j)
     &                               + conjg(eb(i, j)) * xjpz4(i, j) )

            redotj5(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx5(i, j)
     &                               + conjg(ebeta(i, j)) * xjpy5(i, j)
     &                               + conjg(eb(i, j)) * xjpz5(i, j) )

            redotj6(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx6(i, j)
     &                               + conjg(ebeta(i, j)) * xjpy6(i, j)
     &                               + conjg(eb(i, j)) * xjpz6(i, j) )

            redotjs(i, j) = 0.5 *
     &                          real(conjg(ealpha(i, j)) * xj_slox(i, j)
     &                              + conjg(ebeta(i, j)) * xj_sloy(i, j)
     &                              + conjg(eb(i, j)) * xj_sloz(i, j) )

            redotjt(i, j) = redotje(i, j) + redotj1(i, j)
     &                    + redotj2(i, j) + redotj3(i, j)
     &                    + redotj4(i, j) + redotj5(i, j)
     &                    + redotj6(i, j) + redotjs(i, j)

            redotji(i, j) = redotj1(i, j)
     &                    + redotj2(i, j) + redotj3(i, j)
     &                    + redotj4(i, j) + redotj5(i, j)
     &                    + redotj6(i, j)



c           if(nnode_local .le. 0)xjprl(i,j) = redotj_ehst(i,j)

         end do
      end do


*     ----------------------------------------------------------
*     Integrate individual species Estar dot J on solution mesh:
*     ----------------------------------------------------------

      call sgrator(x, y, redotje, 1, nnodex, 1, nnodey, pedotje,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj1, 1, nnodex, 1, nnodey, pedotj1,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj2, 1, nnodex, 1, nnodey, pedotj2,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj3, 1, nnodex, 1, nnodey, pedotj3,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj4, 1, nnodex, 1, nnodey, pedotj4,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj5, 1, nnodex, 1, nnodey, pedotj5,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj6, 1, nnodex, 1, nnodey, pedotj6,
     &   capr, rt, nxmx, nymx)


      call sgrator(x, y, redotjs, 1, nnodex, 1, nnodey, pedotjs,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, redotjt, 1, nnodex, 1, nnodey, pedotjt,
     &   capr, rt, nxmx, nymx)

      call sgrator1(x, y, xjy, 1, nnodex, 1, nnodey, i_ant,
     &   capr, rt, nxmx, nymx, j_cur_max)

      call sgrator(x, y, xjprl, 1, nnodex, 1, nnodey, xjtot,
     &   capr, rt,  nxmx, nymx)

      call sgrator(x, y, wdote, 1, nnodex, 1, nnodey, pe,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti1, 1, nnodex, 1, nnodey, pi1,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti2, 1, nnodex, 1, nnodey, pi2,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti3, 1, nnodex, 1, nnodey, pi3,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti4, 1, nnodex, 1, nnodey, pi4,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti5, 1, nnodex, 1, nnodey, pi5,
     &   capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti6, 1, nnodex, 1, nnodey, pi6,
     &   capr, rt, nxmx, nymx)



      call sgrator(x, y, wdot, 1, nnodex, 1, nnodey, pt,
     &   capr, rt, nxmx, nymx)

      call sgrator(x, y, fyp, 1, nnodex, 1, nnodey, fyptot,
     &   capr, rt, nxmx, nymx)

      call sgrator(x, y, fz0, 1, nnodex, 1, nnodey, fz0tot,
     &   capr, rt, nxmx, nymx)

*     -------------------------
*     Do flux surface averages:
*     -------------------------

      call fluxavg(redotje, redotjeavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj1, redotj1avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj2, redotj2avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj3, redotj3avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj4, redotj4avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj5, redotj5avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj6, redotj6avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotjs, redotjsavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotjt, redotjtavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(xjprl, xjprlavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(wdote, wdoteavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti1, wdoti1avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti2, wdoti2avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti3, wdoti3avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti4, wdoti4avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti5, wdoti5avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti6, wdoti6avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdot, wdotavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)





*     --------------------------------------------
*     Integrate flux averaged quantities over rhon:
*     --------------------------------------------

      call rhograte(rhon, wdoteavg,  1, nnoderho, wdoteavg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti1avg, 1, nnoderho, wdoti1avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti2avg, 1, nnoderho, wdoti2avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti3avg, 1, nnoderho, wdoti3avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti4avg, 1, nnoderho, wdoti4avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti5avg, 1, nnoderho, wdoti5avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, wdoti6avg, 1, nnoderho, wdoti6avg_int,
     &   nrhomax, dvol)


c      call rhograte(rhon, wdote_ql,  1, nnoderho, wdote_ql_int,
c     &   nrhomax, dvol)
c      call rhograte(rhon, wdoti1_ql, 1, nnoderho, wdoti1_ql_int,
c     &   nrhomax, dvol)
c      call rhograte(rhon, wdoti2_ql, 1, nnoderho, wdoti2_ql_int,
c     &   nrhomax, dvol)
c      call rhograte(rhon, wdoti3_ql, 1, nnoderho, wdoti3_ql_int,
c     &   nrhomax, dvol)
c      call rhograte(rhon, wdoti4_ql, 1, nnoderho, wdoti4_ql_int,
c     &   nrhomax, dvol)
c      call rhograte(rhon, wdoti5_ql, 1, nnoderho, wdoti5_ql_int,
c     &   nrhomax, dvol)
c      call rhograte(rhon, wdoti6_ql, 1, nnoderho, wdoti6_ql_int,
c     &   nrhomax, dvol)


      call rhograte(rhon, redotjeavg, 1, nnoderho, redotjeavg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj1avg, 1, nnoderho, redotj1avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj2avg, 1, nnoderho, redotj2avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj3avg, 1, nnoderho, redotj3avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj4avg, 1, nnoderho, redotj4avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj5avg, 1, nnoderho, redotj5avg_int,
     &   nrhomax, dvol)
      call rhograte(rhon, redotj6avg, 1, nnoderho, redotj6avg_int,
     &   nrhomax, dvol)






31109    format(3x, "wdoteavg_int = ", 1p,e12.5, 1x, "Watts")
31110    format(3x, "wdote_ql_int = ", 1p,e12.5, 1x, "Watts")



41109    format(3x, "wdoti1avg_int = ", 1p,e12.5, 1x, "Watts")
41110    format(3x, "wdoti1_ql_int = ", 1p,e12.5, 1x, "Watts")


51109    format(3x, "wdoti2avg_int = ", 1p,e12.5, 1x, "Watts")
51110    format(3x, "wdoti2_ql_int = ", 1p,e12.5, 1x, "Watts")



61109    format(3x, "wdoti3avg_int = ", 1p,e12.5, 1x, "Watts")
61110    format(3x, "wdoti3_ql_int = ", 1p,e12.5, 1x, "Watts")



91109    format(3x, "wdoti4avg_int = ", 1p,e12.5, 1x, "Watts")
91110    format(3x, "fwdoti4_ql_int = ", 1p,e12.5, 1x, "Watts")




      call fluxavg(fyp, fypavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fype, fypeavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi1, fypi1avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi2, fypi2avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi3, fypi3avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi4, fypi4avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi5, fypi5avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi6, fypi6avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)




      call fluxavg(fz0, fz0avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0e, fz0eavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i1, fz0i1avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i2, fz0i2avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i3, fz0i3avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i4, fz0i4avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i5, fz0i5avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i6, fz0i6avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(capr_fzeta, capr_fzeta_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(capr2, capr2_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(rhom1, rhom1avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(muhat, muhat_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(nu_star, nu_star_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(qsafety, qsafety_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(pitch, pitch_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(ipsi, ipsi_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)



      call fluxavg(pressi, pressiavg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(psi_dim, psi_dim_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(fpsi1, fpsi1_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(ftheta1, ftheta1_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(bdotf, bdotf_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(gradprlb, gradprlb_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(bmod2, bmod2_avg, rho, nxmx, nymx, nrhomax,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

c      if (myid .eq. 0) then
c          write(6,*) "gradprlb_avg = ", gradprlb_avg
c      end if


      dpdpsi(1) = 0.0
      do n = 2, nnoderho - 1
         dpdpsi(n) = (pressiavg(n+1) - pressiavg(n-1))
     &              / (psi_dim_avg(n+1) - psi_dim_avg(n-1))
      end do

  300 format (1p,8e11.3)

*     --------------------------------------------
*     Integrate flux averaged quantities over rhon:
*     --------------------------------------------

      call rhograte(rhon, xjprlavg, 1, nnoderho, xjprl_int,
     &   nrhomax, darea)

      call rhograte(rhon, fypavg, 1, nnoderho, fyp_int,
     &   nrhomax, darea)

      call rhograte(rhon, fz0avg, 1, nnoderho, fz0_int,
     &   nrhomax, darea)

c      write(6, 300) (xjprl_int(n), n = 1, nnoderho)

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2843) tmin
         write(15, 2843) tmin
      endif

 2843 format('time to do flux averages =', f9.3, ' min')


      if(muhat_avg(1) .eq. 0.0) go to 9945
c      go to 9945


*     ----------------------------------------
*     Calculate force and flow velocities:
*     ----------------------------------------
      t1 = second1(dummy)


*     ----------------------------------------
*     Calculate coefficients for Khat equation
*     ----------------------------------------

      do n = 1, nnoderho
         if (rhom1avg(n) .gt. 1.0e-15)then
            gamma_avg(n) = muhat_avg(n)
     &                   * bmod2_avg(n) / (ipsi_avg(n) / rt)**2
         endif
         beta_avg(n) = (epsn(n) / qsafety_avg(n))**2
     &                              * (1. + 2.0 * qsafety_avg(n)**2)
      end do


*     -----------------------------
*     Solve for Khat with diffusion
*     -----------------------------

      ank(1) = 0.0
      bnk(1) = - gamma_avg(1) - 4.0 * dn(1) * beta_avg(1) / dr**2
      cnk(1) =   4.0 * dn(1) * beta_avg(2) / dr**2
      ynk(1) = capr_fzeta_avg(1) / rt - rt * bdotf_avg(1) / ipsi_avg(1)

      do n = 2, nnoderho - 1
         ank(n) = rh(n-1) * beta_avg(n-1) * dh(n-1) / (rn(n) * dr**2)
         cnk(n) = rh(n)   * beta_avg(n+1) * dh(n)   / (rn(n) * dr**2)
         bnk(n) = - gamma_avg(n)
     &            - rh(n-1) * beta_avg(n) * dh(n-1) / (rn(n) * dr**2)
     &            - rh(n)   * beta_avg(n) * dh(n)   / (rn(n) * dr**2)
         ynk(n) = capr_fzeta_avg(n) /rt - rt * bdotf_avg(n) /ipsi_avg(n)

      end do

      ank(nnoderho) = 0.0
      bnk(nnoderho) = 1.0
      cnk(nnoderho) = 0.0
      ynk(nnoderho) = 0.0

      call tridag(ank, bnk, cnk, ynk, xkhat, nnoderho)

      do n = 1, nnoderho
         kpsi_avg(n) = 0.0
         if (rhom1avg(n) .gt. 1.0e-10)then
            kpsi_avg(n) = xkhat(n) / (rhom1avg(n) * ipsi_avg(n) / rt)

         endif
      end do


      if(myid .eq. 0) write(6,  *)
      if(myid .eq. 0) write(15,  *)
c      if(myid .eq. 0) write(6,  *)"        n    rhon   pitch angle(deg)"
c      if(myid .eq. 0) write(15, *)"        n    rhon   pitch angle(deg)"

      if(myid .eq. 0)write(6, 76164)
      if(myid .eq. 0)write(15, 76164)
76164 format(" ",5x,"   n  ",3x," rhon ",
     &             6x," pitch",4x,"angle(deg)",4x,"      ",6x,"      ",
     &             6x,"      ",6x,"      ",6x,"      ",6x,"      ",
     &             6x,"      ")
      if(myid .eq. 0) write(6,  *)
      if(myid .eq. 0) write(15,  *)

      do n = 1, nnoderho - 1
         angle = atan2(pitch_avg(n), 1.0)
         angle = angle / pi * 180.0
         if(myid .eq. 0)write(6, 1312)n, rhon(n), pitch_avg(n), angle,
     &                       beta_avg(n), qsafety_avg(n), epsn(n)
         if(myid .eq. 0)write(15,1312)n, rhon(n), pitch_avg(n), angle,
     &                       beta_avg(n), qsafety_avg(n), epsn(n)
      end do




*     --------------
*     Solve for Jhat
*     --------------
      an(1) = 0.0
      bn(1) = - 4.0 * dn(1) / dr**2
      cn(1) = - bn(1)
      yn(1) = - capr_fzeta_avg(1) / rt

      do n = 2, nnoderho - 1
         an(n) = rh(n-1) * dh(n-1) / (rn(n) * dr**2)
         cn(n) = rh(n) * dh(n)     / (rn(n) * dr**2)
         bn(n) = - an(n) - cn(n)
         yn(n) = - capr_fzeta_avg(n) / rt
      end do

      an(nnoderho) = 0.0
      bn(nnoderho) = 1.0
      cn(nnoderho) = 0.0
      yn(nnoderho) = 0.0

      call tridag(an, bn, cn, yn, xjhat, nnoderho)



*     -------------------------
*     Calculate G(psi) and Epsi
*     -------------------------
      if(myid .eq. 0) write(6,  *)
      if(myid .eq. 0) write(15,  *)
      if(myid .eq. 0) write(6,  *)  "        n    xjhat      yn"
      if(myid .eq. 0) write(15, *)  "        n    xjhat      yn"

      do n = 1, nnoderho
         gpsi_avg(n) = 0.0
         epsi_avg(n) = 0.0

         if (rhom1avg(n) .gt. 1.0e-10)then

            gpsi_avg(n) = (xjhat(n) - rhom1avg(n) * kpsi_avg(n) * b0)
     &         / (rhom1avg(n) * rt * (1. + 1.5 * epsn(n)**2))

            term2 = rhom1avg(n) * kpsi_avg(n) * b0

            if(myid .eq. 0) then
               write(6,  1312)n, xjhat(n), yn(n)
               write(15, 1312)n, xjhat(n), yn(n)
            end if

            epsi_avg(n) = gpsi_avg(n)
c     &                 - fpsi1_avg(n)
c     &                 + 1. / (qi2 * xn2avg(n)) * dpdpsi(n)

         endif
      end do



*     --------------------
*     Calculate dE / dpsi
*     --------------------

      dedpsi_avg(1) = 0.0

      do n = 2, nnoderho - 1
         dedpsi_avg(n) = (epsi_avg(n+1) - epsi_avg(n-1))
     &                 / (psi_dim_avg(n+1) - psi_dim_avg(n-1))
      end do

      dedpsi_avg(nnoderho) = (epsi_avg(nnoderho) - epsi_avg(nnoderho-1))
     &               / (psi_dim_avg(nnoderho) - psi_dim_avg(nnoderho-1))



*     -----------------------------------
*     Put G, K, E, Jhat and dE/dpsi on 2D grid
*     -----------------------------------

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1

            gpsi(i,j) = 0.0
            kpsi(i,j) = 0.0
            epsi(i,j) = 0.0
            dedpsi(i,j) = 0.0
            jhat(i,j) = 0.0

            if(n .lt. nnoderho .and. rho(i,j) .le. rholim)then

               gpsi(i,j) = gpsi_avg(n) + (gpsi_avg(n + 1) - gpsi_avg(n))
     &             / drho * (rho(i,j) - rhon(n))
               kpsi(i,j) = kpsi_avg(n) + (kpsi_avg(n + 1) - kpsi_avg(n))
     &             / drho * (rho(i,j) - rhon(n))

               epsi(i,j) = epsi_avg(n) + (epsi_avg(n + 1) - epsi_avg(n))
     &             / drho * (rho(i,j) - rhon(n))
               dedpsi(i,j) = dedpsi_avg(n) + (dedpsi_avg(n + 1)
     &             - dedpsi_avg(n)) / drho * (rho(i,j) - rhon(n))
               jhat(i,j) = xjhat(n) + (xjhat(n + 1)
     &             - xjhat(n)) / drho * (rho(i,j) - rhon(n))

            end if

            if(n .eq. nnoderho .and. rho(i,j) .le. rholim)then

               gpsi(i,j) = gpsi_avg(n)
               kpsi(i,j) = kpsi_avg(n)
               epsi(i,j) = epsi_avg(n)
               dedpsi(i,j) = dedpsi_avg(n)
               jhat(i,j) = xjhat(n)

            end if
         end do
      end do


*     -----------------------------------------------
*     Calculate ExB shearing rate and flow velocities
*     -----------------------------------------------

      do i = 2, nnodex - 1
         do j = 1, nnodey

            omgexb(i,j) = (capr(i) * bpol(i,j))**2 / bmod(i,j)
     &                                           * dedpsi(i,j)
            uzeta(i,j) = gpsi(i,j) * capr(i) +
     &                      kpsi(i,j) * bzeta(i,j) * bmod(i,j)

c            uzeta(i,j) = jhat(i,j) / rhom1(i,j)

            utheta(i,j) = kpsi(i,j) * bpol(i,j)


         end do
      end do


      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2844) tmin
         write(15, 2844) tmin
      endif

 2844 format('time to calculate omg_ExB =', f9.3, ' min')

 9945 continue

      if (myid .eq. 0) then
         write(6,  1309) ptot
         write(15, 1309) ptot

         fedge = ptot_edge / ptot * 100.0
         fcore = ptot_core / ptot * 100.0

         f_wdot_edge = 0.0
         f_wdot_core = 0.0
         if (ptot_wdot .ne. 0.0) then
            f_wdot_edge = wdot_edge / ptot_wdot * 100.0
            f_wdot_core = wdot_core / ptot_wdot * 100.0
         end if


         write(6,  1110) pcrto2, pcito2
         write(15, 1110) pcrto2, pcito2

         write(6, *) "fraction absorbed in edge = ", fedge, " %"
         write(6, *) "fraction absorbed in core = ", fcore, " %"
         write(15, *) "fraction absorbed in edge = ", fedge, " %"
         write(15, *) "fraction absorbed in core = ", fcore, " %"

         write(6, *) "wdot fraction in edge = ", f_wdot_edge, " %"
         write(6, *) "wdot fraction in core = ", f_wdot_core, " %"
         write(15, *) "wdot fraction in edge = ", f_wdot_edge, " %"
         write(15, *) "wdot fraction in core = ", f_wdot_core, " %"


         write(6,  1216)powtot
         write(15, 1216)powtot

         write(6,  1217) xjtot
         write(15, 1217) xjtot

         ntilda_e_max = maxval(ntilda_e_real)
         write(6, *)  "ntilda_e_max = ", ntilda_e_max, "m-3"
         write(15, *) "ntilda_e_max = ", ntilda_e_max, "m-3"

         write(6, 169)
         write(6, 1120) fxtot_ant, fxtot_plasma
         write(6, 1121) fytot_ant, fytot_plasma
         write(6, 1122) fztot_ant, fztot_plasma

         write(6, 1123) fz0tot

         write(15, 169)
         write(15, 1120) fxtot_ant, fxtot_plasma
         write(15, 1121) fytot_ant, fytot_plasma
         write(15, 1122) fztot_ant, fztot_plasma

         write(15, 1123) fz0tot

      endif

*     ---------------------------------------------
*     Calculate species fractions from Estar dot J:
*     ---------------------------------------------
      pcedotje = pedotje / pedotjt * 100.
      pcedotj1 = pedotj1 / pedotjt * 100.
      pcedotj2 = pedotj2 / pedotjt * 100.
      pcedotj3 = pedotj3 / pedotjt * 100.
      pcedotj4 = pedotj4 / pedotjt * 100.
      pcedotj5 = pedotj5 / pedotjt * 100.
      pcedotj6 = pedotj6 / pedotjt * 100.
      pcedotjs = pedotjs / pedotjt * 100.

      pcedotjt = pcedotje + pcedotj1 + pcedotj2 + pcedotj3
     &         + pcedotj4 + pcedotj5 + pcedotj6 + pcedotjs


      pedotje = pedotje * twopi * rt
      pedotj1 = pedotj1 * twopi * rt
      pedotj2 = pedotj2 * twopi * rt
      pedotj3 = pedotj3 * twopi * rt
      pedotj4 = pedotj4 * twopi * rt
      pedotj5 = pedotj5 * twopi * rt
      pedotj6 = pedotj6 * twopi * rt
      pedotjs = pedotjs * twopi * rt
      pedotjt = pedotjt * twopi * rt

      i_ant = i_ant * twopi * rt


71667 format(' Species absorption from Estar dot J:')
      if (myid .eq. 0) then
         write(15, 169)
         write(15, 71667)
         write(15, 1109)  pedotje, pcedotje,
     &                    pedotj1, pcedotj1,
     &                    pedotj2, pcedotj2,
     &                    pedotj3, pcedotj3,
     &                    pedotj4, pcedotj4,
     &                    pedotj5, pcedotj5,
     &                    pedotj6, pcedotj6,
     &                    pedotjs, pcedotjs,
     &                    pedotjt, pcedotjt

         z_load = pedotjt / (i_ant**2)
         write(15, 169)
         write(15, 1148) i_ant
 1148 format(3x, "                  antenna current = ",1e12.5," Amps")
         write(15, 1149) z_load
 1149 format(3x, "                  antenna loading = ",1e12.5," Ohms")


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
     &   3x, "power absorbed by slowing species = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       "/
     &   3x, "             total power absorbed = ",1e12.5,
     &   " Watts ", " = ", f10.4, " %       ")

         write(6, 169)
         write(6, 71667)
         write(6, 1109)   pedotje, pcedotje,
     &                    pedotj1, pcedotj1,
     &                    pedotj2, pcedotj2,
     &                    pedotj3, pcedotj3,
     &                    pedotj4, pcedotj4,
     &                    pedotj5, pcedotj5,
     &                    pedotj6, pcedotj6,
     &                    pedotjs, pcedotjs,
     &                    pedotjt, pcedotjt

         write(6, 169)
         write(6, 1149) z_load


c--      Calculate species fractions from Wdot:
         if (pt .ne. 0.0) then
            pcte =   pe / pt  * 100.
            pcti1 =  pi1 / pt * 100.
            pcti2 =  pi2 / pt * 100.
            pcti3 =  pi3 / pt * 100.
            pcti4 =  pi4 / pt * 100.
            pcti5 =  pi5 / pt * 100.
            pcti6 =  pi6 / pt * 100.
         end if

         pctt = pcte + pcti1 + pcti2 + pcti3
     &               + pcti4 + pcti5 + pcti6


81667 format(' Species absorption from Wdot:')

         pe  = pe  * twopi * rt
         pi1 = pi1 * twopi * rt
         pi2 = pi2 * twopi * rt
         pi3 = pi3 * twopi * rt
         pi4 = pi4 * twopi * rt
         pi5 = pi5 * twopi * rt
         pi6 = pi6 * twopi * rt
         pt  = pt  * twopi * rt

         write(15, 169)
         write(15, 81667)
         write(15, 71109)  pe, pcte,
     &                    pi1, pcti1,
     &                    pi2, pcti2,
     &                    pi3, pcti3,
     &                    pi4, pcti4,
     &                    pi5, pcti5,
     &                    pi6, pcti6,
     &                    pt,  pctt
         write(6, 169)
         write(6, 81667)
         write(6, 71109)  pe,  pcte,
     &                    pi1, pcti1,
     &                    pi2, pcti2,
     &                    pi3, pcti3,
     &                    pi4, pcti4,
     &                    pi5, pcti5,
     &                    pi6, pcti6,
     &                    pt,  pctt

         z_load = pt / (i_ant**2)
         write(15, 169)
         write(15, 1149) z_load


81668 format('Species absorption from quasilinear Wdot:')

c         pe_ql  =  wdote_ql_int(nnoderho)
c        pi1_ql = wdoti1_ql_int(nnoderho)
c        pi2_ql = wdoti2_ql_int(nnoderho)
c        pi3_ql = wdoti3_ql_int(nnoderho)
c        pi4_ql = wdoti4_ql_int(nnoderho)
c        pi5_ql = wdoti5_ql_int(nnoderho)
c        pi6_ql = wdoti6_ql_int(nnoderho)

c        pt_ql  = pe_ql + pi1_ql + pi2_ql + pi3_ql
c     &                   + pi4_ql +pi5_ql + pi6_ql


         if (pt_ql .ne. 0.0) then
            pcte_ql =   pe_ql / pt_ql  * 100.
            pcti1_ql =  pi1_ql / pt_ql * 100.
            pcti2_ql =  pi2_ql / pt_ql * 100.
            pcti3_ql =  pi3_ql / pt_ql * 100.
            pcti4_ql =  pi4_ql / pt_ql * 100.
            pcti5_ql =  pi5_ql / pt_ql * 100.
            pcti6_ql =  pi6_ql / pt_ql * 100.
         end if

         pctt_ql = pcte_ql + pcti1_ql + pcti2_ql + pcti3_ql
     &                     + pcti4_ql + pcti5_ql + pcti6_ql

         write(15, 169)
         write(15, 81668)
         write(15, 71109)  pe_ql, pcte_ql,
     &                    pi1_ql, pcti1_ql,
     &                    pi2_ql, pcti2_ql,
     &                    pi3_ql, pcti3_ql,
     &                    pi4_ql, pcti4_ql,
     &                    pi5_ql, pcti5_ql,
     &                    pi6_ql, pcti6_ql,
     &                    pt_ql,  pctt_ql
         write(6, 169)
         write(6, 81668)
         write(6, 71109)  pe_ql,  pcte_ql,
     &                    pi1_ql, pcti1_ql,
     &                    pi2_ql, pcti2_ql,
     &                    pi3_ql, pcti3_ql,
     &                    pi4_ql, pcti4_ql,
     &                    pi5_ql, pcti5_ql,
     &                    pi6_ql, pcti6_ql,
     &                    pt_ql,  pctt_ql


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

*        ---------------------------------------------------
*        Calculate the differential volume elements, dvol(n):
*        ---------------------------------------------------

c         dvol  = 0.0
c         darea = 0.0

c         do i = 1, nnodex
c            do j = 1, nnodey
c               n = int(rho(i,j) / drho) + 1
c               if(n .lt. nnoderho)then
c                  dvol(n)  =  dvol(n) + dx * dy * 2.0 * pi * capr(i)
c                 darea(n) = darea(n) + dx * dy * capr(i) / rt
c               end if
c            end do
c         end do


         write(38, 309) nmax, mmax
         write(38, 310) (zetai_table(n), n = 1, nmax)
         write(38, 310) (dKdL_table(m), m = 1, mmax)
         write(38, 310) ((z0_table1(n, m), n = 1, nmax),  m = 1, mmax)
         write(38, 310) ((z1_table1(n, m), n = 1, nmax),  m = 1, mmax)
         write(38, 310) ((z2_table1(n, m), n = 1, nmax),  m = 1, mmax)


         write(38, 309) nnodex, nnodey, jmid, idiag, jdiag
         write(38, 310) psilim, prfin, omgrf
         write(38, 310) (x(i), i = 1, nnodex)
         write(38, 310) (y(j), j = 1, nnodey)
         write(38, 310) (capr(i), i = 1, nnodex)
         write(38, 310) (xkphi(i), i = 1, nnodex)



         write(130, 309) nnodex, nnodey, jmid
         write(130, 310) psilim
         write(130, 310) (capr(i), i = 1, nnodex)
         write(130, 310) (y(j), j = 1, nnodey)
         write(130, 310) ((xkperp2_slow(i, j), i = 1, nnodex),
     &      j = 1, nnodey)
         write(130, 310) ((xkperp2_fast(i, j), i = 1, nnodex),
     &      j = 1, nnodey)
         write(130, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
         write(130, 310) ((xkprl_a(i,j), i = 1, nnodex), j = 1, nnodey)
         write(130, 310) ((P_a(i, j), i = 1, nnodex), j = 1, nnodey)



         write(38, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((theta(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((xnea(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((xkte(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((xkti(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((xkperp_cold(i,j), i = 1, nnodex),j = 1,nnodey)
         write(38, 310) ((xkperp_cold2(i,j), i = 1,nnodex),j = 1,nnodey)

         write(38, 310) ((acold(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((bcold(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((ccold(i, j), i = 1, nnodex), j = 1, nnodey)

         write (38, 310) ((reomg1a(i, j), i = 1, nnodex), j = 1, nnodey)
         write (38, 310) ((reomg2a(i, j), i = 1, nnodex), j = 1, nnodey)
         write (38, 310) ((reomg3a(i, j), i = 1, nnodex), j = 1, nnodey)

         write (38, 310) ((reomglha(i, j), i = 1, nnodex),j = 1, nnodey)

         write(38, 310) ((xjx(i, j),  i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((xjy(i, j),  i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((xjz(i, j),  i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)

         write (38, 310) ((xjpxe(i,j), i = 1, nnodex), j = 1, nnodey)
         write (38, 310) ((xjpye(i,j), i = 1, nnodex), j = 1, nnodey)
         write (38, 310) ((xjpze(i,j), i = 1, nnodex), j = 1, nnodey)

         write(38, 309) nkx1, nkx2
         write(38, 309) nky1, nky2

         write(38, 310) (xkxsav(n), n = nkx1, nkx2)
         write(38, 310) (xkysav(m), m = nky1, nky2)

         write(38, 310) ((ealphakmod(n, m), n = nkx1, nkx2),
     &       m = nky1, nky2)
         write(38, 310) ((ebetakmod(n, m), n = nkx1, nkx2),
     &       m = nky1, nky2)
         write(38, 310) ((ebkmod(n, m), n = nkx1, nkx2),
     &       m = nky1, nky2)

         write(38, 310) ((redotje(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj2(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj3(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotjt(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((bmod(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((qsafety(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((ebeta(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((btau(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((bzeta(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((spx(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((spy(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((spz(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((ealphak(n, m), n = nkx1, nkx2), m= nky1, nky2)
         write(38, 310) ((ebetak(n, m), n = nkx1, nkx2), m = nky1, nky2)
         write(38, 310) ((ebk(n, m), n = nkx1, nkx2), m = nky1, nky2)

         write(38, 310) ((wdote(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti2(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdot(i, j), i = 1, nnodex), j = 1, nnodey)

         write (38, 310) ((divq(i, j), i = 1, nnodex), j = 1, nnodey)
         write (38, 310) ((capd(i, n), i = 1, nnodex), n = nkx1, nkx2)
         write (38, 310) ((xkb(i, n), i = 1, nnodex), n = nkx1, nkx2)


         do n = nkx1, nkx2
            write(6, *) "n = ", n, "xkb(64, n) = ", xkb(64, n)
         end do

         write(38, 309) nnoderho
         write(38, 310) (rhon(n), n = 1, nnoderho)
         write(38, 310) (gradprlb_avg(n), n = 1, nnoderho)

         write(38, 309) mnodetheta
         write(38, 310) (thetam(m), m = 1, mnodetheta)

         write(38, 310) (dvol(n), n = 1, nnoderho)
         write(38, 310) (volume(n), n = 1, nnoderho)

         write(38, 310) (xnavg(n), n = 1, nnoderho)
         write(38, 310) (wdoteavg(n), n = 1, nnoderho)
         write(38, 310) (wdoti1avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti2avg(n), n = 1, nnoderho)
         write(38, 310) (vyavg(n), n = 1, nnoderho)
         write(38, 310) (dvydrho(n), n = 1, nnoderho)

         write(38, 310) (redotjeavg(n), n = 1, nnoderho)
         write(38, 310) (redotj1avg(n), n = 1, nnoderho)
         write(38, 310) (redotj2avg(n), n = 1, nnoderho)
         write(38, 310) (redotj3avg(n), n = 1, nnoderho)
         write(38, 310) (redotjtavg(n), n = 1, nnoderho)

         write(38, 310) ((redotjs(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti3(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) (wdoti3avg(n), n = 1, nnoderho)

         write(38, 310) ((fz0e(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i2(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i3(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (fz0eavg(n), n = 1, nnoderho)
         write(38, 310) (fz0i1avg(n), n = 1, nnoderho)
         write(38, 310) (fz0i2avg(n), n = 1, nnoderho)
         write(38, 310) (fz0i3avg(n), n = 1, nnoderho)
         write(38, 310) (fz0avg(n), n = 1, nnoderho)

         write(38, 310) (xjprlavg(n), n = 1, nnoderho)

         write(38, 310) (xjprl_int(n), n = 1, nnoderho)
         write(38, 310) (fz0_int(n), n = 1, nnoderho)

         write(38, 310) (xjhat(n), n = 1, nnoderho)
         write(38, 310) (xkhat(n), n = 1, nnoderho)

         write(38, 310) (gpsi_avg(n), n = 1, nnoderho)
         write(38, 310) ((gpsi(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (kpsi_avg(n), n = 1, nnoderho)
         write(38, 310) ((kpsi(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (muhat_avg(n), n = 1, nnoderho)
         write(38, 310) ((muhat(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (nu_star_avg(n), n = 1, nnoderho)
         write(38, 310) ((nu_star(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (ipsi_avg(n), n = 1, nnoderho)
         write(38, 310) ((ipsi(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((omgexb(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((uzeta(i,j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((utheta(i,j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (xn1avg(n), n = 1, nnoderho)
         write(38, 310) (xn2avg(n), n = 1, nnoderho)
         write(38, 310) (xn3avg(n), n = 1, nnoderho)
         write(38, 310) (xna_sloavg(n), n = 1, nnoderho)
         write(38, 310) (xkteavg(n), n = 1, nnoderho)
         write(38, 310) (xktiavg(n), n = 1, nnoderho)
         write(38, 310) (xkti2avg(n), n = 1, nnoderho)
         write(38, 310) (xkti3avg(n), n = 1, nnoderho)

         write(38, 310) (redotjsavg(n), n = 1, nnoderho)

         write(38, 310) ((bmod_mid(i,j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) (bmod_midavg(n), n = 1, nnoderho)

         write(38, 310) ((capr_bpol_mid2(i, j), i = 1,nnodex),
     &                                                    j = 1,nnodey)
         write(38, 310) (capr_bpol_mid(n), n = 1, nnoderho)

         write(38, 310) ((redotj4(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj5(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj6(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (redotj4avg(n), n = 1, nnoderho)
         write(38, 310) (redotj5avg(n), n = 1, nnoderho)
         write(38, 310) (redotj6avg(n), n = 1, nnoderho)

         write(38, 310) ((wdoti4(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti5(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti6(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (wdoti4avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti5avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti6avg(n), n = 1, nnoderho)

         write(38, 310) ((fz0i4(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i5(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i6(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (fz0i4avg(n), n = 1, nnoderho)
         write(38, 310) (fz0i5avg(n), n = 1, nnoderho)
         write(38, 310) (fz0i6avg(n), n = 1, nnoderho)

         write(38, 310) (xn4avg(n), n = 1, nnoderho)
         write(38, 310) (xn5avg(n), n = 1, nnoderho)
         write(38, 310) (xn6avg(n), n = 1, nnoderho)

         write(38, 310) (xkti4avg(n), n = 1, nnoderho)
         write(38, 310) (xkti5avg(n), n = 1, nnoderho)
         write(38, 310) (xkti6avg(n), n = 1, nnoderho)

         write(38, 310) (wdoteavg_int(n),  n = 1, nnoderho)
         write(38, 310) (wdoti1avg_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti2avg_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti3avg_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti4avg_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti5avg_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti6avg_int(n), n = 1, nnoderho)

         write(38, 310) (wdote_ql_int(n),  n = 1, nnoderho)
         write(38, 310) (wdoti1_ql_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti2_ql_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti3_ql_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti4_ql_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti5_ql_int(n), n = 1, nnoderho)
         write(38, 310) (wdoti6_ql_int(n), n = 1, nnoderho)

         write(38, 310) (redotjeavg_int(n), n = 1, nnoderho)
         write(38, 310) (redotj1avg_int(n), n = 1, nnoderho)
         write(38, 310) (redotj2avg_int(n), n = 1, nnoderho)
         write(38, 310) (redotj3avg_int(n), n = 1, nnoderho)
         write(38, 310) (redotj4avg_int(n), n = 1, nnoderho)
         write(38, 310) (redotj5avg_int(n), n = 1, nnoderho)
         write(38, 310) (redotj6avg_int(n), n = 1, nnoderho)


         write(59,309) nnodex, nnodey
         write(59,310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)

         write(60,309) nnodex, nnodey
         write(60,310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)

         write(51,309) nnodex, nnodey
         write(51,310) psilim
         write(51,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)

         write(69,309) nnodex, nnodey
         write(69,310) ((wdote(i, j), i = 1, nnodex), j = 1, nnodey)
         write(69,310) ((wdoti1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(69,310) ((wdoti2(i, j), i = 1, nnodex), j = 1, nnodey)
         write(69,310) ((wdoti3(i, j), i = 1, nnodex), j = 1, nnodey)


*        -------------------------------
*        write edotj data to edotj.out
*        -------------------------------

         open(unit=39, file='out_edotj', status='unknown',
     &      form='formatted')
         rewind(39)

         write(39, 1314) nnodex, nnodey
         write(39, 310) ((redotje(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotji(i, j), i = 1, nnodex), j = 1, nnodey)

         write(39, 310) ((redotj1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj2(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj3(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj4(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj5(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj6(i, j), i = 1, nnodex), j = 1, nnodey)

         close (39)

*     -----------------------------------------------
*     Write file 46, "naoto_2d" for synthetic diagnostic
*     -----------------------------------------------
      open(unit=46, file='naoto_2d', status='unknown',
     &                                              form= 'formatted')
      write (46, 309) nnodex, nnodey, nphi
      write (46, 310) (capr(i),  i = 1, nnodex)
      write (46, 310) (y(j),  j = 1, nnodey)


      write (46, 310) ((ealpha(i, j), i = 1, nnodex),
     &                                j = 1, nnodey)
      write (46, 310) ((ebeta(i, j), i = 1, nnodex),
     &                               j = 1, nnodey)
      write (46, 310) ((eb(i, j), i = 1, nnodex),
     &                            j = 1, nnodey)


      write (46, 310) ((xjpxe(i,j) , i = 1, nnodex),
     &                                   j = 1, nnodey)
      write (46, 310) ((xjpye(i,j) , i = 1, nnodex),
     &                                   j = 1, nnodey)
      write (46, 310) ((xjpze(i,j) , i = 1, nnodex),
     &                                   j = 1, nnodey)


      write (46, 310) ((xjpxe_lab(i,j) , i = 1, nnodex),
     &                                   j = 1, nnodey)
      write (46, 310) ((xjpye_lab(i,j) , i = 1, nnodex),
     &                                   j = 1, nnodey)
      write (46, 310) ((xjpze_lab(i,j) , i = 1, nnodex),
     &                                   j = 1, nnodey)


      write (46, 310) ((ntilda_e(i, j), i = 1, nnodex),
     &                                   j = 1, nnodey)


      close (46)


      endif


      t1 = second1(dummy)

      nxdim = nxmx
      nydim = nymx

      call run_rf2x(nmodesx, nmodesy, rwleft, rwright,
     &   ytop, ybottom, myid, nxdim, nydim,
     &   rt, b0, rho,
     &   redotje, redotji,
     &   redotj1, redotj2,
     &   redotj3, redotj4,
     &   redotj5, redotj6,
     &   xjprl, wdote,
     &   wdoti1, wdoti2,
     &   wdoti3, wdoti4,
     &   wdoti5, wdoti6, rhomax)

      if (myid.eq.0) then

         write(38, 310) (wdote_ql(n),  n = 1, nnoderho)
         write(38, 310) (wdoti1_ql(n), n = 1, nnoderho)
         write(38, 310) (wdoti2_ql(n), n = 1, nnoderho)
         write(38, 310) (wdoti3_ql(n), n = 1, nnoderho)
         write(38, 310) (wdoti4_ql(n), n = 1, nnoderho)
         write(38, 310) (wdoti5_ql(n), n = 1, nnoderho)
         write(38, 310) (wdoti6_ql(n), n = 1, nnoderho)

         write(38, 310) (dldbavg(n), n = 1, nnoderho)

         write(38, 310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((ntilda_e(i, j), i = 1, nnodex), j = 1,nnodey)

         write(38, 310) ((fpsi0(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((ftheta0(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((pressi(i,j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((dldb_tot12(i, j), i = 1, nnodex),
     &                                      j = 1, nnodey)

         write (38, 309) ndisti1, ndisti2, ndisti3

         write (38, 310) ((eplus_flux_plot(i, j), i = 1, nnodex),
     &                                             j = 1, nnodey)

         write (38, 310) ((eminus_flux_plot(i, j), i = 1, nnodex),
     &                                             j = 1, nnodey)

         write (38, 310) ((xkperp_flux_plot(i, j), i = 1, nnodex),
     &                                             j = 1, nnodey)




         write (38, 310) ((eplus_flux(n, m), n = 1, nnoderho),
     &                                       m = 1, mnodetheta)

         write (38, 310) ((eminus_flux(n, m), n = 1, nnoderho),
     &                                        m = 1, mnodetheta)

         write (38, 310) ((xkperp_flux(n, m), n = 1, nnoderho),
     &                                        m = 1, mnodetheta)

         write (38, 310) ((capr_flux(n,m), n = 1, nnoderho),
     &                                       m = 1, mnodetheta)
         write (38, 310) ((capz_flux(n,m), n = 1, nnoderho),
     &                                       m = 1, mnodetheta)
         write (38, 310) ((bmod_flux(n,m), n = 1, nnoderho),
     &                                       m = 1, mnodetheta)

         write(38, 310) ((dxxuyy(i,j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((dyyuzz(i,j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((gradprlb(i,j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((xkperp2_slow(i, j), i = 1, nnodex),
     &      j = 1, nnodey)
         write(38, 310) ((xkperp2_fast(i, j), i = 1, nnodex),
     &      j = 1, nnodey)
         write(38, 310) ((xkprl_a(i,j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((P_a(i, j), i = 1, nnodex), j = 1, nnodey)


         close (38)

         close (130)



      end if

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2845) tmin
         write(15, 2845) tmin
      endif

 2845 format('time to call run_rf2x =', f9.3, ' min')

      call blacs_barrier(icontxt, 'All')

*     -------------------
*     write fpm (34) file
*     -------------------
      if (myid.eq.0) then

         if(nphi .eq. nphi1 .or. nt .eq. 1)then

            write(34,309) nmodesx, nmodesy
            write(34,310) rwleft, rwright, ytop, ybottom

            write(34,309) nnodex, nnodey
            write(34,310) (capr(i), i = 1, nnodex)
            write(34,310) (y(j), j = 1, nnodey)

            write(166, 309) nnodex, nnodey
            write(166, 310) (capr(i), i = 1, nnodex)
            write(166, 310) (y(j), j = 1, nnodey)
            write(166, 310) psilim
            write(166, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)

            write(166, 309) nphi_number
            write(166, 309) (nphi_array(n), n = 1, nphi_number)

            write(34,310) psilim, rt, b0
            write(34,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
            write(34,310) ((xnea(i, j), i = 1, nnodex), j = 1, nnodey)

            write (34, 309) nuper
            write (34, 309) nupar
            write (34, 309) nnoderho

            if(ndisti1 .eq. 1)write (34, 3310) vc1_cgs
            if(ndisti2 .eq. 1)write (34, 3310) vc2_cgs

            write (34, 3310) UminPara, UmaxPara

            write (34, 3310) (rhon(n), n = 1, nnoderho)
            write (34, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (34, 3310) (upara(i_upara), i_upara = 1, nupar)

         end if

         write(34,309) nphi

         write(34,310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34,310) ((ebeta(i, j),  i = 1, nnodex), j = 1, nnodey)
         write(34,310) ((eb(i, j),     i = 1, nnodex), j = 1, nnodey)

         write(166,310) ((eplus(i, j), i = 1, nnodex), j = 1, nnodey)
         write(166,310) ((eminus(i, j),  i = 1, nnodex), j = 1, nnodey)
         write(166,310) ((eb(i, j),     i = 1, nnodex), j = 1, nnodey)

         write(166,310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
         write(166,310) ((ey(i, j),  i = 1, nnodex), j = 1, nnodey)
         write(166,310) ((ez(i, j),     i = 1, nnodex), j = 1, nnodey)

         write(166,310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
         write(166,310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
         write(166,310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

         write(166,310) ((xjpx(i, j), i = 1, nnodex), j = 1, nnodey)
         write(166,310) ((xjpy(i, j), i = 1, nnodex), j = 1, nnodey)
         write(166,310) ((xjpz(i, j), i = 1, nnodex), j = 1, nnodey)

         write(34,310) ((xjpx(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34,310) ((xjpy(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34,310) ((xjpz(i, j), i = 1, nnodex), j = 1, nnodey)

         write(34,310) ((xjx(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34,310) ((xjy(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34,310) ((xjz(i, j), i = 1, nnodex), j = 1, nnodey)

         write(34,310) ((xjpxe_lab(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34,310) ((xjpye_lab(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34,310) ((xjpze_lab(i, j), i = 1, nnodex), j = 1, nnodey)

         write(34, 310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34, 310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34, 310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

         write(34, 310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34, 310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
         write(34, 310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)

         write(34,310) ((ntilda_e(i, j), i = 1, nnodex), j = 1, nnodey)

         write(34,310) ptot, pcrto2, pcito2, xjtot, spa_fast

         write(34, 309) nnoderho
         write(34, 310) (rhon(n), n = 1, nnoderho)

         write(34, 310) (redotjeavg(n), n = 1, nnoderho)
         write(34, 310) (redotj1avg(n), n = 1, nnoderho)
         write(34, 310) (redotj2avg(n), n = 1, nnoderho)
         write(34, 310) (redotj3avg(n), n = 1, nnoderho)
         write(34, 310) (redotj4avg(n), n = 1, nnoderho)
         write(34, 310) (redotj5avg(n), n = 1, nnoderho)
         write(34, 310) (redotj6avg(n), n = 1, nnoderho)

         write(34, 310) (wdoteavg(n),  n = 1, nnoderho)
         write(34, 310) (wdoti1avg(n), n = 1, nnoderho)
         write(34, 310) (wdoti2avg(n), n = 1, nnoderho)
         write(34, 310) (wdoti3avg(n), n = 1, nnoderho)
         write(34, 310) (wdoti4avg(n), n = 1, nnoderho)
         write(34, 310) (wdoti5avg(n), n = 1, nnoderho)
         write(34, 310) (wdoti6avg(n), n = 1, nnoderho)

         write(34, 310) (wdote_ql(n),  n = 1, nnoderho)
         write(34, 310) (wdoti1_ql(n), n = 1, nnoderho)
         write(34, 310) (wdoti2_ql(n), n = 1, nnoderho)
         write(34, 310) (wdoti3_ql(n), n = 1, nnoderho)
         write(34, 310) (wdoti4_ql(n), n = 1, nnoderho)
         write(34, 310) (wdoti5_ql(n), n = 1, nnoderho)

         write(34, 310) (xjprlavg(n),  n = 1, nnoderho)



         if(ndisti1 .eq. 1)then
            write (34, 3310) (((bqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((cqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((eqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((fqlavg_i1(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) xmi1
         end if


         if(ndisti2 .eq. 1)then
            write (34, 3310) (((bqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((cqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((eqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((fqlavg_i2(i_uperp, i_upara, n),
     &        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) xmi2
         end if

         write(6,*) "finished writing fpm; starting plots"

      end if

 9002 continue


      call blacs_barrier(icontxt, 'All')

*     -----------------------
*     do plotting with pgplot
*     -----------------------
      if (myid.eq.0) then
         t1 = second1(dummy)

         call fieldws( REAL(dfquotient, real32),REAL(rmin_zoom, real32),
     &                 REAL(rmax_zoom, real32) )

         call disp_ono( REAL(yplot, real32) )

         tmin = (second1(dummy) - t1) / 60.
         write(6 , 2846) tmin
         write(15, 2846) tmin
      end if

 2846 format('time to do plots =', f9.3, ' min')
 2847 format('time to sum modes =', f9.3, ' min')

c      write(6,*) "finished plots; starting deallocation"
      call blacs_barrier(icontxt, 'All')


  310 format(1p,6e12.4)
 3310 format(1p,6e18.10)
 8310 format(1p,6e14.6)
  309 format(10i10)
 9310 format(1p,7e12.4)
  311 format(1p,10e12.4)
 3117 format(1p,11e12.4)
 1312 format(i10, 1p,10e12.4)
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
 1021 format(3x,"          npar(R0) = ",1p,e12.4,"       ")
 1022 format(3x,"     npar(antenna) = ",1p,e12.4,"       ")
 1812 format(3x,"                rt = ",1p,e12.4," m     ")
 1882 format(3x,"           rt_init = ",1p,e12.4," m     ")
 1822 format(3x,"            aplasm = ",1p,e12.4," m     ")
 1823 format(3x,"              rant = ",1p,e12.4," m     ")
 1809 format(3x,"                b0 = ",1p,e12.4," T     ")
 1889 format(3x,"           b0_init = ",1p,e12.4," T     ")
 1813 format(3x,"              xn10 = ",1p,e12.4," m-3   ")
 6813 format(3x,"              xne0 = ",1p,e12.4," m-3   ")
 1814 format(3x,"              xn20 = ",1p,e12.4," m-3   ")
 1834 format(3x,"              xn30 = ",1p,e12.4," m-3   ")

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
17014 format(3x,"             lmaxe = ",i12,"       ")
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
71160 format(3x,"        xnue / omg = ",1p,e12.4,"       ")
71161 format(3x,"        xnui / omg = ",1p,e12.4,"       ")

 1309 format(3x,
     1   "             total power absorbed = ",1p,e12.4," watts/m ")
11091 format(3x,
     1   "       total power absorbed (old) = ",1p,e12.4," watts/m ")





81109 format(
     &   3x, "      power absorbed by electrons = ",1e12.5,
     &   " watts/m "/
     1   3x, "  power absorbed by majority ions = ",1e12.5,
     1   " watts/m "/
     1   3x, "  power absorbed by minority ions = ",1e12.5,
     1   " watts/m ")

 1112 format(
     1   3x,"      power absorbed by electrons = ",1f12.4," %      " /
     1   3x,"  power absorbed by majority ions = ",1f12.4," %      " /
     1   3x,"  power absorbed by minority ions = ",1f12.4," %      " /
     1   3x,"power absorbed by 3rd ion species = ",1f12.4," %      " /
     1   3x,"             total power absorbed = ",1f12.4," %      " )
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
     1   3x,"                         x(iedge) = ",1p,e12.4," m"/
     1   3x,"                               bz = ",1p,e12.4," T"/
     1   3x,"              real(eps parallel)  = ",1p,e12.4,"  "/
     1   3x,"    estimate: real(eps parallel)  = ",1p,e12.4,"  "/
     1   3x,"                 real :: (eps left)  = ",1p,e12.4,"  "/
     1   3x,"     estimate:   real :: (eps left)  = ",1p,e12.4,"  "/
     1   3x,"                real :: (eps right)  = ",1p,e12.4,"  "/
     1   3x,"    estimate:   real :: (eps right)  = ",1p,e12.4,"  "/
     1   3x,"                               xn = ",1p,e12.4," m-3"/
     1   3x,"                               LN = ",1p,e12.4," m"/
     1   3x,"                        mod2 (ez) = ",1p,e12.4," (V/m)**2"/
     1   3x,"                            L rf  = ",1p,e12.4," m"/
     1   3x,"                x=omgrf/omgci(i)  = ",1p,e12.4,"    "/
     1   3x,"                             xkti = ",1p,e12.4," deg"/
     1   3x,"                             xkte = ",1p,e12.4," deg")
 1389 format(
     1   3x,"                           dlnndr = ",1p,e12.4," m-1"/
     1   3x,"                            deprl = ",1p,e12.4,
     1   " (V/m)**2/m"/
     1   3x,"                           deleft = ",1p,e12.4,
     1   " (V/m)**2/m"/
     1   3x,"                           derght = ",1p,e12.4,
     1   " (V/m)**2/m"/
     1   3x,"                        real :: (ez) = ",1p,e12.4," V/m"/
     1   3x,"                         imag(ez) = ",1p,e12.4," V/m")
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
     1   3x,"               real power emitted = ",1p,e12.4," watts/m "/
     1   3x,"          imaginary power emitted = ",1p,e12.4," watts/m ")
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

 1123 format(3x,'total toroidal force on plasma = ',1p,e12.4,' Nt/m'/)


  163 format("0")

 2213 format(1i10, 1p,5e12.4)
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


c  164 format(" ",1x," i  ",  3x," R(x) ",
c     &           6x,"  x   ",6x," btau ",6x," bmod ",6x,"re om1",
c     &           6x,"re om2",6x,"re om3",6x,"re om4",6x,"re om5",
c     &           6x," xne  ",6x," Acold")

  164 format(1x ,1x," i  ",  3x," R(x) ",
     &           6x,"  x   ",6x," btau ",6x," bmod ",6x,"re om1",
     &           6x,"re om2",6x,"re om3",6x,"re omgpi",4x,"re omglh",
     &           4x," xne  ",6x," Acold")

c 9164 format(" ",1x," j  ",  3x," R(x) ",
c     &           6x,"  y   ",6x," btau ",6x," bmod ",6x,"re om1",
c     &           6x,"re om2",6x,"re om3",6x,"re om4",6x,"re om5",
c     &           6x," xne  ",6x," Acold")

 9164 format(" ",1x," j  ",  3x," R(x) ",
     &           6x,"  y   ",6x," btau ",6x," bmod ",6x,"re om1",
     &           6x,"re om2",6x,"re om3",6x,"re omgpi",4x,"re omglh",
     &           4x," xne  ",6x," Acold")



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
 1119 format(2x,"ellipticity = ",1p,e12.4,"      ")

 5555 continue

      call blacs_barrier(icontxt, 'All')

      time=second1(dummy)-time0

      ttotal = time/60.

      if (myid .eq. 0) then
         write(15,162)
         write(6, 162)
         write(15,899) ttotal
         write(6, 899) ttotal
c        close (15)
      endif

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
     &                         7x, 'Re Ey', 7x, 'Im Ey',
     &                         7x, 'Re Ez', 7x, 'Im Ez')

 2313 format(1p,8e12.4)
 6834 format(3x,"              eta1 = ",1p,e12.4,"       ")
 6835 format(3x,"              eta2 = ",1p,e12.4,"       ")
 6836 format(3x,"              eta3 = ",1p,e12.4,"       ")



!efd-begin
c      if (myid.eq.0) then
c      call profstat()
c      endif
!efd-end


      call blacs_barrier(icontxt, 'All')


      close(5)


      if (myid.eq.0) then

         close(29)
         close(30)

         close(59)
         close(60)
         close(51)
         close(69)

         close(28)

         close(40)
         close(41)
         close(42)
         close(43)
         close(44)

      endif

!     -----------------
!     deallocate arrays
!     -----------------

      if(ndisti1   .eq. 0 .and.
     &   ndisti2   .eq. 0 .and.
     &   ndisti3   .eq. 0 .and.
     &   ndisti4   .eq. 0 .and.
     &   ndisti5   .eq. 0 .and.
     &   ndisti6   .eq. 0 .and.
     &   ndiste    .eq. 0)   then

         deallocate( UPERP )
         deallocate( UPARA )
         deallocate( VPERP )
         deallocate( VPARA )
         deallocate( UPERP_work )
         deallocate( UPARA_work )

         deallocate( f )

         deallocate( dfdupere )
         deallocate( dfdupare )
         deallocate( fperpe )


         deallocate( dfduper1 )
         deallocate( dfdupar1 )
         deallocate( fperp1 )

         deallocate( dfduper2 )
         deallocate( dfdupar2 )
         deallocate( fperp2 )

         deallocate( dfduper3 )
         deallocate( dfdupar3 )
         deallocate( fperp3 )

         deallocate( dfduper4 )
         deallocate( dfdupar4 )
         deallocate( fperp4 )

         deallocate( dfduper5 )
         deallocate( dfdupar5 )
         deallocate( fperp5 )

         deallocate( dfduper6 )
         deallocate( dfdupar6 )
         deallocate( fperp6 )

         deallocate( bqlavg_e )
         deallocate( cqlavg_e )
         deallocate( eqlavg_e )
         deallocate( fqlavg_e )

         deallocate( bqlavg_work )
         deallocate( cqlavg_work )
         deallocate( eqlavg_work )
         deallocate( fqlavg_work )

         deallocate( bqlavg_i1 )
         deallocate( cqlavg_i1 )
         deallocate( eqlavg_i1 )
         deallocate( fqlavg_i1 )

         if(eta2 .ne. 0.0) then
            deallocate( bqlavg_i2 )
            deallocate( cqlavg_i2 )
            deallocate( eqlavg_i2 )
            deallocate( fqlavg_i2 )
         end if

         if(eta3 .ne. 0.0) then
            deallocate( bqlavg_i3 )
            deallocate( cqlavg_i3 )
            deallocate( eqlavg_i3 )
            deallocate( fqlavg_i3 )
         end if

         if(eta4 .ne. 0.0) then
            deallocate( bqlavg_i4 )
            deallocate( cqlavg_i4 )
            deallocate( eqlavg_i4 )
            deallocate( fqlavg_i4 )
         end if

         if(eta5 .ne. 0.0) then
            deallocate( bqlavg_i5 )
            deallocate( cqlavg_i5 )
            deallocate( eqlavg_i5 )
            deallocate( fqlavg_i5 )
         end if

         if(eta6 .ne. 0.0) then
            deallocate( bqlavg_i6 )
            deallocate( cqlavg_i6 )
            deallocate( eqlavg_i6 )
            deallocate( fqlavg_i6 )
         end if

      end if

      if(ndisti1   .ne. 0 .or.
     &   ndisti2   .ne. 0 .or.
     &   ndisti3   .ne. 0 .or.
     &   ndisti4   .ne. 0 .or.
     &   ndisti5   .ne. 0 .or.
     &   ndisti6   .ne. 0 .or.
     &   ndiste    .ne. 0)   then

         deallocate( UPERP )
         deallocate( UPARA )
         deallocate( VPERP )
         deallocate( VPARA )

         deallocate( UPERP_work )
         deallocate( UPARA_work )

         deallocate( f )

         deallocate( dfdupere )
         deallocate( dfdupare )
         deallocate( fperpe )


         deallocate( dfduper1 )
         deallocate( dfdupar1 )
         deallocate( fperp1 )

         deallocate( dfduper2 )
         deallocate( dfdupar2 )
         deallocate( fperp2 )

         deallocate( dfduper3 )
         deallocate( dfdupar3 )
         deallocate( fperp3 )

         deallocate( dfduper4 )
         deallocate( dfdupar4 )
         deallocate( fperp4 )

         deallocate( dfduper5 )
         deallocate( dfdupar5 )
         deallocate( fperp5 )

         deallocate( dfduper6 )
         deallocate( dfdupar6 )
         deallocate( fperp6 )

         deallocate( dfe_cql_uprp )
         deallocate( dfe_cql_uprl )

         deallocate( df1_cql_uprp )
         deallocate( df1_cql_uprl )

         deallocate( df2_cql_uprp )
         deallocate( df2_cql_uprl )

         deallocate( df3_cql_uprp )
         deallocate( df3_cql_uprl )

         deallocate( df4_cql_uprp )
         deallocate( df4_cql_uprl )

         deallocate( df5_cql_uprp )
         deallocate( df5_cql_uprl )

         deallocate( df6_cql_uprp )
         deallocate( df6_cql_uprl )

         deallocate( bqlavg_e )
         deallocate( cqlavg_e )
         deallocate( eqlavg_e )
         deallocate( fqlavg_e )

         deallocate( bqlavg_work )
         deallocate( cqlavg_work )
         deallocate( eqlavg_work )
         deallocate( fqlavg_work )

         deallocate( bqlavg_i1 )
         deallocate( cqlavg_i1 )
         deallocate( eqlavg_i1 )
         deallocate( fqlavg_i1 )

         if(eta2 .ne. 0.0) then
            deallocate( bqlavg_i2 )
            deallocate( cqlavg_i2 )
            deallocate( eqlavg_i2 )
            deallocate( fqlavg_i2 )
         end if

         if(eta3 .ne. 0.0) then
            deallocate( bqlavg_i3 )
            deallocate( cqlavg_i3 )
            deallocate( eqlavg_i3 )
            deallocate( fqlavg_i3 )
         end if

         if(eta4 .ne. 0.0) then
            deallocate( bqlavg_i4 )
            deallocate( cqlavg_i4 )
            deallocate( eqlavg_i4 )
            deallocate( fqlavg_i4 )
         end if

         if(eta5 .ne. 0.0) then
            deallocate( bqlavg_i5 )
            deallocate( cqlavg_i5 )
            deallocate( eqlavg_i5 )
            deallocate( fqlavg_i5 )
         end if

         if(eta6 .ne. 0.0) then
            deallocate( bqlavg_i6 )
            deallocate( cqlavg_i6 )
            deallocate( eqlavg_i6 )
            deallocate( fqlavg_i6 )
         end if

      end if

      if (read_solution .eqv. .false.) then

c         deallocate( new_to_org )

c         deallocate(niabegin_all )
c         deallocate( isize_all )

c        deallocate(Btmp)
c         deallocate(brhs)

c         deallocate( descBtmp_all )
c         deallocate( descbrhs_all )

         deallocate( p_brhs )
         deallocate( p_ipiv )

c         deallocate( itable )
c         deallocate( jtable )
c         deallocate( mtable )
c         deallocate( ntable )

      end if


 9000 continue
*     -----------------------
*     End of loop over nphi's
*     -----------------------

 9001 continue

      deallocate (xx, yy, xx_inv, yy_inv )

      deallocate ( uxxk, uxyk,
     &  uxzk, uyxk, uyyk, uyzk, uzxk, uzyk, uzzk )

      deallocate (dxuxx, dxuxy, dxuxz,
     &   dxuyx, dxuyy, dxuyz, dxuzx, dxuzy, dxuzz,
     &   dyuxx, dyuxy, dyuxz, dyuyx, dyuyy, dyuyz,
     &   dyuzx, dyuzy, dyuzz, dyyuxx, dyyuxy, dyyuxz,
     &   dyyuyx, dyyuyy, dyyuyz, dyyuzx, dyyuzy, dyyuzz,
     &   dxyuxx, dxyuxy, dxyuxz, dxyuyx, dxyuyy, dxyuyz,
     &   dxyuzx, dxyuzy, dxyuzz, dxxuxx, dxxuxy, dxxuxz,
     &   dxxuyx, dxxuyy, dxxuyz, dxxuzx, dxxuzy, dxxuzz,
     &   spx, spy, spz, uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz,
     &   bmod2, rhom1, bxn, byn, bzn, pressi, muhat, nu_star,
     &   ipsi, xna_slo, gradprlb, omgp2_slo, omgci_slo, zeff,
     &   rho_pol2d, psi_pol2d, drhodr, drhodz,
     &   dthedr, dthedz, reomg1a, reomg2a, reomg3a, xjx, xjy, xjz)

      deallocate (kappa_hatx, kappa_haty, kappa_hatz, rc,
     &   bhatx, bhaty, bhatz)

      deallocate ( acold, bcold, ccold,
     &   xkperp_cold, xkperp_cold2, ex, ey, ez, bxwave, bywave, bzwave,
     &   eplus_flux_plot, eminus_flux_plot, xkperp_flux_plot,
     &   ealpha, ebeta, eb, eplus, eminus, ealphakx, ebetakx,
     &   ealphaky, ebetaky, ebkx, ebky, rho_ant, rho_pla, ntilda_e,
     &   xb, xc, xd )


      deallocate ( wdote, wdoti1,
     &     wdoti2, wdoti3, wdot,
     &     wdoti4, wdoti5, wdoti6, fype, fypi1,
     &     fypi2, fypi3, fyp, fypi4, fypi5, fypi6, fpol0e, fpol0i1,
     &     fpol0i2, fpol0i3, fpol0i4, fpol0i5, fpol0i6, fx0e, fx0i1,
     &     fx0i2, fx0i3,
     &     fx0i4, fx0i5, fx0i6, fy0e , fy0i1, fy0i2, fy0i3,
     &     fy0i4, fy0i5, fy0i6, fz0e , fz0i1,
     &     fz0i2, fz0i3,
     &     fz0i4, fz0i5, fz0i6, bdotf, capr_fzeta, capr2,
     &     uzeta, utheta, gpsi, omgexb, jhat, kpsi, epsi, fprl0e,
     &     fprl0i1,
     &     fprl0i2, fprl0i3, fprl0i4, fprl0i5, fprl0i6, fx0, fy0, fz0,
     &     fpsi0, fpsi1,
     &     ftheta0, ftheta1, fpol1e , fpol1i1,
     &     fpol1i2, fpol1i3, fpol1i4, fpol1i5, fpol1i6, dfpol1e,
     &     dfpol1i1,
     &     dfpol1i2, dfpol1i3, dfpol1i4, dfpol1i5, dfpol1i6, redotj,
     &     pcre, pcim, fx_ant, fy_ant, fz_ant, fx_pla, fy_pla, fz_pla,
     &     rekxdotj, rekydotj,
     &     rekzdotj, xjprl, redotj1, redotj2, redotj3, redotj4, redotj5,
     &     redotj6, redotje, redotjt, redotji, redotjs, redotj_ehst,
     &     divq)

      deallocate ( psi, rho, theta,
     &     rhohatx, rhohaty, theta0, psi_dim, bx, by, bz,
     &     btau, bzeta, dxdth, dzdth, xntau, xkte, xkti,
     &     xkti2, xkti3, xkti4, xkti5, xkti6, xn1a, xnea, xn2a,
     &     xn3a, bmod, omgce, xn4a, xn5a, xn6a,
     &     omgci1, omgci2, omgci3, omgci4, omgci5, omgci6,
     &     omgpe2, bpol, capr_bpol, omgp12, omgp22, omgp32,
     &     omgp42, omgp52, omgp62, xiota, qsafety, bmod_mid,
     &     capr_bpol_mid2, rho_tor2d, psi_tor2d, omglh, reomglha )


      deallocate ( xjpxe, xjpye, xjpze,
     &   xjpx1, xjpy1, xjpz1,xjpx2, xjpy2, xjpz2, xjpx3,
     &   xjpy3, xjpz3, xjpx4, xjpy4, xjpz4, xjpx5, xjpy5,
     &   xjpz5, xjpx6, xjpy6, xjpz6 )

      deallocate ( xjpx_ehst, xjpy_ehst,
     &   xjpz_ehst, xj_slox, xj_sloy, xj_sloz,
     &   xjpx, xjpy, xjpz, xjpx_lab, xjpy_lab, xjpz_lab, pc )

      deallocate (xjpxe_lab, xjpye_lab, xjpze_lab)

      deallocate ( fdksav2d )
      deallocate ( feksav2d )
      deallocate ( ffksav2d )
      deallocate ( fgksav2d )
      deallocate ( faksav2d )
      deallocate ( fpksav2d )
      deallocate ( frksav2d )
      deallocate ( fqksav2d )
      deallocate ( fsksav2d )

      deallocate ( ealphakmod )
      deallocate ( ebetakmod )
      deallocate ( ebkmod )

      deallocate ( ealphak )
      deallocate ( ebetak )
      deallocate ( ebk )

      deallocate ( fdksav )
      deallocate ( feksav )
      deallocate ( ffksav )
      deallocate ( fgksav )
      deallocate ( faksav )
      deallocate ( fpksav )
      deallocate ( frksav )
      deallocate ( fqksav )
      deallocate ( fsksav )

      deallocate (capd)
      deallocate (capd2d)

      deallocate (xkb)
      deallocate (xkb2d)

      deallocate (xkperp2_slow, xkperp2_fast)
      deallocate (xkprl_a, P_a)
      deallocate (work)
      deallocate (workc)

      if (myid .eq. 0) then
         close (963)
         close (34)
         close (166)
      endif

      close (63)

      if(myid .eq. 0 .and. nphi_number .gt. 1)then
         t1 = second1(dummy)

         call aorsa2dSum(myid)

         tmin = (second1(dummy) - t1) / 60.
         write(6 , 2847) tmin
         write(15, 2847) tmin

      end if

      if (myid .eq. 0) then
         close (15)
         close (65)
      endif

      close(6)

*     -------------------------
*     stop parallel environment
*     -------------------------

      call blacs_barrier(icontxt, 'All')

      call blacs_gridexit( icontxt )
      call blacs_exit(0)


      end

!
!*************************************************************************
!

        subroutine pzgecopy( m,n, A,ia,ja,descA, B,ib,jb,descB)
        integer :: m,n,  ia,ja,descA(*), ib,jb,descB(*)
        complex :: A(*), B(*)

        complex :: alpha,beta

        alpha = 1.0
        beta = 0.0
*       ----------------------------------------------------
*       pzgeadd is new capability found in PBLAS V2 or later
*       otherwise, it can be implmented less efficiently by
*       repeated calls to pzcopy
*       ----------------------------------------------------
        call pzgeadd( 'N', m,n, alpha, A,ia,ja,descA,
     &        beta, B,ib,jb,descB )
        return
        end

c
c***************************************************************************
c


      subroutine fluxavg(f, favg, rho, nxdim, nydim, nrhodim,
     &   nnodex, nnodey, nnoderho, drho, dx, dy, capr, r0, vol, fvol)

      implicit none

      integer :: nxdim, nydim, nrhodim, nnodex, nnodey, nnoderho
      integer :: n, i, j

      real :: f(nxdim, nydim), favg(nrhodim), rho(nxdim, nydim), r0,
     &   drho, dx, dy, fvol(nrhodim), vol(nrhodim), capr(nxdim), pi,
     &   twopi

      pi = 3.141592654
      twopi = 2.0 * pi

      fvol = 0.0
      vol = 0.0
      favg = 0.0

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1
            if(n .le. nnoderho)then
               fvol(n) = fvol(n) + dx * dy * twopi * capr(i) * f(i,j)
               vol(n) =  vol(n)  + dx * dy * twopi * capr(i)
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


  100 format (1i10, 1p,8e12.4)
  102 format (2i10)
      return
      end


c
c***************************************************************************
c

c      subroutine rhograte(x, f, nx1, nx2, ans, nxdim, dvol)

c      implicit none

c      integer :: n, nx1, nx2, nxdim

c      real :: f(nxdim), x(nxdim), ans(nxdim), dvol(nxdim)

c      ans(nx1) = 0.0

c      do n = nx1, nx2 - 1

c        ans(n+1) = ans(n) + f(n) * dvol(n)

c      end do


c  100 format (1i10, 1p,8e12.4)
c  102 format (2i10)
c      return
c      end
c
c***************************************************************************
c


      subroutine rhograte(x, f, nx1, nx2, ans, nxdim, dvol)

      implicit none

      integer :: n, nx1, nx2, nxdim

      real :: f(nxdim), x(nxdim), ans(nxdim), dvol(nxdim)

      ans(nx1) = 0.0

      do n = nx1, nx2 - 1

         ans(n+1) = ans(n) + (f(n) * dvol(n) + f(n+1) * dvol(n+1)) / 2.

      end do

  100 format (1i10, 1p,8e12.4)
  102 format (2i10)
      return
      end

c
c***************************************************************************
c
