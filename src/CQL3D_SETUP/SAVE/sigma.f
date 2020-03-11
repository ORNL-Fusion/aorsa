

c
c***************************************************************************
c
      subroutine sigmad_cql3d(i, j, n, m, rho, rho_a,
     .   gradprlb, bmod, bmod0,
     .   xm, q, xn, xnuomg,
     .   xkt, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, nphi, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz,
     .   delta0, ndist, nupar, nuper, n_psi,
     .   n_psi_dim,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, nbessj, ncount)

*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi, ndist, myid, iflag, nproc, ni, mi

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, sgn_kprl, reson, dzetal
      real dakbdkb, xnuomg, gradprlb, bmod, bmod0, nu_coll, gamma_coll
      real akprl, gammab, rho, alpha, eps0, omgrf, v0i, emax
      real a, b, xnurf, pi, delta0, rhol
      real bx, by, bz, bratio, denom
      real time, t1, dummy, tmin, second1

      real xkxsav, xkysav, capr, argd
      real xkphi
      real xkalp, xkbet, xk0, rgamma

      complex zi, zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0_a, sig1_a, sig2_a, sig3_a, sig4_a, sig5_a
      complex sig0_h, sig1_h, sig2_h, sig3_h, sig4_h, sig5_h

      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz

      real uxx, uxy, uxz,
     1     uyx, uyy, uyz,
     2     uzx, uzy, uzz

      parameter (lmaxdim = 99)

      complex xil(0: lmaxdim), xilp(0: lmaxdim)
      complex exil(0: lmaxdim), exilp(0: lmaxdim),
     .                    exilovergam(0: lmaxdim)

      complex zetal, zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      integer  :: n_psi_dim
      integer :: nuper, nupar, n_psi

c      integer, parameter :: NBESSJ = 7
      integer NBESSJ, ncount, nc

      real :: UPERP(NUPER),UPARA(NUPAR)
      real :: UPERP0, UPARA0, dfdupar0, dfduper0

      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: W,K1(3),K2(3),KPER1,KPER2,ENORM,ZSPEC,ASPEC,BMAG,DensSPEC
      integer :: NSBESSJ,IFAIL
      COMPLEX WSPEC(3,3)
      complex :: factor
	real :: UminPara,UmaxPara
      real :: XI1(NUPER),JNXI1(NUPER, NBESSJ)
      real :: XI2(NUPER),JNXI2(NUPER, NBESSJ)


      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: rho_a(n_psi_dim), vc_mks


      integer :: i_uperp, i_upara, i_psi

      common/sigcom/zi, eps0, v0i, omgrf, xk0

      nu_coll =  .01 * omgrf

      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rhol = alpha / omgc
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * xnuomg)



      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi

      xkperp = sqrt(xkalp**2 + xkbet**2)

      akprl = abs(xkprl)



      sgn_kprl = sign(1.0, xkprl)

      if (akprl .lt. 0.01) then

         xkprl = 0.01 * sgn_kprl
         akprl = abs(xkprl)

      end if

*     -----------------------
*     Maxwellian distribution
*     -----------------------
      if(ndist .eq. 0 .or. ndist .eq. 2)then

         gamma = 0.5 * xkperp**2 * rhol**2
         rgamma = real(gamma)


         if(rgamma .ge. 1.0e-08)
     .      call besiexp(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

         if(rgamma .lt. 1.0e-08)
     .      call bes_expand(gamma, lmax, exil, exilp, lmaxdim,
     .                                                      exilovergam)


         sig0 = 0.0
         sig1 = 0.0
         sig2 = 0.0
         sig3 = 0.0
         sig4 = 0.0
         sig5 = 0.0


         do l = lmin, lmax
            labs = abs(l)

            reson = (omgrf - l * real(omgc)) / omgrf
            if (abs(reson) .lt. 0.02)then
               zetal = (omgrfc - l * omgc) / (xkprl * alpha)
               dzetal = omgrf * xnuomg / (xkprl * alpha)
            else
               zetal = (omgrf  - l * omgc) / (xkprl * alpha)
               dzetal = 0.0
            end if


            gammab = l * omgc / (2.0 * alpha * xkprl**2)
     .                                            * gradprlb / bmod
            gamma_coll = nu_coll / (akprl * alpha)


            if(xm .eq. xme)gammab = 0.0
            if(abs(gammab) .gt. 1000.0) gammab = 1000.0
            if(abs(gammab) .lt. .01)gammab = .01



            if (nzfun .eq. 0) call z_approx(sgn_kprl, zetal, 0.0,
     .                                                     z0, z1, z2)
            if (nzfun .eq. 1) call z_approx(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
            if (nzfun .eq. 2) call z_smithe(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
            if (nzfun .eq. 3) call z_table (sgn_kprl, zetal, gammab,
     .                                         gamma_coll, z0, z1, z2)


            al = 1.0 / (xkprl * alpha) * z0
            bl = 1.0 / (xkprl * alpha) * z1
            cl = 1.0 / (xkprl * alpha) * z2


            sig0l = - zieps0 * omgp2 * rhol**2
     .                                 * (exil(labs) - exilp(labs)) * al
            sig1l = - zieps0 * omgp2 * l**2 * exilovergam(labs) * al
            sig2l = - eps0 * omgp2 * l * (exil(labs) - exilp(labs)) * al
            sig3l = - zieps0 * omgp2 * 2.0 * exil(labs) * cl
            sig4l = - zieps0 * omgp2 * rhol * l * exilovergam(labs) * bl
            sig5l = - eps0 * omgp2 * rhol
     .                                 * (exil(labs) - exilp(labs)) * bl


            sig0 = sig0 + sig0l
            sig1 = sig1 + sig1l
            sig2 = sig2 + sig2l
            sig3 = sig3 + sig3l
            sig4 = sig4 + sig4l
            sig5 = sig5 + sig5l

         end do

         sig1 = sig1 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
         sig3 = sig3 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2


         sig0_a = zi * aimag(sig0)
         sig1_a = zi * aimag(sig1)
         sig2_a =       real(sig2)
         sig3_a = zi * aimag(sig3)
         sig4_a = zi * aimag(sig4)
         sig5_a =       real(sig5)



      end if

*     -----------------------------------
*     Funky METS calls for non Maxwellian:
*     -----------------------------------
      if (ndist .eq. 1 .or. ndist .eq. 2) then


         Emax = 0.5 * xm * vc_mks**2
         Enorm = Emax / 1.6e-19

         W = omgrf
         K1(1) = xkperp
         K1(2) = 0.0
         K1(3) = xkprl
         K2(1) = K1(1)
         K2(2) = K1(2)
         K2(3) = K1(3)
         KPER1 = xkperp
         KPER2 = xkperp


         ZSPEC = q / 1.6e-19
         ASPEC = xm / 1.67e-27
         BMAG = omgc * (xm / q)
         NSBESSJ = 2 * NBESSJ + 8
c         NSBESSJ = NBESSJ
         DensSPEC = xn
         bratio = bmod0 / bmod



!         call dummy_dist(NUPAR, NUPER,
!     &                 UminPara, UmaxPara,
!     &                 UPERP, UPARA, DFDUPER, DFDUPAR)



!        ------------------------------------------------
!        get CQL3D distribution function on the midplane
!        ------------------------------------------------
         call cql3d_dist(nupar, nuper, n_psi,
     .                 n_psi_dim, rho_a, rho,
     .                 UminPara,UmaxPara,
     .                 df_cql_uprp, df_cql_uprl,
     .                 UPERP, UPARA, DFDUPER, DFDUPAR)

!        go to 1000

!        ------------------------------------------------
!        map CQL3D distribution function off the midplane
!        ------------------------------------------------
         if(bratio .gt. 0.0)then
            do ni = 1, nuper
               do mi = 1, nupar

                  dfdupar0 = dfdupar(ni, mi)
                  dfduper0 = dfduper(ni, mi)

                  argd = uperp(ni)**2 * (1. - bratio) + upara(mi)**2
                  if (argd .le. 0.0) argd = 1.0e-06

                  upara0  = sign(1.0, upara(mi)) * sqrt(argd)

                  dfduper(ni, mi) = dfduper0 * sqrt(bratio)
     .               + dfdupar0 * uperp(ni) / upara0 * (1.0 - bratio)

                  dfdupar(ni, mi) = dfdupar0 * upara(mi) / upara0


               end do
            end do
         end if

 1000    continue



         t1 = second1(dummy)


         do nc = 1, ncount

            call WMATPRECALC_AORSA(ZSPEC,ASPEC,ENORM,BMAG,KPER1,UPERP,
     &                   NUPER,NBESSJ,NSBESSJ,XI1,JNXI1,IFAIL)

         end do

         time = second1(dummy) - t1

  	   tmin = time / 60.
         if (myid.eq.0) then
	      write(6 ,833) time
            write(15,833) time
         endif
	
  833 format('time taken by WMATPRECALC =',f9.3,4h sec)

         t1 = second1(dummy)

         do nc = 1, ncount

            call GETNONMAXSIGMA_AORSA_NEW(W,
     &                          ZSPEC,ASPEC,DensSPEC,BMAG,
     &                          K1,XI1,JNXI1,
     &                          K1,XI1,JNXI1,NBESSJ,
     &                          Enorm,UminPara,UmaxPara,
     &                          NUPAR,NUPER,UPERP,UPARA,
     &			        DFDUPER,DFDUPAR,
     &                          WSPEC,IFAIL)


         end do

         time = second1(dummy) - t1

  	   tmin = time / 60.
         if (myid.eq.0) then
	      write(6 ,834) time
            write(15,834) time
         endif
	
  834 format('time taken by GETNONMAXSIGMA =',f9.3,4h sec)


         factor = cmplx(0.,-omgrf * eps0)

         sig1 = WSPEC(1,1) * factor
         sig2 = WSPEC(1,2) * factor
         sig3 = WSPEC(3,3) * factor
         sig4 = 0.0
         sig5 = 0.0
         sig0 = 0.0
         if (xkperp .gt. 0.01) then
            sig4 = WSPEC(3,1) / xkperp * factor
            sig5 = WSPEC(3,2) / xkperp * factor
            sig0 = (WSPEC(2,2) * factor - sig1) / xkperp**2
         endif



         sig0_h =       real(sig0)
         sig1_h =       real(sig1)
         sig2_h = zi * aimag(sig2)
         sig3_h =       real(sig3)
         sig4_h =       real(sig4)
         sig5_h = zi * aimag(sig5)




c         write(6, *) "rho = ", rho
c         write(6, *) "xn = ", xn
c         write(6, *) "xkt = ", xkt
c         write(6, *) "xkperp = ", xkperp
c         write(6, *) "xkprl = ", xkprl
c         write(6, *) "bmod = ", bmod
c         write(6, *) "bmod0 = ", bmod0
c         write(6, *) "bratio =", bratio
c         write(6, *) "xm = ", xm
c         write(6, *) "q = ", q
c         write(6, *) "lmax = ", lmax
c         write(6, *) "omgrf = ", omgrf
c         write(6, *) "omgc = ", omgc
         reson = omgrf / omgc
c         write(6, *) "reson = ", reson
c         write(6, *) "vc_mks = ", vc_mks
c         write(6, *) "DFDUPER(10, 30) = ", DFDUPER(10, 30)
c         write(6, *) "DFDUPAR(10, 30) = ", DFDUPAR(10, 30)
c         write(6, *) "UPERP(10) = ", UPERP(10)
c         write(6, *) "UPARA(30) = ", UPARA(30)
c         write(6, *) "Enorm = ", Enorm
c         write(6, *) "UminPara = ", UminPara
c         write(6, *) "UmaxPara = ", UmaxPara
         write(6, *) "NUPER = ", NUPER
         write(6, *) "NUPAR = ", NUPAR

         write(15, *) "NUPER = ", NUPER
         write(15, *) "NUPAR = ", NUPAR




      end if

!      if(ndist .eq. 2) then

!         sig0 = sig0_a + sig0_h
!         sig1 = sig1_a + sig1_h
!         sig2 = sig2_a + sig2_h
!         sig3 = sig3_a + sig3_h
!         sig4 = sig4_a + sig4_h
!         sig5 = sig5_a + sig5_h

!      end if



         write(6, 1312)sig0
         write(6, 1312)sig1
         write(6, 1312)sig2
         write(6, 1312)sig3
         write(6, 1312)sig4
         write(6, 1312)sig5

         write(15, 1312)sig0
         write(15, 1312)sig1
         write(15, 1312)sig2
         write(15, 1312)sig3
         write(15, 1312)sig4
         write(15, 1312)sig5



*     -----------------------------
*     Swanson's rotation (original):
*     -----------------------------
      sigxx = sig1 + sig0 * xkbet**2
      sigxy = sig2 - sig0 * xkbet * xkalp
      sigxz = sig4 * xkalp + sig5 * xkbet

      sigyx = - sig2 - sig0 * xkbet * xkalp
      sigyy =   sig1 + sig0 * xkalp**2
      sigyz =   sig4 * xkbet - sig5 * xkalp

      sigzx = sig4 * xkalp - sig5 * xkbet
      sigzy = sig4 * xkbet + sig5 * xkalp
      sigzz = sig3


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
      end

c
c***************************************************************************
c

      subroutine cql3d_dist(nupar, nuper, n_psi,
     .                       n_psi_dim, rho_a, rho,
     .                       UminPara,UmaxPara,
     .                       df_cql_uprp, df_cql_uprl,
     .                       UPERP, UPARA, DFDUPER, DFDUPAR)

      implicit none
      integer, intent(in) :: nupar, nuper, n_psi
      integer, intent(in) :: n_psi_dim

      real, intent(out) :: UminPara,UmaxPara
      real, intent(out) :: UPERP(NUPER)
      real, intent(out) :: UPARA(NUPAR)
      real, intent(out) :: DFDUPER(NUPER,NUPAR),
     .                      DFDUPAR(NUPER,NUPAR)

      real, intent(in) :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real, intent(in) :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: f, u2, rho, rho_a(n_psi_dim)

      integer n, m, i_psi, i_psin

      i_psi = 1


      do i_psin = 1, n_psi - 1
         if(rho .ge. rho_a(i_psin)     .and.
     .      rho .lt. rho_a(i_psin + 1))
     .      i_psi = i_psin
      end do


      if(rho .ge. rho_a(n_psi)) i_psi = n_psi

c      write(6,  *)
c      write(6,  *) "i_psi = ", i_psi

c      write(15, *)
c      write(15, *) "i_psi = ", i_psi


      if(i_psi .lt. n_psi)then
         do n = 1, nuper
            do m = 1, nupar
               DFDUPER(n, m) = df_cql_uprp(n, m, i_psi)
     .            + (df_cql_uprp(n, m, i_psi + 1)
     .                                     - df_cql_uprp(n, m, i_psi))
     .                    * (rho              - rho_a(i_psi))
     .                    / (rho_a(i_psi + 1) - rho_a(i_psi))
               DFDUPAR(n, m) = df_cql_uprl(n, m, i_psi)
     .            + (df_cql_uprl(n, m, i_psi + 1)
     .                                     - df_cql_uprl(n, m, i_psi))
     .                    * (rho              - rho_a(i_psi))
     .                    / (rho_a(i_psi + 1) - rho_a(i_psi))
            end do
         end do
      end if



      if(i_psi .eq. n_psi)then
         do n = 1, nuper
            do m = 1, nupar
               DFDUPER(n, m) = df_cql_uprp(n, m, i_psi)
               DFDUPAR(n, m) = df_cql_uprl(n, m, i_psi)
            end do
         end do
      end if



  310 format(1p6e12.4)
  311 format(i10, 1p6e12.4)


      return
      end subroutine cql3d_dist

c
c***************************************************************************
c
       subroutine dummy_dist(nupar, nuper,
     &                       UminPara, UmaxPara,
     &                       UPERP, UPARA, DFDUPER, DFDUPAR)
       implicit none

       integer, intent(in) :: nupar, nuper
       real, intent(out) :: UminPara,UmaxPara
       real, intent(out) :: UPERP(NUPER)
       real, intent(out) :: UPARA(NUPAR)
       real, intent(out) :: DFDUPER(NUPER, NUPAR),
     &                      DFDUPAR(NUPER, NUPAR)
       real :: pi,fnorm,f,u2
       integer n,m

       pi = 4.0*atan(1.0)
       fnorm = pi**1.5

       do n = 1, NUPER
          UPERP(n) = 3.5*(real(n-1)/real(NUPER-1))
       end do

       UminPara = -3.5
       UmaxPara =  3.5

       do m = 1, NUPAR
          UPARA(m) = 3.5*(-1.0+2.*(real(m-1)/real(NUPAR-1)))
       end do

       do n = 1, NUPER
          do m = 1, NUPAR
             u2 = UPERP(n)**2 + UPARA(m)**2
             f = exp(-u2) / fnorm
             DFDUPER(n, m) = f *(-2.0 * UPERP(n))
             DFDUPAR(n, m) = f *(-2.0 * UPARA(m))
          end do
       end do

       return
       end subroutine dummy_dist

c
c***************************************************************************
c
      subroutine sigmah_stix(i, j, n, m,
     .   gradprlb, bmod,
     .   xm, q, xn, xnuomg,
     .   xkt, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, nphi, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz,
     .   delta0, iflag_gammab)

*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi, iflag_gammab

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme, delta0
      real xkprl_eff, fgam, y0, sgn_kprl, reson, dzetal
      real dakbdkb, xnuomg, gradprlb, bmod, nu_coll, gamma_coll
      real akprl, gammab, rho, alpha, eps0, omgrf, v0i
      real a, b
      real bx, by, bz

      real xkxsav, xkysav, capr
      real xkphi
      real xkalp, xkbet, xk0, rgamma

      complex zi, zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz

      real uxx, uxy, uxz,
     1     uyx, uyy, uyz,
     2     uzx, uzy, uzz

      parameter (lmaxdim = 99)

      complex xil(0: lmaxdim), xilp(0: lmaxdim)
      complex exil(0: lmaxdim), exilp(0: lmaxdim),
     .                    exilovergam(0: lmaxdim)

      complex zetal, zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff


      common/sigcom/zi, eps0, v0i, omgrf, xk0

      nu_coll =  .01 * omgrf


      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rho = alpha / omgc
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * xnuomg)



      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi

      xkperp = sqrt(xkalp**2 + xkbet**2)

      akprl = abs(xkprl)



      sgn_kprl = sign(1.0, xkprl)

      if (akprl .lt. 0.01) then

         xkprl = 0.01 * sgn_kprl
         akprl = abs(xkprl)

      end if


      gamma = 0.5 * xkperp**2 * rho**2
      rgamma = real(gamma)


      if(rgamma .ge. 1.0e-08)
     .   call besiexp(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

      if(rgamma .lt. 1.0e-08)
     .   call bes_expand(gamma, lmax, exil, exilp, lmaxdim, exilovergam)


      sig0 = 0.0
      sig1 = 0.0
      sig2 = 0.0
      sig3 = 0.0
      sig4 = 0.0
      sig5 = 0.0


      do l = lmin, lmax
         labs = abs(l)


         reson = (omgrf - l * real(omgc)) / omgrf
         if (abs(reson) .lt. 0.02)then
            zetal = (omgrfc - l * omgc) / (xkprl * alpha)
            dzetal = omgrf * xnuomg / (xkprl * alpha)
         else
            zetal = (omgrf  - l * omgc) / (xkprl * alpha)
            dzetal = 0.0
         end if


         gammab = l * omgc / (2.0 * alpha * xkprl**2)
     .                                            * gradprlb / bmod

         gamma_coll = nu_coll / (akprl * alpha)



         if(xm .eq. xme)gammab = 0.0

         if(abs(gammab) .gt. 1000.0) gammab = 1000.0
         if(abs(gammab) .lt. .01)gammab = .01



         if (nzfun .eq. 0) call z_approx(sgn_kprl, zetal, 0.0,
     .                                                     z0, z1, z2)
         if (nzfun .eq. 1) call z_approx(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
         if (nzfun .eq. 2) call z_smithe(sgn_kprl, zetal, gammab,
     .                                                     z0, z1, z2)
         if (nzfun .eq. 3) call z_table (sgn_kprl, zetal, gammab,
     .                                         gamma_coll, z0, z1, z2)





         al = 1.0 / (xkprl * alpha) * z0
         bl = 1.0 / (xkprl * alpha) * z1
         cl = 1.0 / (xkprl * alpha) * z2





         sig0l = - zieps0 * omgp2 * rho**2 * (exil(labs) - exilp(labs))
     .           * al
         sig1l = - zieps0 * omgp2 * l**2 * exilovergam(labs) * al
         sig2l = - eps0 * omgp2 * l * (exil(labs) - exilp(labs)) * al
         sig3l = - zieps0 * omgp2 * 2.0 * exil(labs) * cl
         sig4l = - zieps0 * omgp2 * rho * l * exilovergam(labs)  * bl
         sig5l = - eps0 * omgp2 * rho * (exil(labs) - exilp(labs)) * bl



         sig0 = sig0 + sig0l
         sig1 = sig1 + sig1l
         sig2 = sig2 + sig2l
         sig3 = sig3 + sig3l
         sig4 = sig4 + sig4l
         sig5 = sig5 + sig5l



      end do

      sig1 = sig1 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
      sig3 = sig3 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2

*     -----------------------------
*     Swanson's rotation (original):
*     -----------------------------
      sigxx = sig1 + sig0 * xkbet**2
      sigxy = sig2 - sig0 * xkbet * xkalp
      sigxz = sig4 * xkalp + sig5 * xkbet

      sigyx = - sig2 - sig0 * xkbet * xkalp
      sigyy =   sig1 + sig0 * xkalp**2
      sigyz =   sig4 * xkbet - sig5 * xkalp

      sigzx = sig4 * xkalp - sig5 * xkbet
      sigzy = sig4 * xkbet + sig5 * xkalp
      sigzz = sig3


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
      end

c
c***************************************************************************
c

      subroutine sigmac_stix(i, j, n, m,
     .   xm, q, xn, xnuomg,
     .   xkt, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, nphi, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz,
     .   delta0)

*     ----------------------------------------------------
*     This routine calculates sigma_cold in the Stix frame
*     ----------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xkb, akb
      real dakbdkb
      real akprl, gammab, alpha, eps0, omgrf, v0i
      real cosalp, sinalp, a, b
      real bx, by, bz
      real xkxsav, xkysav, capr
      real xkphi
      real xkalp, xkbet, xnuomg, xk0, delta0

      complex zi, omgrfc


      complex sig0, sig1, sig2, sig3, sig4, sig5

      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz

      real uxx, uxy, uxz,
     1     uyx, uyy, uyz,
     2     uzx, uzy, uzz


      parameter (lmaxdim = 99)

      complex tpl, tml, tplp, tplpp, zetalp, zetalm, zieps0,
     .   bl, gamma



      common/sigcom/zi, eps0, v0i, omgrf, xk0




      zieps0 = zi * eps0
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * xnuomg)


      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi



      sig0 = 0.0
      sig1 = zieps0 * omgrfc * omgp2 / (omgrfc**2 - omgc**2)
      sig2 = - eps0 * omgc   * omgp2 / (omgrfc**2 - omgc**2)
      sig3 = zieps0 * omgp2 / omgrfc
      sig4 = 0.0
      sig5 = 0.0


*     -------------------
*     Swanson's rotation:
*     -------------------
      sigxx = sig1 + sig0 * xkbet**2
      sigxy = sig2 - sig0 * xkbet * xkalp
      sigxz = sig4 * xkalp + sig5 * xkbet

      sigyx = - sig2 - sig0 * xkbet * xkalp
      sigyy =   sig1 + sig0 * xkalp**2
      sigyz =   sig4 * xkbet - sig5 * xkalp

      sigzx = sig4 * xkalp - sig5 * xkbet
      sigzy = sig4 * xkbet + sig5 * xkalp
      sigzz = sig3


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
      end

c
c***************************************************************************
c


      subroutine besiexp(gamma, lmax, expbes, expbesp, lmaxdim,
     .   expbesovergam)

*     ----------------------------------------------------------
*     Calculates exp(-gamma) times the modified bessel functions
*     (and derivatives) of order up to lmax
*     ----------------------------------------------------------

      implicit none

      integer lmax, nmax, ier, l, lmaxdim
      real gammod
      complex gamma, expbes(0: lmaxdim), expbesp(0: lmaxdim),
     .   expbesovergam(0: lmaxdim),
     .   xil(0: lmaxdim), xilp(0: lmaxdim), exgam

      complex b(100)

      exgam = exp(-gamma)
      gammod = cabs(gamma)

      if(gammod .le. 700.)then
         nmax = lmax + 1
         call besic(gamma, nmax, b, ier)
         if(ier .ne. 0)write(6,100) ier

         do l = 0, lmax
            xil(l) = b(l+1)
         end do

         do l = 0, lmax
           if(l .eq. 0) xilp(0) = xil(1)
           if(l .ne. 0) xilp(l) = xil(l-1) - l / gamma * xil(l)
           expbes(l) = exgam * xil(l)
           expbesp(l) = exgam * xilp(l)
         end do
      end if

      if(gammod .gt. 700.)then
         do l = 0, lmax
            call bes_asym(gamma, l, expbes(l), expbesp(l))
         end do
      end if

      do l = 0, lmax
         expbesovergam(l) = expbes(l) / gamma
      end do

  100 format('ier = ', i5, 'besic failed')
      return
      end


c
c***************************************************************************
c


      subroutine bes_expand(gamma, lmax, expbes, expbesp, lmaxdim,
     .   expbesovergam)

*-------------------------------------------------------------------
*     Calculates exp(-gamma) times the modified bessel functions
*     (and derivatives) of order up to lmax using second order
*     expansion for small argument
*-------------------------------------------------------------------

      implicit none

      integer lmax, nmax, ier, l, lmaxdim
      real factrl, factl
      complex gamma, expbes(0:lmaxdim), expbesp(0:lmaxdim),
     .   expbesovergam(0:lmaxdim), xilovergam,
     .   xil(0:lmaxdim), xilp(0:lmaxdim), exgam

      complex b(100)

      exgam = 1.0 - gamma + gamma**2 / 2.0

      do l = 0, lmax
         factl = factrl(l)
         xil(l) = gamma**l / (2**l * factl) *
     .                                 ( 1. + gamma**2 / (4. * (l+1)))
         xilp(l) = gamma**(l-1) / (2**l * factl) *
     .                     (l + (l+2) * gamma**2 / (4. * (l+1)))

         xilovergam = gamma**(l-1) / (2**l * factl) *
     .                                 ( 1. + gamma**2 / (4. * (l+1)))

         expbes(l) = exgam * xil(l)
         expbesp(l) = exgam * xilp(l)
         expbesovergam(l) = exgam * xilovergam

c         write(6, 100)l, expbes(l), expbesp(l)
      end do

  100 format(i10, 1p8e12.4)

      return
      end


c
c***************************************************************************
c

      subroutine bes_asym(z, n, exil, exilp)

      implicit none

      integer mu, n
      real pi
      complex z, exil, exilp
      data pi/3.141592654/

      mu = 4 * n**2
      exil =  1.0 / csqrt(2.0 * pi * z)
     1   * (1.0
     1   - (mu - 1)/(8.0 * z)
     1   + (mu - 1) * (mu - 9) / (2.0 * (8.0 * z)**2)
     1   - (mu - 1) * (mu - 9) * (mu - 25) / (6.0 * (8.0 * z)**3)  )
      exilp = 1.0 / csqrt(2.0 * pi * z)
     1   * (1.0
     1   - (mu + 3)/(8.0 * z)
     1   + (mu - 1) * (mu + 15) / (2.0 * (8.0 * z)**2)
     1   - (mu - 1) * (mu - 9) * (mu + 35) / (6.0 * (8.0 * z)**3)  )

      return
      end

c
c***************************************************************************
c


      function fzeta (arg)
      complex a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, arg, aux0
     1   , aux1, term, z, zz, fzeta
      data d1r/0.0/
      data d1i/1.77245385090551/
      data d2r/0.0/
      data d2i/3.54490770181103/
      data d3r/0.0/
      data d3i/7.08981540362206/
      data d4/0.33333333333333/
      data eps/1.0E-07/

c     data d4/0.33333333333333/
c     common/zetcom/eps
cray  code analysis
cray  optimize
c
      i = 0
      z = arg
      zz = z*z
      x = real(z)
      y = aimag(z)
      d1 = cmplx(d1r,d1i)
      d2 = cmplx(d2r,d2i)
      d3 = cmplx(d3r,d3i)
      ymag = abs(y)
      if (ymag - 1.0 .ge. 0.) then
c
c     continued fraction method: abs(y).ge.1.0
c
         y0 = y
         y = ymag
         aux1 = 1.5 - z*z
         aux2 = 0.0
         del = 1.5
         a1 = 0.0
         a2 = -1.0
         b1 = 1.0
         b2 = aux1
         c1 = a2/b2
c
  100    continue
         aux1 = aux1 + 2.0
         aux2 = aux2 - del
         del = del + 2.0
         a3 = aux1*a2 + aux2*a1
         b3 = aux1*b2 + aux2*b1
         c2 = a3/b3
         c3 = c2 - c1
         c3r = real(c3)
         c3i = aimag(c3)
         if (abs(c3r) + abs(c3i) .lt. eps) go to 110
         a1 = a2
         a2 = a3
         b1 = b2
         b2 = b3
         c1 = c2
         go to 100
  110    continue
         if (y0 .lt. 0.) then
            y = y0
            c2 = conjg(c2) - d3*z*exp(-zz)
         endif
         aux0 = -(0.5*c2 + 1.0)/z
      else
c
c     asymptotic series method: abs(x).ge.4.0 and abs(y).lt.1.0
c
         xmag = abs(x)
         if (xmag - 4.0 .lt. 0.) go to 130
         term = 1.0/z
         aux0 = -term
         aux1 = 0.5*term**2
         p = 1.0
         if (y .le. 0.) then
            if (y .ne. 0.) then
               aux0 = aux0 + d2*exp(-zz)
            else
               aux0 = aux0 + d1*exp(-zz)
            endif
         endif
  120    continue
         term = aux1*term*p
         aux0 = aux0 - term
         p = p + 2.0
         termr = real(term)
         termi = aimag(term)
c     if(abs(termr)+abs(termi).lt.eps)30,18
         if (abs(termr) + abs(termi) .lt. eps) go to 160
         go to 120
c
c     power series method: abs(x).lt.4.0 and abs(y).lt.1.0
c
  130    continue
         aux0 = 1.0
         aux1 = -(zz + zz)
         aux2 = eps/(eps + xmag + ymag)
         term = d4*aux1
         p = 3.0
  140    continue
         aux0 = aux0 + term
         termr = real(term)
         termi = aimag(term)
c     if(abs(termr)+abs(termi).lt.aux2)26,24
         if (abs(termr) + abs(termi) .lt. aux2) go to 150
         p = p + 2.0
         term = aux1*term/p
         go to 140
  150    continue
         aux0 = d1*exp(-zz) - 2.0*z*aux0
      endif
  160 continue
      fzeta = aux0
      if (i .le. 0) return
      fzeta = -2.0*(1.0 + arg*aux0)
      return
      end

c
c***************************************************************************
c
      real function second1_old(dummy)

      implicit none
      
      integer :: v(8)
c      integer mtime, mclock
      real dummy

c*****Fortran 90 standard for wall clock time
c      call date_and_time(values=v)
c      second1_old=(v(5)*3600)+(v(6)*60)+v(7)+0.001d0*v(8)

c*****FALCON:
      double precision mpi_wtime
      external mpi_wtime
      second1_old = MPI_WTIME()


c*****EAGLE:
c      mtime = mclock()
c      second1_old  = 0.01 * mtime

      return
      end
      
c
c***************************************************************************
c
      
      real function second1(dummy)
!     returns a monotonically increasing wall-clock time-stamp in seconds 
!     with an attempt to cross multiple days between calls 
!     (assumes one day if end-of-the-month is spanned between calls)

      implicit none
      
      real dummy

      integer :: v(8)

      integer,save :: ld
      logical,save :: elapsed=.false.

!*****Fortran 90 standard for wall clock time
      call date_and_time(values=v)

      if(elapsed)then  ! second has been called before
        if(v(3).ge.ld)then ! days have advanced monotonically
         second1=(v(3)-ld)*86400+(v(5)*3600)+(v(6)*60)+v(7)+0.001d0*v(8)
        else  ! day counter reset past end of the month
          print *,'WARNING: SECOND1 TRIPPED OVER THE END OF THE MONTH'
          print *,'         SECOND1 IS ASSUMING ONLY ONE ELAPSED DAY'
!         proper fix requires full blown calendar program
          ld=0
         second1=(v(3)-ld)*86400+(v(5)*3600)+(v(6)*60)+v(7)+0.001d0*v(8)
        endif
      else  ! first call
        second1=(v(5)*3600)+(v(6)*60)+v(7)+0.001d0*v(8)
        ld=v(3)
        elapsed=.true.
      endif

      return
      end


