c
c***************************************************************************
c

      subroutine dkw(i, j, n, m, rho, rho_a,
     .   gradprlb, bmod, bmod0,
     .   xm, q, xn, xnuomg,
     .   xkt, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, nphi, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   dkwxx, dkwxy, dkwxz, 
     .   dkwyx, dkwyy, dkwyz, 
     .   dkwzx, dkwzy, dkwzz,
     .   delta0, ndist, nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .   nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .   i_sav, j_sav, upshift, damping, xk_cutoff, rt, 
     .   nkx2, nky2)
        
*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------
c      use zfun_hilbert  

      implicit none
      
      complex dkwxx, dkwxy, dkwxz, 
     .        dkwyx, dkwyy, dkwyz, 
     .        dkwzx, dkwzy, dkwzz
     
      complex dkwxxl, dkwxyl, dkwxzl, 
     .        dkwyxl, dkwyyl, dkwyzl, 
     .        dkwzxl, dkwzyl, dkwzzl           
      
      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi, ndist, myid, iflag, nproc, ni, mi,
     .    i_sav, j_sav, ni0, mi0, upshift, nkx2, nky2
      integer n_upper, n_lower, m_upper, m_lower     
      
     
      real u0, u2, fnorm, f_cql, rt, sgn_xkalp, dkperp_dkalp
      real fnmax, fmmax, nperp
      real fna, fma, fstepn, fstepm, fstep

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, y, sgn_kprl, reson, duperp, dupara
      real xkprl_eff0
      real dzetal(lmin:lmax), descrim
      real dakbdkb, xnuomg, gradprlb, bmod, bmod0, nu_coll
      real akprl,  rho, alpha, eps0, omgrf, v0i, emax, akprl_min
      real gammab(lmin:lmax), gamma_coll(lmin:lmax)
      real a, b, xnurf, pi, delta0, rhol
      real bx, by, bz, bratio, denom
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt

      real xkxsav, xkysav, capr, argd
      real xkphi
      real xkalp, xkbet, xk0, rgamma, xk_cutoff, damping, kr, step

      complex zi, zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0_a, sig1_a, sig2_a, sig3_a, sig4_a, sig5_a
      complex sig0_h, sig1_h, sig2_h, sig3_h, sig4_h, sig5_h

      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
      logical :: l_interp   !new
      logical :: l_first  !new
      integer :: nkperp
      real :: kperp_max !new


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

      complex zetal(lmin:lmax)
      complex  zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      integer  :: n_psi_dim
      integer :: nuper, nupar, n_psi

      integer NBESSJ

      real :: UPERP(NUPER),UPARA(NUPAR)
      real :: UPERP0, UPARA0

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

      integer :: ieer ! fft error flag
      integer :: i_uperp, i_upara, i_psi
      
      common/upcom/akprl_min     
      
      allocate( dfduper0(nuper, nupar) )
      allocate( dfdupar0(nuper, nupar) )

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
      
      dkperp_dkalp = sign(1.0, xkalp)
      
         
      
*     ------------------------------------
*     Optional: leave out upshift in xkprl
*     --------------------------------- --          
      if (upshift .eq. 0)  xkprl = uzz * xkphi
c      if (upshift .eq. 0)   xkprl = nphi / rt
      
      if (upshift .eq. -1) then      
         if (xkperp  .gt. xk_cutoff) xkprl = uzz * xkphi
      end if
      
      if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
      if (xkperp .eq. 0.0) xkperp = 1.0e-08
                        
      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl)                      
                     
!     ----------------------------------------------
!     Optional: Don't allow xkprl to be 0 (upshift = -2)
!     ----------------------------------------------        
      if (upshift .eq. -2) then
         if (akprl .lt. akprl_min) then
            xkprl = akprl_min* sgn_kprl
         end if 
      end if                       
             
      
      
      if(xkperp .gt. kperp_max)then
         write (6, *)"xkperp is gt kperp_max in dkw"
         write (15, *)"xkperp is gt kperp_max in dkw"
      end if
            
      
*     ---------------------------------
*     Calculate zetal(l) and gammab(l)
*     ---------------------------------      
      
      do l = lmin, lmax
         labs = abs(l)

         reson = (omgrf - l * real(omgc)) / omgrf
c        if (abs(reson) .lt. 0.02)then
         if (rho .gt. 1.0) then
            zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
            dzetal(l) = omgrf * xnuomg / (xkprl * alpha)
         else
            zetal(l) = (omgrf  - l * omgc) / (xkprl * alpha)
            dzetal(l) = 0.0
         end if
         
c        zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
c         dzetal(l) = omgrf * xnuomg / (xkprl * alpha)


         gammab(l) = abs(l * omgc / (2.0 * alpha * xkprl**2)
     .                                           * gradprlb / bmod)
         gamma_coll(l) = nu_coll / (akprl * alpha)


         if(xm .eq. xme)gammab(l) = 0.0
c         if(abs(gammab(l)) .gt. 1000.0) gammab(l) = 1000.0
         if(abs(gammab(l)) .lt. .01)gammab(l) = .01


      enddo

      
      
*     ------------------------------------------------
*     Calculate Brambilla's xkrpl_eff using l = 1 only
*     ------------------------------------------------
      y0 = 1.5
      y = y0
      

      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gammab(1) * y) - 1.)
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            descrim = 1. - 4. * gammab(1) * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gammab(1) * y) )
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if
                    

*     -----------------------
*     Maxwellian distribution
*     -----------------------

c      if(ndist .eq. 0)then

         gamma = 0.5 * xkperp**2 * rhol**2
         rgamma = real(gamma)


         if(rgamma .ge. 1.0e-08)
     .      call besiexp(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

         if(rgamma .lt. 1.0e-08)
     .      call bes_expand(gamma, lmax, exil, exilp, lmaxdim,
     .                                                      exilovergam)


         dkwxx = 0.0
         dkwxy = 0.0
         dkwxz = 0.0
         dkwyx = 0.0
         dkwyy = 0.0
         dkwyz = 0.0
         dkwzx = 0.0
         dkwzy = 0.0
         dkwzz = 0.0

         do l = lmin, lmax
            labs = abs(l)

           if(nzfun .eq. 0) call z_approx(sgn_kprl, zetal(l), 0.0,
     .                                                     z0, z1, z2)
           if(nzfun .eq. 1) call z_approx(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 2) call z_smithe(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 3) call z_table(sgn_kprl,zetal(l),gammab(l),
     .                                      gamma_coll(l), z0, z1, z2)


            al = 1.0 / (xkprl * alpha) * z0
            bl = 1.0 / (xkprl * alpha) * z1
            cl = 1.0 / (xkprl * alpha) * z2
            
            
            dkwxxl = -zieps0 * omgp2 * l**2 / xkperp * (exilp(labs)
     .         - exil(labs) * (1. + 1. / gamma) ) * al
     
            dkwxyl =  eps0 * omgp2 * l / xkperp * (
     .         (2. * gamma + l**2 / gamma + 1.) * exil(labs)
     .         - ( 1. + 2. * gamma) * exilp(labs)   ) * al
     
            dkwxzl = -zieps0 * omgp2 * rhol * l * (exilp(labs)
     .         - exil(labs) - exil(labs) / gamma) * bl
     
     
     
            dkwyxl =  -eps0 * omgp2 * l / xkperp * (
     .        (2. * gamma + l**2 / gamma - 1.) * exil(labs)
     .        - ( 1. + 2. * gamma) * exilp(labs)   ) * al
     
            dkwyyl = zieps0 * omgp2 / xkperp * (
     .        exil(labs) * (3. * l**2 - 2. * gamma + l**2 / gamma
     .        + 4. * gamma**2)  - exilp(labs)* (l**2 + 4. * gamma**2) )
     .        * al
     
            dkwyzl =  eps0 * omgp2 * rhol * (
     .        (1. - l**2 / gamma - 2. * gamma) * exil(labs)
     .        + ( 1. + 2. * gamma) * exilp(labs)   ) * bl
     
     
     
            dkwzxl = -zieps0 * omgp2 * rhol * 
     .                l * (exilp(labs) - exil(labs)) * bl
            
            dkwzyl = - eps0 * omgp2 * rhol * (
     .        ( -l**2 / gamma - 2. * gamma) * exil(labs)
     .        + 2. * gamma * exilp(labs)   ) * bl
     
            dkwzzl = -zieps0 * omgp2 * xkperp * rhol**2
     .                  * (exilp(labs) - exil(labs)) * cl

            dkwxx = dkwxx + dkwxxl * dkperp_dkalp
            dkwxy = dkwxy + dkwxyl * dkperp_dkalp
            dkwxz = dkwxz + dkwxzl * dkperp_dkalp
            dkwyx = dkwyx + dkwyxl * dkperp_dkalp
            dkwyy = dkwyy + dkwyyl * dkperp_dkalp
            dkwyz = dkwyz + dkwyzl * dkperp_dkalp
            dkwzx = dkwzx + dkwzxl * dkperp_dkalp
            dkwzy = dkwzy + dkwzyl * dkperp_dkalp
            dkwzz = dkwzz + dkwzzl * dkperp_dkalp           

         end do
         
         

         

c      end if
      
      go to 5600      

*     -----------------------------------
*     Funky METS calls for non Maxwellian:
*     -----------------------------------
      if (ndist .eq. 1) then
      
         if (upshift .ne. 0) xkprl = xkprl_eff
         
      
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
         DensSPEC = xn
         bratio = bmod0 / bmod
         if(bratio .gt. 1.0) bratio = 1.0

         duperp = uperp(nuper) / (nuper - 1)
         dupara = 2.0 * upara(nupar) / (nupar - 1)

         if(i .ne. i_sav .or. j .ne. j_sav)then
         
            dfduper0 = 0.0
            dfdupar0 = 0.0

!           ------------------------------------------------
!           get CQL3D distribution function on the midplane
!           ------------------------------------------------
            call cql3d_dist(nupar, nuper, n_psi,
     .                 n_psi_dim, rho_a, rho,
     .                 UminPara,UmaxPara,
     .                 df_cql_uprp, df_cql_uprl,
     .                 UPERP, UPARA, DFDUPER0, DFDUPAR0)
     

!           ------------------------------------------------
!           map CQL3D distribution function off the midplane
!           ------------------------------------------------

            if(bratio .ge. 0.0)then
            
               dfduper = 0.0
               dfdupar = 0.0
            
               do ni = 1, nuper
                  do mi = 1, nupar

                     argd = uperp(ni)**2 * (1. - bratio) 
     .                                               + upara(mi)**2
                     if (argd .le. 0.0) argd = 1.0e-06
                     
                     uperp0 = uperp(ni) * sqrt(bratio)
                     upara0 = sign(1.0, upara(mi)) * sqrt(argd)
                     
                     dfduper(ni, mi) = 0.0
                     dfdupar(ni, mi) = 0.0
                     
                     if(upara0 .ge. upara(1) .and. 
     .                                   upara0 .le. upara(nupar)) then
     
                        ni0 = int((uperp0 - uperp(1)) / duperp) + 1
                        mi0 = int((upara0 - upara(1)) / dupara) + 1
                        
                        dfduper0_intplt = dfduper0(ni0, mi0)
                        dfdupar0_intplt = dfdupar0(ni0, mi0)
                        
                        if (ni0 .lt. nuper .and. mi0 .lt. nupar) then
                        
                        uperp0_grid = uperp(1) + (ni0 - 1) * duperp
                        upara0_grid = upara(1) + (mi0 - 1) * dupara
                                                
                        zeta = (uperp0 - uperp0_grid) / duperp
                        eta  = (upara0 - upara0_grid) / dupara
                        
                        ai = dfduper0(ni0, mi0)
                        bi = dfduper0(ni0+1 ,mi0) - dfduper0(ni0, mi0)
                        ci = dfduper0(ni0, mi0+1) - dfduper0(ni0, mi0)
                        di = dfduper0(ni0+1, mi0+1)+ dfduper0(ni0, mi0) 
     .                     - dfduper0(ni0+1, mi0) - dfduper0(ni0, mi0+1) 
                           
                        dfduper0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta                         

                        ai = dfdupar0(ni0, mi0)
                        bi = dfdupar0(ni0+1 ,mi0) - dfdupar0(ni0, mi0)
                        ci = dfdupar0(ni0, mi0+1) - dfdupar0(ni0, mi0)
                        di = dfdupar0(ni0+1, mi0+1)+ dfdupar0(ni0, mi0) 
     .                     - dfdupar0(ni0+1, mi0) - dfdupar0(ni0, mi0+1) 
                           
                        dfdupar0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta
     
                        end if                  

                        
                        if (upara0 .ne. 0.0)then
                                             
                           dfdupar(ni, mi) = dfdupar0_intplt * 
     .                        upara(mi) / upara0
     
                           dfduper(ni, mi) = dfduper0_intplt * 
     .                        sqrt(bratio) + dfdupar0_intplt * 
     .                        uperp(ni) / upara0 * (1.0 - bratio)
     
c                          dfdth = upara(mi) * dfduper(ni, mi)
c     .                           - uperp(ni) * dfdupar(ni, mi)


                        end if
                        
                     end if
                     
                     
                     go to 5000
!                    ----------------------------
!                    optional analytic Maxwellian
!                    ----------------------------
                     pi = 3.141592654
                     
                     alpha = sqrt(2.0 * xkt / xm)
!                     vc_mks = 3.5 * alpha
                     u0 = vc_mks / alpha
                     
                     fnorm = u0**3 / pi**1.5 
                           
                     u2 = uperp(ni)**2 + upara(mi)**2
                     
                     f_cql = exp(-u2 * u0**2) * fnorm
                     dfduper(ni, mi) = -f_cql * 2. * uperp(ni) * u0**2
                     dfdupar(ni, mi) = -f_cql * 2. * upara(mi) * u0**2
 5000                continue                


                  end do
               end do
               
            end if
            

!           --------------------------------------
!           Initialize the interpolation in k_perp:
!           --------------------------------------
            if (nkperp .ne. 0) then
               l_first  = .true.
               l_interp = .true.
        
               call GETNONMAXSIGMA_AORSA_NEWi(W,
     .                          ZSPEC,ASPEC,DensSPEC,BMAG,
     .                          K1,XI1,JNXI1,
     .                          K1,XI1,JNXI1,NBESSJ,
     .                          Enorm,UminPara,UmaxPara,
     .                          NUPAR,NUPER,UPERP,UPARA,
     .                          DFDUPER,DFDUPAR,
     .                          WSPEC,IFAIL,
     .                          l_first, l_interp, kperp_max, nkperp,
     .                          xkphi)
            end if



            i_sav = i
            j_sav = j

         end if






!        ------------------------------------
!        Complete integrals; no interpolation
!        ------------------------------------
         if (nkperp .eq. 0)then

            call WMATPRECALC_AORSA(ZSPEC,ASPEC,ENORM,BMAG,KPER1,UPERP,
     &                   NUPER,NBESSJ,NSBESSJ,XI1,JNXI1,IFAIL)

            call GETNONMAX_SIGMA_AORSA_NEW(W,
     .                          ZSPEC, ASPEC, DensSPEC, BMAG,
     .                          K1, XI1, JNXI1,
     .                          K1, XI1, JNXI1, NBESSJ,
     .                          Enorm, UminPara, UmaxPara,
     .                          NUPAR, NUPER, UPERP, UPARA,
     .                          DFDUPER, DFDUPAR,
     .                          WSPEC, IFAIL)

         else
!        -----------------
!        Use interpolation
!        -----------------

            l_first = .false.

              call GETNONMAXSIGMA_AORSA_NEWi(W,
     .                       ZSPEC,ASPEC,DensSPEC,BMAG,
     .                       K1, XI1, JNXI1,
     .                       K1, XI1, JNXI1, NBESSJ,
     .                       Enorm, UminPara, UmaxPara,
     .                       NUPAR, NUPER, UPERP, UPARA,
     .                       DFDUPER, DFDUPAR,
     .                       WSPEC, IFAIL,
     .                       l_first, l_interp, kperp_max, nkperp,
     .                       xkphi)
         end if



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

      end if

      sig1 = sig1 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
      sig3 = sig3 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
            
      
      kr = xkperp / xk_cutoff
      step = damping * kr**16 / (1. + kr**16)       
      if (xm .eq. xme) sig3 = sig3 * (1.0 + step)            
                                      
      

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
      
 5600 continue      
      

      deallocate( dfduper0 )
      deallocate( dfdupar0 )
      



      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
  102 format(2i10, 1p8e12.4)
  103 format(4i10, 1p8e12.4)
      end


c
c***************************************************************************
c

      subroutine delta_(i, j, n, m, rho, rho_a,
     .   delta_x, delta_y, delta_z,
     .   gradprlb, bmod, bmod0,
     .   xm, q, xn, xnuomg,
     .   xkt, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, nphi, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   delta0, ndist, nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .   nkperp, zi, eps0, v0i, omgrf, xk0, 
     .   kperp_max, i_sav, j_sav, upshift, 
     .   damping, xk_cutoff, rt, nkx2, nky2, xkprl, xkperp, 
     .   xkalp, xkbet, 
     .   z0_table, z1_table, z2_table, zetai_table, dKdL_table, 
     .   dKdL_giv, nmax, mmax,  use_new_z2, ntable, mtable)


*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------
c      use zfun_hilbert  
      
      implicit none
      
      integer :: ieer ! fft error flag
      
      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi, ndist, myid, iflag, nproc, ni, mi,
     .    i_sav, j_sav, ni0, mi0, upshift, nkx2, nky2
     
      integer ntable, mtable, nmax, mmax     
      
      complex z0_table(ntable, mtable)
      complex z1_table(ntable, mtable)
      complex z2_table(ntable, mtable)            
      real zetai_table(ntable), dKdL_table(mtable), dKdL_giv     
     
      real u0, u2, fnorm, f_cql, rt

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, y, sgn_kprl, reson, duperp, dupara
      real xkprl_eff0
      real dzetal(lmin:lmax), descrim
      real dakbdkb, xnuomg, gradprlb, bmod, bmod0, nu_coll
      real akprl,  rho, alpha, eps0, omgrf, v0i, emax, akprl_min
      real gammab(lmin:lmax), gamma_coll(lmin:lmax)
      real a, b, xnurf, pi, delta0, rhol
      real bx, by, bz, bratio, denom
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt

      real xkxsav, xkysav, capr, argd
      real xkphi
      real xkalp, xkbet, xk0, rgamma, xk_cutoff, damping, kr, step

      complex zi, zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2
      complex z0_new, z1_new, z2_new

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0_a, sig1_a, sig2_a, sig3_a, sig4_a, sig5_a
      complex sig0_h, sig1_h, sig2_h, sig3_h, sig4_h, sig5_h

      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
      logical :: l_interp   !new
      logical :: l_first  !new
      logical :: use_new_z2
      integer :: nkperp
      real :: kperp_max !new
      
      complex delta_x,  delta_y,  delta_z
      complex delta_xl, delta_yl, delta_zl

      complex sigxx, sigxy, sigxz,
     .        sigyx, sigyy, sigyz,
     .        sigzx, sigzy, sigzz

      real uxx, uxy, uxz,
     .     uyx, uyy, uyz,
     .     uzx, uzy, uzz

      parameter (lmaxdim = 99)

      complex xil(0: lmaxdim), xilp(0: lmaxdim)
      complex exil(0: lmaxdim), exilp(0: lmaxdim),
     .                    exilovergam(0: lmaxdim)

      complex zetal(lmin:lmax)
      complex  zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      integer  :: n_psi_dim
      integer :: nuper, nupar, n_psi

      integer NBESSJ

      real :: UPERP(NUPER),UPARA(NUPAR)
      real :: UPERP0, UPARA0

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
      
      common/upcom/akprl_min
      nu_coll =  .01 * omgrf
      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rhol = alpha / omgc
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * xnuomg)
      
      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl)       

      
*     ---------------------------------
*     Calculate zetal(l) and gammab(l)
*     ---------------------------------      
      
      do l = lmin, lmax
         labs = abs(l)

         reson = (omgrf - l * real(omgc)) / omgrf
c         if (abs(reson) .lt. 0.02)then
         if (rho .gt. 1.0) then
            zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
            dzetal(l) = omgrf * xnuomg / (xkprl * alpha)
         else
            zetal(l) = (omgrf  - l * omgc) / (xkprl * alpha)
            dzetal(l) = 0.0
         end if
         
c        zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
c         dzetal(l) = omgrf * xnuomg / (xkprl * alpha)


         gammab(l) = abs(l * omgc / (2.0 * alpha * xkprl**2)
     .                                           * gradprlb / bmod)
         gamma_coll(l) = nu_coll / (akprl * alpha)


         if(xm .eq. xme)gammab(l) = 0.0
c         if(abs(gammab(l)) .gt. 1000.0) gammab(l) = 1000.0
         if(abs(gammab(l)) .lt. .01)gammab(l) = .01


      enddo

      
      
*     ------------------------------------------------
*     Calculate Brambilla's xkrpl_eff using l = 1 only
*     ------------------------------------------------
      y0 = 1.5
      y = y0
      

      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gammab(1) * y) - 1.)
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            descrim = 1. - 4. * gammab(1) * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gammab(1) * y) )
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if
                    

*     -----------------------
*     Maxwellian electrons
*     -----------------------

         gamma = 0.5 * xkperp**2 * rhol**2
         rgamma = real(gamma)


         if(rgamma .ge. 1.0e-08)
     .      call besiexp(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

         if(rgamma .lt. 1.0e-08)
     .      call bes_expand(gamma, lmax, exil, exilp, lmaxdim,
     .                                                      exilovergam)
         
         delta_x = 0.0 
         delta_y = 0.0 
         delta_z = 0.0

         do l = lmin, lmax
            labs = abs(l)

            if(nzfun .eq. 0) call z_approx(sgn_kprl, zetal(l), 0.0,
     .                                                     z0, z1, z2)
            if(nzfun .eq. 1) then
               call z_approx(sgn_kprl, zetal(l), gammab(l), z0, z1, z2)
               
               if(use_new_z2 .eqv. .true. .and. l .eq. 0) then
                  call z_approx_e(sgn_kprl, zetal(l), gammab(l),
     .               z0_new, z1_new, z2_new, zetai_table, 
     .               dKdL_table, z0_table, z1_table, z2_table, 
     .               dKdL_giv, nmax, mmax, ntable, mtable)
                  z0 = z0_new
                  z1 = z1_new
                  z2 = z2_new
               end if
     
            end if
            
            if(nzfun .eq. 2) call z_smithe(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
            if(nzfun .eq. 3) call z_table(sgn_kprl,zetal(l),gammab(l),
     .                                      gamma_coll(l), z0, z1, z2)
                    
            delta_xl = - zieps0 / q * omgp2 / omgrf**2 * zetal(0)
     .         * xkperp * l * exilovergam(labs)  * omgrf / omgc * z0
            
            delta_yl = - zieps0 / q * omgp2 / omgrf**2 * zetal(0)
     .         * zi * xkperp * omgrf /omgc * (exilp(labs) - exil(labs))
     .         * z0     
         
            delta_zl = - zieps0 / q * omgp2 / omgrf**2 * zetal(0)
     .         * xkprl * 2.0 * zetal(0) * exil(labs) * z1       
                            
            delta_x = delta_x + delta_xl 
            delta_y = delta_y + delta_yl 
            delta_z = delta_z + delta_zl            

         end do


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
  102 format(2i10, 1p8e12.4)
  103 format(4i10, 1p8e12.4)
      end
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
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .   nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .   i_sav, j_sav, upshift, damping, xk_cutoff, 
     .   rt, nkx2, nky2, xkprl, xkperp, xkalp, xkbet)


*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------
c      use zfun_hilbert  
      
      implicit none
      
      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0

      integer :: ieer ! fft error flag 

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi, ndist, myid, iflag, nproc, ni, mi,
     .    i_sav, j_sav, ni0, mi0, upshift, nkx2, nky2
      integer n_upper, n_lower, m_upper, m_lower     
      
     
      real u0, u2, fnorm, f_cql, rt
      real fnmax, fmmax, nperp
      real fna, fma, fstepn, fstepm, fstep

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, y, sgn_kprl, reson, duperp, dupara
      real xkprl_eff0
      real dzetal(lmin:lmax), descrim
      real dakbdkb, xnuomg, gradprlb, bmod, bmod0, nu_coll
      real akprl,  rho, alpha, eps0, omgrf, v0i, emax, akprl_min
      real gammab(lmin:lmax), gamma_coll(lmin:lmax)
      real a, b, xnurf, pi, delta0, rhol
      real bx, by, bz, bratio, denom
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt

      real xkxsav, xkysav, capr, argd
      real xkphi
      real xkalp, xkbet, xk0, rgamma, xk_cutoff, damping, kr, step

      complex zi, zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0_a, sig1_a, sig2_a, sig3_a, sig4_a, sig5_a
      complex sig0_h, sig1_h, sig2_h, sig3_h, sig4_h, sig5_h

      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
      logical :: l_interp   !new
      logical :: l_first  !new
      integer :: nkperp
      real :: kperp_max !new


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

      complex zetal(lmin:lmax)
      complex  zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      integer  :: n_psi_dim
      integer :: nuper, nupar, n_psi

      integer NBESSJ

      real :: UPERP(NUPER),UPARA(NUPAR)
      real :: UPERP0, UPARA0

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

      common/upcom/akprl_min
            
      allocate( dfduper0(nuper, nupar) )
      allocate( dfdupar0(nuper, nupar) )

      nu_coll =  .01 * omgrf
      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rhol = alpha / omgc
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * xnuomg)
      
      
      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl)              


c      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
c      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
c      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi
c      xkperp = sqrt(xkalp**2 + xkbet**2)
      

      
*     ---------------------------------
*     Calculate zetal(l) and gammab(l)
*     ---------------------------------      
      
      do l = lmin, lmax
         labs = abs(l)

         reson = (omgrf - l * real(omgc)) / omgrf
c         if (abs(reson) .lt. 0.02)then
         if (rho .gt. 1.0) then
            zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
            dzetal(l) = omgrf * xnuomg / (xkprl * alpha)
         else
            zetal(l) = (omgrf  - l * omgc) / (xkprl * alpha)
            dzetal(l) = 0.0
         end if
         
c        zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
c         dzetal(l) = omgrf * xnuomg / (xkprl * alpha)


         gammab(l) = abs(l * omgc / (2.0 * alpha * xkprl**2)
     .                                           * gradprlb / bmod)
         gamma_coll(l) = nu_coll / (akprl * alpha)


         if(xm .eq. xme)gammab(l) = 0.0
c         if(abs(gammab(l)) .gt. 1000.0) gammab(l) = 1000.0
         if(abs(gammab(l)) .lt. .01)gammab(l) = .01


      enddo

      
      
*     ------------------------------------------------
*     Calculate Brambilla's xkrpl_eff using l = 1 only
*     ------------------------------------------------
      y0 = 1.5
      y = y0
      

      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gammab(1) * y) - 1.)
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            descrim = 1. - 4. * gammab(1) * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gammab(1) * y) )
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if
                    

*     -----------------------
*     Maxwellian distribution
*     -----------------------

      if(ndist .eq. 0)then

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
c         do l = 0, 0
            labs = abs(l)

           if(nzfun .eq. 0) call z_approx(sgn_kprl, zetal(l), 0.0,
     .                                                     z0, z1, z2)
           if(nzfun .eq. 1) call z_approx(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 2) call z_smithe(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 3) call z_table(sgn_kprl,zetal(l),gammab(l),
     .                                      gamma_coll(l), z0, z1, z2)


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

      end if

*     -----------------------------------
*     Funky METS calls for non Maxwellian:
*     -----------------------------------
      if (ndist .eq. 1) then
      
         if (upshift .ne. 0) xkprl = xkprl_eff
         
      
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
         DensSPEC = xn
         bratio = bmod0 / bmod
         if(bratio .gt. 1.0) bratio = 1.0

         duperp = uperp(nuper) / (nuper - 1)
         dupara = 2.0 * upara(nupar) / (nupar - 1)

         if(i .ne. i_sav .or. j .ne. j_sav)then
         
            dfduper0 = 0.0
            dfdupar0 = 0.0

!           ------------------------------------------------
!           get CQL3D distribution function on the midplane
!           ------------------------------------------------
            call cql3d_dist(nupar, nuper, n_psi,
     .                 n_psi_dim, rho_a, rho,
     .                 UminPara,UmaxPara,
     .                 df_cql_uprp, df_cql_uprl,
     .                 UPERP, UPARA, DFDUPER0, DFDUPAR0)
     

!           ------------------------------------------------
!           map CQL3D distribution function off the midplane
!           ------------------------------------------------

            if(bratio .ge. 0.0)then
            
               dfduper = 0.0
               dfdupar = 0.0
            
               do ni = 1, nuper
                  do mi = 1, nupar

                     argd = uperp(ni)**2 * (1. - bratio) 
     .                                               + upara(mi)**2
                     if (argd .le. 0.0) argd = 1.0e-06
                     
                     uperp0 = uperp(ni) * sqrt(bratio)
                     upara0 = sign(1.0, upara(mi)) * sqrt(argd)
                     
                     dfduper(ni, mi) = 0.0
                     dfdupar(ni, mi) = 0.0
                     
                     if(upara0 .ge. upara(1) .and. 
     .                                   upara0 .le. upara(nupar)) then
     
                        ni0 = int((uperp0 - uperp(1)) / duperp) + 1
                        mi0 = int((upara0 - upara(1)) / dupara) + 1
                        
                        dfduper0_intplt = dfduper0(ni0, mi0)
                        dfdupar0_intplt = dfdupar0(ni0, mi0)
                        
                        if (ni0 .lt. nuper .and. mi0 .lt. nupar) then
                        
                        uperp0_grid = uperp(1) + (ni0 - 1) * duperp
                        upara0_grid = upara(1) + (mi0 - 1) * dupara
                                                
                        zeta = (uperp0 - uperp0_grid) / duperp
                        eta  = (upara0 - upara0_grid) / dupara
                        
                        ai = dfduper0(ni0, mi0)
                        bi = dfduper0(ni0+1 ,mi0) - dfduper0(ni0, mi0)
                        ci = dfduper0(ni0, mi0+1) - dfduper0(ni0, mi0)
                        di = dfduper0(ni0+1, mi0+1)+ dfduper0(ni0, mi0) 
     .                     - dfduper0(ni0+1, mi0) - dfduper0(ni0, mi0+1) 
                           
                        dfduper0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta                         

                        ai = dfdupar0(ni0, mi0)
                        bi = dfdupar0(ni0+1 ,mi0) - dfdupar0(ni0, mi0)
                        ci = dfdupar0(ni0, mi0+1) - dfdupar0(ni0, mi0)
                        di = dfdupar0(ni0+1, mi0+1)+ dfdupar0(ni0, mi0) 
     .                     - dfdupar0(ni0+1, mi0) - dfdupar0(ni0, mi0+1) 
                           
                        dfdupar0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta
     
                        end if                  

                        
                        if (upara0 .ne. 0.0)then
                                             
                           dfdupar(ni, mi) = dfdupar0_intplt * 
     .                        upara(mi) / upara0
     
                           dfduper(ni, mi) = dfduper0_intplt * 
     .                        sqrt(bratio) + dfdupar0_intplt * 
     .                        uperp(ni) / upara0 * (1.0 - bratio)
     
c                          dfdth = upara(mi) * dfduper(ni, mi)
c     .                           - uperp(ni) * dfdupar(ni, mi)


                        end if
                        
                     end if
                     
                     
                     go to 5000
!                    ----------------------------
!                    optional analytic Maxwellian
!                    ----------------------------
                     pi = 3.141592654
                     
                     alpha = sqrt(2.0 * xkt / xm)
!                     vc_mks = 3.5 * alpha
                     u0 = vc_mks / alpha
                     
                     fnorm = u0**3 / pi**1.5 
                           
                     u2 = uperp(ni)**2 + upara(mi)**2
                     
                     f_cql = exp(-u2 * u0**2) * fnorm
                     dfduper(ni, mi) = -f_cql * 2. * uperp(ni) * u0**2
                     dfdupar(ni, mi) = -f_cql * 2. * upara(mi) * u0**2
 5000                continue                


                  end do
               end do
               
            end if
            

!           --------------------------------------
!           Initialize the interpolation in k_perp:
!           --------------------------------------
            if (nkperp .ne. 0) then
               l_first  = .true.
               l_interp = .true.
        
               call GETNONMAXSIGMA_AORSA_NEWi(W,
     .                          ZSPEC,ASPEC,DensSPEC,BMAG,
     .                          K1,XI1,JNXI1,
     .                          K1,XI1,JNXI1,NBESSJ,
     .                          Enorm,UminPara,UmaxPara,
     .                          NUPAR,NUPER,UPERP,UPARA,
     .                          DFDUPER,DFDUPAR,
     .                          WSPEC,IFAIL,
     .                          l_first, l_interp, kperp_max, nkperp,
     .                          xkphi)
            end if



            i_sav = i
            j_sav = j

         end if






!        ------------------------------------
!        Complete integrals; no interpolation
!        ------------------------------------
         if (nkperp .eq. 0)then

            call WMATPRECALC_AORSA(ZSPEC,ASPEC,ENORM,BMAG,KPER1,UPERP,
     &                   NUPER,NBESSJ,NSBESSJ,XI1,JNXI1,IFAIL)

            call GETNONMAX_SIGMA_AORSA_NEW(W,
     .                          ZSPEC, ASPEC, DensSPEC, BMAG,
     .                          K1, XI1, JNXI1,
     .                          K1, XI1, JNXI1, NBESSJ,
     .                          Enorm, UminPara, UmaxPara,
     .                          NUPAR, NUPER, UPERP, UPARA,
     .                          DFDUPER, DFDUPAR,
     .                          WSPEC, IFAIL)

         else
!        -----------------
!        Use interpolation
!        -----------------

            l_first = .false.

              call GETNONMAXSIGMA_AORSA_NEWi(W,
     .                       ZSPEC,ASPEC,DensSPEC,BMAG,
     .                       K1, XI1, JNXI1,
     .                       K1, XI1, JNXI1, NBESSJ,
     .                       Enorm, UminPara, UmaxPara,
     .                       NUPAR, NUPER, UPERP, UPARA,
     .                       DFDUPER, DFDUPAR,
     .                       WSPEC, IFAIL,
     .                       l_first, l_interp, kperp_max, nkperp,
     .                       xkphi)
         end if



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

      end if

      sig1 = sig1 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
      sig3 = sig3 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
            
      
      kr = xkperp / xk_cutoff
      step = damping * kr**16 / (1. + kr**16)       
      if (xm .eq. xme) sig3 = sig3 * (1.0 + step)
                      
      

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
      

      deallocate( dfduper0 )
      deallocate( dfdupar0 )


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
  102 format(2i10, 1p8e12.4)
  103 format(4i10, 1p8e12.4)
      end


c
c***************************************************************************
c

      subroutine sigmad_cql3d_e(i, j, n, m, rho, rho_a,
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
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .   nkperp, zi, eps0, v0i, omgrf, xk0, 
     .   kperp_max, i_sav, j_sav, upshift, damping, 
     .   xk_cutoff, rt, nkx2, nky2, xkprl, xkperp, 
     .   xkalp, xkbet, 
     .   z0_table, z1_table, z2_table, zetai_table, dKdL_table, 
     .   dKdL_giv, nmax, mmax,  use_new_z2, ntable, mtable, sige3_new)

     


*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------
c      use zfun_hilbert  
      
      implicit none
      
      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0

      integer :: ieer ! fft error flag 

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi, ndist, myid, iflag, nproc, ni, mi,
     .    i_sav, j_sav, ni0, mi0, upshift, nkx2, nky2
      integer n_upper, n_lower, m_upper, m_lower, nmax, mmax
      integer ntable, mtable     
      
      complex z0_table(ntable, mtable)
      complex z1_table(ntable, mtable)
      complex z2_table(ntable, mtable)            
      real zetai_table(ntable), dKdL_table(mtable), dKdL_giv
      
      
      real u0, u2, fnorm, f_cql, rt
      real fnmax, fmmax, nperp
      real fna, fma, fstepn, fstepm, fstep

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, y, sgn_kprl, reson, duperp, dupara
      real sgn_dkdl, abs_dkdl
      real xkprl_eff0
      real dzetal(lmin:lmax), descrim
      real dakbdkb, xnuomg, gradprlb, bmod, bmod0, nu_coll
      real akprl,  rho, alpha, eps0, omgrf, v0i, emax, akprl_min
      real gammab(lmin:lmax), gamma_coll(lmin:lmax)
      real a, b, xnurf, pi, delta0, rhol
      real bx, by, bz, bratio, denom
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt

      real xkxsav, xkysav, capr, argd
      real xkphi
      real xkalp, xkbet, xk0, rgamma, xk_cutoff, damping, kr, step

      complex zi, zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2
      complex z0_new, z1_new, z2_new

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0_a, sig1_a, sig2_a, sig3_a, sig4_a, sig5_a
      complex sig0_h, sig1_h, sig2_h, sig3_h, sig4_h, sig5_h

      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
      complex sige3_new
      logical :: l_interp   !new
      logical :: l_first  !new
      logical :: use_new_z2
      integer :: nkperp
      real :: kperp_max !new


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

      complex zetal(lmin:lmax)
      complex  zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      integer  :: n_psi_dim
      integer :: nuper, nupar, n_psi

      integer NBESSJ

      real :: UPERP(NUPER),UPARA(NUPAR)
      real :: UPERP0, UPARA0

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

      common/upcom/akprl_min
            
      allocate( dfduper0(nuper, nupar) )
      allocate( dfdupar0(nuper, nupar) )

      nu_coll =  .01 * omgrf
      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rhol = alpha / omgc
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * xnuomg)

      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl) 
      
      sgn_dKdL = sign(1.0, dKdL_giv)
      abs_dKdL = abs(dKdL_giv)
             
!     -------------------------------------------------------------
!     Optional: Don't allow xkprl = 0 and dKdL = 0 at the same time 
!     -------------------------------------------------------------       
c      if (akprl .lt. akprl_min .and. abs_dKdL .lt. 0.1) then
c         xkprl = akprl_min * sgn_kprl
c        dkdL_giv = 0.1 * sgn_dKdL
c      end if                 
            
      
*     ---------------------------------
*     Calculate zetal(l) and gammab(l)
*     ---------------------------------      
      
      do l = lmin, lmax
         labs = abs(l)

         reson = (omgrf - l * real(omgc)) / omgrf
c         if (abs(reson) .lt. 0.02)then
         if (rho .gt. 1.0) then
            zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
            dzetal(l) = omgrf * xnuomg / (xkprl * alpha)
         else
            zetal(l) = (omgrf  - l * omgc) / (xkprl * alpha)
            dzetal(l) = 0.0
         end if
         
c        zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
c         dzetal(l) = omgrf * xnuomg / (xkprl * alpha)


         gammab(l) = abs(l * omgc / (2.0 * alpha * xkprl**2)
     .                                           * gradprlb / bmod)
         gamma_coll(l) = nu_coll / (akprl * alpha)


         if(xm .eq. xme)gammab(l) = 0.0
c         if(abs(gammab(l)) .gt. 1000.0) gammab(l) = 1000.0
c         if(abs(gammab(l)) .lt. .01)gammab(l) = .01


      enddo

      
      
*     ------------------------------------------------
*     Calculate Brambilla's xkrpl_eff using l = 1 only
*     ------------------------------------------------
      y0 = 1.5
      y = y0
      

      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gammab(1) * y) - 1.)
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            descrim = 1. - 4. * gammab(1) * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gammab(1) * y) )
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if
                    

*     -----------------------
*     Maxwellian distribution
*     -----------------------

      if(ndist .eq. 0)then

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

           if(nzfun .eq. 0) call z_approx(sgn_kprl, zetal(l), 0.0,
     .                                                     z0, z1, z2)
     
           if(nzfun .eq. 1) then
              call z_approx(sgn_kprl, zetal(l), gammab(l), z0, z1, z2)
     
              if(use_new_z2 .eqv. .true. .and. l .eq. 0) then      
                 call z_approx_e(sgn_kprl, zetal(l), gammab(l),
     .              z0_new, z1_new, z2_new, zetai_table, 
     .              dKdL_table, z0_table, z1_table, z2_table, 
     .              dKdL_giv, nmax, mmax, ntable, mtable)
                 z0 = z0_new
                 z1 = z1_new
                 z2 = z2_new
              end if
           end if
                    
         
           if(nzfun .eq. 2) call z_smithe(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 3) call z_table(sgn_kprl,zetal(l),gammab(l),
     .                                      gamma_coll(l), z0, z1, z2)


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
     
*           -----------------------------------------------------------------
*           Replace the l = 0 term with the Fourier expanded term if non-zero
*           -----------------------------------------------------------------           
            if (l .eq. 0 .and. sige3_new .ne. 0.0) sig3l = sige3_new     


            sig0 = sig0 + sig0l
            sig1 = sig1 + sig1l
            sig2 = sig2 + sig2l
            sig3 = sig3 + sig3l
            sig4 = sig4 + sig4l
            sig5 = sig5 + sig5l

         end do

      end if

*     -----------------------------------
*     Funky METS calls for non Maxwellian:
*     -----------------------------------
      if (ndist .eq. 1) then
      
         if (upshift .ne. 0) xkprl = xkprl_eff
         
      
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
         DensSPEC = xn
         bratio = bmod0 / bmod
         if(bratio .gt. 1.0) bratio = 1.0

         duperp = uperp(nuper) / (nuper - 1)
         dupara = 2.0 * upara(nupar) / (nupar - 1)

         if(i .ne. i_sav .or. j .ne. j_sav)then
         
            dfduper0 = 0.0
            dfdupar0 = 0.0

!           ------------------------------------------------
!           get CQL3D distribution function on the midplane
!           ------------------------------------------------
            call cql3d_dist(nupar, nuper, n_psi,
     .                 n_psi_dim, rho_a, rho,
     .                 UminPara,UmaxPara,
     .                 df_cql_uprp, df_cql_uprl,
     .                 UPERP, UPARA, DFDUPER0, DFDUPAR0)
     

!           ------------------------------------------------
!           map CQL3D distribution function off the midplane
!           ------------------------------------------------

            if(bratio .ge. 0.0)then
            
               dfduper = 0.0
               dfdupar = 0.0
            
               do ni = 1, nuper
                  do mi = 1, nupar

                     argd = uperp(ni)**2 * (1. - bratio) 
     .                                               + upara(mi)**2
                     if (argd .le. 0.0) argd = 1.0e-06
                     
                     uperp0 = uperp(ni) * sqrt(bratio)
                     upara0 = sign(1.0, upara(mi)) * sqrt(argd)
                     
                     dfduper(ni, mi) = 0.0
                     dfdupar(ni, mi) = 0.0
                     
                     if(upara0 .ge. upara(1) .and. 
     .                                   upara0 .le. upara(nupar)) then
     
                        ni0 = int((uperp0 - uperp(1)) / duperp) + 1
                        mi0 = int((upara0 - upara(1)) / dupara) + 1
                        
                        dfduper0_intplt = dfduper0(ni0, mi0)
                        dfdupar0_intplt = dfdupar0(ni0, mi0)
                        
                        if (ni0 .lt. nuper .and. mi0 .lt. nupar) then
                        
                        uperp0_grid = uperp(1) + (ni0 - 1) * duperp
                        upara0_grid = upara(1) + (mi0 - 1) * dupara
                                                
                        zeta = (uperp0 - uperp0_grid) / duperp
                        eta  = (upara0 - upara0_grid) / dupara
                        
                        ai = dfduper0(ni0, mi0)
                        bi = dfduper0(ni0+1 ,mi0) - dfduper0(ni0, mi0)
                        ci = dfduper0(ni0, mi0+1) - dfduper0(ni0, mi0)
                        di = dfduper0(ni0+1, mi0+1)+ dfduper0(ni0, mi0) 
     .                     - dfduper0(ni0+1, mi0) - dfduper0(ni0, mi0+1) 
                           
                        dfduper0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta                         

                        ai = dfdupar0(ni0, mi0)
                        bi = dfdupar0(ni0+1 ,mi0) - dfdupar0(ni0, mi0)
                        ci = dfdupar0(ni0, mi0+1) - dfdupar0(ni0, mi0)
                        di = dfdupar0(ni0+1, mi0+1)+ dfdupar0(ni0, mi0) 
     .                     - dfdupar0(ni0+1, mi0) - dfdupar0(ni0, mi0+1) 
                           
                        dfdupar0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta
     
                        end if                  

                        
                        if (upara0 .ne. 0.0)then
                                             
                           dfdupar(ni, mi) = dfdupar0_intplt * 
     .                        upara(mi) / upara0
     
                           dfduper(ni, mi) = dfduper0_intplt * 
     .                        sqrt(bratio) + dfdupar0_intplt * 
     .                        uperp(ni) / upara0 * (1.0 - bratio)
     
c                          dfdth = upara(mi) * dfduper(ni, mi)
c     .                           - uperp(ni) * dfdupar(ni, mi)


                        end if
                        
                     end if
                     
                     
                     go to 5000
!                    ----------------------------
!                    optional analytic Maxwellian
!                    ----------------------------
                     pi = 3.141592654
                     
                     alpha = sqrt(2.0 * xkt / xm)
!                     vc_mks = 3.5 * alpha
                     u0 = vc_mks / alpha
                     
                     fnorm = u0**3 / pi**1.5 
                           
                     u2 = uperp(ni)**2 + upara(mi)**2
                     
                     f_cql = exp(-u2 * u0**2) * fnorm
                     dfduper(ni, mi) = -f_cql * 2. * uperp(ni) * u0**2
                     dfdupar(ni, mi) = -f_cql * 2. * upara(mi) * u0**2
 5000                continue                


                  end do
               end do
               
            end if
            

!           --------------------------------------
!           Initialize the interpolation in k_perp:
!           --------------------------------------
            if (nkperp .ne. 0) then
               l_first  = .true.
               l_interp = .true.
        
               call GETNONMAXSIGMA_AORSA_NEWi(W,
     .                          ZSPEC,ASPEC,DensSPEC,BMAG,
     .                          K1,XI1,JNXI1,
     .                          K1,XI1,JNXI1,NBESSJ,
     .                          Enorm,UminPara,UmaxPara,
     .                          NUPAR,NUPER,UPERP,UPARA,
     .                          DFDUPER,DFDUPAR,
     .                          WSPEC,IFAIL,
     .                          l_first, l_interp, kperp_max, nkperp,
     .                          xkphi)
            end if



            i_sav = i
            j_sav = j

         end if






!        ------------------------------------
!        Complete integrals; no interpolation
!        ------------------------------------
         if (nkperp .eq. 0)then

            call WMATPRECALC_AORSA(ZSPEC,ASPEC,ENORM,BMAG,KPER1,UPERP,
     &                   NUPER,NBESSJ,NSBESSJ,XI1,JNXI1,IFAIL)

            call GETNONMAX_SIGMA_AORSA_NEW(W,
     .                          ZSPEC, ASPEC, DensSPEC, BMAG,
     .                          K1, XI1, JNXI1,
     .                          K1, XI1, JNXI1, NBESSJ,
     .                          Enorm, UminPara, UmaxPara,
     .                          NUPAR, NUPER, UPERP, UPARA,
     .                          DFDUPER, DFDUPAR,
     .                          WSPEC, IFAIL)

         else
!        -----------------
!        Use interpolation
!        -----------------

            l_first = .false.

              call GETNONMAXSIGMA_AORSA_NEWi(W,
     .                       ZSPEC,ASPEC,DensSPEC,BMAG,
     .                       K1, XI1, JNXI1,
     .                       K1, XI1, JNXI1, NBESSJ,
     .                       Enorm, UminPara, UmaxPara,
     .                       NUPAR, NUPER, UPERP, UPARA,
     .                       DFDUPER, DFDUPAR,
     .                       WSPEC, IFAIL,
     .                       l_first, l_interp, kperp_max, nkperp,
     .                       xkphi)
         end if



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

      end if

      sig1 = sig1 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
      sig3 = sig3 + delta0 * eps0 * omgrf * xkperp**2 / xk0**2
            
      
      kr = xkperp / xk_cutoff
      step = damping * kr**16 / (1. + kr**16)       
      if (xm .eq. xme) sig3 = sig3 * (1.0 + step)

              
      

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
      

      deallocate( dfduper0 )
      deallocate( dfdupar0 )


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
  102 format(2i10, 1p8e12.4)
  103 format(4i10, 1p8e12.4)
      end


c
c***************************************************************************
c

      subroutine sigmad_cql3d_1(i, j, n, m, rho, rho_a,
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
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .   nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max, 
     .   i_sav, j_sav, upshift, damping, xk_cutoff,
     .   rt, nkx2, nky2, xkprl, xkperp, xkalp, xkbet,
     .   z0_table, z1_table, z2_table, zetai_table, 
     .   dKdL_table, 
     .   dKdL_giv, nmax, mmax, use_new_z2, ntable, 
     .   mtable)     


*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------
c      use zfun_hilbert  
      
      implicit none
      
      integer :: ieer ! fft error flag     
      
      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi, ndist, myid, iflag, nproc, ni, mi,
     .    i_sav, j_sav, ni0, mi0, upshift, nkx2, nky2
      integer n_upper, n_lower, m_upper, m_lower, nmax, mmax          
      integer ntable, mtable     
      
      complex z0_table(ntable, mtable)
      complex z1_table(ntable, mtable)
      complex z2_table(ntable, mtable)            
      real zetai_table(ntable), dKdL_table(mtable), dKdL_giv
           
      real u0, u2, fnorm, f_cql, rt
      real fnmax, fmmax, nperp
      real fna, fma, fstepn, fstepm, fstep

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, y, sgn_kprl, reson, duperp, dupara
      real xkprl_eff0
      real dzetal(lmin:lmax), descrim
      real dakbdkb, xnuomg, gradprlb, bmod, bmod0, nu_coll
      real akprl,  rho, alpha, eps0, omgrf, v0i, emax, akprl_min
      real gammab(lmin:lmax), gamma_coll(lmin:lmax)
      real a, b, xnurf, pi, delta0, rhol
      real bx, by, bz, bratio, denom
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt

      real xkxsav, xkysav, capr, argd
      real xkphi
      real xkalp, xkbet, xk0, rgamma, xk_cutoff, damping, kr, step

      complex zi, zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2
      complex z0_new, z1_new, z2_new

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0_a, sig1_a, sig2_a, sig3_a, sig4_a, sig5_a
      complex sig0_h, sig1_h, sig2_h, sig3_h, sig4_h, sig5_h

      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
      logical :: l_interp   !new
      logical :: l_first  !new
      logical :: use_new_z2
      integer :: nkperp
      real :: kperp_max !new


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

      complex zetal(lmin:lmax)
      complex  zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      integer  :: n_psi_dim
      integer :: nuper, nupar, n_psi

      integer NBESSJ

      real :: UPERP(NUPER),UPARA(NUPAR)
      real :: UPERP0, UPARA0

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
      
      common/upcom/akprl_min     
      
      allocate( dfduper0(nuper, nupar) )
      allocate( dfdupar0(nuper, nupar) )

      nu_coll =  .01 * omgrf
      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rhol = alpha / omgc
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * xnuomg)

      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl) 

c      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
c      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
c      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi
c      xkperp = sqrt(xkalp**2 + xkbet**2)      

      
*     ------------------------------------
*     Optional: leave out upshift in xkprl
*     --------------------------------- --          
c      if (upshift .eq. 0)  xkprl = uzz * xkphi
      
c      if (upshift .eq. -1) then      
c         if (xkperp  .gt. xk_cutoff) xkprl = uzz * xkphi
c      end if
      
c      if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
c      if (xkperp .eq. 0.0) xkperp = 1.0e-08
                        
                      
                     
!     ----------------------------------------------
!     Optional: Don't allow xkprl to be 0 (upshift = -2)
!     ----------------------------------------------        
c      if (upshift .eq. -2) then
c         if (akprl .lt. akprl_min) then
c            xkprl = akprl_min* sgn_kprl
c         end if 
c      end if   
            
      
c      if(xkperp .gt. kperp_max)then
c         write (6, *)"xkperp is gt kperp_max"
c         write (15, *)"xkperp is gt kperp_max"
c      end if
            
      
*     ---------------------------------
*     Calculate zetal(l) and gammab(l)
*     ---------------------------------      
      
      do l = lmin, lmax
         labs = abs(l)

         reson = (omgrf - l * real(omgc)) / omgrf
c         if (abs(reson) .lt. 0.02)then
         if (rho .gt. 1.0) then
            zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
            dzetal(l) = omgrf * xnuomg / (xkprl * alpha)
         else
            zetal(l) = (omgrf  - l * omgc) / (xkprl * alpha)
            dzetal(l) = 0.0
         end if
         
c        zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
c         dzetal(l) = omgrf * xnuomg / (xkprl * alpha)


         gammab(l) = abs(l * omgc / (2.0 * alpha * xkprl**2)
     .                                           * gradprlb / bmod)
         gamma_coll(l) = nu_coll / (akprl * alpha)


c         if(xm .eq. xme)gammab(l) = 0.0
c         if(abs(gammab(l)) .gt. 1000.0) gammab(l) = 1000.0
         if(abs(gammab(l)) .lt. .01)gammab(l) = .01


      enddo

      
      
*     ------------------------------------------------
*     Calculate Brambilla's xkrpl_eff using l = 1 only
*     ------------------------------------------------
      y0 = 1.5
      y = y0
      

      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gammab(1) * y) - 1.)
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            descrim = 1. - 4. * gammab(1) * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gammab(1) * y) )
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if
                    

*     -----------------------
*     Maxwellian distribution
*     -----------------------

      if(ndist .eq. 0)then

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

           if(nzfun .eq. 0) call z_approx(sgn_kprl, zetal(l), 0.0,
     .                                                     z0, z1, z2)
     
           if(nzfun .eq. 1) then
              call z_approx(sgn_kprl, zetal(l), gammab(l), z0, z1, z2)
     
c              if(use_new_z2 .eq. .true. ) then       
c                 call z_approx_i(sgn_kprl, zetal(l), gammab(l),
c     .              z0_new, z1_new, z2_new, zetai_table, 
c     .              dKdL_table, z0_table, z1_table, z2_table, 
c     .              dKdL_giv, nmax, mmax, ntable, mtable)
c                 z0 = z0_new
c                 z1 = z1_new
c                 z2 = z2_new            
c              end if

           end if                
                 

           if(nzfun .eq. 2) call z_smithe(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 3) call z_table(sgn_kprl,zetal(l),gammab(l),
     .                                      gamma_coll(l), z0, z1, z2)


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

      end if

*     -----------------------------------
*     Funky METS calls for non Maxwellian:
*     -----------------------------------
      if (ndist .eq. 1) then
      
         if (upshift .ne. 0) xkprl = xkprl_eff
         
      
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
         DensSPEC = xn
         bratio = bmod0 / bmod
         if(bratio .gt. 1.0) bratio = 1.0

         duperp = uperp(nuper) / (nuper - 1)
         dupara = 2.0 * upara(nupar) / (nupar - 1)

         if(i .ne. i_sav .or. j .ne. j_sav)then
         
            dfduper0 = 0.0
            dfdupar0 = 0.0

!           ------------------------------------------------
!           get CQL3D distribution function on the midplane
!           ------------------------------------------------
            call cql3d_dist(nupar, nuper, n_psi,
     .                 n_psi_dim, rho_a, rho,
     .                 UminPara,UmaxPara,
     .                 df_cql_uprp, df_cql_uprl,
     .                 UPERP, UPARA, DFDUPER0, DFDUPAR0)
     

!           ------------------------------------------------
!           map CQL3D distribution function off the midplane
!           ------------------------------------------------

            if(bratio .ge. 0.0)then
            
               dfduper = 0.0
               dfdupar = 0.0
            
               do ni = 1, nuper
                  do mi = 1, nupar

                     argd = uperp(ni)**2 * (1. - bratio) 
     .                                               + upara(mi)**2
                     if (argd .le. 0.0) argd = 1.0e-06
                     
                     uperp0 = uperp(ni) * sqrt(bratio)
                     upara0 = sign(1.0, upara(mi)) * sqrt(argd)
                     
                     dfduper(ni, mi) = 0.0
                     dfdupar(ni, mi) = 0.0
                     
                     if(upara0 .ge. upara(1) .and. 
     .                                   upara0 .le. upara(nupar)) then
     
                        ni0 = int((uperp0 - uperp(1)) / duperp) + 1
                        mi0 = int((upara0 - upara(1)) / dupara) + 1
                        
                        dfduper0_intplt = dfduper0(ni0, mi0)
                        dfdupar0_intplt = dfdupar0(ni0, mi0)
                        
                        if (ni0 .lt. nuper .and. mi0 .lt. nupar) then
                        
                        uperp0_grid = uperp(1) + (ni0 - 1) * duperp
                        upara0_grid = upara(1) + (mi0 - 1) * dupara
                                                
                        zeta = (uperp0 - uperp0_grid) / duperp
                        eta  = (upara0 - upara0_grid) / dupara
                        
                        ai = dfduper0(ni0, mi0)
                        bi = dfduper0(ni0+1 ,mi0) - dfduper0(ni0, mi0)
                        ci = dfduper0(ni0, mi0+1) - dfduper0(ni0, mi0)
                        di = dfduper0(ni0+1, mi0+1)+ dfduper0(ni0, mi0) 
     .                     - dfduper0(ni0+1, mi0) - dfduper0(ni0, mi0+1) 
                           
                        dfduper0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta                         

                        ai = dfdupar0(ni0, mi0)
                        bi = dfdupar0(ni0+1 ,mi0) - dfdupar0(ni0, mi0)
                        ci = dfdupar0(ni0, mi0+1) - dfdupar0(ni0, mi0)
                        di = dfdupar0(ni0+1, mi0+1)+ dfdupar0(ni0, mi0) 
     .                     - dfdupar0(ni0+1, mi0) - dfdupar0(ni0, mi0+1) 
                           
                        dfdupar0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta
     
                        end if                  

                        
                        if (upara0 .ne. 0.0)then
                                             
                           dfdupar(ni, mi) = dfdupar0_intplt * 
     .                        upara(mi) / upara0
     
                           dfduper(ni, mi) = dfduper0_intplt * 
     .                        sqrt(bratio) + dfdupar0_intplt * 
     .                        uperp(ni) / upara0 * (1.0 - bratio)
     
c                          dfdth = upara(mi) * dfduper(ni, mi)
c     .                           - uperp(ni) * dfdupar(ni, mi)


                        end if
                        
                     end if
                     
                     
                     go to 5000
!                    ----------------------------
!                    optional analytic Maxwellian
!                    ----------------------------
                     pi = 3.141592654
                     
                     alpha = sqrt(2.0 * xkt / xm)
!                     vc_mks = 3.5 * alpha
                     u0 = vc_mks / alpha
                     
                     fnorm = u0**3 / pi**1.5 
                           
                     u2 = uperp(ni)**2 + upara(mi)**2
                     
                     f_cql = exp(-u2 * u0**2) * fnorm
                     dfduper(ni, mi) = -f_cql * 2. * uperp(ni) * u0**2
                     dfdupar(ni, mi) = -f_cql * 2. * upara(mi) * u0**2
 5000                continue                


                  end do
               end do
               
            end if
            

!           --------------------------------------
!           Initialize the interpolation in k_perp:
!           --------------------------------------
            if (nkperp .ne. 0) then
               l_first  = .true.
               l_interp = .true.
        
               call GETNONMAXSIGMA_AORSA_NEWi_1(W,
     .                          ZSPEC,ASPEC,DensSPEC,BMAG,
     .                          K1,XI1,JNXI1,
     .                          K1,XI1,JNXI1,NBESSJ,
     .                          Enorm,UminPara,UmaxPara,
     .                          NUPAR,NUPER,UPERP,UPARA,
     .                          DFDUPER,DFDUPAR,
     .                          WSPEC,IFAIL,
     .                          l_first, l_interp, kperp_max, nkperp,
     .                          xkphi)
            end if



            i_sav = i
            j_sav = j

         end if






!        ------------------------------------
!        Complete integrals; no interpolation
!        ------------------------------------
         if (nkperp .eq. 0)then

            call WMATPRECALC_AORSA(ZSPEC,ASPEC,ENORM,BMAG,KPER1,UPERP,
     &                   NUPER,NBESSJ,NSBESSJ,XI1,JNXI1,IFAIL)

            call GETNONMAX_SIGMA_AORSA_NEW(W,
     .                          ZSPEC, ASPEC, DensSPEC, BMAG,
     .                          K1, XI1, JNXI1,
     .                          K1, XI1, JNXI1, NBESSJ,
     .                          Enorm, UminPara, UmaxPara,
     .                          NUPAR, NUPER, UPERP, UPARA,
     .                          DFDUPER, DFDUPAR,
     .                          WSPEC, IFAIL)

         else
!        -----------------
!        Use interpolation
!        -----------------

            l_first = .false.

              call GETNONMAXSIGMA_AORSA_NEWi_1(W,
     .                       ZSPEC,ASPEC,DensSPEC,BMAG,
     .                       K1, XI1, JNXI1,
     .                       K1, XI1, JNXI1, NBESSJ,
     .                       Enorm, UminPara, UmaxPara,
     .                       NUPAR, NUPER, UPERP, UPARA,
     .                       DFDUPER, DFDUPAR,
     .                       WSPEC, IFAIL,
     .                       l_first, l_interp, kperp_max, nkperp,
     .                       xkphi)
         end if



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

      end if

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
      

      deallocate( dfduper0 )
      deallocate( dfdupar0 )


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
  102 format(2i10, 1p8e12.4)
  103 format(4i10, 1p8e12.4)
      end


c
c***************************************************************************
c


      subroutine sigmad_cql3d_2(i, j, n, m, rho, rho_a,
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
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .   nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .   i_sav, j_sav, upshift, damping, xk_cutoff, rt, nkx2, nky2,
     .   xkprl, xkperp, xkalp, xkbet)


*     ---------------------------------------------------------
*     This routine uses the modified Z functions Z0, Z1, Z2
*     with the appropriate sign changes for k_parallel < 0.0
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------
c      use zfun_hilbert  
      
      implicit none
      
      integer :: ieer ! fft error flag
      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi, ndist, myid, iflag, nproc, ni, mi,
     .    i_sav, j_sav, ni0, mi0, upshift, nkx2, nky2
      integer n_upper, n_lower, m_upper, m_lower     
     
     
      real u0, u2, fnorm, f_cql, rt
      real fnmax, fmmax, nperp
      real fna, fma, fstepn, fstepm, fstep

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, y, sgn_kprl, reson, duperp, dupara
      real xkprl_eff0
      real dzetal(lmin:lmax), descrim
      real dakbdkb, xnuomg, gradprlb, bmod, bmod0, nu_coll
      real akprl,  rho, alpha, eps0, omgrf, v0i, emax, akprl_min
      real gammab(lmin:lmax), gamma_coll(lmin:lmax)
      real a, b, xnurf, pi, delta0, rhol
      real bx, by, bz, bratio, denom
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt

      real xkxsav, xkysav, capr, argd
      real xkphi
      real xkalp, xkbet, xk0, rgamma, xk_cutoff, damping, kr, step

      complex zi, zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0_a, sig1_a, sig2_a, sig3_a, sig4_a, sig5_a
      complex sig0_h, sig1_h, sig2_h, sig3_h, sig4_h, sig5_h

      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
      logical :: l_interp   !new
      logical :: l_first  !new
      integer :: nkperp
      real :: kperp_max !new


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

      complex zetal(lmin:lmax)
      complex  zieps0, arg,
     .   al, bl, cl,
     .   gamma, zeta_eff

      integer  :: n_psi_dim
      integer :: nuper, nupar, n_psi

      integer NBESSJ

      real :: UPERP(NUPER),UPARA(NUPAR)
      real :: UPERP0, UPARA0

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
      
      common/upcom/akprl_min     
      
      allocate( dfduper0(nuper, nupar) )
      allocate( dfdupar0(nuper, nupar) )

      nu_coll =  .01 * omgrf
      xme = 9.11e-31
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rhol = alpha / omgc
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * xnuomg)

      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl)   


c      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
c      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
c      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi
c      xkperp = sqrt(xkalp**2 + xkbet**2)      
   
      
*     ------------------------------------
*     Optional: leave out upshift in xkprl
*     --------------------------------- --          
c      if (upshift .eq. 0)  xkprl = uzz * xkphi
      
c      if (upshift .eq. -1) then      
c         if (xkperp  .gt. xk_cutoff) xkprl = uzz * xkphi
c      end if
      
c      if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
c      if (xkperp .eq. 0.0) xkperp = 1.0e-08
                        
                      
                     
!     ----------------------------------------------
!     Optional: Don't allow xkprl to be 0 (upshift = -2)
!     ----------------------------------------------        
c      if (upshift .eq. -2) then
c         if (akprl .lt. akprl_min) then
c            xkprl = akprl_min* sgn_kprl
c         end if 
c      end if   
            
      
c      if(xkperp .gt. kperp_max)then
c         write (6, *)"xkperp is gt kperp_max"
c         write (15, *)"xkperp is gt kperp_max"
c      end if
            
      
*     ---------------------------------
*     Calculate zetal(l) and gammab(l)
*     ---------------------------------      
      
      do l = lmin, lmax
         labs = abs(l)

         reson = (omgrf - l * real(omgc)) / omgrf
c         if (abs(reson) .lt. 0.02)then
         if (rho .gt. 1.0) then
            zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
            dzetal(l) = omgrf * xnuomg / (xkprl * alpha)
         else
            zetal(l) = (omgrf  - l * omgc) / (xkprl * alpha)
            dzetal(l) = 0.0
         end if
         
c        zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)
c         dzetal(l) = omgrf * xnuomg / (xkprl * alpha)


         gammab(l) = abs(l * omgc / (2.0 * alpha * xkprl**2)
     .                                           * gradprlb / bmod)
         gamma_coll(l) = nu_coll / (akprl * alpha)


         if(xm .eq. xme)gammab(l) = 0.0
c         if(abs(gammab(l)) .gt. 1000.0) gammab(l) = 1000.0
         if(abs(gammab(l)) .lt. .01)gammab(l) = .01


      enddo

      
      
*     ------------------------------------------------
*     Calculate Brambilla's xkrpl_eff using l = 1 only
*     ------------------------------------------------
      y0 = 1.5
      y = y0
      

      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gammab(1) * y) - 1.)
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            descrim = 1. - 4. * gammab(1) * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gammab(1) * y) )
     .         / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam 

      end if
                    

*     -----------------------
*     Maxwellian distribution
*     -----------------------

      if(ndist .eq. 0)then

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

           if(nzfun .eq. 0) call z_approx(sgn_kprl, zetal(l), 0.0,
     .                                                     z0, z1, z2)
           if(nzfun .eq. 1) call z_approx(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 2) call z_smithe(sgn_kprl,zetal(l),gammab(l),
     .                                                     z0, z1, z2)
           if(nzfun .eq. 3) call z_table(sgn_kprl,zetal(l),gammab(l),
     .                                      gamma_coll(l), z0, z1, z2)


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

      end if

*     -----------------------------------
*     Funky METS calls for non Maxwellian:
*     -----------------------------------
      if (ndist .eq. 1) then
      
         if (upshift .ne. 0) xkprl = xkprl_eff
         
      
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
         DensSPEC = xn
         bratio = bmod0 / bmod
         if(bratio .gt. 1.0) bratio = 1.0

         duperp = uperp(nuper) / (nuper - 1)
         dupara = 2.0 * upara(nupar) / (nupar - 1)

         if(i .ne. i_sav .or. j .ne. j_sav)then
         
            dfduper0 = 0.0
            dfdupar0 = 0.0

!           ------------------------------------------------
!           get CQL3D distribution function on the midplane
!           ------------------------------------------------
            call cql3d_dist(nupar, nuper, n_psi,
     .                 n_psi_dim, rho_a, rho,
     .                 UminPara,UmaxPara,
     .                 df_cql_uprp, df_cql_uprl,
     .                 UPERP, UPARA, DFDUPER0, DFDUPAR0)
     

!           ------------------------------------------------
!           map CQL3D distribution function off the midplane
!           ------------------------------------------------

            if(bratio .ge. 0.0)then
            
               dfduper = 0.0
               dfdupar = 0.0
            
               do ni = 1, nuper
                  do mi = 1, nupar

                     argd = uperp(ni)**2 * (1. - bratio) 
     .                                               + upara(mi)**2
                     if (argd .le. 0.0) argd = 1.0e-06
                     
                     uperp0 = uperp(ni) * sqrt(bratio)
                     upara0 = sign(1.0, upara(mi)) * sqrt(argd)
                     
                     dfduper(ni, mi) = 0.0
                     dfdupar(ni, mi) = 0.0
                     
                     if(upara0 .ge. upara(1) .and. 
     .                                   upara0 .le. upara(nupar)) then
     
                        ni0 = int((uperp0 - uperp(1)) / duperp) + 1
                        mi0 = int((upara0 - upara(1)) / dupara) + 1
                        
                        dfduper0_intplt = dfduper0(ni0, mi0)
                        dfdupar0_intplt = dfdupar0(ni0, mi0)
                        
                        if (ni0 .lt. nuper .and. mi0 .lt. nupar) then
                        
                        uperp0_grid = uperp(1) + (ni0 - 1) * duperp
                        upara0_grid = upara(1) + (mi0 - 1) * dupara
                                                
                        zeta = (uperp0 - uperp0_grid) / duperp
                        eta  = (upara0 - upara0_grid) / dupara
                        
                        ai = dfduper0(ni0, mi0)
                        bi = dfduper0(ni0+1 ,mi0) - dfduper0(ni0, mi0)
                        ci = dfduper0(ni0, mi0+1) - dfduper0(ni0, mi0)
                        di = dfduper0(ni0+1, mi0+1)+ dfduper0(ni0, mi0) 
     .                     - dfduper0(ni0+1, mi0) - dfduper0(ni0, mi0+1) 
                           
                        dfduper0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta                         

                        ai = dfdupar0(ni0, mi0)
                        bi = dfdupar0(ni0+1 ,mi0) - dfdupar0(ni0, mi0)
                        ci = dfdupar0(ni0, mi0+1) - dfdupar0(ni0, mi0)
                        di = dfdupar0(ni0+1, mi0+1)+ dfdupar0(ni0, mi0) 
     .                     - dfdupar0(ni0+1, mi0) - dfdupar0(ni0, mi0+1) 
                           
                        dfdupar0_intplt = ai + bi * zeta 
     .                                     + ci * eta + di * zeta * eta
     
                        end if                  

                        
                        if (upara0 .ne. 0.0)then
                                             
                           dfdupar(ni, mi) = dfdupar0_intplt * 
     .                        upara(mi) / upara0
     
                           dfduper(ni, mi) = dfduper0_intplt * 
     .                        sqrt(bratio) + dfdupar0_intplt * 
     .                        uperp(ni) / upara0 * (1.0 - bratio)
     
c                          dfdth = upara(mi) * dfduper(ni, mi)
c     .                           - uperp(ni) * dfdupar(ni, mi)


                        end if
                        
                     end if
                     
                     
                     go to 5000
!                    ----------------------------
!                    optional analytic Maxwellian
!                    ----------------------------
                     pi = 3.141592654
                     
                     alpha = sqrt(2.0 * xkt / xm)
!                     vc_mks = 3.5 * alpha
                     u0 = vc_mks / alpha
                     
                     fnorm = u0**3 / pi**1.5 
                           
                     u2 = uperp(ni)**2 + upara(mi)**2
                     
                     f_cql = exp(-u2 * u0**2) * fnorm
                     dfduper(ni, mi) = -f_cql * 2. * uperp(ni) * u0**2
                     dfdupar(ni, mi) = -f_cql * 2. * upara(mi) * u0**2
 5000                continue                


                  end do
               end do
               
            end if
            

!           --------------------------------------
!           Initialize the interpolation in k_perp:
!           --------------------------------------
            if (nkperp .ne. 0) then
               l_first  = .true.
               l_interp = .true.
        
               call GETNONMAXSIGMA_AORSA_NEWi_2(W,
     .                          ZSPEC,ASPEC,DensSPEC,BMAG,
     .                          K1,XI1,JNXI1,
     .                          K1,XI1,JNXI1,NBESSJ,
     .                          Enorm,UminPara,UmaxPara,
     .                          NUPAR,NUPER,UPERP,UPARA,
     .                          DFDUPER,DFDUPAR,
     .                          WSPEC,IFAIL,
     .                          l_first, l_interp, kperp_max, nkperp,
     .                          xkphi)
            end if



            i_sav = i
            j_sav = j

         end if






!        ------------------------------------
!        Complete integrals; no interpolation
!        ------------------------------------
         if (nkperp .eq. 0)then

            call WMATPRECALC_AORSA(ZSPEC,ASPEC,ENORM,BMAG,KPER1,UPERP,
     &                   NUPER,NBESSJ,NSBESSJ,XI1,JNXI1,IFAIL)

            call GETNONMAX_SIGMA_AORSA_NEW(W,
     .                          ZSPEC, ASPEC, DensSPEC, BMAG,
     .                          K1, XI1, JNXI1,
     .                          K1, XI1, JNXI1, NBESSJ,
     .                          Enorm, UminPara, UmaxPara,
     .                          NUPAR, NUPER, UPERP, UPARA,
     .                          DFDUPER, DFDUPAR,
     .                          WSPEC, IFAIL)

         else
!        -----------------
!        Use interpolation
!        -----------------

            l_first = .false.

              call GETNONMAXSIGMA_AORSA_NEWi_2(W,
     .                       ZSPEC,ASPEC,DensSPEC,BMAG,
     .                       K1, XI1, JNXI1,
     .                       K1, XI1, JNXI1, NBESSJ,
     .                       Enorm, UminPara, UmaxPara,
     .                       NUPAR, NUPER, UPERP, UPARA,
     .                       DFDUPER, DFDUPAR,
     .                       WSPEC, IFAIL,
     .                       l_first, l_interp, kperp_max, nkperp,
     .                       xkphi)
         end if



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

      end if

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
      

      deallocate( dfduper0 )
      deallocate( dfdupar0 )


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
  102 format(2i10, 1p8e12.4)
  103 format(4i10, 1p8e12.4)
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
      integer   :: nupar, nuper, n_psi
      integer   :: n_psi_dim

      real   :: UminPara,UmaxPara
      real   :: UPERP(NUPER)
      real   :: UPARA(NUPAR)
      real   :: DFDUPER(NUPER,NUPAR),
     .                      DFDUPAR(NUPER,NUPAR)

      real   :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real   :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: f, u2, rho, rho_a(n_psi_dim)

      integer n, m, i_psi, i_psin

      i_psi = 1


      do i_psin = 1, n_psi - 1
         if(rho .ge. rho_a(i_psin)     .and.
     .      rho .lt. rho_a(i_psin + 1))
     .      i_psi = i_psin
      end do


      if(rho .ge. rho_a(n_psi)) i_psi = n_psi


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

      subroutine midplane(igiven, jgiven, f, fmid, rho,
     .   nxdim, nydim, nnodex, nnodey, capr, r0, f0, jmid)

      implicit none

      integer nxdim, nydim, nnodex, nnodey, i, j, jmid, i0
      integer igiven, jgiven, istart

      real f(nxdim, nydim), fmid, rho(nxdim, nydim), r0, capr(nxdim)
      real rhoij, f0

      rhoij = rho(igiven, jgiven)

      fmid =  0.0
      if (rhoij .gt. 1.0) return

c      jmid  = nnodey / 2

c      write(6, *)"jmid = ", jmid

c     DLG: find the first R coord past the axis R=r0
      do i = 1, nnodex
         if(capr(i) .ge. r0)then
            istart = i
            go to 200
         end if
      end do

      stop "ERROR (sigma.f): no capr(i) >= r0?"

  200 continue

c     DLG: for the z slice on the magnetic axis (jmid)
      if(rhoij .ge. rho(istart, jmid)) then

         if(rhoij .lt. rho(istart,jmid)) i0 = istart
         if(rhoij .ge. rho(nnodex,jmid)) then
            write(*,*) 'WARNING (sigma.f) : rhoij > rho(nnodex,jmid)'
            i0 = nnodex-1
         endif

         do i = istart, nnodex - 1
            if(rhoij .ge. rho(i,   jmid) .and.
     .         rhoij .lt. rho(i+1, jmid)) i0 = i
         end do
         fmid = f(i0, jmid) + (f(i0 + 1, jmid) -   f(i0, jmid))
     .                    / (rho(i0 + 1, jmid) - rho(i0, jmid))
     .                    * (rhoij             - rho(i0, jmid))
      end if


      if(rhoij .lt. rho(istart, jmid)) then
         fmid = f0
     .    + (f(istart, jmid) - f0) / (rho(istart, jmid) - 0.0)
     .             * (rhoij - 0.0)
      end if



  100 format (1i10, 1p8e12.4)
  102 format (2i10)
      return
      end

c
c***************************************************************************
c
       subroutine dummy_dist(NUPAR, NUPER,
     .                       UminPara, UmaxPara,
     &                       UPERP, UPARA, DFDUPER, DFDUPAR)
       implicit none

       integer  ::NUPAR,NUPER
       real  :: UminPara,UmaxPara
       real  :: UPERP(NUPER)
       real  :: UPARA(NUPAR)
       real  :: DFDUPER(NUPER,NUPAR),
     &                      DFDUPAR(NUPER,NUPAR)
       real :: pi, fnorm, f, u2
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


       subroutine dummy_dist2(NUPAR,NUPER,UminPara,UmaxPara,
     &                       UPERP,UPARA,DFDUPER,DFDUPAR, fnorm)
       implicit none
       integer  :: NUPAR,NUPER
       real  :: UminPara,UmaxPara
       real  :: UPERP(NUPER)
       real  :: UPARA(NUPAR)
       real  :: DFDUPER(NUPER,NUPAR),
     &                      DFDUPAR(NUPER,NUPAR)
       real :: pi, fnorm, f, u, u2, u3
       integer n,m


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
             u  = sqrt(u2)
             u3 = u**3
             f = fnorm * 1.0 / (1. + u3)
             DFDUPER(n, m) = -3.0 * f / (1. + u3) * UPERP(n) * u
             DFDUPAR(n, m) = -3.0 * f / (1. + u3) * UPARA(n) * u
          end do
       end do

       return
       end subroutine dummy_dist2



c
c***************************************************************************
c

      subroutine sigmah_slow(i, j, n, m,
     .   xm, q, xn, xnuomg,
     .   eslow, omgc, omgp2,
     .   lmin, lmax, nzfun, ibessel,
     .   xkxsav, xkysav, nphi, capr,
     .   bx,    by,    bz,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz,
     .   xkte, zeff, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .   i_sav, j_sav)



*     ---------------------------------------------------------
*     This routine uses call sigslo for the Hermitian part of
*     dispersion tensor for  an arbitrary isotropic
*     distribution function of a single plasma species
*     No rotation is made.  Result is in the Stix frame.
*     ---------------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi,
     .    i_sav, j_sav

      real xkperp, xkprl, xm, q, xn, eslow, omgc, omgp2, xme
      real xkprl_eff, fgam, y0, sgn_kprl
      real dakbdkb, xnuomg
      real akprl, gammab, valpha, eps0, omgrf, v0i, akprl_min
      real a, b, xkte, velect, zeff
      real bx, by, bz

      real xkxsav, xkysav, capr
      real xkphi
      real xkalp, xkbet, xk0

      complex zi, omgrfc
      complex sig0, sig1, sig2, sig3, sig4, sig5


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz

      real uxx, uxy, uxz,
     1     uyx, uyy, uyz,
     2     uzx, uzy, uzz

      real :: kperp_max !new

      common/upcom/akprl_min
      
      xme = 9.11e-31


      valpha = sqrt(2. * eslow / xm)
      velect = sqrt(2.0 * xkte / xme)
      if (velect .eq. 0.0) return

      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * xnuomg)

      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi

      xkperp = sqrt(xkalp**2 + xkbet**2)

      akprl = abs(xkprl)
      sgn_kprl = 1.0
      if(akprl .ne. 0.0)sgn_kprl = xkprl / akprl
    


      call sigslo(q, zeff, velect, valpha, omgc, omgp2, omgrf,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz,
     .   xkprl, xkperp, lmin, lmax)

*     -----------------------------
*     Change to Swanson's notation
*     -----------------------------
      sig1 = sigxx
      sig2 = sigxy
      sig3 = sigzz
      sig4 = sigxz / xkperp
      sig5 = - sigyz / xkperp
      sig0 = (sigyy - sig1) / xkperp**2


*     -----------------------------
*     Swanson's rotation (original)
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
     .   delta0, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .   i_sav, j_sav)



*     ----------------------------------------------------
*     This routine calculates sigma_cold in the Stix frame
*     ----------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel,
     .    i, j, n, m, nphi,
     .    i_sav, j_sav

      real xkperp, xkprl, xm, q, xn, xkt, omgc, omgp2, xkb, akb
      real dakbdkb
      real akprl, gammab, alpha, eps0, omgrf, v0i, akprl_min
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

      real :: kperp_max !new

      common/upcom/akprl_min


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
      real factrl, factl, infinity_check
      complex gamma, expbes(0:lmaxdim), expbesp(0:lmaxdim),
     .   expbesovergam(0:lmaxdim), xilovergam,
     .   xil(0:lmaxdim), xilp(0:lmaxdim), exgam

      complex b(100)

      exgam = 1.0 - gamma + gamma**2 / 2.0

      do l = 0, lmax
         factl = factrl(l)

         infinity_check = 1/(2**l * factl)
         if ( infinity_check .gt. 1e30 ) then

          xil(l) = 0;
          xilp(l) = 0;
          xilovergam = 0;

         else

         xil(l) = gamma**l / (2**l * factl) *
     .                                 ( 1. + gamma**2 / (4. * (l+1)))
         xilp(l) = gamma**(l-1) / (2**l * factl) *
     .                     (l + (l+2) * gamma**2 / (4. * (l+1)))

         xilovergam = gamma**(l-1) / (2**l * factl) *
     .                                 ( 1. + gamma**2 / (4. * (l+1)))

         endif 

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
c***********************************************************************
c



