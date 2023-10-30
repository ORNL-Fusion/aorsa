        module qlsum_myra_mod
         contains
!module METS2AORSA_MYRA

!  private
!  public::  GET_WMAT_MYRA, intplt1d, WINTERP1D_3

!contains


!
!*************************************************************************
!
    subroutine QLSUM_MAXWELLIAN(k_uper, b_sum, c_sum, e_sum, f_sum, &
       & sum_wdot, sum_fx0, sum_fy0, W, ZSPEC, ASPEC, BMAG, &
       & lmax, ENORM, UPARMIN, UPARMAX, &
       & NUPAR, NUPER, UPER, UPAR, DFDUPER, DFDUPAR,  &
       & ealphak, ebetak, ebk, nkdim1, nkdim2, mkdim1, mkdim2,   &
       & nkx1, nkx2, nky1, nky2, &
       & uxx, uxy, uxz, &
       & uyx, uyy, uyz, &
       & uzx, uzy, uzz, &
       & nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i_global, j_global, &
       & lmaxdim, ndist, nzeta, &
       & gradprlb, bmod, omgc, alpha, xm, upshift, xk_cutoff, rt, nphi, rho)

    implicit none

    integer i_global, j_global, ier, nxdim, nydim, k, lmaxdim, ndist
    integer, intent(IN):: NUPAR, NUPER, lmax
    integer nkx1, nkx2, nky1, nky2, j_upar, k_uper, l
    integer nkdim1, nkdim2, mkdim1, mkdim2
    integer:: NHARM, IHARM, M, N, i, nzeta
    integer i_uprl, upshift
    integer ires, iresmax, nphi

    complex, dimension(:,:), allocatable :: zbeta
    complex, dimension(:,:), allocatable :: zbeta_iharm
    complex ttmp_11, cross_13, cross_31, landau_33, ttmp, cross2, cross3, ld
    
    real  y, y0, alpha, xm, akprl, sgn_kprl, omgc, bmod, xkprl_eff,  &
   &   descrim, xme, fgam, gradprlb, gammab, xk_cutoff, rt, akprl_min, rho

    real  uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz
    real  xkphi, sinth, factc, facte, factf, sinth_inv
    real  xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
    real  xkperpn, xkperpni, xkrhon, xketan, xkprln, beta
    real, intent(IN):: W, ZSPEC, ASPEC, BMAG
    real, intent(IN):: ENORM, UPARMIN, UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER, DFDUPAR
    real:: W2, WCW, RRP, RRM, WC, WCI
    real:: MUT0, SQMUT0, PISQMUT0, SQMUT0I, KPARA1
    real:: ISQ2, SQ2, NWCW, DFACTPAR, DFACTPER, U0
    real:: UPAR0, dfdupar0, dfduper0, du, dui, p 
    real:: time, t1, tmsec, second1, dummy, WI, uperpk, uperpk2
    real:: dzeta, dzetai, zetamax, zetamin, zeta0, zetamax1, zetamax2
    real:: A1, A2, A3, u
    real:: temp1, temp2, temp3, factor
    real:: temp1w, temp2w, temp3w
    
    real, dimension(:),     allocatable :: zetai
    real, dimension(:,:),   allocatable :: Jni
    real, dimension(:,:,:), allocatable :: Jn
    real, dimension(:,:),   allocatable :: NPARA_sav
    
    integer, dimension(:),  allocatable :: nres
    integer, dimension(:),  allocatable :: mres
        
    complex, dimension(:),  allocatable :: sumb_11
    complex, dimension(:),  allocatable :: sumb_31
     
    complex, dimension(:),  allocatable :: sumc_11
    complex, dimension(:),  allocatable :: sumc_31
    
    complex, dimension(:),  allocatable :: sume_11
    complex, dimension(:),  allocatable :: sume_31
     
    complex, dimension(:),  allocatable :: sumf_11
    complex, dimension(:),  allocatable :: sumf_31
    
    complex sumf_11_nm
    complex sumf_31_nm

    complex sume_11_nm
    complex sume_31_nm
       
    complex sumc_11_nm
    complex sumc_31_nm
       
    complex sumb_11_nm
    complex sumb_31_nm
 
    logical, dimension(:,:), allocatable :: is_resonance_nm      

    complex epsx, epsy, epsz
    complex ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
    &        ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
    &           ebk(nkdim1 : nkdim2, mkdim1 : mkdim2)
    complex cexpkx, cexpky, zi, zeta

    complex xx(nkdim1 : nkdim2, 1 : nxdim),   &
     &      yy(mkdim1 : mkdim2, 1 : nydim)

    complex cexpn, cexpnp1, cexpnm1, cexp11

    complex cexp1, cexp2, cexp0
    complex sumwdot_11_nm, sumwdot_31_nm
    
    complex sumwdot_11, sumwdot_31
    complex sumwdotkx_11, sumwdotkx_31
    complex sumwdotky_11, sumwdotky_31
        
    complex sum2_1, sum2_2, sum2_3
    complex sum1_1, sum1_2, sum1_3, sum1_4
    complex sumkx2_1, sumkx2_2, sumkx2_3 
    complex sumky2_1, sumky2_2, sumky2_3        
    
    complex b_sum(nupar), c_sum(nupar), e_sum(nupar), f_sum(nupar)
    complex sum_wdot, sum_fx0, sum_fy0
    complex b(100)

    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
    
    real, parameter:: MPC2 = 938271998.38
    real, parameter:: C = 2.99792458e8
    real, parameter:: PI = 3.141592653597932384
    real :: cosbeta_n_m, sinbeta_n_m
    common/upcom/akprl_min    
    
    allocate( zbeta(nkx1:nkx2,nky1:nky2) )
    allocate( zbeta_iharm(nkx1:nkx2,nky1:nky2) )

    allocate(zetai(nzeta + 1) )
    allocate(Jni(-lmaxdim : lmaxdim, nzeta + 1) )
    allocate( Jn(-lmaxdim : lmaxdim, nkdim1 : nkdim2, mkdim1 : mkdim2))    
    allocate(NPARA_sav(nkdim1 : nkdim2, mkdim1 : mkdim2) ) 
    
    allocate(nres(nxdim * nydim) )
    allocate(mres(nxdim * nydim) )
    
    allocate(sumb_11(nupar))
    allocate(sumb_31(nupar))
     
    allocate(sumc_11(nupar))
    allocate(sumc_31(nupar))
    
    allocate(sume_11(nupar))
    allocate(sume_31(nupar))
     
    allocate(sumf_11(nupar))
    allocate(sumf_31(nupar))

    allocate(is_resonance_nm(nkx1:nkx2,nky1:nky2))

!   -------------------------------------
!   initialize allocatable arrays to zero
!   -------------------------------------

    zbeta = 0.0
    zbeta_iharm = 0.0

    sumf_11_nm = 0.0
    sumf_31_nm = 0.0

    sume_11_nm = 0.0
    sume_31_nm = 0.0
    
    sumc_11_nm = 0.0
    sumc_31_nm = 0.0
    
    sumb_11_nm = 0.0
    sumb_31_nm = 0.0

    is_resonance_nm = .false.

    zetai = 0.0
    Jni = 0.0
    NPARA_sav = 0.0

    sumb_11 = 0.0
    sumb_31 = 0.0
      
    sumc_11 = 0.0
    sumc_31 = 0.0

    sume_11 = 0.0
    sume_31 = 0.0

    sumf_11 = 0.0
    sumf_31 = 0.0    
    
    xme = 9.11e-31
    zi = cmplx(0., 1.)
    
    uperpk = uper(k_uper)
    uperpk2 = uperpk**2

    W2 = W * W
    WI = 1.0 / W

    WCW = BMAG * ZSPEC * EOVERMH / ASPEC / W
    WC = WCW * W
    WCI = 1.0 /WC

    MUT0 = 0.5 * MPC2 * ASPEC / ENORM
    SQMUT0 = SQRT(MUT0)
    SQMUT0I = 1.0 / SQMUT0
    PISQMUT0 = SQMUT0 * pi

    ISQ2 = SQRT(0.5)
    SQ2 = SQRT(2.0)
    NHARM = lmax
    
    du = (upar(nupar) - upar(1)) / (nupar - 1)
    dui = 1.0 / du
    
    if(nzeta .eq. 1)then
    
    ! -------------------------------------------------------- !
    ! ---Don't interpolate: precalculate all Bessel functions- !
    ! -------------------------------------------------------- !

       do n = nkx1, nkx2
          do m = nky1, nky2

             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             xkprln   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi
             xkperpn = sqrt(xkrhon**2 + xketan**2)
             
            ! ------------------------------------
            ! Optional: leave out upshift in xkprl
            ! --------------------------------- --          
              if (upshift .eq. 0)xkprln = uzz * xkphi
!             if (upshift .eq. 0)xkprln = nphi /rt
              
             if (upshift .eq. -1) then      
                if (xkperpn  .gt. xk_cutoff) xkprln = uzz * xkphi
             end if           
             
             
             sgn_kprl = sign(1.0, xkprln)
             akprl = abs(xkprln)                      
                     
!            ----------------------------------------------
!            Optional: Don't allow xkprl to be 0 (upshift = -2)
!            ----------------------------------------------        
             if (upshift .eq. -2) then
                if (akprl .lt. akprl_min) then
                    xkprln = akprl_min * sgn_kprl
                end if 
             end if                  
             

             y0 = 1.5
             y = y0
             
      
             l = 1
             if(xkprln .eq. 0)xkprln = 1.0e-06
             gammab = abs(l * omgc / (2.0 * alpha * xkprln**2)  &
     &                                    * gradprlb / bmod)

             if(xm .eq. xme)gammab = 0.0
!             if(abs(gammab) .gt. 1000.0) gammab = 1000.0
             if(abs(gammab) .lt. .01)gammab = .01


             if(sgn_kprl .ge. 0.0)then
                fgam = 1.0

                if(gammab .gt. 1.0e-05)then
                   y = y0
                   fgam = (sqrt(1. +  4. * gammab * y) - 1.)   &
     &               / (2. * gammab * y)
                endif

                xkprl_eff = xkprln / fgam 

             end if


             if(sgn_kprl .lt. 0.0)then
                fgam = 1.0

                if(gammab .gt. 1.0e-05)then
                   descrim = 1. - 4. * gammab * y0
                   if (descrim .ge. 0.0) y =   y0
                   if (descrim .lt. 0.0) y = - y0
                   fgam = (1. - sqrt(1. -  4. * gammab * y) )  &
     &                / (2. * gammab * y)
                endif

                xkprl_eff = xkprln / fgam 

             end if

             if (upshift .ne. 0) xkprln = xkprl_eff          
        
             NPARA_sav(n, m) = xkprln * C / W

             xkperpn = sqrt(xkrhon**2 + xketan**2)
             if(xkperpn .eq. 0.0)xkperpn = 1.0e-08
        
             cosbeta_n_m  = xkrhon / xkperpn
             sinbeta_n_m  = xketan / xkperpn
             zbeta(n,m) = cmplx( cosbeta_n_m , sinbeta_n_m  )

             zeta = xkperpn * uper(k_uper) * c * sqmut0i / wc

             call besjc(zeta, nharm + 2, b, ier)
             if(ier .ne. 0) write(6, *) "ier = ", ier

             do IHARM = 0, NHARM + 1
                Jn(iharm,  n, m) = b(iharm + 1)
                Jn(-iharm, n, m) = (-1.0)**iharm * Jn(iharm, n, m)
             end do

          end do
       end do
    
    else
    
        
       ! -------------------------------------- !
       ! ---Interpolate; calculate zeta mesh -- !
       ! -------------------------------------- !
    
       zetamax = 0.0
       zetamin = 0.0
    
       do n = nkx1, nkx2
          do m = nky1, nky2

             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             xkprln   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi

             xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08
        
             zeta0 = xkperpn * uperpk * c * sqmut0i * wci
          
             if (zeta0 .gt. zetamax) zetamax = zeta0
             if (zeta0 .lt. zetamin) zetamin = zeta0
          end do
       end do
    

        if(zetamax .eq. zetamin)then
          zetamax =  1.0e-06
          zetamin = -1.0e-06
       end if

       dzeta = (zetamax - zetamin) / (nzeta - 1)
       dzetai = 1.0 / dzeta
    
       ! ------------------------------------------------- !
       ! ---Pre-calculate Bessel functions on zeta mesh -- !
       ! ------------------------------------------------- !
          
       do i = 1, nzeta + 1
          zetai(i) = zetamin + (i - 1) * dzeta
          zeta = cmplx(zetai(i), 0.0)
          
          call besjc(zeta, nharm + 2, b, ier)
!          if(ier .ne. 0) write(6, *) "ier = ", ier
          
          do iharm = 0, NHARM + 1
             Jni(iharm,  i) = b(iharm + 1)
             Jni(-iharm, i) = (-1.0)**iharm * b(iharm + 1)
          end do
       end do
     

       ! --------------------------------- !
       ! ---Interpolate Bessel functions-- !
       ! --------------------------------- !

       do n = nkx1, nkx2
          do m = nky1, nky2

             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             xkprln   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi
             xkperpn = sqrt(xkrhon**2 + xketan**2)
             
            ! ------------------------------------
            ! Optional: leave out upshift in xkprl
            ! --------------------------------- --          
             if (upshift .eq. 0)xkprln = uzz * xkphi
!             if (upshift .eq. 0)xkprln = nphi / rt
              
             if (upshift .eq. -1) then      
                if (xkperpn  .gt. xk_cutoff) xkprln = uzz * xkphi
             end if
             
             
             sgn_kprl = sign(1.0, xkprln)
             akprl = abs(xkprln)                      
                     
!            ----------------------------------------------
!            Optional: Don't allow xkprl to be 0 (upshift = -2)
!            ----------------------------------------------        
             if (upshift .eq. -2) then
                if (akprl .lt. akprl_min) then
                    xkprln = akprl_min * sgn_kprl
                end if 
             end if                  
                     

             y0 = 1.5
             y = y0
             
      
             l = 1
             if(xkprln .eq. 0)xkprln = 1.0e-06
             gammab = abs(l * omgc / (2.0 * alpha * xkprln**2)  &
     &                                    * gradprlb / bmod)

             if(xm .eq. xme)gammab = 0.0
!             if(abs(gammab) .gt. 1000.0) gammab = 1000.0
             if(abs(gammab) .lt. .01)gammab = .01


             if(sgn_kprl .ge. 0.0)then
                fgam = 1.0

                if(gammab .gt. 1.0e-05)then
                   y = y0
                   fgam = (sqrt(1. +  4. * gammab * y) - 1.)   &
     &               / (2. * gammab * y)
                endif

                xkprl_eff = xkprln / fgam 

             end if


             if(sgn_kprl .lt. 0.0)then
                fgam = 1.0

                if(gammab .gt. 1.0e-05)then
                   descrim = 1. - 4. * gammab * y0
                   if (descrim .ge. 0.0) y =   y0
                   if (descrim .lt. 0.0) y = - y0
                   fgam = (1. - sqrt(1. -  4. * gammab * y) )  &
     &                / (2. * gammab * y)
                endif

                xkprl_eff = xkprln / fgam 

             end if

             if (upshift .ne. 0) xkprln = xkprl_eff
             
        
             NPARA_sav(n, m) = xkprln * C * WI

             xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08
             xkperpni = 1.0 / xkperpn
        
             cosbeta_n_m = xkrhon * xkperpni
             sinbeta_n_m  = xketan * xkperpni
             zbeta(n,m) = cmplx( cosbeta_n_m , sinbeta_n_m )
          
             zeta0 = xkperpn * uperpk * c * sqmut0i * wci
          
             i = int((zeta0 - zetamin) * dzetai) + 1
             p = (zeta0 - zetai(i)) * dzetai
             A1 = 0.5 * P * (P - 1.)
             A2 = 1. - P * P
             A3 = 0.5 * P * (P + 1.)
          
             do iharm = -NHARM - 1, NHARM + 1
          
                Jn(iharm, n, m) = Jni(iharm, i)    &
     &             + p * (Jni(iharm, i + 1) - Jni(iharm, i))
                if(i .ne. 1 )then
                   Jn(iharm, n, m) = A1 * Jni(iharm, i - 1)     &
     &                             + A2 * Jni(iharm, i)         &
     &                             + A3 * Jni(iharm, i + 1)
                end if
              
             end do
          
          end do
       end do
    

    end if
                    

    ! ------------------------ !
    ! ---Sum over harmonics--- !
    ! ------------------------ !
    
    sum_wdot = 0.0
    sum_fx0  = 0.0
    sum_fy0  = 0.0    
    
    b_sum = 0.0
    c_sum = 0.0
    e_sum = 0.0
    f_sum = 0.0
       
   do IHARM = -NHARM, NHARM   

       NWCW = real(IHARM) * WCW
       
       sumb_11 = 0.0
       sumb_31 = 0.0
          
       sumc_11 = 0.0
       sumc_31 = 0.0    
          
       sume_11 = 0.0
       sume_31 = 0.0
          
       sumf_11 = 0.0
       sumf_31 = 0.0      
                  
       sumwdot_11 = 0.0
       sumwdot_31 = 0.0 
       
       sum1_1 = 0.0
       sum1_2 = 0.0
       sum1_3 = 0.0
       sum1_4 = 0.0      
       
       sumwdotkx_11 = 0.0
       sumwdotkx_31 = 0.0
       
       sumwdotky_11 = 0.0
       sumwdotky_31 = 0.0       

       sum2_1 = 0.0
       sum2_2 = 0.0
       sum2_3 = 0.0
       
       sumkx2_1 = 0.0
       sumkx2_2 = 0.0
       sumkx2_3 = 0.0
       
       sumky2_1 = 0.0
       sumky2_2 = 0.0
       sumky2_3 = 0.0       
       
       call zpow((nkx2-nkx1+1)*(nky2-nky1+1), zbeta, iharm, zbeta_iharm) 
       
       ! ----------------------- !
       ! ---Find resonant modes--!
       ! ----------------------- !
       ires = 0
       do n = nkx1, nkx2
          do m = nky1, nky2     
             ! ------------------------ !
             ! -- Resonance relation -- !
             ! ------------------------ !
             RRP = 1.0 - NWCW - NPARA_sav(n, m) * UPARMAX * SQMUT0i
             RRM = 1.0 - NWCW - NPARA_sav(n, m) * UPARMIN * SQMUT0i

             is_resonance_nm(n,m) = (RRP * RRM .le. 0.0)
             if (is_resonance_nm(n,m)) then
                 ires = ires + 1
                 nres(ires) = n
                 mres(ires) = m
             end if
          end do
       end do
       iresmax = ires
       
       
!       go to 1000
      ! ------------------------------- !
      ! ---Sum over all Fourier modes-- !
      ! --------------------------------!
       do n = nkx1, nkx2
          do m = nky1, nky2
             cexp0 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m)
             cexp1 = cexp0 * zbeta(n,m)
             cexp2 = cexp0 / zbeta(n,m)
                
             epsx = isq2 * (ealphak(n, m) - zi * ebetak(n, m)) * cexp1
             epsy = isq2 * (ealphak(n, m) + zi * ebetak(n, m)) * cexp2
             epsz = ebk(n, m) * cexp0
             
             sum2_1 = sum2_1 + conjg(epsx) * Jn(IHARM + 1, n, m)
             sum2_2 = sum2_2 + conjg(epsy) * Jn(IHARM - 1, n, m)
             sum2_3 = sum2_3 + conjg(epsz) * Jn(IHARM, n, m)         
          end do  !  end sum over ky
       end do    !   end sum over kx
!  1000 continue

      ! ---------------------------- !
      ! ---Sum over resonant modes-- !
      ! ---------------------------- !

       do ires = 1, iresmax
          n = nres(ires)
          m = mres(ires)
                                                         
!          cexp1 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m) * zbeta(n,m)
!          cexp2 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m) / zbeta(n,m)
          cexp0 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m)
          cexp1 = cexp0 * zbeta(n,m)
          cexp2 = cexp0 / zbeta(n,m)
                
          epsx = isq2 * (ealphak(n, m) - zi * ebetak(n, m)) * cexp1
          epsy = isq2 * (ealphak(n, m) + zi * ebetak(n, m)) * cexp2
          epsz = ebk(n, m) * cexp0
             
!          sum2_1 = sum2_1 + conjg(epsx) * Jn(IHARM + 1, n, m)
!          sum2_2 = sum2_2 + conjg(epsy) * Jn(IHARM - 1, n, m)
!          sum2_3 = sum2_3 + conjg(epsz) * Jn(IHARM, n, m)
             
          sumkx2_1 = sumkx2_1 + xkxsav(n) * conjg(epsx) * Jn(IHARM + 1, n, m)
          sumkx2_2 = sumkx2_2 + xkxsav(n) * conjg(epsy) * Jn(IHARM - 1, n, m)
          sumkx2_3 = sumkx2_3 + xkxsav(n) * conjg(epsz) * Jn(IHARM, n, m)
             
          sumky2_1 = sumky2_1 + xkysav(m) * conjg(epsx) * Jn(IHARM + 1, n, m)
          sumky2_2 = sumky2_2 + xkysav(m) * conjg(epsy) * Jn(IHARM - 1, n, m)
          sumky2_3 = sumky2_3 + xkysav(m) * conjg(epsz) * Jn(IHARM, n, m)                           

          UPAR0 = SQMUT0 / NPARA_sav(n, m) * (1. - NWCW)
        
          u = sqrt(upar0**2 + uperpk2) + 1.0e-08
          sinth = uperpk / u + 1.0e-08
          sinth_inv = 1.0 / sinth
                
          facte = (nwcw - sinth**2) /  upar0
        
          i = int((UPAR0 - UPAR(1)) * dui) + 1
          i_uprl = i
          p = (UPAR0 - UPAR(i)) * dui
                
          dfduper0 = dfduper(k_uper, NUPAR)
          if (i .ne. NUPAR) then
             dfduper0 = dfduper(k_uper, i) + (dfduper(k_uper, i+1) - dfduper(k_uper, i)) * p
          end if
                
          U0 = DFDUPER0
                        
          factor = PISQMUT0 / abs(NPARA_sav(n, m)) 
          
          ttmp = UPER(k_uper) **2 * (Jn(IHARM + 1, n, m) * epsx     &
     &                            +  Jn(IHARM - 1, n, m) * epsy )
          cross2 = SQ2 * UPER(k_uper) * UPAR0 * Jn(IHARM, n, m) * epsz
          
          cross3 = SQ2 * UPER(k_uper) * UPAR0 * (Jn(IHARM + 1, n, m) * epsx     &
     &                                         + Jn(IHARM - 1, n, m) * epsy)
      
          ld = 2.0 * UPAR0**2 * Jn(IHARM, n, m) * epsz
          
          sumb_11_nm = ttmp + cross2
          sumb_31_nm = cross3 + ld        
                        
!         sumb_11_nm = UPER(k_uper) * UPER(k_uper)  * Jn(IHARM + 1, n, m) * epsx     &
!    &               + UPER(k_uper) * UPER(k_uper)  * Jn(IHARM - 1, n, m) * epsy     &
!    &               + SQ2 * UPER(k_uper) * UPAR0  * Jn(IHARM, n, m)     * epsz

!          sumb_31_nm = SQ2 * UPER(k_uper) * UPAR0  * Jn(IHARM + 1, n, m) * epsx     &
!     &               + SQ2 * UPER(k_uper) * UPAR0  * Jn(IHARM - 1, n, m) * epsy     &
!     &               + 2.0 * UPAR0 * UPAR0 * Jn(IHARM, n, m)     * epsz
     
          sumb_11_nm = sumb_11_nm * factor 
          sumb_31_nm = sumb_31_nm * factor 
                                
          sume_11_nm = sumb_11_nm * facte
          sume_31_nm = sumb_31_nm * facte
          
          sumc_11_nm = sume_11_nm * sinth_inv
          sumc_31_nm = sume_31_nm * sinth_inv     
     
          sumf_11_nm = sumc_11_nm * facte
          sumf_31_nm = sumc_31_nm * facte
                        
!          sumwdot_11_nm = sumb_11_nm * u0
!          sumwdot_31_nm = sumb_31_nm * u0      
                                                                                             
          sumf_11(i_uprl) = sumf_11(i_uprl) + sumf_11_nm
          sumf_31(i_uprl) = sumf_31(i_uprl) + sumf_31_nm

          sume_11(i_uprl) = sume_11(i_uprl) + sume_11_nm
          sume_31(i_uprl) = sume_31(i_uprl) + sume_31_nm

          sumc_11(i_uprl) = sumc_11(i_uprl) + sumc_11_nm
          sumc_31(i_uprl) = sumc_31(i_uprl) + sumc_31_nm

          sumb_11(i_uprl) = sumb_11(i_uprl) + sumb_11_nm
          sumb_31(i_uprl) = sumb_31(i_uprl) + sumb_31_nm
                

                
          sumwdotkx_11 = sumwdotkx_11 + xkxsav(n) * sumwdot_11_nm
          sumwdotkx_31 = sumwdotkx_31 + xkxsav(n) * sumwdot_31_nm       
                
          sumwdotky_11 = sumwdotky_11 + xkysav(m) * sumwdot_11_nm
          sumwdotky_31 = sumwdotky_31 + xkysav(m) * sumwdot_31_nm
          
          sumwdot_11 = sumwdot_11 + (ttmp + cross2) * factor * u0
          sumwdot_31 = sumwdot_31 + (cross3 + ld  ) * factor * u0         
          
          sum1_1 = sum1_1 + ttmp   * factor * u0   
     
          sum1_2 = sum1_2 + cross2 * factor * u0
     
          sum1_3 = sum1_3 + cross3 * factor * u0         
               
          sum1_4 = sum1_4 + ld     * factor * u0                                                        
                
      enddo  !  End sum over resonant modes
              
      ttmp_11 = 0.0 
      cross_13 = 0.0
      cross_31 = 0.0
      landau_33 = 0.0              
                        
      ttmp_11 = (sum2_1 + sum2_2) * sum1_1 
      cross_13 = (sum2_1 + sum2_2) * sum1_2 
      cross_31 = sum2_3 * sum1_3
      landau_33 = sum2_3 * sum1_4
      
      sum_wdot = sum_wdot + ttmp_11 + cross_13 + cross_31 + landau_33       
              
 
     
     
       sum_fx0 = sum_fx0  + sumkx2_1 * sumwdot_11    &
     &                    + sumkx2_2 * sumwdot_11    &
     &                    + sumkx2_3 * sumwdot_31    &
     &                    + sum2_1 * sumwdotkx_11    &
     &                    + sum2_2 * sumwdotkx_11    &
     &                    + sum2_3 * sumwdotkx_31
     
       sum_fy0 = sum_fy0  + sumky2_1 * sumwdot_11    &
     &                    + sumky2_2 * sumwdot_11    &
     &                    + sumky2_3 * sumwdot_31    &
     &                    + sum2_1 * sumwdotky_11    &
     &                    + sum2_2 * sumwdotky_11    &
     &                    + sum2_3 * sumwdotky_31    
      
       

       do i_uprl = 1, nupar
          b_sum(i_uprl) = b_sum(i_uprl) + sum2_1 * sumb_11(i_uprl)   &
     &                                  + sum2_2 * sumb_11(i_uprl)   &
     &                                  + sum2_3 * sumb_31(i_uprl)
          c_sum(i_uprl) = c_sum(i_uprl) + sum2_1 * sumc_11(i_uprl)   &
     &                                  + sum2_2 * sumc_11(i_uprl)   &
     &                                  + sum2_3 * sumc_31(i_uprl)
          e_sum(i_uprl) = e_sum(i_uprl) + sum2_1 * sume_11(i_uprl)   &
     &                                  + sum2_2 * sume_11(i_uprl)   &
     &                                  + sum2_3 * sume_31(i_uprl)
          f_sum(i_uprl) = f_sum(i_uprl) + sum2_1 * sumf_11(i_uprl)   &
     &                                  + sum2_2 * sumf_11(i_uprl)   &
     &                                  + sum2_3 * sumf_31(i_uprl)     
       end do
       
       

    end do  ! end sum over harmonics
    


    deallocate( zbeta )
    deallocate( zbeta_iharm )
    
    deallocate(nres )
    deallocate(mres )

    deallocate(zetai)
    deallocate(Jni)
    deallocate(Jn)    
    deallocate(NPARA_sav)
    
    
    deallocate(sumb_11)
    deallocate(sumb_31)
    
    deallocate(sumc_11)
    deallocate(sumc_31) 
    
    deallocate(sume_11)
    deallocate(sume_31)
    
    deallocate(sumf_11)
    deallocate(sumf_31)      
    
    deallocate( is_resonance_nm )
    
    return

  end subroutine QLSUM_MAXWELLIAN

!
!*************************************************************************
!

    subroutine QLSUM_NON_MAXWELLIAN(k_uper, b_sum, c_sum, e_sum, f_sum, &
       & sum_wdot, sum_fx0, sum_fy0, W, ZSPEC, ASPEC, BMAG, &
       & lmax, ENORM, UPARMIN, UPARMAX, &
       & NUPAR, NUPER, UPER, UPAR, DFDUPER, DFDUPAR,  &
       & ealphak, ebetak, ebk, nkdim1, nkdim2, mkdim1, mkdim2,   &
       & nkx1, nkx2, nky1, nky2, &
       & uxx, uxy, uxz, &
       & uyx, uyy, uyz, &
       & uzx, uzy, uzz, &
       & nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i_global, j_global, &
       & lmaxdim, ndist, nzeta, &
       & gradprlb, bmod, omgc, alpha, xm, upshift, xk_cutoff, rt, nphi, rho)

    implicit none

    integer i_global, j_global, ier, nxdim, nydim, k, lmaxdim, ndist
    integer, intent(IN):: NUPAR, NUPER, lmax
    integer nkx1, nkx2, nky1, nky2, j_upar, k_uper, l
    integer nkdim1, nkdim2, mkdim1, mkdim2
    integer:: NHARM, IHARM, M, N, i, nzeta
    integer i_uprl, upshift
    integer ires, iresmax, nphi

    complex, dimension(:,:), allocatable :: zbeta
    complex, dimension(:,:), allocatable :: zbeta_iharm
    
    real  y, y0, alpha, xm, akprl, sgn_kprl, omgc, bmod, xkprl_eff,  &
   &   descrim, xme, fgam, gradprlb, gammab, xk_cutoff, rt, akprl_min, rho

    real  uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz
    real  xkphi, sinth, factc, facte, factf, sinth_inv
    real  xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
    real  xkperpn, xkperpni, xkrhon, xketan, xkprln, beta
    real, intent(IN):: W, ZSPEC, ASPEC, BMAG
    real, intent(IN):: ENORM, UPARMIN, UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER, DFDUPAR
    real:: W2, WCW, RRP, RRM, WC, WCI
    real:: MUT0, SQMUT0, PISQMUT0, SQMUT0I, KPARA1
    real:: ISQ2, SQ2, NWCW, DFACTPAR, DFACTPER, U0
    real:: UPAR0, dfdupar0, dfduper0, du, dui, p 
    real:: time, t1, tmsec, second1, dummy, WI, uperpk, uperpk2
    real:: dzeta, dzetai, zetamax, zetamin, zeta0, zetamax1, zetamax2
    real:: A1, A2, A3, u
    real:: temp1, temp2, temp3, factor
    real:: temp1w, temp2w, temp3w
    
    real, dimension(:),     allocatable :: zetai
    real, dimension(:,:),   allocatable :: Jni
    real, dimension(:,:,:), allocatable :: Jn
    real, dimension(:,:),   allocatable :: NPARA_sav
    
    integer, dimension(:),  allocatable :: nres
    integer, dimension(:),  allocatable :: mres
        
    complex, dimension(:),  allocatable :: sumb_11
    complex, dimension(:),  allocatable :: sumb_31
     
    complex, dimension(:),  allocatable :: sumc_11
    complex, dimension(:),  allocatable :: sumc_31
    
    complex, dimension(:),  allocatable :: sume_11
    complex, dimension(:),  allocatable :: sume_31
     
    complex, dimension(:),  allocatable :: sumf_11
    complex, dimension(:),  allocatable :: sumf_31
    
    complex sumf_11_nm
    complex sumf_31_nm

    complex sume_11_nm
    complex sume_31_nm
       
    complex sumc_11_nm
    complex sumc_31_nm
       
    complex sumb_11_nm
    complex sumb_31_nm
 
    logical, dimension(:,:), allocatable :: is_resonance_nm      

    complex epsx, epsy, epsz
    complex ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
    &        ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
    &           ebk(nkdim1 : nkdim2, mkdim1 : mkdim2)
    complex cexpkx, cexpky, zi, zeta

    complex xx(nkdim1 : nkdim2, 1 : nxdim),   &
     &      yy(mkdim1 : mkdim2, 1 : nydim)

    complex cexpn, cexpnp1, cexpnm1, cexp11

    complex cexp1, cexp2, cexp0
    complex sumwdot_11_nm, sumwdot_31_nm
    
    complex sumwdot_11, sumwdot_31
    complex sumwdotkx_11, sumwdotkx_31
    complex sumwdotky_11, sumwdotky_31
        
    complex sum2_1, sum2_2, sum2_3
    complex sumkx2_1, sumkx2_2, sumkx2_3 
    complex sumky2_1, sumky2_2, sumky2_3        
    
    complex b_sum(nupar), c_sum(nupar), e_sum(nupar), f_sum(nupar)
    complex sum_wdot, sum_fx0, sum_fy0
    complex b(100)

    real, parameter:: EOVERAMU=9.64853e7
    real, parameter:: EOVERMH = 9.58084e+07
    
    real, parameter:: MPC2 = 938271998.38
    real, parameter:: C = 2.99792458e8
    real, parameter:: PI = 3.141592653597932384
    real :: cosbeta_n_m, sinbeta_n_m
    common/upcom/akprl_min    
    
    allocate( zbeta(nkx1:nkx2,nky1:nky2) )
    allocate( zbeta_iharm(nkx1:nkx2,nky1:nky2) )

    allocate(zetai(nzeta + 1) )
    allocate(Jni(-lmaxdim : lmaxdim, nzeta + 1) )
    allocate( Jn(-lmaxdim : lmaxdim, nkdim1 : nkdim2, mkdim1 : mkdim2))    
    allocate(NPARA_sav(nkdim1 : nkdim2, mkdim1 : mkdim2) ) 
    
    allocate(nres(nxdim * nydim) )
    allocate(mres(nxdim * nydim) )
    
    allocate(sumb_11(nupar))
    allocate(sumb_31(nupar))
     
    allocate(sumc_11(nupar))
    allocate(sumc_31(nupar))
    
    allocate(sume_11(nupar))
    allocate(sume_31(nupar))
     
    allocate(sumf_11(nupar))
    allocate(sumf_31(nupar))

    allocate(is_resonance_nm(nkx1:nkx2,nky1:nky2))

!   -------------------------------------
!   initialize allocatable arrays to zero
!   -------------------------------------

    zbeta = 0.0
    zbeta_iharm = 0.0

    sumf_11_nm = 0.0
    sumf_31_nm = 0.0

    sume_11_nm = 0.0
    sume_31_nm = 0.0
    
    sumc_11_nm = 0.0
    sumc_31_nm = 0.0
    
    sumb_11_nm = 0.0
    sumb_31_nm = 0.0

    is_resonance_nm = .false.

    zetai = 0.0
    Jni = 0.0
    NPARA_sav = 0.0

    sumb_11 = 0.0
    sumb_31 = 0.0
      
    sumc_11 = 0.0
    sumc_31 = 0.0

    sume_11 = 0.0
    sume_31 = 0.0

    sumf_11 = 0.0
    sumf_31 = 0.0    
    
    xme = 9.11e-31
    zi = cmplx(0., 1.)
    
    uperpk = uper(k_uper)
    uperpk2 = uperpk**2

    W2 = W * W
    WI = 1.0 / W

    WCW = BMAG * ZSPEC * EOVERMH / ASPEC / W
    WC = WCW * W
    WCI = 1.0 /WC

    MUT0 = 0.5 * MPC2 * ASPEC / ENORM
    SQMUT0 = SQRT(MUT0)
    SQMUT0I = 1.0 / SQMUT0
    PISQMUT0 = SQMUT0 * pi

    ISQ2 = SQRT(0.5)
    SQ2 = SQRT(2.0)
    NHARM = lmax
    
    du = (upar(nupar) - upar(1)) / (nupar - 1)
    dui = 1.0 / du
    
    if(nzeta .eq. 1)then
    
    ! -------------------------------------------------------- !
    ! ---Don't interpolate: precalculate all Bessel functions- !
    ! -------------------------------------------------------- !

       do n = nkx1, nkx2
          do m = nky1, nky2

             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             xkprln   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi
             xkperpn = sqrt(xkrhon**2 + xketan**2)
             
            ! ------------------------------------
            ! Optional: leave out upshift in xkprl
            ! --------------------------------- --          
             if (upshift .eq. 0)xkprln = uzz * xkphi
!            if (upshift .eq. 0)xkprln = nphi / rt
              
             if (upshift .eq. -1) then      
                if (xkperpn  .gt. xk_cutoff) xkprln = uzz * xkphi
             end if
             
                     
             
             sgn_kprl = sign(1.0, xkprln)
             akprl = abs(xkprln)                      
                     
!            ----------------------------------------------
!            Optional: Don't allow xkprl to be 0 (upshift = -2)
!            ----------------------------------------------        
             if (upshift .eq. -2) then
                if (akprl .lt. akprl_min) then
                    xkprln = akprl_min * sgn_kprl
                end if 
             end if                  
                     

             y0 = 1.5
             y = y0
             
      
             l = 1
             if(xkprln .eq. 0)xkprln = 1.0e-06
             gammab = abs(l * omgc / (2.0 * alpha * xkprln**2)  &
     &                                    * gradprlb / bmod)

             if(xm .eq. xme)gammab = 0.0
!             if(abs(gammab) .gt. 1000.0) gammab = 1000.0
             if(abs(gammab) .lt. .01)gammab = .01


             if(sgn_kprl .ge. 0.0)then
                fgam = 1.0

                if(gammab .gt. 1.0e-05)then
                   y = y0
                   fgam = (sqrt(1. +  4. * gammab * y) - 1.)   &
     &               / (2. * gammab * y)
                endif

                xkprl_eff = xkprln / fgam 

             end if


             if(sgn_kprl .lt. 0.0)then
                fgam = 1.0

                if(gammab .gt. 1.0e-05)then
                   descrim = 1. - 4. * gammab * y0
                   if (descrim .ge. 0.0) y =   y0
                   if (descrim .lt. 0.0) y = - y0
                   fgam = (1. - sqrt(1. -  4. * gammab * y) )  &
     &                / (2. * gammab * y)
                endif

                xkprl_eff = xkprln / fgam 

             end if

             if (upshift .ne. 0) xkprln = xkprl_eff          
        
             NPARA_sav(n, m) = xkprln * C / W

             xkperpn = sqrt(xkrhon**2 + xketan**2)
             if(xkperpn .eq. 0.0)xkperpn = 1.0e-08
        
             cosbeta_n_m  = xkrhon / xkperpn
             sinbeta_n_m  = xketan / xkperpn
             zbeta(n,m) = cmplx( cosbeta_n_m , sinbeta_n_m  )

             zeta = xkperpn * uper(k_uper) * c * sqmut0i / wc

             call besjc(zeta, nharm + 2, b, ier)
             if(ier .ne. 0) write(6, *) "ier = ", ier

             do IHARM = 0, NHARM + 1
                Jn(iharm,  n, m) = b(iharm + 1)
                Jn(-iharm, n, m) = (-1.0)**iharm * Jn(iharm, n, m)
             end do

          end do
       end do
    
    else
    
        
       ! -------------------------------------- !
       ! ---Interpolate; calculate zeta mesh -- !
       ! -------------------------------------- !
    
       zetamax = 0.0
       zetamin = 0.0
    
       do n = nkx1, nkx2
          do m = nky1, nky2

             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             xkprln   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi

             xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08
        
             zeta0 = xkperpn * uperpk * c * sqmut0i * wci
          
             if (zeta0 .gt. zetamax) zetamax = zeta0
             if (zeta0 .lt. zetamin) zetamin = zeta0
          end do
       end do
    

        if(zetamax .eq. zetamin)then
          zetamax =  1.0e-06
          zetamin = -1.0e-06
       end if

       dzeta = (zetamax - zetamin) / (nzeta - 1)
       dzetai = 1.0 / dzeta
    
       ! ------------------------------------------------- !
       ! ---Pre-calculate Bessel functions on zeta mesh -- !
       ! ------------------------------------------------- !
          
       do i = 1, nzeta + 1
          zetai(i) = zetamin + (i - 1) * dzeta
          zeta = cmplx(zetai(i), 0.0)
          
          call besjc(zeta, nharm + 2, b, ier)
!          if(ier .ne. 0) write(6, *) "ier = ", ier
          
          do iharm = 0, NHARM + 1
             Jni(iharm,  i) = b(iharm + 1)
             Jni(-iharm, i) = (-1.0)**iharm * b(iharm + 1)
          end do
       end do
     

       ! --------------------------------- !
       ! ---Interpolate Bessel functions-- !
       ! --------------------------------- !

       do n = nkx1, nkx2
          do m = nky1, nky2

             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             xkprln   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi
             xkperpn = sqrt(xkrhon**2 + xketan**2)
             
            ! ------------------------------------
            ! Optional: leave out upshift in xkprl
            ! --------------------------------- --          
             if (upshift .eq. 0)xkprln = uzz * xkphi
!            if (upshift .eq. 0)xkprln = nphi / rt
             
             if (upshift .eq. -1) then      
                if (xkperpn  .gt. xk_cutoff) xkprln = uzz * xkphi
             end if
             
             
             
             sgn_kprl = sign(1.0, xkprln)
             akprl = abs(xkprln)
                      
                     
!            ----------------------------------------------
!            Optional: Don't allow xkprl to be 0 (upshift = -2)
!            ----------------------------------------------        
             if (upshift .eq. -2) then
                if (akprl .lt. akprl_min) then
                    xkprln = akprl_min * sgn_kprl
                end if 
             end if                  
                     

             y0 = 1.5
             y = y0
             
      
             l = 1
             if(xkprln .eq. 0)xkprln = 1.0e-06
             gammab = abs(l * omgc / (2.0 * alpha * xkprln**2)  &
     &                                    * gradprlb / bmod)

             if(xm .eq. xme)gammab = 0.0
!             if(abs(gammab) .gt. 1000.0) gammab = 1000.0
             if(abs(gammab) .lt. .01)gammab = .01


             if(sgn_kprl .ge. 0.0)then
                fgam = 1.0

                if(gammab .gt. 1.0e-05)then
                   y = y0
                   fgam = (sqrt(1. +  4. * gammab * y) - 1.)   &
     &               / (2. * gammab * y)
                endif

                xkprl_eff = xkprln / fgam 

             end if


             if(sgn_kprl .lt. 0.0)then
                fgam = 1.0

                if(gammab .gt. 1.0e-05)then
                   descrim = 1. - 4. * gammab * y0
                   if (descrim .ge. 0.0) y =   y0
                   if (descrim .lt. 0.0) y = - y0
                   fgam = (1. - sqrt(1. -  4. * gammab * y) )  &
     &                / (2. * gammab * y)
                endif

                xkprl_eff = xkprln / fgam 

             end if

             if (upshift .ne. 0) xkprln = xkprl_eff
             
        
             NPARA_sav(n, m) = xkprln * C * WI

             xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08
             xkperpni = 1.0 / xkperpn
        
             cosbeta_n_m = xkrhon * xkperpni
             sinbeta_n_m  = xketan * xkperpni
             zbeta(n,m) = cmplx( cosbeta_n_m , sinbeta_n_m )
          
             zeta0 = xkperpn * uperpk * c * sqmut0i * wci
          
             i = int((zeta0 - zetamin) * dzetai) + 1
             p = (zeta0 - zetai(i)) * dzetai
             A1 = 0.5 * P * (P - 1.)
             A2 = 1. - P * P
             A3 = 0.5 * P * (P + 1.)
          
             do iharm = -NHARM - 1, NHARM + 1
          
                Jn(iharm, n, m) = Jni(iharm, i)    &
     &             + p * (Jni(iharm, i + 1) - Jni(iharm, i))
                if(i .ne. 1 )then
                   Jn(iharm, n, m) = A1 * Jni(iharm, i - 1)     &
     &                             + A2 * Jni(iharm, i)         &
     &                             + A3 * Jni(iharm, i + 1)
                end if
              
             end do
          
          end do
       end do
    

    end if
                    

    ! ------------------------ !
    ! ---Sum over harmonics--- !
    ! ------------------------ !
    
    sum_wdot = 0.0
    sum_fx0  = 0.0
    sum_fy0  = 0.0    
    
    b_sum = 0.0
    c_sum = 0.0
    e_sum = 0.0
    f_sum = 0.0
       
    do IHARM = -NHARM, NHARM

       NWCW = real(IHARM) * WCW
       
       sumb_11 = 0.0
       sumb_31 = 0.0
          
       sumc_11 = 0.0
       sumc_31 = 0.0    
          
       sume_11 = 0.0
       sume_31 = 0.0
          
       sumf_11 = 0.0
       sumf_31 = 0.0      
                  
       sumwdot_11 = 0.0
       sumwdot_31 = 0.0       
       
       sumwdotkx_11 = 0.0
       sumwdotkx_31 = 0.0
       
       sumwdotky_11 = 0.0
       sumwdotky_31 = 0.0       

       sum2_1 = 0.0
       sum2_2 = 0.0
       sum2_3 = 0.0
       
       sumkx2_1 = 0.0
       sumkx2_2 = 0.0
       sumkx2_3 = 0.0
       
       sumky2_1 = 0.0
       sumky2_2 = 0.0
       sumky2_3 = 0.0       
       
       call zpow((nkx2-nkx1+1)*(nky2-nky1+1), zbeta, iharm, zbeta_iharm) 
       
       ! ----------------------- !
       ! ---Find resonant modes--!
       ! ----------------------- !
       ires = 0
       do n = nkx1, nkx2
          do m = nky1, nky2     
             ! ------------------------ !
             ! -- Resonance relation -- !
             ! ------------------------ !
             RRP = 1.0 - NWCW - NPARA_sav(n, m) * UPARMAX * SQMUT0i
             RRM = 1.0 - NWCW - NPARA_sav(n, m) * UPARMIN * SQMUT0i

             is_resonance_nm(n,m) = (RRP * RRM .le. 0.0)
             if (is_resonance_nm(n,m)) then
                 ires = ires + 1
                 nres(ires) = n
                 mres(ires) = m
             end if
          end do
       end do
       iresmax = ires 
       
!        go to 1000
      ! ------------------------------- !
      ! ---Sum over all Fourier modes-- !
      ! --------------------------------!
       do n = nkx1, nkx2
          do m = nky1, nky2
!            cexp1 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m) * zbeta(n,m)
!            cexp2 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m) / zbeta(n,m)
             cexp0 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m)
             cexp1 = cexp0 * zbeta(n,m)
             cexp2 = cexp0 / zbeta(n,m)
                
             epsx = isq2 * (ealphak(n, m) - zi * ebetak(n, m)) * cexp1
             epsy = isq2 * (ealphak(n, m) + zi * ebetak(n, m)) * cexp2
             epsz = ebk(n, m) * cexp0
             
             sum2_1 = sum2_1 + conjg(epsx) * Jn(IHARM + 1, n, m)
             sum2_2 = sum2_2 + conjg(epsy) * Jn(IHARM - 1, n, m)
             sum2_3 = sum2_3 + conjg(epsz) * Jn(IHARM, n, m)         
          end do
       end do
! 1000  continue

      ! ---------------------------- !
      ! ---Sum over resonant modes-- !
      ! ---------------------------- !

       do ires = 1, iresmax
          n = nres(ires)
          m = mres(ires)
                                                         
          cexp1 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m) * zbeta(n,m)
          cexp2 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m) / zbeta(n,m)
          cexp0 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m)
                
          epsx = isq2 * (ealphak(n, m) - zi * ebetak(n, m)) * cexp1
          epsy = isq2 * (ealphak(n, m) + zi * ebetak(n, m)) * cexp2
          epsz = ebk(n, m) * cexp0
             
!          sum2_1 = sum2_1 + conjg(epsx) * Jn(IHARM + 1, n, m)
!          sum2_2 = sum2_2 + conjg(epsy) * Jn(IHARM - 1, n, m)
!          sum2_3 = sum2_3 + conjg(epsz) * Jn(IHARM, n, m)
             
          sumkx2_1 = sumkx2_1 + xkxsav(n) * conjg(epsx) * Jn(IHARM + 1, n, m)
          sumkx2_2 = sumkx2_2 + xkxsav(n) * conjg(epsy) * Jn(IHARM - 1, n, m)
          sumkx2_3 = sumkx2_3 + xkxsav(n) * conjg(epsz) * Jn(IHARM, n, m)
             
          sumky2_1 = sumky2_1 + xkysav(m) * conjg(epsx) * Jn(IHARM + 1, n, m)
          sumky2_2 = sumky2_2 + xkysav(m) * conjg(epsy) * Jn(IHARM - 1, n, m)
          sumky2_3 = sumky2_3 + xkysav(m) * conjg(epsz) * Jn(IHARM, n, m)                           

          UPAR0 = SQMUT0 / NPARA_sav(n, m) * (1. - NWCW)
        
          u = sqrt(upar0**2 + uperpk2) + 1.0e-08
          sinth = uperpk / u + 1.0e-08
          sinth_inv = 1.0 / sinth
                
          facte = (nwcw - sinth**2) /  upar0
        
          i = int((UPAR0 - UPAR(1)) * dui) + 1
          i_uprl = i
          p = (UPAR0 - UPAR(i)) * dui
                
          dfduper0 = dfduper(k_uper, NUPAR)
          if (i .ne. NUPAR) then
            dfduper0 = dfduper(k_uper, i) + (dfduper(k_uper, i+1) - dfduper(k_uper, i)) * p
          end if
          
                
          DFACTPAR = NPARA_sav(n, m) * UPAR0 * SQMUT0I
          DFACTPER = NPARA_sav(n, m) * UPER(k_uper) * SQMUT0I
          
          dfdupar0 = dfdupar(k_uper, NUPAR)
          if(i .ne. NUPAR)then
             dfdupar0 = dfdupar(k_uper, i) +                    &
     &       (dfdupar(k_uper, i+1) - dfdupar(k_uper, i)) * p
          end if          
        
          U0 = (1. - DFACTPAR) * DFDUPER0 + DFACTPER * DFDUPAR0   
                        
          factor = PISQMUT0 / abs(NPARA_sav(n, m)) 
                        
          sumb_11_nm = UPER(k_uper) * UPER(k_uper)  * Jn(IHARM + 1, n, m) * epsx     &
     &               + UPER(k_uper) * UPER(k_uper)  * Jn(IHARM - 1, n, m) * epsy     &
     &               + SQ2 * UPER(k_uper) * UPAR0  * Jn(IHARM, n, m)     * epsz

          sumb_31_nm = SQ2 * UPER(k_uper) * UPAR0  * Jn(IHARM + 1, n, m) * epsx     &
     &               + SQ2 * UPER(k_uper) * UPAR0  * Jn(IHARM - 1, n, m) * epsy     &
     &               + 2.0 * UPAR0 * UPAR0 * Jn(IHARM, n, m)     * epsz
     
          sumb_11_nm = sumb_11_nm * factor 
          sumb_31_nm = sumb_31_nm * factor 
                                
          sume_11_nm = sumb_11_nm * facte
          sume_31_nm = sumb_31_nm * facte
          
          sumc_11_nm = sume_11_nm * sinth_inv
          sumc_31_nm = sume_31_nm * sinth_inv     
     
          sumf_11_nm = sumc_11_nm * facte
          sumf_31_nm = sumc_31_nm * facte
                        
          sumwdot_11_nm = sumb_11_nm * u0
          sumwdot_31_nm = sumb_31_nm * u0       
                                                                                             
          sumf_11(i_uprl) = sumf_11(i_uprl) + sumf_11_nm
          sumf_31(i_uprl) = sumf_31(i_uprl) + sumf_31_nm

          sume_11(i_uprl) = sume_11(i_uprl) + sume_11_nm
          sume_31(i_uprl) = sume_31(i_uprl) + sume_31_nm

          sumc_11(i_uprl) = sumc_11(i_uprl) + sumc_11_nm
          sumc_31(i_uprl) = sumc_31(i_uprl) + sumc_31_nm

          sumb_11(i_uprl) = sumb_11(i_uprl) + sumb_11_nm
          sumb_31(i_uprl) = sumb_31(i_uprl) + sumb_31_nm
                
          sumwdot_11 = sumwdot_11 + sumwdot_11_nm
          sumwdot_31 = sumwdot_31 + sumwdot_31_nm
                
          sumwdotkx_11 = sumwdotkx_11 + xkxsav(n) * sumwdot_11_nm
          sumwdotkx_31 = sumwdotkx_31 + xkxsav(n) * sumwdot_31_nm       
                
          sumwdotky_11 = sumwdotky_11 + xkysav(m) * sumwdot_11_nm
          sumwdotky_31 = sumwdotky_31 + xkysav(m) * sumwdot_31_nm                                                       
                
       enddo
       
       
       sum_wdot = sum_wdot + sum2_1 * sumwdot_11    &
     &                     + sum2_2 * sumwdot_11    &
     &                     + sum2_3 * sumwdot_31 
     
     
       sum_fx0 = sum_fx0  + sumkx2_1 * sumwdot_11    &
     &                    + sumkx2_2 * sumwdot_11    &
     &                    + sumkx2_3 * sumwdot_31    &
     &                    + sum2_1 * sumwdotkx_11    &
     &                    + sum2_2 * sumwdotkx_11    &
     &                    + sum2_3 * sumwdotkx_31
     
       sum_fy0 = sum_fy0  + sumky2_1 * sumwdot_11    &
     &                    + sumky2_2 * sumwdot_11    &
     &                    + sumky2_3 * sumwdot_31    &
     &                    + sum2_1 * sumwdotky_11    &
     &                    + sum2_2 * sumwdotky_11    &
     &                    + sum2_3 * sumwdotky_31    
      
       

       do i_uprl = 1, nupar
          b_sum(i_uprl) = b_sum(i_uprl) + sum2_1 * sumb_11(i_uprl)   &
     &                                  + sum2_2 * sumb_11(i_uprl)   &
     &                                  + sum2_3 * sumb_31(i_uprl)
          c_sum(i_uprl) = c_sum(i_uprl) + sum2_1 * sumc_11(i_uprl)   &
     &                                  + sum2_2 * sumc_11(i_uprl)   &
     &                                  + sum2_3 * sumc_31(i_uprl)
          e_sum(i_uprl) = e_sum(i_uprl) + sum2_1 * sume_11(i_uprl)   &
     &                                  + sum2_2 * sume_11(i_uprl)   &
     &                                  + sum2_3 * sume_31(i_uprl)
          f_sum(i_uprl) = f_sum(i_uprl) + sum2_1 * sumf_11(i_uprl)   &
     &                                  + sum2_2 * sumf_11(i_uprl)   &
     &                                  + sum2_3 * sumf_31(i_uprl)     
       end do
       
       

    end do
    


    deallocate( zbeta )
    deallocate( zbeta_iharm )
    
    deallocate(nres )
    deallocate(mres )

    deallocate(zetai)
    deallocate(Jni)
    deallocate(Jn)    
    deallocate(NPARA_sav)
    
    
    deallocate(sumb_11)
    deallocate(sumb_31)
    
    deallocate(sumc_11)
    deallocate(sumc_31) 
    
    deallocate(sume_11)
    deallocate(sume_31)
    
    deallocate(sumf_11)
    deallocate(sumf_31)      
    
    deallocate( is_resonance_nm )
    
    return

  end subroutine QLSUM_NON_MAXWELLIAN


!
!*************************************************************************
!

        subroutine zpow(n, z, iharm, zout )
        implicit none
        integer n, iharm
        complex z(n), zout(n)

        integer i,ipow
        complex one, zero
        complex zin

        integer nharm
        logical isodd
        intrinsic mod

        logical use_zdiv
        parameter(use_zdiv=.true.)

        integer nb
        parameter(nb=1024*4*4)
        complex zk(nb)
        integer istart,iend,isize

        one = 1.0d0
        zero = 0.0d0

        if (iharm.eq.0) then
             do i=1,n
                zout(i) = one
             enddo
             return
        endif


        do istart=1,n,nb

           iend = min(n,istart+nb-1)
           isize = iend-istart+1

          do i=1,isize
            zout(istart-1+i) = one
          enddo


         do i=1,isize
           zk(i) = z(istart-1+i)
         enddo

        nharm = abs(iharm)
        do while (nharm .gt. 0)
           isodd = (mod(nharm,2).eq.1)
           if (isodd) then
               do i=1,isize
                 zout(istart-1+i) = zout(istart-1+i) * zk(i)
               enddo
           endif
           do i=1,isize
              zk(i) = zk(i) * zk(i)
           enddo
           nharm = int( nharm/2 )
        enddo



        if (iharm.lt.0) then
           if (use_zdiv) then
             do i=1,isize
                zin = zout(istart-1+i)
                call zdiv( zin, zout(istart-1+i) )
            enddo
           else
             do i=1,isize
                zin = zout(istart-1+i)
                zout(istart-1+i) = one/zin
             enddo
           endif
        endif

        enddo

        return
        end subroutine 
        
!
!*************************************************************************
!

        subroutine zdiv( zin, zout )
        implicit none
        complex zin, zout
        real a, b
        real d

        real one
        parameter(one=1.0d0)
        real rd, a_over_b, b_over_a


        a = real(zin)
        b = aimag(zin)

!       z = (a + i * b)
!       1/z =  a/(a^2 + b^2) - i * b/(a^2 + b^2)
!
!       or    1/(a + (b/a)*b) - i * (b/a) / (a + (b/a)*b)
!       or    (a/b)/( (a/b)*a + b ) - i * 1/( (a/b)*a + b )
!        
        if (abs(a).gt.abs(b)) then
            b_over_a = b/a
            d = a + (b_over_a)*b
            rd = one/d
            zout = cmplx( rd, -(b_over_a)*rd )
        else
            a_over_b = a/b
            d = (a_over_b)*a + b
            rd = one/d
            zout = cmplx( (a_over_b)*rd, -rd )
        endif

        return
        end subroutine
        


!
!*************************************************************************
!

         end module qlsum_myra_mod
         

