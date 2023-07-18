      module ql_myra_mod
      use qlsum_myra_mod
!     include 'mpif.h'
      use mpi
      contains
      
      
c
c***************************************************************************
c

      subroutine ql_myra_write(bqlavg, cqlavg, eqlavg, fqlavg, 
     &   wdot_inout, fx0_inout, fy0_inout, fz0_inout, 
     &   vol, nrhodim, nnoderho, drho,
     &   dx, dy, r0, nxdim, nydim, nnodex, nnodey,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     &   xm, q, xn, xkt, omgc, omgp2, lmax,
     &   xkxsav, xkysav, nzfun, ibessel,
     &   exk, eyk, ezk, nphi, capr,
     &   bxn, byn, bzn,
     &   uxx, uxy, uxz,
     &   uyx, uyy, uyz,
     &   uzx, uzy, uzz,
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     &   xx, yy, isigma, xnuomga, psi, psilim, nboundary,
     &   myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper, dfdupar, fperp, 
     &   UminPara_cql, UmaxPara_cql, UPERP_cql, UPARA_cql, UPERP, UPARA,
     &   vc_mks_cql, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj,
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     &   nkperp, lmaxdim, nzeta_wdot, theta_,
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, 
     &   n_bin, upshift, i_write, xk_cutoff)

*----------------------------------------------------------------------
*     This subroutine calculates wdot and the quasi-linear operator for 
*     a single species
*----------------------------------------------------------------------

      implicit none

      integer splitrank, nloops, partition, start, finish
      integer pstart, pfinish, status(MPI_STATUS_SIZE), ierr
      integer recv_size, remainder
c      real, allocatable :: bql_store(:,:,:)
c      real, allocatable :: cql_store(:,:,:)
c      real, allocatable :: eql_store(:,:,:)
c      real, allocatable :: fql_store(:,:,:)
      logical iam_root
      integer left_neighbor, right_neighbor
      real token
      real upara_test
      integer mi_max, mi_min

      logical ismine
      
      real u2, fnorm, f_cql
      
      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0
      
      integer n_theta_max, n_psi_max, i_psi_eq, n_bin, upshift
      integer n_theta_(n_psi_max), n_theta, ntheta_giv, i_write

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp, lmaxdim
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     &   icontxt, nboundary, nzeta_wdot
      integer dlen_, desc_amat(dlen_)
      integer nxdim, nydim, nnodex, nnodey, nrhodim, nnoderho,
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma
      integer ftrap, k, ni0, mi0
      
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt, xk_cutoff
      
      real theta_(n_theta_max, n_psi_max), thegiv, dthetag
      real deriv, dtheta, dtau_ratio_giv, tau_bounce_giv
      real uperp0, upara0, xlamda, derivb, dtheta0, factor, upara_mi
      real duperp, dupara

      real zeffcd, clight, rt, b0, damp, r1, xm1, ceta, epsa, xmut2,
     &   eta0, xnexp, xjtild, signkz, cfit, c1, afit, yt, bmaxa, vphase,
     &   omgrf, c_ehst, akprl, xkprl, xkphi, a, rmaxa, rmina,
     &   wphase, xlnlam, vth
      real x(nxdim), y(nydim), argd, bratio, drho


      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim),
     &     psi(nxdim, nydim), psilim, alpha
      real capr(nxdim), xnuomg, delta0, dx, dy, r0
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim),
     &     bmod_mid(nxdim, nydim),xnuomga(nxdim, nydim) 

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     &        yy(mkdim1 : mkdim2, 1 : nydim)
     
      complex wdoti, fx0i, fy0i
      
      integer  :: n_psi_dim, nuper, nupar, n_psi, mi, ni

      real bqlavg(nuper, nupar, nnoderho)
      real cqlavg(nuper, nupar, nnoderho)
      real eqlavg(nuper, nupar, nnoderho)
      real fqlavg(nuper, nupar, nnoderho)
      
      real wdot_inout(nxdim, nydim)
      real  fx0_inout(nxdim, nydim)
      real  fy0_inout(nxdim, nydim)
      real  fz0_inout(nxdim, nydim)
      
      real wdot(nnodex, nnodey)
      real fx0(nnodex, nnodey)      
      real fy0(nnodex, nnodey)      
      real fz0(nnodex, nnodey)
      
!      real count(0 : 10000, 1), sum_count

      real vol(nrhodim), dldbavg(nrhodim)
        


      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz


      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        eyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real q, xn(nxdim, nydim), xkt(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim),
     &     uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim),
     &     uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)


      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: UPERP_cql(NUPER), UPARA_cql(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: FPERP(NUPER)

      real bql, cql, eql, fql

      complex, dimension(:),  allocatable :: b_sum
      complex, dimension(:),  allocatable :: c_sum
      complex, dimension(:),  allocatable :: e_sum
      complex, dimension(:),  allocatable :: f_sum
      
      complex, dimension(:),  allocatable :: wdot_sum
      complex, dimension(:),  allocatable :: sum_fx0
      complex, dimension(:),  allocatable :: sum_fy0
      
c      real, dimension(:,:,:), allocatable :: factvol
c      real, dimension(:,:), allocatable :: factvol2d

      real, dimension(:,:,:), allocatable :: bqlvol
      real, dimension(:,:), allocatable :: bqlvol2d

      real, dimension(:,:,:), allocatable :: cqlvol
      real, dimension(:,:), allocatable :: cqlvol2d

      real, dimension(:,:,:), allocatable :: eqlvol
      real, dimension(:,:), allocatable :: eqlvol2d

      real, dimension(:,:,:), allocatable :: fqlvol
      real, dimension(:,:), allocatable :: fqlvol2d

      real :: W, ENORM, ZSPEC, ASPEC, BMAG
      real :: UminPara,UmaxPara
      real :: UminPara_cql,UmaxPara_cql
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, vc_mks_cql, rho_a(n_psi_dim)
      real :: eps0, pi, emax, u0, u_0, u_, costh, costh0, psic

      parameter (eps0 = 8.85e-12)
      parameter (PI = 3.141592653597932384)


      integer :: nwork, ip
      logical :: has_work
      integer, dimension(nnodex*nnodey) :: i_table, j_table
      
      allocate( dfduper0(nuper, nupar) )
      allocate( dfdupar0(nuper, nupar) )
      
      allocate(b_sum(NUPAR) )
      allocate(c_sum(NUPAR) )
      allocate(e_sum(NUPAR) )
      allocate(f_sum(NUPAR) )
      
      allocate(wdot_sum(NUPER) )
      allocate(sum_fx0(NUPER) )
      allocate(sum_fy0(NUPER) )      
      
c      allocate(factvol(nuper, nupar, nnoderho) )
c      allocate(factvol2d(nuper, nupar) )

      allocate(bqlvol(nuper, nupar, nnoderho) )
      allocate(bqlvol2d(nuper, nupar) )

      allocate(cqlvol(nuper, nupar, nnoderho) )
      allocate(cqlvol2d(nuper, nupar) )

      allocate(eqlvol(nuper, nupar, nnoderho) )
      allocate(eqlvol2d(nuper, nupar) )

      allocate(fqlvol(nuper, nupar, nnoderho) )
      allocate(fqlvol2d(nuper, nupar) )

!     -------------------------------------
!efd  initialize allocatable arrays to zero
!     -------------------------------------
      wdot = 0.0
      fx0 = 0.0
      fy0 = 0.0
      fz0 = 0.0

      b_sum = 0.0
      c_sum = 0.0
      e_sum = 0.0
      f_sum = 0.0

      wdot_sum = 0.0
      sum_fx0 = 0.0
      sum_fy0 = 0.0

!      count = 0.0
      
c      factvol = 0.0
c      factvol2d = 0.0      

      bqlvol = 0.0
      bqlvol2d = 0.0

      cqlvol = 0.0
      cqlvol2d = 0.0

      eqlvol = 0.0
      eqlvol2d = 0.0

      fqlvol = 0.0
      fqlvol2d = 0.0




      W = omgrf
      ZSPEC = q / 1.6e-19

      if(ndist .eq. 0)then   !--Maxwellian--!

         do ni = 1, NUPER
            UPERP(ni) = (real(ni-1)/real(NUPER-1))
         end do

         UminPara = -1.0
         UmaxPara =  1.0

         do mi = 1, NUPAR
            UPARA(mi) = (-1.0 + 2. * (real(mi-1) / real(NUPAR-1)))
         end do

      else   !--non-Maxwellian--!

            vc_mks = vc_mks_cql
            UminPara = UminPara_cql
            UmaxPara = UmaxPara_cql

            do ni = 1, nuper
               uperp(ni) = uperp_cql(ni)
            end do

            do mi = 1, nupar
               upara(mi) = upara_cql(mi)
            end do

      end if


      do n = 1, nnoderho

         do mi = 1, nupar
            do ni = 1, nuper
            

               bqlvol(ni, mi, n) = 0.0
               cqlvol(ni, mi, n) = 0.0
               eqlvol(ni, mi, n) = 0.0
               fqlvol(ni, mi, n) = 0.0


               bqlavg(ni, mi, n) = 0.0
               cqlavg(ni, mi, n) = 0.0
               eqlavg(ni, mi, n) = 0.0
               fqlavg(ni, mi, n) = 0.0

            end do
         end do
      end do



      nwork = 0
      do j=1,nnodey
         do i=1,nnodex
           has_work = (psi(i,j) .le. psilim .and. nboundary .eq. 1
     &                                      .or. nboundary .eq. 0)
           if (has_work) then
              nwork = nwork + 1
              i_table(nwork) = i
              j_table(nwork) = j
           endif
         enddo
      enddo
      

*     -----------------------
*     Loop over spatial mesh:
*     -----------------------

      nloops = nwork
      partition=int(nloops/nproc)
      remainder=mod(nloops,nproc)
      splitrank=nproc-remainder
      if(myid.lt.splitrank)then
         start=1+myid*partition
         finish=start+partition-1
      else
         start=1+myid*partition+myid-splitrank
         finish=start+partition
      endif
      
c      if(ndist.eq.1) then
c             allocate(bql_store(start:finish,1:nuper,1:nupar),
c     &         cql_store(start:finish,1:nuper,1:nupar),
c     &         eql_store(start:finish,1:nuper,1:nupar),
c     &         fql_store(start:finish,1:nuper,1:nupar))
c      else
c      endif


      do ip = start,finish
      
            i = i_table(ip)
            j = j_table(ip)

        
            xkphi = nphi / capr(i)
            if(xkphi .eq. 0.0)xkphi = 1.0e-05
        
            alpha = sqrt(2.0 * xkt(i, j) / xm)
            n = int(rho(i,j) / drho) + 1
        
            if (ndist .eq. 0) then 
               vc_mks =  3.0 * alpha    !--Maxwellian only--!
               vc_mks_cql = vc_mks
               UminPara_cql = -1.0 
               UmaxPara_cql = 1.0
            end if
            
            u0 = vc_mks / alpha
                    
            Emax = 0.5 * xm * vc_mks**2
            Enorm = Emax / 1.6e-19
            ASPEC = xm / 1.67e-27
            BMAG = omgc(i,j) * (xm / q)
               
            bratio = bmod_mid(i,j) / bmod(i,j)
            if (bratio .gt. 1.0) bratio = 1.0
            
            duperp = (uperp(nuper) - uperp(1)) / (nuper - 1)
            dupara = (upara(nupar) - upara(1)) / (nupar - 1)
            
            psic = 1.0 / bratio

!           -----------------------------------------------------
!           get CQL3D distribution function on the midplane:
!           used for Wdot only - not the quasilinear coefficients
!           -----------------------------------------------------

            if(ndist .eq. 0)then   !--Maxwellian--!
        
               call maxwell_dist(u0, NUPAR, NUPER,
     &              UminPara, UmaxPara,
     &              UPERP, UPARA, DFDUPER, DFDUPAR, FPERP)

            else   !--non-Maxwellian--!
        
               call cql3d_dist(nupar, nuper, n_psi,
     &              n_psi_dim, rho_a, rho(i,j),
     &              UminPara,UmaxPara,
     &              df_cql_uprp, df_cql_uprl,
     &              UPERP, UPARA, DFDUPER0, DFDUPAR0)

!              ---------------------------------------------------------
!              map CQL3D distribution function off the midplane for Wdot
!              ---------------------------------------------------------
               if(bratio .gt. 0.0)then
               
                  dfduper = 0.0
                  dfdupar = 0.0
                       
                  do ni = 1, nuper
                     do mi = 1, nupar

                        argd =  uperp(ni)**2 * (1. - bratio)
     &                                                   + upara(mi)**2
                        if (argd .le. 0.0) argd = 1.0e-06
                                                
                        uperp0 = uperp(ni) * sqrt(bratio)
                        upara0  = sign(1.0, upara(mi)) * sqrt(argd)
                        
                        dfduper(ni, mi) = 0.0
                        dfdupar(ni, mi) = 0.0
                                
                        if(upara0 .ge. upara(1) .and. 
     &                                   upara0 .le. upara(nupar)) then                 
                                        
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
                           bi = dfduper0(ni0+1,mi0) - dfduper0(ni0,mi0)
                           ci = dfduper0(ni0,mi0+1) - dfduper0(ni0,mi0)
                           di = dfduper0(ni0+1,mi0+1)+ dfduper0(ni0,mi0) 
     &                        - dfduper0(ni0+1,mi0)- dfduper0(ni0,mi0+1) 
                           
                           dfduper0_intplt = ai + bi * zeta 
     &                                     + ci * eta + di * zeta * eta                         

                           ai = dfdupar0(ni0, mi0)
                           bi = dfdupar0(ni0+1,mi0) - dfdupar0(ni0,mi0)
                           ci = dfdupar0(ni0,mi0+1) - dfdupar0(ni0,mi0)
                           di = dfdupar0(ni0+1,mi0+1)+ dfdupar0(ni0,mi0) 
     &                        - dfdupar0(ni0+1,mi0)- dfdupar0(ni0,mi0+1) 
                           
                           dfdupar0_intplt = ai + bi * zeta 
     &                                     + ci * eta + di * zeta * eta 
     
                           end if                       
                                                        
                                        
                           if (upara0 .ne. 0.0)then
                        
                              dfdupar(ni, mi) = dfdupar0_intplt * 
     &                           upara(mi) / upara0
     
                              dfduper(ni, mi) = dfduper0_intplt * 
     &                           sqrt(bratio) + dfdupar0_intplt * 
     &                           uperp(ni) / upara0 * (1.0 - bratio)
     
c                             dfdth = upara(mi) * dfduper(ni, mi)
c     &                             - uperp(ni) * dfdupar(ni, mi)
     
                           end if
                           
                        end if
                        
                        
                        go to 5000                      
!                       ----------------------------
!                       optional analytic Maxwellian
!                       ----------------------------
                                        
                        alpha = sqrt(2.0 * xkt(i, j) / xm)
!                        vc_mks = 3.5 * alpha
                        u0 = vc_mks / alpha
                     
                        fnorm = u0**3 / pi**1.5 
                           
                        u2 = uperp(ni)**2 + upara(mi)**2
                     
                        f_cql = exp(-u2 * u0**2) * fnorm
                        dfduper(ni, mi) = -f_cql * 2.* uperp(ni) * u0**2
                        dfdupar(ni, mi) = -f_cql * 2.* upara(mi) * u0**2
 5000                   continue                        
                                                
                        
                     end do
                  end do
               end if

            end if
        
!           ----------------------------------
!           Loop over perpendicular velocities
!           ----------------------------------

            do ni = 1, nuper
            
            if (ndist .eq. 0)
     &          call QLSUM_MAXWELLIAN(ni, b_sum, c_sum, e_sum, f_sum,
     &               wdot_sum(ni), sum_fx0(ni), sum_fy0(ni), W, ZSPEC, 
     &               ASPEC, BMAG, lmax, ENORM, UminPara, UmaxPara,
     &               NUPAR, NUPER, UPERP, UPARA, DFDUPER, DFDUPAR,
     &               exk, eyk, ezk, nkdim1, nkdim2, mkdim1, mkdim2,
     &               nkx1, nkx2, nky1, nky2,
     &               uxx(i,j), uxy(i,j), uxz(i,j),
     &               uyx(i,j), uyy(i,j), uyz(i,j),
     &               uzx(i,j), uzy(i,j), uzz(i,j),
     &               nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i, j,
     &               lmaxdim, ndist, nzeta_wdot,
     &               gradprlb(i,j), bmod(i,j), omgc(i,j), alpha, xm,
     &               upshift, xk_cutoff, rt, nphi, rho(i,j))
     
            if (ndist .eq. 1)
     &         call QLSUM_NON_MAXWELLIAN(ni, b_sum, c_sum, e_sum, f_sum,
     &               wdot_sum(ni), sum_fx0(ni), sum_fy0(ni), W, ZSPEC, 
     &               ASPEC, BMAG, lmax, ENORM, UminPara, UmaxPara,
     &               NUPAR, NUPER, UPERP, UPARA, DFDUPER, DFDUPAR,
     &               exk, eyk, ezk, nkdim1, nkdim2, mkdim1, mkdim2,
     &               nkx1, nkx2, nky1, nky2,
     &               uxx(i,j), uxy(i,j), uxz(i,j),
     &               uyx(i,j), uyy(i,j), uyz(i,j),
     &               uzx(i,j), uzy(i,j), uzz(i,j),
     &               nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i, j,
     &               lmaxdim, ndist, nzeta_wdot,
     &               gradprlb(i,j), bmod(i,j), omgc(i,j), alpha, xm,
     &               upshift, xk_cutoff, rt, nphi, rho(i,j))    
                                                               
               
!              -----------------------------
!              Loop over parallel velocities
!              -----------------------------

               upara_test = sqrt(1.0 - (UPERP(ni))**2)
               mi_max = ceiling(upara_test*(nupar-1)/2 + (nupar+1)/2)
               mi_min = floor(-1.0*upara_test*(nupar-1)/2
     &                        + (nupar+1)/2)

               do mi = 1, nupar
               
                  bql = 0.0 
                  cql = 0.0
                  eql = 0.0
                  fql = 0.0

                  if(mi.ge.mi_min .and. mi.le.mi_max)then

                  bql = 1.0 / (8. * emax * dupara)
     &                  * eps0 * omgp2(i,j) / omgrf * real(b_sum(mi))

                  cql = 1.0 / (8. * emax * dupara)
     &                  * eps0 * omgp2(i,j) / omgrf * real(c_sum(mi))

                  eql = 1.0 / (8. * emax * dupara)
     &                  * eps0 * omgp2(i,j) / omgrf * real(e_sum(mi))

                  fql = 1.0 / (8. * emax * dupara)
     &                  * eps0 * omgp2(i,j) / omgrf * real(f_sum(mi))



!                  if(bql .ne. 0.0)count(myid, 1) = count(myid, 1) + 1.0
                     
                     
                  if(n .le. nnoderho)then
                     
*                    ----------------------
*                    calculate midplane u's
*                    ----------------------
                     argd = uperp(ni)**2 * (1.-bratio) + upara(mi)**2
                     if (argd .le. 0.0) argd = 1.0e-06

                     uperp0 = uperp(ni) * sqrt(bratio)
                     upara0 = sign(1.0, upara(mi)) * sqrt(argd)
                         
                     u_  = sqrt(uperp(ni)**2 + upara(mi)**2)
                     if (u_  .eq. 0.0) u_  = 1.0e-08
                     u_0 = u_
                         

                     ni0 = int((uperp0 - uperp(1)) / duperp) + 1
                     mi0 = int((upara0 - upara(1)) / dupara) + 1

*                    --------------------------------------
*                    bounce average and map to midplane u's
*                    --------------------------------------
                     if(ni0 .ge. 1 .and. ni0 .le. nuper .and. 
     &                     mi0 .ge. 1 .and. mi0 .le. nupar)then
     
     
     
                        costh0 = upara0 / u_0
                        costh = upara(mi) / u_
                        
                        if(costh0 .eq. 0.0)costh0 = 1.0e-08
                        
                        upara_mi = upara(mi)
                        if(upara_mi .eq. 0.0)
     &                        upara_mi = (upara(mi) + upara(mi+1)) / 2.0 
        

c                        factor = abs(psic * upara0 / upara_mi)
c                        factor = abs(sqrt(psic))
                        factor = 1.0
 
                        
c                        factvol(ni0, mi0, n) = factvol(ni0, mi0, n)
c     &                      + dx * dy * capr(i) / r0 * factor                  

                        bqlvol(ni0, mi0, n) = bqlvol(ni0, mi0, n)
     &                      + dx * dy * capr(i) / r0 * bql * factor 

                        cqlvol(ni0, mi0, n) = cqlvol(ni0, mi0, n)
     &                      + dx * dy * capr(i) / r0 * cql * factor
     &                      * costh / costh0 / sqrt(psic)  

                        eqlvol(ni0, mi0, n) = eqlvol(ni0, mi0, n)
     &                      + dx * dy * capr(i) / r0 * eql * factor 
     &                      * costh / costh0 / psic

                        fqlvol(ni0, mi0, n) = fqlvol(ni0, mi0, n)
     &                      + dx * dy * capr(i) / r0 * fql * factor 
     &                      * (costh / costh0)**2 / psic**1.5

                     end if
                        


                  end if

                  else
                  endif

c                  if(ndist.eq.1)then
                  
c                    bql_store(ip,ni,mi) = bql
c                   cql_store(ip,ni,mi) = cql
c                   eql_store(ip,ni,mi) = eql
c                   fql_store(ip,ni,mi) = fql
c                 endif

   
               end do

            end do

               
            wdoti = 0.0
            fx0i = 0.0
            fy0i = 0.0
               
*           --------------------------------------------
*           integrate wdot over perpendicular velocities
*           --------------------------------------------                       

            do ni = 1, nuper - 1
               wdoti = wdoti + 0.5 * duperp * 
     &            (wdot_sum(ni) + wdot_sum(ni + 1) )
               fx0i = fx0i + 0.5 * duperp * 
     &            (sum_fx0(ni)  + sum_fx0(ni + 1)  )
               fy0i = fy0i + 0.5 * duperp * 
     &            (sum_fy0(ni)  + sum_fy0(ni + 1)  )        
            end do     
              
              
     
            wdot(i,j) = - pi / 2.0 * eps0 * omgp2(i,j)
     &                                     / omgrf * real(wdoti)
               
            fx0(i,j)  = - pi / 2.0 * eps0 * omgp2(i,j)
     &                      / omgrf * real(fx0i) / (2.0 * omgrf)
            fy0(i,j)  = - pi / 2.0 * eps0 * omgp2(i,j)
     &                       / omgrf * real(fy0i)/ (2.0 * omgrf)
                    
            fz0(i,j)  = xkphi / omgrf * wdot(i,j)

        
!     --------------------------------
!     end loop over i,j spatial points
!     --------------------------------
      end do


c      if(ndist.eq.1 .and. i_write .ne. 0)then
      
c         iam_root = (myid .eq. 0)
c         left_neighbor = mod( myid - 1 + nproc, nproc)
c         right_neighbor = mod( myid + 1, nproc )
c         token = 1.0

c         if(iam_root)then
         
c           open(unit=43,file='out_orbitrf.coef',  status='replace',
c     &                                             form='formatted')
         
c          do ip=start,finish
c             i=i_table(ip)
c             j=j_table(ip)
c             do ni=1,nuper
c               do mi=1,nupar
c                 if(abs(bql_store(ip,ni,mi)) .gt. 10**(-10)) then

c                    write(43,102) i, j, ni, mi,
c     &                 bql_store(ip,ni,mi), cql_store(ip,ni,mi),
c     &                 eql_store(ip,ni,mi), fql_store(ip,ni,mi)
c                 endif
c               enddo
c             enddo
c           enddo
           
c          close(43)
           
c          call MPI_SEND(token,1, MPI_REAL,right_neighbor,
c     &             2, MPI_COMM_WORLD, ierr)
c           call MPI_RECV(token,1, MPI_REAL,left_neighbor,
c     &             2, MPI_COMM_WORLD, status, ierr)
c         else
         
c           call MPI_RECV(token,1, MPI_REAL,left_neighbor,
c     &             2, MPI_COMM_WORLD, status, ierr)
           
c           open(unit=43,file='out_orbitrf.coef',  status='old',
c     &                    form='formatted', position='append')
         
c          do ip=start,finish
c             i=i_table(ip)
c             j=j_table(ip)
c             do ni=1,nuper
c               do mi=1,nupar
c                 if(abs(bql_store(ip,ni,mi)) .gt. 10**(-10)) then

c                    write(43,102) i, j, ni, mi,
c     &                 bql_store(ip,ni,mi), cql_store(ip,ni,mi),
c     &                 eql_store(ip,ni,mi), fql_store(ip,ni,mi)
c                 endif
c               enddo
c             enddo
c           enddo
           
c          close(43)
           
c          call MPI_SEND(token,1, MPI_REAL,right_neighbor,
c     &             2, MPI_COMM_WORLD, ierr)
       
c         endif  
c      endif

      call blacs_barrier(icontxt, 'All')

!     -------------------
!     Sum over processors
!     -------------------

!      call dgsum2d(icontxt, 'All', ' ', 5001, 1, count,
!     &      5001, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex, nnodey, wdot,
     &      nnodex, -1, -1)
     
      call dgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fx0,
     &      nnodex, -1, -1)
     
      call dgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fy0,
     &      nnodex, -1, -1)
     
      call dgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fz0,
     &      nnodex, -1, -1)
     
      do n = 1, nnoderho

         do mi = 1, nupar
            do ni = 1, nuper
c               factvol2d(ni, mi) = factvol(ni, mi, n)
               bqlvol2d(ni, mi) = bqlvol(ni, mi, n)
               cqlvol2d(ni, mi) = cqlvol(ni, mi, n)
               eqlvol2d(ni, mi) = eqlvol(ni, mi, n)
               fqlvol2d(ni, mi) = fqlvol(ni, mi, n)
            end do
         end do

c         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, factvol2d,
c     &      nuper, -1, -1)
         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, bqlvol2d,
     &      nuper, -1, -1)
         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, cqlvol2d,
     &      nuper, -1, -1)
         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, eqlvol2d,
     &      nuper, -1, -1)
         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, fqlvol2d,
     &      nuper, -1, -1)
     



         do mi = 1, nupar
            do ni = 1, nuper
c              factvol(ni, mi, n) = factvol2d(ni, mi)
               bqlvol(ni, mi, n) = bqlvol2d(ni, mi)
               cqlvol(ni, mi, n) = cqlvol2d(ni, mi)
               eqlvol(ni, mi, n) = eqlvol2d(ni, mi)
               fqlvol(ni, mi, n) = fqlvol2d(ni, mi)
            end do
         end do


      end do
      
      
*     ------------------------------------
*     Divide by volume element (passed in)
*     ------------------------------------
      do n = 1, nnoderho
         do mi = 1, nupar
            do ni = 1, nuper

*           --------------
*           bounce average
*           --------------
            if (vol(n) .ne. 0.0)then
            bqlavg(ni, mi, n) = bqlvol(ni, mi, n) / vol(n) * dldbavg(n) 
            cqlavg(ni, mi, n) = cqlvol(ni, mi, n) / vol(n) * dldbavg(n) 
            eqlavg(ni, mi, n) = eqlvol(ni, mi, n) / vol(n) * dldbavg(n) 
            fqlavg(ni, mi, n) = fqlvol(ni, mi, n) / vol(n) * dldbavg(n)
            end if
            
            if (bqlavg(ni, mi, n) .lt. 0.0)bqlavg(ni, mi, n) = 0.0
c           if (cqlavg(ni, mi, n) .lt. 0.0)cqlavg(ni, mi, n) = 0.0
c           if (eqlavg(ni, mi, n) .lt. 0.0)eqlavg(ni, mi, n) = 0.0
            if (fqlavg(ni, mi, n) .lt. 0.0)fqlavg(ni, mi, n) = 0.0          
            
            end do
         end do
      end do
      
      deallocate( dfduper0 )
      deallocate( dfdupar0 )


      deallocate(b_sum)
      deallocate(c_sum)
      deallocate(e_sum)
      deallocate(f_sum)
      
      deallocate(wdot_sum)
      deallocate(sum_fx0)
      deallocate(sum_fy0)      
      
c      deallocate(factvol)
c      deallocate(factvol2d)      
      
      deallocate(bqlvol)
      deallocate(bqlvol2d)

      deallocate(cqlvol)
      deallocate(cqlvol2d)

      deallocate(eqlvol)
      deallocate(eqlvol2d)

      deallocate(fqlvol)
      deallocate(fqlvol2d)
      
c      if(ndist.eq.1) then
c         deallocate(bql_store, cql_store, eql_store,  fql_store)
c      endif

      wdot_inout(1:nnodex,1:nnodey) = wdot(1:nnodex,1:nnodey)
      fx0_inout(1:nnodex,1:nnodey) = fx0(1:nnodex,1:nnodey)      
      fy0_inout(1:nnodex,1:nnodey) = fy0(1:nnodex,1:nnodey)      
      fz0_inout(1:nnodex,1:nnodey) = fz0(1:nnodex,1:nnodey)
            
!      if (myid .eq. 0)then
!         sum_count = 0.0
!         do n = 0, 5000
!            write(43, 101)n, count(n, 1)
!           sum_count = sum_count + count(n, 1)
!         end do
!        write(43, *)"sum_count = ", sum_count
!      end if

      return

 1311 format(1p,9e12.4)
  100 format (1p,8e12.4)
  101 format (1i6, 1p,8e12.4)
  102 format (4i6, 1p,8e12.4)

      end subroutine ql_myra_write

      
c
c***************************************************************************
c



      subroutine wdot_qlcheck(wdot_check, 
     &   nnoderho, nrhodim,
     &   bqlavg, cqlavg, xm, omgrf, xktavg, ndist,
     &   nupar, nuper, n_psi,
     &   n_psi_dim, dfduper, dfdupar, fperp, 
     &   UminPara_cql, UmaxPara_cql, UPERP_cql, UPARA_cql, UPERP, UPARA,
     &   vc_mks_cql, df_cql_uprp, df_cql_uprl, rhon, rho_a, myid,
     &   dldbavg)

*     ------------------------------------------------------------------
*     This subroutine calculates the flux averaged wdot for checking QL
*     ------------------------------------------------------------------

      implicit none


      integer nnoderho, nrhodim, n, m, ndist, myid
      integer  :: n_psi_dim, nuper, nupar, n_psi, mi0, ni0

      real bqlavg(nuper, nupar, nnoderho)
      real cqlavg(nuper, nupar, nnoderho)
      real xktavg(nnoderho), alpha, e
      real ans, dfdth, dfdu, u, xm, omgrf, jacobian, upara_mi
      real, dimension(:,:), allocatable :: wdot_int
      real wdot_check(nrhodim), rhon(nrhodim)
      real dldbavg(nrhodim)

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: UPERP_cql(NUPER), UPARA_cql(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: FPERP(NUPER)     
      
      real :: DFDUPER0, DFDUPAR0, UPARA0
      real :: W, ENORM
      real :: UminPara,UmaxPara
      real :: UminPara_cql,UmaxPara_cql

      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)


      real :: vc_mks, vc_mks_cql, rho_a(n_psi_dim)
      real :: eps0, pi, emax, u0, dfdu0, dfdth0, u_0

      parameter (eps0 = 8.85e-12)
      parameter (PI = 3.141592653597932384)

      allocate(wdot_int(nuper, nupar) )

!     ------------------------------------
!efd  initialize allocatable array to zero
!     ------------------------------------


      wdot_int = 0.0

      e = 1.6e-19
      W = omgrf

      if(ndist .eq. 0)then   !--Maxwellian--!

         do n = 1, NUPER
            UPERP(n) = (real(n-1)/real(NUPER-1))
         end do

         UminPara = -1.0
         UmaxPara =  1.0

         do m = 1, NUPAR
            UPARA(m) = (-1.0 + 2. * (real(m-1) / real(NUPAR-1)))
         end do

      else   !--non-Maxwellian--!

         vc_mks = vc_mks_cql
         UminPara = UminPara_cql
         UmaxPara = UmaxPara_cql

         do n = 1, nuper
            uperp(n) = uperp_cql(n)
         end do

         do m = 1, nupar
            upara(m) = upara_cql(m)
         end do

      end if


*     -------------------
*     Loop over rho mesh:
*     -------------------

      do n = 1, nnoderho

         alpha = sqrt(2.0 * xktavg(n) / xm)
         if (ndist .eq. 0) vc_mks =  3.0 * alpha    !--Maxwellian only--!
         u0 = vc_mks / alpha

         Emax = 0.5 * xm * vc_mks**2
         Enorm = Emax / 1.6e-19

!        ------------------------------------------------
!        get CQL3D distribution function on the midplane
!        ------------------------------------------------

         if(ndist .eq. 0)then   !--Maxwellian--!
        
            call maxwell_dist(u0, NUPAR, NUPER,
     &                 UminPara, UmaxPara,
     &                 UPERP, UPARA, DFDUPER, DFDUPAR, FPERP)

         else   !--non-Maxwellian--!
        
            call cql3d_dist(nupar, nuper, n_psi,
     &                 n_psi_dim, rho_a, rhon(n),
     &                 UminPara,UmaxPara,
     &                 df_cql_uprp, df_cql_uprl,
     &                 UPERP, UPARA, DFDUPER, DFDUPAR)



         end if
         
!        -----------------------------
!        loop over MIDPLANE velocities
!        -----------------------------   

         do ni0 = 1, nuper
            do mi0 = 1, nupar

               dfdupar0 = dfdupar(ni0, mi0)
               dfduper0 = dfduper(ni0, mi0)

               u_0 = sqrt(uperp(ni0)**2 + upara(mi0)**2)
               if (u_0 .eq. 0.0) u_0 = 1.0e-08
                
               dfdu0 = (uperp(ni0) * dfduper0 + upara(mi0) * dfdupar0)
     &                / u_0
               dfdth0 = upara(mi0) * dfduper0 - uperp(ni0) * dfdupar0
               
               wdot_int(ni0, mi0) = (bqlavg(ni0, mi0, n) * dfdu0
     &                             + cqlavg(ni0, mi0, n) * dfdth0) / u_0 

            end do

         end do


        
!        ---------------------------------------------------
!        Do velocity space integral over midplane velocities
!        ---------------------------------------------------

         wdot_check(n) = 0.0


         call ugrate(wdot_int, uperp, upara, nuper, nupar, ans, myid,xm)


         wdot_check(n) = - 4.0 * pi * e * enorm / dldbavg(n) * ans


      end do


      deallocate(wdot_int)


      return

 1311 format(1p,9e12.4)
 1312 format(i10, 1p,9e12.4)
 1313 format(2i10, 1p,9e12.4)
  100 format (1p,8e12.4)
  101 format (2i10, 1p,8e12.4)

      end subroutine 
      

c
c***************************************************************************
c



      subroutine ugrate(f, uperp, upara, nuper, nupar, fint, myid, xm)

      implicit none

      integer nuper, nupar, ni, mi, myid

      real uperp(nuper), upara(nupar), f(nuper, nupar)
      real favg, fint, duperp, dupara, xm

      duperp = (uperp(nuper) - uperp(1)) / (nuper - 1)
      dupara = (upara(nupar) - upara(1)) / (nupar - 1)

      fint = 0.0

      do ni = 1, nuper - 1
         do mi = 1, nupar - 1
            favg = (f(ni, mi)   + f(ni+1, mi)
     &            + f(ni, mi+1) + f(ni+1, mi+1)) / 4.0

            fint = fint + favg * uperp(ni) * duperp * dupara

c            if (myid.eq.0 .and. ni .eq. 32 .and. mi .eq. 95) then
c                write(6 ,1313)ni, mi, xm, favg, fint
c                write(15,1313)ni, mi, xm, favg, fint
c           end if

         end do
      end do

 1313 format(2i10, 1p,9e12.4)

      return
      end subroutine  ugrate

c
c***************************************************************************
c

       end module ql_myra_mod
