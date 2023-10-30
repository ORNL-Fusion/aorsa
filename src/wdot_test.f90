      module wdot_mod
      use qlsum_myra_mod
      use wdot_sum_mod
      !      include 'mpif.h'
      use mpi
      contains

!
!***************************************************************************
!

      subroutine wdot_new_maxwellian(z2_electron,       &
     &   wdot_inout,                                    &
     &   vol, nrhodim, nnoderho, drho,                  &
     &   dx, dy, r0, nxdim, nydim, nnodex, nnodey,                   &
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,    &
     &   xm, q, xn, xkt, omgc, omgp2, lmax,    &
     &   xkxsav, xkysav, nzfun, ibessel,    &
     &   exk, eyk, ezk, nphi, capr,    &
     &   bxn, byn, bzn,    &
     &   uxx, uxy, uxz,    &
     &   uyx, uyy, uyz,    &
     &   uzx, uzy, uzz,    &
     &   dxuzx, dyuzx, dxuzy, dyuzy, dxuzz, dyuzz, drdx,     &
     &   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,    &
     &   xx, yy, isigma, xnuomga, psi, psilim, nboundary,    &
     &   myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid,    &
     &   nupar, nuper, n_psi,    &
     &   n_psi_dim, dfduper, dfdupar, fperp,    &
     &   UminPara_cql, UmaxPara_cql, UPERP_cql, UPARA_cql, UPERP, UPARA,    &
     &   vc_mks_cql, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj,    &
     &   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,    &
     &   nkperp, lmaxdim, nzeta_wdot, theta_,    &
     &   n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg,     &
     &   n_bin, upshift, i_write, xk_cutoff, z0_table, z1_table,     &
     &   z2_table, ntable, mtable, zetai_table, dKdL_table,     &
     &   nmax, mmax)

!----------------------------------------------------------------------
!     This subroutine calculates wdot and the quasi-linear operator for 
!     a single species
!----------------------------------------------------------------------

      implicit none

      integer z2_electron
      integer splitrank, nloops, partition, start, finish
      integer pstart, pfinish, status(MPI_STATUS_SIZE), ierr
      integer recv_size, remainder

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
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,    &
     &   icontxt, nboundary, nzeta_wdot
      integer dlen_, desc_amat(dlen_)
      integer nxdim, nydim, nnodex, nnodey, nrhodim, nnoderho,    &
     &   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma
      integer ftrap, k, ni0, mi0
      
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt, xk_cutoff
      
      real theta_(n_theta_max, n_psi_max), thegiv, dthetag
      real deriv, dtheta, dtau_ratio_giv, tau_bounce_giv
      real uperp0, upara0, xlamda, derivb, dtheta0, factor, upara_mi
      real duperp, dupara

      real zeffcd, clight, rt, b0, damp, r1, xm1, ceta, epsa, xmut2,    &
     &   eta0, xnexp, xjtild, signkz, cfit, c1, afit, yt, bmaxa, vphase,    &
     &   omgrf, c_ehst, akprl, xkprl, xkphi, a, rmaxa, rmina,    &
     &   wphase, xlnlam, vth
      real x(nxdim), y(nydim), argd, bratio, drho


      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim),    &
     &     psi(nxdim, nydim), psilim, alpha    
      real capr(nxdim), xnuomg, delta0, dx, dy, r0
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim),    &
     &   bmod_mid(nxdim, nydim), xnuomga(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim),    &
     &        yy(mkdim1 : mkdim2, 1 : nydim)
     
      complex wdoti, fx0i, fy0i

      
      real wdot_inout(nxdim, nydim)

      
      real wdot(nnodex, nnodey)

      
!      real count(0 : 10000, 1), sum_count

      real vol(nrhodim), dldbavg(nrhodim)
        


      complex sigxx, sigxy, sigxz,    &
     &        sigyx, sigyy, sigyz,    &
     &        sigzx, sigzy, sigzz
     
      integer ntable, mtable, nmax, mmax
      
      complex z0_table(ntable, mtable)
      complex z1_table(ntable, mtable)
      complex z2_table(ntable, mtable)
      real zetai_table(ntable), dKdL_table(mtable)              

      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2),    &
     &        eyk(nkdim1 : nkdim2, mkdim1 : mkdim2),    &
     &        ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real q, xn(nxdim, nydim), xkt(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim),    &
     &     uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim),    &
     &     uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)
     
      real dxuzx(nxdim, nydim), dyuzx(nxdim, nydim), dxuzy(nxdim,nydim),   &
     &     dyuzy(nxdim, nydim), dxuzz(nxdim, nydim), dyuzz(nxdim,nydim) 
     
      real drdx(nxdim)
     
      real dbxdx, dbxdy, dbydx, dbydy, dbphidx, dbphidy, dbphirdx
      real bhatgradbx, bhatgradby, bhatgradbz  

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi, mi, ni

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: UPERP_cql(NUPER), UPARA_cql(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: FPERP(NUPER)


      complex, dimension(:),  allocatable :: b_sum
      complex, dimension(:),  allocatable :: c_sum
      complex, dimension(:),  allocatable :: e_sum
      complex, dimension(:),  allocatable :: f_sum
      
      complex, dimension(:),  allocatable :: wdot_sum
      complex, dimension(:),  allocatable :: sum_fx0
      complex, dimension(:),  allocatable :: sum_fy0
      

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

!     -------------------------------------
!efd  initialize allocatable arrays to zero
!     -------------------------------------
      wdot = 0.0


      b_sum = 0.0
      c_sum = 0.0
      e_sum = 0.0
      f_sum = 0.0

      wdot_sum = 0.0
!      sum_fx0 = 0.0
!      sum_fy0 = 0.0


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



      nwork = 0
      do j=1,nnodey
         do i=1,nnodex
           has_work = (psi(i,j) .le. psilim .and. nboundary .eq. 1    &
     &                                       .or. nboundary .eq. 0)
           if (has_work) then
              nwork = nwork + 1
              i_table(nwork) = i
              j_table(nwork) = j
           endif
         enddo
      enddo
      

!     -----------------------
!     Loop over spatial mesh:
!     -----------------------

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
      

!     ----------------------------
!     loop over i,j spatial points
!     ----------------------------
      do ip = start,finish
      
            i = i_table(ip)
            j = j_table(ip)
            
            xnuomg = xnuomga(i,j)
            
            dbxdx = dxuzx(i,j)
            dbxdy = dyuzx(i,j)
                     
            dbydx = dxuzy(i,j)
            dbydy = dyuzy(i,j)
                     
            dbphidx = dxuzz(i,j)
            dbphidy = dyuzz(i,j)            
            
            bhatgradbx =  bxn(i,j) * dbxdx     &
     &                  + byn(i,j) * dbxdy   
            bhatgradby =  bxn(i,j) * dbydx     &
     &                  + byn(i,j) * dbydy 
            bhatgradbz =  bxn(i,j) * (dbphidx      &
     &                  - bzn(i,j) / capr(i) * drdx(i))    &
     &                  + byn(i,j) * dbphidy        

        
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
        
               call maxwell_dist(u0, NUPAR, NUPER,    &
     &              UminPara, UmaxPara,    &
     &              UPERP, UPARA, DFDUPER, DFDUPAR, FPERP)

            else   !--non-Maxwellian--!
        
               call cql3d_dist(nupar, nuper, n_psi,    &
     &              n_psi_dim, rho_a, rho(i,j),    &
     &              UminPara,UmaxPara,    &
     &              df_cql_uprp, df_cql_uprl,    &
     &              UPERP, UPARA, DFDUPER0, DFDUPAR0)

!              ---------------------------------------------------------
!              map CQL3D distribution function off the midplane for Wdot
!              ---------------------------------------------------------
               if(bratio .gt. 0.0)then
               
                  dfduper = 0.0
                  dfdupar = 0.0
                       
                  do ni = 1, nuper
                     do mi = 1, nupar

                        argd =  uperp(ni)**2 * (1. - bratio)    &
     &                                                   + upara(mi)**2
                        if (argd .le. 0.0) argd = 1.0e-06
                                                
                        uperp0 = uperp(ni) * sqrt(bratio)
                        upara0  = sign(1.0, upara(mi)) * sqrt(argd)
                        
                        dfduper(ni, mi) = 0.0
                        dfdupar(ni, mi) = 0.0
                                
                        if(upara0 .ge. upara(1) .and.     &
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
                           di = dfduper0(ni0+1,mi0+1)+ dfduper0(ni0,mi0)     &
     &                        - dfduper0(ni0+1,mi0)- dfduper0(ni0,mi0+1) 
                           
                           dfduper0_intplt = ai + bi * zeta     &
     &                                     + ci * eta + di * zeta * eta                         

                           ai = dfdupar0(ni0, mi0)
                           bi = dfdupar0(ni0+1,mi0) - dfdupar0(ni0,mi0)
                           ci = dfdupar0(ni0,mi0+1) - dfdupar0(ni0,mi0)
                           di = dfdupar0(ni0+1,mi0+1)+ dfdupar0(ni0,mi0)     &
     &                        - dfdupar0(ni0+1,mi0)- dfdupar0(ni0,mi0+1) 
                           
                           dfdupar0_intplt = ai + bi * zeta     &
     &                                     + ci * eta + di * zeta * eta 
     
                           end if                       
                                                        
                                        
                           if (upara0 .ne. 0.0)then
                        
                              dfdupar(ni, mi) = dfdupar0_intplt *     &
     &                           upara(mi) / upara0
     
                              dfduper(ni, mi) = dfduper0_intplt *     &
     &                           sqrt(bratio) + dfdupar0_intplt *     &
     &                           uperp(ni) / upara0 * (1.0 - bratio)
     
!                             dfdth = upara(mi) * dfduper(ni, mi)    &
!     &                             - uperp(ni) * dfdupar(ni, mi)
     
                           end if
                           
                        end if
                        
                                                                                                                        
                     end do
                  end do
               end if

            end if
        
!           ----------------------------------
!           Loop over perpendicular velocities
!           ----------------------------------

            do ni = 1, nuper
            
            if (ndist .eq. 0)    &
     &          call wdot_new_maxwellian_sum(i,j, ni, b_sum, c_sum,    & 
     &               e_sum, f_sum,    &
     &               wdot_sum(ni), sum_fx0(ni), sum_fy0(ni), W, ZSPEC,     &
     &               ASPEC, BMAG, lmax, ENORM, UminPara, UmaxPara,    &
     &               NUPAR, NUPER, UPERP, UPARA, DFDUPER, DFDUPAR,FPERP,    &
     &               exk, eyk, ezk, nkdim1, nkdim2, mkdim1, mkdim2,    &
     &               nkx1, nkx2, nky1, nky2,    &
     &               uxx(i,j), uxy(i,j), uxz(i,j),    &
     &               uyx(i,j), uyy(i,j), uyz(i,j),    &
     &               uzx(i,j), uzy(i,j), uzz(i,j),    &
     &               nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i, j,    &
     &               lmaxdim, ndist, nzeta_wdot,    &
     &               gradprlb(i,j), bmod(i,j), omgc(i,j), alpha, xm,    &
     &               upshift, xk_cutoff, rt, nphi, rho(i,j),     &
     &               z0_table, z1_table, z2_table, ntable, mtable,    &
     &               zetai_table, dKdL_table,  nmax, mmax,    &
     &               bhatgradbx, bhatgradby, bhatgradbz, z2_electron,    &
     &               xnuomg)
     
     
            if (ndist .eq. 1)    &
     &         call QLSUM_NON_MAXWELLIAN(ni, b_sum, c_sum, e_sum, f_sum,    &
     &               wdot_sum(ni), sum_fx0(ni), sum_fy0(ni), W, ZSPEC,     &
     &               ASPEC, BMAG, lmax, ENORM, UminPara, UmaxPara,    &
     &               NUPAR, NUPER, UPERP, UPARA, DFDUPER, DFDUPAR,    &
     &               exk, eyk, ezk, nkdim1, nkdim2, mkdim1, mkdim2,    &
     &               nkx1, nkx2, nky1, nky2,    &
     &               uxx(i,j), uxy(i,j), uxz(i,j),    &
     &               uyx(i,j), uyy(i,j), uyz(i,j),    &
     &               uzx(i,j), uzy(i,j), uzz(i,j),    &
     &               nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i, j,    &
     &               lmaxdim, ndist, nzeta_wdot,    &
     &               gradprlb(i,j), bmod(i,j), omgc(i,j), alpha, xm,    &
     &               upshift, xk_cutoff, rt, nphi, rho(i,j))    
                                                              

            end do   ! end of loop over uperp

               

               
!           --------------------------------------------
!           integrate wdot over perpendicular velocities
!           --------------------------------------------
            wdoti = 0.0
            do ni = 1, nuper - 1
               wdoti = wdoti + 0.5 * duperp *     &
     &            (wdot_sum(ni) + wdot_sum(ni + 1) )
            end do
            
                                                              
            wdot(i,j) = pi**1.5 * omgp2(i,j) * eps0 / omgrf     &
     &                                          * aimag(wdoti)
            
              
        
!     --------------------------------
!     end loop over i,j spatial points
!     --------------------------------
      end do



      call blacs_barrier(icontxt, 'All')

!     -------------------
!     Sum over processors
!     -------------------

      call dgsum2d(icontxt, 'All', ' ', nnodex, nnodey, wdot,    &
     &      nnodex, -1, -1)
     

      
      deallocate( dfduper0 )
      deallocate( dfdupar0 )

      deallocate(b_sum)
      deallocate(c_sum)
      deallocate(e_sum)
      deallocate(f_sum)
      
      deallocate(wdot_sum)
      deallocate(sum_fx0)
      deallocate(sum_fy0)
                    
      wdot_inout(1:nnodex,1:nnodey) = wdot(1:nnodex,1:nnodey)

            

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (1i6, 1p8e12.4)
  102 format (4i6, 1p8e12.4)

      end subroutine wdot_new_maxwellian

      
!    
!***************************************************************************
!

     end module wdot_mod

