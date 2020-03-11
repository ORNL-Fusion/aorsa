c
c***************************************************************************
c

      subroutine odd_order_derivs(i, j, n, m, nxdim, nydim, 
     .   nnodex, nnodey, mkdim1, mkdim2, capr,
     .   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, 
     .   xm, q, xn, xkt, omgc, omgp2, lmax,
     .   xkxsav, xkysav, xkphi, nzfun, ibessel, nphi, 
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz,    
     .   myrow, mycol, nprow, npcol, icontxt, 
     .   desc_amat, dlen_, nboundary,
     .   xx, yy, isigma, xnuomg, psi, psilim, 
     .   myid, nproc, delta0, gradprlb, bmod, 
     .   bmod_mid, ndist, nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, rho, 
     .   rho_a, nbessj, nkperp,
     .   zi, eps0, v0i, omgrf, xk0, kperp_max, 
     .   i_sav, j_sav, upshift, 
     .   damping, xk_cutoff, rt, dx, dy)

*-----------------------------------------------------------------
*     This subroutine adds the odd order derivative terms to sigma
*-----------------------------------------------------------------

      implicit none

      logical ismine

      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max, dx, dy, q
      integer i_sav, j_sav, upshift 
      real damping, xk_cutoff, rt

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt, nboundary
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim),
     .     psi(nxdim, nydim), psilim
      real capr(nxdim), xnuomg, delta0
      real xkphi(nxdim)
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim),
     .                         bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim)

      complex sigxx, sigxy, sigxz,
     .        sigyx, sigyy, sigyz,
     .        sigzx, sigzy, sigzz
     
      complex delta_x, delta_y, delta_z
      
      complex cexpkxky

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real  xn(nxdim, nydim), xkt(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim),
     .     uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim),
     .     uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)
     
      complex dkwxxp, dkwxyp, dkwxzp,
     .        dkwyxp, dkwyyp, dkwyzp, 
     .        dkwzxp, dkwzyp, dkwzzp
     
      complex dkwxxm, dkwxym, dkwxzm,
     .        dkwyxm, dkwyym, dkwyzm, 
     .        dkwzxm, dkwzym, dkwzzm                 
     
      complex dxdkwxx, dxdkwxy, dxdkwxz,
     .        dxdkwyx, dxdkwyy, dxdkwyz,
     .        dxdkwzx, dxdkwzy, dxdkwzz 
                   
      complex dydkwxx, dydkwxy, dydkwxz,
     .        dydkwyx, dydkwyy, dydkwyz,
     .        dydkwzx, dydkwzy, dydkwzz 
          
      complex wxx, wxy, wxz,
     .        wyx, wyy, wyz,
     .        wzx, wzy, wzz
     
      complex wxxp, wxyp, wxzp,
     .        wyxp, wyyp, wyzp,
     .        wzxp, wzyp, wzzp              

      real xkalp, xkbet, xkprl, xkperp2, xkperp, cosa, sina, 
     .   cosa2, sina2

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)
      
      xkalp = uxx(i, j) * xkxsav(n) 
     .      + uxy(i, j) * xkysav(m) 
     .      + uxz(i, j) * xkphi(i)    
      xkbet = uyx(i, j) * xkxsav(n) 
     .      + uyy(i, j) * xkysav(m) 
     .      + uyz(i, j) * xkphi(i)
      xkprl = uzx(i, j) * xkxsav(n) 
     .      + uzy(i, j) * xkysav(m) 
     .      + uzz(i, j) * xkphi(i)
	   		  
      xkperp2 = xkalp**2 + xkbet**2
      xkperp = sqrt(xkperp2) 
        
      cosa = xkalp / xkperp
      sina = xkbet / xkperp 
      
      cosa2 = cosa**2
      sina2 = sina**2 
                  

!     ------------------------------
!     take kperp derivatives in +-x
!     ------------------------------     
      if (i .ne. 1 .and. i .ne. nnodex)then
                             
         call dkw(i+1, j, n, m, rho(i+1,j), rho_a,
     .      gradprlb(i+1,j), bmod(i+1,j), bmod_mid(i+1,j),
     .      xm, q, xn(i+1,j), xnuomg,
     .      xkt(i+1,j), omgc(i+1,j), omgp2(i+1,j),
     .      -lmax, lmax, nzfun, ibessel,
     .      xkxsav(n), xkysav(m), nphi, capr(i+1),
     .      bxn(i+1,j), byn(i+1,j), bzn(i+1,j),
     .      uxx(i+1,j), uxy(i+1,j), uxz(i+1,j),
     .      uyx(i+1,j), uyy(i+1,j), uyz(i+1,j),
     .      uzx(i+1,j), uzy(i+1,j), uzz(i+1,j),
     .      dkwxxp, dkwxyp, dkwxzp,
     .      dkwyxp, dkwyyp, dkwyzp, 
     .      dkwzxp, dkwzyp, dkwzzp,
     .      delta0, ndist, nupar, nuper, n_psi,
     .      n_psi_dim, dfduper, dfdupar,
     .      UminPara, UmaxPara, UPERP, UPARA,
     .      vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .      nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .      i_sav, j_sav, upshift, damping, xk_cutoff,rt,
     .      nkx2, nky2) 
     
         call dkw(i-1, j, n, m, rho(i-1,j), rho_a,
     .      gradprlb(i-1,j), bmod(i-1,j), bmod_mid(i-1,j),
     .      xm, q, xn(i-1,j), xnuomg,
     .      xkt(i-1,j), omgc(i-1,j), omgp2(i-1,j),
     .      -lmax, lmax, nzfun, ibessel,
     .      xkxsav(n), xkysav(m), nphi, capr(i-1),
     .      bxn(i-1,j), byn(i-1,j), bzn(i-1,j),
     .      uxx(i-1,j), uxy(i-1,j), uxz(i-1,j),
     .      uyx(i-1,j), uyy(i-1,j), uyz(i-1,j),
     .      uzx(i-1,j), uzy(i-1,j), uzz(i-1,j),
     .      dkwxxm, dkwxym, dkwxzm,
     .      dkwyxm, dkwyym, dkwyzm, 
     .      dkwzxm, dkwzym, dkwzzm,
     .      delta0, ndist, nupar, nuper, n_psi,
     .      n_psi_dim, dfduper, dfdupar,
     .      UminPara, UmaxPara, UPERP, UPARA,
     .      vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .      nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .      i_sav, j_sav, upshift, damping, xk_cutoff, rt,
     .      nkx2, nky2)
     

!        ------------------
!        take x derivatives
!        ------------------ 
         dxdkwxx = (dkwxxp - dkwxxm) / (2. * dx)
         dxdkwxy = (dkwxyp - dkwxym) / (2. * dx)
         dxdkwxz = (dkwxzp - dkwxzm) / (2. * dx)
         dxdkwyx = (dkwyxp - dkwyxm) / (2. * dx)
         dxdkwyy = (dkwyyp - dkwyym) / (2. * dx)
         dxdkwyz = (dkwyzp - dkwyzm) / (2. * dx)
         dxdkwzx = (dkwzxp - dkwzxm) / (2. * dx)
         dxdkwzy = (dkwzyp - dkwzym) / (2. * dx)
         dxdkwzz = (dkwzzp - dkwzzm) / (2. * dx)
	 
      end if
      
      
!     ------------------------------
!     take kperp derivatives in +-y
!     ------------------------------	    
      if (j .ne. 1. and. j .ne. nnodey)then
         call dkw(i, j+1, n, m, rho(i,j+1), rho_a,
     .      gradprlb(i,j+1), bmod(i,j+1), bmod_mid(i,j+1),
     .      xm, q, xn(i,j+1), xnuomg,
     .      xkt(i,j+1), omgc(i,j+1), omgp2(i,j+1),
     .      -lmax, lmax, nzfun, ibessel,
     .      xkxsav(n), xkysav(m), nphi, capr(i),
     .      bxn(i,j+1), byn(i,j+1), bzn(i,j+1),
     .      uxx(i,j+1), uxy(i,j+1), uxz(i,j+1),
     .      uyx(i,j+1), uyy(i,j+1), uyz(i,j+1),
     .      uzx(i,j+1), uzy(i,j+1), uzz(i,j+1),
     .      dkwxxp, dkwxyp, dkwxzp,
     .      dkwyxp, dkwyyp, dkwyzp, 
     .      dkwzxp, dkwzyp, dkwzzp,
     .      delta0, ndist, nupar, nuper, n_psi,
     .      n_psi_dim, dfduper, dfdupar,
     .      UminPara, UmaxPara, UPERP, UPARA,
     .      vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .      nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .      i_sav, j_sav, upshift, damping, xk_cutoff,rt,
     .      nkx2, nky2) 
       				 
			
         call dkw(i, j-1, n, m, rho(i,j-1), rho_a,
     .      gradprlb(i,j-1), bmod(i,j-1), bmod_mid(i,j-1),
     .      xm, q, xn(i,j-1), xnuomg,
     .      xkt(i,j-1), omgc(i,j-1), omgp2(i,j-1),
     .      -lmax, lmax, nzfun, ibessel,
     .      xkxsav(n), xkysav(m), nphi, capr(i),
     .      bxn(i,j-1), byn(i,j-1), bzn(i,j-1),
     .      uxx(i,j-1), uxy(i,j-1), uxz(i,j-1),
     .      uyx(i,j-1), uyy(i,j-1), uyz(i,j-1),
     .      uzx(i,j-1), uzy(i,j-1), uzz(i,j-1),
     .      dkwxxm, dkwxym, dkwxzm,
     .      dkwyxm, dkwyym, dkwyzm, 
     .      dkwzxm, dkwzym, dkwzzm,
     .      delta0, ndist, nupar, nuper, n_psi,
     .      n_psi_dim, dfduper, dfdupar,
     .      UminPara, UmaxPara, UPERP, UPARA,
     .      vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .      nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max,
     .      i_sav, j_sav, upshift, damping, xk_cutoff,rt,
     .      nkx2, nky2) 
          
			
!        ------------------
!        take y derivatives
!        ------------------
         dydkwxx = (dkwxxp - dkwxxm) / (2. * dy)
         dydkwxy = (dkwxyp - dkwxym) / (2. * dy)
         dydkwxz = (dkwxzp - dkwxzm) / (2. * dy)
         dydkwyx = (dkwyxp - dkwyxm) / (2. * dy)
         dydkwyy = (dkwyyp - dkwyym) / (2. * dy)
         dydkwyz = (dkwyzp - dkwyzm) / (2. * dy)
         dydkwzx = (dkwzxp - dkwzxm) / (2. * dy)
         dydkwzy = (dkwzyp - dkwzym) / (2. * dy)
         dydkwzz = (dkwzzp - dkwzzm) / (2. * dy)
			     			   			   			   
      end if
      
      wxx = uxx(i,j) * dxdkwxx + uxy(i,j) * dydkwxx 
      wxy = uxx(i,j) * dxdkwxy + uxy(i,j) * dydkwxy     
      wxz = uxx(i,j) * dxdkwxz + uxy(i,j) * dydkwxz 
          
      wyx = uxx(i,j) * dxdkwyx + uxy(i,j) * dydkwyx     
      wyy = uxx(i,j) * dxdkwyy + uxy(i,j) * dydkwyy      
      wyz = uxx(i,j) * dxdkwyz + uxy(i,j) * dydkwyz 
           
      wzx = uxx(i,j) * dxdkwzx + uxy(i,j) * dydkwzx     
      wzy = uxx(i,j) * dxdkwzy + uxy(i,j) * dydkwzy      
      wzz = uxx(i,j) * dxdkwzz + uxy(i,j) * dydkwzz 
                                  			

!     -------------------------------------------------------
!     rotate odd order derivatives for kbeta .ne. 0 (Swanson)
!     -------------------------------------------------------      
      wxxp = cosa2 * wxx + sina2 * wyy - cosa * sina * (wxy + wyx) 
      wxyp = cosa2 * wxy - sina2 * wyx + cosa * sina * (wxx - wyy)
      wxzp = cosa * wxz - sina * wyz	      
	           
      wyxp = cosa2 * wyx - sina2 * wxy + cosa * sina * (wxx - wyy)
      wyyp = cosa2 * wyy + sina2 * wxx + cosa * sina * (wxy + wyx)
      wyzp = wxz * sina + wyz * cosa	 	 
     
      wzxp = wzx * cosa - wzy * sina
      wzyp = wzx * sina + wzy * cosa
      wzzp = wzz

            
!     ------------------------------------------
!     add rotated odd-order derivatives to sigma
!     ------------------------------------------				   
      sigxx = sigxx - zi * wxxp 
      sigxy = sigxy - zi * wxyp      
      sigxz = sigxz - zi * wxzp 
          
      sigyx = sigyx - zi * wyxp    
      sigyy = sigyy - zi * wyyp      
      sigyz = sigyz - zi * wyzp 
           
      sigzx = sigzx - zi * wzxp     
      sigzy = sigzy - zi * wzyp     
      sigzz = sigzz - zi * wzzp
                                                    			      
      return
      end


c
c***************************************************************************
c

      subroutine ntilda_(ntilda, nxdim, nydim, nnodex, nnodey,
     .   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     .   xm, q, xn, xkt, omgc, omgp2, lmax,
     .   xkxsav, xkysav, nzfun, ibessel,
     .   exk, eyk, ezk, nphi, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     .   xx, yy, isigma, xnuomg, psi, psilim, nboundary,
     .   myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj, 
     .   nkperp,
     .   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, 
     .   damping, xk_cutoff, rt)

*-----------------------------------------------------------------------
*     This subroutine calculates the plasma current for a single species
*-----------------------------------------------------------------------

      implicit none

      logical ismine

      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max
      integer i_sav, j_sav, upshift 
      real damping, xk_cutoff, rt

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt, nboundary
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim),
     .     psi(nxdim, nydim), psilim
      real capr(nxdim), xnuomg, delta0
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim),
     .                         bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim)

      complex ntilda(nxdim, nydim)

      complex sigxx, sigxy, sigxz,
     .        sigyx, sigyy, sigyz,
     .        sigzx, sigzy, sigzz
     
      complex delta_x, delta_y, delta_z

      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     .        eyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     .        ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real q, xn(nxdim, nydim), xkt(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim),
     .     uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim),
     .     uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
	real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)
      
      ntilda(:,:) = 0.0

c--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (i - 1) * nnodey + j
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

	    ntilda(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     .                              .or.  nboundary .eq. 0)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     call delta_(i, j, n, m, rho(i,j),rho_a,
     .                      delta_x, delta_y, delta_z,
     .                      gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     .                      xm, q, xn(i,j), xnuomg,
     .                      xkt(i,j), omgc(i,j), omgp2(i,j),
     .                      -lmax, lmax, nzfun, ibessel,
     .                      xkxsav(n), xkysav(m), nphi, capr(i),
     .                      bxn(i,j), byn(i,j), bzn(i,j),
     .                      uxx(i,j), uxy(i,j), uxz(i,j),
     .                      uyx(i,j), uyy(i,j), uyz(i,j),
     .                      uzx(i,j), uzy(i,j), uzz(i,j),
     .                      delta0, ndist, nupar, nuper, n_psi,
     .                      n_psi_dim, dfduper, dfdupar,
     .                      UminPara, UmaxPara, UPERP, UPARA,
     .                      vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .                      nkperp, zi, eps0, v0i, omgrf, xk0,
     .                      kperp_max, i_sav, j_sav, upshift, 
     .                      damping, xk_cutoff, rt)


                     cexpkxky = xx(n, i) * yy(m, j)

                     ntilda(i,j) = ntilda(i,j) 
     .                           + (delta_x * exk(n,m)
     .                           +  delta_y * eyk(n,m)
     .                           +  delta_z * ezk(n,m)) * cexpkxky

                  end do
               end do

            end if

            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')
     
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, ntilda,
     .   nxdim, -1, -1)

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end

c
c***************************************************************************
c

      subroutine current_cd(xjpx, xjpy, xjpz,
     .   nxdim, nydim, nnodex, nnodey,
     .   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     .   xm, q, xn, xkt, omgc, omgp2, lmax,
     .   xkxsav, xkysav, nzfun, ibessel,
     .   exk, eyk, ezk, nphi, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     .   xx, yy, isigma, xnuomg, psi, psilim, nboundary,
     .   myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj,
     .   zeffcd, clight, x, y, rt, b0, ftrap, omgrf,
     .   xjpx_ehst, xjpy_ehst, xjpz_ehst, nkperp,  zi, eps0, 
     .   v0i, xk0, kperp_max, i_sav, j_sav, upshift, damping, 
     .   xk_cutoff, odd_order, dx, dy, xkphi)

*-----------------------------------------------------------------------
*     This subroutine calculates the plasma current for a single species
*-----------------------------------------------------------------------

      implicit none

      logical ismine

      complex zi
      real eps0, v0i, xk0, kperp_max
      integer i_sav, j_sav, upshift 
      real damping, xk_cutoff

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt, nboundary
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey,
     .   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      integer ftrap, odd_order
      real dx, dy
      real zeffcd, clight, rt, b0, damp, r1, xm1, ceta, epsa, xmut2,
     .   eta0, xnexp, xjtild, signkz, cfit, c1, afit, yt, bmaxa, vphase,
     .   omgrf, c_ehst, akprl, xkprl, a, rmaxa, rmina,
     .   wphase, xlnlam, vth
      real x(nxdim), y(nydim), xkphi(nxdim)


      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim),
     .     psi(nxdim, nydim), psilim
      real capr(nxdim), xnuomg, delta0
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim),
     .                         bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

      complex xjpx_ehst(nxdim, nydim),
     .        xjpy_ehst(nxdim, nydim),
     .        xjpz_ehst(nxdim, nydim)


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
     .     uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim),
     .     uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
	real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)



      if(nphi .gt. 0) signkz =  1.0
      if(nphi .lt. 0) signkz = -1.0
      if(nphi .eq. 0) signkz =  1.0


      ceta = 8.27
      xnexp = 2.48
      cfit = .0987
      afit = 12.3
      damp = 23.83/(0.678 + zeffcd)


      xjpx(:,:) = 0.0
      xjpy(:,:) = 0.0
      xjpz(:,:) = 0.0

      xjpx_ehst(:,:) = 0.0
      xjpy_ehst(:,:) = 0.0
      xjpz_ehst(:,:) = 0.0


c--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (i - 1) * nnodey + j
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            xjpx_ehst(i,j) = 0.0
            xjpy_ehst(i,j) = 0.0
            xjpz_ehst(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     .                              .or.  nboundary .eq. 0)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

*                    ----------------
*                    Calculate C_ehst
*                    ----------------

c                     xkprl = uzx(i, j) * xkxsav(n)
c     .                    + uzy(i, j) * xkysav(m)
c     .                    + uzz(i, j) * xkphi(i)

                     xkprl = uzz(i, j) * xkphi(i)


                     akprl = abs(xkprl)
                     if(akprl .lt. 1.0e-05)akprl = 1.0e-05

                     vphase = omgrf / akprl
                     if (vphase .ge. clight)vphase = clight

                     c_ehst = 0.0

                     if(xkt(i, j) .ne. 0.0) then

                        xlnlam = 24. - alog(sqrt(xn(i, j)
     .                     /1.e+06) / (xkt(i, j)/ abs(q)))
                        vth = sqrt(xkt(i, j) / xm)
                        wphase = vphase / vth

                        if (wphase .lt. 10.0 .and. wphase .gt. 0.1)then

*                       -----------------------------
*                       Assume circular flux surfaces
*                       -----------------------------
                        a = sqrt(x(i)**2 + y(j)**2)
                        rmaxa = rt + a
                        rmina = rt - a
                        bmaxa = b0 * rt / rmina
*                       ------------------------------------
*                       End circular flux surface assumption
*                       ------------------------------------


                        epsa = (rmaxa - rmina)/(rmaxa + rmina)


                        xmut2 = 1. - bmod(i, j) / bmaxa
                        eta0 = 8. * wphase**2 / (5. + zeffcd) +
     .                     ceta / zeffcd**.707 + damp / wphase
                        r1 = 1. - epsa**.77 * sqrt(3.5**2 + wphase**2)
     .                     /(3.5 * epsa**.77 + wphase)
                        xm1 = 1.0
                        c1 = 1.0

                        if(ftrap .ne. 0)then
                           if (xmut2 .gt. 1.0E-06) then
                              xm1 = 1.0 + afit *
     .                                        (xmut2**1.5 / wphase**3)
                              yt = (1. - xmut2) * wphase**2 / xmut2
                              if(yt.ge.0.0)c1 = 1. -
     .                                          exp(-(cfit*yt)**xnexp)
                           endif
                        endif

                        xjtild = c1 * xm1 * eta0 * r1


                        c_ehst = -19.19E+15 / xlnlam *
     .                     (xkt(i, j)
     .                     / abs(q))/ xn(i, j) *
     .                     xjtild * signkz
c                        c_ehst = 1.0
                        end if
                     end if

                     if (vphase .ge. clight .or. akprl .lt. 1.0e-05)
     .                  c_ehst = 0.0



                     if (isigma .eq. 1)

     .                  call sigmad_cql3d(i, j, n, m, rho(i,j),rho_a,
     .                      gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     .                      xm, q, xn(i,j), xnuomg,
     .                      xkt(i,j), omgc(i,j), omgp2(i,j),
     .                      -lmax, lmax, nzfun, ibessel,
     .                      xkxsav(n), xkysav(m), nphi, capr(i),
     .                      bxn(i,j), byn(i,j), bzn(i,j),
     .                      uxx(i,j), uxy(i,j), uxz(i,j),
     .                      uyx(i,j), uyy(i,j), uyz(i,j),
     .                      uzx(i,j), uzy(i,j), uzz(i,j),
     .                      sigxx, sigxy, sigxz,
     .                      sigyx, sigyy, sigyz,
     .                      sigzx, sigzy, sigzz,
     .                      delta0, ndist,
     .                      nupar, nuper, n_psi,
     .                      n_psi_dim, dfduper, dfdupar,
     .                      UminPara, UmaxPara, UPERP, UPARA,
     .                      vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .                      nkperp, zi, eps0, v0i, omgrf, xk0,
     .                      kperp_max, i_sav, j_sav, upshift, 
     .                      damping, xk_cutoff, rt, nkx2, nky2)
     
                        if(odd_order .ne. 0) 
     .                     call odd_order_derivs(i, j, n, m, nxdim, 
     .                     nydim, 
     .                     nnodex, nnodey, mkdim1, mkdim2, capr,
     .                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, 
     .                     xm, q, xn, xkt, omgc, omgp2, lmax,
     .                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi, 
     .                     bxn, byn, bzn,
     .                     uxx, uxy, uxz,
     .                     uyx, uyy, uyz,
     .                     uzx, uzy, uzz,
     .                     sigxx, sigxy, sigxz,
     .                     sigyx, sigyy, sigyz,
     .                     sigzx, sigzy, sigzz,        
     .                     myrow, mycol, nprow, npcol, icontxt, 
     .                     desc_amat, dlen_, nboundary,
     .                     xx, yy, isigma, xnuomg, psi, psilim, 
     .                     myid, nproc, delta0, gradprlb, bmod,
     .                     bmod_mid, ndist, nupar, nuper, n_psi,
     .                     n_psi_dim, dfduper, dfdupar,
     .                     UminPara, UmaxPara, UPERP, UPARA,
     .                     vc_mks, df_cql_uprp, df_cql_uprl, rho, 
     .                     rho_a, nbessj, nkperp,
     .                     zi, eps0, v0i, omgrf, xk0, kperp_max, 
     .                     i_sav, j_sav, upshift, 
     .                     damping, xk_cutoff, rt, dx, dy)     
          

                     if (isigma .eq. 0)

     .                  call sigmac_stix(i, j, n, m,
     .                      xm, q, xn(i,j), xnuomg,
     .                      xkt(i,j), omgc(i,j), omgp2(i,j),
     .                      -lmax, lmax, nzfun, ibessel,
     .                      xkxsav(n), xkysav(m), nphi, capr(i),
     .                      bxn(i,j), byn(i,j), bzn(i,j),
     .                      uxx(i,j), uxy(i,j), uxz(i,j),
     .                      uyx(i,j), uyy(i,j), uyz(i,j),
     .                      uzx(i,j), uzy(i,j), uzz(i,j),
     .                      sigxx, sigxy, sigxz,
     .                      sigyx, sigyy, sigyz,
     .                      sigzx, sigzy, sigzz,
     .                      delta0, zi, eps0, v0i, omgrf, xk0,
     .                      kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m)
     .                               + sigxy * eyk(n,m)
     .                               + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m)
     .                               + sigyy * eyk(n,m)
     .                               + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m)
     .                               + sigzy * eyk(n,m)
     .                               + sigzz * ezk(n,m)) * cexpkxky


                     xjpx_ehst(i,j) = xjpx_ehst(i,j) + (sigxx * exk(n,m)
     .                           + sigxy * eyk(n,m)
     .                           + sigxz * ezk(n,m)) * cexpkxky * c_ehst
                     xjpy_ehst(i,j) = xjpy_ehst(i,j) + (sigyx * exk(n,m)
     .                           + sigyy * eyk(n,m)
     .                           + sigyz * ezk(n,m)) * cexpkxky * c_ehst
                     xjpz_ehst(i,j) = xjpz_ehst(i,j) + (sigzx * exk(n,m)
     .                           + sigzy * eyk(n,m)
     .                           + sigzz * ezk(n,m)) * cexpkxky * c_ehst

                  end do
               end do

            end if

            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz,
     .   nxdim, -1, -1)


      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx_ehst,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy_ehst,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz_ehst,
     .   nxdim, -1, -1)


      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end



c
c***************************************************************************
c

      subroutine current(xjpx, xjpy, xjpz, nxdim, nydim, nnodex, nnodey,
     .   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     .   xm, q, xn, xkt, omgc, omgp2, lmax,
     .   xkxsav, xkysav, nzfun, ibessel,
     .   exk, eyk, ezk, nphi, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     .   xx, yy, isigma, xnuomg, psi, psilim, nboundary,
     .   myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj, nkperp,
     .   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, 
     .   damping, xk_cutoff, rt, odd_order, dx, dy, xkphi)

*-----------------------------------------------------------------------
*     This subroutine calculates the plasma current for a single species
*-----------------------------------------------------------------------

      implicit none

      logical ismine

      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max
      integer i_sav, j_sav, upshift 
      real damping, xk_cutoff, rt

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt, nboundary, odd_order
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim),
     .     psi(nxdim, nydim), psilim, xkphi(nxdim)
      real capr(nxdim), xnuomg, delta0, dx, dy
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim),
     .                         bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

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
     .     uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim),
     .     uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
	real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)


      xjpx(:,:) = 0.0
      xjpy(:,:) = 0.0
      xjpz(:,:) = 0.0


c--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (i - 1) * nnodey + j
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     .                              .or.  nboundary .eq. 0)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     if (isigma .eq. 1)

     .                  call sigmad_cql3d(i, j, n, m, rho(i,j),rho_a,
     .                      gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     .                      xm, q, xn(i,j), xnuomg,
     .                      xkt(i,j), omgc(i,j), omgp2(i,j),
     .                      -lmax, lmax, nzfun, ibessel,
     .                      xkxsav(n), xkysav(m), nphi, capr(i),
     .                      bxn(i,j), byn(i,j), bzn(i,j),
     .                      uxx(i,j), uxy(i,j), uxz(i,j),
     .                      uyx(i,j), uyy(i,j), uyz(i,j),
     .                      uzx(i,j), uzy(i,j), uzz(i,j),
     .                      sigxx, sigxy, sigxz,
     .                      sigyx, sigyy, sigyz,
     .                      sigzx, sigzy, sigzz,
     .                      delta0, ndist,
     .                      nupar, nuper, n_psi,
     .                      n_psi_dim, dfduper, dfdupar,
     .                      UminPara, UmaxPara, UPERP, UPARA,
     .                      vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .                      nkperp, zi, eps0, v0i, omgrf, xk0,
     .                      kperp_max, i_sav, j_sav, upshift, 
     .                      damping, xk_cutoff, rt, nkx2, nky2)
     
                        if(odd_order .ne. 0) 
     .                     call odd_order_derivs(i, j, n, m, nxdim, 
     .                     nydim, 
     .                     nnodex, nnodey, mkdim1, mkdim2, capr,
     .                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, 
     .                     xm, q, xn, xkt, omgc, omgp2, lmax,
     .                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi, 
     .                     bxn, byn, bzn,
     .                     uxx, uxy, uxz,
     .                     uyx, uyy, uyz,
     .                     uzx, uzy, uzz,
     .                     sigxx, sigxy, sigxz,
     .                     sigyx, sigyy, sigyz,
     .                     sigzx, sigzy, sigzz,        
     .                     myrow, mycol, nprow, npcol, icontxt, 
     .                     desc_amat, dlen_, nboundary,
     .                     xx, yy, isigma, xnuomg, psi, psilim, 
     .                     myid, nproc, delta0, gradprlb, bmod,
     .                     bmod_mid, ndist, nupar, nuper, n_psi,
     .                     n_psi_dim, dfduper, dfdupar,
     .                     UminPara, UmaxPara, UPERP, UPARA,
     .                     vc_mks, df_cql_uprp, df_cql_uprl, rho, 
     .                     rho_a, nbessj, nkperp,
     .                     zi, eps0, v0i, omgrf, xk0, kperp_max, 
     .                     i_sav, j_sav, upshift, 
     .                     damping, xk_cutoff, rt, dx, dy)      

                     if (isigma .eq. 0)

     .                  call sigmac_stix(i, j, n, m,
     .                      xm, q, xn(i,j), xnuomg,
     .                      xkt(i,j), omgc(i,j), omgp2(i,j),
     .                      -lmax, lmax, nzfun, ibessel,
     .                      xkxsav(n), xkysav(m), nphi, capr(i),
     .                      bxn(i,j), byn(i,j), bzn(i,j),
     .                      uxx(i,j), uxy(i,j), uxz(i,j),
     .                      uyx(i,j), uyy(i,j), uyz(i,j),
     .                      uzx(i,j), uzy(i,j), uzz(i,j),
     .                      sigxx, sigxy, sigxz,
     .                      sigyx, sigyy, sigyz,
     .                      sigzx, sigzy, sigzz,
     .                      delta0, zi, eps0, v0i, omgrf, xk0,
     .                      kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m)
     .                               + sigxy * eyk(n,m)
     .                               + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m)
     .                               + sigyy * eyk(n,m)
     .                               + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m)
     .                               + sigzy * eyk(n,m)
     .                               + sigzz * ezk(n,m)) * cexpkxky

                  end do
               end do

            end if

            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz,
     .   nxdim, -1, -1)

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end


c
c***************************************************************************
c


      subroutine current_1(xjpx, xjpy, xjpz, nxdim, nydim, 
     .   nnodex, nnodey,
     .   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     .   xm, q, xn, xkt, omgc, omgp2, lmax,
     .   xkxsav, xkysav, nzfun, ibessel,
     .   exk, eyk, ezk, nphi, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     .   xx, yy, isigma, xnuomg, psi, psilim, nboundary,
     .   myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj, nkperp,
     .   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, 
     .   damping, xk_cutoff, rt, odd_order, dx, dy, xkphi)

*-----------------------------------------------------------------------
*     This subroutine calculates the plasma current for a single species
*-----------------------------------------------------------------------

      implicit none

      logical ismine

      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max
      integer i_sav, j_sav, upshift 
      real damping, xk_cutoff, rt

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt, nboundary, odd_order
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim),
     .     psi(nxdim, nydim), psilim, dx, dy
     
      real capr(nxdim), xnuomg, delta0, xkphi(nxdim)
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim),
     .                         bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

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
     .     uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim),
     .     uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
	real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)


      xjpx(:,:) = 0.0
      xjpy(:,:) = 0.0
      xjpz(:,:) = 0.0


c--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (i - 1) * nnodey + j
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     .                              .or.  nboundary .eq. 0)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     if (isigma .eq. 1)

     .                  call sigmad_cql3d_1(i, j, n, m, rho(i,j),rho_a,
     .                      gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     .                      xm, q, xn(i,j), xnuomg,
     .                      xkt(i,j), omgc(i,j), omgp2(i,j),
     .                      -lmax, lmax, nzfun, ibessel,
     .                      xkxsav(n), xkysav(m), nphi, capr(i),
     .                      bxn(i,j), byn(i,j), bzn(i,j),
     .                      uxx(i,j), uxy(i,j), uxz(i,j),
     .                      uyx(i,j), uyy(i,j), uyz(i,j),
     .                      uzx(i,j), uzy(i,j), uzz(i,j),
     .                      sigxx, sigxy, sigxz,
     .                      sigyx, sigyy, sigyz,
     .                      sigzx, sigzy, sigzz,
     .                      delta0, ndist,
     .                      nupar, nuper, n_psi,
     .                      n_psi_dim, dfduper, dfdupar,
     .                      UminPara, UmaxPara, UPERP, UPARA,
     .                      vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .                      nkperp, zi, eps0, v0i, omgrf, xk0,
     .                      kperp_max, i_sav, j_sav, upshift, 
     .                      damping, xk_cutoff, rt, nkx2, nky2)
     
                        if(odd_order .ne. 0) 
     .                     call odd_order_derivs(i, j, n, m, nxdim, 
     .                     nydim, 
     .                     nnodex, nnodey, mkdim1, mkdim2, capr,
     .                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, 
     .                     xm, q, xn, xkt, omgc, omgp2, lmax,
     .                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi, 
     .                     bxn, byn, bzn,
     .                     uxx, uxy, uxz,
     .                     uyx, uyy, uyz,
     .                     uzx, uzy, uzz,
     .                     sigxx, sigxy, sigxz,
     .                     sigyx, sigyy, sigyz,
     .                     sigzx, sigzy, sigzz,        
     .                     myrow, mycol, nprow, npcol, icontxt, 
     .                     desc_amat, dlen_, nboundary,
     .                     xx, yy, isigma, xnuomg, psi, psilim, 
     .                     myid, nproc, delta0, gradprlb, bmod,
     .                     bmod_mid, ndist, nupar, nuper, n_psi,
     .                     n_psi_dim, dfduper, dfdupar,
     .                     UminPara, UmaxPara, UPERP, UPARA,
     .                     vc_mks, df_cql_uprp, df_cql_uprl, rho, 
     .                     rho_a, nbessj, nkperp,
     .                     zi, eps0, v0i, omgrf, xk0, kperp_max, 
     .                     i_sav, j_sav, upshift, 
     .                     damping, xk_cutoff, rt, dx, dy)      

                     if (isigma .eq. 0)

     .                  call sigmac_stix(i, j, n, m,
     .                      xm, q, xn(i,j), xnuomg,
     .                      xkt(i,j), omgc(i,j), omgp2(i,j),
     .                      -lmax, lmax, nzfun, ibessel,
     .                      xkxsav(n), xkysav(m), nphi, capr(i),
     .                      bxn(i,j), byn(i,j), bzn(i,j),
     .                      uxx(i,j), uxy(i,j), uxz(i,j),
     .                      uyx(i,j), uyy(i,j), uyz(i,j),
     .                      uzx(i,j), uzy(i,j), uzz(i,j),
     .                      sigxx, sigxy, sigxz,
     .                      sigyx, sigyy, sigyz,
     .                      sigzx, sigzy, sigzz,
     .                      delta0, zi, eps0, v0i, omgrf, xk0,
     .                      kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m)
     .                               + sigxy * eyk(n,m)
     .                               + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m)
     .                               + sigyy * eyk(n,m)
     .                               + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m)
     .                               + sigzy * eyk(n,m)
     .                               + sigzz * ezk(n,m)) * cexpkxky

                  end do
               end do

            end if

            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz,
     .   nxdim, -1, -1)

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end


c
c***************************************************************************
c


      subroutine current_2(xjpx, xjpy, xjpz, nxdim, nydim, 
     .   nnodex, nnodey,
     .   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     .   xm, q, xn, xkt, omgc, omgp2, lmax,
     .   xkxsav, xkysav, nzfun, ibessel,
     .   exk, eyk, ezk, nphi, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     .   xx, yy, isigma, xnuomg, psi, psilim, nboundary,
     .   myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper, dfdupar,
     .   UminPara, UmaxPara, UPERP, UPARA,
     .   vc_mks, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj, nkperp,
     .   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, 
     .   damping, xk_cutoff, rt, odd_order, dx, dy, xkphi)

*-----------------------------------------------------------------------
*     This subroutine calculates the plasma current for a single species
*-----------------------------------------------------------------------

      implicit none

      logical ismine

      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max
      integer i_sav, j_sav, upshift 
      real damping, xk_cutoff, rt

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt, nboundary, odd_order
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim),
     .     psi(nxdim, nydim), psilim, xkphi(nxdim)
      real capr(nxdim), xnuomg, delta0, dx, dy
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim),
     .                         bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

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
     .     uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim),
     .     uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
	real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)


      xjpx(:,:) = 0.0
      xjpy(:,:) = 0.0
      xjpz(:,:) = 0.0


c--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (i - 1) * nnodey + j
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     .                              .or.  nboundary .eq. 0)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     if (isigma .eq. 1)

     .                  call sigmad_cql3d_2(i, j, n, m, rho(i,j),rho_a,
     .                      gradprlb(i,j), bmod(i,j), bmod_mid(i,j),
     .                      xm, q, xn(i,j), xnuomg,
     .                      xkt(i,j), omgc(i,j), omgp2(i,j),
     .                      -lmax, lmax, nzfun, ibessel,
     .                      xkxsav(n), xkysav(m), nphi, capr(i),
     .                      bxn(i,j), byn(i,j), bzn(i,j),
     .                      uxx(i,j), uxy(i,j), uxz(i,j),
     .                      uyx(i,j), uyy(i,j), uyz(i,j),
     .                      uzx(i,j), uzy(i,j), uzz(i,j),
     .                      sigxx, sigxy, sigxz,
     .                      sigyx, sigyy, sigyz,
     .                      sigzx, sigzy, sigzz,
     .                      delta0, ndist,
     .                      nupar, nuper, n_psi,
     .                      n_psi_dim, dfduper, dfdupar,
     .                      UminPara, UmaxPara, UPERP, UPARA,
     .                      vc_mks, df_cql_uprp, df_cql_uprl, nbessj,
     .                      nkperp, zi, eps0, v0i, omgrf, xk0,
     .                      kperp_max, i_sav, j_sav, upshift, 
     .                      damping, xk_cutoff, rt, nkx2, nky2)
     
                        if(odd_order .ne. 0) 
     .                     call odd_order_derivs(i, j, n, m, nxdim, 
     .                     nydim, 
     .                     nnodex, nnodey, mkdim1, mkdim2, capr,
     .                     nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, 
     .                     xm, q, xn, xkt, omgc, omgp2, lmax,
     .                     xkxsav, xkysav, xkphi, nzfun, ibessel, nphi, 
     .                     bxn, byn, bzn,
     .                     uxx, uxy, uxz,
     .                     uyx, uyy, uyz,
     .                     uzx, uzy, uzz,
     .                     sigxx, sigxy, sigxz,
     .                     sigyx, sigyy, sigyz,
     .                     sigzx, sigzy, sigzz,        
     .                     myrow, mycol, nprow, npcol, icontxt, 
     .                     desc_amat, dlen_, nboundary,
     .                     xx, yy, isigma, xnuomg, psi, psilim, 
     .                     myid, nproc, delta0, gradprlb, bmod,
     .                     bmod_mid, ndist, nupar, nuper, n_psi,
     .                     n_psi_dim, dfduper, dfdupar,
     .                     UminPara, UmaxPara, UPERP, UPARA,
     .                     vc_mks, df_cql_uprp, df_cql_uprl, rho, 
     .                     rho_a, nbessj, nkperp,
     .                     zi, eps0, v0i, omgrf, xk0, kperp_max, 
     .                     i_sav, j_sav, upshift, 
     .                     damping, xk_cutoff, rt, dx, dy)      

                     if (isigma .eq. 0)

     .                  call sigmac_stix(i, j, n, m,
     .                      xm, q, xn(i,j), xnuomg,
     .                      xkt(i,j), omgc(i,j), omgp2(i,j),
     .                      -lmax, lmax, nzfun, ibessel,
     .                      xkxsav(n), xkysav(m), nphi, capr(i),
     .                      bxn(i,j), byn(i,j), bzn(i,j),
     .                      uxx(i,j), uxy(i,j), uxz(i,j),
     .                      uyx(i,j), uyy(i,j), uyz(i,j),
     .                      uzx(i,j), uzy(i,j), uzz(i,j),
     .                      sigxx, sigxy, sigxz,
     .                      sigyx, sigyy, sigyz,
     .                      sigzx, sigzy, sigzz,
     .                      delta0, zi, eps0, v0i, omgrf, xk0,
     .                      kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m)
     .                               + sigxy * eyk(n,m)
     .                               + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m)
     .                               + sigyy * eyk(n,m)
     .                               + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m)
     .                               + sigzy * eyk(n,m)
     .                               + sigzz * ezk(n,m)) * cexpkxky

                  end do
               end do

            end if

            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz,
     .   nxdim, -1, -1)

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end


c
c***************************************************************************
c


      subroutine cur_slo(xjpx, xjpy, xjpz, nxdim, nydim,
     .   nnodex, nnodey,
     .   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     .   xm, q, xn, eslow, omgc, omgp2, lmax,
     .   xkxsav, xkysav, nzfun, ibessel,
     .   exk, eyk, ezk, nphi, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_,
     .   xx, yy, isigma, xnuomg, psi, psilim, nboundary,
     .   xkte, zeffcd, myid, nproc,
     .   zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav)

      implicit none

      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max
      integer i_sav, j_sav

      logical ismine


      integer nproc, myid, ngrid, id
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx,
     .   icontxt, nboundary
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey,
     1   nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim),
     .     psi(nxdim, nydim), psilim, zeffcd
      real capr(nxdim), xnuomg

      complex xx(nkdim1 : nkdim2, 1 : nxdim),
     .        yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

      complex sigxx, sigxy, sigxz,
     1        sigyx, sigyy, sigyz,
     1        sigzx, sigzy, sigzz


      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        eyk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real q, xn(nxdim, nydim), eslow, xkte(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim),
     .     uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim),
     .     uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)


c--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (j - 1) * nnodex + i
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1
     .                              .or.  nboundary .eq. 0)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     call sigmah_slow(i, j, n, m,
     .                   xm, q, xn(i,j), xnuomg,
     .                   eslow, omgc(i,j), omgp2(i,j),
     .                   -lmax, lmax, nzfun, ibessel,
     .                   xkxsav(n), xkysav(m), nphi, capr(i),
     .                   bxn(i,j), byn(i,j), bzn(i,j),
     .                   uxx(i,j), uxy(i,j), uxz(i,j),
     .                   uyx(i,j), uyy(i,j), uyz(i,j),
     .                   uzx(i,j), uzy(i,j), uzz(i,j),
     .                   sigxx, sigxy, sigxz,
     .                   sigyx, sigyy, sigyz,
     .                   sigzx, sigzy, sigzz,
     .                   xkte(i,j), zeffcd, zi, eps0, v0i, omgrf, xk0,
     .                   kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m)
     1                               + sigxy * eyk(n,m)
     1                               + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m)
     1                               + sigyy * eyk(n,m)
     1                               + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m)
     1                               + sigzy * eyk(n,m)
     1                               + sigzz * ezk(n,m)) * cexpkxky

                  end do
               end do

            end if
            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy,
     .   nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz,
     .   nxdim, -1, -1)

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end

c
c***************************************************************************
c

