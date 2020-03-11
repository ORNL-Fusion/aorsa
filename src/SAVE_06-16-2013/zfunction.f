


      subroutine z_table(sgn_kprl, zeta, gamma, xnu, z0, z1, z2)

*     ----------------------------------------------------
*     z_table is a wrapper routine for David Smithe's table
*     lookup Z function
*     ----------------------------------------------------

      implicit none

      complex, intent(inout)   :: z0, z1, z2
      complex, intent(in)    :: zeta
      real, intent(in)       :: sgn_kprl, gamma, xnu
      complex, dimension(0:4)  :: zfunc

      real alpha_s
      complex zi

      zi = cmplx(0.0, 1.0)

      alpha_s = -2.0 * gamma


      if (sgn_kprl .ge. 0.0) then

         call ztable_zfunc( real(zeta), alpha_s, xnu, zfunc )
         z0 = zfunc(0)
         z1 = - 0.5 * zfunc(1) - gamma / (2.0 * zi) * zfunc(2)
         z2 = 0.5 * zfunc(0) + 0.25 * zfunc(2)
     .      + gamma / (2.0 * zi) * zfunc(3)
     .      - gamma**2 / 4. * zfunc(4)


      else if (sgn_kprl .lt. 0.0) then

         call ztable_zfunc( - real(zeta), - alpha_s, - xnu, zfunc )
         z0 = - zfunc(0)
         z1 = - 0.5 * zfunc(1) + gamma / (2.0 * zi) * zfunc(2)
         z2 = - 0.5 * zfunc(0) - 0.25 * zfunc(2)
     .      + gamma / (2.0 * zi) * zfunc(3)
     .      + gamma**2 / 4. * zfunc(4)


      end if

      return

      end subroutine z_table



c
c*******************************************************************************
c

      subroutine z_exact(sgn_kprl, zeta, gamma, xnu, z0, z1, z2,
     .   dz0, dz1, dz2)

*     ----------------------------------------------------
*     z_exact is a wrapper routine for David Smithe's own
*     version of "z_smithe"
*     ----------------------------------------------------

      implicit none

      complex, intent(inout)   :: z0, z1, z2, dz0, dz1, dz2
      complex, intent(in)    :: zeta
      real, intent(in)       :: sgn_kprl, gamma, xnu
      complex, dimension(0:4)  :: zfunc

      real alpha_s
      complex zi

      zi = cmplx(0.0, 1.0)

      alpha_s = -2.0 * gamma

      if ( sgn_kprl .ge. 0.0 ) then

         call ztable_zfunc( real(zeta), alpha_s, xnu, zfunc )
         z0 = zfunc(0)
         z1 = - 0.5 * zfunc(1) - gamma / (2.0 * zi) * zfunc(2)
         z2 = 0.5 * zfunc(0) + 0.25 * zfunc(2)
     .      + gamma / (2.0 * zi) * zfunc(3)
     .      - gamma**2 / 4. * zfunc(4)

         dz0 = zfunc(1)
         dz1 = - 0.5 * zfunc(2) - gamma / (2.0 * zi) * zfunc(3)
c         dz2 = 0.5 * zfunc(1) + 0.25 * zfunc(3)
c     .      + gamma / (2.0 * zi) * zfunc(4)



      else if ( sgn_kprl .lt. 0.0 ) then

         call ztable_zfunc( - real(zeta), - alpha_s, - xnu, zfunc )
         z0 = - zfunc(0)
         z1 = - 0.5 * zfunc(1) + gamma / (2.0 * zi) * zfunc(2)
         z2 = - 0.5 * zfunc(0) - 0.25 * zfunc(2)
     .      + gamma / (2.0 * zi) * zfunc(3)
     .      + gamma**2 / 4. * zfunc(4)

         dz0 =  zfunc(1)
         dz1 =  0.5 * zfunc(2) - gamma / (2.0 * zi) * zfunc(3)
c         dz2 =  0.5 * zfunc(1) + 0.25 * zfunc(3)
c     .      - gamma / (2.0 * zi) * zfunc(4)

      end if

      return

      end subroutine z_exact


c
c***************************************************************************
c


      subroutine z_approx(sgn_kprl, zeta, gamma, z0, z1, z2)
      use z_erfc
c      use zfun_hilbert
      implicit none

*     ------------------
*     Finite k_parallel:
*     ------------------
      real gamma, fgam, y0, y, sgn_kprl, descrim
      complex zeta, z0, z1, z2, zfunct, zetat, fzeta
      logical :: flag
      integer :: ieer ! fft error flag 
           

      
      y0 = 1.5
      y = y0            


      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gamma .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gamma * y) - 1.)
     .         / (2. * gamma * y)
         endif

         zetat = fgam * zeta
	 
c        zfunct = fzeta(zetat)
         call zfun (zetat, zfunct)
c        call zfun_erfc (zetat, zfunct)
c	 call zfun_hil(zetat, zfunct, flag)
	 
	 
         z0 = fgam * zfunct
         z1 = fgam * (1.0 + fgam * zeta * zfunct)
         z2 =  fgam**2 * zeta * (1.0 + fgam * zeta * zfunct)
      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gamma .gt. 1.0e-05)then
            descrim = 1. - 4. * gamma * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gamma * y) )
     .         / (2. * gamma * y)
         endif


         zetat = - fgam * zeta
	 
c         zfunct = fzeta( zetat)
         call zfun (zetat, zfunct)
c        call zfun_erfc (zetat, zfunct)
c	 call zfun_hil(zetat, zfunct, flag)	 
	 	 
         z0 = - fgam * zfunct
         z1 = y / abs(y) * fgam * (1.0 - fgam * zeta * zfunct)
         z2 =  fgam**2 * zeta * (1.0 - fgam * zeta * zfunct)
      end if


      return
      end

c
c******************************************************************************
c


      subroutine z_approx_im(sgn_kprl, zeta, gamma, z0, z1, z2)

      implicit none

*     ------------------
*     Finite k_parallel:
*     ------------------
      real gamma, fgam, y0, y, sgn_kprl, descrim
      complex zeta, z0, z1, z2, zfunct, zetat, fzeta


      y0 = 1.5
      y = y0


      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gamma .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gamma * y) - 1.)
     .         / (2. * gamma * y)
         endif

         zetat = fgam * zeta
c        zfunct = fzeta(zetat)
         call zfun_im (zetat, zfunct)

         z0 = fgam * zfunct
         z1 = fgam * (1.0 + fgam * zeta * zfunct)
         z2 =  fgam**2 * zeta * (1.0 + fgam * zeta * zfunct)
      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gamma .gt. 1.0e-05)then
            descrim = 1. - 4. * gamma * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gamma * y) )
     .         / (2. * gamma * y)
         endif


         zetat = - fgam * zeta
c         zfunct = fzeta( zetat)
         call zfun_im (zetat, zfunct)

         z0 = - fgam * zfunct
         z1 = y / abs(y) * fgam * (1.0 - fgam * zeta * zfunct)
         z2 =  fgam**2 * zeta * (1.0 - fgam * zeta * zfunct)
      end if


      return
      end

c
c******************************************************************************
c

      subroutine z_approx0(sgn_kprl, alpha, xkprl_eff, zeta_eff,
     .   al, bl, cl)
      use z_erfc
c      use zfun_hilbert      
      implicit none

*     -----------------------
*     Small k_parallel limit:
*     -----------------------
      real sgn_kprl, alpha, xkprl_eff
      complex zeta_eff, al, bl, cl, zfunct, zetat, fzeta
      logical :: flag
      integer :: ieer ! fft error flag 
            

      if(sgn_kprl .ge. 0.0)then
      
c        zfunct = fzeta(zeta_eff)
         call zfun (zeta_eff, zfunct)
c        call zfun_erfc (zeta_eff, zfunct)
c	 call zfun_hil(zeta_eff, zfunct, flag)	 
	 
         al = 1. /(xkprl_eff * alpha) * zfunct
         bl = 1. /(xkprl_eff * alpha) * (1. + zeta_eff * zfunct)
         cl = zeta_eff /(xkprl_eff * alpha) *(1. + zeta_eff * zfunct)
      end if


      if(sgn_kprl .lt. 0.0)then
c        zfunct = fzeta(-zeta_eff)
         call zfun (-zeta_eff, zfunct)
c        call zfun_erfc (-zeta_eff, zfunct)
c	 call zfun_hil(-zeta_eff, zfunct, flag)	 	 

         al = - 1. /(xkprl_eff * alpha) * zfunct
         bl = - 1. /(xkprl_eff * alpha) * (1. - zeta_eff * zfunct)
         cl = zeta_eff /(xkprl_eff * alpha) * (1. - zeta_eff * zfunct)
            end if


      return
      end


c
c***************************************************************************
c
      subroutine z2pole(zeta, zfunct)

      complex zeta, zfunct, zi, znum
      data pi/3.141592654/
      data zi/(0.0,1.0)/
      data g/1.141592654/
      data sqpi/1.772453851/
      znum = zi*sqpi + g*zeta
      zfunct = znum/(1.0 - zeta*znum)
      return
      end
c
c******************************************************************************
c

      subroutine z_smithe(sgn_kprl, zeta, gamma, z0, z1, z2)

      implicit none

      real dx, xmax, xnueff, gamma, sgn_kprl, xmax0
      integer j, npts, nmaxdim
      real x(10000), zeta_mod
      complex z0, z1, z2, zi, f(10000), zeta
      complex cexpx(10000)

      parameter (nmaxdim = 10000)

      zi = cmplx(0.,1.)
      xnueff = 0.0
      zeta_mod = sqrt(conjg(zeta) * zeta)


      if(zeta_mod .gt. 1.0e+02)then

         z0 = - 1.0 / zeta - 0.5 / zeta**3
     .                          + 3.0 * gamma / (zi * zeta**4)
         z1 = - 0.5 / zeta**2 + gamma / (zi * zeta**3)
     .                                       - 0.75 / zeta**4
         z2 = - 0.5 / zeta - 0.75 / zeta**3
     .                          + 4.5 * gamma / (zi * zeta**4)
         return

      else

         xmax0 = 7.0
         npts = 10000
         xmax = sgn_kprl * xmax0
         dx = xmax / (npts - 1)


*        ------------
*        calculate Z0
*        ------------
         do  j = 1, npts

            x(j) = (j-1) * dx
            cexpx(j) = exp(zi * zeta * x(j)
     .         - x(j)**2 / 4. * (1. + gamma * x(j))**2
     .         - xnueff / 8. * x(j)**3   )

            f(j) = zi * cexpx(j)
         end do

         call dgrat1(x, f, 1, npts, z0, nmaxdim)


*        ------------
*        calculate Z1
*        ------------
         do  j = 1, npts
            f(j) = 0.5 * x(j)* (1.0 + gamma * x(j)) * cexpx(j)
         end do

         call dgrat1(x, f, 1, npts, z1, nmaxdim)


*        ------------
*        calculate Z2
*        ------------
         do  j = 1, npts
            f(j) = zi / 2.0 * (1.0 -  x(j)**2 / 2.
     .           * (1. + gamma * x(j))**2) * cexpx(j)
         end do

         call dgrat1(x, f, 1, npts, z2, nmaxdim)



         return

      end if

      end

c
c*******************************************************************************
c

      subroutine dgrat1(x, f, nx1, nx2, ans, nmaxdim)

      implicit none

      integer n, nx1, nx2, nmaxdim
      real x(nmaxdim), dx
      complex f(nmaxdim), ans

      ans = 0.0

      do n = nx1, nx2 - 1
         dx = x(n+1) - x(n)
         ans = ans + dx * (f(n) + f(n+1)) / 2.0
      end do

      return
      end


c
c*******************************************************************************
c
