!
!********************************************************************
!

      subroutine znew_integral1(zetai, dKdL, z0, z1, z2)

      implicit none

      real dz, zmax
      integer j, npts, nmaxdim, j_filt
      
      parameter (nmaxdim = 10000) 
         
      real z(nmaxdim), dKdL, zetai, mult
      complex z0, z1, z2, zi, f2(nmaxdim), beta, cexpx
      complex f1(nmaxdim), f0(nmaxdim)

      zi = cmplx(0.,1.)

      npts = nmaxdim
c      zmax = 10000.
c      dz = zmax / (npts - 1)
      dz = 0.0628            
      if(dz .gt. 1.0 / (5. * abs(zetai))) dz = 1.0 / (5. * abs(zetai))            
      mult = 1.0     

*     --------------------
*     calculate Z0, Z1, Z2
*     --------------------
      j_filt = npts / 2
      
      f0 = 0.0
      f1 = 0.0
      f2 = 0.0           
      
      do  j = 1, npts
         z(j) = (j - 1) * dz     
         beta = 1. / sqrt(1. - 0.5 * zi * dKdL * z(j)**2 )
         if (j .ge. npts - j_filt) then
             mult = cos((float(j - j_filt) / float(j_filt))
     .            * 2.0 * atan(1.0))   
         end if             
         cexpx = exp(zi * z(j) - 0.25 * (z(j) * beta * zetai)**2) * mult         
         f0(j) = zi * zetai * beta * cexpx
         
         f1(j) = 0.5 * zetai**2 * beta**3 * z(j) * cexpx
         
         f2(j) = 0.5 * zi * zetai * beta**3
     .      * (1.0 -  0.5 * (beta * zetai * z(j))**2) * cexpx              
     
      end do      

      z0 = sum(f0) * dz
      z1 = sum(f1) * dz
      z2 = sum(f2) * dz             

      return
      end



c
c*******************************************************************************
c

      subroutine z_approx_e(sgn_kprl, zeta, gamma, 
     .   z0, z1, z2, zetai_table, 
     .   dKdL_table, z0_table, z1_table, z2_table, 
     .   dKdL_giv, nmax, mmax, ntable, mtable)
     
      use z_erfc
c      use zfun_hilbert
      implicit none

*     ------------------
*     Finite k_parallel:
*     ------------------
      real gamma, fgam, y0, y, sgn_kprl, descrim
      complex zeta, z0, z1, z2, zfunct, zetat, fzeta
      complex z0_out, z1_out, z2_out
      logical :: flag
      integer :: ieer ! fft error flag 
      integer nmax, mmax, ntable, mtable
      
      complex z0_table(ntable, mtable)
      complex z1_table(ntable, mtable)      
      complex z2_table(ntable, mtable)      
      
      real zetai_table(ntable), zetai_giv
      real dKdL_table(mtable), dKdL_giv
           
      
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
c        call zfun_hil(zetat, zfunct, flag)
         
         
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
c        call zfun_hil(zetat, zfunct, flag)      
                 
         z0 = - fgam * zfunct
         z1 = y / abs(y) * fgam * (1.0 - fgam * zeta * zfunct)
         z2 =  fgam**2 * zeta * (1.0 - fgam * zeta * zfunct)                     
         
      end if
      
      zetai_giv = 1.0 / real(zeta) 
                 
      call intplt_z(zetai_giv, dKdL_giv, z0_out, nmax, mmax, z0_table, 
     .   ntable, mtable, zetai_table, dKdL_table)
     
      call intplt_z(zetai_giv, dKdL_giv, z1_out, nmax, mmax, z1_table, 
     .   ntable, mtable, zetai_table, dKdL_table)
     
      call intplt_z(zetai_giv, dKdL_giv, z2_out, nmax, mmax, z2_table, 
     .   ntable, mtable, zetai_table, dKdL_table)                     
     
      if (z0_out .ne. 0.0) z0 = z0_out
      if (z1_out .ne. 0.0) z1 = z1_out
      if (z2_out .ne. 0.0) z2 = z2_out
                  
      return
      end

!
!********************************************************************
!

      subroutine z_approx_i(sgn_kprl, zeta, gamma, 
     .   z0, z1, z2, zetai_table, 
     .   dKdL_table, z0_table, z1_table, z2_table, 
     .   dKdL_giv, nmax, mmax, ntable, mtable)
     
      use z_erfc
c      use zfun_hilbert
      implicit none

*     ------------------
*     Finite k_parallel:
*     ------------------
      real gamma, fgam, y0, y, sgn_kprl, descrim, gammad
      complex zeta, z0, z1, z2, zfunct, zetat, fzeta
      complex z0_out, z1_out, z2_out
      logical :: flag
      integer :: ieer ! fft error flag 
      integer nmax, mmax, ntable, mtable
      
      complex z0_table(ntable, mtable)
      complex z1_table(ntable, mtable)      
      complex z2_table(ntable, mtable)      
      
      real zetai_table(ntable), zetai_giv
      real dKdL_table(mtable), dKdL_giv
              
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
c        call zfun_hil(zetat, zfunct, flag)
 
 
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
c        call zfun_hil(zetat, zfunct, flag)      
                 
         z0 = - fgam * zfunct
         z1 = y / abs(y) * fgam * (1.0 - fgam * zeta * zfunct)
         z2 =  fgam**2 * zeta * (1.0 - fgam * zeta * zfunct)                     
         
      end if
      
      return
      
      zetai_giv = 1.0 / real(zeta) 
                 
      call intplt_z(zetai_giv, dKdL_giv, z0_out, nmax, mmax, z0_table, 
     .   ntable, mtable, zetai_table, dKdL_table)
     
      call intplt_z(zetai_giv, dKdL_giv, z1_out, nmax, mmax, z1_table, 
     .   ntable, mtable, zetai_table, dKdL_table)
     
      call intplt_z(zetai_giv, dKdL_giv, z2_out, nmax, mmax, z2_table, 
     .   ntable, mtable, zetai_table, dKdL_table)                     
     
      if (z0_out .ne. 0.0) z0 = z0_out
      if (z1_out .ne. 0.0) z1 = z1_out
      if (z2_out .ne. 0.0) z2 = z2_out
                  
      return
      end

!
!********************************************************************
!
      subroutine z_integral(sgn_kprl, zeta, delta, z2)

      implicit none

      real dx, xmax, xnueff, sgn_kprl, xmax0
      integer j, npts, nmaxdim
      real x(10000), zeta_mod, delta
      complex z0, z1, z2, zi, f(10000), zeta, beta
      complex cexpx(10000)

      parameter (nmaxdim = 10000)

      zi = cmplx(0.,1.)
      xnueff = 0.0
      zeta_mod = sqrt(conjg(zeta) * zeta)


      if(zeta_mod .gt. 1.0e+03)then
         z2 = - 0.5 / zeta - 0.75 / zeta**3
         return
      else

         xmax0 = 7.0
         npts = 10000
         xmax = sgn_kprl * xmax0
         dx = xmax / (npts - 1)

*        ------------
*        calculate Z2
*        ------------    

         do  j = 1, npts
            x(j) = (j-1) * dx     
            beta = 1. / sqrt(1. - zi * delta * x(j)**2 )            
            cexpx(j) = exp(zi * zeta * x(j) - x(j)**2 / 4. * beta**2)
            f(j) = zi / 2.0 * (1.0 -  x(j)**2 / 2. * beta**2)
     .           * cexpx(j) * beta**3

         end do

         call dgrat1(x, f, 1, npts, z2, nmaxdim)         



         return

      end if

      end


c
c*******************************************************************************
c



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
c        call zfun_hil(zetat, zfunct, flag)
         
         
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
c        call zfun_hil(zetat, zfunct, flag)      
                 
         z0 = - fgam * zfunct
         z1 = y / abs(y) * fgam * (1.0 - fgam * zeta * zfunct)
         z2 =  fgam**2 * zeta * (1.0 - fgam * zeta * zfunct)
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




