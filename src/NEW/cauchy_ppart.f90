        subroutine cauchy_ppart2(x,nx,xres,fx,pint)
	
!     -----------------------------------------------
!     Returns principal part of int(x,fx/(x-xres))
!     -----------------------------------------------

        implicit none
        integer, parameter :: idebug = 0

        integer   nx
        real xres 
        real x(nx), fx(nx)
        real pint
!       ---------------
!       local variables
!       ---------------
        logical, parameter :: use_linear = .true.

        integer :: i
        real, dimension(nx-1) :: va0,va1,va2
        real :: a0,a1,a2,  x1,x2,x3,  f1,f2,f3
        real :: tmp1,tmp2,  d1,d2,d3,  d31,d21
        real :: value0,value1,value2, xleft,xright
        real :: v0,v1,v2
        intrinsic :: log


!       ------------------------------------------------------------
!       compute piecewise quadratic over 
!       [x(i-1), x(i), x(i+1)] as    a0 + a1*(x-xres)+ a2*(x-xres)^2
!
!       assume xres is not exactly one of x(1:n)
!       ------------------------------------------------------------
        do i=2,nx-1
           x1 = x(i-1)
           x2 = x(i)
           x3 = x(i+1)
           f1 = fx(i-1)
           f2 = fx(i)
           f3 = fx(i+1)
!          --------------------------------------
!          find a0,a1,a2 such that
!
!          a0 + a1*(x1-xres) + a2*(x1-xres)^2 = f1
!          a0 + a1*(x2-xres) + a2*(x2-xres)^2 = f2
!          a0 + a1*(x3-xres) + a2*(x3-xres)^2 = f3
!          --------------------------------------
           d1 = x1-xres
           d2 = x2-xres
           d3 = x3-xres
           d31 = (d3-d1)*(d3+d1)
           d21 = (d2-d1)*(d2+d1)
!          ----------------------------------
!          a1 * ( d3-d1 ) + a2*( d31 ) = f3-f1
!          a1 * ( d2-d1 ) + a2*( d21 ) = f2-f1
!          -----------------------------------
           tmp1 = ((f3-f1)*(d2-d1) - (f2-f1)*(d3-d1))
           tmp2 = ( d31*(d2-d1) - d21*(d3-d1) )

           a2 = tmp1/tmp2
           a1 = ((f3-f1) - a2*d31 )/(d3-d1)
           a0 = f1 - a2*d1*d1 - a1*d1

           va0(i) = a0
           va1(i) = a1
           va2(i) = a2
       enddo

!       --------------------------------------------
!       special handling for first interval
!       --------------------------------------------
       if (use_linear) then
!         ------------------------------
!         backward compatible with df_dv
!
!         a0 + a1*(x1-xres) = f1
!         a0 + a1*(x2-xres) = f2
!         ------------------------------
               x1 = x(1)
               x2 = x(2)
               f1 = fx(1)
               f2 = fx(2)

               a1 = (f2-f1)/(x2-x1)
               a0 = f1 - a1*(x1-xres)
               a2 = 0.0

               va0(1) = a0
               va1(1) = a1
               va2(1) = a2
       else
!        ------------------------------------
!        use quadratic even in first interval
!        ------------------------------------
               va0(1) = va0(2)
               va1(1) = va1(2)
               va2(1) = va2(2)

       endif

       if (idebug.ge.1) then
         do i=1,nx-1
          write(*,*) 'i,a0,a1,a2 ',i,va0(i),va1(i),va2(i)
         enddo
       endif


!      ----------------------------------------------
!      sum contribtions from a0
!
!      int( 1/(x-xres), x=xleft..xright ) = log( (xright-xres)/(xleft-xres) )
!
!      sum contributions from a1
!
!      int( a1, x=xleft..xright) = (xright-xleft)*a1
!
!      sum contributions from a2
!
!      int( a2*(x-xres), x=xleft..xright) = 
!      = (a2/2) * ( (xright-xres)**2 - (xleft-xres)**2)
!      = (a2/2) * ( (xright+xleft - 2*xres) * (xright-xleft) )
!      ----------------------------------------------
       value0 = 0.0d0
       value1 = 0.0d0
       value2 = 0.0d0

       do i=1,nx-1
          a0 = va0(i)
          a1 = va1(i)
          a2 = va2(i)

          xleft = x(i)
          xright = x(i+1)

          v0 = a0*log( abs( (xright-xres)/(xleft-xres) ) )
          v1 = a1*(xright-xleft)
          v2 = a2*(xright+xleft-2.0d0*xres)*(xright-xleft)

          value0 = value0 + v0
          value1 = value1 + v1
          value2 = value2 + v2

          if (idebug.ge.1) then
            write(*,*) 'i,v0,v1,v2 ',i,v0,v1,v2
          endif

        enddo 
        value2 = value2*0.5d0


        pint = value0 + value1 + value2
          
        return
        end subroutine cauchy_ppart2
	
	
!
!*************************************************************************
!


        subroutine cauchy_ppart6(x,nx,xres,nfx,fx,pint,is_uniform)
!     -----------------------------------------------
!     Returns principal part of int(x,fx/(x-xres))
!     -----------------------------------------------
        implicit none
        integer, parameter :: idebug = 0

        integer   nx,nfx
        real*8 xres 
        real*8 x(nx), fx(nx,nfx)
        real*8 pint(nfx)
        logical  :: is_uniform
!       ---------------
!       local variables
!       ---------------

        integer :: i,j,ih,nlog
        real*8, dimension(nx) :: vx, vlogt
        real*8 :: a0,a1,a2,  x1,x2,x3,  f1,f2,f3
        real*8 :: tmp1,tmp2,  d1,d2,d3,  d31,d21
        real*8 :: d3md1, d2md1
        real*8, dimension(nfx) :: value0,value1,value2 
        real*8 :: xleft,xright
        real*8 :: v0,v1,v2
        real*8 :: xc,a,b

        intrinsic :: log,abs,mod

        logical, parameter ::use_approx = .false.


        logical, parameter :: use_pade = .true.
        
        real*8 ::  h,t,dx,t2, logt
        real*8, parameter :: zero=0.0d0
        real*8, parameter :: one=1.0d0
        real*8, parameter :: two=2.0d0
        real*8, parameter :: three=3.0d0
        real*8, parameter :: four=4.0d0
        real*8, parameter :: five=5.0d0
        real*8, parameter :: six=6.0d0
        real*8, parameter :: seven=7.0d0
        real*8, parameter :: eight=8.0d0
        logical :: isok


!       ------------------------------------------------------------
!       compute piecewise quadratic over 
!       [x(i-1), x(i), x(i+1)] as    a0 + a1*(x-xres)+ a2*(x-xres)^2
!
!       assume xres is not exactly one of x(1:n)
!       ------------------------------------------------------------

        value0 = zero
        value1 = zero
        value2 = zero
        pint = zero

!       ------------------
!       precompute the logs
!       ------------------
        nlog = (nx-1)/2
        do ih=1,nlog
           i = 2*ih
           x1 = x(i-1)
           x2 = x(i)
           x3 = x(i+1)
           vx(ih) = abs( (x3-xres)/(x1-xres) )
        enddo
        call vlog( vlogt, vx, nlog )
           

        if (is_uniform) then

        dx = (x(nx)-x(1))/dble(nx-1)
        h = dx


     
        do j=1,nfx
        do ih=1,nlog
           i = 2*ih

           x1 = x(i-1)
           x2 = x(i)
           x3 = x(i+1)

             d2 = x2-xres
             d1 = d2 - h
             d3 = d2 + h



           f1 = fx(i-1,j)
           f2 = fx(i,j)
           f3 = fx(i+1,j)

            
!          --------------------------------------
!          find a0,a1,a2 such that
!
!          a0 + a1*(x1-xres) + a2*(x1-xres)^2 = f1
!          a0 + a1*(x2-xres) + a2*(x2-xres)^2 = f2
!          a0 + a1*(x3-xres) + a2*(x3-xres)^2 = f3
!
!          a0 + a1*(x2-xres) + a2*(x2-xres)^2 = f2
!          a0 + a1*((x2-xres)-h) + a2*((x2-xres)-h)^2 = f1
!          a0 + a1*((x2-xres)+h) + a2*((x2-xres)+h)^2 = f3
!
!          f2 + a1*h + a2*(h*h + 2*h*(x2-xres)) = f3
!          f2 - a1*h + a2*(h*h - 2*h*(x2-xres)) = f1
!          --------------------------------------



!             ---------------------
!             d31 = (d3-d1)*(d3+d1)
!             d21 = (d2-d1)*(d2+d1)
!
!          a1 * ( d3-d1 ) + a2*( d31 ) = f3-f1
!          a1 * ( d2-d1 ) + a2*( d21 ) = f2-f1
!          -----------------------------------

             a2 = (f3-two*f2+f1)/(two*h*h)
             a1 = (f3-f1)/(two*h) - two*a2*d2
             a0 = f2 - (a1*d2 + a2*d2*d2)


             logt = vlogt(ih)
             value0(j) = value0(j) + a0*logt
             value1(j) = value1(j) + a1 
             value2(j) = value2(j) + a2*d2

            enddo
            enddo

            value1 = value1 * (two*h)
            value2 = value2 * (four*h)


        else


        
        do j=1,nfx
        do ih=1,nlog
           i = 2*ih
           x1 = x(i-1)
           x2 = x(i)
           x3 = x(i+1)


            
!          --------------------------------------
!          find a0,a1,a2 such that
!
!          a0 + a1*(x1-xres) + a2*(x1-xres)^2 = f1
!          a0 + a1*(x2-xres) + a2*(x2-xres)^2 = f2
!          a0 + a1*(x3-xres) + a2*(x3-xres)^2 = f3
!          --------------------------------------
             d1 = x1-xres
             d2 = x2-xres
             d3 = x3-xres
             d31 = (d3-d1)*(d3+d1)
             d21 = (d2-d1)*(d2+d1)
!          ----------------------------------
!          a1 * ( d3-d1 ) + a2*( d31 ) = f3-f1
!          a1 * ( d2-d1 ) + a2*( d21 ) = f2-f1
!          -----------------------------------
             logt = log( abs( (x3-xres)/(x1-xres) ))

             f1 = fx(i-1,j)
             f2 = fx(i,j)
             f3 = fx(i+1,j)

             tmp1 = ((f3-f1)*(d2-d1) - (f2-f1)*(d3-d1))
             tmp2 = ( d31*(d2-d1) - d21*(d3-d1) )

             a2 = tmp1/tmp2
             a1 = ((f3-f1) - a2*d31 )/(d3-d1)
             a0 = f1 - a2*d1*d1 - a1*d1
!
!         f(x) = a0 + a1*(x-xres) + a2*(x-xres)^2           
!
!         int( f(x)/(x-xres), x=x1..x3 ) = 
!         a0 * int( 1/(x-xres), x=x1..x3 ) + 
!         a1 * int( 1, x=x1..x3)
!         a2 * int( (x-xres), x=x1..x3 )
!

             logt = vlogt(ih)
             value0(j) = value0(j) + a0 * logt
             value1(j) = value1(j) + a1 * (x3-x1)
             value2(j) = value2(j) + a2 * (x3-x1)*((x3-xres)+(x1-xres))

              enddo

          enddo
          endif


!         ----------------------------------
!         do we have an odd number of intervals?
!         if so, need to account for last interval
!         ----------------------------------
          if (mod(nx,2).eq.0) then

           x1 = x(nx-2)
           x2 = x(nx-1)
           x3 = x(nx)

             d1 = x1-xres
             d2 = x2-xres
             d3 = x3-xres
             d31 = (d3-d1)*(d3+d1)
             d21 = (d2-d1)*(d2+d1)

             logt = log( abs( (x3-xres)/(x2-xres) ))

           do j=1,nfx

           f1 = fx(nx-2,j)
           f2 = fx(nx-1,j)
           f3 = fx(nx,j)

            
!          --------------------------------------
!          find a0,a1,a2 such that
!
!          a0 + a1*(x1-xres) + a2*(x1-xres)^2 = f1
!          a0 + a1*(x2-xres) + a2*(x2-xres)^2 = f2
!          a0 + a1*(x3-xres) + a2*(x3-xres)^2 = f3
!          --------------------------------------
!          ----------------------------------
!          a1 * ( d3-d1 ) + a2*( d31 ) = f3-f1
!          a1 * ( d2-d1 ) + a2*( d21 ) = f2-f1
!          -----------------------------------
             tmp1 = ((f3-f1)*(d2-d1) - (f2-f1)*(d3-d1))
             tmp2 = ( d31*(d2-d1) - d21*(d3-d1) )

             a2 = tmp1/tmp2
             a1 = ((f3-f1) - a2*d31 )/(d3-d1)
             a0 = f1 - a2*d1*d1 - a1*d1
!
!         f(x) = a0 + a1*(x-xres) + a2*(x-xres)^2           
!
!         int( f(x)/(x-xres), x=x2..x3 ) = 
!         a0 * int( 1/(x-xres), x=x2..x3 ) + 
!         a1 * int( 1, x=x1..x3)
!         a2 * int( (x-xres), x=x2..x3 )
!
             value1(j) = value1(j) + a1 * (x3-x2)
             value2(j) = value2(j) + a2 * (x3-x2)*(x3+x2-two*xres)
             value0(j) = value0(j) + a0 * logt
           enddo

           endif
            
          value2 = value2 / two
          pint = value0 + value1 + value2


          return

        end subroutine cauchy_ppart6



