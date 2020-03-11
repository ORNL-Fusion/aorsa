
c
c***************************************************************************
c
      subroutine invert(rho_in, theta_in, capr_out, capz_out,
     . nnodex, nnodey, nxmx, nymx, rho, theta, capr, capz, myid,
     . drhodr, drhodz, dthedr, dthedz, xj)
     
      implicit none
      integer myid, nnodex, nnodey, nxmx, nymx, n, m, nmin, mmin
      real rho_in, theta_in, capr_out, capz_out, diff_rho, diff_theta, 
     .   diff, diffmin
      real  capr(nxmx), capz(nymx), rk, zk,
     .     rho(nxmx, nymx), theta(nxmx, nymx), 
     .  drhodr(nxmx, nymx), drhodz(nxmx, nymx), 
     .  dthedr(nxmx, nymx), dthedz(nxmx, nymx)
      real  xj, thetak, rhok, dthedrk, drhodzk, drhodrk, dthedzk
    
      
      diffmin = 1.0
      
      
      do n = 1, nnodex
         do m = 1, nnodey
	    diff_rho = rho_in - rho(n,m)
	    diff_theta = theta_in - theta(n,m)
	    diff = sqrt(diff_rho**2 + diff_theta**2)
	    
	    if(diff .le. diffmin)then
	       diffmin = diff
	       nmin = n
	       mmin = m
	    end if 
	    	 	 
	 end do 
      end do
            
      n = nmin
      m = mmin
      
      rk = capr(n)
      zk = capz(m)
      rhok = rho(n,m)
      thetak = theta(n,m)
      drhodrk = drhodr(n,m)
      drhodzk = drhodz(n,m)
      dthedrk = dthedr(n,m)
      dthedzk = dthedz(n,m)
      
      xj = dthedrk * drhodzk - drhodrk * dthedzk
      
      if(xj .eq. 0)then
         capr_out = rk
	 capz_out = zk
	 return
      end if
      
      capr_out = rk + 1. / xj * (drhodzk * (theta_in - thetak)
     .                         - dthedzk * (rho_in   - rhok)  )
      capz_out = zk + 1. / xj * (dthedrk * (rho_in   - rhok) 
     .                         - drhodrk * (theta_in - thetak))
      
      return
      end


c
c***************************************************************************
c
      subroutine checkx(myid, i, j, nkx1, nkx2, nky1, nky2,
     1   fdk, fek, ffk, u, v, w,
     1   xb, xc, xd, nkdim1, nkdim2, mkdim1, mkdim2,
     1   nnodex, nnodey, nxmx, nymx, xprime, yprime,
     1   idiag, jdiag, isweep)

      implicit none

      integer nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     1        nnodex, nnodey, nxmx, nymx, idiag, jdiag, isweep,
     1        i, j, n, m, myid

      complex fdk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        fek(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        ffk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        xb(nxmx, nymx), term1,
     1        xc(nxmx, nymx), term2,
     1        xd(nxmx, nymx), term3,
     1        xlhs, xrhs
      complex u(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        v(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        w(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real    xprime(nxmx), yprime(nymx)

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



         term1 = 0.0
         term2 = 0.0
         term3 = 0.0

         do n = nkx1, nkx2
            do m = nky1, nky2

            term1 = term1 + fdk(n, m) * u(n, m)
            term2 = term2 + fek(n, m) * v(n, m)
            term3 = term3 + ffk(n, m) * w(n, m)

            end do
         end do

         xlhs = term1 + term2 + term3
         xrhs = xb(i, j)
         write(15, 900) i, j, xprime(i), yprime(j), xlhs, xrhs
         write(6,  900) i, j, xprime(i), yprime(j), xlhs, xrhs
	 
	 write(15, 901) term1, term2, term3
         write(6,  901) term1, term2, term3
	 
	 write(15, 901) u(64,64), v(64,64), w(64,64)
         write(6,  901) u(64,64), v(64,64), w(64,64)
	 
	 write(15, 901) fdk(64,64), fek(64,64), ffk(64,64)
         write(6,  901) fdk(64,64), fek(64,64), ffk(64,64)	 
	 
	 
      endif


      return
  900 format(i5, i5, 1p9e12.3)
  901 format(10x, 1p9e12.3)
  940 format(1h , 4h   i, 5h    j, 4x, 6h x(i) , 6x, 6h y(j) ,
     1                    6x,  6hre lhs, 6x, 6him lhs, 6x,
     1                         6hre rhs, 6x, 6him rhs)
  910 format (1h1)
  920 format (1h0)
  930 format(3x, 14h u diagnostics)
      end

c
c***************************************************************************
c
      subroutine checky(myid, i, j, nkx1, nkx2, nky1, nky2,
     1   fgk, fak, fpk, u, v, w,
     1   xb, xc, xd, nkdim1, nkdim2, mkdim1, mkdim2,
     1   nnodex, nnodey, nxmx, nymx, xprime, yprime,
     1   idiag, jdiag, isweep)

      implicit none

      integer nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     1        nnodex, nnodey, nxmx, nymx, idiag, jdiag, isweep,
     1        i, j, n, m, myid

      complex fgk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        fak(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        fpk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        xb(nxmx, nymx), term1,
     1        xc(nxmx, nymx), term2,
     1        xd(nxmx, nymx), term3,
     1        xlhs, xrhs
      complex u(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        v(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        w(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real    xprime(nxmx), yprime(nymx)

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


         term1 = 0.0
         term2 = 0.0
         term3 = 0.0

         do n = nkx1, nkx2
            do m = nky1, nky2

            term1 = term1 + fgk(n, m) * u(n, m)
            term2 = term2 + fak(n, m) * v(n, m)
            term3 = term3 + fpk(n, m) * w(n, m)

            end do
         end do

         xlhs = term1 + term2 + term3
         xrhs = xc(i, j)
         write(15, 900) i, j, xprime(i), yprime(j), xlhs, xrhs
         write(6,  900) i, j, xprime(i), yprime(j), xlhs, xrhs
	 
         write(15, 901) term1, term2, term3
         write(6,  901) term1, term2, term3
	 
	 write(15, 901) u(64,64), v(64,64), w(64,64)
         write(6,  901) u(64,64), v(64,64), w(64,64)
	 
	 write(15, 901) fgk(64,64), fak(64,64), fpk(64,64)
         write(6,  901) fgk(64,64), fak(64,64), fpk(64,64)
      endif


      return
  900 format(i5, i5, 1p9e12.3)
  901 format(10x, 1p9e12.3)
  940 format(1h , 4h   i, 5h    j, 4x, 6h x(i) , 6x, 6h y(j) ,
     1                    6x,  6hre lhs, 6x, 6him lhs, 6x,
     1                         6hre rhs, 6x, 6him rhs)
  910 format (1h1)
  920 format (1h0)
  930 format(3x, 14h v diagnostics)
      end

c
c***************************************************************************
c
      subroutine checkz(myid, i, j, nkx1, nkx2, nky1, nky2,
     1   frk, fqk, fsk, u, v, w,
     1   xb, xc, xd, nkdim1, nkdim2, mkdim1, mkdim2,
     1   nnodex, nnodey, nxmx, nymx, xprime, yprime,
     1   idiag, jdiag, isweep)

      implicit none

      integer nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2,
     1        nnodex, nnodey, nxmx, nymx, idiag, jdiag, isweep,
     1        i, j, n, m, myid

      complex frk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        fqk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        fsk(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        xb(nxmx, nymx), term1,
     1        xc(nxmx, nymx), term2,
     1        xd(nxmx, nymx), term3,
     1        xlhs, xrhs
      complex u(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        v(nkdim1 : nkdim2, mkdim1 : mkdim2),
     1        w(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real    xprime(nxmx), yprime(nymx)

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



         term1 = 0.0
         term2 = 0.0
         term3 = 0.0

         do n = nkx1, nkx2
            do m = nky1, nky2

            term1 = term1 + frk(n, m) * u(n, m)
            term2 = term2 + fqk(n, m) * v(n, m)
            term3 = term3 + fsk(n, m) * w(n, m)

            end do
         end do

         xlhs = term1 + term2 + term3
         xrhs = xd(i, j)
         write(15, 900) i, j, xprime(i), yprime(j), xlhs, xrhs
         write(6,  900) i, j, xprime(i), yprime(j), xlhs, xrhs
	 
	 write(15, 901) term1, term2, term3
         write(6,  901) term1, term2, term3
	 
	 write(15, 901) u(64,64), v(64,64), w(64,64)
         write(6,  901) u(64,64), v(64,64), w(64,64)
	 
	 write(15, 901) frk(64,64), fqk(64,64), fsk(64,64)
         write(6,  901) frk(64,64), fqk(64,64), fsk(64,64)	 
	 
      endif

      return
  900 format(i5, i5, 1p9e12.3)
  901 format(10x, 1p9e12.3)
  940 format(1h , 4h   i, 5h    j, 4x, 6h x(i) , 6x, 6h y(j) ,
     1                    6x,  6hre lhs, 6x, 6him lhs, 6x,
     1                         6hre rhs, 6x, 6him rhs)
  910 format (1h1)
  920 format (1h0)
  930 format(3x, 14h w diagnostics)
      end

c
c***************************************************************************
c
