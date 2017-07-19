c
c**********************************************************************
c
      subroutine sigslo(q, zions, velect, valpha, omgc, omgp2, omgrf,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz,
     .   xkprl, xkperp, lmin, lmax)

      implicit none

      integer lmin, lmax

      real q, zeff, velect, valpha, omgc, omgp2, omgrf, xkprl, xkperp
      real beta, alfa, xk0, chrg, xkz, rdz, vm, vte, zions, qe,
     .     clight, xmu0, eps0

      complex zi, xkx2,
     .   sigxx, sigxy, sigxz,
     .   sigyx, sigyy, sigyz,
     .   sigzx, sigzy, sigzz

      complex sig(3,3), rdx

      common/slowdist/zeff, vte

      xkz = xkprl
      xkx2 = xkperp**2
      zeff = zions

      qe = -1.6e-19
      eps0 = 8.85e-12
      xmu0 = 1.26e-06
      clight = 1.0 / sqrt(eps0 * xmu0)
      zi = cmplx(0.0, 1.0)
      xk0 = omgrf / clight

      alfa = omgp2 / omgrf**2
      beta = (omgc / omgrf)**2

      chrg = q / abs(qe)

      vm =  valpha / clight
      vte = velect / clight

      rdx = csqrt(xkx2) / xk0
c      rdx = xkperp / xk0
      rdz = xkz / xk0


c      write(6, 11)lmin, lmax
c      write(6, 10)alfa, beta, chrg, vm, rdx, rdz


      call isosig(alfa, beta, chrg, vm, rdx, rdz, lmin, lmax, sig)

      sigxx = sig(1,1) * eps0 * omgrf / zi
      sigxy = sig(1,2) * eps0 * omgrf / zi
      sigxz = sig(1,3) * eps0 * omgrf / zi

      sigyx = sig(2,1) * eps0 * omgrf / zi
      sigyy = sig(2,2) * eps0 * omgrf / zi
      sigyz = sig(2,3) * eps0 * omgrf / zi

      sigzx = sig(3,1) * eps0 * omgrf / zi
      sigzy = sig(3,2) * eps0 * omgrf / zi
      sigzz = sig(3,3) * eps0 * omgrf / zi

c      write(6, 10)sigxx, sigxy, sigxz
c      write(6, 10)sigyx, sigyy, sigyz
c      write(6, 10)sigzx, sigzy, sigzz


      return

   10 format(8e12.4)
   11 format(8i10)

      end




c
c**********************************************************************
c

      subroutine isosig(alfa,beta,chrg,vm,rdx,rdz,nhmin,nhmax,sig)


*     ------------------------------------------------------------------------
*     calculates Hermitian part of dispersion tensor for arbitrary isotropic
*     distribution function of a single plasma species
*
*     sig(i,j) = hermitian part of 4*pi*i/omega*sigma
*
*     alfa=(plasma freq/wave freq)**2 for the given species
*     beta=(cyclotron freq/wave freq)**2
*     chrg=species charge/abs(electron charge) only sgn of chrg is used
*     vmax=maximum velocity of nonzero distribution/speed of light
*     rdx=complex perpendicular refractive index
*     rdz=real parallel refractive index (Plemelj relation assumes kz real)
*     nhmin=minimum harmonic number allowed
*     nhmax=maximum harmonic number allowed
*     note: nmin and nmax are calculated from presence of resonant particles.
*     however, nhmin and nhmax place limits on allowed values
*------------------------------------------------------------------------


      dimension bk(3)
      complex rdx,rx,rx2
      complex sig(3,3),sv(6),v(6),vc(6)
      common /int/rx,rx2,b,bet,vmax,rz,rz2,u,n
c
      external intz1
      data ncall/0/,pi/3.141592654/
c
c
      ncall=ncall+1
      vmax=vm
      bet=beta
      b=chrg*sqrt(beta)/abs(chrg)
      rx=rdx
      rx2=rx**2
      rz=abs(rdz)
      rz2=rz**2
c
      do 10 i=1,6
10    sv(i)=cmplx(0.,0.)
c
c do special case nz=0
c
      if (rz .eq.0.) go to 300
c
c
calculate contribution from each cyclotron harmonic
c
c
c
      do 200 n=nhmin,nhmax
c
c determine if resonant particles exist for this harmonic.  If not
c there is no contribution.
c
      if (abs(1.-n*b).ge. rz*vmax) go to 200
c
      x0=0.
      x1=1.
      np=30
      ndim=6
      call simpc (x0,x1,np,intz1,v,ndim)
c
calculate contribution from velocity cutoff
c

c
calculate u resonant
c
      u=(1.-n*b)/vmax/rz
c
c if u is too big there are no resonant particles, u-parallel integral is 0
c
      x=sqrt(1.-u**2)
      call mofx (x,vc)
      call dgdv2(x,dg,gcut)
c
c

      do 100 i=1,6
100   sv(i)=sv(i)+v(i)-gcut*vc(i)/rz/2.

200   continue
c
c
calculate normalizaton
c
      anorm=-4.*pi**2*alfa
c
240   do 250 i=1,6
250   sv(i)=anorm*cmplx(0.,1.)*sv(i)
c
300   continue
      sig(1,1)=sv(1)
      sig(1,2)=sv(2)
      sig(1,3)=sv(3)
      sig(2,1)=-sig(1,2)
      sig(2,2)=sv(4)
      sig(2,3)=sv(5)
      sig(3,1)=sig(1,3)
      sig(3,2)=-sig(2,3)
      sig(3,3)=sv(6)
c
c
c
c      call matout (sig,'ISOSIG',beta,'BETA',10)
c
      return
      end


c
c**********************************************************************
c

      subroutine intz1 (x,v)

*     ---------------------
*     integrand fo nz=0 case 1
*     ---------------------
      complex rx,rx2
      complex v(6),vm(6)
      common /int/rx,rx2,b,bet,vmax,rz,rz2,u,n
c
calculate u resonant
c
      u=(1.-n*b)/vmax/rz
c
c if u is too big there are no resonant particles, u-parallel integral is 0
c
      if (abs(u).lt.sqrt(1.-x**2)) then
      ucut=1.
      call mofx (x,vm)
      call dgdv2(x,dg,gcut)
      else
      ucut=0.
      end if
c
50    do 100 i=1,6
100   v(i)=ucut*x*dg*vm(i)/rz
c
      return
      end


c
c***************************************************************************
c
      subroutine mofx (x,um)

*     ------------------
*     generate m matrix
*     ------------------

      complex rx,rx2,bj,bjp,besc(50),arg
      dimension bes(3)
      complex um(6),ci
      common /int/ rx,rx2,b,bet,vmax,rz,rz2,u,n
      ci=cmplx(0.,1.)
      arg=rx*x*vmax/b

*     ----------------------------
*     generate j(n,x) and jprime(n,x)
*     ----------------------------

      nabs=iabs(n)
      if (nabs.ge.1) go to 220
      a1=0.
      nj=2
        nmax=1
        call besjc(arg,nmax,besc,ier)
        bj=besc(1)
        bjp=-besc(2)
        if (ier .ne. 0) write(59,10) ier
10      format( ' ier in besjc is ',i2)
      go to 225
220   a1=nabs-1
      nj=3
        nmax=nabs+1
                                call besjc(arg,nmax,besc,ier)
        bjp=(besc(nabs)-besc(nabs+2))/2.
        bj=besc(nabs+1)
        if(ier .ne. 0) write(59,10)ier
225   continue
c
c evaluate um vectors, consisting of the 6 independant
c elements of matrix m
c
c
      um(1)=n**2*b**2/rx2*bj**2
      um(2)=ci*n*b/rx*x*vmax*bj*bjp
      um(3)=u*n*b/rx*vmax*bj**2
      um(4)=x**2*vmax**2*bjp**2
      um(5)=-ci*x*u*vmax**2*bj*bjp
      um(6)=u**2*vmax**2*bj**2
      return
      end


c
c*********************************************************************
c


      subroutine dgdv2 (x,dg,gcut)

*     ---------------------------------------------------------
*     c**3*vmax**2* dg/dv**2 for a slowing down distribution
*     vte=electron thermal speed /c
*     rmass=species mass (probably alpha particles)/electron mass
*     vmax=cutoff speed/c
*     zeff= z effective
*     ---------------------------------------------------------

      complex rx,rx2
      common /int/ rx,rx2,b,bet,vmax,rz,rz2,u,n
      common /slowdist/ zeff, vte
c     data pi/3.1415926/, pi32/5.568328/,rmass/7296./,vmax/0.137/
      data pi/3.1415926/, pi32/5.568328/,rmass/7296./

*     ---------------------------------------------------------
*     calculate vcrit,delt and normalization factor anorm
*     recalculate only if vte or zeff has changed since last time
*     ---------------------------------------------------------

      if ((vte.eq.vte0).and.(zeff.eq.zeff0)) go to 50
      vte0=vte
      zeff0=zeff
      vc=(3.*sqrt(pi)/rmass*zeff)**.3333*vte
      delt=vc/vmax
      a=3./4./pi/alog(1.+1./delt**3)/vc**3
      anorm=1.5*a/delt**3
      gcut0=a/(1.+1./delt**3)
c
50    vu=sqrt(u**2+x**2)
      dg=-anorm*vu/(1.+(vu/delt)**3)**2
      gcut=gcut0
c
      return
      end


c
c***************************************************************************
c

      subroutine dgdv2max (x,dg)


*     -----------------------------------------------------
*     c**3*vmax**2* dg/dv**2 for a Maxwellian distribution

*     t=temperature in kev
*     rmass=species mass/electron mass
*     vmax=cutoff speed/c
*     v=species thermal speed (sqrt(2T/M)) divided by speed of light c
*     -----------------------------------------------------------------


      complex rx,rx2
      common /int/ rx,rx2,b,bet,vmax,rz,rz2,u,n
      common /maxdist/ v
      data pi32/5.568328/
c
      delt=v/vmax
      dg1=-exp(-(u**2+x**2)/delt**2)/v**3/delt**2/pi32
      dg=dg1
      x1=x
c
      return
      end




c
c***************************************************************************
c


      subroutine simpc (x0,x1,np,fname,fint,ndim)

*     -------------------------------------------
*     complex vector integration using simpsons rule
*     ------------------------------------------

      complex fint(ndim),veven(100),vodd(100)

      n=2*np
      h=(x1-x0)/n

      call fname (x0,fint)

      do 50 i=1,np
      x=x0+2.*h*i
      call fname (x,veven)
      x=x-h
      call fname (x,vodd)
      wt=2.
      if (i.eq.np) wt=1.
      do 25 j=1,ndim
25    fint(j)=fint(j)+wt*veven(j)+4.*vodd(j)
50    continue
c
      do 100 i=1,ndim
      fint(i)=h*fint(i)/3.
100   continue
      return
      end


