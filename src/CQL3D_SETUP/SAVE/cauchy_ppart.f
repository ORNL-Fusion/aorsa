      subroutine CAUCHY_PPART2(X, NX, XRES, FX, PINT)

!     -----------------------------------------------
!     Returns principal part of int(x,fx/(x-xres))
!     -----------------------------------------------

      implicit none

      integer, intent(IN):: NX
      real, dimension(NX), intent(IN):: X, FX
      real, intent(IN):: XRES
      real, intent(OUT):: PINT

      real:: CXP,CXM,FXP,FXM,LLOGF,XMAX,XMIN
      real:: CXF,DXF,DX
      real, dimension(NX):: LXF,INTGD,FXPP
      integer:: IX, IXM1, IXP1

      integer lx
      real fxa(1000), xa(1000)

!     ----------------
!     input parameters
!     ----------------
	external df_dv !function
      real :: a,b !uper, lower limits
	real :: c !location of singular value
      real :: epsabs = 0.1 !absolute accuracy
      real :: epsrel = 0.1  !relative accuracy\
      integer, parameter :: limit = 2000  !max number of function evaluations

!     ------------------
!     return parameters	
!     ------------------
      real :: result !approximation to integral
      real :: abserr !estimate of modulus of the absolute error
	integer :: neval !number of integrand evaluations
      integer :: ier  !error code on return 0 = normal
	real, dimension(limit) ::  alist, blist, rlist, elist  !integral stuff
      integer, dimension(limit) :: iord !pointer to error stuff
	integer :: last !actual number of subintervals

      real cauchy_cb(2002)
      

      lx = nx
      cauchy_cb(2002) = nx
!dir$ concurrent
      do ix = 1, nx
         cauchy_cb(1000+ix) = x(ix)
         cauchy_cb(ix) = fx(ix)
      end do
      
      a = x(1)
      b = x(nx)
      c = xres

      dx = (b - a) / (lx - 1)
      cauchy_cb(2001) = dx
      call qawce(cauchy_cb, a, b, c, epsabs, epsrel, limit, 
     &           pint, abserr, neval, ier, alist, blist, rlist, 
     &           elist, iord, last)



      end subroutine CAUCHY_PPART2

!
!*************************************************************************
!

      function df_dv(x,cauchy_cb)

      implicit none

      real :: df_dv
      real :: x

      integer lx, i
      real  dx
      real a1, a2, a3, p, p2
      real :: cauchy_cb(2002)
      integer,parameter:: ifxa=0, ixa=1000, idx=2001, ilx=2002


      dx=cauchy_cb(idx)
      lx = cauchy_cb(ilx)

      i = int((x-cauchy_cb(ixa+1)) / dx) + 1
      p = (x - cauchy_cb(ixa+i)) / dx
c       go to 2
*      -----------------------
*      quadradic interpolation
*      -----------------------
c       call winterp1d_2(xa, lx, fxa, x, df_dv)

      p2 = p * p

      if(i .ne. 1 .and. i .ne. lx)then
        a1 = 0.5 * (p2 - p)
        a2 = 1. - p2
        a3 = 0.5 * (p2 + p)
        df_dv=  a1*cauchy_cb(ifxa+i-1) +
     &          a2*cauchy_cb(ifxa+i) +
     &          a3*cauchy_cb(ifxa+i+1)

      else
        if(i .eq. lx) then
           df_dv = cauchy_cb(ifxa+i)
        endif
        if(i .eq. 1) then 
           df_dv = cauchy_cb(ifxa+i)+
     &            (cauchy_cb(ifxa+i+1) -
     &             cauchy_cb(ifxa+i))*p
        endif
      end if

      return


 2    continue
*      --------------------
*      linear interpolation
*      --------------------
      if(i .eq. lx)then
        df_dv = cauchy_cb(ifxa+lx)
      else
        df_dv = cauchy_cb(ifxa+i)+
     &         (cauchy_cb(ifxa+i+1)-
     &          cauchy_cb(ifxa+i))*p
      end if

      return

      end function df_dv
*
*-----------------------------------------------------------------------
*
      subroutine qawce(f_cb,a,b,c,epsabs,epsrel,limit,result,
     &   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  qawce
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a1,j4
c***keywords  automatic integrator, special-purpose,
c             cauchy principal value, clenshaw-curtis method
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***  purpose  the routine calculates an approximation result to a
c              cauchy principal value i = integral of f*w over (a,b)
c              (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying
c              following claim for accuracy
c              abs(i-result).le.max(epsabs,epsrel*abs(i))
c***description
c
c        computation of a cauchy principal value
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            c      - real
c                     parameter in the weight function, c.ne.a, c.ne.b
c                     if c = a or c = b, the routine will end with
c                     ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of
c                             limit. however, if this yields no
c                             improvement it is advised to analyze the
c                             the integrand, in order to determine the
c                             the integration difficulties. if the
c                             position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some interior points of
c                             the integration interval.
c                         = 6 the input is invalid, because
c                             c = a or c = b or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1.
c                             result, abserr, neval, rlist(1), elist(1),
c                             iord(1) and last are set to zero. alist(1)
c                             and blist(1) are set to a and b
c                             respectively.
c
c            alist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the left
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            blist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the right
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            rlist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the integral
c                      approximations on the subintervals
c
c            elist   - real
c                      vector of dimension limit, the first  last
c                      elements of which are the moduli of the absolute
c                      error estimates on the subintervals
c
c            iord    - integer
c                      vector of dimension at least limit, the first k
c                      elements of which are pointers to the error
c                      estimates over the subintervals, so that
c                      elist(iord(1)), ..., elist(iord(k)) with k = last
c                      if last.le.(limit/2+2), and k = limit+1-last
c                      otherwise, form a decreasing sequence
c
c            last    - integer
c                      number of subintervals actually produced in
c                      the subdivision process
c
c***references  (none)
c***routines called  qc25c,qpsrt,r1mach
c***end prologue  qawce
c
      real a,aa,abserr,alist,area,area1,area12,area2,a1,a2,b,bb,blist,
     *  b1,b2,c,r1mach,elist,epmach,epsabs,epsrel,errbnd,errmax,error1,
     *  error2,errsum,f,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,krule,last,limit,maxerr,nev,
     *  neval,nrmax
c
      dimension alist(limit),blist(limit),rlist(limit),elist(limit),
     *  iord(limit)
c
      real f_cb(2002)
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qawce
c      epmach = r1mach(4)
      epmach = 1.0e-30
c      uflow = r1mach(1)
      uflow = 1.0e-30
c
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 6
      neval = 0
      last = 0
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0e+00
      elist(1) = 0.0e+00
      iord(1) = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(c.eq.a.or.c.eq.b.or.(epsabs.le.0.0e+00.and
     *  .epsrel.lt.amax1(0.5e+02*epmach,0.5e-14))) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      aa=a
      bb=b
      if (a.le.b) go to 10
      aa=b
      bb=a
10    ier=0
      krule = 1
      call qc25c(f_cb,aa,bb,c,result,abserr,krule,neval)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      alist(1) = a
      blist(1) = b
c
c           test on accuracy
c
      errbnd = amax1(epsabs,epsrel*abs(result))
      if(limit.eq.1) ier = 1
      if(abserr.lt.amin1(0.1e-01*abs(result),errbnd)
     *  .or.ier.eq.1) go to 70
c
c           initialization
c           --------------
c
      alist(1) = aa
      blist(1) = bb
      rlist(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
c
c           main do-loop
c           ------------
c
      do 40 last = 2,limit
c
c           bisect the subinterval with nrmax-th largest
c           error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5e+00*(alist(maxerr)+blist(maxerr))
        b2 = blist(maxerr)
        if(c.le.b1.and.c.gt.a1) b1 = 0.5e+00*(c+b2)
        if(c.gt.b1.and.c.lt.b2) b1 = 0.5e+00*(a1+c)
        a2 = b1
        krule = 2
        call qc25c(f_cb,a1,b1,c,area1,error1,krule,nev)
        neval = neval+nev
        call qc25c(f_cb,a2,b2,c,area2,error2,krule,nev)
        neval = neval+nev
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(abs(rlist(maxerr)-area12).lt.0.1e-04*abs(area12)
     *    .and.erro12.ge.0.99e+00*errmax.and.krule.eq.0)
     *    iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax.and.krule.eq.0)
     *    iroff2 = iroff2+1
        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = amax1(epsabs,epsrel*abs(area))
        if(errsum.le.errbnd) go to 15
c
c           test for roundoff error and eventually
c           set error flag.
c
        if(iroff1.ge.6.and.iroff2.gt.20) ier = 2
c
c           set error flag in the case that number of interval
c           bisections exceeds limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(amax1(abs(a1),abs(b2)).le.(0.1e+01+0.1e+03*epmach)
     *    *(abs(a2)+0.1e+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
   15   if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine qpsrt to maintain the descending ordering
c           in the list of error estimates and select the
c           subinterval with nrmax-th largest error estimate (to be
c           bisected next).
c
   30    call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 50
   40 continue
c
c           compute final result.
c           ---------------------
c
   50 result = 0.0e+00
      do 60 k=1,last
        result = result+rlist(k)
   60 continue
      abserr = errsum
   70 if (aa.eq.b) result=-result
  999 return
      end

*
*-----------------------------------------------------------------------
*

      subroutine qc25c(f_cb,a,b,c,result,abserr,krul,neval)
c***begin prologue  qc25c
c***date written   810101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a2a2,j4
c***keywords  25-point clenshaw-curtis integration
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (a,b) with
c            error estimate, where w(x) = 1/(x-c)
c***description
c
c        integration rules for the computation of cauchy
c        principal value integrals
c        standard fortran subroutine
c        real version
c
c        parameters
c           f      - real
c                    function subprogram defining the integrand function
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l  in the driver program.
c
c           a      - real
c                    left end point of the integration interval
c
c           b      - real
c                    right end point of the integration interval, b.gt.a
c
c           c      - real
c                    parameter in the weight function
c
c           result - real
c                    approximation to the integral
c                    result is computed by using a generalized
c                    clenshaw-curtis method if c lies within ten percent
c                    of the integration interval. in the other case the
c                    15-point kronrod rule obtained by optimal addition
c                    of abscissae to the 7-point gauss rule, is applied.
c
c           abserr - real
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           krul   - integer
c                    key which is decreased by 1 if the 15-point
c                    gauss-kronrod scheme has been used
c
c           neval  - integer
c                    number of integrand evaluations
c
c***references  (none)
c***routines called  qcheb,qk15w,qwgtc
c***end prologue  qc25c
c
      real a,abserr,ak22,amom0,amom1,amom2,b,c,cc,
     *  centr,cheb12,cheb24,qwgtc,f,fval,hlgth,p2,p3,p4,
     *  resabs,resasc,result,res12,res24,u,x
      integer i,isym,k,kp,krul,neval
c
      dimension x(11),fval(25),cheb12(13),cheb24(25)
c
      external qwgtc
      real f_cb(2002)
c
c           the vector x contains the values cos(k*pi/24),
c           k = 1, ..., 11, to be used for the chebyshev series
c           expansion of f
c
      data x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),
     *  x(11)/
     *     0.9914448613738104e+00,     0.9659258262890683e+00,
     *     0.9238795325112868e+00,     0.8660254037844386e+00,
     *     0.7933533402912352e+00,     0.7071067811865475e+00,
     *     0.6087614290087206e+00,     0.5000000000000000e+00,
     *     0.3826834323650898e+00,     0.2588190451025208e+00,
     *     0.1305261922200516e+00/
c
c           list of major variables
c           ----------------------
c           fval   - value of the function f at the points
c                    cos(k*pi/24),  k = 0, ..., 24
c           cheb12 - chebyshev series expansion coefficients,
c                    for the function f, of degree 12
c           cheb24 - chebyshev series expansion coefficients,
c                    for the function f, of degree 24
c           res12  - approximation to the integral corresponding
c                    to the use of cheb12
c           res24  - approximation to the integral corresponding
c                    to the use of cheb24
c           qwgtc - external function subprogram defining
c                    the weight function
c           hlgth  - half-length of the interval
c           centr  - mid point of the interval
c
c
c           check the position of c.
c
c***first executable statement  qc25c
      cc = (0.2e+01*c-b-a)/(b-a)
      if(abs(cc).lt.0.11e+01) go to 10
c
c           apply the 15-point gauss-kronrod scheme.
c
      krul = krul-1
      call qk15w(f_cb,c,p2,p3,p4,kp,a,b,result,abserr,
     *  resabs,resasc)
      neval = 15
      if (resasc.eq.abserr) krul = krul+1
      go to 50
c
c           use the generalized clenshaw-curtis method.
c
   10 hlgth = 0.5e+00*(b-a)
      centr = 0.5e+00*(b+a)
      neval = 25
      fval(1) = 0.5e+00*df_dv(hlgth+centr,f_cb)
      fval(13) = df_dv(centr,f_cb)
      fval(25) = 0.5e+00*df_dv(centr-hlgth,f_cb)
!dir$ concurrent
      do 20 i=2,12
        u = hlgth*x(i-1)
        isym = 26-i
!        fval(i) = f(u+centr,f_cb)
!        fval(isym) = f(centr-u,f_cb)
        fval(i) = df_dv(u+centr,f_cb)
        fval(isym) = df_dv(centr-u,f_cb)
   20 continue
c
c           compute the chebyshev series expansion.
c
      call qcheb(x,fval,cheb12,cheb24)
c
c           the modified chebyshev moments are computed
c           by forward recursion, using amom0 and amom1
c           as starting values.
c
      amom0 = alog(abs((0.1e+01-cc)/(0.1e+01+cc)))
      amom1 = 0.2e+01+cc*amom0
      res12 = cheb12(1)*amom0+cheb12(2)*amom1
      res24 = cheb24(1)*amom0+cheb24(2)*amom1
      do 30 k=3,13
        amom2 = 0.2e+01*cc*amom1-amom0
        ak22 = (k-2)*(k-2)
        if((k/2)*2.eq.k) amom2 = amom2-0.4e+01/(ak22-0.1e+01)
        res12 = res12+cheb12(k)*amom2
        res24 = res24+cheb24(k)*amom2
        amom0 = amom1
        amom1 = amom2
   30 continue
      do 40 k=14,25
        amom2 = 0.2e+01*cc*amom1-amom0
        ak22 = (k-2)*(k-2)
        if((k/2)*2.eq.k) amom2 = amom2-0.4e+01/
     *  (ak22-0.1e+01)
        res24 = res24+cheb24(k)*amom2
        amom0 = amom1
        amom1 = amom2
   40 continue
      result = res24
      abserr = abs(res24-res12)
   50 return
      end


*
*-----------------------------------------------------------------------
*

      subroutine qcheb(x,fval,cheb12,cheb24)
c***begin prologue  qcheb
c***refer to  qc25c,qc25f,qc25s
c***routines called  (none)
c***revision date  830518   (yymmdd)
c***keywords  chebyshev series expansion, fast fourier transform
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this routine computes the chebyshev series expansion
c            of degrees 12 and 24 of a function using a
c            fast fourier transform method
c            f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
c            f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
c            where t(k,x) is the chebyshev polynomial of degree k.
c***description
c
c        chebyshev series expansion
c        standard fortran subroutine
c        real version
c
c        parameters
c          on entry
c           x      - real
c                    vector of dimension 11 containing the
c                    values cos(k*pi/24), k = 1, ..., 11
c
c           fval   - real
c                    vector of dimension 25 containing the
c                    function values at the points
c                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
c                    where (a,b) is the approximation interval.
c                    fval(1) and fval(25) are divided by two
c                    (these values are destroyed at output).
c
c          on return
c           cheb12 - real
c                    vector of dimension 13 containing the
c                    chebyshev coefficients for degree 12
c
c           cheb24 - real
c                    vector of dimension 25 containing the
c                    chebyshev coefficients for degree 24
c
c***end prologue  qcheb
c
      real alam,alam1,alam2,cheb12,cheb24,
     *  fval,part1,part2,part3,v,x
      integer i,j
c
      dimension cheb12(13),cheb24(25),fval(25),v(12),x(11)
c
c***first executable statement  qcheb
      do 10 i=1,12
        j = 26-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   10 continue
      alam1 = v(1)-v(9)
      alam2 = x(6)*(v(3)-v(7)-v(11))
      cheb12(4) = alam1+alam2
      cheb12(10) = alam1-alam2
      alam1 = v(2)-v(8)-v(10)
      alam2 = v(4)-v(6)-v(12)
      alam = x(3)*alam1+x(9)*alam2
      cheb24(4) = cheb12(4)+alam
      cheb24(22) = cheb12(4)-alam
      alam = x(9)*alam1-x(3)*alam2
      cheb24(10) = cheb12(10)+alam
      cheb24(16) = cheb12(10)-alam
      part1 = x(4)*v(5)
      part2 = x(8)*v(9)
      part3 = x(6)*v(7)
      alam1 = v(1)+part1+part2
      alam2 = x(2)*v(3)+part3+x(10)*v(11)
      cheb12(2) = alam1+alam2
      cheb12(12) = alam1-alam2
      alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8)
     *  +x(9)*v(10)+x(11)*v(12)
      cheb24(2) = cheb12(2)+alam
      cheb24(24) = cheb12(2)-alam
      alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8)
     *  +x(3)*v(10)-x(1)*v(12)
      cheb24(12) = cheb12(12)+alam
      cheb24(14) = cheb12(12)-alam
      alam1 = v(1)-part1+part2
      alam2 = x(10)*v(3)-part3+x(2)*v(11)
      cheb12(6) = alam1+alam2
      cheb12(8) = alam1-alam2
      alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6)
     *  -x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
      cheb24(6) = cheb12(6)+alam
      cheb24(20) = cheb12(6)-alam
      alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8)
     *  -x(9)*v(10)-x(5)*v(12)
      cheb24(8) = cheb12(8)+alam
      cheb24(18) = cheb12(8)-alam
      do 20 i=1,6
        j = 14-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   20 continue
      alam1 = v(1)+x(8)*v(5)
      alam2 = x(4)*v(3)
      cheb12(3) = alam1+alam2
      cheb12(11) = alam1-alam2
      cheb12(7) = v(1)-v(5)
      alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
      cheb24(3) = cheb12(3)+alam
      cheb24(23) = cheb12(3)-alam
      alam = x(6)*(v(2)-v(4)-v(6))
      cheb24(7) = cheb12(7)+alam
      cheb24(19) = cheb12(7)-alam
      alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
      cheb24(11) = cheb12(11)+alam
      cheb24(15) = cheb12(11)-alam
      do 30 i=1,3
        j = 8-i
        v(i) = fval(i)-fval(j)
        fval(i) = fval(i)+fval(j)
   30 continue
      cheb12(5) = v(1)+x(8)*v(3)
      cheb12(9) = fval(1)-x(8)*fval(3)
      alam = x(4)*v(2)
      cheb24(5) = cheb12(5)+alam
      cheb24(21) = cheb12(5)-alam
      alam = x(8)*fval(2)-fval(4)
      cheb24(9) = cheb12(9)+alam
      cheb24(17) = cheb12(9)-alam
      cheb12(1) = fval(1)+fval(3)
      alam = fval(2)+fval(4)
      cheb24(1) = cheb12(1)+alam
      cheb24(25) = cheb12(1)-alam
      cheb12(13) = v(1)-v(3)
      cheb24(13) = cheb12(13)
      alam = 0.1e+01/0.6e+01
      do 40 i=2,12
        cheb12(i) = cheb12(i)*alam
   40 continue
      alam = 0.5e+00*alam
      cheb12(1) = cheb12(1)*alam
      cheb12(13) = cheb12(13)*alam
      do 50 i=2,24
        cheb24(i) = cheb24(i)*alam
   50 continue
      cheb24(1) = 0.5e+00*alam*cheb24(1)
      cheb24(25) = 0.5e+00*alam*cheb24(25)
      return
      end

*
*-----------------------------------------------------------------------
*

      subroutine qk15w(f_cb,p1,p2,p3,p4,kp,a,b,result,abserr,
     *   resabs,resasc)
c***begin prologue  qk15w
c***date written   810101   (yymmdd)
c***revision date  830518   (mmddyy)
c***category no.  h2a2a2
c***keywords  15-point gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  to compute i = integral of f*w over (a,b), with error
c                           estimate
c                       j = integral of abs(f*w) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real version
c
c           parameters
c             on entry
c              f      - real
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              w      - real
c                       function subprogram defining the integrand
c                       weight function w(x). the actual name for w
c                       needs to be declared e x t e r n a l in the
c                       calling program.
c
c              p1, p2, p3, p4 - real
c                       parameters in the weight function
c
c              kp     - integer
c                       key for indicating the type of weight function
c
c              a      - real
c                       lower limit of integration
c
c              b      - real
c                       upper limit of integration
c
c            on return
c              result - real
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule (resg).
c
c              abserr - real
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - real
c                       approximation to the integral of abs(f)
c
c              resasc - real
c                       approximation to the integral of abs(f-i/(b-a))
c
c***references  (none)
c***routines called  r1mach
c***end prologue  qk15w
c
      real a,absc,absc1,absc2,abserr,b,centr,dhlgth,
     *  r1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,
     *  hlgth,p1,p2,p3,p4,resabs,resasc,resg,resk,reskh,result,uflow,
     *  w,wg,wgk,xgk
      integer j,jtw,jtwm1,kp
      external df_dv,qwgtc
      real f_cb(2002)
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(4)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 15-point gauss-kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ... abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point gauss-kronrod rule
c
c           wg     - weights of the 7-point gauss rule
c
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),
     *  xgk(8)/
     *     0.9914553711208126e+00,     0.9491079123427585e+00,
     *     0.8648644233597691e+00,     0.7415311855993944e+00,
     *     0.5860872354676911e+00,     0.4058451513773972e+00,
     *     0.2077849550078985e+00,     0.0000000000000000e+00/
c
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),
     *  wgk(8)/
     *     0.2293532201052922e-01,     0.6309209262997855e-01,
     *     0.1047900103222502e+00,     0.1406532597155259e+00,
     *     0.1690047266392679e+00,     0.1903505780647854e+00,
     *     0.2044329400752989e+00,     0.2094821410847278e+00/
c
      data wg(1),wg(2),wg(3),wg(4)/
     *     0.1294849661688697e+00,    0.2797053914892767e+00,
     *     0.3818300505051889e+00,    0.4179591836734694e+00/
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of f*w over (a,b),
c                    i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qk15w
c      epmach = r1mach(4)
c      uflow = r1mach(1)
      epmach = 1.0e-30
      uflow = 1.0e-30
c
      centr = 0.5e+00*(a+b)
      hlgth = 0.5e+00*(b-a)
      dhlgth = abs(hlgth)
c
c           compute the 15-point kronrod approximation to the
c           integral, and estimate the error.
c
      fc = df_dv(centr,f_cb)*qwgtc(centr,p1,p2,p3,p4,kp)
      resg = wg(4)*fc
      resk = wgk(8)*fc
      resabs = abs(resk)
!dir$ concurrent
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = df_dv(absc1,f_cb)*qwgtc(absc1,p1,p2,p3,p4,kp)
        fval2 = df_dv(absc2,f_cb)*qwgtc(absc2,p1,p2,p3,p4,kp)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
!dir$ concurrent
      do 15 j=1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        absc1 = centr-absc
        absc2 = centr+absc
        fval1 = df_dv(absc1,f_cb)*qwgtc(absc1,p1,p2,p3,p4,kp)
        fval2 = df_dv(absc2,f_cb)*qwgtc(absc2,p1,p2,p3,p4,kp)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   15 continue
      reskh = resk*0.5e+00
      resasc = wgk(8)*abs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.0.0e+00.and.abserr.ne.0.0e+00)
     *  abserr = resasc*amin1(0.1e+01,
     *  (0.2e+03*abserr/resasc)**1.5e+00)
      if(resabs.gt.uflow/(0.5e+02*epmach)) abserr = amax1((epmach*
     *  0.5e+02)*resabs,abserr)
      return
      end


*
*-----------------------------------------------------------------------
*

      subroutine qpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
c***begin prologue  qpsrt
c***refer to  qage,qagie,qagpe,qagse,qawce,qawse,qawoe
c***routines called  (none)
c***keywords  sequential sorting
c***description
c
c 1.        qpsrt
c           ordering routine
c              standard fortran subroutine
c              real version
c
c 2.        purpose
c              this routine maintains the descending ordering
c              in the list of the local error estimates resulting from
c              the interval subdivision process. at each call two error
c              estimates are inserted using the sequential search
c              method, top-down for the largest error estimate
c              and bottom-up for the smallest error estimate.
c
c 3.        calling sequence
c              call qpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
c
c           parameters (meaning at output)
c              limit  - integer
c                       maximum number of error estimates the list
c                       can contain
c
c              last   - integer
c                       number of error estimates currently
c                       in the list
c
c              maxerr - integer
c                       maxerr points to the nrmax-th largest error
c                       estimate currently in the list
c
c              ermax  - real
c                       nrmax-th largest error estimate
c                       ermax = elist(maxerr)
c
c              elist  - real
c                       vector of dimension last containing
c                       the error estimates
c
c              iord   - integer
c                       vector of dimension last, the first k
c                       elements of which contain pointers
c                       to the error estimates, such that
c                       elist(iord(1)),... , elist(iord(k))
c                       form a decreasing sequence, with
c                       k = last if last.le.(limit/2+2), and
c                       k = limit+1-last otherwise
c
c              nrmax  - integer
c                       maxerr = iord(nrmax)
c
c 4.        no subroutines or functions needed
c***end prologue  qpsrt
c
      real elist,ermax,errmax,errmin
      integer i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,
     *  nrmax
      dimension elist(last),iord(last)
c
c           check whether the list contains more than
c           two error estimates.
c
c***first executable statement  qpsrt
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
c
c           this part of the routine is only executed
c           if, due to a difficult integrand, subdivision
c           increased the error estimate. in the normal case
c           the insert procedure should start after the
c           nrmax-th largest error estimate.
c
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
c ***jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20    continue
c
c           compute the number of elements in the list to
c           be maintained in descending order. this number
c           depends on the number of subdivisions still
c           allowed.
c
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
c
c           insert errmax by traversing the list top-down,
c           starting comparison from the element elist(iord(nrmax+1)).
c
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
c ***jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
c
c           insert errmin by traversing the list bottom-up.
c
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
c ***jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
c
c           set maxerr and ermax.
c
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end


*
*-----------------------------------------------------------------------
*

      real function qwgtc(x,c,p2,p3,p4,kp)
c***begin prologue  qwgtc
c***refer to qk15w
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  weight function, cauchy principal value
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this function subprogram is used together with the
c            routine qawc and defines the weight function.
c***end prologue  qwgtc
c
      real c,p2,p3,p4,x
      integer kp
c***first executable statement
      qwgtc = 0.1e+01/(x-c)
      return
      end

*
*-----------------------------------------------------------------------
*



