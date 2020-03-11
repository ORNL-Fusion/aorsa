SUBROUTINE cubic_B_spline_basis(t, n_sp, sp, dsp_dt)

! Evaluates n_sp cubic B-splines on the interval 0 ² t ² 1.
! Returns zero if t is outside this interval.

! Generates its own knot sequence with 4 knots at t = 0 and 4 knots at t = 1
! With N uniformly spaced interior knots, i.e. at 1/(N+1), 2/(N+1) .. , N/(N+1) you get:
! one spline non-zero at t = 0 with support on [0,1/(N+1)),
! one spline non-zero at t = 1 with support on (N/(N+1), 1]
! one spline that vanishes at t = 0 and t = 1/(N+1) and peaks inside this interval
! one spline that vanishes at t = N/(N+1) and t = 1 and peaks inside this interval
! N splines that peak at n/(N+1) and support on ((n-1)/(N+1), (n+1)/(N+1)).
! Which makes the total spline count n_sp = N + 4. So n_sp needs to be at least 5.
! Total knot count is N + 8 = n_sp + 4

! Note the present implementation is far from optimum.  It is better described as
! quick and dirty.  In particular the knot sequence is regenerated at each call.

IMPLICIT none


	integer, parameter :: rprec = kind(1.0)
	integer, parameter :: iprec = kind(1)
	integer, parameter :: dp = kind(1.0d0)

	REAL(rprec), INTENT(IN) :: t
	INTEGER(iprec),INTENT(IN) :: n_sp
	REAL(rprec), INTENT(OUT) :: sp(n_sp)
!	REAL(rprec), INTENT(OUT),  OPTIONAL :: dsp_dt(n_sp)
	REAL(rprec), INTENT(OUT) :: dsp_dt(n_sp)

	
	  integer :: i, n_interior, ifail
      real(rprec), dimension(n_sp+4) :: xv
      real(rprec) :: dx

	sp = 0
	IF (t < 0.0 .or. t > 1.0) RETURN
	
! Set up knot sequence

! interior values for knots
		n_interior = n_sp - 4
		dx = 1.0/(n_interior + 1)
		
		do i = 1, n_interior
			xv(i+4) = dx*i
		end do
		
! End values of knots
		  do i=1,4
		  xv(i)   =  0.0
		  xv( i+ n_interior + 4) = 1.0
		  end do
		
! Generate spline values
		
		call speval_v(n_sp, xv, t, ifail, sp, dsp_dt)
!      write(6,*) "dsp_dt = "
!      write(6, *) dsp_dt			
	
	RETURN
	
END SUBROUTINE cubic_B_spline_basis



! ******************************************************************************************



subroutine speval_v(n, xk, t, ifail, f, df_dt)

      implicit none

      integer, parameter :: rprec = kind(1.0)
      integer, parameter :: iprec = kind(1)
      integer, parameter :: dp = kind(1.0d0)

      integer :: i, l, n, ip1, ifail
      real(rprec):: xk(n+4), f(n)

!      real(rprec), optional :: df_dt(n)
      real(rprec) :: df_dt(n)

      real(rprec) :: t
!      real(rprec) :: t, a10, a11, a12, a20, a21	<-Original code
      real(rprec) :: b02, b12, b03, b13, b23, b04, b14, b24, b34
      real(rprec) :: d10, d11, d12, d20, d21, d30
      real(rprec) :: tm0, tm1, tm2, tp1, tp2, tp3

!     evaluate the cubic spline with knots xk and b-spline
!     coefficients c at the point t. return function value
!     in f(1) and derivatives in f(2)-f(4).

      ifail=1

	f(1:n) = 0.
	
	  if((t.lt.xk(4)) .or. (t.gt.xk(n+1))) return
	
      i=4
      ip1=n+1
   20 l=(i+ip1)/2
      if(ip1-i.le.1) go to 40
      if(t.lt.xk(l)) go to 30
      i=l
      go to 20
   30 ip1=l
      go to 20
   40 tm2=t-xk(i-2)
      tm1=t-xk(i-1)
      tm0=t-xk(i)
      tp1=xk(i+1)-t
      tp2=xk(i+2)-t
      tp3=xk(i+3)-t
      d10=tp1+tm0
      d11=tp1+tm1
      d20=tp2+tm0
      d12=tp1+tm2
      d21=tp2+tm1
      d30=tp3+tm0
      b12=tp1/d10
      b02=tm0/d10
      b23=tp1*b12/d11
      b13=tm1*b12/d11+tp2*b02/d20
      b03=tm0*b02/d20
      b34=tp1*b23/d12
      b24=tm2*b23/d12+tp2*b13/d21
      b14=tm1*b13/d21+tp3*b03/d30
      b04=tm0*b03/d30

!      f(i)=c(i-3)*b34+c(i-2)*b24+c(i-1)*b14+c(i)*b04	Original code

      f(i) = b04
      f(i-1) = b14
      f(i-2) = b24
      f(i-3) = b34

 ! first derivatives

!      a12=(c(i-2)-c(i-3))/d12	Original code
!      a11=(c(i-1)-c(i-2))/d21
!      a10=(c(i)-c(i-1))/d30
!      f(2)=3.*(a12*b23+a11*b13+a10*b03)


      ifail=0	! normal return

!	IF (.not.PRESENT(df_dt)) RETURN
		
	df_dt(1:n) = 0.

	df_dt(i) = 3.0 * b03/d30
	df_dt(i-1) = 3.0 * (b13/d21 - b03/d30)
	df_dt(i-2) = 3.0 * (b23/d12 - b13/d21)
	df_dt(i-3) = -3.0 * b23/d12

! Don't need higher derivatives but I could get them from f(3), f(4) below
!      a21=(a11-a12)/d11
!      a20=(a10-a11)/d20
!      f(3)=6.*(a21*b12+a20*b02)
!      f(4)=6.*(a20-a21)/d10


      return

end subroutine speval_v
