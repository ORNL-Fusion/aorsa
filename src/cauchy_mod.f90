
        module toeplitz_mod
        implicit none
!       ------------------------------------
!       module to perform matrix multiply
!       with circulant and toeplitz matrices
!       ------------------------------------
        interface
!       -------------------------------------------
!       use FFTPACK routines available from
!       netlib http://netlib.org/fftpack/index.html
!       -------------------------------------------

        subroutine zffti(n,wsave)
        integer n
        integer, parameter :: r15 = selected_real_kind(15,307)
        complex(r15):: wsave(*)
        end subroutine zffti

        subroutine zfftb(n,c,wsave)
        integer n
        integer, parameter :: r15 = selected_real_kind(15,307)
        complex(r15):: c(*), wsave(*)
        end subroutine zfftb

        subroutine zfftf(n,c,wsave)
        integer n
        integer, parameter :: r15 = selected_real_kind(15,307)
        complex(r15):: c(*), wsave(*)
        end subroutine zfftf

        subroutine dffti(n,wsave)
        integer n
        double precision wsave(*)
        end subroutine dffti

        subroutine dfftf(n,r,wsave)
        integer n
        double precision r(*), wsave(*)
        end subroutine dfftf

        subroutine dfftb(n,r,wsave)
        integer n
        double precision r(*), wsave(*)
        end subroutine dfftb

        subroutine dzffti(n,wsave)
        integer n
        double precision wsave(*)
        end subroutine dzffti

        subroutine dzfftf(n,r,azero,a,b,wsave)
        integer n
        double precision r(*),a(*),b(*),wsave(*)
        double precision azero
        end subroutine dzfftf

        subroutine dzfftb(n,r,azero,a,b,wsave)
        integer n
        double precision r(*),a(*),b(*),wsave(*)
        double precision azero
        end subroutine dzfftb



        end interface

        private
        integer, parameter :: r15 = selected_real_kind(15)         
        public :: circulant, circulant_mult
        public :: toeplitz, toeplitz_mult

        contains

        subroutine assert( lcond, msg, ival)
        implicit none
        logical, intent(in) ::  lcond
        character*(*), intent(in) ::  msg
        integer, intent(in) ::  ival

        if (.not.lcond) then
                write(*,*) msg,ival
                stop '** assertion error ** '
        endif
        return
        end subroutine assert
        recursive function padfft(n) result(m)
        integer, intent(in) :: n
        integer :: m

        integer, parameter :: ntry = 4
        integer, dimension(ntry) :: htry 

        integer, parameter :: idebug = 0
        integer :: itry, i


        if (n.le.1) then
                m = 1
                return
        endif

        if (idebug.ge.1) then
                write(*,*) 'padfft(',n,') '
        endif

         htry(1) = 4
         htry(2) = 2
         htry(3) = 3
         htry(4) = 5

        do i=1,ntry
          itry = htry(i)
          if (mod(n,itry).eq.0) then
                  m = itry * padfft( int(n/itry) )
                  return
          endif
        enddo

!       --------------------
!       no match, odd number
!       --------------------
        m = 2* padfft( int( (n+1)/2 ) )
        return
        end function padfft
        subroutine circulant( n, c, Cmat, ldc )
        implicit none
!       ----------------------
!       form circulant matrix
!       given the first column
!       ----------------------
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: c
        integer, intent(in) :: ldc
        double precision, dimension(ldc,*) :: Cmat

        integer :: j

        Cmat(1:n,1) = c(1:n)
        do j=2,n
          Cmat(2:n,j) = Cmat(1:(n-1),j-1)
          Cmat(1,j) = Cmat(n,j-1)
        enddo
        
        return
        end subroutine circulant

        subroutine circulant_mult(n, c, x, y )
        implicit none
!       --------------------------------------
!       compute y = C*x, C is circulant matrix
!       fast algorithm based on FFT
!
!       The columns of Fourier matrix are the
!       eigenvectors of the circulant matrix
!       --------------------------------------
        logical, parameter :: use_zfft = .false.

        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: c
        double precision, dimension(n), intent(in) :: x
        double precision, dimension(n), intent(inout) :: y

        integer, parameter :: idebug = 0
        complex(r15), dimension(n) :: zx, zc, zy
        complex(r15), dimension(4*n+15) :: wsave
        double precision :: t1,t2,  time1, time2, diff

        double precision, dimension(n) :: dx,dc,dy
        double precision, dimension(4*n+15) :: dwsave

        integer :: i,ldc,ierr
        double precision, allocatable, dimension(:) :: y2
        double precision, allocatable, dimension(:,:) :: Cmat
        double precision :: rn

!       ------------------------------
!       use complex fft for simplicity
!       y = ifft( fft(x) .* fft(c) )
!       ------------------------------

                

        call cpu_time(t1)

        if (use_zfft) then
          wsave(1) = 0
          call zffti(n, wsave)
        else
          dwsave(1) = 0
          call dffti(n, dwsave)
        endif
        

        if (use_zfft) then
!        ----------------
!        zc(1:n) = c(1:n)
!        ----------------
         do i=1,n
          zc(i) = c(i)
         enddo
         call zfftf(n,zc,wsave)

        else
!       ---------------
!       dc(1:n) = c(1:n)
!       ---------------
          do i=1,n
            dc(i) = c(i)
          enddo
          call dfftf(n,dc,dwsave)
        endif


        if (use_zfft) then
!        ----------------
!        zx(1:n) = x(1:n)
!        ----------------
         do i=1,n
          zx(i) = x(i)
         enddo
         call zfftf(n,zx,wsave)
        else
         do i=1,n
          dx(i) = x(i)
         enddo
         call dfftf(n,dx,dwsave)
        endif


        if (use_zfft) then
!        ---------------------------
!        zy(1:n) = zx(1:n) * zc(1:n)
!        ---------------------------
         do i=1,n
          zy(i) = zx(i) * zc(i)
         enddo
        else
!        ---------------------------
!        simulate complex arithmetic
!        but do only 1/2 the work
!        ---------------------------
         dy(1) = dx(1)*dc(1)
         do i=2,n-1,2
           dy(i) = dx(i)*dc(i) - dx(i+1)*dc(i+1)
           dy(i+1) = dx(i)*dc(i+1) + dx(i+1)*dc(i)
         enddo
         if (mod(n,2).eq.0) then
           dy(n) = dx(n)*dc(n)
         endif
        endif

        if (use_zfft) then
         call zfftb(n, zy, wsave )
        else
         call dfftb(n, dy, dwsave)
        endif

!        ------------------------------
!        y(1:n) = dble(zy(1:n))/dble(n)
!        ------------------------------

         rn = dble(1)/dble(n)
         if (use_zfft) then
          do i=1,n
           y(i) = dble(zy(i))*rn
          enddo
         else
          do i=1,n
           y(i) = dy(i) * rn
          enddo
         endif
       


        call cpu_time(t2)
        time1 = t2-t1

        if (idebug.ge.1) then
         if (use_zfft) then
           write(*,*) 'i,zc,zx,zy '
           do i=1,n
             write(*,*) i,zc(i),zx(i),zy(i)
           enddo
         else
           write(*,*) 'i,dc,dx,dy '
           do i=1,n
             write(*,*) i,dc(i),dx(i),dy(i)
           enddo
         endif
        endif
        

        if (idebug .ge.1) then
           ldc = n + mod(n+1,2)
           allocate( Cmat(ldc,n), y2(n), stat=ierr)
           call assert(ierr.eq.0,'circulant_mult: alloc',ierr)

           Cmat(1:n,1:n) = 0.0d0
           y2(1:n) = 0.0d0

           call cpu_time(t1)
           call circulant(n,c,Cmat,ldc)
           y2(1:n) = matmul( Cmat(1:n,1:n), x(1:n) )
           call cpu_time(t2)
           time2 = t2-t1

           diff = maxval(abs(y(1:n)-y2(1:n)) )
           write(*,9020) time1,time2,diff
 9020      format(' circulant_mult:fast time1,time2,diff',3(1x,1pe12.4))

           if (idebug.ge.2) then
             do i=1,n
               write(*,*) 'i,y(i),y2(i) ',i,y(i),y2(i)
             enddo
           endif

           deallocate( Cmat, y2, stat=ierr)
           call assert(ierr.eq.0,'circulant_mult:dealloc',ierr)
        endif


        return
        end subroutine circulant_mult

          
          

        subroutine toeplitz(nc, c,nr, r, T,ldt )
        implicit none
!       ---------------------------------------------
!       create Toeplitz matrix, c(1:nc) is 1st column
!       r(1:nr) is 1st row
!       ---------------------------------------------
        integer, intent(in) :: nc, nr, ldt
        double precision, intent(in) :: c(nc), r(nr)
        double precision, intent(inout) :: T(ldt,*)

        integer :: idiag, j
        double precision :: dval
!       ------------------------------
!       Toeplitz matrix T is nc by nr
!       ------------------------------

!       --------------------------
!       construct lower triangular
!       --------------------------
        do idiag=0,(nc-1)
          dval = c(idiag+1)
          do j=1,(nc-idiag)
           T(j+idiag,j) = dval
          enddo
        enddo

!       --------------------------
!       construct upper triangular
!       --------------------------
        do idiag=1,(nr-1)
          dval = r(idiag+1)
          do j=1,(nr-idiag)
            T(j,j+idiag) = dval
          enddo
        enddo

        return
        end subroutine toeplitz

        subroutine toeplitz_mult(nc_in,c,nr,r, x, y )
        implicit none
!        --------------------------------
!        compute y = T*x, T is toeplitz
!        use fast algorithm for Circulant matrix
!        by embeding toeplitz matrix in
!        a larger circulant matrix
!        --------------------------------
        integer, intent(in) :: nc_in, nr
        double precision, intent(in) :: c(nc_in), r(nr)
        double precision, intent(in) :: x(nr)
        double precision, intent(inout) :: y(nc_in)

!       -------------------------------------
!       use fast multiply by circulant matrix
!       -------------------------------------
        integer, parameter :: idebug = 0
        double precision, dimension(:,:), allocatable :: Tmat
        double precision, dimension(:), allocatable :: y2

        integer :: nc, ncc
        double precision, dimension(2*(nc_in+nr)+1024) :: cc, xx, yy
        integer :: i, ip, ldt, ierr
        double precision :: t1, t2, time1, time2, diff

        logical :: isok
        logical, parameter :: use_padfft = .true.


        nc = nc_in
        if (use_padfft) then
          ncc = padfft( nc + (nr-1) )
          nc = ncc - (nr-1)
        else
          ncc = nc + (nr-1)
        endif
        if (idebug.ge.1) then
          write(*,*) 'toepltiz_mult: nc_in,nr,ncc ',nc_in,nr,ncc
        endif

        isok = (ncc .le. size(cc) )
        call assert(isok,'toeplitz_mult: ncc > size(cc) ',ncc)

        call cpu_time(t1)

        do i=1,nc_in
          cc(i) = c(i)
        enddo
        do i=nc_in+1,nc
          cc(i) = 0.0d0
        enddo

        ip = nc + 1
        do i=nr,2,-1
          cc( ip ) = r(i)
          ip = ip + 1
        enddo

        do i=1,nr
          xx(i) = x(i)
        enddo
        do i=nr+1,ncc
          xx(i) = 0.0d0
        enddo
        do i=1,ncc
          yy(i) = 0.0d0
        enddo

        call circulant_mult( ncc, cc, xx, yy )
        do i=1,nc_in
          y(i) = yy(i)
        enddo

        call cpu_time(t2)
        time1 = t2-t1
        
        if (idebug.ge.1) then

          nc = nc_in
          call cpu_time(t1)
          ldt = nc + mod(nc+1,2)
          allocate( Tmat(ldt,nr), y2(nc), stat=ierr)
          call assert(ierr.eq.0,'toeplitz_mult: alloc ',ierr)

          call toeplitz( nc,c,nr,r, Tmat,ldt )
          y2(1:nc) = matmul( Tmat(1:nc, 1:nr), x(1:nr) )
          call cpu_time(t2)
          time2 = t2 - t1

          diff = maxval( abs(y(1:nc) - y2(1:nc)) )
          write(*,9020) time1, time2, diff
 9020     format(' toeplitz_mult: time1, time2, diff ',3(1x,1pe14.4))

          deallocate( Tmat,y2, stat=ierr)
          call assert(ierr.eq.0,'toeplitz_mult:dealloc',ierr)
        endif

        
        return
        end subroutine toeplitz_mult

        end module toeplitz_mod


        module cauchy_mod
        use toeplitz_mod
        implicit none
!       ---------------------------------------------
!       module to compute principal value integral of
!       g(z) = integrate( f(x)/(x-z), x=a..b )
!       where f(x) is given as discrete values
!       (x(i),y(i)) i=1..n
!
!
!       hilb0 : computes g(z) assuming piecewise constant value
!               on each interval
!       hilb1 : computes g(z) assuming piecewise linear value
!               on each interval
!       hilb3 : computes g(z) assuming cubic spline coefficients
!               are precomputed
!
!       setup_cauchy:  preevaluates g(z) at mid points, then
!                      uses a spline approximation
!       eval_cauchy:   performs interpolation if z is within
!                      the interior of grid
!
!       ----------------------------------------------
!       Auxiliary routines:
!       spline/eval:   cubic spline routines obtained from 
!                      FMM library on netlib.org
!
!       test_cauchy_mod: simple tester
!       test_cauchy_driver: simple driver to test code
!       ---------------------------------------------




!       ------------------------------------------------------
!       data structure to save precomputed spline coefficients
!       ------------------------------------------------------
        type cauchy_t 
          integer :: n
          double precision, pointer, dimension(:) :: x,y,b,c,d
          double precision, pointer, dimension(:) :: z,gz, zb,zc,zd
        end type cauchy_t

        private :: fcn1,fcn2,fcn3, spline, seval, assert



        interface eval_cauchy
          module procedure eval_cauchy_1, eval_cauchy_vec
        end interface

        contains
      subroutine spline (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
!
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
      integer nm1, ib, i
      double precision t
!
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.0d0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0d0
      c(n) = 0.0d0
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
!
!  back substitution
!
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.0d0*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0d0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.0d0*c(i)
   40 continue
      c(n) = 3.0d0*c(n)
      d(n) = d(n-1)
      return
!
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.0d0
      d(1) = 0.0d0
      b(2) = b(1)
      c(2) = 0.0d0
      d(2) = 0.0d0
      return
      end subroutine spline
      double precision function seval(n, u, x, y, b, c, d)
      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n)
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
      integer i, j, k
      double precision dx

!     ------------------------------------
!     initial guess assumes a uniform mesh
!     ------------------------------------
      i = int(  ((u-x(1))/(x(n)-x(1)))*dble(n) )
      i = max(1, min(n-1, i) )
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
!
!  evaluate spline
!
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end function seval

        subroutine assert( lcond, msg, ival)
        implicit none
        logical, intent(in) ::  lcond
        character*(*), intent(in) ::  msg
        integer, intent(in) ::  ival

        if (.not.lcond) then
                write(*,*) msg,ival
                stop '** assertion error ** '
        endif
        return
        end subroutine assert

        subroutine vlog(logw, w, n)
        implicit none
!       ------------------------------------------------
!       vector log evaluation
!       similar to IBM MASSV
!       use the real MASSV vlog if running on IBM system
!       ------------------------------------------------
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: w
        double precision, dimension(n), intent(inout) :: logw

        integer :: i
        intrinsic :: log

        do i=1,n
          logw(i) = log( w(i) )
        enddo

        return
        end subroutine vlog


        subroutine hilb3(n,x,y,b,c,d,  z,gzi)
        implicit none

!       ------------------------------------------------
!       evaluate  g(z) = integral( f(x)/(x-z), x=x1..xn )
!       assume piecewise cubic polynomial
!       over interval  [x(i),x(i+1)]
!       ------------------------------------------------
        integer, intent(in) ::  n
        double precision, intent(in) ::  x(n),y(n),b(n),c(n),d(n)
        double precision, intent(in) ::  z 
        double precision, intent(inout) :: gzi

        integer idebug
        parameter(idebug = 0)


        integer i
        double precision gz0, gz1, gz2, gz3
        double precision dz, dz0, dz1
        double precision a0,a1,a2,a3, ans0
        intrinsic :: log, abs

        logical, parameter :: use_vlog = .true.
        double precision, dimension(n-1) :: w, logw
        double precision :: t1, t2

!       ---------------------------------------
!       over [x(i),x(i+1)], the function is
!       approximated by the cubic
!
!       f(x) = y(i) + (x-x(i))*b(i) + (x-x(i))^2*c(i) + 
!              (x-x(i))^3*d(i)
!
!       f(x) = y(i) + ((x-z) + (z-x(i)))*b(i) + 
!              ((x-z)+(z-x(i)))^2*c(i) + 
!              ((x-z)+(z-x(i)))^3*d(i)
!            = a0 + a1*(x-z) + a2*(x-z)^2 + a3*(x-z)^3
!
!       ---------------------------------------

        if (use_vlog) then
          do i=1,n-1
            dz1 = x(i+1)-z
            dz0 = x(i)-z
            w(i) = abs( dz1 / dz0 )
          enddo

          call cpu_time(t1)
          call vlog( logw, w, n-1 )
          call cpu_time(t2)

          if (idebug.ge.1) then
            write(*,9010) n,t2-t1
 9010       format(' hilb3,vlog: n,time ',i5,1x,1pe12.4)
          endif

        endif


        gz0 = 0.0d0
        gz1 = 0.0d0
        gz2 = 0.0d0
        gz3 = 0.0d0
        do i=1,n-1
          dz = z-x(i)

          a0 = y(i) + dz*( b(i) + dz*(c(i) + dz*d(i)))
          a1 = b(i) + dz*(2.0d0*c(i) + 3.0d0*dz*d(i))
          a2 = c(i) + 3.0d0*dz*d(i)
          a3 = d(i)

!         --------------------------------------
!         integral( a0/(x-z), x(i)..x(i+1) )  = 
!         a0*ln( abs(  (x(i+1)-z)/(x(i)-z) ) )
!
!         other integrals are simple since the (x-z) term
!         is cancelled out
!         --------------------------------------
          dz1 = x(i+1)-z
          dz0 = x(i)-z

          if (use_vlog) then
            gz0 = gz0 + a0*logw(i)
          else
            gz0 = gz0 + a0*log( abs(  dz1 / dz0 ) )
          endif

          gz1 = gz1 + a1*(x(i+1)-x(i))
          gz2 = gz2 + (a2/2.0d0)*(x(i+1)-x(i))*(x(i+1)+x(i)-2.0d0*z)
          gz3 = gz3 + (a3/3.0d0)*(x(i+1)-x(i))*         &                      
     &                (dz1*dz1 + dz1*dz0 + dz0*dz0)

        enddo

        gzi = gz0 + gz1 + gz2 + gz3
        if (idebug.ge.2) then
           write(*,9020) z, gzi
 9020      format(' hilb3: z,gzi ', 2(1x,1pe12.4))

           write(*,9030) gz0, gz1, gz2, gz3
 9030      format(' hilb3: gz0,gz1,gz2,gz3 ',4(1x,1pe12.4))

!        -------------------------------
!        compute the low accuracy answer
!        -------------------------------
           ans0 = 0.0d0
           do i=1,n-1
             a0 = (y(i) + y(i+1))/2.0d0
             dz1 = x(i+1)-z
             dz0 = x(i)-z
             if (dz0.eq.0.0d0) then
                write(*,9040) n,i,z,x(i),x(i+1)
 9040           format(' hilb3:n,i,z,x(i),x(i+1) ',     &                
     &                   2(1x,i5),3(1x,1pe12.4))
             endif

             ans0 = ans0 + a0*log( abs(  dz1 / dz0 ) )
           enddo

           write(*,9050) ans0
 9050      format(' hilb3: ans0 = ',1pe12.4)
        endif

        return
        end  subroutine hilb3

        subroutine hilb1(n,x,y,  zi, gzi )
        implicit none
!       ------------------------------------------------
!       evaluate  g(z) = integral( f(x)/(x-z), x=x1..xn )
!       assume piecewise linear value f(x) = a1*(x-z) + a0
!       over interval  [x(i),x(i+1)]
!       ------------------------------------------------
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: x,y
        double precision, intent(in) :: zi
        double precision, intent(inout) :: gzi

        integer :: i
        double precision :: a0, a1, gz0, gz1
        intrinsic :: log, abs

        integer, parameter :: idebug = 0
        logical, parameter :: use_vlog = .true.
        double precision, dimension(n-1) :: w, logw
        double precision :: t1, t2

        if (use_vlog) then
          do i=1,n-1
            w(i) = abs( (x(i+1)-zi)/(x(i)-zi) )
          enddo

          call cpu_time(t1)
          call vlog(logw, w, n-1 )
          call cpu_time(t2)

          if (idebug.ge.1) then
          write(*,*) 'hilb1,vlog: n, time ',n,t2-t1
          endif

        endif

        gz0 = 0.0d0
        gz1 = 0.0d0
        do i=1,n-1
!         -----------------------------------------
!         assume y = a1 * (x-z) + a0
!         a1 = slope of line
!         
!         
!         integral( a0/(x-z) + a1, x=x(i)..x(i+1))
!         = a0*log( abs((x(i+1)-z)/(x(i)-z)) ) + a1*(x(i+1)-x(i))
!         -----------------------------------------
          a1 = (y(i+1)-y(i))/(x(i+1)-x(i))
          a0 = (y(i)*(x(i+1)-zi) - y(i+1)*(x(i)-zi)  )/(x(i+1)-x(i))
          gz1 = gz1 + a1*(x(i+1)-x(i))
          if (use_vlog) then
            gz0 = gz0 + a0*logw(i)
          else
            gz0 = gz0 + a0*log(abs( (x(i+1)-zi)/(x(i)-zi) ))
          endif
        enddo

        gzi = gz0 + gz1

        return
        end subroutine hilb1
        subroutine hilb0(n,x,y,  zi, gzi )
        implicit none
!       ------------------------------------------------
!       evaluate  g(z) = integral( f(x)/(x-z), x=x1..xn )
!       assume piecewise constant value (f(x(i)) + f(x(i+1)))/2
!       over interval  [x(i),x(i+1)]
!       ------------------------------------------------        
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: x,y
        double precision, intent(in) :: zi
        double precision, intent(inout) :: gzi

        integer :: i
        double precision :: a0
        intrinsic :: log, abs

        integer, parameter :: idebug = 0
        logical, parameter :: use_vlog = .true.
        double precision, dimension(n-1) :: w, logw
        double precision :: t1, t2

        if (use_vlog) then
          do i=1,n-1
            w(i) = abs( (x(i+1)-zi)/(x(i)-zi) )
          enddo

          call cpu_time(t1)
          call vlog(logw, w, n-1 )
          call cpu_time(t2)

          if (idebug.ge.1) then
          write(*,*) 'hilb0,vlog: n, time ',n,t2-t1
          endif

        endif

        gzi = 0.0d0
        do i=1,n-1
          a0 = (y(i) + y(i+1))/2.0d0
          if (use_vlog) then
            gzi = gzi + a0*logw(i)
          else
            gzi = gzi + a0*log(abs( (x(i+1)-zi)/(x(i)-zi) ))
          endif
        enddo

        return
        end subroutine hilb0



!23456789!23456789!23456789!23456789!23456789!23456789!23456789!23456789
        subroutine setup_cauchy(n,xin,yin,  p)
        implicit none
!       -----------------------------------------------
!       precompute g(z) = integral(f(x)/(x-z),x=x1..xn)
!       and use faster interpolation
!       -----------------------------------------------
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) ::  xin,yin
        type(cauchy_t), intent(inout) :: p

        integer, parameter :: idebug = 0
        integer :: i, m,k
        double precision :: h, zi, gzi, t2, t1
        double precision, pointer, dimension(:) :: x,y,b,c,d
        double precision, pointer, dimension(:) :: z,gz,zb,zc,zd

        logical, parameter :: use_hilb1 = .false.
!       ---------------------------------
!       generate spline approximation
!       
!       over [x(i),x(i+1)], the function is
!       approximated by the cubic
!
!       f(x) = y(i) + (x-x(i))*b(i) + (x-x(i))^2*c(i) + 
!              (x-x(i))^3*d(i)
!
!       ---------------------------------
        integer ierr 
        logical isok

        nullify(x)
        nullify(y)
        nullify(b)
        nullify(c)
        nullify(d)

        nullify(z)
        nullify(gz)
        nullify(zb)
        nullify(zc)
        nullify(zd)

        allocate( x(n),y(n),  b(n),c(n),d(n), stat=ierr)
        isok = (ierr.eq.0)
        call assert( isok,       &                                       
     &      'setup_cauchy: allocate x(n) failed, ierr=',ierr)

        x(:) = 0.0d0
        y(:) = 0.0d0
        b(:) = 0.0d0
        c(:) = 0.0d0
        d(:) = 0.0d0


        do i=1,n
          x(i) = xin(i)
          y(i) = yin(i)
        enddo

!       -----------------------------------------
!       z's located at mid point of [x(i),x(i+1)]
!       k extra points on both ends
!       -----------------------------------------
        k = 8
        m = (n-1) + 2*k
        allocate( z(m),gz(m), zb(m),zc(m),zd(m), stat=ierr )

        isok = (ierr.eq.0)
        call assert( isok,     &                                         
     &       'setup_cauchy: allocate z(m) failed, ierr=',ierr)

        z(:) = 0.0d0
        gz(:) = 0.0d0
        zb(:) = 0.0d0
        zc(:) = 0.0d0
        zd(:) = 0.0d0

!       ---------------
!       setup interior points
!       ---------------
        do i=1,n-1
          z(k+i) = (x(i) + x(i+1))/2.0d0;
        enddo
        h = (x(2)-x(1))
        h = (h/2.0d0)/dble(k+1)
        do i=k,1,-1
          z(i) = z(i+1)-h
        enddo
        h = x(n)-x(n-1)
        h = h/2.0d0/dble(k+1)
        do i=1,k
          z(k+(n-1)+i) = z(k+(n-1)+(i-1)) + h
        enddo

        
        call cpu_time(t1)
        call spline( n, x,y,   b,c,d)
        call cpu_time(t2)
        if (idebug .ge. 1) then
          write(*,9010) n,t2-t1
 9010     format(' setup_cauchy,spline1: n, time ',   &                  
     &             i5,1x,1pe12.4)
        endif
!
!       --------------------------------------------------
!       evaluate g(z) = integral( f(x)/(x-z), x=x(1)..x(n)
!       --------------------------------------------------


!       ------------------
!       precompute answer 
!       ------------------
        call cpu_time(t1)
        do i=1,m
           zi = z(i)
           gzi = 0.0d0

           if (use_hilb1) then
             call hilb1(n,x,y, zi,gzi)
           else
             call hilb3( n, x, y, b, c,d, zi, gzi )
           endif

           gz(i) = gzi
        enddo
        call cpu_time(t2)
        if (idebug.ge.1) then
          write(*,9020) n, t2-t1
 9020     format(' setup_cauchy,hilb3: n, time ',     &                  
     &          i5,1x,1pe12.4)
        endif

        call cpu_time(t1)
        call spline(size(z), z, gz, zb, zc, zd )
        call cpu_time(t2)
        if (idebug.ge.1) then
          write(*,9030) n,t2-t1
 9030     format(' setup_cauchy,spline2: n, time ',    &                 
     &           i5,1x,1pe12.4)
        endif

!       --------------------------------
!       save arrays for later evaluation
!       --------------------------------

        p%n = n

        p%x => x
        p%y => y
        p%b => b
        p%c => c
        p%d => d

        p%z => z
        p%gz => gz
        p%zb => zb
        p%zc => zc
        p%zd => zd


        return
        end subroutine setup_cauchy

        subroutine eval_cauchy_1( zi, gzi, p )
        implicit none
        double precision, intent(in) :: zi
        double precision, intent(inout) :: gzi
        type(cauchy_t), intent(in) :: p

        double precision, dimension(1) :: zt, gzt
        integer :: m

        m = 1
        zt(1) = zi
        gzt(1) = 0.0d0

        call eval_cauchy_vec(m,zt,gzt,p)

        gzi = gzt(1)
        return
        end subroutine eval_cauchy_1

        subroutine eval_cauchy_vec( m, zt, gzt, p )
        implicit none
        integer , intent(in) :: m
        double precision, dimension(m), intent(in) :: zt
        double precision, dimension(m), intent(inout) :: gzt
        type( cauchy_t), intent(in)  :: p

        integer, parameter :: idebug = 0

        integer :: j,n, jmax
        double precision, pointer, dimension(:) :: x,y,b,c,d
        double precision, pointer, dimension(:) :: z,gz,zb,zc,zd
        double precision :: zlo, zhi
        double precision :: abserr, maxerr, gztj
        logical :: in_range



!       -------------------------------------------
!       evaluate integral( f(x)/(x-zt(j)), x=x1..xn) 
!       for j = 1:m
!       -------------------------------------------

        n = p%n
        x => p%x
        y => p%y
        b => p%b
        c => p%c
        d => p%d

        z => p%z
        gz => p%gz
        zb => p%zb
        zc => p%zc
        zd => p%zd


        zlo = x(lbound(x,1)+2)
        zhi = x(ubound(x,1)-2)

        if (idebug.ge.1) then
          write(*,*) 'eval_cauchy:lbound(z,1), ubound(z,1) ',   &        
     &                  lbound(z,1),ubound(z,1)
        endif

        maxerr = 0.0d0
        jmax = -1
        do j=1,m
          in_range = (zlo.le.zt(j)).and.(zt(j).le.zhi)
          if (in_range) then
!           ---------------------
!           use spline evaluation
!           ---------------------
            gzt(j) = seval(size(z), zt(j), z, gz, zb, zc, zd)
            gztj = 0.0d0
            if (idebug.ge.1) then
              call hilb3( n, x, y, b,c, d,   zt(j), gztj )
              abserr = abs(gztj - gzt(j))
              if (abserr .gt. maxerr) then
                  maxerr = abserr
                  jmax = j
              endif
            endif

          else
!           --------------------
!           use basic evaluation
!           --------------------
            call hilb3( n, x, y, b,c, d,   zt(j), gzt(j) )
          endif
        enddo


        if (idebug.ge.2) then
          if ( (1.le.jmax).and.(jmax.le.(m-1)) ) then
             write(*,*) 'eval_cauchy: jmax, maxerr ',jmax,maxerr
          endif
        endif

        return
        end subroutine eval_cauchy_vec



        subroutine reset_cauchy( p )
        implicit none
!       --------------
!       release memory
!       --------------
        type(cauchy_t), intent(inout) :: p

        integer :: ierr


        if (associated(p%x)) then
          deallocate( p%x, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%x ',ierr)
        endif

        if (associated(p%y)) then
          deallocate( p%y, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%y ',ierr)
        endif

        if (associated(p%b)) then
          deallocate( p%b, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%b ',ierr)
        endif

        if (associated(p%c)) then
          deallocate( p%c, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%c ',ierr)
        endif

        if (associated(p%d)) then
          deallocate( p%d, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%d ',ierr)
        endif

        if (associated(p%z)) then
          deallocate( p%z, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%z ',ierr)
        endif

        if (associated(p%gz)) then
          deallocate( p%gz, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%gz ',ierr)
        endif

        if (associated(p%zb)) then
          deallocate( p%zb, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%zb ',ierr)
        endif

        if (associated(p%zc)) then
          deallocate( p%zc, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%zc ',ierr)
        endif

        if (associated(p%zd)) then
          deallocate( p%zd, stat=ierr )
          call assert(ierr.eq.0,'reset_cauchy: deallocate p%zd ',ierr)
        endif




        nullify(p%x)
        nullify(p%y)
        nullify(p%z)
        nullify(p%b)
        nullify(p%c)
        nullify(p%d)

        nullify(p%z)
        nullify(p%gz)
        nullify(p%zb)
        nullify(p%zc)
        nullify(p%zd)

        p%n = -1
        return
        end subroutine reset_cauchy

        double precision function fcn1(x)
        double precision, intent(in) :: x
        fcn1 = sin(x)
        return
        end function fcn1


        double precision function fcn2(x)
        double precision, intent(in) :: x
        fcn2 = exp(-x*x)
        return
        end function fcn2

        double precision function fcn3(x)
        double precision, intent(in) :: x
        fcn3 = cos(x)
        return
        end function fcn3

!234567890!234567890!234567890!234567890!234567890!234567890!234567890!234567890

        subroutine test_cauchy_mod(n, xlo, xhi, fcn )
        implicit none
!       -------------------------------------
!       simple code for testing and debugging
!       -------------------------------------
        integer, intent(in) :: n
        double precision, intent(in) :: xlo, xhi
        double precision, external :: fcn


        integer, parameter :: idebug = 0
        double precision, dimension(n) :: x, y

        integer :: m, i, j, ierr
        double precision, dimension(:), allocatable :: zt,gzt
        double precision, dimension(:), allocatable :: gztm,gztm1,gztm3
        double precision :: dx,gzt0,gzt1,gzt3, t, t1,t2

        type (cauchy_t) :: p,p3



        dx = (xhi-xlo)/dble(n-1)
        do i=1,n
          x(i) = xlo + dble(i-1)*dx
          y(i) = fcn( x(i) )
        enddo

        call cpu_time(t1)
        call setup_cauchy0(n, x, y, p )
        call cpu_time(t2)
        write(*,9010) n, t2-t1
 9010   format(' test_cauchy_mod,setup_cauchy0: n,time ',       &        
     &         i5,(1x,1pe12.4))
        call reset_cauchy(p)


        call cpu_time(t1)
        call setup_cauchy1(n, x, y, p )
        call cpu_time(t2)
        write(*,9020) n, t2-t1
 9020   format(' test_cauchy_mod,setup_cauchy1: n,time ',      &         
     &         i5,(1x,1pe12.4))
        call reset_cauchy(p)

        call cpu_time(t1)
        call setup_cauchy3(n, x, y, p3 )
        call cpu_time(t2)
        write(*,9025) n, t2-t1
 9025   format(' test_cauchy_mod,setup_cauchy3: n,time ',    &           
     &         i5,(1x,1pe12.4))


        call cpu_time(t1)
        call setup_cauchy( n, x, y,  p )
        call cpu_time(t2)
        write(*,9030) n, t2-t1
 9030   format(' test_cauchy_mod,setup_cauchy: n, time ',    &           
     &           i5,(1x,1pe12.4))


!       ----------------------------------------------------
!       compare results from setup_cauchy3 and setup_cauchy
!       ----------------------------------------------------
        write(*,*) 'p3 vs p: diff z ', maxval(abs(p%z-p3%z))
        write(*,*) 'p3 vs p: diff gz ', maxval(abs(p%gz-p3%gz))
        if (idebug.ge.1) then
          do i=1,10
           write(*,*)'i,p3%gz,p%gz',i,p3%gz(i),p%gz(i)
          enddo
          do j=1,10
           i = size(p%gz)-j+1
           
           write(*,*)'i,p3%gz,p%gz',i,p3%gz(i),p%gz(i)
          enddo
        endif


        write(*,7010) maxval(abs(p%zb-p3%zb))
 7010   format(' p3 vs p: diff zb ',1pe14.4)

        write(*,7020) maxval(abs(p%zc-p3%zc))
 7020   format(' p3 vs p: diff zc ',1pe14.4)

        write(*,7030) maxval(abs(p%zd-p3%zd))
 7030   format(' p3 vs p: diff zd ',1pe14.4)



        dx = dble(xhi-xlo)/dble(m-1)

        t = 0.1d0
        m = n-1
        allocate(zt(m),gzt(m),gztm(m),gztm1(m),gztm3(m),stat=ierr)
        call assert(ierr.eq.0,                           &               
     &          'test_cauchy_mod: allocate zt',ierr)
        do i=1,m
          zt(i) = t*x(i) + (1.0d0-t)*x(i+1)
        enddo

        call cpu_time(t1)
        call eval_cauchy( m, zt, gzt, p3 )
        call cpu_time(t2)
        write(*,9040) m, t2-t1
 9040   format(' eval_cauchy: m, time ',i5,(1x,1pe12.4))

        call cpu_time(t1)
        do i=1,m
          call hilb0( n,x,y,  zt(i), gzt0 )
          gztm(i) = gzt0
        enddo
        call cpu_time(t2)
        write(*,9050) m, t2-t1
 9050   format(' test_cauchy_mod,hilb0: m,time ',i5,(1x,1pe12.4))


        call cpu_time(t1)
        do i=1,m
          call hilb1(n,x,y,   zt(i), gzt1 )
          gztm1(i) = gzt1
        enddo
        call cpu_time(t2)
        write(*,9060) m, t2-t1
 9060   format(' test_cauchy_mod,hilb1: m,time ',i5,(1x,1pe12.4))

        call cpu_time(t1)
        do i=1,m
          call hilb3(n,x,y,   p%b, p%c, p%d, zt(i), gzt3 )
          gztm3(i) = gzt3
        enddo
        call cpu_time(t2)
        write(*,9070) m, t2-t1
 9070   format(' test_cauchy_mod,hilb3: m,time ',i5,(1x,1pe12.4))



        write(*,*) 'i,   zt(i),   gzt(i),  diff 0, diff 1, diff 3 '
        do i=1,m
          gzt0 = gztm(i)
          gzt1 = gztm1(i)
          gzt3 = gztm3(i)
          write(*,9080) i,zt(i),gzt(i),           &                      
     &          abs(gzt(i)-gzt0),abs(gzt(i)-gzt1),abs(gzt(i)-gzt3)
 9080     format(1x, i5, 5(1x,1pe12.4))
        enddo


        call reset_cauchy( p )

        deallocate( zt, gzt, gztm, gztm1, gztm3, stat=ierr)
        call assert(ierr.eq.0,                     &                     
     &          'test_cauchy_mod: deallocate zt',ierr)
        return
        end subroutine test_cauchy_mod
!234567890!234567890!234567890!234567890!234567890!234567890!234567890!234567890

        subroutine test_cauchy_driver(indev,filename)
        implicit none
!       ---------------------------------------
!       simple driver to read in parameters and
!       run tester code
!       ---------------------------------------
        integer, intent(in) :: indev
        character(len=*), intent(in) :: filename


        integer, parameter :: idebug = 0

        integer :: icase, ierr
        double precision :: xlo, xhi
        integer :: n

! ---------------------------
! input file is something like
! 64  -3.14  3.14  1
! 128  -3.14  3.14  1
! 256  -3.14  3.14  1
! 1024  -3.14  3.14  1
! 2048  -3.14  3.14  1
! ---------------------------
        if (idebug.ge.1) then
          write(*,*) 'test_cauchy_driver: indev ',indev
          write(*,*) 'filename:',trim(filename),':'
        endif
        

        open(indev,file=trim(filename),form='formatted',    &            
     &       access='sequential',iostat=ierr)
        call assert(ierr.eq.0,                           &               
     &          'test_cauchy_driver: open failed ',ierr)
        rewind(indev,iostat=ierr)
        call assert(ierr.eq.0,                         &                 
     &          'test_cauchy_driver: rewind failed ',ierr)
        

        ierr = 0
        do while (ierr.eq.0)

          read(indev,*,iostat=ierr) n,xlo,xhi,icase
          if (ierr.ne.0) exit

          if (idebug.ge.1) then
            write(*,*) 'n,xlo,xhi,icase',n,xlo,xhi,icase
          endif


          call test_cal_log(n,xlo,xhi)

          if (icase.eq.1) then
            call test_cauchy_mod(n,xlo,xhi, fcn1 )
          else if (icase.eq.2) then
            call test_cauchy_mod(n,xlo,xhi, fcn2 )
          else
            call test_cauchy_mod(n,xlo,xhi, fcn3 )
          endif

        enddo


        close(indev,iostat=ierr)
        call assert(ierr.eq.0,             &                             
     &       'test_cauchy_driver: clsoe failed ',ierr)

        return
        end subroutine test_cauchy_driver



        subroutine cal_log(n,xlo,xhi, a0, y)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: xlo, xhi, a0(n-1)
        double precision, intent(inout) :: y(n-1)

        double precision :: h, zi, yi, xi, xip1, xj,xjp1
        integer :: i,j
        
        h = (xhi-xlo)/dble(n-1)

        do i=1,n-1
          xi = xlo + dble(i-1)*h
          xip1 = xi + h

          zi = (xi+xip1)/2.0d0

          yi = 0.0d0
          do j=1,n-1
            xj = xlo + dble(j-1)*h 
            xjp1 = xj + h
            yi = yi + a0(j)*log( abs( (xjp1-zi)/(xj-zi) ) )
          enddo

          y(i) = yi
        enddo

        return
        end subroutine cal_log
        subroutine fast_cal_log(n,xlo,xhi, a0, y)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: xlo,xhi
        double precision, dimension(n-1), intent(in) :: a0
        double precision, dimension(n-1), intent(inout) :: y

        integer :: j,nc, nr
        double precision, dimension(n-1) :: c1,c2,r1,r2, cin,rin
        double precision :: h


        h = (xhi-xlo)/dble(n-1)

        nc = n-1
        nr = n-1

        c1(1) = log(h/2.0d0)
        c1(2) = log(h/2.0d0)
        do j=3,nc
          c1(j) = log( h/2.0d0 + dble(j-2)*h )
        enddo

        r1(1) = c1(1)
        do j=2,nr
          r1(j) = log( h/2.0d0 + (j-1)*h )
        enddo

        c2(1) = log(h/2.0d0)
        c2(2) = log(h + h/2.0d0)
        do j=3,nc
          c2(j) = log( h/2.0d0 + dble(j-1)*h)
        enddo

        r2(1) = c2(1)
        r2(2) = log(h/2.0d0)
        do j=3,nr
          r2(j) = log( dble(j-1)*h - h/2.0d0)
        enddo

        do j=1,nc
          cin(j) = c1(j) - c2(j)
        enddo
        do j=1,nr
          rin(j) = r1(j) - r2(j)
        enddo
        call toeplitz_mult(nc,cin,nr,rin, a0, y)


        return
        end subroutine fast_cal_log
        subroutine test_cal_log(n,xlo,xhi)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: xlo, xhi
!       --------------------------------------------
!       check implementation of fast toeplitz method
!       --------------------------------------------
        double precision, dimension(n-1) :: y,y2, a0
        double precision :: t1, t2, time1, time2, diff

        integer :: i
        integer, parameter :: idebug = 0


        call random_number(a0)

        call cpu_time(t1)
        call fast_cal_log(n,xlo,xhi,a0,y2)
        call cpu_time(t2)
        time1 = t2-t1

        call cpu_time(t1)
        call cal_log(n,xlo,xhi, a0, y )
        call cpu_time(t2)
        time2 = t2-t1



        diff = maxval( abs(y(1:(n-1)) - y2(1:(n-1))) )
        write(*,9020) time1, time2, diff
 9020   format(1x,'test_cal_log: fast time1, time2, diff ',  &           
     &         3(1x,1pe14.4))

        if (idebug .ge. 1) then
          write(*,*) 'test_cal_log: '
          write(*,*) 'i, y(i), y2(i) '
          do i=1,n-1
           write(*,9030) i,y(i),y2(i)
 9030      format(1x,i6,2(1x,1pe14.4))
          enddo
        endif

        return
        end subroutine test_cal_log

        subroutine setup_cauchy0(n,xin,yin, p)
        implicit none
!       -----------------------------------------------
!       precompute g(z) = integral(f(x)/(x-z),x=x1..xn)
!       and use faster interpolation
!       -----------------------------------------------
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: xin,yin
        type(cauchy_t), intent(inout) :: p

        double precision, pointer, dimension(:) :: x,y,b,c,d
        double precision, pointer, dimension(:) :: z,gz,zb,zc,zd


        integer, parameter :: idebug = 0
        logical :: isok
        integer :: i, m, ierr
        double precision, dimension(n-1) :: gz2, a0vec
        double precision :: xlo, xhi, zi, gzi
        double precision :: t1, t2, time1, time2, diff


        nullify(x)
        nullify(y)
        nullify(b)
        nullify(c)
        nullify(d)

        nullify(z)
        nullify(gz)
        nullify(zb)
        nullify(zc)
        nullify(zd)

        allocate( x(n),y(n),  b(n),c(n),d(n), stat=ierr)
        isok = (ierr.eq.0)
        call assert( isok,                                   &           
     &      'setup_cauchy0: allocate x(n) failed, ierr=',ierr)

        x(:) = 0.0d0
        y(:) = 0.0d0
        b(:) = 0.0d0
        c(:) = 0.0d0
        d(:) = 0.0d0


        do i=1,n
          x(i) = xin(i)
          y(i) = yin(i)
        enddo

!       -----------------------------------------
!       z's located at mid point of [x(i),x(i+1)]
!       -----------------------------------------
        m = n-1
        allocate( z(m),gz(m), zb(m),zc(m),zd(m), stat=ierr )

        isok = (ierr.eq.0)
        call assert( isok,                                   &           
     &       'setup_cauchy0: allocate z(m) failed, ierr=',ierr)

        z(:) = 0.0d0
        gz(:) = 0.0d0
        zb(:) = 0.0d0
        zc(:) = 0.0d0
        zd(:) = 0.0d0


        do i=1,m
         z(i) = (x(i)+x(i+1))*0.5d0
        enddo


          call cpu_time(t1)
          do i=1,n-1
             a0vec(i) = (y(i+1)+y(i))*0.5d0
          enddo

          xlo = x(1)
          xhi = x(n)
          call fast_cal_log( n, xlo, xhi, a0vec, gz )
          call cpu_time(t2)
          time1 = t2-t1

       if (idebug .ge. 1) then

        call cpu_time(t1)
        do i=1,m
          zi = z(i)
          call hilb0( n,x,y,  zi,gzi )
          gz2(i) = gzi
        enddo
        call cpu_time(t2)
        time2 = t2-t1


        diff = maxval( abs( gz - gz2 ) )
        write(*,9010) time1, time2, diff
 9010   format(1x,'setup_cauchy0: fast time1, time2, diff ',  &          
     &              3(1x,1pe12.4))

        endif



        call cpu_time(t1)
        call spline(size(z), z, gz, zb, zc, zd )
        call cpu_time(t2)
        if (idebug.ge.1) then
          write(*,9020) n,t2-t1
 9020     format(' setup_cauchy0,spline2: n, time ',i5,1x,1pe12.4)
        endif

!       --------------------------------
!       save arrays for later evaluation
!       --------------------------------

        p%n = n

        p%x => x
        p%y => y
        p%b => b
        p%c => c
        p%d => d

        p%z => z
        p%gz => gz
        p%zb => zb
        p%zc => zc
        p%zd => zd


        return
        end subroutine setup_cauchy0
        subroutine setup_cauchy1(n,xin,yin, p)
        implicit none
!       -----------------------------------------------
!       precompute g(z) = integral(f(x)/(x-z),x=x1..xn)
!       and use faster interpolation
!       -----------------------------------------------
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: xin,yin
        type(cauchy_t), intent(inout) :: p

        double precision, pointer, dimension(:) :: x,y,b,c,d
        double precision, pointer, dimension(:) :: z,gz,zb,zc,zd


        integer, parameter :: idebug = 0
        logical :: isok
        integer :: i, m, ierr
        double precision, dimension(n-1) :: gz2
        double precision :: zi, gzi
        double precision :: t1, t2, time1, time2, diff

        integer :: j, nr, nc
        double precision, dimension(n-1) :: xx, xx2, y1, y2
        double precision, dimension(n-1) :: c1,r1, c2, r2 
        double precision, dimension(n-1) :: cin,rin, cin2,rin2
        double precision :: logc12, logr12, h, hinv, gz0, gz1 


        nullify(x)
        nullify(y)
        nullify(b)
        nullify(c)
        nullify(d)

        nullify(z)
        nullify(gz)
        nullify(zb)
        nullify(zc)
        nullify(zd)

        allocate( x(n),y(n),  b(n),c(n),d(n), stat=ierr)
        isok = (ierr.eq.0)
        call assert( isok,                                    &          
     &      'setup_cauchy1: allocate x(n) failed, ierr=',ierr)

        x(:) = 0.0d0
        y(:) = 0.0d0
        b(:) = 0.0d0
        c(:) = 0.0d0
        d(:) = 0.0d0


        do i=1,n
          x(i) = xin(i)
          y(i) = yin(i)
        enddo


!       -----------------------------------------
!       z's located at mid point of [x(i),x(i+1)]
!       -----------------------------------------
        m = n-1
        allocate( z(m),gz(m), zb(m),zc(m),zd(m), stat=ierr )

        isok = (ierr.eq.0)
        call assert( isok,                                    &          
     &       'setup_cauchy1: allocate z(m) failed, ierr=',ierr)

        z(:) = 0.0d0
        gz(:) = 0.0d0
        zb(:) = 0.0d0
        zc(:) = 0.0d0
        zd(:) = 0.0d0


        do i=1,m
         z(i) = (x(i)+x(i+1))*0.5d0
        enddo


          call cpu_time(t1)
!         ----------------------------------------------------------
!         a0 = (y(j)*(x(j+1)-z(i)) - y(j+1)*(x(j)-z(i))  )/(x(j+1)-x(j))
!         gz0 = gz0 + a0*log( abs( (x(j+1)-zi)/(x(j)-zi) ) )
!         sum with 4 Toeplitz matrix multiply
!         ----------------------------------------------------------
          j = 1
          do i=1,n-1
            c1(i) = x(j+1) - z(i)
            c2(i) = x(j) - z(i)
!            ------------------------------
!            logc1(i) = log( abs(c1(i)) )
!            logc2(i) = log( abs(c2(i)) )
!            ------------------------------
          enddo
          i = 1
          do j=1,n-1
            r1(j) = x(j+1) - z(i)
            r2(j) = x(j) - z(i)
!            ----------------------------
!            logr1(j) = log( abs(r1(j)) )
!            logr2(j) = log( abs(r2(j)) )
!            ----------------------------
          enddo
!         ------------------------------------------
!         fast multiply of Toeplitz matrices
!
!         T1 = toepltiz(c1,r1) 
!         T2 = toeplitz(c2,r2)
!         then T1 .* T2 = toeplitz( c1.*c2, r1.*r2 )
!         T1 - T2 = toeplitz(c1-c2, r1-r2)
!         ------------------------------------------
          h = x(2)-x(1)
          hinv = dble(1)/dble(h)
          nr = n-1
          nc = n-1
          do i=1,n-1
!             ----------------------------
!             logc12 = logc1(i) - logc2(i)
!             logr12 = logr1(i) - logr2(i)
!             ----------------------------
             logc12 = log( abs(c1(i)/c2(i)) )
             logr12 = log( abs(r1(i)/r2(i)) )
             cin(i) = c1(i) * logc12
             rin(i) = r1(i) * logr12
             cin2(i) = c2(i) * logc12
             rin2(i) = r2(i) * logr12
             xx(i) = y(i)*hinv
             xx2(i) = y(i+1)*hinv
          enddo

          call toeplitz_mult(nc,cin,  nr, rin,  xx,  y1 )
          call toeplitz_mult(nc,cin2, nr, rin2, xx2, y2 )

          gz1 = y(n) - y(1)
          do i=1,n-1
            gz0 = y1(i) - y2(i)
            gz(i) = gz1 + gz0
          enddo

          call cpu_time(t2)
          time1 = t2-t1

       if (idebug .ge. 1) then

        call cpu_time(t1)
        do i=1,m
          zi = z(i)
          call hilb1( n,x,y,  zi,gzi )
          gz2(i) = gzi
        enddo
        call cpu_time(t2)
        time2 = t2-t1

        if (idebug.ge.2) then
        do i=1,m
          write(*,9010) i,gz(i),gz2(i),gz(i)-gz2(i)
 9010     format('setup_cauchy1:i,gz,gz2,diff ',i5,3(1x,1pe12.4))
        enddo
        endif


        diff = maxval( abs( gz - gz2 ) )
        write(*,9020) time1, time2, diff
 9020   format(1x,'setup_cauchy1: fast time1, time2, diff ',   &         
     &         3(1x,1pe12.4))

        endif



        call cpu_time(t1)
        call spline(size(z), z, gz, zb, zc, zd )
        call cpu_time(t2)
        if (idebug.ge.1) then
          write(*,9030) n, t2-t1
 9030     format(' setup_cauchy1,spline2:n, time ',i5,(1x,1pe12.4))
        endif

!       --------------------------------
!       save arrays for later evaluation
!       --------------------------------

        p%n = n

        p%x => x
        p%y => y
        p%b => b
        p%c => c
        p%d => d

        p%z => z
        p%gz => gz
        p%zb => zb
        p%zc => zc
        p%zd => zd


        return
        end subroutine setup_cauchy1
!23456789!23456789!23456789!23456789!23456789!23456789!23456789!23456789
        subroutine setup_cauchy3(n,xin,yin,  p)
        implicit none
!       -----------------------------------------------
!       precompute g(z) = integral(f(x)/(x-z),x=x1..xn)
!       and use faster interpolation
!       -----------------------------------------------
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) ::  xin,yin
        type(cauchy_t), intent(inout) :: p

        integer, parameter :: idebug = 0
        integer :: i, m,k
        double precision :: h, zi, gzi, t2, t1
        double precision, pointer, dimension(:) :: x,y,b,c,d
        double precision, pointer, dimension(:) :: z,gz,zb,zc,zd


        double precision, dimension(n-1) :: gz2,g0,g1,g2,g3,gtmp
        double precision, dimension(n-1) :: cin,rin,c1,r1,c2,r2
        double precision, dimension(n-1) :: c12,r12, logc12,logr12,yout

        double precision :: dz, a0,a1,a2,a3, dz0,dz1
        double precision :: diff, time1, time2, rc2, rc3
        integer :: j, nr, nc
!       ---------------------------------
!       generate spline approximation
!       
!       over [x(i),x(i+1)], the function is
!       approximated by the cubic
!
!       f(x) = y(i) + (x-x(i))*b(i) + (x-x(i))^2*c(i) + 
!              (x-x(i))^3*d(i)
!
!       ---------------------------------
        logical, parameter :: use_vlog = .true.

        integer :: icase, ilo,ihi, ierr 
        logical isok

        nullify(x)
        nullify(y)
        nullify(b)
        nullify(c)
        nullify(d)

        nullify(z)
        nullify(gz)
        nullify(zb)
        nullify(zc)
        nullify(zd)

        allocate( x(n),y(n),  b(n),c(n),d(n), stat=ierr)
        isok = (ierr.eq.0)
        call assert( isok,                                     &         
     &      'setup_cauchy3: allocate x(n) failed, ierr=',ierr)

        x(:) = 0.0d0
        y(:) = 0.0d0
        b(:) = 0.0d0
        c(:) = 0.0d0
        d(:) = 0.0d0


        do i=1,n
          x(i) = xin(i)
          y(i) = yin(i)
        enddo

!       -----------------------------------------
!       z's located at mid point of [x(i),x(i+1)]
!       k extra points on both ends
!       -----------------------------------------
        k = 0
        m = n-1
        allocate( z(m),gz(m),  stat=ierr )

        isok = (ierr.eq.0)
        call assert( isok,                                     &         
     &       'setup_cauchy3: allocate z(m) failed, ierr=',ierr)

        z(:) = 0.0d0
        gz(:) = 0.0d0

!       ---------------
!       setup interior points
!       ---------------
        do i=1,n-1
          z(k+i) = (x(i) + x(i+1))/2.0d0;
        enddo


        
        call cpu_time(t1)
        call spline( n, x,y,   b,c,d)
        call cpu_time(t2)
        if (idebug .ge. 1) then
          write(*,9010) n,t2-t1
 9010     format(' setup_cauchy3,spline1: n,time ',i5,1x,1pe12.4)
        endif
!
!       --------------------------------------------------
!       evaluate g(z) = integral( f(x)/(x-z), x=x(1)..x(n)
!       --------------------------------------------------


!       ----------------------------------
!       fast evaluation at interior points
!       ----------------------------------

           
!          --------------------------------------
!          dz = z - x(j)
!          a0 = y(j) + dz*( b(j) + dz*(c(j) + dz*d(j)))
!          a1 = b(j) + dz*(2.0d0*c(j) + 3.0d0*dz*d(j))
!          a2 = c(j) + 3.0d0*dz*d(j)
!          a3 = d(j)
!
!          a0 = y(j) + dz*b(j) + dz^2*c(j) + dz^3*d(j)
!          a1 = b(j) + 2*dz*c(j) + 3*dz^3*d(j)
!          a2 = c(j) + 3*dz*d(j)
!          a3 = d(j)
!
!          f(x) = a0 + a1*(x-z) + a2*(x-z)^2 + a3*(x-z)^3
!          g(z) = a0*log( (x(j+1)-z)/(x(j)-z) ) + 
!                    a1*(x(j+1)-x(j)) +
!                    (a2/2)*((x(j+1)-z)^2 - (x(j)-z)^2) +
!                    (a3/3)*((x(j+1)-z)^3 - (x(j)-z)^3)          
!          ---------------------------------------

        call cpu_time(t1)

        nr = n-1
        nc = n-1
        j = 1
        do i=1,nc
          c1(i) = z(i) - x(j)
          c2(i) = z(i) - x(j+1)
          if (use_vlog) then
            c12(i) = abs(c2(i)/c1(i))
          else
             logc12(i) = log( abs(c2(i)/c1(i)) )
          endif
        enddo

        if (use_vlog) then
          call vlog( logc12, c12, nc) 
        endif

        i = 1
        do j=1,nr
          r1(j) = z(i) - x(j)
          r2(j) = z(i) - x(j+1)
          if (use_vlog) then
            r12(j) = abs(r2(j)/r1(j))
          else
            logr12(j) = log( abs(r2(j)/r1(j)) )
          endif
        enddo

        if (use_vlog) then
          call vlog( logr12, r12, nr )
        endif

     

        do i=1,nc
           cin(i) = logc12(i)
        enddo
        do j=1,nr
           rin(j) = logr12(j)
        enddo
        call toeplitz_mult(nc,cin,nr,rin,y,yout)
!  g0(1:nc) = yout(1:nc)
        do i=1,nc
           g0(i) = yout(i)
        enddo

        do i=1,nc
          cin(i) = c1(i) * cin(i)
        enddo
        do j=1,nr
          rin(j) = r1(j) * rin(j)
        enddo
        call toeplitz_mult( nc, cin, nr, rin, b, yout)
!  g0(1:nc) = g0(1:nc) + yout(1:nc)
        do i=1,nc
          g0(i) = g0(i) + yout(i)
        enddo

        do i=1,nc
          cin(i) = c1(i)*cin(i)
        enddo
        do j=1,nr
          rin(j) = r1(j)*rin(j)
        enddo
        call toeplitz_mult(nc,cin,nr,rin, c, yout )
!  g0(1:nc) = g0(1:nc) + yout(1:nc)
        do i=1,nc
          g0(i) = g0(i) + yout(i)
        enddo


        do i=1,nc
          cin(i) = c1(i)*cin(i)
        enddo
        do j=1,nr
          rin(j) = r1(j)*rin(j)
        enddo
        call toeplitz_mult(nc,cin,nr,rin, d, yout)
!  g0(1:nc) = g0(1:nc) + yout(1:nc)
        do i=1,nc
           g0(i) = g0(i) + yout(i)
        enddo

        if (idebug.ge.2) then
!       --------
!       check g0
!       --------
        gtmp(1:nc) = 0.0d0
        do i=1,nc
          do j=1,nr
             dz = z(i) - x(j)
             dz1 = x(j+1)-z(i)
             dz0 = x(j)-z(i)


             a0 = y(j) + dz*( b(j) + dz*(c(j) + dz*d(j)))
             gtmp(i) = gtmp(i)+a0*log( abs(dz1/dz0) )
          enddo
        enddo

        write(*,*) 'check g0: diff ',maxval( abs(g0 - gtmp) )
        g0(1:nc) = gtmp(1:nc)
        endif


!          -----------------------------------
!          a1 = b(j) + 2*dz*c(j) + 3*dz^3*d(j)
!          contribution is a1*(x(j+1)-x(j))
!          -----------------------------------
        do j=1,nr
           rin(j) = x(j+1)-x(j)
        enddo
        do i=1,nc
           cin(i) = rin(i);
        enddo
        call toeplitz_mult(nc,cin,nr,rin,b,yout)
!  g1(1:nc) = yout(1:nc)
        do i=1,nc
          g1(i) = yout(i)
        enddo


        do i=1,nc
           cin(i) = c1(i)*cin(i)
        enddo
        do j=1,nr
           rin(j) = r1(j)*rin(j)
        enddo
        call toeplitz_mult(nc,cin,nr,rin,c,yout)
!  g1(1:nc) = g1(1:nc) + 2.0d0*yout(1:nc)
        do i=1,nc
          g1(i) = g1(i) + 2.0d0*yout(i)
        enddo


        do i=1,nc
           cin(i) = c1(i)*cin(i)
        enddo
        do j=1,nc
           rin(j) = r1(j)*rin(j)
        enddo
        call toeplitz_mult(nc,cin,nr,rin,d,yout)
!  g1(1:nc) = g1(1:nc) + 3.0d0*yout(1:nc)
        do i=1,nc
          g1(i) = g1(i) + 3.0d0*yout(i)
        enddo


        if (idebug.ge.2) then
!       --------
!       check g1
!       --------
        gtmp(:) = 0.0d0
        do i=1,nc
        do j=1,nr
            dz = z(i) - x(j)
            a1 = b(j) + dz*(2.0d0*c(j) + 3.0d0*dz*d(j))
            gtmp(i) = gtmp(i) + a1*(x(j+1)-x(j))
        enddo
        enddo

        write(*,*) 'check g1: diff ',maxval(abs(g1-gtmp))
        g1(1:nc) = gtmp(1:nc)
        endif


!          ---------------------------------------
!          a2 = c(j) + 3*dz*d(j)
!          contribution is
!          (a2/2.0d0)*( (x(j+1)-z)^2 - (x(j)-z)^2 )
!          (a2/2.0d0)*(x(j+1)-x(j))*(x(j+1)+x(j)-2.0d0*z)
!          ---------------------------------------

        do j=1,nr
! rin(j) = r2(j)*r2(j) - r1(j)*r1(j)
          rin(j) = (r2(j)-r1(j))*(r2(j)+r1(j))
        enddo
        do i=1,nc
! cin(i) = c2(i)*c2(i) - c1(i)*c1(i)
          cin(i) = (c2(i)-c1(i))*(c2(i)+c1(i))
        enddo
        call toeplitz_mult(nc,cin,nr,rin,c,yout)
!  g2(1:nc) = yout(1:nc)
        do i=1,nc
           g2(i) = yout(i)
        enddo

!  cin(1:nc) = cin(1:nc) * c1(1:nc)
        do i=1,nc
          cin(i) = cin(i) * c1(i)
        enddo

!  rin(1:nr) = rin(1:nr) * r1(1:nr)
        do j=1,nr
           rin(j) = rin(j) * r1(j)
        enddo

        call toeplitz_mult(nc,cin,nr,rin,d,yout)
!  g2(1:nc) = g2(1:nc) + 3.0d0*yout(1:nc)
        do i=1,nc
           g2(i) = g2(i) + 3.0d0*yout(i)
        enddo

        if (idebug.ge.2) then
!       --------
!       check g2
!       --------
        gtmp(1:nc) = 0.0d0
        do i=1,nc
        do j=1,nr
          dz = z(i) - x(j)
          dz1 = x(j+1)-z(i)
          dz0 = x(j)-z(i)

          a2 = c(j) + 3.0d0*dz*d(j)
          gtmp(i) = gtmp(i) +                                 &          
     &         a2*(x(j+1)-x(j))*(x(j+1)+x(j)-2.0d0*z(i))



        enddo
        enddo

        write(*,*) 'check g2: diff ', maxval(abs(gtmp-g2))
        do i=1,min(5,size(gtmp))
          write(*,*) 'i,gtmp,g2 ',i,gtmp(i),g2(i)
        enddo
        g2(1:nc) = gtmp(1:nc)

        endif



!          -----------------------------------
!          a3 = d(j)
!          dz1 = x(j+1) - z(i)
!          dz0 = x(j) - z(i)
!          contribution (a3/3)*(dz1^3 - dz0^3)
!          -----------------------------------

        j = 1
        do i=1,nc
          dz1 = x(j+1)-z(i)
          dz0 = x(j) - z(i)
! cin(i) = dz1*dz1*dz1 - dz0*dz0*dz0
          cin(i) = (dz1-dz0)*(dz1*dz1 + dz1*dz0 + dz0*dz0)
        enddo

        i = 1
        do j=1,nr
          dz1 = x(j+1)-z(i)
          dz0 = x(j) - z(i)
! rin(j) = dz1*dz1*dz1 - dz0*dz0*dz0
          rin(j) = (dz1-dz0)*(dz1*dz1 + dz1*dz0 + dz0*dz0)
        enddo
        call toeplitz_mult(nc,cin,nr,rin,d,yout)
!  g3(1:nc) = yout(1:nc)
        do i=1,nc
           g3(i) = yout(i)
        enddo


        if (idebug.ge.2) then
!       --------
!       check g3
!       --------
        gtmp(:) = 0.0d0
        do i=1,nc
        do j=1,nr
          dz = z(i) - x(j)
          dz1 = x(j+1)-z(i)
          dz0 = x(j)-z(i)
          a3 = d(j)
          gtmp(i) = gtmp(i) + a3*(x(j+1)-x(j))*               &          
     &                (dz1*dz1 + dz1*dz0 + dz0*dz0)
        enddo
        enddo

        write(*,*) 'check g3: diff ', maxval(abs(g3-gtmp))
        g3(1:nc) = gtmp(1:nc)
        endif

        
        rc2 = dble(1)/dble(2)
        rc3 = dble(1)/dble(3)
!        gz(1:nc) = g0(1:nc) + g1(1:nc) + rc2*g2(1:nc) +                   &
!     &                   rc3*g3(1:nc)
        do i=1,nc
         gz(i) = g0(i) + g1(i) + rc2*g2(i) +            &                
     &                   rc3*g3(i)
        enddo
          

        call cpu_time(t2)
        time1 = t2-t1

!       --------------------
!       evaluate at interior points
!       --------------------
        if (idebug.ge.1) then
        ilo = k+1
        ihi = ilo + (n-1) -1

        call cpu_time(t1)
        do i=ilo,ihi
           zi = z(i)
           call hilb3(n,x,y,b,c,d, zi, gzi )
           gz2(i) = gzi
        enddo
        call cpu_time(t2)
        time2 = t2-t1

        diff = maxval( abs(gz - gz2) )
        write(*,9020) n,time1, time2,diff
 9020   format(' setup_cauchy3:n, time1, time2, diff ',     &            
     &              i5,3(1x,1pe12.4))

        if (idebug.ge.2) then
          do i=1,min(5,size(gz))
            write(*,*) 'i,gz,gz2 ',i,gz(i),gz2(i)
          enddo
          gz(1:nc) = gz2(1:nc)
        endif


        endif

!       ----------------------------------------------
!       save results, extra evaluation near end points
!       ----------------------------------------------
        gtmp(1:nc) = gz(1:nc)
        deallocate( z, gz, stat=ierr)
        call assert(ierr.eq.0,'setup_cauchy3:deallocate ',ierr)
        call cpu_time(t2)
        time1 = t2-t1

        k = 8
        m = (n-1) + 2*k
        allocate( z(m),gz(m), zb(m),zc(m),zd(m), stat=ierr )

        z(:) = 0.0d0
        gz(:) = -987654321.0
        zb(:) = 0.0d0
        zc(:) = 0.0d0
        zd(:) = 0.0d0


!       ---------------
!       setup interior points
!       ---------------
        do i=1,n-1
          z(k+i) = (x(i) + x(i+1))/2.0d0;
        enddo
        h = (x(2)-x(1))
        h = (h/2.0d0)/dble(k+1)
        do i=k,1,-1
          z(i) = z(i+1)-h
        enddo
        h = x(n)-x(n-1)
        h = h/2.0d0/dble(k+1)
        do i=1,k
          z(k+(n-1)+i) = z(k+(n-1)+(i-1)) + h
        enddo

!       ------------------------
!       copy precomputed results
!       ------------------------
        ilo = k+1
        ihi = ilo + nc - 1
        gz( ilo:ihi ) = gtmp(1:nc)


!       ------------------
!       precompute answer 
!       ------------------
        call cpu_time(t1)
        do icase=1,2
          if (icase.eq.1) then
             ilo = lbound(z,1)
             ihi = ilo + k-1
          else
             ilo = ((k+1) + (n-1) - 1) + 1
             ihi = ubound(z,1)
          endif

          do i=ilo,ihi
           zi = z(i)
           gzi = 0.0d0

           call hilb3( n, x, y, b, c,d, zi, gzi )

           gz(i) = gzi
          enddo
        enddo
        call cpu_time(t2)
        if (idebug.ge.1) then
          write(*,9025) n, t2-t1
 9025     format(' setup_cauchy3,hilb3: n, time ',i5,1x,1pe12.4)
        endif

        call cpu_time(t1)
        call spline(size(z), z, gz, zb, zc, zd )
        call cpu_time(t2)
        if (idebug.ge.1) then
          write(*,9030) n,t2-t1
 9030     format(' setup_cauchy3,spline2: n, time ',i5,1x,1pe12.4)
        endif

!       --------------------------------
!       save arrays for later evaluation
!       --------------------------------

        p%n = n

        p%x => x
        p%y => y
        p%b => b
        p%c => c
        p%d => d

        p%z => z
        p%gz => gz
        p%zb => zb
        p%zc => zc
        p%zd => zd


        return
        end subroutine setup_cauchy3

        end module cauchy_mod
