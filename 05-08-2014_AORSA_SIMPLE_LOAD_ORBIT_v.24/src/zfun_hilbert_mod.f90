module zfun_hilbert
use types, only : rprec => dp
private
public :: zfun_hil_init, zfun_hil

! to get a double precision version, set n = 30 (or so)
! residuals are ~ 1.0e-14
integer, parameter ::  n = 32 ! number of eigenfunctions used in expansion
integer :: ierr ! error flag for fftpack
integer, parameter :: mult = 2! sample multiplyer--default is 2
integer, parameter :: m = mult * n! number of velocity function sample
integer, parameter :: lensav = 2 * m + int(log(real(2 * m))/log(2.0)) + 4
real(kind = rprec), parameter :: inv_m = 1.0_rprec / real(m, kind = rprec)
real(kind = rprec), dimension(2 * m) :: work
real(kind = rprec), dimension(lensav) ::  wsave
real(kind = rprec) :: length, & ! renormalization parameter--set by n
  d_theta ! theta grid increment
logical, save :: first_flag ! use for initialization
real(kind = rprec), parameter :: infinity = huge(1.0_rprec)
real(kind = rprec), parameter :: pi = 4.0_rprec * atan(1.0_rprec)
real(kind = rprec), dimension(:), allocatable :: &
  f_vel, & ! array of f values on velocity grid
  theta_grid, & ! theta grid-- -pi to pi - dtheta
  veloc_grid    ! grid of velocity values from vel = length * tan(theta/2.0)
                ! first value is -infinity, not actuall used--value of function will
                ! be set to zero.
complex(kind = rprec), dimension(:,:), allocatable ::   exp_m ! precomputed exponentials
complex(kind = rprec), dimension(:), allocatable :: an !array for eigenfunction expansion coefficents
real(kind = rprec), dimension(:), allocatable :: re_an
complex(kind = rprec), parameter :: zi = cmplx(0.0_rprec, 1.0_rprec)
contains
subroutine zfun_hil_init(ieer)
integer :: i
!print*, 'lensav = ', lensav, 'fft len = ', m
call rfft1i(2*m, wsave, lensav, ieer)
!print *, 'fft init '
!print *,  2*m, lensav
!do i = 1, lensav
!    print*, wsave(i)
!end do

!if(ieer .eq. 2) then
!    print *, 'input parameter lensav not big enough'
!    return
!end if
! move initial calculations to init
allocate(theta_grid(0:2 * m))
! a standard [0, 2*pi] grid will be used with an
! even number of intervals and an odd number to
! grid points.  all will be indexed [0,2*m]
! but the last on won't be used because of periodicity.
!print *, 'theta grid', 'points = ', 2*m + 1
d_theta = pi * inv_m
do i = 0, 2 * m
    theta_grid(i) = real(i, rprec) * d_theta
    !print *, i, theta_grid(i)
end do
! set the renormaliztion parameter
length = sqrt(real(n, rprec) * 0.7_rprec)
! fill the velocity grid
allocate(veloc_grid(0:2 * m ))
veloc_grid(0) = -infinity
veloc_grid(2 * m) = infinity
veloc_grid(1:2 * m - 1 ) = length * tan((theta_grid(1:2 * m - 1)  - pi) * 0.5_rprec)

! fill the function grid
allocate(f_vel(0:2 * m ))
! may need to truncate the fill for practical cases
! put zeros at the ends
f_vel(0) = 0.0_rprec
f_vel(2 * m) = 0.0_rprec
! and the rest
f_vel(1:2 * m -1) = exp(- veloc_grid(1:2 * m-1)**2) *  &
    (length**2 + veloc_grid(1:2 * m-1)**2)
!last term is takes out the weighting function in
!the hilbert transorm
!f_vel = 1.0_rprec
!do i = 0, 2 * m
!    print *, i, theta_grid(i), veloc_grid(i), f_vel(i)
!end do
!end do
!print *, ''
! ft f_vel
! compute the exponentials (not needed in general when using fft)
!allocate(exp_m(0:2 * m - 1, 0:2 * m - 1))
!do i = 0, 2 * m -1
!    do j = 0, 2 * m -1
!        exp_m(i, j) = exp(- zi * j * theta_grid(i))
!        !if(j .eq. 0) print*, exp_m(i, j)
!    end do
!end do


! get the coefficients--will use fft of real function ultimately
allocate(an(0:2 * m))
allocate(re_an(0:m))


an = cmplx(0.0_rprec, 0.0_rprec)
!do j = 0, 2 * m - 1
!    do i = 0, 2 * m - 1
!        an(j) = an(j) + exp_m(i, j) * f_vel(i)
!        ! for odd harmonics need to give the fft negative values
!        !if(j.eq.0) print*, an(j)
!    end do
!end do
! try and make the fftpack call
!print *, 'ier for fftf call = ', ieer!, sum(f_vel(0:2*m -1))/2.0_rprec/m
!print *, 'size of work ', size(work), 2 * m, size(f_vel(0:2 * m -1))
call rfft1f(2 * m, 1, f_vel(0:2 * m -1), 2 * m, wsave, lensav, &
    work, 2 * m, ieer)
!print *, 'ier for fftf call = ', ieer, sum(f_vel)
!print *, 'long, real fft ', j, an(0)* 0.5_rprec * inv_m, f_vel(0)
do j = 1, m
!    print *, 'long, real fft ', j, an(j) * inv_m, f_vel(2*j-1), f_vel(2*j)
!    an(j) = (-1)**j * an(j) * 0.5_rprec * inv_m
    !an(j) = (-1)**j * f_vel(2*j-1) * 0.5_rprec
    re_an(j) = (-1)**j * f_vel(2*j-1) * 0.5_rprec
    !  need to adjust signs because of coordinate shift from
    !  [-pi, pi] to [0, 2*pi] and remove the 2*m dtf normalization
    !  don't need to higher j's since they are aliased to negative
    !  frequencies
end do
!an(0) = an(0) * 0.5_rprec * inv_m
!an(m) = an(m) * 2.0_rprec
re_an(m) = 2.0_rprec * re_an(m)
!an(0) = cmplx(f_vel(0), kind = rprec)
re_an(0) = f_vel(0)
!re_an = real(an(0:m))
!print *, 'fourier coefficients'

return
end subroutine zfun_hil_init

subroutine zfun_hil(z_in, z_out, flag)
!  This routine computes the 'Z' function using an expansion of the
!  in the distribution function in eigenfunctions of the
!  eigenfunctions of the Hilbert transform.  These eigen functions,
!  when the integration variable is changed from +/1 infinity via a
!  y = tan(theta/2.0), are Fourier basis functions on the interval
!  [-pi, pi].  With this renormalized integration variable, the
!  coefficients in the eigenvalue expansions can be computed with
!  FFTs.  In addition, with the further normalizaton z = L*tan(theta/2),
!  the convergence can be made competitive with the classical three-region
!  algogrithm presented by GPM Poppe, CMJ Wigers:  More efficient
!  computation of the complex error function ACM Trans. Math. Software
!  Vol. 16, No. 1, pp.47 presented by Hammett, netlib, and used (until replaced
!  because of license restrictions) in scipy.  While the advantages maybe small
!  for CPU Maxwellian calcualtions, there are no if tests, except perhaps
!  for faults.  The present routine is only good in the upper half plane
!  including the real axis
!  version2--convert to ffts, first use fftpack
implicit none
complex(kind = rprec), intent(in) :: z_in
complex(kind = rprec), intent(out) :: z_out
logical, intent(out) :: flag
integer :: i   ! dummy index
complex(kind = rprec) :: z_poly, z_term, z_out2
real(kind = rprec) :: x_term, y_term, x_poly, y_poly, x_in, y_in, x_out, y_out, &
    x_sq, y_sq, len_sq, denom_in, denom, x_out_sav, z_mod_sq, denom_l_in
! check for allocations
z_term = (length + zi * z_in) / (length - zi * z_in)
!x_in = real(z_in); y_in = aimag(z_in)
!x_sq = x_in **2; y_sq = y_in**2; len_sq = length**2
!z_mod_sq = x_sq + y_sq
!denom = z_mod_sq + len_sq + 2.0_rprec * length * y_in
!denom_l_in = 1.0_rprec /( z_mod_sq * (x_sq - y_sq) - 2.0_rprec * len_sq * x_in * y_in)

!denom_in = 1.0_rprec / denom
!x_term = (len_sq - z_mod_sq) * denom_in
!y_term = 2.0_rprec * x_in * length * denom_in
!print *, z_term, x_term, y_term
z_out = re_an(m)
!x_out = re_an(m)
!y_out = 0.0_rprec

do i = m, 2, -1
    z_out = z_term * z_out + re_an(i -1)
    !x_out_sav = x_out
    !x_out = re_an(i-1) + x_term * x_out_sav - y_term * y_out
    !y_out = x_term * y_out + y_term * x_out_sav
end do
!print *, 'long/short ', z_out, x_out, y_out
!z_out = cmplx(x_out, y_out, kind = rprec)
z_out = 2.0_rprec * z_term * z_out / (length**2 + z_in**2)
z_out = z_out + re_an(0) / length / (length - zi * z_in)
!x_out_sav = x_out
!x_out = 2.0_rprec * (x_term * x_out - y_term * y_out) * denom_l_in
!y_out = 2.0_rprec * (x_term * y_out + y_term * x_out_sav) * denom_l_in
!x_out_sav = x_out
!x_out = x_out * (len_sq + x_sq + y_sq)
!y_out = -y_out * (x_in * y_in)
!x_out_sav = x_out
!x_out = x_out + re_an(0) / length /
!y_out = y_out * re_an(0) / length



!print *, an(0) / length * sqrt(pi)
!x_out_sav = x_out
!x_out = -sqrt(pi) * y_out
!Y_out = sqrt(pi) * x_out_sav
z_out = zi * sqrt(pi) * z_out
! set the flag for success
flag = .true.
if(aimag(z_in) .lt. 0.0_rprec) stop 'negative im(zeta)'
return


end subroutine zfun_hil

end module zfun_hilbert

