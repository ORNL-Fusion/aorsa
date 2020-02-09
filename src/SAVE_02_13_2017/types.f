module types
use, intrinsic :: iso_fortran_env
implicit none
integer, parameter :: sp = REAL(32)
integer, parameter :: dp = REAL(64)
! there is a bug in gfortran that
! for now, make the REAL128 returns extended
! precsion, REAL80
!integer, parameter :: qp = REAL128
integer, parameter :: qpp = selected_real_kind(32)
end module types
