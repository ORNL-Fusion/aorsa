module types
use, intrinsic :: iso_fortran_env
implicit none
!integer, parameter :: sp = REAL32
!integer, parameter :: dp = REAL64
integer, parameter :: dp = selected_real_kind(p=13,r=200)
! there is a bug in gfortran that
! for now, make the REAL128 returns extended
! precsion, REAL80
!integer, parameter :: qp = REAL128
integer, parameter :: qpp = selected_real_kind(32)
end module types
