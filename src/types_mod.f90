module types
!jcwright: this is not used in aorsa, requires fortran 2003
!use, intrinsic :: iso_fortran_env !conflicts with -r8
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
