       MODULE CQL_kinds_m

!----------------------------------------------------------------------
!  Kind specifications
!-------------------------------------------------------------------------------
!For REAL variables:
!  Use p=6, r=35   for single precision on 32 bit machine
!  Use p=12,r=100  for double precision on 32 bit machine
!                  for single precision on 64 bit machine
!Parameters for SELECTED_REAL_KIND:
!  p                   -number of digits
!  r                   -range of exponent from -r to r
!
! For INTEGER variables:
!	Use r=8 for single precision on 32 bit machine
!-------------------------------------------------------------------------------

	IMPLICIT NONE

	PRIVATE

! Defaults
      INTEGER, PARAMETER :: default_real = kind(1.0)
      INTEGER, PARAMETER :: default_int = kind(1)
      INTEGER, PARAMETER :: default_complex = kind(1.0)

! *8 kinds

      REAL*8, PARAMETER :: xdbl_real = 1.
      INTEGER*8, PARAMETER :: xdbl_int = 1
      INTEGER, PARAMETER :: dbl_real = kind(xdbl_real)
      INTEGER, PARAMETER :: dbl_int = kind(xdbl_int)
      INTEGER, PARAMETER :: dbl_complex = kind(xdbl_real)

! Selected kinds

       INTEGER, PARAMETER :: rspec = SELECTED_REAL_KIND(p=12,r=100)
       INTEGER, PARAMETER :: ispec = SELECTED_INT_KIND(r=100)

! Exported kinds

      INTEGER, PUBLIC, PARAMETER :: real_kind = default_real
      INTEGER, PUBLIC, PARAMETER :: int_kind = default_int
      INTEGER, PUBLIC, PARAMETER :: comp_kind = default_complex

       END MODULE CQL_kinds_m
