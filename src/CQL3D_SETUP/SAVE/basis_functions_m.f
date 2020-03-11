

MODULE basis_functions_m

	IMPLICIT NONE

	CHARACTER(len = 20) :: u_fns_name		! kind of u basis functions to use.  select one of:
											! "polynomial_u", "Max_poly_u", or "B_spline_u"

	REAL :: u0_Max_poly	! Width of Maxwellian used in Max_poly_u_basis
	
CONTAINS

! *******************************************************************	



SUBROUTINE u_basis(u, n_g, g, dg_du)

! Evaluates the g(u) basis functions (n_g of them) at argument u

	IMPLICIT NONE
	
	REAL, INTENT(in) :: u
	INTEGER, INTENT(in) :: n_g
	REAL, INTENT(out) :: g(n_g)
	REAL, OPTIONAL, INTENT(out) :: dg_du(n_g)

	SELECT CASE (TRIM(u_fns_name))
		
		CASE ("polynomial_u")
			CALL polynomial_u_basis(u, n_g, g, dg_du)
		CASE ("Max_poly_u")
			CALL Max_poly_u_basis(u, n_g, g)
		CASE ("B_spline_u")
			CALL cubic_B_spline_basis(u, n_g, g, dg_du)
		CASE DEFAULT
			WRITE(*,*) "Unimplemented u basis = ", u_fns_name
			STOP
			
	END SELECT
!              write(6, *) "dg_du = "
!              write(6, *) dg_du
	RETURN
	
END SUBROUTINE u_basis	


! *******************************************************************	


SUBROUTINE polynomial_u_basis(u, n_g, g, dg_du)

! Generates powers of u, u**n where n ranges from 0 to n_g - 1

	IMPLICIT NONE
	
	REAL, INTENT(in) :: u
	INTEGER, INTENT(in) :: n_g
	REAL, INTENT(out) :: g(n_g)
	REAL, OPTIONAL, INTENT(out) :: dg_du(n_g)
	
	INTEGER :: i
	
	g(1) = 1.
	
	DO i = 2, n_g
		g(i) = u * g(i-1)
	END DO
	
	IF (.not.PRESENT(dg_du))	RETURN
	
! Derivatives
	
	DO i = 1, n_g
		dg_du(i) = (i-1) * g(i-1)
	END DO

END SUBROUTINE polynomial_u_basis	


! *******************************************************************	


SUBROUTINE Max_poly_u_basis(u, n_g, g)

! Generates powers of u times a Gaussian of width u_th

	IMPLICIT NONE
	
	REAL, INTENT(in) :: u
	INTEGER, INTENT(in) :: n_g
	REAL, INTENT(out) :: g(n_g)
	REAL :: Max_factor
	
	INTEGER :: i
	
	g(1) = 1.
	
	DO i = 2, n_g
		g(i) = u * g(i-1)
	END DO
	
	Max_factor = EXP(-u**2/u0_Max_poly**2)
	
	g = g * Max_factor
	
	RETURN
	
END SUBROUTINE Max_poly_u_basis	


! *******************************************************************	






SUBROUTINE cos_theta_basis(theta, n_cosn, cosn, dcosn_dtheta)

! Generates cos(n * theta) where n ranges from 0 to n_g-1

	IMPLICIT NONE
	
	REAL, INTENT(in) :: theta
	INTEGER, INTENT(in) :: n_cosn
	REAL, INTENT(out) :: cosn(n_cosn)
	REAL, OPTIONAL, INTENT(out) :: dcosn_dtheta(n_cosn)
	
	INTEGER :: i
	
	cosn(1) = 1.
	
	DO i = 2, n_cosn
		cosn(i) = COS((i-1) * theta)
	END DO
	
	IF (.not.PRESENT(dcosn_dtheta))	RETURN
	
	dcosn_dtheta(1) = 0.
	
	DO i = 2, n_cosn
		dcosn_dtheta(i) = -(i-1)*SIN((i-1) * theta)
	END DO
	
	RETURN
	
END SUBROUTINE cos_theta_basis	


! *******************************************************************	

END MODULE basis_functions_m

