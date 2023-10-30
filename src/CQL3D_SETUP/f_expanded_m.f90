MODULE f_expanded_m

        IMPLICIT NONE
                
        CHARACTER (LEN = 34) :: C_matrix_file

! Data read in from C_matrix_file

        character(31) :: in_file_CQL     ! name of Harvey's netCDF file

        INTEGER :: n_psi        ! Number of flux surfaces
        INTEGER :: n_u_coeffs, n_theta_coeffs
        REAL :: u_norm
        REAL :: u0_Max_poly     ! Width of Maxwellian used in Max_poly_u_basis
        REAL, ALLOCATABLE :: C_matrix(:,:,:)    ! Expansion coeffs C(i_psi, i_u, j_theta)
        REAL, ALLOCATABLE :: dens_Maxwell_psi(:)        ! density of Maxwellian component
        REAL, ALLOCATABLE :: T_Maxwell_psi(:)   ! Temperature in keV of Maxwellian component
        
CONTAINS
                
! ********************************************************************

SUBROUTINE read_C_matrix_m
        
        USE basis_functions_m, only : u_fns_name
        
        IMPLICIT NONE
        
        INTEGER, PARAMETER :: C_matrix_unit = 55
        INTEGER :: istat
                
        OPEN (UNIT=C_matrix_unit, FILE=C_matrix_file, FORM='FORMATTED', status = 'old')
        
        READ (C_matrix_unit, *) in_file_CQL
        READ (C_matrix_unit, *) u_fns_name
        
        READ (C_matrix_unit, *) n_psi
        READ (C_matrix_unit, *) n_u_coeffs
        READ (C_matrix_unit, *) n_theta_coeffs
                        
        ALLOCATE( C_matrix(n_psi, n_u_coeffs, n_theta_coeffs) , stat=istat )
        IF (istat /= 0 ) THEN
        WRITE (*,'("read_C_matrix: allocate failed for C_matrix")')
        ! PAUSE
        END IF
        
        ALLOCATE( dens_Maxwell_psi(n_psi), stat=istat )
    IF (istat /= 0 ) THEN
    WRITE (*,'("read_C_matrix: allocate failed for dens_Maxwell_psi")')
    ! PAUSE
    END IF
        
        ALLOCATE( T_Maxwell_psi(n_psi), stat=istat )
    IF (istat /= 0 ) THEN
    WRITE (*,'("read_C_matrix: allocate failed for T_Maxwell_psi")')
    ! PAUSE
    END IF
        
        READ (C_matrix_unit, '(6(1pe16.8))') C_matrix
        
        READ (C_matrix_unit, *) u_norm
        READ (C_matrix_unit, *) u0_Max_poly
        READ (C_matrix_unit, *) dens_Maxwell_psi
        READ (C_matrix_unit, *) T_Maxwell_psi

END SUBROUTINE read_C_matrix_m

                
! ********************************************************************

SUBROUTINE eval_f_exp_m(u, theta, i_psi, n_u_basis, n_theta_basis, n0, T, u_norm, f, &
                                                &       df_du, df_dtheta)

! Evaluates f_u_theta for arbitrary u, theta on flux surface i_psi by first calculating
! delta_f (summing expansion coefficients times the basis functions) then adding
! in the isotropic Maxwellian component for this flux surface.

! n_u_basis, and n_theta_basis are the number of terms in each basis set actually summed.
! This has to be ² number of coefficents in the C_matrix. n_u_basis ² n_u_coeffs
! The user is responsible for making sure this is the case. He should check outside this routine.

        USE basis_functions_m, only : u_basis, cos_theta_basis

        IMPLICIT NONE
        
! Declare arguments
        REAL, INTENT(in) :: u, theta, n0, T, u_norm
        INTEGER, INTENT(in) :: n_u_basis, n_theta_basis
        INTEGER, INTENT(in) :: i_psi
        REAL, INTENT(out) :: f
        REAL, OPTIONAL, INTENT(out) :: df_du
        REAL, OPTIONAL, INTENT(out) :: df_dtheta

! Local variables

        REAL :: u_fns(n_u_basis), theta_fns(n_theta_basis)
        
        REAL :: du_fns_du(n_u_basis), dtheta_fns_dtheta(n_theta_basis)

        REAL :: delta_f, ddelta_f_du, f_Maxwell, df_Maxwell_du
        
        CALL u_basis(u, n_u_basis, u_fns, du_fns_du)
                        
        CALL cos_theta_basis(theta, n_theta_basis, theta_fns, dtheta_fns_dtheta)
        
        delta_f = DOT_PRODUCT( u_fns, MATMUL(C_matrix(i_psi, 1:n_u_basis, 1:n_theta_basis), &
                        &        theta_fns) )
                        
        CALL f_Maxwell_3D(u, u_norm, n0, T, f_Maxwell, df_Maxwell_du)
        
        f = u_norm**3*f_Maxwell + delta_f

        IF ( .NOT. PRESENT(df_du) ) RETURN
        
                ddelta_f_du = DOT_PRODUCT( du_fns_du, &
                                        & MATMUL(C_matrix(i_psi, 1:n_u_basis, 1:n_theta_basis), theta_fns) )
                                        
                df_du = u_norm**3*df_Maxwell_du + ddelta_f_du

        IF ( .NOT. PRESENT(df_dtheta) ) RETURN
                
                df_dtheta = DOT_PRODUCT( u_fns, MATMUL(C_matrix(i_psi, 1:n_u_basis, 1:n_theta_basis), &
                                &        dtheta_fns_dtheta) )
                                
        RETURN

END SUBROUTINE eval_f_exp_m

! *******************************************************************   

                
! ********************************************************************

SUBROUTINE eval_delta_f_exp_m(u, theta, i_psi, n_u_basis, n_theta_basis, f, &
                                                &       df_du, df_dtheta)

! Evaluates delta_f_u_theta by summing expansion coefficients times the basis functions
! for arbitrary u, theta

! n_u_basis, and n_theta_basis are the number of terms in each basis set actually summed.
! This has to be ² number of coefficents in the C_matrix. n_u_basis ² n_u_coeffs
! The user is responsible for making sure this is the case. He should check outside this routine.

        USE basis_functions_m, only : u_basis, cos_theta_basis

        IMPLICIT NONE
        
! Declare arguments
        REAL, INTENT(in) :: u, theta
        INTEGER, INTENT(in) :: n_u_basis, n_theta_basis
        INTEGER, INTENT(in) :: i_psi
        REAL, INTENT(out) :: f
        REAL, OPTIONAL, INTENT(out) :: df_du
        REAL, OPTIONAL, INTENT(out) :: df_dtheta

! Local variables

        REAL :: u_fns(n_u_basis), theta_fns(n_theta_basis), delta_f
        
        REAL :: du_fns_du(n_u_basis), dtheta_fns_dtheta(n_theta_basis)
        
        CALL u_basis(u, n_u_basis, u_fns, du_fns_du)
                        
        CALL cos_theta_basis(theta, n_theta_basis, theta_fns, dtheta_fns_dtheta)
        
        delta_f = DOT_PRODUCT( u_fns, MATMUL(C_matrix(i_psi, 1:n_u_basis, 1:n_theta_basis), &
                        &        theta_fns) )
        
        f = delta_f

        IF ( .NOT. PRESENT(df_du) ) RETURN
        
                df_du = DOT_PRODUCT( du_fns_du, &
                                        & MATMUL(C_matrix(i_psi, 1:n_u_basis, 1:n_theta_basis), theta_fns) )
                                        
        IF ( .NOT. PRESENT(df_dtheta) ) RETURN
                
                df_dtheta = DOT_PRODUCT( u_fns, MATMUL(C_matrix(i_psi, 1:n_u_basis, 1:n_theta_basis), &
                                &        dtheta_fns_dtheta) )
                                
        RETURN

END SUBROUTINE eval_delta_f_exp_m

! *******************************************************************   

END MODULE f_expanded_m

! *******************************************************************   

        SUBROUTINE f_Maxwell_3D(v, v_norm, n0, T, f, df_dv)
! Gives 3D Maxwellian distribution of n0 and T
! v is velocity in units of vnorm which is in cm/sec
! n0 is in cm^-3
! T is in keV

        USE global_data_m, only : pi, kev_per_erg, mass_D
                
        IMPLICIT none
        
        REAL, INTENT(IN) :: v, v_norm, n0, T
        REAL, INTENT(OUT) :: f
!       REAL, OPTIONAL, INTENT(OUT) :: df_dv
        REAL, INTENT(OUT) :: df_dv
        
        REAL :: v0
        
        v0 = SQRT( 2.*T/kev_per_erg/mass_D )
        
        f = n0/pi/SQRT(pi)/v0**3 * exp( -(v*v_norm)**2/v0**2 )
        
!       IF (.not.PRESENT(df_dv)) RETURN
        
        df_dv = -2.0*v*v_norm**2/v0**2*f
        
        END SUBROUTINE f_Maxwell_3D



