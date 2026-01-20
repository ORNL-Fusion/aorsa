module read_4df_mod
  use netcdf
  implicit none
contains

  !-----------------------------------------------------------------
  ! Read full 4D distribution (allocates outputs)
  ! f_rzvv(R, Z, vPer, vPar)
  !-----------------------------------------------------------------
  subroutine netcdfr4d_single_gen(netcdfnm, R_binCenters, z_binCenters, &
                                  vPer_binCenters, vPar_binCenters, f_rzvv)
    implicit none

    character(len=*), intent(in)  :: netcdfnm
    real,    allocatable, intent(out) :: R_binCenters(:)
    real,    allocatable, intent(out) :: z_binCenters(:)
    real,    allocatable, intent(out) :: vPer_binCenters(:)
    real,    allocatable, intent(out) :: vPar_binCenters(:)
    real,    intent(out) :: f_rzvv(:,:,:,:)

    integer :: ncid, vid, istatus
    integer :: n_vpar, n_vper, n_z, n_r

    ! open
    istatus = nf90_open(trim(netcdfnm), NF90_NOWRITE, ncid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Failed to open file: ', trim(netcdfnm)
       stop
    end if

    ! get dims
    call get_dimlen(ncid, 'vPar_binCenters', n_vpar)
    call get_dimlen(ncid, 'vPer_binCenters', n_vper)
    call get_dimlen(ncid, 'z_binCenters',    n_z)
    call get_dimlen(ncid, 'R_binCenters',    n_r)

    ! allocate outputs
    if (.not. allocated(vPar_binCenters)) allocate(vPar_binCenters(n_vpar))
    if (.not. allocated(vPer_binCenters)) allocate(vPer_binCenters(n_vper))
    if (.not. allocated(z_binCenters))   allocate(z_binCenters(n_z))
    if (.not. allocated(R_binCenters))   allocate(R_binCenters(n_r))
   !  if (.not. allocated(f_rzvv)) allocate(f_rzvv(n_r, n_z, n_vper, n_vpar))

    ! read coordinate variables
    call get_var1d(ncid, 'vPar_binCenters', vPar_binCenters)
    call get_var1d(ncid, 'vPer_binCenters', vPer_binCenters)
    call get_var1d(ncid, 'z_binCenters',    z_binCenters)
    call get_var1d(ncid, 'R_binCenters',    R_binCenters)

    ! read distribution
    istatus = nf90_inq_varid(ncid, 'f_rzvv', vid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Could not get varid for f_rzvv'
       istatus = nf90_close(ncid)
       if (istatus /= nf90_noerr) then
       print *, 'ERROR: nf90_close failed for ncid=', ncid
       stop
       end if
       stop
    end if

    istatus = nf90_get_var(ncid, vid, f_rzvv)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Could not read f_rzvv'
       istatus = nf90_close(ncid)
       if (istatus /= nf90_noerr) then
       print *, 'ERROR: nf90_close failed for ncid=', ncid
       stop
       end if
       stop
    end if

    istatus = nf90_close(ncid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: nf90_close failed'
       stop
    end if

  end subroutine netcdfr4d_single_gen


  !-----------------------------------------------------------------
  ! Get shapes: returns n_r, n_z, n_vper, n_vpar
  !-----------------------------------------------------------------
  subroutine get_dims_from_file(netcdfnm, n_r, n_z, n_vper, n_vpar)
    implicit none

    character(len=*), intent(in)  :: netcdfnm
    integer, intent(out) :: n_r, n_z, n_vper, n_vpar
    integer :: ncid, istatus

    istatus = nf90_open(trim(netcdfnm), NF90_NOWRITE, ncid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Failed to open NetCDF file: ', trim(netcdfnm)
       stop
    end if

    call get_dimlen(ncid, 'R_binCenters',    n_r)
    call get_dimlen(ncid, 'z_binCenters',    n_z)
    call get_dimlen(ncid, 'vPer_binCenters', n_vper)
    call get_dimlen(ncid, 'vPar_binCenters', n_vpar)

    istatus = nf90_close(ncid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: nf90_close failed in get_dims_from_file'
       stop
    end if

  end subroutine get_dims_from_file


  !-----------------------------------------------------------------
  ! Read a 1D variable into an allocatable array (caller supplies array)
  ! intent(inout) used so caller may pre-allocate or expect it to be filled
  !-----------------------------------------------------------------
  subroutine get_var1d(ncid, varname, array)
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varname
    real, intent(out) :: array(:)
    integer :: vid, istatus

    istatus = nf90_inq_varid(ncid, varname, vid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Could not get varid for '//trim(varname)
       stop
    end if

    istatus = nf90_get_var(ncid, vid, array)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Failed to read '//trim(varname)
       stop
    end if
  end subroutine get_var1d

  ! working on a 2D density read input 
  subroutine get_density_RZ(netcdfnm, n_r, n_z, density2d) 
    implicit none 
    integer, intent(in) :: n_r, n_z
    integer :: n_rcheck, n_zcheck
    character(len=*), intent(in)  :: netcdfnm
    real,  intent(out) :: density2d(:,:)
    integer :: ncid, vid, istatus
    ! open
    istatus = nf90_open(trim(netcdfnm), NF90_NOWRITE, ncid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Failed to open file: ', trim(netcdfnm)
       stop
    end if   
    ! check the dims of the density array 
    call get_dimlen(ncid, 'z_binCenters',    n_zcheck)
    call get_dimlen(ncid, 'R_binCenters',    n_rcheck)

    if (n_zcheck .ne. n_z) then
      error stop "density_rz z dim is incorrect"
    end if

    if (n_rcheck .ne. n_r) then
      error stop "density_rz r dim is incorrect"
    end if

    ! raise an error if these do not match the expected dimensions 
    ! read density 
    istatus = nf90_inq_varid(ncid, 'density_rz', vid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Could not get varid for density_rz'
       istatus = nf90_close(ncid)
       if (istatus /= nf90_noerr) then
       print *, 'ERROR: nf90_close failed for ncid=', ncid
       stop
       end if
       stop
    end if

    istatus = nf90_get_var(ncid, vid, density2d)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Could not read density_rz'
       istatus = nf90_close(ncid)
       if (istatus /= nf90_noerr) then
       print *, 'ERROR: nf90_close failed for ncid=', ncid
       stop
       end if
       stop
    end if

    istatus = nf90_close(ncid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: nf90_close failed'
       stop
    end if
  end subroutine get_density_RZ




   
  !-----------------------------------------------------------------
  ! Query variable to get its first dimension length (works for 1D var)
  !-----------------------------------------------------------------
  subroutine get_dimlen(ncid, varname, len_out)
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varname
    integer, intent(out) :: len_out

    integer :: vid, istatus, ndims
    integer, allocatable :: dimids(:)

    istatus = nf90_inq_varid(ncid, varname, vid)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Variable not found - ', trim(varname)
       stop
    end if

    ! ask number of dimensions
    istatus = nf90_inquire_variable(ncid, vid, ndims=ndims)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Could not get ndims for '//trim(varname)
       stop
    end if

    if (ndims < 1) then
       print *, 'ERROR: variable '//trim(varname)//' has ndims < 1'
       stop
    end if

    allocate(dimids(ndims))

    istatus = nf90_inquire_variable(ncid, vid, dimids=dimids)
    if (istatus /= nf90_noerr) then
       print *, 'ERROR: Could not get dim IDs for '//trim(varname)
       deallocate(dimids)
       stop
    end if

    ! Take the first dimid (for 1D coordinate variables)
    istatus = nf90_inquire_dimension(ncid, dimids(1), len=len_out)
    if (istatus /= nf90_noerr) then
    print *, 'ERROR: nf90_inquire_dimension failed for dimid=', dimids(1)
    stop
    end if

    deallocate(dimids)
  end subroutine get_dimlen


  subroutine deriv_dfdvperp(f_rzvv, dfdvper_rzvv, i_r, i_z, &
                            vPer_rzvv, n_r, n_z, n_vper, n_vpar)

   ! JVDL calculates the df/dvper derivatives at R,Z. 
    implicit none 

    integer, intent(in) :: i_r, i_z, n_vper, n_vpar, n_r, n_z
    real, intent(in) :: vPer_rzvv(n_vper)
    real, intent(in) :: f_rzvv(n_r, n_z, n_vper, n_vpar)

    real, intent(out) ::  dfdvper_rzvv(n_r, n_z, n_vper, n_vpar)

    ! local vars 
    integer :: ivpar, ivper

    
    do ivpar = 1, n_vpar
      do ivper = 1, n_vper
         ! interior vperp derivatives centered difference 
         if(ivper .ne. 1 .and. ivper .ne. n_vper)then
            dfdvper_rzvv(i_r, i_z, ivper, ivpar) = (f_rzvv(i_r, i_z, ivper+1, ivpar) - &
                        f_rzvv(i_r, i_z, ivper-1, ivpar)) / &
                        (vPer_rzvv(ivper+1) -  vPer_rzvv(ivper-1))

   ! not sure i even need to calculate the second derivative? 
   !          d2fdr2 = (f_rzvv(i_r, i_z, ivper+1, i_vpar) - 2.0 * f_rzvv(i_r, i_z, ivper,i_vpar) + f_rzvv(i_r, i_z, ivper-1, i_vpar)) / &
   !   &              (vPer_rzvv(i_vper+1) - vPer_rzvv(i_vper))**2
         end if

         ! BC at vperp = 0
         if(ivper .eq. 1)then
             dfdvper_rzvv(i_r, i_z, ivper, ivpar) = 0.0
         end if 

         ! BC at vperp = n_vperp, the largest vperp on the velocity grid. TODO confirm 2nd derivs not needed. 
         if(ivper .eq. n_vper)then
            dfdvper_rzvv(i_r, i_z, ivper, ivpar) = 0.0
         end if

         end do
      end do
   end subroutine deriv_dfdvperp 

   subroutine deriv_dfdvpar(f_rzvv, dfdvpar_rzvv, i_r, i_z, &
                            vPar_rzvv, n_r, n_z, n_vper, n_vpar)

   ! JVDL calculates the df/dvpar derivatives at R,Z. 
    implicit none 

    integer, intent(in) :: i_r, i_z, n_vper, n_vpar, n_r, n_z
    real, intent(in) :: vPar_rzvv(n_vpar)
    real, intent(in) :: f_rzvv(n_r, n_z, n_vper, n_vpar)

    real, intent(out) ::  dfdvpar_rzvv(n_r, n_z, n_vper, n_vpar)

    ! local vars 
    integer :: ivpar, ivper

    
    do ivper = 1, n_vper
      do ivpar = 1, n_vpar
         ! interior vpar derivatives centered difference 
         if(ivpar .ne. 1 .and. ivpar .ne. n_vpar)then
            dfdvpar_rzvv(i_r, i_z, ivper, ivpar) = (f_rzvv(i_r, i_z, ivper, ivpar+1) - &
                        f_rzvv(i_r, i_z, ivper, ivpar-1)) / &
                        (vPar_rzvv(ivpar+1) -  vPar_rzvv(ivpar-1))

   ! not sure i even need to calculate the second derivative? 
   !          d2fdr2 = (f_rzvv(i_r, i_z, ivper+1, ivpar) - 2.0 * f_rzvv(i_r, i_z, ivper,ivpar) + f_rzvv(i_r, i_z, ivper-1, ivpar)) / &
   !   &              (vPer_rzvv(i_vper+1) - vPer_rzvv(i_vper))**2
         end if

         ! BC at vpar = - bound on velocity grid 
         if(ivpar .eq. 1)then
             dfdvpar_rzvv(i_r, i_z, ivper, ivpar) = 0.0
         end if 

         ! BC at vpar = n_vpar, the largest vpar on the velocity grid. TODO confirm 2nd derivs not needed. 
         if(ivpar .eq. n_vpar)then
            dfdvpar_rzvv(i_r, i_z, ivper, ivpar) = 0.0
         end if

         end do
      end do
   end subroutine deriv_dfdvpar

   ! Example usage
   !       call  sgrate_cyl(uperp, upara, f, 1, nuper, 1, nupar, &
   !   &      dens(i_psi), NUPER, NUPAR)
   
      subroutine sgrate_cyl(x, y, f, nx1, nx2, ny1, ny2, ans, nxmax, nymax)

      implicit none

      integer i, j, nx1, nx2, ny1, ny2, nxmax, nymax
      real f(nxmax, nymax), x(nxmax), y(nymax), ans, dx, dy
      real :: pi

      data pi/3.141592654/

      ans = 0.0

      do i = nx1, nx2 - 1
         dx = x(i+1) - x(i)

         do j = ny1, ny2 - 1
            dy = y(j+1) - y(j)

            ans = ans + dx * dy * (x(i)   * f(i, j)  &
     &                           + x(i+1) * f(i+1, j) &
     &                           + x(i)   * f(i, j+1) &
     &                           + x(i+1) * f(i+1, j+1) ) / 4.0
         end do

      end do

      ans = ans * 2.0 * pi

      return
      end

end module read_4df_mod