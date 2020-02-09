      program process_cql3d_output

!- -------------------------------------------------------------------
! " processing cql3d output
!Lee A. Berry 12/6/06  based on D. McMcune's test program
!BobH, 03/28/07        based on Berry  process_aorsa_output
!--------------------------------------------------------------------
!
 
! State elements are traditional f77 integer, floating point, and character
! string scalars and arrays all contained within a large container data 
! type "plasma_state".  Two instances of this container data type are
! declared in the module "plasma_state_mod":
!
!   ps -- the current state (timestep now being computed).
!   psp -- the prior or "committed" state (completed, prior timestep)
!
! Elements of the state can be referenced directly using the standard
! f95 syntax-- e.g. ps%nrho is the size of one of the radial flux 
! coordinate "rho" grids in the state.
!
! State elements can be directly modified by codes that use the plasma
! state-- although this should be done carefully.  Generally, items in 
! the state are meant to be shared between multiple components of the 
! IPS; conventions for use of the state data will need to evolve.
!
! States will be mapped to NetCDF files.  The module "plasma_state_mod"
! defines two variables for these filenames:
!
!   CHARACTER*256 :: state_file = 'cur_state.cdf'
!   CHARACTER*256 :: prior_file = 'prev_state.cdf'
!
! The assigned default values can of course be modified.
!-----------------------------------
! In addition to direct access to state data elements, the following
! subroutines are available (defined in the f95 module plasma_state_mod)
! (this is the module's public interface):
!
!    integer :: ierr -- status code returned by many routines (0=OK).
!

!
!    SUBROUTINE ps_clear_profs(ierr)
!	Set all profile data to zero-- this might be desirable when starting
!	to build a state at a new time.  It will be easier to tell what has
!	been added, and what not, if quantities not added are zero, rather
!	than from the prior timestep.  All the prior timestep data is still
!	accessible in the prior state object psp-- i.e. psp%rho_eq(...), etc.
!	Scalar data and grids are not affected by this call-- just profiles.
!
!    SUBROUTINE ps_update_equilibrium(<g-filename>,ierr)
!	Update state MHD equilibrium from G-eqdsk file <g-filename>
!	   (could also be an MDSplus G-eqdsk timeslice from an experiment).
!	Compute state quantities derived from the equilibrium
!	Arrays such as enclosed volumes (m^3) ps%vol(1:ps%nrho_eq) are 
!	filled in.
!
!    SUBROUTINE ps_store_plasma_state(ierr)
!	Update interpolation information and store the state (ps) in
!	the file (state_file).
!
!    SUBROUTINE ps_update_plasma_state(ierr)
!	Update interpolation information but do not write a file.
!
!    SUBROUTINE ps_commit_plasma_state(ierr)
!	Copy current state (ps) to prior state (psp).  Save prior state
!	in file (prior_state).  The file (state_file) is not modified--
!	use ps_store_plasma_state for this.
!
!    SUBROUTINE ps_get_plasma_state(ierr)
!	Read the current state (ps) with all interpolation information
!	from (state_file) --AND-- read the prior state (psp) with all
!	its interpolation information from the file (prior_state).
!
! Profile IDs:  each state array that is defined over one or more coordinate
!   grids is assigned an ID variable or array with name ID_<name>.  These IDs
!   are needed to refer to specific profiles for interpolation or rezoning.
!
!   Examples mapping from swim_state_spec.dat to plasma state object "ps":
!
!      R|pclin  chi_e(nrho)	     ! electron thermal conductivity
!
!      -> ps%chi_e(...)  (1d allocated REAL(KIND=rspec) array (1:ps%nrho)) &
!      -> ps%id_chi_e	 INTEGER scalar
!
!      R|units=m^-3|step  ns(~nrho,0:nspec_th)       ! thermal specie density
!
!      -> ps%ns(...)	 (2d allocated REAL(KIND=rspec) array
!      -> ps%id_ns(...)  (1d allocated INTEGER array)
!
!	       ps%ns(1:ps%nrho-1,0:ps%nspec_th)
!	       ps%id_ns(0:ps%nspec_th)
!
! Direct interpolation:
!
!    SUBROUTINE PS_INTRP_1D(...)  for 1D profiles  LAB based on steps
!
!	CALL PS_INTRP_1D( &
!	     x,  &	! target of interpolation: scalar or 1d vector
!	     id, &	! profile(s) interpolated: scalar 1d or 2d INT array
!	     ans,&	! result, dimensioned to match x(...) and id(...)
!	     ierr,   &  ! completion code 0=OK
!	     icur,   &  ! OPTIONAL: "ps_previous" or "ps_current" (default)
!	     ideriv, &  ! OPTIONAL: derivative control for all IDs
!	     ideriv1s,& ! OPTIONAL: derivative control for each ID separately
!	     iccw_th)	! OPTIONAL: .FALSE. for clockwise poloidal angle
!			  (theta) coordinate
!
!	If x is a vector, the 1st dimension size of ans matches size(x);
!	subsequent dimensions match sizes of dimensions of id(...).
!
!	icur can be used to select the current or prior state; current state
!	is the default.
!
!	ideriv1s is only available if id(...) is an array.  If so, ideriv1s
!	takes precedence over ideriv, if both are present.  The dimensioning
!	of INT array ideriv1s must match dimensioning of ID exactly.
!
!	If neither ideriv nor ideriv1s are specified, the interpolating 
!	function value is evaluated with no derivatives.
!
!	iccw_th would only be needed if a profile f(theta) is ever defined.
!
!    SUBROUTINE PS_INTRP_2D(...)  for 2D profiles
!
!	(interface is like PS_INTRP_1D, except that:
!
!	   x -> x1,x2 -- two interpolation target scalars or vectors must
!			 be supplied; the coordinate to which they belong
!			 will match the declaration
!
!	   similarly, optional derivative control is available separately
!	   for each coordinate:
!	     ideriv -> ideriv1,ideriv2
!	     ideriv1s -> ideriv1s,ideriv2s
!
! Profile rezoning integration:
!
!    SUBROUTINE PS_RHO_REZONE(...) for "conservative" rezoning
!      of profiles f(rho) 
!      (rezoning to 1d of profiles f(rho,theta) will be done but is
!      not yet implemented -- DMC 16 Oct 2006).
!
!    SUBROUTINE PS_RHOTH_REZONE(...) for "conservative" rezoning
!      of profiels f(rho,theta) -- not yet implemented DMC 16 Oct 2006
! 
!-----------------------------------
! Implementation of interpolation and rezoning-- the "visible" state
! elements ps%* and psp%* are not the entire state.  When interpolation
! and rezoning operations are carried out, additional hidden information
! is accessed which defines the profile interpolation methods.
!
! Interpolation methods for profiles are part of the state specification
! "swim_state_spec.dat" and are handled by generated code.
!-----------------------------------

! define state object and interface...

      USE plasma_state_mod
!--------------------------------------------------------------------------

      use swim_global_data_mod, only :
     1 rspec, ispec,                  ! int: kind specification for real and integer
!!            & swim_string_length, & ! length of strings for names, files, etc.
     1 swim_error                     ! error routine
    

    
!--------------------------------------------------------------------------
!
!   Data declarations
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------------
!
!   AORSA data that will be given to the state via the swim_out file
!   the aorsa output data is in a file called swim_out
!
!--------------------------------------------------------------------------
   

!BH070328:  Begin by reading in toroidal current density from cql3d
!           output netcdf file (mnemonic.nc) and putting it into the state

 
!      IMPLICIT NONE
      implicit none

      include 'netcdf.inc'
    
      integer, parameter :: swim_string_length = 256 !for compatibility LAB
      integer :: nnoderho, r0dim, vid_
      integer :: iarg

!     Storage for netcdf file elements and retrieval
      real*8, allocatable, dimension(:)   :: rya   ! radial grid--bin centers
      ! note that mesh in not uniform
      real*8, allocatable, dimension(:)   :: dvol   ! bin volumes   
      real*8, allocatable, dimension(:,:) :: wperp   !perp energy/particle
                                                      !tdim, rdim
      real*8, allocatable, dimension(:,:) :: wpar    !par energy/particle
                                                      !tdim, rdim
      real*8, allocatable, dimension(:,:,:) :: density !density of all species
                                                      !tdim, r0dim, species dim
      real*8, allocatable, dimension(:,:,:,:) :: powers !collection of power flows
      ! rdim, power_dim, general_species dim, rdim
      ! do ncdump on an output file to get documentation

      real*8, allocatable, dimension(:) ::   temp_1, temp_2, temp_3 !work arrays
      integer :: ncid,vid,istatus
      integer :: start(2),count(2)
      integer :: start_4(4), count_4(4)
      integer :: nt_id,lrz_id   !Time step dim id, radial fp bins id
      integer :: ngen_id, ntotal_id        !number of general species id
      integer :: nt,lrz, ngen, ntotal     !values for number of time bins.radial bins
                                         !general species , total number of species
      integer :: r0dim_id, rdim
      character*256 ::  cur_state_file, cql_out_file



     


!------------------------------------
!  local
      INTEGER :: ierr
      INTEGER :: iout = 6

 


! this call retrives the plasma state at the present ps%, and
! prvious time steps psp%, plus sets storage for the variables
! get plasma state from command line


      call get_arg_count(iarg)
      if(iarg .ne.2) then 
         print*, 'two command line arguments needed'
         print*, 'usage:  process_fp_rfmin_clq3d_output:   cur_state_file, cql_out_file'
 	 stop 1
      end if
      call get_arg(1, cur_state_file)
      write(*,*)  'current state file = ', trim(cur_state_file)
      call get_arg(2, cql_out_file)
      write(*,*) 'cql output file = ', trim(cql_out_file)


      CALL ps_get_plasma_state(ierr, trim(cur_state_file))

      if(ierr.ne.0) then
         write(iout,*) ' process cql_output: ps_get_plasma_state: ierr=',ierr
         stop 1
      end if
      
   

!.......................................................................
!     Open cql3d netcdf file
!.......................................................................

      istatus = nf_open(trim(cql_out_file),nf_nowrite,ncid)
      write(*,*)'after nf_open ncid=',ncid,'istatus',istatus

!     read in dimension IDs
!     we need the time, radial, and species dimensions
      istatus = nf_inq_dimid(ncid,'rdim',lrz_id)
      write(*,*)'proc_cql3d_op: after ncdid lrz_id',lrz_id,'istatus',istatus
      
      istatus = nf_inq_dimid(ncid,'tdim',nt_id)
      write(*,*)'proc_cql3d_op: after ncdid nt_id',nt_id,'istatus',istatus
      
      istatus = nf_inq_dimid(ncid,'gen_species_dim',ngen_id)
      write(*,*)'proc_cql3d_op:after ncdid ngen_id',ngen_id,'istatus',istatus
      
      istatus = nf_inq_dimid(ncid,'species_dim',ntotal_id)
      write(*,*)'proc_cql3d_op:after ncdid ntotal_id',ntotal_id,'istatus',istatus
      
      istatus =  nf_inq_dimid(ncid,'r0dim',r0dim_id)
      write(*,*)'proc_cql3d_op:after ncdid r0dim_id',r0dim_id,'istatus',istatus
      
      
      ! first the radial dimension
      
      istatus = nf_inq_dimlen(ncid, lrz_id, lrz)
      !call ncdinq(ncid,lrz_id,'rdim',lrz,istatus)
      write(*,*)'proc_cql3d_op: after ncdinq, # of rad bins= ',lrz,
     1     '  istatus=',istatus
           
      
      ! and the second radial dimension
      istatus = nf_inq_dimlen(ncid, r0dim_id, r0dim)
      !call ncdinq(ncid,r0dim_id,'r0dim',r0dim,istatus)
      write(*,*)'proc_cql3d_op: after ncdinq, # of rad bins= ',lrz,
     1     '  istatus=',istatus
      ! check if the are the same, if not more code is needed--grid
      ! subset us used for  solutions
      
      if(lrz .ne. r0dim) then
         print*, 'grid subset used for solutions'
         print*, 'another day, another time'
      stop 1
       end if 

     
      ! second the time dimension (will always use the last)
      istatus = nf_inq_dimlen(ncid, nt_id, nt)
      !call ncdinq(ncid, nt_id,'tdim', nt, istatus)
      write(*,*)'proc_cql3d_op: after ncdinq, # of t steps = ',nt,
     1     '  istatus=',istatus
     


      
      ! third the number of general species (will be one for now)

      !call ncdinq(ncid, ngen_id,'gen_species_dim', ngen, istatus)
      istatus = nf_inq_dimlen(ncid, ngen_id, ngen)
      write(*,*)'proc_cql3d_op: after ncdinq, #  gen specs = ', ngen,
     1     '  istatus=',istatus
!.......................................................................
!     inquire about dimension sizes:time steps,grid size,#species---
!.......................................................................




      ! fourth the number of total species (will be 1 (general) + 1 
      ! duplicate max
      ! + thermal + electron is the CMOD order for now)
      istatus = nf_inq_dimlen(ncid, ntotal_id, ntotal)
      !call ncdinq(ncid, ntotal_id,'species_dim',ntotal,istatus)
      write(*,*)'proc_cql3d_op: after ncdinq, #  species total = ',ngen,
     1     '  istatus=',istatus
     

     

!     allocate space for arrays to be read
      allocate (rya(lrz))  !the radial grid
      allocate (dvol(lrz))
      allocate (powers(lrz, 13, 1, nt))
      allocate (temp_1(lrz+1), temp_2(lrz+1), temp_3(lrz+1))
      allocate (wperp(lrz,nt),wpar(lrz,nt))
      
      
      


!     Specify reading ranges
      start(1:2)=1
      count(1)=lrz
      count(2)=nt
      
      
      start_4 = 1
      count_4(1) = lrz
      count_4(2) = 13
      count_4(3) = 1
      count_4(4) = nt

!.......................................................................
!     read netcdf variables
!.......................................................................
!     will need grid, volumes, powers, wpar, wperp, density


      ! the grid
      istatus = nf_inq_varid(ncid, 'rya', vid)
      istatus = nf_get_var_double(ncid, vid, rya)
      write(*,*)'proc_cql3d_op: after ncvgt, rya = ',rya
      
      ! the bin volumes
      istatus = nf_inq_varid(ncid, 'dvol', vid)
      istatus = nf_get_var_double(ncid, vid, dvol)
      write(*,*)'proc_cql3d_op: after ncvgt, dvol = ',dvol
      
      ! the powers 
      istatus = nf_inq_varid(ncid, 'powers', vid)            
      istatus = nf_get_var_double(ncid, vid, powers)
      write(*,*)'proc_cql3d_op: after ncvgt, powers(e,i) = ',
     1    powers(:, 1:2,  1, nt)
      
      
      ! the energies wperp--wpar
      print*, 'shape of wperp ', shape(wperp)  
      istatus = nf_inq_varid(ncid, 'wperp', vid)
      istatus = nf_get_var_double(ncid, vid, wperp)
      print*, 'shape of wperp ', shape(wperp)            
      write(*,*)'proc_cql3d_op: after ncvgt, wperp = ',
     1    wperp(:, nt)






      istatus = nf_inq_varid(ncid, 'wpar', vid)
      istatus = nf_get_var_double(ncid, vid, wpar)           
      write(*,*)'proc_cql3d_op: after ncvgt, wpar = ',
     1    wpar(:, nt)
      
!     now need to put into state
!     the perp, par, power goes into the same grid as they 
!     were created on--no interpolation needed--need multiply 
!     by dvol no conversion needed--watts=watts, keV = keV
      
      ! double eperp_mini(dim_nspec_rfmin, dm1_nrho_icrf) 
      !double epll_mini(dim_nspec_rfmin, dm1_nrho_icrf)
      !double pmini(dm1_nrho_icrf)
      !double pmine(dm1_nrho_icrf)
      
      
      ! this will work with state grids 
      !ps%pmine = dvol * powers(:,1,1,nt)
      !ps%pmini = dvol * powers(:,2,1,nt)
      !ps%eperp_mini(1,:) = wperp
      !ps%epll_mini(1,:) = wpar
      
      ! this will work for the test problem where the
      ! grids are not compatible 
      ! powers(*,1,k,t)=due to collisions with Maxw electrons
      ! k is the number of the general species
      ! nt is the last time slice of the cql3d run
      ps%pmine = -dvol * powers(1:ps%nrho_icrf-1,1,1,nt)
      ! powers(*,2,k,t)=due to collisions with Maxw ions
      ps%pmini = -dvol * powers(1:ps%nrho_icrf-1,2,1,nt)
      
      ps%eperp_mini(1:ps%nrho_icrf-1,1) = wperp(1:ps%nrho_icrf-1,nt)
      ps%epll_mini(1:ps%nrho_icrf-1,1) = wpar(1:ps%nrho_icrf-1,nt)

      call ncclos(ncid,istatus)

      call ps_store_plasma_state(ierr, trim(cur_state_file))
      if(ierr .ne. 0) then
         write(iout,*) 'plasma state not stored in process_cql3d_output'
         stop 1
      end if
      print*, 'power check on cql'
      print*, 'minority to electron power = ', sum(ps%pmine)
      print*, 'minority to ion power = ', sum(ps%pmini)

      print*, 'end of process cql3d output'
      
      end program process_cql3d_output
