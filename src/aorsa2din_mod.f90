      module aorsa2din_mod
      
      implicit none
      
      logical :: write_for_azzam = .false. ! write out matrix for Azzam
      logical :: use_new_z2                                               
      logical :: use_fourier_z2       

!     --------------------------------------------------------
!     Declarations and defaults for aorsa2d.in input variables
!     --------------------------------------------------------

      CHARACTER(128) :: eqdsk = 'g1080408021.txt'               ! eqdsk name
      CHARACTER(128) :: netCDF_file1 = 'phillips_nstx3.5.2.nc'  !cql3d distribution file name 1
      CHARACTER(128) :: netCDF_file2 = 'phillips_nstx3.5.2.nc'  !cql3d distribution file name 2
      
      CHARACTER(128) :: root =  'none'  ! if (root = 'none') sum all modes
                                       ! if (root = 'fast') sum only modes within 20% of fast wave 
                                       ! if (root = 'slow') sum only modes within 20% of slow wave                                                                          
      
      real :: xkx_ono = -30.          !kx (m-1) for which Ono's dispersion relations is plotted
      real :: xky_ono = 0.0           !ky (m-1) for which Ono's dispersion relations is plotted      
      
      integer :: nstrap = 4           ! number of current straps in antenna
      
      real :: scrape = 100.           ! scrape-off length in rho for iprofile = 5 (numerical profiles) 
      real :: nmin = 1.0e+18          ! minimum density allowed in scrape-off layer for iprofile = 5 (numerical profiles)
      real :: Temin = .015          ! minimum temperature (keV) allowed in scrape-off layer for iprofile = 5
      real :: rhomax = 1.0            ! maximum value of the rho grid

            
      real :: xlt = 20.0E-02          ! width of current strap (m)
      real :: wd = 50.0E-02           ! separation between centers of antenna elements
      
      real :: enorm_factor = 0.0      ! if (enorm_factor .eq. 0.0) AORSA & CQL3D use same enorm (default)
                                      ! if (enorm_factor .gt. 0.0) AORSA enorm = enorm_factor x the maximum energy 
                                       
      real :: xkperp_cutoff = 0.75    ! fraction of xkperp above which the electron conductivity (sig3) 
                                      !  is enhanced to short out noise in E_parallel (default = 0.75)
                                      
      real :: damping = 0.0           !  enhancement factor (default = 0.0) for the electron conductivity (sig3) 
                                      !  applied above the fractional value of xkperp (xkperp_cutoff) to 
                                      !  short out noise in E_parallel 
      
      real :: taue = 50.0e-03         ! energy diffusion time used in flow drive calculation
      real :: theta_ant = 0.0
      real :: eslowev = 3.5e+06
      real :: amu_slo = 4.0
      real :: z_slo = 2.0
      real :: eta_slo = 0.0
      real :: xnuomi = 0.0             !-----xnuomi is the ion collison rate used in hot and cold plasma dielectrics
      real :: xnuome = 0.0             !-----xnuome is the electron collison rate used in hot and cold plasma dielectrics
      real :: xnuead = 0.0000E+00      !-----xnuead = ad hoc collision frequency for electron in sec-1
      real :: xnu1ad = 0.0000E+00      !-----nu1ad=ad hoc collision frequency for majority ions in sec-1                                      
      real :: xnu2ad = 0.0000E+00      !-----nu2ad=ad hoc collision frequency for minority ions in sec-1
      real :: rant = 0.0               !-----rant = major radius of antenna in meters (default is 0.0 in which case, psiant = .95)
      real :: antang = 0.0             !-----antang = angle that the antenna straps make with the vertical 
                                       !-----in degrees (default is 0 deg. - vertical straps)

      real :: dthetant0 = 40.
      real :: dpsiant0 = .05
      real :: antlen = 1.0
      real :: antlc = .0001            !-----antlc = propagation constant along the antenna = c / vphase
      real :: psilim = .99
      real :: psimask = 1.05
      real :: psiant = 0.95
      real :: psimol = 1.00            !-----psimol: if (psimol .eq. 1.0) molifier is NOT used for plasma profiles (default)
                                       !-----        if (psimol .ne. 1.0) molifier is used for profiles, centered at psimol
      real :: psipne = 0.50
      real :: psipte = .30
      real :: psipti1 = .30
      real :: psipti2 = .30
      real :: psipti3 = .30
      real :: psipti4 = .30
      real :: psipti5 = .30
      real :: psipti6 = .30
      real :: yant = 0.0               !-----yant = location of antenna center in Z (m)
      real :: te0  = 4.2900E+03        !-----te0=central value of eletron temperature in eV
      real :: ti0  = 7.0700E+03        !-----ti0=central value of ion temperature in eV
      real :: ti02 = 7.0700E+03
      real :: ti03 = 7.0700E+03
      real :: ti04 = 7.0700E+03
      real :: ti05 = 7.0700E+03
      real :: ti06 = 7.0700E+03
      real :: epszet = 1.0000E-07      !-----epszet=error criterion for Z function calculation if disp is used
      real :: delta0 = 0.0000E+00      !-----delta0=numerical damping for Bernstein wave:  about 1.e-04 (dimensionless)
      real :: xwall = 0.0000E+00       !-----xwall = not used
      real :: xnwall = 0.0000E+00      !-----xnwall = density of metal put on last mesh point
      real :: amu1 = 2.0000E+00        !-----amu1= ratio of majority ion to hydrogen ion mass
      real :: amu2 = 1.0000E+00        !-----amu2= ratio of minority ion to hydrogen ion mass
      real :: z1 = 1.0000E+00          !-----z1=ratio of majority ion charge to hydrogen ion charge
      real :: z2 = 1.0000E+00          !-----z2=ratio of minority ion charge to hydrogen ion charge
      real :: eta = 0.0                !-----eta=ratio of minority ion density to electron density
      real :: b0 = 2.08                !-----b0=value of magnetic field at x=0 in Tesla
      real :: q0 = 1.0                 !-----q0=value of inverse rotational transform on axis
      real :: rt = 1.68                !-----rt= major radius of torus
      real :: ekappa = 2.0             !-----ekappa = elongation
      real :: rwleft  = 0.0            !-----rwleft = major radius of the left conducting wall
      real :: rwright = 0.0            !-----rwright = major radius of the right conducting wall
      real :: ytop    = 0.0            !-----ytop = location of the upper conducting wall
      real :: ybottom = 0.0            !-----ybottom = location of the lower conducting wall
      real :: ymax = 0.0               !-----ymax = radius in vertical (y) direction- in default it is set to awallx
      real :: freqcy = 3.2000E+07      !-----freqcy= rf frequency in Hertz
      real :: aplasm = 7.0000E-01      !-----aplasm = not used
      real :: alim = 100.0             !-----alim = location of limiter
      real :: grad = 0.0               !-----grad = 0.0 ignors gradients in Wdot (default)
                                       !-----grad = 1.0 includes gradients in Wdot      
      real :: signbz = 1.0000E+00
      real :: xn0 = 3.1100E+19         !-----xn0=electron density at x=0
      real :: xn2 = 0.0
      real :: xn3 = 0.0
      real :: xn4 = 0.0
      real :: xn5 = 0.0
      real :: xn6 = 0.0
      real :: xnslo = 0.0
      real :: flat = 0.0000E+00         !-----flat = not used
      real :: b1rat = 7.0000E-01        !-----b1rat = not used
      real :: b2rat = 1.3000E+00        !-----b2rat = not used
      real :: curdnx = 1.0              !-----curdnx=Amps/meter of toroidal length of antenna in the x direction
      real :: curdny = 1.0              !-----curdny=Amps/meter of toroidal length of antenna in the y direction
      real :: curdnz = 1.0              !-----curdnz=Amps/meter of toroidal length of antenna in the z direction
      real :: prfin = 0.0               !-----prfin = total applied RF power
      real :: xnuabs = 0.0000E+00       !-----xnuabs not used
      real :: xbnch = 0.0000E+00        !-----xbnch not used
      real :: xleft = -7.000E-01        !-----xleft=left boundary for energy integrals and outgoing energy flux
      real :: xright = 7.0000E-01       !-----xright=right boundary for energy integrals and incoming energy flux
      real :: qavg0 = 1.0               !-----qavg0 is the rotational transform on axis
      real :: xnlim  = 0.0              !-----xnlim=electron density in scrape-off region (x>aplasm)
      real :: xn2lim = 0.0
      real :: xn3lim = 0.0
      real :: xn4lim = 0.0
      real :: xn5lim = 0.0
      real :: xn6lim = 0.0
      real :: xnslolim = 0.0
      real :: telim  = 0.0000E+00       !-----telim = electron temperature in scrape-off region (x>aplasm)
      real :: tilim  = 0.0000E+00       !-----tilim = ion temperature in scrape-off region (x>aplasm)
      real :: ti2lim = 0.0000E+00
      real :: ti3lim = 0.0000E+00
      real :: ti4lim  = 0.0000E+00
      real :: ti5lim = 0.0000E+00
      real :: ti6lim = 0.0000E+00      
      real :: dfreq = 0.0000E+00        !-----dfreq not used   
      real :: dkz = 0.0000E+00          !-----dkz not used
      real :: xnudip = 2.5000E+00       !-----xnudip is not used
      real :: adip = 0.0000E+00         !-----adip is not used                          
      real :: efold = 0.0000E+00        !-----efold is not used
      real :: amu3 = 12.0               !-----amu3=ratio of third ion mass to hydrogen ion mass
      real :: z3 = 6.0                  !-----z3=ratio of third ion charge to hydrogen ion charge
      real :: eta3 = 0.0                !-----eta3=ratio of third ion density to electron density
      real :: xnu3ad = 0.0              !-----xnu3ad not used
      real :: amu4 = 12.0               !-----amu4=ratio of 4th ion mass to hydrogen ion mass
      real :: z4 = 6.0                  !-----z4=ratio of 4th ion charge to hydrogen ion charge
      real :: eta4 = 0.0                !-----eta4=ratio of 4th ion density to electron density
      real :: xnu4ad = 0.0              !-----xnu4ad not used
      real :: amu5 = 12.0               !-----amu5=ratio of 5th ion mass to hydrogen ion mass
      real :: z5 = 6.0                  !-----z5=ratio of 5th ion charge to hydrogen ion charge
      real :: eta5 = 0.0                !-----eta5=ratio of 5th ion density to electron density
      real :: xnu5ad = 0.0              !-----xnu5ad not used
      real :: amu6 = 12.0               !-----amu6=ratio of 6th ion mass to hydrogen ion mass
      real :: z6 = 6.0                  !-----z6=ratio of 6th ion charge to hydrogen ion charge
      real :: eta6 = 0.0                !-----eta6=ratio of 6th ion density to electron density
      real :: xnu6ad = 0.0              !-----xnu6ad not used
      real :: xdelta = 5.5000E-01       !-----xdelta not used
      real :: wdelta = 0.0000E+00       !-----wdelta not used
      real :: xdelt2 = -7.000E-02       !-----xdelt2 not used
      real :: wdelt2 = 0.0000E+00       !-----wdelt2 not used
      real :: zeffcd = 2.5000E+00       !-----zeffcd=Zeff for Ehst-Karney current drive calculation
      real :: rzoom1 = 0.0              !-----rzoom1 = R on left side of wdot calculation
      real :: rzoom2 = 0.0              !-----rzoom2 = R on right side of wdot calculation
      real :: yzoom1 = 0.0              !-----yzoom1 = Y on bottom of wdot calculation
      real :: yzoom2 = 0.0              !-----yzoom2 = Y on top of wdot calculation
      
      real :: ydisp = 0.0               !-----ydisp = y value (m) at which contours of D(x, kx) are drawn
      real :: xkydisp = 0.0             !-----xkydisp = ky value (m-1) for which contours of D(x, kx) are drawn
      
      real :: dfquotient = 0.0        !-----dfquotient = divisor of df for which a single contour of D(x, kx) is drawn
                                        !     if (dfquotient .eq.  0.0) all 20 contours are drawn (default)
                                        !     if (dfquotient .ne.  0.0) only one coutour is drawn with df/dfquotient
                                        
      real :: rmin_zoom = 2.15        !-----R_min for edge density plot
      real :: rmax_zoom = 2.35        !-----R_max for edge density plot
      
      real :: yplot = .05             !-----yplot = cut in y for slow wave dipersion plots
      
      real :: alphan = 1.0
      real :: alphan2 = 1.0
      real :: alphan3 = 1.0
      real :: alphan4 = 1.0
      real :: alphan5 = 1.0
      real :: alphan6 = 1.0
      real :: alphan_slo = 1.0
      real :: alphate = 1.0
      real :: alphati = 1.0
      real :: alphati2 = 1.0
      real :: alphati3 = 1.0
      real :: alphati4 = 1.0
      real :: alphati5 = 1.0
      real :: alphati6 = 1.0
      real :: betan = 2.0
      real :: betan2 = 2.0
      real :: betan3 = 2.0
      real :: betan4 = 2.0
      real :: betan5 = 2.0
      real :: betan6 = 2.0
      real :: betan_slo = 2.0
      real :: betate = 2.0
      real :: betati = 2.0
      real :: betati2 = 2.0
      real :: betati3 = 2.0
      real :: betati4 = 2.0
      real :: betati5 = 2.0
      real :: betati6 = 2.0

           
      
      integer :: version_number = 24      !-----version_number: version of AORSA         
      integer :: n_prof_flux = 0          !-----n_prof_flux = flag determining whether profiles are wrt toroidal or poloidal flux
                                          !           if(n_prof_flux .eq. 0) profiles are wrt sqrt(poloidal) flux (default)
                                          !           if(n_prof_flux .ne. 0) profiles are wrt sqrt(toroidal) flux 
                                       
      integer :: upshift = 1              !-----upshift: if (upshift .ne.  0) upshift is turned on (default)
                                          !-----if (upshift .eq. -1) upshift is turned off for xkperp > xkperp_cutoff
                                          !-----if (upshift .eq. -2) don't allow k_parallel too close to zero
                                          !-----if (upshift .eq.  0) upshift is turned off always
                                          
      integer :: odd_order = 0            !-----odd_order: 
                                          !     if (odd_order .eq.  0) odd_order derivatives are NOT included (default)
                                          !     if (odd_order .ne.  0) odd_order derivatives are included (Maxwellians ONLY) 
                                          
      integer :: nsmooth = 4              !-----nsmooth: 
                                          !     if (nsmooth .eq. 0) derivatives are NOT smoothed
                                          !     if (nsmooth .ne. 0) derivatives are smoothed nsmooth times (default: 4 times)                                     
                                          
      integer :: eqdsk_special = 0        !     if (eqdsk_special .eq.  0) regular geqdsk (default)
                                          !     if (eqdsk_special .ne.  0) special eqdsk for Alcator-C lower hybrid                                                                               
                                          
      logical :: nphi_sum_only = .false.  !-----if (nphi_sum_only .eq. .true.) skip aorsa and just sum the modes for nphi_number > 1                                 
                                          !-----if (nphi_sum_only .eq. .false.) run aorsa and sum modes if nphi_number > 1 (default)

      logical :: anti_alias = .false.     !-----if (anti_alias .eq. .false.) no filter:  keep all modes (default)                                 
                                          !-----if (anti_alias .eq. .true. ) anti-aliasing filter (two-thirds rule)
                                          
      logical :: use_new_wdot = .false.   !-----if (use_new_wdot .eq. .false.) use original wdote - resonant terms only (default)                                 
                                          !-----if (use_new_wdot .eq. .true. ) use new wdote - both resonant and non-resonant terms                                      
                                         
      logical :: use_no_damp = .false.    !-----if (use_no_damp .eq. .true.) turn off artificial damping for power calculation                                  
                                          !-----if (use_no_damp .eq. .false.) keep artificial damping for power calculation calculation (default)                                         
                                                                                          
      logical :: read_solution = .false.  !-----if (read_solution .eq. .true.) skip scalapack solution and read fields from out_fields                                 
                                          !-----if (read_solution .eq. .false.) solve for field with scalapack (default)                                                                                   
                                          

      logical :: plot_only = .false.      !-----if (plot_only .eq. .true.) skip aorsa and just do plots                                 
                                          !-----if (plot_only .eq. .false.) run aorsa and do plots (default)
                                          
      logical :: curved_ant = .true.      !-----if (curved_ant .eq. .true.) antenna follows flux surface (default)                                 
                                          !-----if (curved_ant .eq. .false.) antenna is vertical
                                          
      integer :: norm = 1                 !-----if (norm .eq. 1) normalize to jdote (default)   
                                          !-----if (norm .eq. 2) normalize to wdot 
                                          !-----if (norm .eq. 3) normalize to wdot_ql                                                                                      
                                          

      integer :: z2_electron = 0          !-----if (z2_electron .eq. 0) use the original Z2 function for electrons (default)    
                                          !-----if (z2_electron .eq. 1) use the Z2 table for electrons with l = 0 (Taylor expansion along field line)
                                          !-----if (z2_electron .eq. 2) use Fourier expansion along field line for electrons with l = 0 (full orbits)
                                                                                                                                                                                                                                                                                                                   
                                            
      integer :: i_write = 0              !-----i_write: if (i_write .eq. 0) 4-D ORBIT_RF file is NOT written (default)
                                          !-----         if (i_write .ne. 0) 4-D ORBIT_RF file IS written

      integer :: n_bin = 2
      integer :: iql = 1
      integer :: i_antenna = 1            ! i_antenna = flag determining which antenna model is used
                                          ! if(i_antenna .eq. 0) antenna current is Gaussian 
                                          ! if(i_antenna .eq. 1) antenna current is cos(ky * y)  (default)
                                          ! where ky = omgrf / vphase = (omgrf / clight) * antlc = k0 * antlc
                                          ! For constant current, set antlc = 0.0
      integer :: nuper = 65
      integer :: nupar = 129
      integer :: nkperp = 201             !-----nkperp: number of kperp values used in Lee's interpolation version of the 
                                          !     non-Maxwellian sigma (default = 201: interpolates on 201 points)
                                          !     if (nkperp .eq. 0) there is no interpolation
                                     
      integer :: nzeta_wdote = 51         !-----nzeta_wdote: if (nzeta_wdote .eq. 0) no wdot calculation for electrons
                                          !-----             if (nzeta_wdote .eq. 1) wdote is calculated without interpolation
                                          !-----             if (nzeta_wdote .ge. 2) wdote is calculated with interpolation
                                          !-----                 over nzeta_wdote grid points (default is 51) 
                                             
      integer :: nzeta_wdoti = 51         !-----nzeta_wdoti: if (nzeta_wdoti .eq. 0) no wdot calculation for ions
                                          !-----             if (nzeta_wdoti .eq. 1) wdoti is calculated without interpolation
                                          !-----             if (nzeta_wdoti .ge. 2) wdoti is calculated with interpolation
                                          !-----                 over nzeta_wdoti grid points (default is 51)                                             
                                          
      integer :: i_sav = 0
      integer :: j_sav = 0      
      integer :: i_sav1 = 0
      integer :: j_sav1 = 0      
      integer :: i_sav2 = 0
      integer :: j_sav2 = 0      
      integer :: isolve = 1               !-----isolve = flag determining whether wave solution is calculated or not
                                          !     if(isolve .ne. 0) wave solution is calculated (default)
                                          !     if(isolve .eq. 0) wave solutions is read from file='fields_fourier'

      integer :: ftrap = 1                !-----ftrap = integer flag determining whether trapped particles effect current drive
                                          !        if(ftrap.eq.0) no trapped particles
                                          !        if(ftrap.ne.0) include trapped particles (default)
                                     
      integer :: nnode_local = 0          !-----nnode_local = number of local Fourier modes used to calculate local Wdot
                                          !-----      if (nnode_local .le. 0)Wdot is NOT calculated (default)
                                          !-----      if (nnode_local .gt. 0)Wdot is calculated

      integer :: nnode_overlap = 0        !-----nnode_overlap is the number of overlapped points
      integer :: iprofile = 3             !-----iprofile:  if (iprofile .eq. 1) generic profiles (Gaussian)
                                          !-----           if (iprofile .eq. 2) generic profiles (parabolas)
                                          !-----           if (iprofile .eq. 3) fits of form (1 - rho**beta)**alpha (default)
                                          !-----           if (iprofile .eq. 5) numerical profiles from namelist
      integer :: nboundary = 1            !-----nboundary: if(nboundary .eq. 1)flux surface boundary (default)
                                          !-----           if(nboundary .eq. 0)square boundary
      integer :: nprow = 8
      integer :: npcol = 8

      integer :: nwdot = 0                !-----nwdot=number of radial modes used in wdot and flow (fy) calculation

      integer :: lmax = 5                 !-----lmax = highest order Bessel function kept in ion plasma conductivity wdot
      integer :: lmaxe = 1                !-----lmaxe = highest order Bessel function kept in electron plasma conductivity and wdot     
      
      
      integer :: ibessel = 1              !-----ibessel = flag determining which Bessel functions are used in Wdot
                                          !        if(ibessel.eq.0) Complex Bessel functions are used from besic
                                          !        if(ibessel.eq.1) Real Bessel functions used from ribesl (default)
                                          !        if(ibessel.eq.2) 2nd order expanded Bessel functions are used in Wdot only
      integer :: inu = 0                  !-----inu:   if(inu.eq.0)real collisions are left out
      
      integer :: iprint = 1               !-----iprint:  if (iprint .eq. 28) print fields_local
                                          !-----         if (iprint .eq. 1) don't print fields_local (default)
                                     
      integer :: iexact = 1               !-----iexact:  not used
      integer :: iroot = 2                !-----iroot: not used
      integer :: iequat = 1
      integer :: igeom = 5                !-----igeom: if(igeom.eq.2)Solovev flux surfaces
                                          !-----       if(igeom.eq.5)EQDSK (GA) flux surfaces (default)
      integer :: iqx = 4                  !-----iqx:  not used
      
      integer :: iqprof = 1               !-----iqprof: if(iqprof.eq.1) q proportional to density (default)
                                          !-----        if(iqprof.eq.2) q proportional to sqrt(density) for TAE runs
      integer :: iez = 0                  !-----iez:  not used
      integer :: nstep  = 16              !-----nstep: not used
      integer :: nabs = 2                 !-----nabs not used
      integer :: isigma = 1               !-----if(isigma.eq.0) cold plasma conductivity is used.
                                          !-----if(isigma.eq.1) hot  plasma conductivity is used (default)
      integer :: nzfun = 1                !-----nzfun:  if(nzfun.eq.0) Simple Z function is used from ZFUN 
                                          !-----        if(nzfun.eq.1) Generalized Z function of Brambilla is used (default)
                                          !-----        if(nzfun.eq.2) Z function of Smithe is used by doing numerical integrals.
                                          !-----        if(nzfun.eq.3) Z function table lookup of Smithe is used                                  
                                          
                                     
      integer :: iabsorb = 2              !-----iabsorb not used
      integer :: itemp = 0                !-----itemp not used
      integer :: nfreqm = 1               !-----nfreqm not used
      integer :: nkzm = 1                 !-----nkzm not used
      integer :: idens = 0                !-----idens not used
      integer :: ibackground = 1          !-----ibackground is not used
      integer :: idiag = 5
      integer :: jdiag = 4
      integer :: ndiste = 0               !-----ndist:  if (ndist .eq. 0) Maxwellian is used in sigmad_stix
                                          !-----        if (ndist .eq. 1) non-Maxwellian is used in sigmad_stix
      integer :: ndisti1 = 0
      integer :: ndisti2 = 0
      integer :: ndisti3 = 0
      integer :: ndisti4 = 0
      integer :: ndisti5 = 0
      integer :: ndisti6 = 0      

      integer :: nphimx, nphasemx 
      integer :: nmodesx = 128            !-----nmodesx=number of modes used in the x direction
      integer :: nmodesy = 128            !-----nmodesy=number of modes used in the y direction
      integer :: nmodeslb = 256           !-----nmodeslb =number of modes used along the field line
      
      integer :: izfunc                   !-----not used 

      integer :: nphi_number = 10         !-----number of nphi's in the nphi_array 
           
      parameter (nphimx = 200)            !-----dimension of nphi_array 
      integer :: nphi_array(nphimx) = 0   !-----list of nphi's in the nphi_array
      
      parameter (nphasemx = 12) 
      real :: phase_array(nphasemx) = 0.0   !-----list of phase on each antenna strap     
      
      real phase, zmin, zmax, phi0, amplt(20) 
                   
      namelist/aorsa2din/nmodesx, nmodesy, nmodeslb, nwdot, lmax, lmaxe, ibessel, &
     &    ti0, xnuead, xnu1ad, xnu2ad, rant, te0, yant,  &
     &    ti02, ti03, ti2lim, ti3lim, nuper, nupar, &
     &    ti04, ti05, ti06, ti4lim, ti5lim, ti6lim,  &
     &    inu, iprint, iexact, delta0, xwall, xnwall,  &
     &    iroot, iequat, igeom, epszet,  &
     &    dthetant0, dpsiant0, psilim, psiant, psimol, psipne, psipte,  &
     &    psipti1, psipti2, psipti3,  &
     &    psipti4, psipti5, psipti6,  &
     &    iqx, iqprof, iez, nprow, npcol, &
     &    amu1, amu2, z1, z2, eta, izfunc,  &
     &    b0, rt, ytop, ybottom, freqcy, aplasm,  &
     &    xnlim, xn2lim, xn3lim, xnslolim, signbz,  &
     &    xn4lim, xn5lim, xn6lim,  &
     &    xn0, xn2, xn3, xnslo, flat, b1rat, b2rat, curdnx, curdny,  &
     &    xn4, xn5, xn6, curdnz,  &
     &    nstep, nabs, xnuabs, xbnch, xleft, xright,  &
     &    isigma, itemp, telim, tilim,  &
     &    nfreqm, dfreq, nkzm, dkz,  &
     &    idens, xnudip, adip, efold,  &
     &    amu3, z3, eta3, xnu3ad,  &
     &    amu4, z4, eta4, xnu4ad,  &
     &    amu5, z5, eta5, xnu5ad,  &
     &    amu6, z6, eta6, xnu6ad,  &
     &    xdelta, wdelta, xdelt2, wdelt2, zeffcd, &
     &    rzoom1, rzoom2, yzoom1, yzoom2, ibackground, iabsorb, q0,  &
     &    prfin, nzfun, alim, grad, qavg0, ymax,  &
     &    alphan,  alphan2, alphan3, alphan_slo,  &
     &    alphan4, alphan5, alphan6,  &
     &    alphate,  alphati, alphati2, alphati3,  &
     &    alphati4, alphati5, alphati6,  &
     &    ekappa, rwleft, rwright, xnuome, xnuomi,&
     &    nboundary, eta_slo, amu_slo, z_slo, eslowev, nnode_local,  &
     &    nnode_overlap, iprofile, ftrap, isolve, &
     &    betan, betan2, betan3, betan_slo, betate, betati, betati2,  &
     &    betan4, betan5, betan6, betati4, betati5, betati6,  &
     &    betati3, taue, theta_ant, Temin,  &
     &    ndiste, ndisti1, ndisti2, ndisti3,  &
     &    ndisti4, ndisti5, ndisti6, nkperp, nzeta_wdote, n_bin, antlen,  &
     &    eqdsk, iql, i_antenna, antlc, n_prof_flux, netcdf_file1,  &
     &    netcdf_file2, upshift, xkperp_cutoff, damping, i_write,  &
     &    nstrap, xlt, wd, nphi_number, nphi_array, phase_array,  &
     &    enorm_factor, nphi_sum_only, odd_order, read_solution, &
     &    eqdsk_special, nsmooth, plot_only, ydisp, xkydisp, dfquotient, &
     &    anti_alias, scrape, nmin, rmin_zoom, rmax_zoom, root, rhomax, &
     &    xkx_ono, xky_ono, curved_ant, yplot, z2_electron, use_new_wdot,  &
     &    use_no_damp, nzeta_wdoti, norm, antang, write_for_azzam               

      end module aorsa2din_mod

