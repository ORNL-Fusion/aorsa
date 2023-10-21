
      subroutine cql3d_setup1(netcdf_file, nuper, nupar, xm, enorm_factor, vc_mks_cql3d)
      use size_mod   

      USE read_CQL3D
      USE quick_vector_write_m
      USE f_expanded_m, only : C_matrix_file, in_file_CQL, &
     &  n_u_coeffs, n_theta_coeffs, u_norm, u0_Max_poly, &
     &  C_matrix, dens_Maxwell_psi, T_Maxwell_psi,              &
     &  read_C_matrix_m, eval_delta_f_exp_m, eval_f_exp_m
        
      USE basis_functions_m, only : u_fns_name
        

      IMPLICIT NONE
      
      INTEGER :: j, k, n, n_expand, ncount, number_points
      integer i, i_psi, i_theta, i_u, iflag, nproc, myid
      integer n_u_basis, n_theta_basis, ni, mi,kspeci
      real enorm_factor

      integer, parameter :: n_u_dim = 300
      integer, parameter :: n_theta_dim = 200
      integer, parameter :: n_psi_dim = 64

!     ----------------
!     450 x 450 modes:
!     ----------------
!      integer, parameter :: nmodesmax = 450
!      integer, parameter :: mmodesmax = 450

      integer, parameter :: nxmx = nmodesmax
      integer, parameter :: nymx = mmodesmax
      integer, parameter :: nrhomax = nmodesmax * 2

      integer, parameter :: nkdim1 = - nmodesmax / 2
      integer, parameter :: nkdim2 =   nmodesmax / 2

      integer, parameter :: mkdim1 = - mmodesmax / 2
      integer, parameter :: mkdim2 =   mmodesmax / 2

      real :: UminPara, UmaxPara, UmaxPerp, dupara, duperp, drho, drho10
      real fb, dfb_du, dfb_dtheta, xn_alpha, xne, dfdth
      real betat, alphat, shapen, shapet, xkt0, xnlim, xktlim, betan, alphan
      real wperp_max, wpar_max, enorm_aorsa_kev
      real enorm_cql3d_joules, enorm_cql3d_kev
      real enorm_aorsa_joules, vc_mks_aorsa, vc_ratio, vc_mks_cql3d  

      integer :: nuper, nupar

      CHARACTER(128) :: netCDF_file

!      namelist/cql3din/nuper, nupar, rho, UmaxPara, UmaxPerp, n_distribution, &
!     &   n_expand, lmax, xnurf, ncount, xktev, xn0, amu, &
!     &   C_matrix_file, n_u_basis, n_theta_basis


      real, dimension(:),   allocatable :: UPERP, UPARA
      real, dimension(:,:), allocatable :: DFDUPER, DFDUPAR
      real, dimension(:,:), allocatable :: f
      real, dimension(:,:,:), allocatable :: f_cql_cart
      real, dimension(:,:,:), allocatable :: df_cql_uprp
      real, dimension(:,:,:), allocatable :: df_cql_uprl


      real :: dens(n_psi_dim)
      real :: xn(n_psi_dim), xkt(n_psi_dim), xkt_prt(n_psi_dim)

      integer :: i_uperp, i_upara
      real :: dfdu, dfdtheta
      real :: d2fdu2, d2fdtheta2

!      double precision :: f_cql_n(0 : n_theta_dim, n_u_dim, n_psi_dim)      
!      double precision :: f_cql_sum(n_theta_dim, n_u_dim, n_psi_dim)      
!      double precision :: theta_dp(n_theta_dim, n_psi_dim)
!      double precision :: pi_dp, piinv_dp
!      double precision :: f_cql_dp(n_theta_dim, n_u_dim, n_psi_dim)

!      real :: dfdu_cql(n_theta_dim, n_u_dim, n_psi_dim)
!      real :: dfdtheta_cql(n_theta_dim, n_u_dim, n_psi_dim)
!      real :: f_aorsa(n_theta_dim, n_u_dim, n_psi_dim)

      real, dimension(:,:,:), allocatable :: dfdu_cql
      real, dimension(:,:,:), allocatable :: dfdtheta_cql
      real, dimension(:,:,:), allocatable :: f_aorsa


      real :: f1(n_u_dim), fout
      
!      real :: f2(n_u_dim, n_theta_dim)
!      real :: dfdu_(n_u_dim, n_theta_dim)
!      real :: dfdth_(n_u_dim, n_theta_dim)
!      real :: f_(n_u_dim, n_theta_dim)
      
      real, dimension(:,:), allocatable :: f2
      real, dimension(:,:), allocatable :: dfdu_
      real, dimension(:,:), allocatable :: dfdth_
      real, dimension(:,:), allocatable :: f_      
      

!      real uxx(nxmx, nymx), uxy(nxmx, nymx), uxz(nxmx, nymx), &
!     &     uyx(nxmx, nymx), uyy(nxmx, nymx), uyz(nxmx, nymx), &
!     &     uzx(nxmx, nymx), uzy(nxmx, nymx), uzz(nxmx, nymx)

!      complex sig2xx, sig2xy, sig2xz, &
!     &        sig2yx, sig2yy, sig2yz, &
!     &        sig2zx, sig2zy, sig2zz

!      real bxn(nxmx, nymx), byn(nxmx, nymx), bzn(nxmx, nymx)
!      real gradprlb(nxmx, nymx), bmod(nxmx, nymx), bmod_mid(nxmx, nymx)
!      real xkti2(nxmx, nymx), xn2a(nxmx, nymx)
!      real capr(nxmx), omgci2(nxmx, nymx), omgp22(nxmx, nymx)

!      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)

      integer nzfun, ibessel, lmax, lmin, ndisti2, nphi, m, nbessj
      integer n_distribution
      real delta0, xnuomg, xmi, qi2, omgrf, xk0, v0i

      real z2, amu2, eta2, a, umax, ve3, ve, alphae, xkte, e0
      real z3, amu3, eta3
      real z1, amu1, eta1, xme, zeff, vmax
      real z_alpha, amu_alpha, eta_alpha

      real xktev, xmih, amu, xn0, xmi_alpha, xm
      real eps0, clight, xmu0
      real q, xnurf, sqx, rho, xlogf, rho10
      complex zi


      real :: theta_(n_theta_dim), u_(n_u_dim)
      real :: u_cql, u_aorsa

      real :: pi, piinv, fl, ugiv, thegiv, fnorm, fmin, fmax, &
     &    vc_mks, u0, alpha, e, u_dbb, u_norm_mks



      common/sigcom/zi, eps0, v0i, omgrf, xk0


      open(unit=237,file='out237',status='unknown', form='formatted')

      open(unit=140, file='f_cql_3D.vtk', status='unknown', &
     &                                              form='formatted')       
      
!      open(unit=115,file='out115',status='unknown', form='formatted')
!      open(unit=143,file='ahern',status='unknown', form='formatted')
      
      open(unit=50,file='cql3d.out',status='unknown', form='formatted')



!      open(unit=163, file='cql3d.in', status='old', form='formatted')
!      rewind(163)

!     --------------------------------
!     set default values of input data:
!     --------------------------------

      C_matrix_file = 'phillips_nstx3.5._C_matrix.dat'
!      netCDF_file = 'phillips_nstx3.5.2.nc'

      n_distribution = 1
      xktev = 17.1e+03
      xn0 = 1.9403e+18
      amu = 1.0

      n_u_basis = 13
      n_theta_basis = 10

      rho = .10
      ncount = 1

      UmaxPara = 1.0
      UmaxPerp = 1.0

      n_expand = 0
      lmax = 4

      xnurf = 60.00E+06


!-----n_distribution: if (n_distribution .eq. 0) Maxwellian distribution is used
!-----                if (n_distribution .eq. 1) CQL3D distribution is used from netCDF_file (default)
!-----                if (n_distribution .eq. 2) Smoothed CQL3D distribution is used from C_matrix_file
!-----                if (n_distribution .eq. 3) Slowing down distribution is used

!-----xktev   : temperature of Maxwellian in eV when n_distribution = 0
!-----xn0     : density of Maxwellian in m-3 when n_distribution = 0
!-----amu     : atomic mass for the Maxwellian species when n_distribution = 0

!-----n_u_basis     : number of u basis functions in smoothed CQL3D distribution
!-----n_theta_basis : number of theta basis functions in smoothed CQL3D distribution

!-----rho:      radial location at which sigma is calculated

!-----nuper:    number of perpendicular velocities in METS integration grid
!-----nupar:    number of parallel velocities in METS integration grid

!-----UmaxPerp: maximum perpendicular velocity on METS integration grid
!-----UmaxPara: maximum parallel velocity on METS integration grid

!-----n_expand: number of terms in cosine(theta) expansion
!-----          if n_expand = 0 (default) no expansion is done


!     ---------------
!     read input data
!     ---------------
!      read (163, cql3din)


!       -----------------
!       initialize arrays
!       -----------------

      allocate( UPERP(nuper) )
      allocate( UPARA(nupar) )
      allocate( dfduper(nuper, nupar) )
      allocate( dfdupar(nuper, nupar) )
      allocate( f(nuper, nupar) )

      allocate( f_cql_cart(nuper, nupar, n_psi_dim) )
      allocate( df_cql_uprp(nuper, nupar, n_psi_dim) )
      allocate( df_cql_uprl(nuper, nupar, n_psi_dim) )
      
      allocate( dfdu_cql(n_theta_dim, n_u_dim, n_psi_dim) )
      allocate( dfdtheta_cql(n_theta_dim, n_u_dim, n_psi_dim) )
      allocate( f_aorsa(n_theta_dim, n_u_dim, n_psi_dim) )

      allocate( f2(n_u_dim, n_theta_dim) )
      allocate( dfdu_(n_u_dim, n_theta_dim) )
      allocate( dfdth_(n_u_dim, n_theta_dim) )
      allocate( f_(n_u_dim, n_theta_dim) )


      e = 1.6e-19
      xmih = 1.67e-27
      pi = 4.0 * atan(1.0)
      piinv = 1.0 / pi

!      pi_dp = 4.0 * datan(1.0d+00)
!      piinv_dp = 1.0 / pi_dp


!     -------------------------------------------------
!     Read CQL3D distribution function from netcdf file
      !     -------------------------------------------------
      kspeci = 1
      CALL netcdfr3d_multigen(netCDF_file,kspeci)

      vc_mks = vc * 1.0e-02
      vc_mks_cql3d = vc_mks

!     ---------------
!     Write some data
!     ---------------

        WRITE(6,*)
        WRITE (6,*) "n_theta_max = ", n_theta_max
        WRITE (6,*) "n_u = ", n_u
        WRITE (6,*) "n_psi = ", n_psi
        WRITE (6,*) "n_theta_(i) = ", n_theta_
        WRITE (6,*) "vc_cql3d_cgs = ", vc, "cm/sec"
        WRITE (6,*) "vc_cql3d_mks = ", vc_mks_cql3d, "m/sec"
        WRITE (6,*) "n_t = ", n_t



      WRITE(6,*)

      write(6, *)
      write(6, *) "rho/a(i_psi)"
      write(6, 310) (rho_a(i_psi), i_psi = 1, n_psi)
      
      write(6, *)
      write(6, *) "wperp_cql(i_psi, n_t) = "
      write(6, 310) (wperp_cql(i_psi, n_t), i_psi = 1, n_psi)
      
      wperp_max = wperp_cql(1, n_t)
      do i_psi = 1, n_psi
         wperp_max = amax1(wperp_max, wperp_cql(i_psi, n_t))
      end do
      
      
      
      write(6, *)
      write(6, *) "wpar_cql(i_psi, n_t) = "
      write(6, 310) (wpar_cql(i_psi, n_t), i_psi = 1, n_psi) 
      
      wpar_max = wpar_cql(1, n_t)
      do i_psi = 1, n_psi
         wpar_max = amax1(wpar_max, wpar_cql(i_psi, n_t))
      end do
      
      write(15, *)"wperp_max = ", wperp_max, " keV"
      write(15, *)"wpar_max = ", wpar_max, " keV"       
            
      write(6, *)"wperp_max = ", wperp_max, " keV"
      write(6, *)"wpar_max = ", wpar_max, " keV" 
      
      enorm_cql3d_joules = .5 * xm * vc_mks_cql3d**2
      enorm_cql3d_kev = enorm_cql3d_joules / e / 1000.
      
!     ----------------------------------------------------
!     Set AORSA enorm to enorm_factor x the maximum energy
!     ----------------------------------------------------      
      enorm_aorsa_kev =  (wperp_max + wpar_max) * enorm_factor      

      
      enorm_aorsa_joules = 1000. * enorm_aorsa_kev * e             
      vc_mks_aorsa = sqrt(2.0 * enorm_aorsa_joules / xm)
                   
      if (enorm_factor .eq. 0.0 .or. enorm_aorsa_kev .gt. enorm_cql3d_kev) then
         vc_mks_aorsa = vc_mks_cql3d 
         enorm_aorsa_joules = .5 * xm * vc_mks_aorsa**2
         enorm_aorsa_kev = enorm_aorsa_joules / e / 1000.
      end if      
      
      vc_ratio = vc_mks_aorsa / vc_mks_cql3d
      
      write(6, *)"enorm_cql3d_kev = ", enorm_cql3d_kev, " keV"       
      write(6, *)"enorm_aorsa_kev = ", enorm_aorsa_kev, " keV" 
      
      write(15, *)"enorm_cql3d_kev = ", enorm_cql3d_kev, " keV"       
      write(15, *)"enorm_aorsa_kev = ", enorm_aorsa_kev, " keV"       
               

      write(6, *)
      write(6, *) "u(i_u) ="
      write(6, 310) (u(i_u), i_u = 1, n_u)


      i_psi = 5
      write(6, *)
      write(6, *) "theta(i_theta, 5) = "
      write(6, 310) (theta(i_theta, i_psi), i_theta = 1, n_theta_(i_psi))
      

      if (n_distribution .eq. 0)then


!        -------------------------------------------------------------------
!        Isotropic Maxwellian with thermal velocity = sqrt(2.0 * xkt / xmi)
!        -------------------------------------------------------------------
         amu = 1.0
         xmi =  xmih * amu
      
         xn0 = 1.022e+19
         xnlim = 4.84e+18
      
         xkt0 = 2.0e+03  * e
         xktlim = .285e+03  * e
      
         alphan = 1.1 
         betan = 1.5
             
         alphat = 1.9
         betat = 1.9
      
!        -----------------------------------------
!        n and T profiles for isotropic Maxwellian
!        -----------------------------------------
         do i_psi = 1, n_psi
         
            shapen = 0.0
            shapet = 0.0

            if(rho_a(i_psi) .le. 1.0)then
                  shapen  = 1.0 - rho_a(i_psi)**betan
                  shapet  = 1.0 - rho_a(i_psi)**betat
            end if

            xn(i_psi) = xnlim  + (xn0 - xnlim)   * shapen**alphan 
            xkt(i_psi)= xktlim + (xkt0 - xktlim) * shapet**alphat 
            xkt_prt(i_psi) = xkt(i_psi) / e
            
         end do
         

         do i_psi = 1, n_psi
         
            alpha = sqrt(2.0 * xkt(i_psi) / xmi)
!           vc_mks = 3.5 * alpha
            u0 = vc_mks / alpha
            fnorm = xn(i_psi) * u0**3 / pi**1.5
            
            
            do  i_u = 1, n_u
               do i_theta = 1, n_theta_(i_psi)
                  f_cql(i_theta, i_u, i_psi) = exp(-u(i_u)**2 * u0**2) &
     &                * fnorm
               end do
            end do
         end do

      end if

!     -------------------------
!     Slowing down distribution
!     -------------------------

      if (n_distribution .eq. 3)then

         z1 = 1.0
         amu1 = 2.0


         z2 = 1.0
         amu2 = 3.0
         eta2 =  .50

         z3 = 2.0
         amu3 = 3.0
         eta3 =  .0


         z_alpha = 2.0
         amu_alpha = 4.0

         xme = 9.11e-31
         xktev = 33.7E+03
         xne = 0.7e+20
         xkte = xktev * e

         eta_alpha =  .02
         xn_alpha = eta_alpha * xne
         xmi_alpha = xmih * amu_alpha
         e0 = 3.5e+06 * e

         eta1 = 1.0 / z1 * (1.0 - z2 * eta2 - z3 * eta3  &
     &                                      - z_alpha * eta_alpha)


        
         vmax = sqrt(2.0 * e0 / xmi_alpha)
         vc_mks = 1.1 * vmax
         alphae = sqrt(2.0 * xkte / xme)
         zeff = z1 / amu1 * eta1 + z2 / amu2 * eta2  &
     &            + z3 / amu3 * eta3 + z_alpha / amu_alpha * eta_alpha
        
         ve3 = 3.0 * sqrt(pi) * xme/xmi_alpha * zeff * alphae**3
         ve = ve3**(1./3.)
         A = 3.0 / (4.0 * pi * alog(1.0 + vmax**3 / ve3))
         umax = vmax / vc_mks

         u0 = vc_mks / ve

!         write(6, *)
!         write(6, *) "u0 = vc_mks / ve =", u0

         fnorm = xn_alpha * u0**3 * A


         do i_psi = 1, n_psi
            do  i_u = 1, n_u
               do i_theta = 1, n_theta_(i_psi)
                  if(u(i_u) .gt. umax) then
                     f_cql(i_theta, i_u, i_psi) = 0.0
                  else
                     f_cql(i_theta, i_u, i_psi) = fnorm / &
     &                (1.0 + u(i_u)**3 * u0**3)
                  end if
               end do
            end do
         end do

      end if



!     ------------------------
!     integrate to get density
!     ------------------------

       do i_psi = 1, n_psi


          do i_theta = 1, n_theta_(i_psi)
             theta_(i_theta) =  theta(i_theta, i_psi)
          end do

          do  i_u = 1, n_u
             do i_theta = 1, n_theta_(i_psi)
                f2(i_u, i_theta) = f_cql(i_theta, i_u, i_psi)
             end do
          end do

          if (i_psi .eq. 1) then
             call a2dmnmx(f2, n_u_dim, n_theta_dim, n_u, n_theta_(i_psi), fmin, fmax)
          end if

          call sgrate_sphere(u, theta, f2, 1, n_u, 1, n_theta_(i_psi), &
     &           dens(i_psi), n_u_dim, n_theta_dim)


       end do

      write(6, *)
      write(6, *) "density ="
      write(6, 310) (dens(i_psi), i_psi = 1, n_psi)

!      write(6, *)
!      write(6, *) "f_cql(1, 2, 1) = ", f_cql(1, 2, 1)

      write(6, *)
      write(6, *) "f_cql_max = ", fmax
      
           
      
!      write(6, *)
!      write(6, *) "temperature ="
!      write(6, 310) (xkt_prt(i_psi), i_psi = 1, n_psi)




!     ----------------------
!     normalize f_cql to 1.0
!     ----------------------

       do i_psi = 1, n_psi

          do  i_u = 1, n_u
             do i_theta = 1, n_theta_(i_psi)
                f_cql(i_theta, i_u, i_psi) = f_cql(i_theta, i_u, i_psi) / dens(i_psi)
             end do
          end do

       end do
       

!     ----------------------
!     calculate derivatives:
!     ----------------------

      do i_psi = 1, n_psi

         do i_theta = 1, n_theta_(i_psi)
            theta_(i_theta) =  theta(i_theta, i_psi)
         end do

         do  i_u = 1, n_u
            do i_theta = 1, n_theta_(i_psi)
               f2(i_u, i_theta) = f_cql(i_theta, i_u, i_psi)
            end do
         end do

         do  i_u = 1, n_u
            do i_theta = 1, n_theta_(i_psi)

               call deriv_1(f2, n_u_dim, n_theta_dim, i_u, i_theta, n_u, n_theta_(i_psi), &
     &                u, dfdu, d2fdu2)
               call deriv_2(f2, n_u_dim, n_theta_dim, i_u, i_theta, n_u, n_theta_(i_psi), &
     &                theta_, dfdtheta, d2fdtheta2)

               dfdu_cql(i_theta, i_u, i_psi) = dfdu
               dfdtheta_cql(i_theta, i_u, i_psi) = dfdtheta

            end do
         end do

      end do
       
       
!     ------------------------------
!     interpolate onto AORSA u grid: 
!     ------------------------------       
       
      do i_psi = 1, n_psi
          do i_theta = 1, n_theta_(i_psi)
       
             do  i_u = 1, n_u             
                f1(i_u) = f_cql(i_theta, i_u, i_psi)
             end do
                     
             do  i_u = 1, n_u 
                u_aorsa = u(i_u)             
                u_cql = u_aorsa * vc_ratio           
                ugiv = u_cql
          
                call intplt_1d(ugiv, fout, n_u, f1, n_u_dim, u)
                f_aorsa(i_theta, i_u, i_psi) = fout * vc_ratio**3
             end do
          
          end do
       end do
              

!     ----------------------------------------------
!     Write data for plotting f_aorsa in u and theta
!     ----------------------------------------------

      write (237, 311) n_u
      write (237, 311) n_psi
      write (237, 310) vc

      write (237, 310) (u(i_u), i_u = 1, n_u)
      write (237, 311) (n_theta_(i_psi), i_psi = 1, n_psi)
      write (237, 310) ((theta(i_theta, i_psi), &
     &          i_theta = 1, n_theta_(i_psi)), i_psi = 1, n_psi)     
           
      write (237, 310) (((f_aorsa(i_theta, i_u, i_psi), &
     &     i_theta = 1, n_theta_(i_psi)), i_u = 1, n_u), i_psi = 1, n_psi)
     
           
!     ----------------------------------------------
!     interpolate onto METS grid: u_perp, u_parallel
!     ----------------------------------------------
      do i_u = 1, n_u
         u_(i_u) =  u(i_u)
      end do

      do i_psi = 1, n_psi

         do i_theta = 1, n_theta_(i_psi)
            theta_(i_theta) =  theta(i_theta, i_psi)
         end do

         do i_u = 1, n_u
             do i_theta = 1, n_theta_(i_psi)
               f_(i_u, i_theta) =  f_cql(i_theta, i_u, i_psi)
               dfdu_(i_u, i_theta) = dfdu_cql(i_theta, i_u, i_psi)
               dfdth_(i_u, i_theta) = dfdtheta_cql(i_theta, i_u, i_psi)
            end do
         end do  


         call mets_grid(nupar, nuper, UminPara, UmaxPara, UmaxPerp, &
     &                UPERP, UPARA, DFDUPER, DFDUPAR, f, fnorm, &
     &                u_, theta_, f_, dfdu_, dfdth_, n_u_dim, n_theta_dim,  &
     &                n_u, n_theta_(i_psi), vc_ratio)
 

         do i_uperp = 1, NUPER
            do i_upara = 1, NUPAR

               f_cql_cart(i_uperp, i_upara, i_psi) = f(i_uperp, i_upara)
               df_cql_uprp(i_uperp, i_upara, i_psi)= DFDUPER(i_uperp, i_upara)
               df_cql_uprl(i_uperp, i_upara, i_psi)= DFDUPAR(i_uperp, i_upara)

               if(abs(f_cql_cart(i_uperp, i_upara, i_psi)) .lt. 1.0e-99) &
     &                               f_cql_cart(i_uperp, i_upara, i_psi) = 0.0
               if(abs(df_cql_uprp(i_uperp, i_upara, i_psi)) .lt. 1.0e-99) &
     &                               df_cql_uprp(i_uperp, i_upara, i_psi) = 0.0
               if(abs(df_cql_uprl(i_uperp, i_upara, i_psi)) .lt. 1.0e-99) &
     &                               df_cql_uprl(i_uperp, i_upara, i_psi) = 0.0

            end do
         end do

      end do

!     ------------------------
!     integrate to get density
!     ------------------------

      do i_psi = 1, n_psi

         do i_uperp = 1, NUPER
            do i_upara = 1, NUPAR
               f(i_uperp, i_upara)= f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

          if (i_psi .eq. 1) then
             call a2dmnmx(f, NUPER, NUPAR, nuper, nupar, fmin, fmax)
          end if

         call  sgrate_cyl(uperp, upara, f, 1, nuper, 1, nupar, &
     &      dens(i_psi), NUPER, NUPAR)

      end do


!      write(6, *)
!      write(6, *) "uperp ="
!      write(6, 310) (uperp(i_uperp), i_uperp = 1, nuper)

!      write(6, *)
!      write(6, *) "upara ="
!      write(6, 310) (upara(i_upara), i_upara = 1, nupar)
      
      

      write(6, *)
      write(6, *) "density ="
      write(6, 310) (dens(i_psi), i_psi = 1, n_psi)
      
!      call exit      

      write(6, *)
      write(6, *) "f_cql_cart_max = ", fmax

      write(6, *)
      write(6, *) "rho/a(i_psi)"
      write(6, 310) (rho_a(i_psi), i_psi = 1, n_psi)
      write(6, *)




!     ------------------------------------------------
!     Write data for AORSA2D in u_perp and u_parallel
!     ------------------------------------------------
      rewind (50)
      
      write (50, 309) nuper
      write (50, 309) nupar      
      write (50, 309) n_psi

      write (50, 3310) vc_mks_aorsa
      write (50, 3310) UminPara, UmaxPara

      write (50, 3310) (rho_a(i_psi), i_psi = 1, n_psi)
      write (50, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
      write (50, 3310) (upara(i_upara), i_upara = 1, nupar)

      write (50, 3310) (((f_cql_cart(i_uperp, i_upara, i_psi), &
     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

      write (50, 3310) (((df_cql_uprp(i_uperp, i_upara, i_psi), &
     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

      write (50, 3310) (((df_cql_uprl(i_uperp, i_upara, i_psi), &
     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)
     

      rewind (50)

!        read (50, 309) nuper
!        read (50, 309) nupar

!      read (50, 309) n_psi

!      read (50, 3310) vc_mks
!      read (50, 3310) UminPara, UmaxPara

!      read (50, 3310) (rho_a(i_psi), i_psi = 1, n_psi)
!      read (50, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
!      read (50, 3310) (upara(i_upara), i_upara = 1, nupar)

!      read (50, 3310) (((f_cql_cart(i_uperp, i_upara, i_psi), &
!     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

!      read (50, 3310) (((df_cql_uprp(i_uperp, i_upara, i_psi), &
!     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

!      read (50, 3310) (((df_cql_uprl(i_uperp, i_upara, i_psi), &
!     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)


!     ------------------------------------------------
!     get CQL3D distribution function on the midplane
!     ------------------------------------------------
!      call cql3d_dist(nupar, nuper, n_psi,  &
!     &              n_psi_dim, rho_a, rho,  &
!     &              UminPara,UmaxPara,      &
!     &              df_cql_uprp, df_cql_uprl,  &
!     &              UPERP, UPARA, DFDUPER, DFDUPAR)
     
!      do ni = 1, nuper
!          do mi = 1, nupar
!     
!             dfdth = upara(mi) * dfduper(ni, mi)  &
!     &                           - uperp(ni) * dfdupar(ni, mi)
!     
!            write(6,  1314)ni, mi, uperp(ni),  upara(mi), dfdth
!            write(115, 1314)ni, mi, uperp(ni),  upara(mi), dfdth
!                 
!         end do
!      end do
               
! 1314  format(2i10, 1p,9e12.4)  

!     ------------------------------------------------
!     Write data for plotting in u_perp and u_parallel
!     ------------------------------------------------

      write (237, 309) NUPER
      write (237, 309) NUPAR
      write (237, 309) n_psi

      write (237, 310) (uperp(i_uperp), i_uperp = 1, nuper)
      write (237, 310) (upara(i_upara), i_upara = 1, nupar)
      write (237, 310) (((f_cql_cart(i_uperp, i_upara, i_psi), &
     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)
     
      write(237, 310) (wperp_cql(i_psi, n_t), i_psi = 1, n_psi)
      write(237, 310) (wpar_cql(i_psi, n_t), i_psi = 1, n_psi)
     
      close (237)


 
!     ----------------------
!     take log of f_cql_cart
!     ----------------------
      do i_psi = 1, n_psi      
         do i_uperp = 1, nuper
            do i_upara = 1, nupar
               
               if(f_cql_cart(i_uperp, i_upara, i_psi) .le. 0.0) then
                  f_cql_cart(i_uperp, i_upara, i_psi) = -100.0
               end if
               
               if(f_cql_cart(i_uperp, i_upara, i_psi) .gt. 0.0) then
                  f_cql_cart(i_uperp, i_upara, i_psi) =  &
     &                    log10(f_cql_cart(i_uperp, i_upara, i_psi))
               end if

            end do
         end do
      end do 
     
      
!     ---------------
!     expand rho grid
!     ----------------
      do i_psi = 1, n_psi      
          rho_a(i_psi) = rho_a(i_psi) * 5.0      
      end do 
 
!     ---------------------------------------
!     Write f_cql_3D.vtk file "rectilinear_grid"
!     ---------------------------------------
      number_points = nupar * nuper * n_psi
      
      write(140, 2840)
 2840 format('# vtk DataFile Version 2.0')
 
      write(140, 2845)
 2845 format('Cql3d distribution function')      
      
      write(140, 2846)
 2846 format('ASCII')       
           
      write(140, 2847)
 2847 format('DATASET RECTILINEAR_GRID')  
      
      write(140, 2841)nupar, nuper, n_psi
 2841 format('DIMENSIONS', 3i8)
           
      write(140, 2842) nupar
 2842 format('X_COORDINATES', i8, ' float')
      write (140, 3411) (upara(i_upara), i_upara = 1, nupar) 
 
      write(140, 2843) nuper
 2843 format('Y_COORDINATES', i8, ' float')
      write (140, 3411) (uperp(i_uperp), i_uperp = 1, nuper) 
      
      write(140, 2851) n_psi
 2851 format('Z_COORDINATES', i8, ' float')
      write (140, 3411) (rho_a(i_psi), i_psi = 1, n_psi)      
      
      write(140, 2844) number_points
 2844 format('POINT_DATA', i10)      
      
      write(140, 2848)
 2848 format('SCALARS f_cql_cart float 1')       
      
      write(140, 2849)
 2849 format('LOOKUP_TABLE default')        
                  
      write (140, 3411) (((f_cql_cart(i_uperp, i_upara, i_psi), &
     &   i_upara = 1, nupar), i_uperp = 1, nuper), i_psi = 1, n_psi) 
     
     
 3410 format(1p,4e10.2)  
 3411 format(6f12.4)     
     
 
                 
!     ----------------
!     Write ahern file
!     ----------------          
!      write (143, 311) nuper
!      write (143, 311) nupar
!      write (143, 311) n_psi

!      write (143, 310) (uperp(i_uperp), i_uperp = 1, nuper)
!      write (143, 310) (upara(i_upara), i_upara = 1, nupar)
!      write (143, 310) (rho_a(i_psi), i_psi = 1, n_psi)
      
      
!      write (143, 310) (((f_cql_cart(i_uperp, i_upara, i_psi), &
!     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)



!      go to 5000

!      n = 0
!      m = 0
!      i = 5
!      j = 5

!      zi = cmplx(0.0, 1.0)
!      eps0 = 8.85e-12
!      xmu0 = 1.26e-06
!      clight = 1.0 / sqrt(eps0 * xmu0)
!      pi = 3.141592654

!      omgrf = 2.0 * pi * xnurf
!      xk0 = omgrf / clight

!      gradprlb(i,j) = 0.0


!     -----------------
!     DIII-D parameters
!     -----------------

!      qi2 = 1.6e-19
!      xn2a(i,j) = 0.112920789750647245E+19
!      xnuomg = 0.0
!      xkti2(i,j) = xkt0


!      omgp22(i,j) =  xn2a(i,j) * qi2**2 / (eps0 * xmi)


!      lmin = -lmax
!      nbessj = lmax + 2

!      capr(i) =  1.56400000000000006
!      nphi = -15

!      xkxsav(n) =  0.000000000000000000E+00
!      xkysav(m) =  0.000000000000000000E+00
!      bxn(i,j) =  0.978610000000000070E-02
!      byn(i,j) =  -0.837220000000000047E-01
!      bzn(i,j) =  0.996439999999999992
!      bmod(i,j) =  2.01040000000000019
!      bmod_mid(i,j) = 1.62991393393751482


!      omgci2(i,j) = qi2 * bmod(i,j) / xmi

!     ---------------------------
!     Calculate rotation matrix U
!     ---------------------------

!      sqx = sqrt(1.0 - bxn(i,j)**2)

!      uxx(i, j) =   sqx
!      uxy(i, j) = - bxn(i, j) * byn(i, j) / sqx
!      uxz(i, j) = - bxn(i, j) * bzn(i, j) / sqx
!      uyx(i, j) =   0.0
!      uyy(i, j) =   bzn(i, j) / sqx
!      uyz(i, j) = - byn(i, j) / sqx
!      uzx(i, j) =   bxn(i, j)
!      uzy(i, j) =   byn(i, j)
!      uzz(i, j) =   bzn(i, j)

!      nzfun = 0

!      myid = 0
!      nproc = 1
!     -------------------------------------
!     setup Z function table when nzfun = 3
!     -------------------------------------
!      if (nzfun .eq. 3) then
!         iflag = 0
!         call ztable_setup(myid, iflag)
!         if(iflag .eq. 1)call ztable_setup(myid, nproc)
!      end if


!      ndisti2 = 0

!      call sigmad_cql3d(i, j, n, m, rho, rho_a,  &
!     &   gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
!     &   xmi, qi2, xn2a(i,j), xnuomg,  &
!     &   xkti2(i,j), omgci2(i,j), omgp22(i,j), &
!     &   -lmax, lmax, nzfun, ibessel,          &
!     &   xkxsav(n), xkysav(m), nphi, capr(i),  &
!     &   bxn(i,j), byn(i,j), bzn(i,j), &
!     &   uxx(i,j), uxy(i,j), uxz(i,j), &
!     &   uyx(i,j), uyy(i,j), uyz(i,j), &
!     &   uzx(i,j), uzy(i,j), uzz(i,j), &
!     &   sig2xx, sig2xy, sig2xz, &
!     &   sig2yx, sig2yy, sig2yz, &
!     &   sig2zx, sig2zy, sig2zz, &
!     &   delta0, ndisti2, nupar, nuper, n_psi, &
!     &   n_psi_dim, &
!     &   UminPara, UmaxPara, UPERP, UPARA, &
!     &   vc_mks, df_cql_uprp, df_cql_uprl, nbessj, ncount)

!       ndisti2 = 1

!      call sigmad_cql3d(i, j, n, m, rho, rho_a,  &
!     &   gradprlb(i,j), bmod(i,j), bmod_mid(i,j),  &
!     &   xmi, qi2, xn2a(i,j), xnuomg,  &
!     &   xkti2(i,j), omgci2(i,j), omgp22(i,j), &
!     &   -lmax, lmax, nzfun, ibessel,          &
!     &   xkxsav(n), xkysav(m), nphi, capr(i),  &
!     &   bxn(i,j), byn(i,j), bzn(i,j), &
!     &   uxx(i,j), uxy(i,j), uxz(i,j), &
!     &   uyx(i,j), uyy(i,j), uyz(i,j), &
!     &   uzx(i,j), uzy(i,j), uzz(i,j), &
!     &   sig2xx, sig2xy, sig2xz, &
!     &   sig2yx, sig2yy, sig2yz, &
!     &   sig2zx, sig2zy, sig2zz, &
!     &   delta0, ndisti2, nupar, nuper, n_psi, &
!     &   n_psi_dim, &
!     &   UminPara, UmaxPara, UPERP, UPARA, &
!     &   vc_mks, df_cql_uprp, df_cql_uprl, nbessj, ncount)


  310 format(1p,6e12.4)
 3310 format(1p,6e18.10)
  311 format(10i10)
  309 format(10i10)

  312 format(i10, 1p,6e12.4)
  
  
! 5000 continue


!     -----------------
!     deallocate arrays
!     -----------------

      deallocate( UPERP )
      deallocate( UPARA )
      deallocate( dfduper )
      deallocate( dfdupar )
      deallocate( f )

      deallocate( f_cql_cart )
      deallocate( df_cql_uprp )
      deallocate( df_cql_uprl )
      
      deallocate( dfdu_cql )
      deallocate( dfdtheta_cql )
      deallocate( f_aorsa ) 
      
      deallocate( f2 )
      deallocate( dfdu_ )
      deallocate( dfdth_ )
      deallocate( f_ )           
      
      close(50)
      close (140)


      return
      END subroutine cql3d_setup1

!
!***********************************************************************
!
      subroutine cql3d_setup2(netcdf_file, nuper, nupar, xm, enorm_factor, vc_mks_cql3d)
      use size_mod   

      USE read_CQL3D
      USE quick_vector_write_m
      USE f_expanded_m, only : C_matrix_file, in_file_CQL, &
     &  n_u_coeffs, n_theta_coeffs, u_norm, u0_Max_poly, &
     &  C_matrix, dens_Maxwell_psi, T_Maxwell_psi,              &
     &  read_C_matrix_m, eval_delta_f_exp_m, eval_f_exp_m
        
      USE basis_functions_m, only : u_fns_name
        

      IMPLICIT NONE
      
      INTEGER :: j, k, n, n_expand, ncount, number_points
      integer i, i_psi, i_theta, i_u, iflag, nproc, myid
      integer n_u_basis, n_theta_basis, ni, mi, kspeci
      real enorm_factor

      integer, parameter :: n_u_dim = 300
      integer, parameter :: n_theta_dim = 200
      integer, parameter :: n_psi_dim = 64

!     ----------------
!     450 x 450 modes:
!     ----------------
!      integer, parameter :: nmodesmax = 450
!      integer, parameter :: mmodesmax = 450

      integer, parameter :: nxmx = nmodesmax
      integer, parameter :: nymx = mmodesmax
      integer, parameter :: nrhomax = nmodesmax * 2

      integer, parameter :: nkdim1 = - nmodesmax / 2
      integer, parameter :: nkdim2 =   nmodesmax / 2

      integer, parameter :: mkdim1 = - mmodesmax / 2
      integer, parameter :: mkdim2 =   mmodesmax / 2

      real :: UminPara, UmaxPara, UmaxPerp, dupara, duperp, drho, drho10
      real fb, dfb_du, dfb_dtheta, xn_alpha, xne, dfdth
      real betat, alphat, shapen, shapet, xkt0, xnlim, xktlim, betan, alphan
      real wperp_max, wpar_max, enorm_aorsa_kev
      real enorm_cql3d_joules, enorm_cql3d_kev
      real enorm_aorsa_joules, vc_mks_aorsa, vc_ratio, vc_mks_cql3d  

      integer :: nuper, nupar

      CHARACTER(128) :: netCDF_file

!      namelist/cql3din/nuper, nupar, rho, UmaxPara, UmaxPerp, n_distribution, &
!     &   n_expand, lmax, xnurf, ncount, xktev, xn0, amu, &
!     &   C_matrix_file, n_u_basis, n_theta_basis


      real, dimension(:),   allocatable :: UPERP, UPARA
      real, dimension(:,:), allocatable :: DFDUPER, DFDUPAR
      real, dimension(:,:), allocatable :: f
      real, dimension(:,:,:), allocatable :: f_cql_cart
      real, dimension(:,:,:), allocatable :: df_cql_uprp
      real, dimension(:,:,:), allocatable :: df_cql_uprl


      real :: dens(n_psi_dim)
      real :: xn(n_psi_dim), xkt(n_psi_dim), xkt_prt(n_psi_dim)

      integer :: i_uperp, i_upara


      real :: dfdu, dfdtheta
      real :: d2fdu2, d2fdtheta2

!      double precision :: f_cql_n(0 : n_theta_dim, n_u_dim, n_psi_dim)      
!      double precision :: f_cql_sum(n_theta_dim, n_u_dim, n_psi_dim)      
!      double precision :: theta_dp(n_theta_dim, n_psi_dim)
!      double precision :: pi_dp, piinv_dp
!      double precision :: f_cql_dp(n_theta_dim, n_u_dim, n_psi_dim)

!      real :: dfdu_cql(n_theta_dim, n_u_dim, n_psi_dim)
!      real :: dfdtheta_cql(n_theta_dim, n_u_dim, n_psi_dim)
!      real :: f_aorsa(n_theta_dim, n_u_dim, n_psi_dim)

      real, dimension(:,:,:), allocatable :: dfdu_cql
      real, dimension(:,:,:), allocatable :: dfdtheta_cql
      real, dimension(:,:,:), allocatable :: f_aorsa


      real :: f1(n_u_dim), fout
      real :: f2(n_u_dim, n_theta_dim)

      real :: dfdu_(n_u_dim, n_theta_dim)
      real :: dfdth_(n_u_dim, n_theta_dim)
      real :: f_(n_u_dim, n_theta_dim)

!      real uxx(nxmx, nymx), uxy(nxmx, nymx), uxz(nxmx, nymx), &
!     &     uyx(nxmx, nymx), uyy(nxmx, nymx), uyz(nxmx, nymx), &
!     &     uzx(nxmx, nymx), uzy(nxmx, nymx), uzz(nxmx, nymx)

!      complex sig2xx, sig2xy, sig2xz, &
!     &        sig2yx, sig2yy, sig2yz, &
!     &        sig2zx, sig2zy, sig2zz

!      real bxn(nxmx, nymx), byn(nxmx, nymx), bzn(nxmx, nymx)
!      real gradprlb(nxmx, nymx), bmod(nxmx, nymx), bmod_mid(nxmx, nymx)
!      real xkti2(nxmx, nymx), xn2a(nxmx, nymx)
!      real capr(nxmx), omgci2(nxmx, nymx), omgp22(nxmx, nymx)

!      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)

      integer nzfun, ibessel, lmax, lmin, ndisti2, nphi, m, nbessj
      integer n_distribution
      real delta0, xnuomg, xmi, qi2, omgrf, xk0, v0i

      real z2, amu2, eta2, a, umax, ve3, ve, alphae, xkte, e0
      real z3, amu3, eta3
      real z1, amu1, eta1, xme, zeff, vmax
      real z_alpha, amu_alpha, eta_alpha

      real xktev, xmih, amu, xn0, xmi_alpha, xm
      real eps0, clight, xmu0
      real q, xnurf, sqx, rho, xlogf, rho10
      complex zi


      real :: theta_(n_theta_dim), u_(n_u_dim)
      real :: u_cql, u_aorsa

      real :: pi, piinv, fl, ugiv, thegiv, fnorm, fmin, fmax, &
     &    vc_mks, u0, alpha, e, u_dbb, u_norm_mks



      common/sigcom/zi, eps0, v0i, omgrf, xk0


      open(unit=238,file='out238',status='unknown', form='formatted')
!      open(unit=139,file='Point.3D',status='unknown', form='formatted')
      open(unit=140, file='f_cql_3D.vtk', status='unknown', &
     &                                              form='formatted')       
      
!      open(unit=115,file='out115',status='unknown', form='formatted')
!      open(unit=143,file='ahern',status='unknown', form='formatted')
      
      open(unit=50,file='cql3d.out',status='unknown', form='formatted')



!      open(unit=163, file='cql3d.in', status='old', form='formatted')
!      rewind(163)

!     --------------------------------
!     set default values of input data:
!     --------------------------------

      C_matrix_file = 'phillips_nstx3.5._C_matrix.dat'
!      netCDF_file = 'phillips_nstx3.5.2.nc'

      n_distribution = 1
      xktev = 17.1e+03
      xn0 = 1.9403e+18
      amu = 1.0

      n_u_basis = 13
      n_theta_basis = 10

      rho = .10
      ncount = 1

      UmaxPara = 1.0
      UmaxPerp = 1.0

      n_expand = 0
      lmax = 4

      xnurf = 60.00E+06


!-----n_distribution: if (n_distribution .eq. 0) Maxwellian distribution is used
!-----                if (n_distribution .eq. 1) CQL3D distribution is used from netCDF_file (default)
!-----                if (n_distribution .eq. 2) Smoothed CQL3D distribution is used from C_matrix_file
!-----                if (n_distribution .eq. 3) Slowing down distribution is used

!-----xktev   : temperature of Maxwellian in eV when n_distribution = 0
!-----xn0     : density of Maxwellian in m-3 when n_distribution = 0
!-----amu     : atomic mass for the Maxwellian species when n_distribution = 0

!-----n_u_basis     : number of u basis functions in smoothed CQL3D distribution
!-----n_theta_basis : number of theta basis functions in smoothed CQL3D distribution

!-----rho:      radial location at which sigma is calculated

!-----nuper:    number of perpendicular velocities in METS integration grid
!-----nupar:    number of parallel velocities in METS integration grid

!-----UmaxPerp: maximum perpendicular velocity on METS integration grid
!-----UmaxPara: maximum parallel velocity on METS integration grid

!-----n_expand: number of terms in cosine(theta) expansion
!-----          if n_expand = 0 (default) no expansion is done


!     ---------------
!     read input data
!     ---------------
!      read (163, cql3din)


!       -----------------
!       initialize arrays
!       -----------------

      allocate( UPERP(nuper) )
      allocate( UPARA(nupar) )
      allocate( dfduper(nuper, nupar) )
      allocate( dfdupar(nuper, nupar) )
      allocate( f(nuper, nupar) )

      allocate( f_cql_cart(nuper, nupar, n_psi_dim) )
      allocate( df_cql_uprp(nuper, nupar, n_psi_dim) )
      allocate( df_cql_uprl(nuper, nupar, n_psi_dim) )

      allocate( dfdu_cql(n_theta_dim, n_u_dim, n_psi_dim) )
      allocate( dfdtheta_cql(n_theta_dim, n_u_dim, n_psi_dim) )
      allocate( f_aorsa(n_theta_dim, n_u_dim, n_psi_dim) )


      e = 1.6e-19
      xmih = 1.67e-27
      pi = 4.0 * atan(1.0)
      piinv = 1.0 / pi

!      pi_dp = 4.0 * datan(1.0d+00)
!      piinv_dp = 1.0 / pi_dp


!     -------------------------------------------------
!     Read CQL3D distribution function from netcdf file
!     -------------------------------------------------
      kspeci=2
      CALL netcdfr3d_multigen(netCDF_file,kspeci)

      vc_mks = vc * 1.0e-02
      vc_mks_cql3d = vc_mks

!     ---------------
!     Write some data
!     ---------------

        WRITE(6,*)
        WRITE (6,*) "n_theta_max = ", n_theta_max
        WRITE (6,*) "n_u = ", n_u
        WRITE (6,*) "n_psi = ", n_psi
        WRITE (6,*) "n_theta_(i) = ", n_theta_
        WRITE (6,*) "vc_cql3d_cgs = ", vc, "cm/sec"
        WRITE (6,*) "vc_cql3d_mks = ", vc_mks_cql3d, "m/sec"
        WRITE (6,*) "n_t = ", n_t



      WRITE(6,*)

      write(6, *)
      write(6, *) "rho/a(i_psi)"
      write(6, 310) (rho_a(i_psi), i_psi = 1, n_psi)
      
      write(6, *)
      write(6, *) "wperp_cql(i_psi, n_t) = "
      write(6, 310) (wperp_cql(i_psi, n_t), i_psi = 1, n_psi)
      
      wperp_max = wperp_cql(1, n_t)
      do i_psi = 1, n_psi
         wperp_max = amax1(wperp_max, wperp_cql(i_psi, n_t))
      end do
      
      
      
      write(6, *)
      write(6, *) "wpar_cql(i_psi, n_t) = "
      write(6, 310) (wpar_cql(i_psi, n_t), i_psi = 1, n_psi) 
      
      wpar_max = wpar_cql(1, n_t)
      do i_psi = 1, n_psi
         wpar_max = amax1(wpar_max, wpar_cql(i_psi, n_t))
      end do
      
      write(15, *)"wperp_max = ", wperp_max, " keV"
      write(15, *)"wpar_max = ", wpar_max, " keV"       
            
      write(6, *)"wperp_max = ", wperp_max, " keV"
      write(6, *)"wpar_max = ", wpar_max, " keV" 
      
      enorm_cql3d_joules = .5 * xm * vc_mks_cql3d**2
      enorm_cql3d_kev = enorm_cql3d_joules / e / 1000.
      
!     ----------------------------------------------------
!     Set AORSA enorm to enorm_factor x the maximum energy
!     ----------------------------------------------------      
      enorm_aorsa_kev =  (wperp_max + wpar_max) * enorm_factor      

      
      enorm_aorsa_joules = 1000. * enorm_aorsa_kev * e             
      vc_mks_aorsa = sqrt(2.0 * enorm_aorsa_joules / xm)
                   
      if (enorm_factor .eq. 0.0 .or. enorm_aorsa_kev .gt. enorm_cql3d_kev) then
         vc_mks_aorsa = vc_mks_cql3d 
         enorm_aorsa_joules = .5 * xm * vc_mks_aorsa**2
         enorm_aorsa_kev = enorm_aorsa_joules / e / 1000.
      end if      
      
      vc_ratio = vc_mks_aorsa / vc_mks_cql3d
      
      write(6, *)"enorm_cql3d_kev = ", enorm_cql3d_kev, " keV"       
      write(6, *)"enorm_aorsa_kev = ", enorm_aorsa_kev, " keV" 
      
      write(15, *)"enorm_cql3d_kev = ", enorm_cql3d_kev, " keV"       
      write(15, *)"enorm_aorsa_kev = ", enorm_aorsa_kev, " keV"       
               

      write(6, *)
      write(6, *) "u(i_u) ="
      write(6, 310) (u(i_u), i_u = 1, n_u)


      i_psi = 5
      write(6, *)
      write(6, *) "theta(i_theta, 5) = "
      write(6, 310) (theta(i_theta, i_psi), i_theta = 1, n_theta_(i_psi))
      

      if (n_distribution .eq. 0)then


!        -------------------------------------------------------------------
!        Isotropic Maxwellian with thermal velocity = sqrt(2.0 * xkt / xmi)
!        -------------------------------------------------------------------
         amu = 1.0
         xmi =  xmih * amu
      
         xn0 = 1.022e+19
         xnlim = 4.84e+18
      
         xkt0 = 2.0e+03  * e
         xktlim = .285e+03  * e
      
         alphan = 1.1 
         betan = 1.5
             
         alphat = 1.9
         betat = 1.9
      
!        -----------------------------------------
!        n and T profiles for isotropic Maxwellian
!        -----------------------------------------
         do i_psi = 1, n_psi
         
            shapen = 0.0
            shapet = 0.0

            if(rho_a(i_psi) .le. 1.0)then
                  shapen  = 1.0 - rho_a(i_psi)**betan
                  shapet  = 1.0 - rho_a(i_psi)**betat
            end if

            xn(i_psi) = xnlim  + (xn0 - xnlim)   * shapen**alphan 
            xkt(i_psi)= xktlim + (xkt0 - xktlim) * shapet**alphat 
            xkt_prt(i_psi) = xkt(i_psi) / e
            
         end do
         

         do i_psi = 1, n_psi
         
            alpha = sqrt(2.0 * xkt(i_psi) / xmi)
!           vc_mks = 3.5 * alpha
            u0 = vc_mks / alpha
            fnorm = xn(i_psi) * u0**3 / pi**1.5
            
            
            do  i_u = 1, n_u
               do i_theta = 1, n_theta_(i_psi)
                  f_cql(i_theta, i_u, i_psi) = exp(-u(i_u)**2 * u0**2) &
     &                * fnorm
               end do
            end do
         end do

      end if

!     -------------------------
!     Slowing down distribution
!     -------------------------

      if (n_distribution .eq. 3)then

         z1 = 1.0
         amu1 = 2.0


         z2 = 1.0
         amu2 = 3.0
         eta2 =  .50

         z3 = 2.0
         amu3 = 3.0
         eta3 =  .0


         z_alpha = 2.0
         amu_alpha = 4.0

         xme = 9.11e-31
         xktev = 33.7E+03
         xne = 0.7e+20
         xkte = xktev * e

         eta_alpha =  .02
         xn_alpha = eta_alpha * xne
         xmi_alpha = xmih * amu_alpha
         e0 = 3.5e+06 * e

         eta1 = 1.0 / z1 * (1.0 - z2 * eta2 - z3 * eta3  &
     &                                      - z_alpha * eta_alpha)


        
         vmax = sqrt(2.0 * e0 / xmi_alpha)
         vc_mks = 1.1 * vmax
         alphae = sqrt(2.0 * xkte / xme)
         zeff = z1 / amu1 * eta1 + z2 / amu2 * eta2  &
     &            + z3 / amu3 * eta3 + z_alpha / amu_alpha * eta_alpha
        
         ve3 = 3.0 * sqrt(pi) * xme/xmi_alpha * zeff * alphae**3
         ve = ve3**(1./3.)
         A = 3.0 / (4.0 * pi * alog(1.0 + vmax**3 / ve3))
         umax = vmax / vc_mks

         u0 = vc_mks / ve

!         write(6, *)
!         write(6, *) "u0 = vc_mks / ve =", u0

         fnorm = xn_alpha * u0**3 * A


         do i_psi = 1, n_psi
            do  i_u = 1, n_u
               do i_theta = 1, n_theta_(i_psi)
                  if(u(i_u) .gt. umax) then
                     f_cql(i_theta, i_u, i_psi) = 0.0
                  else
                     f_cql(i_theta, i_u, i_psi) = fnorm / &
     &                (1.0 + u(i_u)**3 * u0**3)
                  end if
               end do
            end do
         end do

      end if



!     ------------------------
!     integrate to get density
!     ------------------------

       do i_psi = 1, n_psi


          do i_theta = 1, n_theta_(i_psi)
             theta_(i_theta) =  theta(i_theta, i_psi)
          end do

          do  i_u = 1, n_u
             do i_theta = 1, n_theta_(i_psi)
                f2(i_u, i_theta) = f_cql(i_theta, i_u, i_psi)
             end do
          end do

          if (i_psi .eq. 1) then
             call a2dmnmx(f2, n_u_dim, n_theta_dim, n_u, n_theta_(i_psi), fmin, fmax)
          end if

          call sgrate_sphere(u, theta, f2, 1, n_u, 1, n_theta_(i_psi), &
     &           dens(i_psi), n_u_dim, n_theta_dim)


       end do

      write(6, *)
      write(6, *) "density ="
      write(6, 310) (dens(i_psi), i_psi = 1, n_psi)

!      write(6, *)
!      write(6, *) "f_cql(1, 2, 1) = ", f_cql(1, 2, 1)

      write(6, *)
      write(6, *) "f_cql_max = ", fmax
      
           
      
!      write(6, *)
!      write(6, *) "temperature ="
!      write(6, 310) (xkt_prt(i_psi), i_psi = 1, n_psi)




!     ----------------------
!     normalize f_cql to 1.0
!     ----------------------

       do i_psi = 1, n_psi

          do  i_u = 1, n_u
             do i_theta = 1, n_theta_(i_psi)
                f_cql(i_theta, i_u, i_psi) = f_cql(i_theta, i_u, i_psi) / dens(i_psi)
             end do
          end do

       end do
       

!     ----------------------
!     calculate derivatives:
!     ----------------------

      do i_psi = 1, n_psi

         do i_theta = 1, n_theta_(i_psi)
            theta_(i_theta) =  theta(i_theta, i_psi)
         end do

         do  i_u = 1, n_u
            do i_theta = 1, n_theta_(i_psi)
               f2(i_u, i_theta) = f_cql(i_theta, i_u, i_psi)
            end do
         end do

         do  i_u = 1, n_u
            do i_theta = 1, n_theta_(i_psi)

               call deriv_1(f2, n_u_dim, n_theta_dim, i_u, i_theta, n_u, n_theta_(i_psi), &
     &                u, dfdu, d2fdu2)
               call deriv_2(f2, n_u_dim, n_theta_dim, i_u, i_theta, n_u, n_theta_(i_psi), &
     &                theta_, dfdtheta, d2fdtheta2)

               dfdu_cql(i_theta, i_u, i_psi) = dfdu
               dfdtheta_cql(i_theta, i_u, i_psi) = dfdtheta

            end do
         end do

      end do
       
       
!     ------------------------------
!     interpolate onto AORSA u grid: 
!     ------------------------------       
       
      do i_psi = 1, n_psi
          do i_theta = 1, n_theta_(i_psi)
       
             do  i_u = 1, n_u             
                f1(i_u) = f_cql(i_theta, i_u, i_psi)
             end do
                     
             do  i_u = 1, n_u 
                u_aorsa = u(i_u)             
                u_cql = u_aorsa * vc_ratio           
                ugiv = u_cql
          
                call intplt_1d(ugiv, fout, n_u, f1, n_u_dim, u)
                f_aorsa(i_theta, i_u, i_psi) = fout * vc_ratio**3
             end do
          
          end do
       end do
              

!     ----------------------------------------------
!     Write data for plotting f_aorsa in u and theta
!     ----------------------------------------------

      write (238, 311) n_u
      write (238, 311) n_psi
      write (238, 310) vc

      write (238, 310) (u(i_u), i_u = 1, n_u)
      write (238, 311) (n_theta_(i_psi), i_psi = 1, n_psi)
      write (238, 310) ((theta(i_theta, i_psi), &
     &          i_theta = 1, n_theta_(i_psi)), i_psi = 1, n_psi)     
           
      write (238, 310) (((f_aorsa(i_theta, i_u, i_psi), &
     &     i_theta = 1, n_theta_(i_psi)), i_u = 1, n_u), i_psi = 1, n_psi)
     
           
!     ----------------------------------------------
!     interpolate onto METS grid: u_perp, u_parallel
!     ----------------------------------------------
      do i_u = 1, n_u
         u_(i_u) =  u(i_u)
      end do

      do i_psi = 1, n_psi

         do i_theta = 1, n_theta_(i_psi)
            theta_(i_theta) =  theta(i_theta, i_psi)
         end do

         do i_u = 1, n_u
             do i_theta = 1, n_theta_(i_psi)
               f_(i_u, i_theta) =  f_cql(i_theta, i_u, i_psi)
               dfdu_(i_u, i_theta) = dfdu_cql(i_theta, i_u, i_psi)
               dfdth_(i_u, i_theta) = dfdtheta_cql(i_theta, i_u, i_psi)
            end do
         end do  


         call mets_grid(nupar, nuper, UminPara, UmaxPara, UmaxPerp, &
     &                UPERP, UPARA, DFDUPER, DFDUPAR, f, fnorm, &
     &                u_, theta_, f_, dfdu_, dfdth_, n_u_dim, n_theta_dim,  &
     &                n_u, n_theta_(i_psi), vc_ratio)
 

         do i_uperp = 1, NUPER
            do i_upara = 1, NUPAR

               f_cql_cart(i_uperp, i_upara, i_psi) = f(i_uperp, i_upara)
               df_cql_uprp(i_uperp, i_upara, i_psi)= DFDUPER(i_uperp, i_upara)
               df_cql_uprl(i_uperp, i_upara, i_psi)= DFDUPAR(i_uperp, i_upara)

               if(abs(f_cql_cart(i_uperp, i_upara, i_psi)) .lt. 1.0e-99) &
     &                               f_cql_cart(i_uperp, i_upara, i_psi) = 0.0
               if(abs(df_cql_uprp(i_uperp, i_upara, i_psi)) .lt. 1.0e-99) &
     &                               df_cql_uprp(i_uperp, i_upara, i_psi) = 0.0
               if(abs(df_cql_uprl(i_uperp, i_upara, i_psi)) .lt. 1.0e-99) &
     &                               df_cql_uprl(i_uperp, i_upara, i_psi) = 0.0

            end do
         end do

      end do

!     ------------------------
!     integrate to get density
!     ------------------------

      do i_psi = 1, n_psi

         do i_uperp = 1, NUPER
            do i_upara = 1, NUPAR
               f(i_uperp, i_upara)= f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

          if (i_psi .eq. 1) then
             call a2dmnmx(f, NUPER, NUPAR, nuper, nupar, fmin, fmax)
          end if

         call  sgrate_cyl(uperp, upara, f, 1, nuper, 1, nupar, &
     &      dens(i_psi), NUPER, NUPAR)

      end do


!      write(6, *)
!      write(6, *) "uperp ="
!      write(6, 310) (uperp(i_uperp), i_uperp = 1, nuper)

!      write(6, *)
!      write(6, *) "upara ="
!      write(6, 310) (upara(i_upara), i_upara = 1, nupar)
      
      

      write(6, *)
      write(6, *) "density ="
      write(6, 310) (dens(i_psi), i_psi = 1, n_psi)
      
!      call exit(-1)

      write(6, *)
      write(6, *) "f_cql_cart_max = ", fmax

      write(6, *)
      write(6, *) "rho/a(i_psi)"
      write(6, 310) (rho_a(i_psi), i_psi = 1, n_psi)
      write(6, *)




!     ------------------------------------------------
!     Write data for AORSA2D in u_perp and u_parallel
!     ------------------------------------------------
      rewind (50)
      
      write (50, 309) nuper
      write (50, 309) nupar      
      write (50, 309) n_psi

      write (50, 3310) vc_mks_aorsa
      write (50, 3310) UminPara, UmaxPara

      write (50, 3310) (rho_a(i_psi), i_psi = 1, n_psi)
      write (50, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
      write (50, 3310) (upara(i_upara), i_upara = 1, nupar)

      write (50, 3310) (((f_cql_cart(i_uperp, i_upara, i_psi), &
     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

      write (50, 3310) (((df_cql_uprp(i_uperp, i_upara, i_psi), &
     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

      write (50, 3310) (((df_cql_uprl(i_uperp, i_upara, i_psi), &
     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)
     

      rewind (50)

!        read (50, 309) nuper
!        read (50, 309) nupar

!      read (50, 309) n_psi

!      read (50, 3310) vc_mks
!      read (50, 3310) UminPara, UmaxPara

!      read (50, 3310) (rho_a(i_psi), i_psi = 1, n_psi)
!      read (50, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
!      read (50, 3310) (upara(i_upara), i_upara = 1, nupar)

!      read (50, 3310) (((f_cql_cart(i_uperp, i_upara, i_psi), &
!     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

!      read (50, 3310) (((df_cql_uprp(i_uperp, i_upara, i_psi), &
!     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

!      read (50, 3310) (((df_cql_uprl(i_uperp, i_upara, i_psi), &
!     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)


!     ------------------------------------------------
!     get CQL3D distribution function on the midplane
!     ------------------------------------------------
!      call cql3d_dist(nupar, nuper, n_psi,  &
!     &              n_psi_dim, rho_a, rho,  &
!     &              UminPara,UmaxPara,      &
!     &              df_cql_uprp, df_cql_uprl,  &
!     &              UPERP, UPARA, DFDUPER, DFDUPAR)
     
!      do ni = 1, nuper
!          do mi = 1, nupar
!     
!             dfdth = upara(mi) * dfduper(ni, mi)  &
!     &                           - uperp(ni) * dfdupar(ni, mi)
!     
!            write(6,  1314)ni, mi, uperp(ni),  upara(mi), dfdth
!            write(115, 1314)ni, mi, uperp(ni),  upara(mi), dfdth
!                 
!         end do
!      end do
               
! 1314  format(2i10, 1p,9e12.4)  

!     ------------------------------------------------
!     Write data for plotting in u_perp and u_parallel
!     ------------------------------------------------

      write (238, 309) NUPER
      write (238, 309) NUPAR
      write (238, 309) n_psi

      write (238, 310) (uperp(i_uperp), i_uperp = 1, nuper)
      write (238, 310) (upara(i_upara), i_upara = 1, nupar)
      write (238, 310) (((f_cql_cart(i_uperp, i_upara, i_psi), &
     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)
     
      write(238, 310) (wperp_cql(i_psi, n_t), i_psi = 1, n_psi)
      write(238, 310) (wpar_cql(i_psi, n_t), i_psi = 1, n_psi)
     
      close (238)
 
 
!     ----------------------
!     take log of f_cql_cart
!     ----------------------
      do i_psi = 1, n_psi      
         do i_uperp = 1, nuper
            do i_upara = 1, nupar
               
               if(f_cql_cart(i_uperp, i_upara, i_psi) .le. 0.0) then
                  f_cql_cart(i_uperp, i_upara, i_psi) = -100.0
               end if
               
               if(f_cql_cart(i_uperp, i_upara, i_psi) .gt. 0.0) then
                  f_cql_cart(i_uperp, i_upara, i_psi) =  &
     &                    log10(f_cql_cart(i_uperp, i_upara, i_psi))
               end if

            end do
         end do
      end do 
     
      
!     ---------------
!     expand rho grid
!     ----------------
      do i_psi = 1, n_psi      
          rho_a(i_psi) = rho_a(i_psi) * 5.0      
      end do 
 
!     ---------------------------------------
!     Write f_cql_3D.vtk file "rectilinear_grid"
!     ---------------------------------------
      number_points = nupar * nuper * n_psi
      
      write(140, 2840)
 2840 format('# vtk DataFile Version 2.0')
 
      write(140, 2845)
 2845 format('Cql3d distribution function')      
      
      write(140, 2846)
 2846 format('ASCII')       
           
      write(140, 2847)
 2847 format('DATASET RECTILINEAR_GRID')  
      
      write(140, 2841)nupar, nuper, n_psi
 2841 format('DIMENSIONS', 3i8)
           
      write(140, 2842) nupar
 2842 format('X_COORDINATES', i8, ' float')
      write (140, 3411) (upara(i_upara), i_upara = 1, nupar) 
 
      write(140, 2843) nuper
 2843 format('Y_COORDINATES', i8, ' float')
      write (140, 3411) (uperp(i_uperp), i_uperp = 1, nuper) 
      
      write(140, 2851) n_psi
 2851 format('Z_COORDINATES', i8, ' float')
      write (140, 3411) (rho_a(i_psi), i_psi = 1, n_psi)      
      
      write(140, 2844) number_points
 2844 format('POINT_DATA', i10)      
      
      write(140, 2848)
 2848 format('SCALARS f_cql_cart float 1')       
      
      write(140, 2849)
 2849 format('LOOKUP_TABLE default')        
                  
      write (140, 3411) (((f_cql_cart(i_uperp, i_upara, i_psi), &
     &   i_upara = 1, nupar), i_uperp = 1, nuper), i_psi = 1, n_psi) 
     
     
 3410 format(1p,4e10.2)  
 3411 format(6f12.4)     
     
 
                 
!     ----------------
!     Write ahern file
!     ----------------          
!      write (143, 311) nuper
!      write (143, 311) nupar
!      write (143, 311) n_psi

!      write (143, 310) (uperp(i_uperp), i_uperp = 1, nuper)
!      write (143, 310) (upara(i_upara), i_upara = 1, nupar)
!      write (143, 310) (rho_a(i_psi), i_psi = 1, n_psi)
      
      
!      write (143, 310) (((f_cql_cart(i_uperp, i_upara, i_psi), &
!     &  i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)



!      go to 5000

!      n = 0
!      m = 0
!      i = 5
!      j = 5

!      zi = cmplx(0.0, 1.0)
!      eps0 = 8.85e-12
!      xmu0 = 1.26e-06
!      clight = 1.0 / sqrt(eps0 * xmu0)
!      pi = 3.141592654

!      omgrf = 2.0 * pi * xnurf
!      xk0 = omgrf / clight

!      gradprlb(i,j) = 0.0


!     -----------------
!     DIII-D parameters
!     -----------------

!      qi2 = 1.6e-19
!      xn2a(i,j) = 0.112920789750647245E+19
!      xnuomg = 0.0
!      xkti2(i,j) = xkt0


!      omgp22(i,j) =  xn2a(i,j) * qi2**2 / (eps0 * xmi)


!      lmin = -lmax
!      nbessj = lmax + 2

!      capr(i) =  1.56400000000000006
!      nphi = -15

!      xkxsav(n) =  0.000000000000000000E+00
!      xkysav(m) =  0.000000000000000000E+00
!      bxn(i,j) =  0.978610000000000070E-02
!      byn(i,j) =  -0.837220000000000047E-01
!      bzn(i,j) =  0.996439999999999992
!      bmod(i,j) =  2.01040000000000019
!      bmod_mid(i,j) = 1.62991393393751482


!      omgci2(i,j) = qi2 * bmod(i,j) / xmi

!     ---------------------------
!     Calculate rotation matrix U
!     ---------------------------

!      sqx = sqrt(1.0 - bxn(i,j)**2)

!      uxx(i, j) =   sqx
!      uxy(i, j) = - bxn(i, j) * byn(i, j) / sqx
!      uxz(i, j) = - bxn(i, j) * bzn(i, j) / sqx
!      uyx(i, j) =   0.0
!      uyy(i, j) =   bzn(i, j) / sqx
!      uyz(i, j) = - byn(i, j) / sqx
!      uzx(i, j) =   bxn(i, j)
!      uzy(i, j) =   byn(i, j)
!      uzz(i, j) =   bzn(i, j)

!      nzfun = 0

!      myid = 0
!      nproc = 1
!     -------------------------------------
!     setup Z function table when nzfun = 3
!     -------------------------------------
!      if (nzfun .eq. 3) then
!         iflag = 0
!         call ztable_setup(myid, iflag)
!         if(iflag .eq. 1)call ztable_setup(myid, nproc)
!      end if


!      ndisti2 = 0

!      call sigmad_cql3d(i, j, n, m, rho, rho_a,  &
!     &   gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
!     &   xmi, qi2, xn2a(i,j), xnuomg,  &
!     &   xkti2(i,j), omgci2(i,j), omgp22(i,j), &
!     &   -lmax, lmax, nzfun, ibessel,          &
!     &   xkxsav(n), xkysav(m), nphi, capr(i),  &
!     &   bxn(i,j), byn(i,j), bzn(i,j), &
!     &   uxx(i,j), uxy(i,j), uxz(i,j), &
!     &   uyx(i,j), uyy(i,j), uyz(i,j), &
!     &   uzx(i,j), uzy(i,j), uzz(i,j), &
!     &   sig2xx, sig2xy, sig2xz, &
!     &   sig2yx, sig2yy, sig2yz, &
!     &   sig2zx, sig2zy, sig2zz, &
!     &   delta0, ndisti2, nupar, nuper, n_psi, &
!     &   n_psi_dim, &
!     &   UminPara, UmaxPara, UPERP, UPARA, &
!     &   vc_mks, df_cql_uprp, df_cql_uprl, nbessj, ncount)

!       ndisti2 = 1

!      call sigmad_cql3d(i, j, n, m, rho, rho_a,  &
!     &   gradprlb(i,j), bmod(i,j), bmod_mid(i,j),  &
!     &   xmi, qi2, xn2a(i,j), xnuomg,  &
!     &   xkti2(i,j), omgci2(i,j), omgp22(i,j), &
!     &   -lmax, lmax, nzfun, ibessel,          &
!     &   xkxsav(n), xkysav(m), nphi, capr(i),  &
!     &   bxn(i,j), byn(i,j), bzn(i,j), &
!     &   uxx(i,j), uxy(i,j), uxz(i,j), &
!     &   uyx(i,j), uyy(i,j), uyz(i,j), &
!     &   uzx(i,j), uzy(i,j), uzz(i,j), &
!     &   sig2xx, sig2xy, sig2xz, &
!     &   sig2yx, sig2yy, sig2yz, &
!     &   sig2zx, sig2zy, sig2zz, &
!     &   delta0, ndisti2, nupar, nuper, n_psi, &
!     &   n_psi_dim, &
!     &   UminPara, UmaxPara, UPERP, UPARA, &
!     &   vc_mks, df_cql_uprp, df_cql_uprl, nbessj, ncount)


  310 format(1p,6e12.4)
 3310 format(1p,6e18.10)
  311 format(10i10)
  309 format(10i10)

  312 format(i10, 1p,6e12.4)
  
  
! 5000 continue


!     -----------------
!     deallocate arrays
!     -----------------

      deallocate( UPERP )
      deallocate( UPARA )
      deallocate( dfduper )
      deallocate( dfdupar )
      deallocate( f )

      deallocate( f_cql_cart )
      deallocate( df_cql_uprp )
      deallocate( df_cql_uprl )
      
      deallocate( dfdu_cql )
      deallocate( dfdtheta_cql )
      deallocate( f_aorsa )      
      
      close(50)
      close (140)


      return
      END subroutine cql3d_setup2

!
!***********************************************************************
!



       subroutine mets_grid(nupar, nuper, UminPara, UmaxPara, UmaxPerp, &
     &                       UPERP, UPARA, DFDUPER, DFDUPAR, f, fnorm, &
     &                       u_, theta_, f_, dfdu_, dfdth_, n_u_dim, n_theta_dim,  &
     &                       n_u, n_theta, vc_ratio)
       implicit none

       integer, intent(in) :: nupar, nuper
       integer, intent(in) :: n_u_dim, n_theta_dim
       real,    intent(in) :: u_(n_u_dim), theta_(n_theta_dim)
       real,    intent(in) :: f_(n_u_dim, n_theta_dim)
       real,    intent(in) :: dfdu_(n_u_dim, n_theta_dim)
       real,    intent(in) :: dfdth_(n_u_dim, n_theta_dim)
     
       integer, intent(in) :: n_u, n_theta

       real, intent(out) :: UminPara, UmaxPara, UmaxPerp
       real, intent(out) :: UPERP(NUPER)
       real, intent(out) :: UPARA(NUPAR)
       real, intent(out) :: DFDUPER(NUPER, NUPAR),  &
     &                      DFDUPAR(NUPER, NUPAR)
       real, intent(out) :: f(NUPER, NUPAR)
       real :: dfdu(NUPER, NUPAR)
       real :: dfdth(NUPER, NUPAR)

       real :: dfduprp, d2fduprp2
       real :: dfduprl, d2fduprl2

       real surf2
       real :: zp(n_u_dim, n_theta_dim, 3)
       real :: zpu(n_u_dim, n_theta_dim, 3)
       real :: zpth(n_u_dim, n_theta_dim, 3)

       real :: zx1(n_theta_dim), zxm(n_theta_dim)
       real :: zy1(n_u_dim), zyn(n_u_dim)
       real :: pi, fnorm, u, u2,u3, theta
       real :: u2_cql3d, u2_aorsa
       real :: ugiv, thegiv, costh, sinth, sigma
       integer n, m, islpsw, islpsw1, ierr
       real :: temp(2 *(n_u_dim + n_theta_dim) )
       real :: zxy11, zxym1, zxy1n, zxymn
       
       real ::  vc_ratio

!      ---------------------------------------------------
!      sigma = 0.0 for tensor product cubic splines
!      sigma = 50 for bi-linear interpolation
!      documentation recommends sigma=1.0 as standard value
!        ---------------------------------------------------
       sigma = 1.0
       islpsw = 255
       islpsw1 = 3


       do n = 1, NUPER
          UPERP(n) = UmaxPerp * (real(n-1) / real(NUPER-1))
       end do


!       UmaxPara =  0.6
       UminPara = -UmaxPara

       do m = 1, NUPAR
          UPARA(m) = UmaxPara * (-1.0 + 2. * (real(m-1)/real(NUPAR-1)))
       end do




       do n = 1, NUPER
          do m = 1, NUPAR
             u2_aorsa = UPERP(n)**2 + UPARA(m)**2
             
             u2_cql3d = u2_aorsa * vc_ratio**2
                     
             u2 = u2_cql3d


             ugiv  = sqrt(u2)
             thegiv = atan2(UPERP(n), UPARA(m))

             costh = cos(thegiv)
             sinth = sin(thegiv)

!            -----------------------
!            linear interpolation
!            -----------------------

             call intplt_2d(ugiv, thegiv, f(n, m), &
     &             n_u, n_theta, f_, n_u_dim, n_theta_dim, &
     &             u_, theta_)

             call intplt_2d(ugiv, thegiv, dfdu(n, m), &
     &             n_u, n_theta, dfdu_, n_u_dim, n_theta_dim, &
     &             u_, theta_)

             call intplt_2d(ugiv, thegiv, dfdth(n, m), &
     &             n_u, n_theta, dfdth_, n_u_dim, n_theta_dim, &
     &             u_, theta_)

!            -----------------------
!            quadradic interpolation
!            -----------------------

!             call WINTERP2D(u_, n_u, theta_, n_theta, f_, &
!     &                ugiv, thegiv, f(n,m) )

!             call WINTERP2D(u_, n_u, theta_, n_theta, dfdu_, &
!     &                ugiv, thegiv, dfdu(n,m) )

!             call WINTERP2D(u_, n_u, theta_, n_theta, dfdth_, &
!     &                ugiv, thegiv, dfdth(n,m) )


!            -----------------------
!            cubic spline interpolation
!            -----------------------

!             f(n, m) = surf2(ugiv, thegiv, n_u, n_theta, u_, theta_, f_, &
!     &             n_u_dim, zp, sigma)

!             dfdu(n, m) = surf2(ugiv, thegiv, n_u, n_theta, u_, theta_, dfdu_, &
!     &             n_u_dim, zpu, sigma)

!             dfdth(n, m) = surf2(ugiv, thegiv, n_u, n_theta, u_, theta_, dfdth_, &
!     &             n_u_dim, zpth, sigma)
    
                                             
             f(n,m) = f(n,m) * vc_ratio**3

             DFDUPER(n, m) = 0.0
             DFDUPAR(n, m) = 0.0

             if(ugiv .ne. 0)then
                DFDUPER(n, m) =  (costh / ugiv * dfdth(n, m) + sinth * dfdu(n, m)) * vc_ratio**4
                DFDUPAR(n, m) = (-sinth / ugiv * dfdth(n, m) + costh * dfdu(n, m)) * vc_ratio**4
             end if

             if(n .eq. 1) DFDUPER(n, m) = 0.0

          end do
       end do





  312 format(i10, 1p,6e12.4)

       return
       end subroutine mets_grid

!
!***********************************************************************
!
      subroutine intplt_2d(rgiv, thegiv, fl, &
     &                nr, nth, f, nrmax, nthmax, &
     &                r, theta)

      implicit none

      integer m, n, nr, nth, nrmax, nthmax, i, j, mp1, np1
      real f(nrmax, nthmax), r(nrmax), theta(nthmax)
      real zeta, c, d, b, eta, a, fl, rgiv, thegiv

      if (rgiv .gt. r(nr))then
         fl = 0.0
         return
      end if

      fl = 0.0

      do i = 1, nr - 1
         if(rgiv .ge. r(i) .and. rgiv .lt. r(i+1)) n = i
      end do
      if (rgiv .eq. r(nr)) n = nr



      do j = 1, nth - 1
         if(thegiv .ge. theta(j) .and. thegiv .lt. theta(j + 1)) m = j
      end do
      if(thegiv .ge. theta(nth)) m = nth



      if (n .gt. nr) then
         WRITE (6,*) "n .gt. nr", n
         !         call exit(-1)
         stop 1
      end if

      if (m .gt. nth) then
         WRITE (6,*) "m .gt. nth", m
         !         call exit(-1)
         stop 1
      end if


      if (m .lt. nth) mp1 = m + 1
      if (m .eq. nth) mp1 = m - 1

      if (n .lt. nr) np1 = n + 1
      if (n .eq. nr) np1 = n - 1


         zeta = (rgiv - r(n)) / (r(np1) - r(n))
         eta = (thegiv - theta(m)) / (theta(mp1) - theta(m))

         a = f(n, m)
         b = f(np1, m) - f(n, m)
         c = f(n, mp1) - f(n, m)
         d = f(np1, mp1) + f(n, m) - f(np1, m) - f(n, mp1)



      fl = a + b * zeta + c * eta + d * zeta * eta

      return
      end subroutine intplt_2d
!
!***********************************************************************
!

      subroutine cosft(f, theta, a, n_theta_max, n_theta, n_expand)

      implicit none

      integer, intent(in) :: n_theta_max, n_theta, n_expand

      real, intent(in) :: f(n_theta_max), theta(n_theta_max)

      real, intent(out) :: a(0: n_theta_max)

      real :: pi
      integer n, i

      data pi/3.141592654/


      do n = 0, n_expand
         a(n) = 0.0
         do i = 1, n_theta - 1
            a(n) = a(n) + 1.0 / pi * (f(i) * cos(n * theta(i)) &
     &                            + f(i+1) * cos(n * theta(i+1))) &
     &                            * (theta(i+1) - theta(i))
         end do
      end do

      a(0) = a(0) / 2.0

      return

      end subroutine cosft
!
!***********************************************************************
!

      subroutine cosft_inv(f, theta, a, n_theta_max, n_theta)

      implicit none

      integer, intent(in) :: n_theta_max, n_theta
      real, intent(in) :: a(0: n_theta_max), theta(n_theta_max)

      real, intent(out) :: f(n_theta_max)

      integer n, i


      do i = 1, n_theta
         f(i) = 0.0
         do n = 0, n_theta - 1
            f(i) = f(i) + a(n) * cos(n * theta(i))
         end do
      end do

      return

      end subroutine cosft_inv

!
!***********************************************************************
!

      subroutine intplt_1d(rgiv, fout, nr, fin, nrmax, r)

      implicit none

      integer n, nr, nrmax, i, np1
      real fin(nrmax), r(nrmax)
      real zeta, b,  a, fout, rgiv
      
      fout = 0.0

      if (rgiv .gt. r(nr))return
         
      do i = 1, nr - 1
         if(rgiv .ge. r(i) .and. rgiv .lt. r(i+1)) n = i
      end do
      if (rgiv .eq. r(nr)) n = nr

      if (n .gt. nr) then
         WRITE (6,*) "n .gt. nr", n
         !         call exit(-1)
         stop 1
      end if

      if (n .lt. nr) np1 = n + 1
      if (n .eq. nr) np1 = n - 1

      zeta = (rgiv - r(n)) / (r(np1) - r(n))

      a = fin(n)
      b = fin(np1) - fin(n)
         
      fout = a + b * zeta 

      return
      end subroutine intplt_1d
      
!
!***********************************************************************
!

      subroutine deriv_1(f, id, jd, i, j, imax, jmax, rg, dfdr, d2fdr2)

      implicit none

      integer id, jd, i, j, imax, jmax
      real f(id, jd), rg(id), dfdr, d2fdr2

      if(i .ne. 1 .and. i .ne. imax)then
         dfdr = (f(i+1, j) - f(i-1, j)) / (rg(i+1) - rg(i-1))
         d2fdr2 = (f(i+1, j) - 2.0 * f(i,j) + f(i-1, j)) / &
     &              (rg(i+1) - rg(i))**2
      end if

      if(i .eq. 1)then
!         dfdr = (f(i+1, j) - f(i, j)) / (rg(i+1) - rg(i))
         dfdr = 0.00
         d2fdr2 = (f(i+1, j) - 2.0 * f(i,j) + f(i+1, j)) / &
     &              (rg(i+1) - rg(i))**2
      end if

      if(i .eq. imax)then
!         dfdr = (f(i, j) - f(i-1, j)) /  (rg(i) - rg(i-1))
         dfdr = 0.00
         d2fdr2 = (f(i, j) - 2.0 * f(i-1,j) + f(i-2, j)) / &
     &              (rg(i) - rg(i-1))**2
      end if

      return
      end

!
!***************************************************************************
!
      subroutine deriv_2(f, id, jd, i, j, imax, jmax, zg, dfdz, d2fdz2)

      implicit none

      integer id, jd, i, j, imax, jmax
      real f(id, jd), zg(jd), dfdz, d2fdz2

      if(j .ne. 1 .and. j .ne. jmax)then
         dfdz = (f(i, j+1) - f(i, j-1)) / (zg(j+1) - zg(j-1))
         d2fdz2 = (f(i, j+1) - 2.0 * f(i,j) + f(i, j-1)) / &
     &      (zg(j+1) - zg(j))**2
      end if


      if(j .eq. 1)then
!         dfdz = (f(i, j+1) - f(i, j)) / (zg(j+1) - zg(j))
         dfdz = 0.0
         d2fdz2 = (f(i, j+2) - 2.0 * f(i,j+1) + f(i, j)) / &
     &      (zg(j+1) - zg(j))**2
      end if


      if(j .eq. jmax)then
!         dfdz = (f(i, j) - f(i, j-1)) /  (zg(j) - zg(j-1))
         dfdz = 0.0
         d2fdz2 = (f(i, j) - 2.0 * f(i,j-1) + f(i, j-2)) / &
     &      (zg(j) - zg(j-1))**2
      end if

  312 format(i10, 1p,8e12.4)

      return
      end

!
!***************************************************************************
!

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
!
!***************************************************************************
!

      subroutine sgrate_sphere(u, theta, f, nu1, nu2, ntheta1, ntheta2, &
     &   ans, nu_max, ntheta_max)

      implicit none

      integer i, j, nu1, nu2, ntheta1, ntheta2, nu_max, ntheta_max

      real f(nu_max, ntheta_max), u(nu_max), theta(ntheta_max), ans, du, dtheta
      real :: pi

      data pi/3.141592654/

      ans = 0.0

      do i = nu1, nu2 - 1
         du = u(i+1) - u(i)

         do j = ntheta1, ntheta2 - 1
            dtheta = theta(j+1) - theta(j)

            ans = ans + du * dtheta * (u(i)**2 * sin(theta(j))   * f(i, j)  &
     &                             + u(i+1)**2 * sin(theta(j))   * f(i+1, j) &
     &                               + u(i)**2 * sin(theta(j+1)) * f(i, j+1) &
     &                             + u(i+1)**2 * sin(theta(j+1)) * f(i+1, j+1) ) / 4.0
         end do

      end do

      ans = ans * 2.0 * pi

      return
      end
!
!***************************************************************************
!

!----------------------------------------------------------------------------!
! Subroutine a2dmnmx
!----------------------------------------------------------------------------!
! Minimum and the maximum elements of a 2-d array.
!----------------------------------------------------------------------------!

      subroutine a2dmnmx(f, nxmax, nymax, nx, ny, fmin, fmax)

      implicit none

      integer nx, ny, i, j, nxmax, nymax
      real f(nxmax, nymax), fmin, fmax

      fmax=f(2, 1)
      fmin=fmax

      do 23000 j=1, ny
!         do 23002 i=1, nx

         do 23002 i=2, nx-1
            fmax=amax1(fmax, f(i, j))
            fmin=amin1(fmin, f(i, j))
23002    continue

23000 continue

      return
 2201 format(2i5,1p,8e12.4)
      end
!
!***************************************************************************
!





