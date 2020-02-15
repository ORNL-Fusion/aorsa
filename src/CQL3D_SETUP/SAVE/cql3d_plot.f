      program cql3d_plot

      implicit none


      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titx
      character*32 tity
      character*32 titz
      real logmax, ycut, dy
      integer nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, nlevmax
      integer ibackground, nkx1, nkx2, nky1, nky2, n, m,
     .   nkpltdim, mkpltdim, nkxplt, nkyplt, ipage

      integer i_psi1, i_psi2, i_psi3, i_psi4

      integer  ncolln10, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, ncollin, ncollab, ncolion,
     1   ncolelec

      integer nrhomax, nnoderho, iflag
      integer n_theta_max, n_u_max, n_psi_max

      real vc_cgs, duminpara, dumaxpara
      integer dnuper, dnupar, dnnoderho, n1, n2, n3, n4

      parameter (nlevmax = 101)


      parameter (n_theta_max = 150)
      parameter (n_u_max = 150)
      parameter (n_psi_max = 150) 

      integer, parameter :: LNUPER_dim = 129
      integer, parameter :: LNUPAR_dim = 257

      real :: UPERP(LNUPER_dim),UPARA(LNUPAR_dim)
      real :: DFDUPER(LNUPER_dim, LNUPAR_dim),
     .        DFDUPAR(LNUPER_dim, LNUPAR_dim)

      real :: f(LNUPER_dim, LNUPAR_dim)
      real :: dUPERP(LNUPER_dim), dUPARA(LNUPAR_dim) 


      real :: drhon(n_psi_max)  

      real  y1max_giv

      real :: f_cql_cart(LNUPER_dim, LNUPAR_dim, n_psi_max)
      real :: df_cql_uprp(LNUPER_dim, LNUPAR_dim, n_psi_max)
      real :: df_cql_uprl(LNUPER_dim, LNUPAR_dim, n_psi_max)
      real ::   bqlavg_i1(LNUPER_dim, LNUPAR_dim, n_psi_max) 

      real :: f_cql_cart_2d(LNUPAR_dim, LNUPER_dim)
      real ::  bqlavg_i1_2d(LNUPAR_dim, LNUPER_dim) 

      integer :: i_uperp, i_upara, lnuper, lnupar

      real u(n_u_max), theta(n_theta_max)
      real f_cql_2d(n_u_max, n_theta_max)
      real f_cql_1d(n_u_max)

	real xg(n_u_max, n_theta_max), yg(n_u_max, n_theta_max)

      real f_cql(n_theta_max, n_u_max, n_psi_max)
      real dfdu_cql(n_theta_max, n_u_max, n_psi_max)
      real dfdtheta_cql(n_theta_max, n_u_max, n_psi_max)


      integer n_theta_(n_psi_max)
      real theta_(n_theta_max, n_psi_max)



      integer n_theta, n_u, n_psi
      integer i_theta, i_u, i_psi
      real vc, r0


      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1   ncolion,ncolelec,ncollin,ncollab,ncolln2
      common/boundcom/rhoplasm
      common/zoom/ xmaxz, xminz, ymaxz

      real xmaxz, xminz, ymaxz


      real fmin, fmax, fminre, fmaxre, fminim, fmaxim, fmin1,
     .   fmax1, fmax2, fmax3, fmaxs, fmaxt, fmin2, fmin3, fmins,fmint


      real ff(nlevmax)
      real q, omgrf, xk0, n0, clight, xmu0, eps0, rhoplasm
      real temax, temin, timin, tmin, tmax, timax, caprmaxp,
     1   caprmin, caprminp, caprmax, xnmax, xnmin, qmin, qmax
      real bpmin, bpmax
      integer nnodex, j, i, nnodey, numb, jmid, jcut
      complex zi

      namelist/cql3d_plotin/ibackground, xminz, xmaxz, ymaxz, logmax,
     .   ipage, numb, ycut, y1max_giv, i_psi1, i_psi2, i_psi3, i_psi4

      open(unit=63, file='cql3d_plot.in', status='old',
     .   form='formatted')
      rewind(63)


!     --------------------------------
!     set default values of input data:
!     --------------------------------
      ibackground = 1
      ipage = 1
      numb = 20
      xminz = 0.51
      xmaxz = 0.62
      ymaxz = 0.24
      logmax = 2.0
      ycut = 0.0
      y1max_giv = 0.0

      i_psi1 = 3
      i_psi2 = 5
      i_psi3 = 9
      i_psi4 = 15



!     ----------------
!     read input data:
!     ----------------
      read (63, cql3d_plotin)

      open(unit=138,file='out138',status='old',form='formatted')
    

      zi = cmplx(0.0, 1.0)
      eps0 = 8.85e-12
      xmu0 = 1.26e-06
      clight = 1./sqrt(eps0 * xmu0)
      xk0 = omgrf / clight
      q = 1.6e-19

      r0 = 0.0


      nblack=0
      nred=1
      nyellow=2
      ngreen=3
      naqua=4
      npink=5
      nwheat=6
      ngrey=7
      nbrown=8
      nblue=9
      nblueviolet=10
      ncyan=11
      nturquoise=12
      nmagenta=13
      nsalmon=14
      nwhite=15

      ncollin=ncyan
      ncolln2=nblueviolet
      ncolln3=ngreen
      ncolln4=naqua
      ncolln5=nyellow
      ncolbrd=nwheat
      ncollab=ncyan

      ncolbox=nred
      ncolbrd=nwheat
      ncolion=naqua
      ncolelec=nyellow

      if (ibackground.eq.1)then
         ncolbox=nred
         ncolbrd=nwheat
      end if

      if (ibackground.eq.0)then
         ncolbox=nblack
         ncolbrd=nblack
      end if

      call plsfnam('cql3d_plot.plp')
      if (ipage .eq. 2)call plstart('plmeta',2, 2)
      if (ipage .eq. 1)call plstart('plmeta',1, 1)


!     -----------------------
!	read data for plotting
!     -----------------------

	read (138, 311) n_u
	read (138, 311) n_psi
	read (138, 310) vc

      read (138, 310) (u(i_u), i_u = 1, n_u)
	read (138, 311) (n_theta_(i_psi), i_psi = 1, n_psi)
      read (138, 310) ((theta_(i_theta, i_psi),
     .                        i_theta = 1, n_theta_(i_psi)),
     .                        i_psi = 1, n_psi)

      read (138, 310) (((f_cql(i_theta, i_u, i_psi),
     .                         i_theta = 1, n_theta_(i_psi)),
     .                         i_u = 1, n_u),
     .                         i_psi = 1, n_psi)

      read (138, 310) (((dfdu_cql(i_theta, i_u, i_psi),
     .                         i_theta = 1, n_theta_(i_psi)),
     .                         i_u = 1, n_u),
     .                         i_psi = 1, n_psi)
      read (138, 310) (((dfdtheta_cql(i_theta, i_u, i_psi),
     .                         i_theta = 1, n_theta_(i_psi)),
     .                         i_u = 1, n_u),
     .                         i_psi = 1, n_psi)

      titx = 'u_parallel'
      tity = 'u_perp'



!     -----------------------
!	2D plots of f(u, theta)
!     -----------------------
      title = 'f(u, theta)'

      i_psi = i_psi1

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = f_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)



      i_psi = i_psi2

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = f_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)



      i_psi = i_psi3

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = f_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)



      i_psi = i_psi4

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = f_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)



!     -----------------------------
!	1D plots of f(u, theta = 0)
!     -----------------------------
      open(unit=36,file='Swain',status='unknown',form='formatted')
      
      title = 'log f(u, theta = 90 deg.)'
      titll= 'log f(u)'
      titlr='       '

      i_psi = i_psi1
      i_theta = n_theta_(i_psi) / 2.0

      do i_u = 1, n_u
         f_cql_1d(i_u) = -40.0
	 
         if (f_cql(i_theta, i_u, i_psi) .gt. 0.0) then
	    f_cql_1d(i_u) = alog10(f_cql(i_theta, i_u, i_psi))
	 end if
	 
	 write(6, 312)i_u, u(i_u), f_cql(i_theta, i_u, i_psi), 
     .      f_cql_1d(i_u)
      end do
      call ezplot1(title, titll, titlr, u, f_cql_1d,
     .    n_u, n_u_max, y1max_giv)


      i_psi = i_psi2
      i_theta = n_theta_(i_psi) / 2.0
      
      do i_u = 1, n_u
         f_cql_1d(i_u) = -40.0
         if (f_cql(i_theta, i_u, i_psi) .gt. 0.0) then
	    f_cql_1d(i_u) = alog10(f_cql(i_theta, i_u, i_psi))
	 end if
      end do
      call ezplot1(title, titll, titlr, u, f_cql_1d,
     .    n_u, n_u_max, y1max_giv)



      i_psi = i_psi3
      i_theta = n_theta_(i_psi) / 2.0

      do i_u = 1, n_u
         f_cql_1d(i_u) = -40.0
         if (f_cql(i_theta, i_u, i_psi) .gt. 0.0) then
	    f_cql_1d(i_u) = alog10(f_cql(i_theta, i_u, i_psi))
	 end if
      end do
      call ezplot1(title, titll, titlr, u, f_cql_1d,
     .    n_u, n_u_max, y1max_giv)



      i_psi = i_psi4
      i_theta = n_theta_(i_psi) / 2.0

      do i_u = 1, n_u
         f_cql_1d(i_u) = -40.0
         if (f_cql(i_theta, i_u, i_psi) .gt. 0.0) then
	    f_cql_1d(i_u) = alog10(f_cql(i_theta, i_u, i_psi))
	 end if
      end do
      call ezplot1(title, titll, titlr, u, f_cql_1d,
     .    n_u, n_u_max, y1max_giv)




!      i_psi = 18
!      do i_theta = 1, n_theta_(i_psi)
!         theta(i_theta) = theta_(i_theta, i_psi)
!         do i_u = 1, n_u
!            f_cql_2d(i_u, i_theta) = f_cql(i_theta, i_u, i_psi)
!         end do
!      end do
!      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
!     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)


      title = 'dfdu(u, theta)'

      i_psi = i_psi1

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = dfdu_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)


      i_psi = i_psi2

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = dfdu_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)

      i_psi = i_psi3

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = dfdu_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)

      i_psi = i_psi4

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = dfdu_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)





      title = 'dfdtheta(u, theta)'

      i_psi = i_psi1

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = dfdtheta_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)


      i_psi = i_psi2

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = dfdtheta_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)

      i_psi = i_psi3

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = dfdtheta_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)

      i_psi = i_psi4

      do i_theta = 1, n_theta_(i_psi)
         theta(i_theta) = theta_(i_theta, i_psi)
         do i_u = 1, n_u
            f_cql_2d(i_u, i_theta) = dfdtheta_cql(i_theta, i_u, i_psi)
         end do
      end do

      call ezconpx(u, theta, f_cql_2d, ff, n_u, n_theta_(i_psi), numb,
     .   n_u_max, n_theta_max, nlevmax, xg, yg, title, r0)



!     --------------------------------------------
!	read data for plotting f(u_perp, u_parallel)
!     --------------------------------------------

	read (138, 311) lnuper
	read (138, 311) lnupar

      read (138, 310) (uperp(i_uperp), i_uperp = 1, lnuper)
      read (138, 310)  (upara(i_upara), i_upara = 1, lnupar)
      read (138, 310) (((f_cql_cart(i_uperp, i_upara, i_psi),
     .  i_uperp = 1, lnuper), i_upara = 1, lnupar), i_psi = 1, n_psi)

      read (138, 310) (((df_cql_uprp(i_uperp, i_upara, i_psi),
     .  i_uperp = 1, lnuper), i_upara = 1, lnupar), i_psi = 1, n_psi)

      read (138, 310) (((df_cql_uprl(i_uperp, i_upara, i_psi),
     .  i_uperp = 1, lnuper), i_upara = 1, lnupar), i_psi = 1, n_psi)

      read (138, 310) ((DFDUPER(i_uperp, i_upara),
     &  i_uperp = 1, lnuper), i_upara = 1, lnupar)

      read (138, 310) ((DFDUPAR(i_uperp, i_upara),
     &  i_uperp = 1, lnuper), i_upara = 1, lnupar)


      title = 'f(u_perp, u_parallel)'
      i_psi = i_psi1
      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             f_cql_cart(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)
     


 3410 format(1p4e10.2)      
 
 



      i_psi = i_psi2

      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             f_cql_cart(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)

      i_psi = i_psi3

      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             f_cql_cart(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)


      i_psi = i_psi4

      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             f_cql_cart(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)
     
      write(6, *)"i_psi1 = ", i_psi1
      write(6, *)"i_psi2 = ", i_psi2
      write(6, *)"i_psi3 = ", i_psi3
      write(6, *)"i_psi4 = ", i_psi4

      go to 5004

*     -----------------------------------------------------------
*     Read quasilinear diffusion coefficients from file cql3d.coef
*     -----------------------------------------------------------
      open(unit=42,file='cql3d.coef',  status='unknown',
     .                                             form='formatted')

      read (42, 309) dnuper
      read (42, 309) dnupar
      read (42, 309) dnnoderho    
      

      read (42, 3310) vc_cgs
      read (42, 3310) dUminPara, dUmaxPara

      read (42, 3310) (drhon(n), n = 1, dnnoderho)
      read (42, 3310) (duperp(i_uperp), i_uperp = 1, dnuper)
      read (42, 3310) (dupara(i_upara), i_upara = 1, dnupar)

      read (42, 3310) (((bqlavg_i1(i_uperp, i_upara, n),
     &   i_uperp = 1, dnuper), i_upara = 1, dnupar), n = 1, dnnoderho)
     
     
      close (42)
      
!     --------------------------------------
!     2D plots of bqlavg(u_perp, u_parallel)
!     --------------------------------------


      
      write(6, *)"bqlavg_i1(25, 50, 8) = ", bqlavg_i1(25, 50, 8)

      numb = 23
      
      title = 'bqlavg(u_perp, u_parallel)'
      i_psi = i_psi1
      
      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            bqlavg_i1_2d(i_upara, i_uperp) =
     .                             bqlavg_i1(i_uperp, i_upara, i_psi)
         end do
      end do
      
      write(6, *)"bqlavg_i1_2d(50, 25) = ", bqlavg_i1_2d(50, 25)
      write(6, *)"upara(50) = ", upara(50)
      write(6, *)"uperp(25) = ", uperp(25)
      write(6, *)"numb = ", numb
      write(6, *)"nlevmax = ", nlevmax
      write(6, *)"i_psi = ", i_psi
      write(6, *)"lnupar = ", lnupar
      write(6, *)"lnuper = ", lnuper
      
      call ezconcx(upara, uperp, bqlavg_i1_2d, ff, lnupar, lnuper, numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)
     
     
      i_psi = i_psi2

      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            bqlavg_i1_2d(i_upara, i_uperp) =
     .                             bqlavg_i1(i_uperp, i_upara, i_psi)
         end do
      end do
      
      call ezconcx(upara, uperp, bqlavg_i1_2d, ff, lnupar, lnuper, numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)
     
     
      i_psi = i_psi3

      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            bqlavg_i1_2d(i_upara, i_uperp) =
     .                             bqlavg_i1(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, bqlavg_i1_2d, ff, lnupar, lnuper, numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)
     
     
      i_psi = i_psi4

      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            bqlavg_i1_2d(i_upara, i_uperp) =
     .                             bqlavg_i1(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, bqlavg_i1_2d, ff, lnupar, lnuper, numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)
     
      numb = 20

      
 5004 continue
 
 
!      i_psi = 18
!      do  i_upara = 1, lnupar
!         do i_uperp = 1, lnuper
!            f_cql_cart_2d(i_upara, i_uperp) =
!     .                             f_cql_cart(i_uperp, i_upara, i_psi)
!         end do
!      end do
!      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
!     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)


      title = 'df/duperp'
      i_psi = i_psi1
      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             df_cql_uprp(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)


      i_psi = i_psi2
      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             df_cql_uprp(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)


      i_psi = i_psi3
      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             df_cql_uprp(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)


      i_psi = i_psi4
      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             df_cql_uprp(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)


!      i_psi = 18
!      do  i_upara = 1, lnupar
!         do i_uperp = 1, lnuper
!            f_cql_cart_2d(i_upara, i_uperp) =
!     .                             df_cql_uprp(i_uperp, i_upara, i_psi)
!         end do
!      end do
!      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper, numb,
!     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)




      title = 'df/du_parallel'
      i_psi = i_psi1
!      do  i_upara = 1, lnupar
!         do i_uperp = 1, lnuper
!            f_cql_cart_2d(i_upara, i_uperp) =
!     .                             df_cql_uprl(i_uperp, i_upara, i_psi)
!         end do
!      end do


      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             df_cql_uprl(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)


      i_psi = i_psi2
      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             df_cql_uprl(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,
     .   numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)

      i_psi = i_psi3
      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             df_cql_uprl(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)

      i_psi = i_psi4
      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             df_cql_uprl(i_uperp, i_upara, i_psi)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)



!     ------------------------------------------
!     plot derivatives at rho given in data file
!     ------------------------------------------

      title = 'df/duperp'

      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             dfduper(i_uperp, i_upara)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)

      title = 'df/duparallel'


      do  i_upara = 1, lnupar
         do i_uperp = 1, lnuper
            f_cql_cart_2d(i_upara, i_uperp) =
     .                             dfdupar(i_uperp, i_upara)
         end do
      end do
      call ezconcx(upara, uperp, f_cql_cart_2d, ff, lnupar, lnuper,numb,
     .   LNUPAR_dim, LNUPER_dim, nlevmax, title, titx, tity, iflag)


      call plend


  310 format(1p6e12.4)
  311 format(10i10)
  312 format(i10, 1p6e12.4)
  309 format(10i10)
 3310 format(1p6e18.10)

 1312 format(i10,1p8e12.4)
 1313 format(2i10,1p8e12.4)
 9310 format(1p7e12.4)
   10 format(i10,1p4e10.3,i10,1pe10.3)
      stop
      end
c
c********************************************************************
c

c Subroutine a1mnmx
c----------------------------------------------------------------------------!
c Minimum and the maximum elements of a 1-d array.
c----------------------------------------------------------------------------!

      subroutine a1mnmx(f, nxmax, nx, fmin, fmax)

      implicit none

      integer nx, i, nxmax
      real f(nxmax), fmin, fmax

      fmax=f(1)
      fmin=fmax
         do 23002 i=1, nx
            fmax=amax1(fmax, f(i))
            fmin=amin1(fmin, f(i))
23002    continue

      return
 2201 format(2i5,1p8e12.4)
      end

c
c********************************************************************
c
c Subroutine a2mnmx
c----------------------------------------------------------------------------!
c Minimum and the maximum elements of a 1-d array between x1 and x2
c----------------------------------------------------------------------------!

      subroutine a2mnmx(f, nxmax, nx, x, x1, x2, fmin, fmax)

      implicit none

      integer nx, i, nxmax
      real f(nxmax), fmin, fmax, x1, x2
      real x(nxmax)

      fmax=0.0
      fmin=fmax
      do i = 1, nx
         if(x(i).gt.x1.and.x(i).lt.x2)then
            fmax=amax1(fmax, f(i))
            fmin=amin1(fmin, f(i))
         endif
      end do

      return

      end

c
c********************************************************************
c
c----------------------------------------------------------------------------!
c Subroutine a2dmnmx
c----------------------------------------------------------------------------!
c Minimum and the maximum elements of a 2-d array.
c----------------------------------------------------------------------------!

      subroutine a2dmnmx(f, nxmax, nymax, nx, ny, fmin, fmax)

      implicit none

      integer nx, ny, i, j, nxmax, nymax
      real f(nxmax, nymax), fmin, fmax

      fmax=f(2, 1)
      fmin=fmax
      do 23000 j=1, ny
c         do 23002 i=1, nx
         do 23002 i=2, nx-1
            fmax=amax1(fmax, f(i, j))
            fmin=amin1(fmin, f(i, j))
23002    continue
23000 continue

      return
 2201 format(2i5,1p8e12.4)
      end
c
c********************************************************************
c

      subroutine ezconp(r, theta, f, flevel, nr, nth, nlevel,
     1   nrmax, nthmax, nlevmax, xg, yg, title, r0)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1   ncolelec, j, n, m, nth1
      real xmaxz, xminz, ymaxz, r0

      real r(nrmax),theta(nthmax)
      real f(nrmax,nthmax)
      real xg(nrmax,nthmax)
      real yg(nrmax,nthmax)
      real flevel(nlevmax)

      character*8 xopt,yopt
      character*32 title


      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1   ncolion,ncolelec,ncollin,ncollab,ncolln2



c--set up transformation
      do n = 1, nr
         do m = 1, nth
            xg(n,m) = r(n) * cos(theta(m)) + r0
            yg(n,m) = r (n) * sin(theta(m))
         end do
         xg(n, nth+1) = xg(n,1)
         yg(n, nth+1) = yg(n,1)
      end do
      nth1 = nth + 1
      call a2dmnmx(xg, nrmax, nthmax, nr, nth, xmin, xmax)
      call a2dmnmx(yg, nrmax, nthmax, nr, nth, ymin, ymax)



c--set up contour levels
      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)
      df=abs(fmax-fmin)/float(nlevel)

c        write (6, *) "df = ", df

      if(df .eq. 0.0)return

      do i = 1, nlevel
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel.eq.1)flevel(1)=1.0e-03
      if(nlevel.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif


c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0

      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do

  100 continue

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt



c--Advance graphics frame and get ready to plot
      call plcol(ncolbox)
      call plenv(xmin, xmax, ymin, ymax, 1, 0)
c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps 
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps


      call pllsty(1)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then
         call plwid(1)
         call pllsty(2)
         call plcol(npink)
         call plcon2(f,nrmax,nthmax,1,nr,1,nth1,flevel(ilevlt),
     &   nlevlt, xg, yg)
      endif

      if(nlevgt .gt. 0) then
         call plwid(1)
         call pllsty(1)
         call plcol(ncollin)
         call plcon2(f,nrmax,nthmax,1,nr,1,nth1,flevel(ilevgt),
     &   nlevgt, xg, yg)
      endif
      call plwid(1)


      call plcol(ncollab)
      call pllab('u_parallel','u_perp', title)


  310 format(1p6e12.4)
  311 format(10i10)

      return
      end
c
c*********************************************************************
c

      subroutine ezconpx(r, theta, f, flevel, nr, nth, nlevel0,
     1   nrmax, nthmax, nlevmax, xg, yg, title, r0)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1   ncolelec, j, n, m, nth1, nlevel0
      real xmaxz, xminz, ymaxz, r0, fact

      real r(nrmax),theta(nthmax)
      real f(nrmax,nthmax)
      real xg(nrmax,nthmax)
      real yg(nrmax,nthmax)
      real flevel(nlevmax)

      character*8 xopt,yopt
      character*32 title


      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1   ncolion,ncolelec,ncollin,ncollab,ncolln2



c--set up transformation
      do n = 1, nr
         do m = 1, nth
            xg(n,m) = r(n) * cos(theta(m)) + r0
            yg(n,m) = r (n) * sin(theta(m))
         end do
         xg(n, nth+1) = xg(n,1)
         yg(n, nth+1) = yg(n,1)
      end do
      nth1 = nth + 1
      call a2dmnmx(xg, nrmax, nthmax, nr, nth, xmin, xmax)
      call a2dmnmx(yg, nrmax, nthmax, nr, nth, ymin, ymax)



c--set up contour levels
      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)
      df=abs(fmax-fmin)/float(nlevel0)

c        write (6, *) "df = ", df

      if(df .eq. 0.0)return

      do i = 1, nlevel0
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel0.eq.1)flevel(1)=1.0e-03
      if(nlevel0.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif

      if (nlevel0 .eq. 20) then

         nlevel = nlevel0 * 2.0

         fact = 2.0
         i = 1
         flevel(i) = fmin

         do i = 2, nlevel / 2
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do


         fact = 2.0
         i = nlevel / 2 + 1
         flevel(i) = fmax

         do i = nlevel / 2 + 2, nlevel
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do

      else
         nlevel = nlevel0
      end if


c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0

      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do

  100 continue

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt



c--Advance graphics frame and get ready to plot
      call plcol(ncolbox)
      call plenv(xmin, xmax, ymin, ymax, 0, 0)
c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps


      call pllsty(1)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then
         call plwid(1)
         call pllsty(2)
c         call plcol(ncolln2)
         call plcol(nyellow)
         call plcon2(f,nrmax,nthmax,1,nr,1,nth1,flevel(ilevlt),
     &   nlevlt, xg, yg)
      endif

      if(nlevgt .gt. 0) then
         call plwid(1)
         call pllsty(1)
c         call plcol(ncollin)
         call plcol(nblue)
         call plcon2(f,nrmax,nthmax,1,nr,1,nth1,flevel(ilevgt),
     &   nlevgt, xg, yg)
      endif
      call plwid(1)


c      call plcol(ncollab)
      call plcol(nblue)
      call pllab('u_parallel','u_perp', title)


  310 format(1p6e12.4)
  311 format(10i10)

      return
      end
c
c*********************************************************************
c
      subroutine ezconc(r, theta, f, flevel, nr, nth, nlevel,
     1   nrmax, nthmax, nlevmax, title, titx, tity, iflag)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, theta, r, flevel, f, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1    ncolelec, iflag, j

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1    nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1    nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1    ncolion,ncolelec,ncollin,ncollab,ncolln2



  310 format(1p6e12.4)
  312 format(i10, 1p6e12.4)


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

c--set up contour levels

      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)


      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
         iflag = 1
         return
      end if



      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
c        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel.eq.1)flevel(1)=1.0e-03

      if(nlevel.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif



c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do
  100 continue
      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


c--Advance graphics frame and get ready to plot
c     call pladv(0)
      call plcol(ncolbox)
      call plenv(xmin, xmax, ymin, ymax, 1, 0)
c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps
         call pllsty(1)
c      call plvpas(0.15,0.85,0.15,0.85,0.0)
c      call plwind(xpmin,xpmax,ypmin,ypmax)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0
c     call plbox(xopt,xtick,nxsub,yopt,ytick,nysub)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

c         call plwid(0.1)
         call pllsty(4)
         call plcol(ncolln2)
c         call plcol(nyellow)

         call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &      nlevlt, r, theta)
      endif

      if(nlevgt .gt. 0) then

c         call plwid(10)
         call pllsty(1)
         call plcol(ncollin)
c         call plcol(nblue)

         call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &      nlevgt, r, theta)
      endif

      call pllsty(1)
      call plcol(ncollab)
      call pllab(titx,tity,title)
      return
      end

c
c*********************************************************************
c
      subroutine ezconcx(r, theta, f, flevel, nr, nth, nlevel0,
     1   nrmax, nthmax, nlevmax, title, titx, tity, iflag)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, theta, r, flevel, f, xmin, ymin, xmax, fact
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1    ncolelec, iflag, j, nlevel0

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1    nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1    nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1    ncolion,ncolelec,ncollin,ncollab,ncolln2



  310 format(1p6e12.4)
  312 format(i10, 1p6e12.4)


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

c--set up contour levels

      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)


      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
         iflag = 1
         return
      end if



      df = abs(fmax - fmin) / float(nlevel0)
      do i = 1, nlevel0
c        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel0.eq.1)flevel(1)=1.0e-03

      if(nlevel0.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif

      if (nlevel0 .eq. 20) then

         nlevel = nlevel0 * 2.0

         fact = 2.0
         i = 1
         flevel(i) = fmin

         do i = 2, nlevel / 2
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do



         fact = 2.0
         i = nlevel / 2 + 1
         flevel(i) = fmax

         do i = nlevel / 2 + 2, nlevel
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do

      else

         nlevel = nlevel0

      end if


c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do
  100 continue
      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


c--Advance graphics frame and get ready to plot
c     call pladv(0)
      call plcol(ncolbox)
      call plenv(xmin, xmax, ymin, ymax, 0, 0)
c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps
         call pllsty(1)
c      call plvpas(0.15,0.85,0.15,0.85,0.0)
c      call plwind(xpmin,xpmax,ypmin,ypmax)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0
c     call plbox(xopt,xtick,nxsub,yopt,ytick,nysub)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

         call plwid(1)
         call pllsty(2)
         call plcol(nyellow)

         call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &      nlevlt, r, theta)
      endif

      if(nlevgt .gt. 0) then

         call plwid(1)
         call pllsty(1)
         call plcol(nblue)

         call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &      nlevgt, r, theta)
      endif

      call pllsty(1)

      call plcol(nblue)
      call pllab(titx,tity,title)
      
      
      return
      end
c
c*********************************************************************
c
      subroutine ezconz(r, theta, f, flevel, nr, nth, nlevel,
     1   nrmax, nthmax, nlevmax, title, titx, tity)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, theta, r, flevel, f, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1   ncolelec, j
      real xmaxz, xminz, ymaxz

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1    nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1    nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1    ncolion,ncolelec,ncollin,ncollab,ncolln2
      common/zoom/ xmaxz, xminz, ymaxz


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

c--set up contour levels

c      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)


      fmax=f(2, 1)
      fmin=fmax
      do j = 1, nth
         do i = 2, nr - 1
            if (r(i) .gt. xminz .and. r(i) .le. xmaxz .and.
     .         theta(j) .gt. -ymaxz .and. theta(j) .le. ymaxz) then
               fmax = amax1(fmax, f(i, j))
               fmin = amin1(fmin, f(i, j))
            end if
         end do
      end do

      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)return


      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
c        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel .eq. 1)flevel(1) = 1.0e-03

      if(nlevel.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif

c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do
  100 continue
      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


c--Advance graphics frame and get ready to plot
c     call pladv(0)
      call plcol(ncolbox)

      xmin = xminz
      xmax = xmaxz
      ymax = ymaxz
      ymin = -ymax

      call plenv(xmin, xmax, ymin, ymax, 1, 0)
c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps
         call pllsty(1)
c      call plvpas(0.15,0.85,0.15,0.85,0.0)
c      call plwind(xpmin,xpmax,ypmin,ypmax)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0
c     call plbox(xopt,xtick,nxsub,yopt,ytick,nysub)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

c         call plwid(0.1)
         call pllsty(4)
         call plcol(ncolln2)
c         call plcol(nyellow)

         call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &      nlevlt, r, theta)
      endif

      if(nlevgt .gt. 0) then

c         call plwid(10)
         call pllsty(1)
         call plcol(ncollin)
c         call plcol(nblue)

         call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &      nlevgt, r, theta)
      endif

      call plcol(ncollab)
      call pllab(titx,tity,title)
      return
      end
c
c*********************************************************************
c
      subroutine ezcon3d(r,theta,f,flevel,nr,nth,nlevel,
     1   nrmax,nthmax,nlevmax,title,titx,tity,titz)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1    ncolelec

      real r(nrmax), theta(nthmax)
      real f(nrmax, nthmax)
      real flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity
      character*32 titz

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1    nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1    nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1    ncolion,ncolelec,ncollin,ncollab,ncolln2

      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

c--set up contour levels
      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)

      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)return
      df = abs(fmax - fmin) / float(nlevel)


  100 format(1p8e12.4)


c--Advance graphics frame and get ready to plot
         call plcol(ncolbox)
         call plenv(-2.5, 2.5, -2.5, 4.0, 0, -2)
c         call plenv(-2.3, 1.8, -1.5, 3.5, 0, -2)
         call plw3d(2.,4.,3.,xmin,xmax,ymin,ymax,fmin,fmax,30.,30.)
         call plschr(0.,.50)
         call plbox3('binstu',titx,0.0,0,'binstu',
     *     tity,0.0,0,'binstu',titz,0.0,0)
         call plcol(ngreen)
         call plmesh(r,theta,f,nr,nth,3,nrmax)
         call plschr(0.,1.0)
         call plcol(ncollab)
         call plmtex('t',1.0,0.5,0.5,title)
      return
      end

c
c********************************************************************
c
      subroutine boundary(r, theta, f, flevel, nr, nth, nlevel,
     1   nrmax, nthmax, nlevmax, title, titx, tity)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, theta, r, flevel, f, xmin, ymin, xmax, rhoplasm
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1    ncolelec, nlevelb

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1    nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1    nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1    ncolion,ncolelec,ncollin,ncollab,ncolln2
      common/boundcom/rhoplasm

      nlevelb = 1

      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

c--set up contour levels

      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)

      if(nlevelb.eq.1)then
         flevel(1)= rhoplasm
      endif

c Split contours into two parts, f > 0, and f < 0.
c Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
c Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevelb
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do
  100 continue
      ilevgt = ilevlt + nlevlt
      nlevgt = nlevelb - nlevlt


c--Advance graphics frame and get ready to plot
c     call pladv(0)
      call plcol(ncolbox)
c      call plenv(xmin, xmax, ymin, ymax, 1, 0)
c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps
         call pllsty(1)
c      call plvpas(0.15,0.85,0.15,0.85,0.0)
c      call plwind(xpmin,xpmax,ypmin,ypmax)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0
c     call plbox(xopt,xtick,nxsub,yopt,ytick,nysub)

c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then
         call pllsty(4)
         call plcol(ncolbox)
         call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &      nlevlt, r, theta)
      endif

      if(nlevgt .gt. 0) then
         call pllsty(1)
         call plcol(ncolbox)
         call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &      nlevgt, r, theta)
      endif

      call plcol(ncollab)
      call pllab(titx,tity,title)
      return
      end
c
c*********************************************************************
c

      subroutine ezcon3dq(r,theta,f,flevel,nr,nth,nlevel,
     1   nrmax,nthmax,nlevmax,title,titx,tity,titz)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1    ncolelec

      real r(nrmax), theta(nthmax)
      real f(nrmax, nthmax)
      real flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity
      character*32 titz

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1    nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1    nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1    ncolion,ncolelec,ncollin,ncollab,ncolln2

      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

c--set up contour levels
      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)
      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)return
      fmax = 10.0
      df = abs(fmax - fmin) / float(nlevel)

c      write(6, 100)fmin, fmax
  100 format(1p8e12.4)


c--Advance graphics frame and get ready to plot
         call plcol(ncolbox)
         call plenv(-2.5, 2.5, -2.5, 4.0, 0, -2)
c         call plenv(-2.3, 1.8, -1.5, 3.5, 0, -2)
         call plw3d(2.,4.,3.,xmin,xmax,ymin,ymax,fmin,fmax,30.,30.)
         call plschr(0.,.50)
         call plbox3('binstu',titx,0.0,0,'binstu',
     *     tity,0.0,0,'binstu',titz,0.0,0)
         call plcol(ngreen)
         call plmesh(r,theta,f,nr,nth,3,nrmax)
         call plschr(0.,1.0)
         call plcol(ncollab)
         call plmtex('t',1.0,0.5,0.5,title)
      return
      end
c
c***************************************************************************
c
      subroutine ezplot3(title, titll, titlr, x1, y1, y2, y3, y4,
     .    nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax)
      real y1max, y2max, y3max, y4max
      real y1min, y2min, y3min, y4min
      real ymin, ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1 nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1 nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1 ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)

      ymax = amax1(y1max, y2max, y3max, y4max)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1
      ymin=0.0



      call pladv(0)
      call pllsty(1)
      call plcol(ncolbox)
      call plvpor(0.15,0.85,0.15,0.85)
      call plwind(rhomin,rhomax,ymin,ymax)
      call plbox('bcnst',0.0,0,'bnstv',0.0,0)

      call plcol(nblue)
      call plline(nr,x1,y1)
      call plptex(0.65*rhomax,.95*ymax,1.0,0.0,0.0,'minority ions')

      call plcol(nred)
      call plline(nr,x1,y2)
      call plptex(0.65*rhomax,.89*ymax,1.0,0.0,0.0,'electrons')

      call plcol(ncyan)
      call plline(nr,x1,y3)
      call plptex(0.65*rhomax,.83*ymax,1.0,0.0,0.0,'majority ions')

      call plcol(ngreen)
      call plline(nr, x1, y4)
      call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'fast ions')


      call plcol(ncollin)
      call plmtex('l',5.0,0.5,0.5,titll)
      call plmtex('b',3.2,0.5,0.5,'rho')
      call plmtex('t',2.0,0.5,0.5,title)

      call pllsty(1)
      call plcol(ncyan)
c     call plpoin(n12,x2,y3,4)
      call pllsty(1)
      call plcol(ncolbox)
      call plwind(rhomin,rhomax,ymin,ymax)
      call plbox(' ',0.0,0,'cstv',0.0,0)
c     call plcol(3)
c     call pllsty(1)
c     call plline(nr,x1,y2)
c     call plpoin(n12,x2,y4,4)

c     call pllsty(1)
c     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p9e11.3)
      return
      end

c
c***************************************************************************
c

      subroutine ezplot5p(title, titll, titlr, x1, y1, y2, y3, y4, y5,
     .    nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax),
     .     y5(nrmax)
      real y1max, y2max, y3max, y4max, y5max
      real y1min, y2min, y3min, y4min, y5min
      real ymin, ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1 nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1 nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1 ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)

      ymax = amax1(y1max, y2max, y3max, y4max, y5max)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1
      ymin=0.0



      call pladv(0)
      call pllsty(1)
      call plcol(ncolbox)
      call plvpor(0.15,0.85,0.15,0.85)
      call plwind(rhomin,rhomax,ymin,ymax)
      call plbox('bcnst',0.0,0,'bnstv',0.0,0)

      call plcol(nblue)
      call plline(nr,x1,y1)
      call plptex(0.65*rhomax,.95*ymax,1.0,0.0,0.0,'minority ions')

      call plcol(nred)
      call plline(nr,x1,y2)
      call plptex(0.65*rhomax,.89*ymax,1.0,0.0,0.0,'electrons')

      call plcol(ncyan)
      call plline(nr,x1,y3)
      call plptex(0.65*rhomax,.83*ymax,1.0,0.0,0.0,'majority ions')

      call plcol(ngreen)
      call plline(nr, x1, y4)
      call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'fast ions')

      call plcol(ngreen)
      call plline(nr, x1, y5)
      call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'fast ions')


      call plcol(ncollin)
      call plmtex('l',5.0,0.5,0.5,titll)
      call plmtex('b',3.2,0.5,0.5,'rho')
      call plmtex('t',2.0,0.5,0.5,title)

      call pllsty(1)
      call plcol(ncyan)
c     call plpoin(n12,x2,y3,4)
      call pllsty(1)
      call plcol(ncolbox)
      call plwind(rhomin,rhomax,ymin,ymax)
      call plbox(' ',0.0,0,'cstv',0.0,0)
c     call plcol(3)
c     call pllsty(1)
c     call plline(nr,x1,y2)
c     call plpoin(n12,x2,y4,4)

c     call pllsty(1)
c     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p9e11.3)
      return
      end

c
c***************************************************************************
c
      subroutine ezplot2(title, titll, titlr, x1, y1, y2, nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax), y2(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr, nrmax, n

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1 nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1 nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1 ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen



      call a1mnmx(y1, nrmax, nr, y1min, y1max)


c      write(6, 300) y2min, y2max
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return


      ymax = y1max
      ymin = y1min

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1



      call pladv(0)
      call plvpor(0.15,0.85,0.15,0.85)
      call plwind(rhomin, rhomax, ymin, ymax)
      call plcol(ncolbox)
      call pllsty(1)
      call plbox('bcnst', 0.0, 0, 'bnstv', 0.0, 0)
      call plcol(ncollin)
      call plmtex('l',5.0,0.5,0.5,titll)

      call plline(nr, x1, y1)

      call plmtex('b',3.2,0.5,0.5,'rho')
      call plmtex('t',2.0,0.5,0.5,title)


      call pllsty(1)
      call plcol(ncolbox)

      call a1mnmx(y2, nrmax, nr, y2min, y2max)
c      write(6, 300) y2min, y2max
      y2min = y2min * 1.1
      y2max = y2max * 1.1
c      write(6, 300) y2min, y2max


      call plwind(rhomin, rhomax, y2min, y2max)

      call plbox(' ', 0.0, 0, 'cmstv', 0.0, 0)

      call plcol(nblue)
      call pllsty(1)
      call plline(nr, x1, y2)





      call plcol(nblue)
      call pllsty(1)
      call plmtex('r', 5.0, 0.5, 0.5, titlr)



  300 format (1p9e11.3)
      return
      end
c
c***************************************************************************
c



      subroutine ezplot4(title, titll, titlr, x1, y1, y2, y3, y4,
     .   nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax),y2(nrmax),y3(nrmax),y4(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real y4max, y4min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1 nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1 nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1 ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      call a1mnmx(y2,nrmax,nr,y2min,y2max)
      call a1mnmx(y3,nrmax,nr,y3min,y3max)
      call a1mnmx(y4,nrmax,nr,y4min,y4max)

      ymax = amax1(y1max,y2max,y3max,y4max)
      ymin = amin1(y1min,y2min,y3min,y4min)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin * 1.1



      call pladv(0)
      call pllsty(1)
      call plcol(ncolbox)
      call plvpor(0.15,0.85,0.15,0.85)
      call plwind(rhomin,rhomax,ymin,ymax)
      call plbox('bcnst',0.0,0,'bnstv',0.0,0)

      call plcol(nblue)
      call plline(nr,x1,y1)
      call plptex(0.65*rhomax,.95*ymax,1.0,0.0,0.0,'minority ions')

      call plcol(nred)
      call plline(nr,x1,y2)
      call plptex(0.65*rhomax,.89*ymax,1.0,0.0,0.0,'electrons')

      call plcol(ncyan)
      call plline(nr,x1,y3)
      call plptex(0.65*rhomax,.83*ymax,1.0,0.0,0.0,'majority ions')

      call plcol(ngreen)
      call plline(nr,x1,y4)
      call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'fast ions')

      call plcol(ncollin)
      call plmtex('l',5.0,0.5,0.5,titll)
      call plmtex('b',3.2,0.5,0.5,'rho')
      call plmtex('t',2.0,0.5,0.5,title)

      call pllsty(1)
      call plcol(ncyan)
c     call plpoin(n12,x2,y3,4)
      call pllsty(1)
      call plcol(ncolbox)
      call plwind(rhomin,rhomax,ymin,ymax)
      call plbox(' ',0.0,0,'cstv',0.0,0)
c     call plcol(3)
c     call pllsty(1)
c     call plline(nr,x1,y2)
c     call plpoin(n12,x2,y4,4)

c     call pllsty(1)
c     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p9e11.3)
      return
      end
c
c***************************************************************************
c
      subroutine ezplot5(title,titll,titlr,x1,
     1   y1,y2,y3,y4,y5,y6,nr,nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax),y2(nrmax),y3(nrmax)
      real y4(nrmax),y5(nrmax),y6(nrmax)
      real y1max,y2max,y3max,y4max,y5max,y6max
      real y1min,y2min,y3min,y4min,y5min,y6min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
      integer ncolln4,ncolln5,ncolln6
      integer nblack,nred,nyellow, ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1 nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1 nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1 ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen
      ncolln4=nred
      ncolln5=ngrey
      ncolln6=nbrown


      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      call a1mnmx(y2,nrmax,nr,y2min,y2max)
      call a1mnmx(y3,nrmax,nr,y3min,y3max)
      call a1mnmx(y4,nrmax,nr,y4min,y4max)
      call a1mnmx(y5,nrmax,nr,y5min,y5max)
      call a1mnmx(y6,nrmax,nr,y6min,y6max)
      ymax = amax1(y1max,y2max,y3max,y4max,y5max,y6max)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1
      ymin=0.0



      call pladv(0)
      call pllsty(1)
      call plcol(ncolbox)
      call plvpor(0.15,0.85,0.15,0.85)
      call plwind(rhomin,rhomax,ymin,ymax)
      call plbox('bcnst',0.0,0,'bnstv',0.0,0)

      call plcol(ncollin)
      call pllsty(3)
      call plline(nr,x1,y1)
      call plptex(0.6*rhomax,.95*ymax,1.0,0.0,0.0,'Mode -2')

      call plcol(ncolln2)
      call pllsty(2)
      call plline(nr,x1,y2)
      call plptex(0.6*rhomax,.90*ymax,1.0,0.0,0.0,'Mode -1')

      call plcol(ncolln3)
      call pllsty(1)
      call plline(nr,x1,y3)
      call plptex(0.6*rhomax,.85*ymax,1.0,0.0,0.0,'Mode 0')

      call plcol(ncolln4)
      call pllsty(2)
      call plline(nr,x1,y4)
      call plptex(0.6*rhomax,.80*ymax,1.0,0.0,0.0,'Mode 1')

      call plcol(ncolln5)
      call pllsty(3)
      call plline(nr,x1,y5)
      call plptex(0.6*rhomax,.75*ymax,1.0,0.0,0.0,'Mode 2')



      call plcol(ncollin)
      call plmtex('l',5.0,0.5,0.5,titll)
      call plmtex('b',3.2,0.5,0.5,'rho')
      call plmtex('t',2.0,0.5,0.5,title)

      call pllsty(1)
      call plcol(ncyan)
c     call plpoin(n12,x2,y3,4)
      call pllsty(1)
      call plcol(ncolbox)
      call plwind(rhomin,rhomax,ymin,ymax)
      call plbox(' ',0.0,0,'cstv',0.0,0)
c     call plcol(3)
c     call pllsty(1)
c     call plline(nr,x1,y2)
c     call plpoin(n12,x2,y4,4)

c     call pllsty(1)
c     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p9e11.3)
      return
      end

c
c***************************************************************************
c
      subroutine ezplot1(title, titll, titlr, x1, y1, nr, nrmax,
     .   y1max_giv)

      implicit none

      real xzmax, xzmin, xnmin, xnmax, rhomin, rhomax
      real x1(nrmax),y1(nrmax)
      real y1max, y2max, y3max, y1min, y2min, y3min
      real ymin, ymax
      real y1max_giv

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink,
     1 nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan,
     1 nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd,
     1 ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return




      ymax = y1max
      ymin = y1min
      ymin = y1max - 12

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1

      if(y1max_giv .ne. 0.0)ymax = y1max_giv



      call pladv(0)
      call plvpor(0.15, 0.85, 0.15, 0.85)
      
      call plwind(rhomin, rhomax, ymin, ymax)
      call plcol(ncolbox)
      call pllsty(1)
      call plbox('bcnst',0.0,0,'lbnstv',0.0,0)
      call plcol(ncollin)
      call plmtex('l', 5.0, 0.5, 0.5, titll)

      call plcol(ncollin)
      
      call plline(nr, x1, y1)

      call plmtex('b', 3.2, 0.5, 0.5, 'u')
      call plmtex('t', 2.0, 0.5, 0.5, title)
      call plcol(ncolbox)
      call plwind(rhomin, rhomax, ymin, ymax)
      call plbox(' ', 0.0, 0, 'lcstv', 0.0, 0)

      call pllsty(1)
      call plcol(ncyan)
      call pllsty(1)


  300 format (1p9e11.3)
      return
      end
c
c***************************************************************************
c


