c
c--version 1.0 (10/16/2006)
c
      subroutine plot(ndisti2) 
      
      use size_mod             
      
      implicit none
      
      integer i_psi, i_psi1, i_psi2, i_psi3, i_psi4, i_psi5, i_psi6
      integer nphimx, nphmax, ntmax, ntmin, nchmax, jhalf, khalf, kshift
      integer ftrap, iezrant, nr, k, nt1, nt2, jmid, 
     .   nant, nphi3d, n, nstrap, nphi, nstrpol, nphi_fpm, nt_max
      integer nphi3d_1quarter, nphi3d_half, nphi3d_3quarter 
      
      integer nxmx, nymx, i, j
      integer nnodex, nnodey, nphi1, nphi2, nphi1_dum, nphi2_dum
      integer nxplot_dim, nyplot_dim
      integer nxplot, nyplot, ir, nt, imax, jmax, kmax
      
      integer pgopen, pgbeg, ier, numb, iflag, nlevmax
      integer number_points, nrhomax
      integer nnoderho
      
c      integer nmodesmax, mmodesmax
            
      parameter (nlevmax = 101)
      
c      parameter (nxmx   = 500)
c      parameter (nymx   = 500)
      parameter (nxmx = nmodesmax)
      parameter (nymx = mmodesmax)
      
      parameter (nrhomax = nxmx * 2)
      parameter (nphimx = 200)
      parameter (nxplot_dim   = 1000)
      parameter (nyplot_dim   = 1000)
      
     
      parameter (nphmax = 2048)
      parameter (ntmax = nphmax / 2, ntmin = -ntmax)
      
      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb
      character*32 titx
      character*32 tity
      character*32 titz
      character*32 tityl
      character*32 tityr
      
      real, dimension(:),   allocatable :: UPERP, UPARA
      
      real, dimension(:,:), allocatable :: bqlsum_i2_2d 
      real, dimension(:,:), allocatable :: cqlsum_i2_2d     
      real, dimension(:,:), allocatable :: eqlsum_i2_2d     
      real, dimension(:,:), allocatable :: fqlsum_i2_2d          
      
      real*8, dimension(:,:,:), allocatable :: bqlsum_i2
      real*8, dimension(:,:,:), allocatable :: cqlsum_i2
      real*8, dimension(:,:,:), allocatable :: eqlsum_i2
      real*8, dimension(:,:,:), allocatable :: fqlsum_i2
      
      integer :: i_uperp, i_upara, i_psi_eq
      real vce_mks, vc1_mks, vc2_mks, vc3_mks, vc4_mks, vc5_mks, vc6_mks 
      real vce_cgs, vc1_cgs, vc2_cgs, vc3_cgs, vc4_cgs, vc5_cgs, vc6_cgs
      real UminPara, UmaxPara, xmi2
      
      real fmod(nphimx)

      
      real xnurf, pi, omgrf, q, width, eps0, xmu0, clight, qe, absqe, 
     .     xk0, phimin, phimax, dphi, dphi3d, rt, r0, zmin, zmax, xltcm
      real xlt, phase, wd, xi0, phi0, separ0, rhoplasm, rant, rmax, rmin
      real dphase0, phase0, phi(nphmax), capr_min, capr_max
      real phi3d(nphimx), phi3d_shift(nphimx), zk3d(nphimx), dxplot
      real zk(nphmax), amplt(20), dyplot, dx, dy, emax
      real workr(nphmax), worki(nphmax), capr(nxmx), capz(nymx)
      real rho(nxmx, nymx), drho, phase_deg
      real ff(101), dummy, tmin, t1, second1, ptot, pcito2, pcrto2
      real capr_giv, phi3d_giv, fout, dcapr, x_giv, ptotal, xjtot
      
      real pcte, pcti1, pcti2, pcti3, pcti4, pcti5, pcti6, pctt, pt
      
      real  pcedotje, pcedotj1, pcedotj2, pcedotj3, pcedotj4,
     .      pcedotj5, pcedotj6, pcedotjs, pcedotjt, pedotjt
    
      
      complex work(nphmax), zi, csum, curent(nphmax),cn(ntmin : ntmax)
      complex cur3d(nphimx), jfun, cursum(nphmax)
      complex cexpkz      
      
      real rhon(nrhomax)
            
      real wdoti1_dvol(nrhomax), wdoti2_dvol(nrhomax),
     .     wdoti3_dvol(nrhomax), wdoti4_dvol(nrhomax), 
     .     wdoti5_dvol(nrhomax), wdoti6_dvol(nrhomax),
     .     wdote_dvol(nrhomax)     
     
      real redotj1_dvol(nrhomax), redotj2_dvol(nrhomax),
     .     redotj3_dvol(nrhomax), redotj4_dvol(nrhomax),
     .     redotj5_dvol(nrhomax), redotj6_dvol(nrhomax),
     .     redotje_dvol(nrhomax)
          
      real xjprlavg(nrhomax), xjprl_sum(nrhomax)
          
      real wdotesum(nrhomax),
     .     wdot1sum(nrhomax), wdot2sum(nrhomax),
     .     wdot3sum(nrhomax), wdot4sum(nrhomax),
     .     wdot5sum(nrhomax), wdot6sum(nrhomax) 
     
      real wdotesum_ql(nrhomax),
     .     wdot1sum_ql(nrhomax), wdot2sum_ql(nrhomax),
     .     wdot3sum_ql(nrhomax), wdot4sum_ql(nrhomax),
     .     wdot5sum_ql(nrhomax), wdot6sum_ql(nrhomax)      
          
      real redotjesum(nrhomax),
     .     redotj1sum(nrhomax), redotj2sum(nrhomax),
     .     redotj3sum(nrhomax), redotj4sum(nrhomax),
     .     redotj5sum(nrhomax), redotj6sum(nrhomax)
     
      real redotje_int(nrhomax),
     .     redotji1_int(nrhomax), redotji2_int(nrhomax),
     .     redotji3_int(nrhomax), redotji4_int(nrhomax), 
     .     redotji5_int(nrhomax), redotji6_int(nrhomax)
     
      real wdote_int(nrhomax),
     .     wdoti1_int(nrhomax), wdoti2_int(nrhomax), 
     .     wdoti3_int(nrhomax), wdoti4_int(nrhomax),
     .     wdoti5_int(nrhomax), wdoti6_int(nrhomax)
     
      real xjprl_int(nrhomax)
     
      real  dvol(nrhomax), darea(nrhomax)        
      

      real, dimension (:, :), allocatable :: freal
      real, dimension (:, :), allocatable :: fimag
      real, dimension (:, :), allocatable :: power
      real, dimension (:, :), allocatable :: capr_plot
      
      real, dimension (:, :), allocatable :: freal_phi
      real, dimension (:, :), allocatable :: fimag_phi
      real, dimension (:, :), allocatable :: power_phi    
      real, dimension (:, :), allocatable :: capr_plot_phi 
      
      real fmidre(nxmx), fmidim(nxmx)
      
      real cnmod2(nphimx), xnphi(nphimx), pabs(nphimx), jdriven(nphimx)
      real pabs_weight, pabs_sum, prfin, pscale, j_driven_weight
      real jdriven_weight, jdriven_sum
      
      real xplot(nxplot_dim), yplot(nyplot_dim)
      real xplotm(nxplot_dim), yplotm(nyplot_dim)
      
      real, dimension (:, :), allocatable :: frealp
      real, dimension (:, :), allocatable :: fimagp
      real, dimension (:, :), allocatable :: powerp
      real, dimension (:, :), allocatable :: capr_plotp                  
      
      real, dimension (:, :, :), allocatable :: ealpha_sum_mod  
      real, dimension (:, :, :), allocatable :: ebeta_sum_mod           
      complex, dimension(:,:,:), allocatable :: ealpha_sum
      complex, dimension(:,:,:), allocatable :: ebeta_sum
      complex, dimension(:,:,:), allocatable :: eb_sum            
      real, dimension (:, :, :), allocatable :: redotj 
      
      complex, dimension(:, :, :), allocatable :: ntilda       
            
      common/boundcom/rhoplasm      
      
      integer nmodesx, nmodesy, nwdot, lmax, ibessel,
     .    nuper, nupar, 
     .    inu, iprint, iexact,
     .    iroot, iequat, igeom, 
     .    iqx, iqprof, iez, nprow, npcol, 
     .    izfunc,
     .    nstep, nabs, 
     .    isigma, itemp, 
     .    nfreqm,  nkzm, 
     .    ibackground, iabsorb, nzfun, 
     .    nboundary, nnode_local,
     .    nnode_overlap, iprofile, isolve, 
     .    ndiste, ndisti1, ndisti2, ndisti3,
     .    ndisti4, ndisti5, ndisti6, nkperp, nzeta_wdot, n_bin, 
     .    iql, i_antenna, n_prof_flux, i_write, idens
     
      CHARACTER*128 :: eqdsk
      CHARACTER*128 :: netCDF_file1
      CHARACTER*128 :: netCDF_file2      
      
      real ti0, xnuead, xnu1ad, xnu2ad, te0, yant,
     .    ti02, ti03, ti2lim, ti3lim,  
     .    ti04, ti05, ti06, ti4lim, ti5lim, ti6lim,
     .    delta0, xwall, xnwall,
     .     epszet,
     .    dthetant0, dpsiant0, psilim, psiant, psimol, psipne, psipte,
     .    psipti1, psipti2, psipti3,
     .    psipti4, psipti5, psipti6,
     .    amu1, amu2, z1, z2, eta, 
     .    b0, ytop, ybottom, freqcy, aplasm,
     .    xnlim, xn2lim, xn3lim, xnslolim, signbz,
     .    xn4lim, xn5lim, xn6lim,
     .    xn0, xn2, xn3, xnslo, flat, b1rat, b2rat, curdnx, curdny,
     .    xn4, xn5, xn6,
     .    curdnz,
     .    xnuabs, xbnch, xleft, xright,
     .    telim, tilim,
     .    dfreq, dkz,
     .    xnudip, adip, efold,
     .    amu3, z3, eta3, xnu3ad,
     .    amu4, z4, eta4, xnu4ad,
     .    amu5, z5, eta5, xnu5ad,
     .    amu6, z6, eta6, xnu6ad,
     .    xdelta, wdelta, xdelt2, wdelt2, zeffcd,
     .    rzoom1, rzoom2, yzoom1, yzoom2,  q0,
     .    alim, grad, qavg0, ymax,
     .    alphan,  alphan2, alphan3, alphan_slo,
     .    alphan4, alphan5, alphan6,
     .    alphate,  alphati, alphati2, alphati3,
     .    alphati4, alphati5, alphati6,
     .    ekappa, rwleft, rwright, xnuomg,
     .    eta_slo, amu_slo, z_slo, eslowev, 
     .    betan, betan2, betan3, betan_slo, betate, betati, betati2,
     .    betan4, betan5, betan6, betati4, betati5, betati6,
     .    betati3, taue, theta_ant,
     .    antlen,
     .    antlc, upshift, xkperp_cutoff, damping
     
      allocate( ealpha_sum_mod(nxmx, nymx, nphimx)  )
      allocate( ebeta_sum_mod(nxmx, nymx, nphimx)  )
      allocate( ealpha_sum(nxmx, nymx, nphimx) )
      allocate( ebeta_sum(nxmx, nymx, nphimx) )
      allocate( eb_sum(nxmx, nymx, nphimx) ) 
      allocate( redotj(nxmx, nymx, nphimx)  )
      
      allocate( ntilda(nxmx, nymx, nphimx) )         
      
      
      allocate( frealp(nxplot_dim, nyplot_dim) )
      allocate( fimagp(nxplot_dim, nyplot_dim) )
      allocate( powerp(nxplot_dim, nyplot_dim) )
      allocate( capr_plotp(nxplot_dim, nyplot_dim) ) 
      
      allocate( freal(nxmx, nymx) )
      allocate( fimag(nxmx, nymx) )
      allocate( capr_plot(nxmx, nymx) )
      allocate( power(nxmx, nymx) ) 
      
      allocate( freal_phi(nxmx, nphimx) )
      allocate( fimag_phi(nxmx, nphimx) )
      allocate( capr_plot_phi(nxmx, nphimx) )
      allocate( power_phi(nxmx, nphimx) ) 
   
      
            
c      open(unit=63, file='aorsa2d.in', status='old', form='formatted')
c      rewind(63)
      
     
      open(unit=242,file='Bql_sum_2D.vtk',status='unknown',
     .                                                 form='formatted') 
      open(unit=342,file='Cql_sum_2D.vtk',status='unknown',
     .                                                 form='formatted')      
      open(unit=362,file='murakami3d',status='unknown',form='formatted') 
      
*     ---------------------------------------
*     Open and read file "out36" for plotting
*     ---------------------------------------
      open(unit=36, file='out36', status='unknown', form= 'formatted')
      open(unit=37, file='ryan',  status='unknown', form= 'formatted')
      
c      read (36, 309) nnodex, nnodey, nphi3d
c      read (36, 310) (capr(i),  i = 1, nnodex)
c      read (36, 310) (capz(j),  j = 1, nnodey)
c      read (36, 310) (phi3d(k), k = 1, nphi3d)
      
c      read (36, 310) (((ealpha_sum(i, j, k), i = 1, nnodex), 
c     .                                       j = 1, nnodey), 
c     .                                       k = 1, nphi3d)

      read (36, 309) nphi1, nphi2, nt_max, nnoderho
      read (36, 310) prfin      
      
      read (36, 310) (xnphi(nt),  nt = 1, nt_max)
      read (36, 310) (cnmod2(nt), nt = 1, nt_max)
      read (36, 310) (pabs(nt),   nt = 1, nt_max) 
      read (36, 310) (jdriven(nt),   nt = 1, nt_max)
      
      
      read (36, 310) (rhon(n),    n = 1, nnoderho) 
      read (36, 310) (redotj2sum(n), n = 1, nnoderho) 
      read (36, 310) (redotjesum(n), n = 1, nnoderho) 
      read (36, 310) (redotj1sum(n), n = 1, nnoderho) 
      read (36, 310) (redotj3sum(n), n = 1, nnoderho)  
      read (36, 310) (redotj4sum(n), n = 1, nnoderho)  
      read (36, 310) (redotj5sum(n), n = 1, nnoderho) 
      
      write (37, 310) (rhon(n),    n = 1, nnoderho)
             
      read (36, 310) (redotji2_int(n), n = 1, nnoderho) 
      read (36, 310) (redotje_int(n),  n = 1, nnoderho) 
      read (36, 310) (redotji1_int(n), n = 1, nnoderho) 
      read (36, 310) (redotji3_int(n), n = 1, nnoderho)  
      read (36, 310) (redotji4_int(n), n = 1, nnoderho)  
      read (36, 310) (redotji5_int(n), n = 1, nnoderho) 
      
      read (36, 310) (xjprl_sum(n), n = 1, nnoderho)  
      read (36, 310) (xjprl_int(n), n = 1, nnoderho)
      
      read (36, 310) (wdot2sum(n), n = 1, nnoderho) 
      read (36, 310) (wdotesum(n), n = 1, nnoderho) 
      read (36, 310) (wdot1sum(n), n = 1, nnoderho) 
      read (36, 310) (wdot3sum(n), n = 1, nnoderho)  
      read (36, 310) (wdot4sum(n), n = 1, nnoderho)  
      read (36, 310) (wdot5sum(n), n = 1, nnoderho)
      

      read (36, 310) (wdot2sum_ql(n), n = 1, nnoderho) 
      read (36, 310) (wdotesum_ql(n), n = 1, nnoderho) 
      read (36, 310) (wdot1sum_ql(n), n = 1, nnoderho) 
      read (36, 310) (wdot3sum_ql(n), n = 1, nnoderho)  
      read (36, 310) (wdot4sum_ql(n), n = 1, nnoderho)  
      read (36, 310) (wdot5sum_ql(n), n = 1, nnoderho)
      
      write (37, 310) (wdot2sum_ql(n), n = 1, nnoderho) 
      write (37, 310) (wdotesum_ql(n), n = 1, nnoderho) 
      write (37, 310) (wdot1sum_ql(n), n = 1, nnoderho) 
      write (37, 310) (wdot3sum_ql(n), n = 1, nnoderho)  
      write (37, 310) (wdot4sum_ql(n), n = 1, nnoderho)  
      write (37, 310) (wdot5sum_ql(n), n = 1, nnoderho)        
               
       
      read (36, 310) (wdoti2_int(n), n = 1, nnoderho) 
      read (36, 310) (wdote_int(n),  n = 1, nnoderho) 
      read (36, 310) (wdoti1_int(n), n = 1, nnoderho) 
      read (36, 310) (wdoti3_int(n), n = 1, nnoderho)  
      read (36, 310) (wdoti4_int(n), n = 1, nnoderho)  
      read (36, 310) (wdoti5_int(n), n = 1, nnoderho)       

      read (36, 310) (dvol(n), n = 1, nnoderho) 
      
      read (36, 309) nnodex, nnodey, nphi3d
      read (36, 310) (capr(i),  i = 1, nnodex)
      read (36, 310) (capz(j),  j = 1, nnodey)
      read (36, 310) (phi3d(k), k = 1, nphi3d)
      read (36, 310) (phi3d_shift(k), k = 1, nphi3d)   
      read (36, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      read (36, 310) rhoplasm, rt
            
      read (36, 310) (((ealpha_sum(i, j, k), i = 1, nnodex), 
     .                                       j = 1, nnodey), 
     .                                       k = 1, nphi3d)
     
      read (36, 310) (((ealpha_sum_mod(i, j, k), i = 1, nnodex), 
     .                                           j = 1, nnodey), 
     .                                           k = 1, nphi3d)     
      
      read (36, 310) (((ebeta_sum(i, j, k), i = 1, nnodex), 
     .                                      j = 1, nnodey), 
     .                                      k = 1, nphi3d) 
     
      read (36, 310) (((ebeta_sum_mod(i, j, k), i = 1, nnodex), 
     .                                          j = 1, nnodey), 
     .                                          k = 1, nphi3d) 
     
      read (36, 310) (((redotj(i, j, k), i = 1, nnodex), 
     .                                   j = 1, nnodey), 
     .                                   k = 1, nphi3d) 
     
      read (36, 309) nxplot, nyplot
      read (36, 310) capr_max, capr_min
      read (36, 310) (xplotm(i),  i = 1, nxplot)      
      read (36, 310) (yplotm(j),  j = 1, nyplot)       
      read (36, 310) ((frealp(i, j), i = 1, nxplot), j = 1, nyplot)
      read (36, 310) ((fimagp(i, j), i = 1, nxplot), j = 1, nyplot) 
      read (36, 310) ((capr_plotp(i, j), i = 1, nxplot), j = 1, nyplot)
      read (36, 310) ((powerp(i, j), i = 1, nxplot), j = 1, nyplot) 
      
      read (36, 310) (((ntilda(i, j, k), i = 1, nnodex), 
     .                                   j = 1, nnodey), 
     .                                   k = 1, nphi3d)      
      
c      write(6, *)"xplotm(i) = "
c      write (6, 310) (xplotm(i),  i = 1, nxplot)       
      

*     ----------------------------------
*     Open and read the SWIM output file:
*     ----------------------------------
         open(unit=99, file='out_swim', status='unknown',
     .                                              form='formatted')                  
         read(99, 309) nnoderho
         read(99, 310) (rhon(n), n = 1, nnoderho)
	 
	 read(99, 310) (dvol(n), n = 1, nnoderho - 1)
	 
	 read(99, 310) (redotje_int(n), n = 1, nnoderho)
         read(99, 310) (redotji1_int(n), n = 1, nnoderho) 
	 read(99, 310) (redotji2_int(n), n = 1, nnoderho) 
         read(99, 310) (redotji3_int(n), n = 1, nnoderho)  
         read(99, 310) (redotji4_int(n), n = 1, nnoderho)  
         read(99, 310) (redotji5_int(n), n = 1, nnoderho)
         read(99, 310) (redotji6_int(n), n = 1, nnoderho)
	 
	 
	 read(99, 310) (wdote_int(n), n = 1, nnoderho)
         read(99, 310) (wdoti1_int(n), n = 1, nnoderho)
         read(99, 310) (wdoti2_int(n), n = 1, nnoderho)
         read(99, 310) (wdoti3_int(n), n = 1, nnoderho)
         read(99, 310) (wdoti4_int(n), n = 1, nnoderho)
	 read(99, 310) (wdoti5_int(n), n = 1, nnoderho)
	 read(99, 310) (wdoti6_int(n), n = 1, nnoderho)
	 
	 read(99, 310) (redotje_dvol(n), n = 1, nnoderho)
         read(99, 310) (redotj1_dvol(n), n = 1, nnoderho)
         read(99, 310) (redotj2_dvol(n), n = 1, nnoderho)
         read(99, 310) (redotj3_dvol(n), n = 1, nnoderho)
         read(99, 310) (redotj4_dvol(n), n = 1, nnoderho)
         read(99, 310) (redotj5_dvol(n), n = 1, nnoderho)
         read(99, 310) (redotj6_dvol(n), n = 1, nnoderho)
	 
	 
         read(99, 310) (wdote_dvol(n), n = 1, nnoderho)
         read(99, 310) (wdoti1_dvol(n), n = 1, nnoderho)
         read(99, 310) (wdoti2_dvol(n), n = 1, nnoderho)
         read(99, 310) (wdoti3_dvol(n), n = 1, nnoderho)
         read(99, 310) (wdoti4_dvol(n), n = 1, nnoderho)
	 read(99, 310) (wdoti5_dvol(n), n = 1, nnoderho)
	 read(99, 310) (wdoti6_dvol(n), n = 1, nnoderho)   
      close (99)       
      
                              
                      
      close (36)
      close (37)      
      


      i_psi1 = 2
      i_psi2 = 10      
      i_psi3 = 20
      i_psi4 = 30
      i_psi5 = 40
      i_psi6 = 50
	 
	       
*     ----------------------------------------------------------------
*     read quasilinear diffusion coefficients from file out_cql3d.coef2
*     ----------------------------------------------------------------	


      if(ndisti2  .eq. 1) then

         open(unit=42, file='out_cql3d.coef2',  
     .                               status='unknown', form='formatted')
      	 
         read (42, 309) nuper
         read (42, 309) nupar
         read (42, 309) nnoderho
      
         allocate( UPERP(nuper) )
         allocate( UPARA(nupar) ) 
      
         allocate( bqlsum_i2(nuper, nupar, nnoderho) )
         allocate( cqlsum_i2(nuper, nupar, nnoderho) )
         allocate( eqlsum_i2(nuper, nupar, nnoderho) )
         allocate( fqlsum_i2(nuper, nupar, nnoderho) ) 
      
         allocate( bqlsum_i2_2d(nupar, nuper) )
         allocate( cqlsum_i2_2d(nupar, nuper) )
	 allocate( eqlsum_i2_2d(nupar, nuper) )
         allocate( fqlsum_i2_2d(nupar, nuper) )      
              	       

         read (42, 3310) vc2_cgs
         read (42, 3310) UminPara, UmaxPara

         read (42, 3310) (rhon(n), n = 1, nnoderho)
         read (42, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
         read (42, 3310) (upara(i_upara), i_upara = 1, nupar)
      
         read (42, 3310) (((bqlsum_i2(i_uperp, i_upara, n),
     .        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) (((cqlsum_i2(i_uperp, i_upara, n),
     .        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) (((eqlsum_i2(i_uperp, i_upara, n),
     .        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) (((fqlsum_i2(i_uperp, i_upara, n),
     .        i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
     
         read (42, 3310) xmi2	         

         close (42)
	 
      end if 
      
      write(362, 309) nnoderho
      do n = 1, nnoderho
         write(362,1312)n, rhon(n), wdotesum(n), wdot1sum(n),
     .              wdot2sum(n), wdot3sum(n), xjprl_sum(n)
      end do      
      
      
c      call exit               

*     ----------------
*     Plot with pgplot
*     ---------------- 

      t1 = second1(dummy)      
      
*     --------------------
*     Open graphics device
*     --------------------        

c      IER = PGBEG(0, 'aorsa2dSum.ps/cps', 2, 2)
      IER = PGBEG(0, 'aorsa2dSum.ps/vcps', 1, 1)
      
      
      IF (IER.NE.1) STOP
      
c      if (pgopen('aorsa2d.cps/cps') .lt. 1) stop

      call PGSCH (1.5)

      title = 'Toroidal spectrum (power)'
      titx  = 'toroidal mode number'
      tityl = 'weight'
      tityr = 'pabs'
              
      call ezplot2_sum(title, tityl, tityr, titx, xnphi, cnmod2, 
     .    pabs, nt_max, nphimx)
          
     
      tityr = ' '
     
      call ezplot1_sum(title, tityl, tityl, tityr, xnphi, cnmod2, 
     .    nt, nphimx)
     
     
     
c      title = 'Mod(E_alpha)'
c      titx  = 'phi'
c      tityl = 'Mod(E_alpha)'
c      tityr = ' '

c      i = nnodex / 2
c      j = nnodey / 2 
c      khalf = nphi3d / 2   
      
c      do k = 1, nphi3d	 
c	 if (k .lt. khalf) kshift = k + khalf
c	 if (k .ge. khalf) kshift = k - khalf + 1	    
c         fmod(k) = real(ealpha_sum_mod(i, j, kshift))	    
c      end do      
     
c      call ezplot1_sum(title, titx, tityl, tityr, phi3d_shift, fmod, 
c     .    nphi3d, nphimx)
     
          
     
     
      title = 'Toroidal spectrum (current)'
      titx  = 'toroidal mode number'
      tityl = 'weight'
      tityr = 'current'
      
      
              
      call ezplot2_sum(title, tityl, tityr, titx, xnphi, cnmod2, 
     .    jdriven, nt_max, nphimx)
     
     
      if(ndisti2  .eq. 1) then     
     
!        --------------------------------------
!        2D plots of bqlsum(u_perp, u_parallel)
!        --------------------------------------

         titx = 'u_parallel'
         tity = 'u_perp'

         write(6, *)"i_psi1 = ", i_psi1
         write(6, *)"i_psi2 = ", i_psi2
         write(6, *)"i_psi3 = ", i_psi3
         write(6, *)"i_psi4 = ", i_psi4
         write(6, *)"i_psi5 = ", i_psi5
         write(6, *)"i_psi6 = ", i_psi6
      
      
         write(16, *)"i_psi1 = ", i_psi1
         write(16, *)"i_psi2 = ", i_psi2
         write(16, *)"i_psi3 = ", i_psi3
         write(16, *)"i_psi4 = ", i_psi4
         write(16, *)"i_psi5 = ", i_psi5
         write(16, *)"i_psi6 = ", i_psi6      
	    
         numb = 15      
      
      
      
         i_psi = i_psi1      
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlsum_i2_2d(i_upara, i_uperp) =
     .                     bqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do
	 
         title = 'bqlsum_psi1(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, bqlsum_i2_2d, ff, nupar, nuper, 
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)
     
     
!        ---------------------------------------------
!        Write Bql_sum_2D.vtk file "structured points"
!        ---------------------------------------------
         number_points = nupar * nuper
         dx = upara(2) - upara(1)
         dy = uperp(2) - uperp(1)
	 
         write(242, 2840) 
 2840    format('# vtk DataFile Version 2.0')
 
         write(242, 4850)
 4850    format('Quasilinear diffusion coefficients') 
            
         write(242, 2846)
 2846    format('ASCII') 
                    
         write(242, 2847) 
 2847    format('DATASET STRUCTURED_POINTS')
          
         write(242, 2841) nupar,  nuper, 1 
 2841    format('DIMENSIONS', 3i8) 
          
         write(242, 2842) upara(1), uperp(1), 0 
 2842    format('ORIGIN', 2f9.3, 1i8)
 	 
         write(242, 2843) dx, dy, 1 
 2843    format('SPACING', 2f9.3, 1i8) 
            
         write(242, 2844) number_points 
 2844    format('POINT_DATA', i10) 
   
	                 
         write(242, 4851)
 4851    format('SCALARS bqlsum_psi1 float 1') 
              
         write(242, 2849)
 2849    format('LOOKUP_TABLE default')
                      
         write(242, 3411) ((bqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper)
                    
 3410    format(1p4e10.2)  
 3411    format(6f16.4)   
 
 
         i_psi = i_psi2      
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlsum_i2_2d(i_upara, i_uperp) =
     .                     bqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do
	 
         title = 'bqlsum_psi2(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, bqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(242, 4852)
 4852    format('SCALARS bqlsum_psi2 float 1')              
         write(242, 2849)                   
         write(242, 3411) ((bqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper)  
     
     
         i_psi = i_psi3      
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlsum_i2_2d(i_upara, i_uperp) =
     .                     bqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do
	 
         title = 'bqlsum_psi3(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, bqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(242, 4853)
 4853    format('SCALARS bqlsum_psi3 float 1')              
         write(242, 2849)                   
         write(242, 3411) ((bqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper)  
     
     
         i_psi = i_psi4     
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlsum_i2_2d(i_upara, i_uperp) =
     .                     bqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do
	 
         title = 'bqlsum_psi4(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, bqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(242, 4854)
 4854    format('SCALARS bqlsum_psi4 float 1')              
         write(242, 2849)                   
         write(242, 3411) ((bqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper) 
     
     
         i_psi = i_psi5     
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlsum_i2_2d(i_upara, i_uperp) =
     .                     bqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do
	 
         title = 'bqlsum_psi5(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, bqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(242, 4855)
 4855    format('SCALARS bqlsum_psi5 float 1')              
         write(242, 2849)                   
         write(242, 3411) ((bqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper) 
     
     
     
         i_psi = i_psi6    
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlsum_i2_2d(i_upara, i_uperp) =
     .                     bqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do
	 
         title = 'bqlsum_psi6(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, bqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(242, 4856)
 4856    format('SCALARS bqlsum_psi6 float 1')              
         write(242, 2849)                   
         write(242, 3411) ((bqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper) 
     
      end if                             
     
     
      if(ndisti2  .eq. 1) then     
     
!        --------------------------------------
!        2D plots of cqlsum(u_perp, u_parallel)
!        --------------------------------------  
	    
         numb = 15            
      
      
         i_psi = i_psi1      
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlsum_i2_2d(i_upara, i_uperp) =
     .                     cqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+26
c               write(6,*)i_upara, i_uperp, cqlsum_i2_2d(i_upara,i_uperp)
            end do
         end do
	 
         title = 'cqlsum_psi1(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, cqlsum_i2_2d, ff, nupar, nuper, 
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)
     
     
!        ---------------------------------------------
!        Write Cql_sum_2D.vtk file "structured points"
!        ---------------------------------------------
         number_points = nupar * nuper
         dx = upara(2) - upara(1)
         dy = uperp(2) - uperp(1)
	 
         write(342, 2840)  
         write(342, 4850)            
         write(342, 2846)
                    
         write(342, 2847) 
          
         write(342, 2841) nupar,  nuper, 1 
          
         write(342, 2842) upara(1), uperp(1), 0 
 	 
         write(342, 2843) dx, dy, 1 
            
         write(342, 2844) number_points 
   
	                 
         write(342, 5851)
 5851    format('SCALARS cqlsum_psi1 float 1') 
              
         write(342, 2849)
                      
         write(342, 3411) ((cqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper)
                    

 
 
         i_psi = i_psi2      
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlsum_i2_2d(i_upara, i_uperp) =
     .                     cqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do
	 
         title = 'cqlsum_psi2(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, cqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(342, 5852)
 5852    format('SCALARS cqlsum_psi2 float 1')              
         write(342, 2849)                   
         write(342, 3411) ((cqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper)  
     
     
         i_psi = i_psi3      
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlsum_i2_2d(i_upara, i_uperp) =
     .                     cqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do
	 
         title = 'cqlsum_psi3(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, cqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(342, 5853)
 5853    format('SCALARS cqlsum_psi3 float 1')              
         write(342, 2849)                   
         write(342, 3411) ((cqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper)  
     
     
         i_psi = i_psi4     
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlsum_i2_2d(i_upara, i_uperp) =
     .                     cqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do
	 
         title = 'cqlsum_psi4(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, cqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(342, 5854)
 5854    format('SCALARS cqlsum_psi4 float 1')              
         write(342, 2849)                   
         write(342, 3411) ((cqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper) 
     
     
         i_psi = i_psi5     
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlsum_i2_2d(i_upara, i_uperp) =
     .                     cqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do
	 
         title = 'cqlsum_psi5(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, cqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(342, 5855)
 5855    format('SCALARS cqlsum_psi5 float 1')              
         write(342, 2849)                   
         write(342, 3411) ((cqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper) 
     
     
     
         i_psi = i_psi6    
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               cqlsum_i2_2d(i_upara, i_uperp) =
     .                     cqlsum_i2(i_uperp, i_upara, i_psi) / 1.0e+26
            end do
         end do
	 
         title = 'cqlsum_psi6(u_perp, u_parallel)'        
         call ezconc_sum(upara, uperp, cqlsum_i2_2d, ff, nupar, nuper,
     .      numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)      
     
         write(342, 5856)
 5856    format('SCALARS cqlsum_psi6 float 1')              
         write(342, 2849)                   
         write(342, 3411) ((cqlsum_i2_2d(i,j), i = 1, nupar),  
     .                                         j = 1, nuper) 
     
      end if                             
          
     
     
      title= 'Flux surface average redotj'
      titll= 'redotj (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot6(title, titll, titlr, titlb, rhon, redotj2sum,
     .    redotjesum, redotj1sum, redotj3sum, redotj4sum, redotj5sum,
     .    nnoderho, nrhomax) 
     
      title= 'Flux surface average redotj1'
     
      call ezplot1_sum(title, titlb, titll, titlr, rhon, redotj1sum,
     .      nnoderho, nrhomax) 
     
      title= 'Flux surface average redotj2'
     
      call ezplot1_sum(title, titlb, titll, titlr, rhon, redotj2sum,
     .      nnoderho, nrhomax)      
                      

      title= 'Integrated redotj'
      titll= 'P (watts)'
      titlr= '       '
      titlb= 'rho'
     
      call ezplot6(title, titll, titlr, titlb, rhon, redotji2_int,
     .    redotje_int, redotji1_int, redotji3_int, redotji4_int, 
     .    redotji5_int, nnoderho, nrhomax)
          
      write(16, *)
      write(16, *) 'Flux surface driven current'
      write(16, *)
      write(16, *) '        n      rho       J (A/m2)     I (A)'
      write(16, *)
      
      do n = 1, nnoderho       
	 write (16, 1312) n, rhon(n), xjprl_sum(n), xjprl_int(n)
     	 write (6, 1312) n, rhon(n), xjprl_sum(n), xjprl_int(n)
      end do 
         
                                                    
      title= 'Flux surface driven current'
      titll= 'xjprl (Amps/m2)'
      titlr= 'I (Amps)'
      titlb= 'rho'

      call ezplot2_sum(title, titll, titlr, titlb, rhon, xjprl_sum,
     .    xjprl_int, nnoderho, nrhomax)
                         
      if(xjprl_int(nnoderho) .lt. 0.0) then 
      
         titll= '-xjprl (Amps/m2)'
         titlr= '-I (Amps)' 
	    
         xjprl_sum = - xjprl_sum
         xjprl_int = - xjprl_int
      
         call ezplot2_sum(title, titll, titlr, titlb, rhon, xjprl_sum,
     .      xjprl_int, nnoderho, nrhomax) 
     
         titll= 'xjprl (MA/m2/MW)'
	 xjprl_sum = xjprl_sum  / prfin	 
     
         call ezplot1_sum(title, titlb, titll, titlr, rhon, xjprl_sum,
     .      nnoderho, nrhomax)
     
         titll= 'I (kA)'
	 xjprl_int  = xjprl_int / 1.0e+03     
         call ezplot1_sum(title, titlb, titll, titlr, rhon, xjprl_int,
     .      nnoderho, nrhomax)           
                                       
      end if
            
      if(xjprl_int(nnoderho) .gt. 0.0) then 
      
         titll= 'xjprl (Amps/m2)'
         titlr= 'I (Amps)' 
	    
         xjprl_sum = xjprl_sum
         xjprl_int = xjprl_int
           
         titll= 'xjprl (MA/m2/MW)'
	 xjprl_sum = xjprl_sum  / prfin	 
     
         call ezplot1_sum(title, titlb, titll, titlr, rhon, xjprl_sum,
     .      nnoderho, nrhomax)
     
         titll= 'I (kA)'
	 xjprl_int  = xjprl_int / 1.0e+03     
         call ezplot1_sum(title, titlb, titll, titlr, rhon, xjprl_int,
     .      nnoderho, nrhomax)           
                                       
      end if      
    
      title= 'Flux surface average wdot'
      titll= 'wdot (watts/m3)'
      titlr= '       '
      titlb= 'rho'
     
      call ezplot6(title, titll, titlr, titlb, rhon, wdot2sum,
     .    wdotesum, wdot1sum, wdot3sum, wdot4sum, wdot5sum,
     .    nnoderho, nrhomax)
     
     
      title= 'Integrated wdot'
      titll= 'P (watts)'
      titlr= '       '
      titlb= 'rho'
     
      call ezplot6(title, titll, titlr, titlb, rhon, wdoti2_int,
     .    wdote_int, wdoti1_int, wdoti3_int, wdoti4_int, wdoti5_int,
     .    nnoderho, nrhomax)
     
     
      title= 'Quasilinear wdot'
      titll= 'wdot (watts/m3)'
      titlr= '       '
      titlb= 'rho'
     
      call ezplot6(title, titll, titlr, titlb, rhon, wdot2sum_ql,
     .    wdotesum_ql, wdot1sum_ql, wdot3sum_ql, wdot4sum_ql, 
     .    wdot5sum_ql, nnoderho, nrhomax) 
     
      title= 'Quasilinear wdot1'
     
      call ezplot1_0_sum(title, titlb, titll, titlr, rhon, wdot1sum_ql,
     .      nnoderho, nrhomax)
     
      title= 'Quasilinear wdot2'
     
      call ezplot1_0_sum(title, titlb, titll, titlr, rhon, wdot2sum_ql,
     .      nnoderho, nrhomax)         
     
     
      title= 'Volume element, dvol'
      titll= 'dvol (m3)'
      titlr='       '
      titlb= 'rho'
	    
      call ezplot1_sum(title, titlb, titll, titlr, rhon, dvol, 
     .    nnoderho, nrhomax)
     
     
      numb = 19

      titx = 'R (m)'
      tity = 'Z (m)'

*     -------------------------------------------------------------------
*     plot electric field in poloidal cross-section (midplane of antenna)
*     ------------------------------------------------------------------- 

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ealpha_sum(i, j, 1))
            fimag(i,j) = imag(ealpha_sum(i, j, 1))
	    power(i,j) = real(redotj(i, j, 1)) 
         end do
      end do 

      title = 'Real E alpha'
      
      call ezconc_sum(capr, capz, freal, ff, nnodex, nnodey, numb,
     .   nxmx, nymx, nlevmax, title, titx, tity, iflag)
     
      write(6, *) "iflag = ", iflag
     
      if (iflag .eq. 0) call boundary (capr, capz, rho, ff, nnodex,
     .   nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity)
     
     
*     -------------------------------------------------------------------
*     plot fluctuation level in poloidal cross-section (midplane of antenna)
*     ------------------------------------------------------------------- 

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ntilda(i, j, 1))
            fimag(i,j) = imag(ntilda(i, j, 1))
         end do
      end do 

      title = 'Real ntilda'
      
      call ezconc_sum(capr, capz, freal, ff, nnodex, nnodey, numb,
     .   nxmx, nymx, nlevmax, title, titx, tity, iflag)
     
      write(6, *) "iflag = ", iflag
     
      if (iflag .eq. 0) call boundary (capr, capz, rho, ff, nnodex,
     .   nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity)
          
                      
*     ------------------------
*     plot E alpha in midplane
*     ------------------------

      j = 64 
      k = 1    

      do i = 1, nnodex
         fmidre(i) = real(ealpha_sum(i, j, k))
         fmidim(i) = imag(ealpha_sum(i, j, k))
      end do
      
      title = 'Ealpha'
      titll = 'Re Ealpha (V/m)'
      titlr = 'Im Ealpha (V/m)'
      titlb = 'R (m)'
	    
      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, 
     .   nnodex, nxmx)      


               
*     -----------------------------------------
*     plot mod E in midplane (plane of antenna)
*     -----------------------------------------

      j = 64 
      k = 1    
 
      do i = 1, nnodex
         fmidre(i) = ealpha_sum_mod(i, j, k)
         fmidim(i) = ealpha_sum_mod(i, j, k)
      end do
      
      title = 'Mod ealpha'
      titll = 'Mod ealpha (V/m)'
      titlr = ' '
      titlb = 'R (m)'
	    
      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, 
     .   nnodex, nxmx)
     
        
      do i = 1, nnodex
         fmidre(i) = ebeta_sum_mod(i, j, k)
         fmidim(i) = ebeta_sum_mod(i, j, k)
      end do
      
      title = 'Mod ebeta'
      titll = 'Mod ebeta (V/m)'
      titlr = ' '
      titlb = 'R (m)'
	    
      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, 
     .   nnodex, nxmx)     
     
     
      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ealpha_sum(i, j, 1))
            fimag(i,j) = imag(ealpha_sum(i, j, 1))
	    power(i,j) = real(redotj(i, j, 1)) 
         end do
      end do      
          
     
      title = 'Imag E alpha'
     
       call ezconc_sum(capr, capz, fimag, ff, nnodex, nnodey, numb,
     .   nxmx, nymx, nlevmax, title, titx, tity, iflag)
     
      if (iflag .eq. 0) call boundary (capr, capz, rho, ff, nnodex,
     .   nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity)
     
     
*     ------------------------------------
*     plot power in poloidal cross section
*     ------------------------------------
            
      title = 'Real E dot J'
      
      call ezconc_sum(capr, capz, power, ff, nnodex, nnodey, numb,
     .   nxmx, nymx, nlevmax, title, titx, tity, iflag)
     
      if (iflag .eq. 0) call boundary (capr, capz, rho, ff, nnodex,
     .   nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity)
     
     
      nphi3d_1quarter = nphi3d / 4 * 1
      nphi3d_half     = nphi3d / 4 * 2
      nphi3d_3quarter = nphi3d / 4 * 3   
          
         
*     -----------------------------------------------------------
*     plot electric field in poloidal cross-section (one quarter)
*     ----------------------------------------------------------- 

      do i = 1, nnodex          
         do j = 1, nnodey
            freal(i,j) = real(ealpha_sum(i, j, nphi3d_1quarter))
            fimag(i,j) = imag(ealpha_sum(i, j, nphi3d_1quarter))
	    power(i,j) = real(redotj(i, j, nphi3d_1quarter)) 
         end do
      end do 

      title = 'Real E alpha (1/4)'
      
      call ezconc_sum(capr, capz, freal, ff, nnodex, nnodey, numb,
     .   nxmx, nymx, nlevmax, title, titx, tity, iflag)
     
      if (iflag .eq. 0) call boundary (capr, capz, rho, ff, nnodex,
     .   nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity)
     
     
*     --------------------------------------------------------
*     plot electric field in poloidal cross-section (one half)
*     -------------------------------------------------------- 

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ealpha_sum(i, j, nphi3d_half))
            fimag(i,j) = imag(ealpha_sum(i, j, nphi3d_half))
	    power(i,j) = real(redotj(i, j, nphi3d_half)) 
         end do
      end do 

      title = 'Real E alpha (1/2)'
      
      call ezconc_sum(capr, capz, freal, ff, nnodex, nnodey, numb,
     .   nxmx, nymx, nlevmax, title, titx, tity, iflag)
     
      if (iflag .eq. 0) call boundary (capr, capz, rho, ff, nnodex,
     .   nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity)
     
*     -----------------------------------------------------------
*     plot electric field in poloidal cross-section (three quarter)
*     ----------------------------------------------------------- 

      do i = 1, nnodex          
         do j = 1, nnodey
            freal(i,j) = real(ealpha_sum(i, j, nphi3d_3quarter))
            fimag(i,j) = imag(ealpha_sum(i, j, nphi3d_3quarter))
	    power(i,j) = real(redotj(i, j, nphi3d_3quarter)) 
         end do
      end do 

      title = 'Real E alpha (3/4)'
      
      call ezconc_sum(capr, capz, freal, ff, nnodex, nnodey, numb,
     .   nxmx, nymx, nlevmax, title, titx, tity, iflag)
     
      if (iflag .eq. 0) call boundary (capr, capz, rho, ff, nnodex,
     .   nnodey, numb, nxmx, nymx, nlevmax, title, titx, tity)  
     
     
*     --------------------------------------------------------------
*     plot E field in toroidal cross-section (rectangular - shifted)
*     --------------------------------------------------------------

      titx  = 'R (m)'
      tity  = 'phi (radians)'
     
      jhalf = nnodey / 2
      khalf = nphi3d / 2
      
      write(6, *) "jhalf = ", jhalf
      write(6, *) "nnodex = ", nnodex
      write(6, *) "nphi3d = ", nphi3d
      
             
      do i = 1, nnodex
         do k = 1, nphi3d
	 
	    if (k .lt. khalf) kshift = k + khalf
	    if (k .ge. khalf) kshift = k - khalf + 1
	    
c           freal(i, k) = real(ealpha_sum(i, jhalf, kshift))
c           fimag(i, k) = imag(ealpha_sum(i, jhalf, kshift))
c	    power(i, k) = real(redotj(i, jhalf, kshift))	    
c	    capr_plot(i, k) = capr(i)
	    
            freal_phi(i, k) = real(ealpha_sum(i, jhalf, kshift))
            fimag_phi(i, k) = imag(ealpha_sum(i, jhalf, kshift))
	    power_phi(i, k) = real(redotj(i, jhalf, kshift))	    
	    capr_plot_phi(i, k) = capr(i)
	    
	    
         end do
      end do
      
     

      
      title = 'Real E alpha'
      call ezconc_sum(capr, phi3d_shift, freal_phi, ff, nnodex, nphi3d,
     .           numb, nxmx, nphimx, nlevmax, title, titx, tity, iflag)
     
      title = 'Imag E alpha'
      call ezconc_sum(capr, phi3d_shift, fimag_phi, ff, nnodex, nphi3d,
     .           numb, nxmx, nphimx, nlevmax, title, titx, tity, iflag) 
     
      title = 'Real E dot J'
      call ezconc_sum(capr, phi3d_shift, power_phi, ff, nnodex, nphi3d,
     .           numb, nxmx, nphimx, nlevmax, title, titx, tity, iflag) 
     
*     --------------------
*     plot cut in R vs phi
*     --------------------      
      title = 'Mod(E_alpha)'
      titx  = 'phi'
      tityl = 'Mod(E_alpha)'
      tityr = ' '

      i = nnodex / 2
      j = nnodey / 2 
      khalf = nphi3d / 2   
      
      do k = 1, nphi3d	 
	 if (k .lt. khalf) kshift = k + khalf
	 if (k .ge. khalf) kshift = k - khalf + 1	    
         fmod(k) = real(ealpha_sum_mod(i, j, kshift))	    
      end do      
     
      call ezplot1_sum(title, titx, tityl, tityr, phi3d_shift, fmod, 
     .    nphi3d, nphimx)
     
          
     
*     --------------------------------------
*     plot in toroidal cross-section (polar)
*     --------------------------------------    
     
      titx  = ' '
      tity  = ' '
      numb = 10
            
c      write(6, *)"xplotm(i) = "
c      write (6, 310) (xplotm(i),  i = 1, nxplot) 
      
      title = 'Real E alpha'           
      call ezconc_sum(xplotm, yplotm, frealp, ff, nxplot, nyplot, numb,
     .   nxplot_dim, nyplot_dim, nlevmax, title, titx, tity, iflag) 
     
     
              
      if (iflag .eq. 0) call boundary_tor(xplotm, yplotm, capr_plotp, 
     .   ff, nxplot, nyplot, numb,
     .   nxplot_dim, nyplot_dim, nlevmax, title, titx, tity, iflag,
     .   capr_min, capr_max)
               
     
      title = 'Imag E alpha'           
      call ezconc_sum(xplotm, yplotm, fimagp, ff, nxplot, nyplot, numb,
     .   nxplot_dim, nyplot_dim, nlevmax, title, titx, tity, iflag) 
              
      if (iflag .eq. 0) call boundary_tor(xplotm, yplotm, capr_plotp, 
     .   ff, nxplot, nyplot, numb,
     .   nxplot_dim, nyplot_dim, nlevmax, title, titx, tity, iflag,
     .   capr_min, capr_max)
     
     
      title = 'Real EdotJ'           
      call ezconc_sum(xplotm, yplotm, powerp, ff, nxplot, nyplot, numb,
     .   nxplot_dim, nyplot_dim, nlevmax, title, titx, tity, iflag) 
              
      if (iflag .eq. 0) call boundary_tor(xplotm, yplotm, capr_plotp, 
     .   ff, nxplot, nyplot, numb,
     .   nxplot_dim, nyplot_dim, nlevmax, title, titx, tity, iflag,
     .   capr_min, capr_max)
     
     
     
*     -------------------------
*     Plot the SWIM output file:
*     ------------------------- 

      title= 'Volume element, dvol'
      titll= 'dvol (m3)'
      titlr='       '
      titlb= 'rho'
	    
      call ezplot1_sum(title, titlb, titll, titlr, rhon, dvol, 
     .    nnoderho, nrhomax)
     
    
      title= 'Integrated redotj'
      titll= 'P (watts)'
      titlr= '       '
      titlb= 'rho'
     
      call ezplot6(title, titll, titlr, titlb, rhon, redotji2_int,
     .    redotje_int, redotji1_int, redotji3_int, redotji4_int, 
     .    redotji5_int, nnoderho, nrhomax) 
     
      title= 'Integrated wdot'
      titll= 'P (watts)'
      titlr= '       '
      titlb= 'rho'
     
      call ezplot6(title, titll, titlr, titlb, rhon, wdoti2_int,
     .    wdote_int, wdoti1_int, wdoti3_int, wdoti4_int, wdoti5_int,
     .    nnoderho, nrhomax)
     
     
      title= 'Real(EdotJ) * dvol'
      titll= 'P(watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot6(title, titll, titlr, titlb, rhon, redotj1_dvol, 
     .   redotj2_dvol, redotj3_dvol, redotj4_dvol, 
     .    redotj5_dvol, redotje_dvol, nnoderho, nrhomax)  
     
     
      title= 'Wdot * dvol'
      titll= 'P (watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot6(title, titll, titlr, titlb, rhon, 
     .    wdoti1_dvol, wdoti2_dvol, wdoti3_dvol, wdoti4_dvol, 
     .    wdoti5_dvol, wdote_dvol, nnoderho, nrhomax)              
      
             
 6000 continue                  

*     -------------------------
*     Close the graphics device
*     -------------------------		  
      call pgclos 
      
      if(ndisti2  .eq. 1) then      
         deallocate( UPERP )
         deallocate( UPARA )
      
         deallocate( bqlsum_i2 )
         deallocate( cqlsum_i2 )
         deallocate( eqlsum_i2 )
         deallocate( fqlsum_i2 )
      
         deallocate( bqlsum_i2_2d ) 
	 deallocate( cqlsum_i2_2d )  
         deallocate( eqlsum_i2_2d )  
	 deallocate( fqlsum_i2_2d )       
      end if
            
      close (242)
      close (342)      
      close (362)      
      

      deallocate( ealpha_sum_mod ) 
      deallocate( ebeta_sum_mod ) 
      deallocate (ealpha_sum)  
      deallocate( ebeta_sum )
      deallocate( eb_sum )
      deallocate( redotj ) 
      
      deallocate( ntilda )        
      
      deallocate( frealp )
      deallocate( fimagp )
      deallocate( powerp )
      deallocate( capr_plotp ) 
      
      deallocate( freal )
      deallocate( fimag )
      deallocate( capr_plot )
      deallocate( power ) 
      
      deallocate( freal_phi )
      deallocate( fimag_phi )
      deallocate( capr_plot_phi )
      deallocate( power_phi )                     
              

  309 format(10i10) 
 3310 format(1p6e18.10)     
11001 format(1e16.8,2e15.8)
 1312 format(i10,1p8e12.4)
99998 format(1h ,f12.4, 1p9e12.4)
99997 format(i10)
88887 format(2i10)
  899 format(1h ,5x,12h      cpu = ,1pe12.4,4h min/
     1       1h ,5x,12h       io = ,1pe12.4,4h min/
     1       1h ,5x,12h      sys = ,1pe12.4,4h min/
     1       1h ,5x,12h      mem = ,1pe12.4,4h min/
     1       1h ,5x,12h    total = ,1pe12.4,4h min)
 2000 format(1i10,1e10.3,1i10,1e10.3)
 2001 format(1i10,1p9e12.4)
 2517 format(1i10,1e10.3,1i10,1e10.3)
 2211 format(8i10)
 2220 format(3i10,3e10.3,2i10)
 2230 format(2i10,7e10.3)
 2240 format(1i10,7e10.3)
 2260 format(1i10,1e10.3,1i10,1e10.3)
 2280 format(2i10,5e10.3,1i10)
 1009 format(3x,21h   power absorbed  = ,1pe12.4,9h watts   )
  300 format(i10/i10,4e10.3,i10)
  310 format(1p6e12.4)
  311 format(1p10e12.4)
  122 format(1i10,1p9e12.4)
  123 format(2i10,1p9e12.4)  
 1051 format(6i10)
 1052 format(1i10,6e10.4)
 1053 format(6e10.4)
 1013 format(3x,15h         lt  = ,1pe12.4,9h m       )
 1014 format(3x,15h         wd  = ,1pe12.4,9h m       )
 1015 format(3x,15h       dphi  = ,1pe12.4,9h radians )
 3015 format(3x,15h       phi0  = ,1pe12.4,9h radians )
 3016 format(3x,15h        xi0  = ,1pe12.4,9h Amps    )
c     skip a line
  163 format(1h0)
 1163 format(1h ,5x,6h  k   ,2x,6h phi  ,6x,6hre cur,4x,8h im cur ,
     1       4x,8h re jsum,4x,8h im jsum,
     1       4x,8h re csum,4x,8h im csum,
     1       4x,8h        ,4x,8h        )
 1164 format(1h ,5x,6hnphi  ,2x,6h  k// ,6x,6h pabs ,5x,8h weight ,
     1       5x,8h re pc  ,4x,8h im pc  ,
     1       4x,8h pabs wt,4x,8h gammacd,
     1       4x,8h  damp1 ,4x,8h        )
 1165 format(1h ,5x,6h  k   ,2x,6h phi  ,6x,6hre cur,4x,8h im cur ,
     1       4x,8hre jthsm,4x,8him jthsm,
     1       4x,8h re epsi,4x,8h im epsi,
     1       4x,8h        ,4x,8h        )
  162 format(1h1)
 2011 format(3x,33htotal fast wave current driven = ,1pe12.4,8h Amps   )
 2016 format(3x,22htotal ohmic current = ,1pe12.4,9h Amps    )
 2017 format(3x,27htotal boot strap current = ,1pe12.4,9h Amps    )
 7014 format(3x,26hcurrent driven per watt = ,1pe12.4,10h Amps/watt)
 7015 format(3x,27hcurrent drive efficiency = ,1pe12.4,
     1   15h Amps/watt/m**2)
 7016 format(3x,31h average single pass damping = ,1pe12.4,
     1   15h               )
 2009 format(3x,30hpower absorbed by electrons = ,1pe12.4,9h watts   )
 2012 format(3x,27hpower absorbed by alphas = ,1pe12.4,9h watts   )
 2013 format(3x,29hpower absorbed by majority = ,1pe12.4,9h watts   )
 2014 format(3x,29hpower absorbed by minority = ,1pe12.4,9h watts   )
 2015 format(3x,29hpower absorbed by impurity = ,1pe12.4,9h watts   )
 1008 format(3x,21h      phi antenna  = ,1pe12.4,9h pi      )
 1010 format(3x,42h   real power transferred from antenna  = ,
     1   1pe12.4,9h watts   /
     1      ,3x,47h   imaginary power transferred from antenna  = ,
     1   1pe12.4,9h watts   )
 1011 format(3x,46hreal part of antenna impedance (resistance) = ,
     1   1pe12.4,9h ohms    /
     1      ,3x,50himaginary part of antenna impedance (reactance) = ,
     1   1pe12.4,9h ohms    )
 1016 format(3x,33htotal power coupled at antenna = ,
     1   1pe12.4,9h watts   )
 1002 format(2i10,7e10.3)
 2165 format(3i10,5e10.3)
 1101 format(2i10,5e10.3,1i10)
 1001 format(8e10.3)
 1000 format(1i10,7e10.3)
11009 format(3x,21h        frequency  = ,1pe12.4,7h hertz )
81822 format(3x,21h            xwall  = ,1pe12.4,7h m     )
81823 format(3x,21h           xnwall  = ,1pe12.4,7h m-3   )
81824 format(3x,21h             amu1  = ,1pe12.4,7h amu   )
81825 format(3x,21h             amu2  = ,1pe12.4,7h amu   )
81826 format(3x,21h             amu3  = ,1pe12.4,7h amu   )
98824 format(3x,21h             amu4  = ,1pe12.4,7h amu   )
98825 format(3x,21h             amu5  = ,1pe12.4,7h amu   )
98826 format(3x,21h             amu6  = ,1pe12.4,7h amu   )
98827 format(3x,21h             amu7  = ,1pe12.4,7h amu   )
81827 format(3x,21h               z1  = ,1pe12.4,9h q proton)
81828 format(3x,21h               z2  = ,1pe12.4,9h q proton)
81829 format(3x,21h               z3  = ,1pe12.4,9h q proton)
98834 format(3x,21h               z4  = ,1pe12.4,9h q proton)
98835 format(3x,21h               z5  = ,1pe12.4,9h q proton)
98836 format(3x,21h               z6  = ,1pe12.4,9h q proton)
98837 format(3x,21h               z7  = ,1pe12.4,9h q proton)
81830 format(3x,21h              eta  = ,1pe12.4,7h       )
81831 format(3x,21h             eta3  = ,1pe12.4,7h       )
91834 format(3x,21h             eta4  = ,1pe12.4,7h       )
91835 format(3x,21h             eta5  = ,1pe12.4,7h       )
91836 format(3x,21h             eta6  = ,1pe12.4,7h       )
91837 format(3x,21h             eta7  = ,1pe12.4,7h       )
81832 format(3x,21h            xnlim  = ,1pe12.4,7h m-3   )
81833 format(3x,21h            qavg0  = ,1pe12.4,7h       )
81834 format(3x,21h           scrape  = ,1pe12.4,7h m     )
81835 format(3x,21h           alphac  = ,1pe12.4,7h       )
81836 format(3x,21h           alphan  = ,1pe12.4,7h       )
81837 format(3x,21h           alphat  = ,1pe12.4,7h       )
81838 format(3x,21h           zeffcd  = ,1pe12.4,9h q proton)
81839 format(3x,21h            dpsi0  = ,1pe12.4,7h       )
81840 format(3x,21h            telim  = ,1pe12.4,7h eV    )
81841 format(3x,21h            tilim  = ,1pe12.4,7h eV    )
91841 format(3x,21h           ti2lim  = ,1pe12.4,7h eV    )
91842 format(3x,21h           ti3lim  = ,1pe12.4,7h eV    )
91844 format(3x,21h           ti4lim  = ,1pe12.4,7h eV    )
91845 format(3x,21h           ti5lim  = ,1pe12.4,7h eV    )
91846 format(3x,21h           ti6lim  = ,1pe12.4,7h eV    )
91847 format(3x,21h           ti7lim  = ,1pe12.4,7h eV    )
81842 format(3x,21h           deltap  = ,1pe12.4,7h       )
81843 format(3x,21h                p  = ,1pe12.4,7h       )
81845 format(3x,21h            q07qa  = ,1pe12.4,7h       )
 1012 format(3x,21h             omgrf = ,1pe12.4,7h hertz )
 1714 format(3x,21h               xk0 = ,1pe12.4,7h m-1   )
 1021 format(3x,21h                nz = ,1pe12.4,7h       )
 1921 format(3x,21h           w phase = ,1pe12.4,7h       )
 1812 format(3x,21h                rt = ,1pe12.4,7h m     )
 1835 format(3x,21h            ekappa = ,1pe12.4,7h       )
 1822 format(3x,21h            aplasm = ,1pe12.4,7h m     )
 1823 format(3x,21h              xant = ,1pe12.4,7h m     )
 1809 format(3x,21h                b0 = ,1pe12.4,7h T     )
81844 format(3x,21h              xne0 = ,1pe12.4,7h m-3   )
 1813 format(3x,21h              xn10 = ,1pe12.4,7h m-3   )
 1814 format(3x,21h              xn20 = ,1pe12.4,7h m-3   )
 1834 format(3x,21h              xn30 = ,1pe12.4,7h m-3   )
 9834 format(3x,21h              xn40 = ,1pe12.4,7h m-3   )
 9835 format(3x,21h              xn50 = ,1pe12.4,7h m-3   )
 9836 format(3x,21h              xn60 = ,1pe12.4,7h m-3   )
 9837 format(3x,21h              xn70 = ,1pe12.4,7h m-3   )
 8834 format(3x,21h              xn40 = ,1pe12.4,7h m-3   )
 8835 format(3x,21h              xn50 = ,1pe12.5,7h m-3   )
 8836 format(3x,21h              xn60 = ,1pe12.6,7h m-3   )
 8837 format(3x,21h              xn70 = ,1pe12.7,7h m-3   )
 1815 format(3x,21h               te0 = ,1pe12.4,7h eV    )
 1821 format(3x,21h               ti0 = ,1pe12.4,7h eV    )
 1922 format(3x,21h              ti02 = ,1pe12.4,7h eV    )
 1923 format(3x,21h              ti03 = ,1pe12.4,7h eV    )
 1924 format(3x,21h              ti04 = ,1pe12.4,7h eV    )
 1925 format(3x,21h              ti05 = ,1pe12.4,7h eV    )
 1926 format(3x,21h              ti06 = ,1pe12.4,7h eV    )
 1927 format(3x,21h              ti07 = ,1pe12.4,7h eV    )
 1017 format(3x,21h        xnu1/omgrf = ,1pe12.4,7h       )
 1018 format(3x,21h        xnu2/omgrf = ,1pe12.4,7h       )
 1019 format(3x,21h        xnu3/omgrf = ,1pe12.4,7h       )
 1054 format(3x,21h           flandau = ,1pe12.4,7h       )
 1022 format(3x,21h                nl = ,i12,7h       )
71022 format(3x,21h                nr = ,i12,7h       )
71023 format(3x,21h            ntheta = ,i12,7h       )
71024 format(3x,21h             mdiag = ,i12,7h       )
71034 format(3x,21h            mpdiag = ,i12,7h       )
73031 format(3x,21h             jharm = ,i12,7h       )
73032 format(3x,21h              lmax = ,i12,7h       )
73033 format(3x,21h               iez = ,i12,7h       )
73034 format(3x,21h             isort = ,i12,7h       )
71025 format(3x,21h             ndiag = ,i12,7h       )
71026 format(3x,21h             ikprl = ,i12,7h       )
71027 format(3x,21h             igeom = ,i12,7h       )
71028 format(3x,21h             nphi1 = ,i12,7h       )
71029 format(3x,21h             nphi2 = ,i12,7h       )
71030 format(3x,21h             ncoll = ,i12,7h       )
91031 format(3x,21h            ncoll7 = ,i12,7h       )
71031 format(3x,21h             imode = ,i12,7h       )
31012 format(3x,34h    power absorbed by electrons = ,
     1   1pe12.4,9h  %      )
31013 format(3x,34h    power absorbed by majority  = ,
     1   1pe12.4,9h  %      )
31014 format(3x,34h    power absorbed by minority  = ,
     1   1pe12.4,9h  %      )
31015 format(3x,34h    power absorbed by species 3 = ,
     1   1pe12.4,9h  %      )
31016 format(3x,34h    power absorbed by species 4 = ,
     1   1pe12.4,9h  %      )
31017 format(3x,34h    power absorbed by species 5 = ,
     1   1pe12.4,9h  %      )
31018 format(3x,34h    power absorbed by species 6 = ,
     1   1pe12.4,9h  %      )
31019 format(3x,34h    power absorbed by species 7 = ,
     1   1pe12.4,9h  %      )
     
 5000 continue
 
      return
      end      

c
c*********************************************************************
c
      subroutine boundary_tor(r, theta, f, flevel, nr, nth, nlevel,
     .   nrmax, nthmax, nlevmax, title, titx, tity, iflag, 
     .   capr_min, capr_max)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1, norange

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, theta, r, flevel, f, xmin, ymin, xmax 
     
      realcapr_min, capr_max
      
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1    ncolelec, iflag, nlevelb
     
      real tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8
      
      dx = r(2) - r(1)
      dy = theta(2) - theta(1)

      tr(1) = r(1) - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

c--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)
      
c      write(6, *)"fmax = ", fmax, "   fmin = ", fmin
c      write(16,*)"fmax = ", fmax, "   fmin = ", fmin

      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
c         write(6, *)"fmax = ", fmax
c	  write(6, *)"fmin = ", fmin
         iflag = 1
         return
      end if

      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
c        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel.eq.1)flevel(1)=1.0e-03

       nlevelb = 2
       
       flevel(1) = capr_min
       flevel(2) = capr_max

      

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

      
c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps

      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0
      
c-- get ready to plot

      call pgsci(nblack)
c      call pgenv(xmin, xmax, ymin, ymax, 1, 0)
      CALL PGVSTD
      CALL PGWNAD (xmin, xmax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
      
c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

c	 call pgsci(nyellow)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &       nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then
      
c	 call pgsci(nblue)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &       nlevgt, tr)
      endif

      call pgsci(nblack)
C      call pglab(titx, tity, title)
      
  310 format(1p6e12.4)
  312 format(i10, 1p6e12.4)

      
      return
      end
      
c
c*********************************************************************
c
     

      subroutine ezplot6(title, titll, titlr, titlb, x1, y1,
     .   y2, y3, y4, y5, y6,
     .   nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax),
     .     y5(nrmax), y6(nrmax)
      real y1max, y2max, y3max, y4max, y5max, y6max
      real y1min, y2min, y3min, y4min, y5min, y6min
      real ymin, ymax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb
      
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
     
      integer nblack,nred,nyellow, ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)
      call a1mnmx(y6, nrmax, nr, y6min, y6max)

      ymax = amax1(y1max, y2max, y3max, y4max, y5max, y6max)
      ymin = amin1(y1min, y2min, y3min, y4min, y5min, y6min)
c      ymin=0.0

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1
      
c      ymax = 2.8e+06
      
c Advance plotter to a new page, define coordinate range of graph and draw axes
c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)      
      
c Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)
      
c Plot the line graph.

      call pgsci(nblue)      
      call pgline(nr, x1, y1)

      call pgsci(nred)       
      call pgline(nr, x1, y2)
      
      call pgsci(ncyan) 
      call pgline(nr, x1, y3)
      
      call pgsci(ngreen) 
      call pgline(nr, x1, y4)
      
      call pgsci(nyellow) 
      call pgline(nr, x1, y5)
      
      call pgsci(nmagenta) 
      call pgline(nr, x1, y6)
      
      call pgsci(nblack) 
	    
  300 format (1p9e11.3)

      return
      end

c
c***************************************************************************
c

      subroutine ezplot1_sum(title, titx, titll, titlr, x1, y1, nr, 
     .    nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titx
      character*32 titll
      character*32 titlr
      
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
     
      integer nblack, nred, nyellow, ngreen, naqua, npink,
     1   nwheat, ngrey, nbrown, nblue, nblueviolet, ncyan1,
     1   nturquoise, nmagenta, nsalmon, nwhite, norange, ncolln3
     
     
      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8
     
      ncolln3 = ngreen

      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = x1(nr)
      rhomin = x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1
      
c Advance plotter to a new page, define coordinate range of graph and draw axes

c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)
      
      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
      
     
c Label the axes (note use of \u and \d for raising exponent).

c      call pglab(titx, titll, title)
      call pgmtxt('b', 3.2, 0.5, 0.5, titx)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)

      CALL PGSCI(nblue)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)
      
c Plot the line graph.

      
      call pgline(nr, x1, y1)
      
      CALL PGSCI(nblack)
	    
  300 format (1p9e11.3)
  
      return
      end
      
c
c***************************************************************************
c


      subroutine ezplot1_0_sum(title, titx, titll, titlr, x1, y1, 
     .   nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titx
      character*32 titll
      character*32 titlr
      
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
     
      integer nblack, nred, nyellow, ngreen, naqua, npink,
     1   nwheat, ngrey, nbrown, nblue, nblueviolet, ncyan1,
     1   nturquoise, nmagenta, nsalmon, nwhite, norange, ncolln3
     
     
      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8
     
      ncolln3 = ngreen

      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = x1(nr)
      rhomin = x1(1)
      ymax = ymax * 1.1
      ymin = 0.0

      if (ymin .le. 0.0) ymin = ymin * 1.1
      
c Advance plotter to a new page, define coordinate range of graph and draw axes

c      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)
      
      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
      
     
c Label the axes (note use of \u and \d for raising exponent).

c      call pglab(titx, titll, title)
      call pgmtxt('b', 3.2, 0.5, 0.5, titx)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)

      CALL PGSCI(nblue)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)
      
c Plot the line graph.

      
      call pgline(nr, x1, y1)
      
      CALL PGSCI(nblack)
	    
  300 format (1p9e11.3)
  
      return
      end
      
c
c***************************************************************************
c

      subroutine ezplot2_sum(title, titll, titlr, titlb, 
     .                                       x1, y1, y2, nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y2(nrmax)
      real y1max, y2max, y3max, y1min, y2min, y3min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb
      
      integer nr, nrmax, n

      integer nplot1,ncollab, ncolion,ncolbox, ncyan,
     1    ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink,
     1   nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1,
     1   nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange
     
      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)

c      write(6, 300) y2min, y2max
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin
      
c      write(*,*)"ymin = ", ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BNST', 0.0, 0)
      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)

      CALL PGSCI(nblue)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)

      
      call pgline(nr, x1, y1)  
      
c      write(6, *)y2      

      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      y2min = y2min 
      y2max = y2max * 1.1
      
c      write(6, *) y2min, y2max
c      call exit


      CALL PGSWIN (rhomin, rhomax, y2min, y2max)
      CALL PGSCI(nblack)
      CALL PGBOX  (' ', 0.0, 0, 'CMST', 0.0, 0)
      
c      write(*,*) "y2 = ", y2
      
      CALL PGSCI(ngreen)
      call pgline(nr, x1, y2)  
      call PGMTXT ('r', 2.0, 0.5, 0.5, titlr)

      CALL PGSCI(nblack)
      

  300 format (1p9e11.3)
      return
      end
c
c***************************************************************************
c

      subroutine ezconc_sum(r, theta, f, flevel, nr, nth, nlevel,
     1   nrmax, nthmax, nlevmax, title, titx, tity, iflag)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua,
     1   npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen,
     1   nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5,
     1   ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2,
     1   ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant,
     1   nndm1, norange

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax,
     1   ypmin, theta, r, flevel, f, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr,
     1   nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin,
     1    ncolelec, iflag
     
      real tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8
      
      dx = r(2) - r(1)
      dy = theta(2) - theta(1)

      tr(1) = r(1) - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)
      
c      write(6, *) "xmin = ", xmin
c      write(6, *) "xmax = ", xmax
c      write(6, *) "ymin = ", ymin
c      write(6, *) "ymax = ", ymax

c--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)
      
c      write(6, *)"fmax = ", fmax, "   fmin = ", fmin
c      write(16,*)"fmax = ", fmax, "   fmin = ", fmin

      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
c         write(6, *)"fmax = ", fmax
c	  write(6, *)"fmin = ", fmin
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

      
c--Scale window to user coordinates
c--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps

      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0
      
c--Advance graphics frame and get ready to plot

      call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 1, 0)
      
c Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

	 call pgsci(nyellow)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
     &       nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then
      
	 call pgsci(nblue)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
     &       nlevgt, tr)
      endif

      call pgsci(nblack)
      call pglab(titx, tity, title)
      
  310 format(1p6e12.4)
  312 format(i10, 1p6e12.4)

      
      return
      end
c
c*********************************************************************
c
